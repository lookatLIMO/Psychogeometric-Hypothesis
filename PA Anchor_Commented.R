# ============================================================
# MONTE CARLO POWER SIMULATION FOR A LINEAR MIXED MODEL
# IN THE ANCHOR DESIGN
# ------------------------------------------------------------
#
# Purpose:
# This script estimates statistical power for detecting a focal
# pairwise-difference predictor (for example, dH) in a linear
# mixed-effects model with a standardized singular MDS dimension
# as the outcome.
#
# Standardization and calibration strategy:
# - All pairwise-difference predictors are z-standardized within each replicate.
# - The outcome is z-standardized within each replicate.
# - The random-intercept standard deviations and residual standard deviation
#   are fixed to 1 in the data-generating model.
# - For each n_targets condition, the script first estimates the expected
#   non-focal variance of the outcome and then calibrates the raw focal
#   coefficient so that beta_true corresponds approximately to the target
#   standardized beta in the fitted model.
#
# Simulation workflow:
# 1) Build all unordered target pairs for each sample size.
# 2) Precompute calibration constants for each n_targets value.
# 3) Simulate one dataset at a time under the requested focal effect.
# 4) Fit the mixed model and extract the focal estimate and p-value.
# 5) Repeat adaptively until either the Monte Carlo precision target is met
#    or the maximum number of simulations is reached.
# 6) Save checkpoints after each design condition and summarize the minimum
#    detectable standardized effect size.
# ============================================================

## --------------------------------------------------
## 0. Packages
## --------------------------------------------------

# Install lme4 if it is not already installed.
if (!requireNamespace("lme4", quietly = TRUE)) install.packages("lme4")

# Install lmerTest if it is not already installed.
if (!requireNamespace("lmerTest", quietly = TRUE)) install.packages("lmerTest")

# Load lme4 for fitting linear mixed-effects models.
library(lme4)

# Load lmerTest so p-values are available in model summaries.
library(lmerTest)

# Set a global random seed for reproducibility.
set.seed(123)

## --------------------------------------------------
## 1. Global settings
## --------------------------------------------------

# Number of simulation replications to run in each batch.
batch_size <- 1

# Minimum number of usable simulations required before early stopping is allowed.
min_nsim <- 100

# Maximum number of simulations allowed for each design condition.
max_nsim <- 1000

# Target Monte Carlo standard error for the estimated power.
mc_se_target <- 0.01

# Grid of sample sizes to evaluate.
n_grid <- c(1505)

# Grid of target standardized beta values to evaluate.
beta_grid <- seq(0.001, 0.005, by = 0.001)

# Desired target power.
target_power <- 0.80

# Familywise alpha level.
alpha_family <- 0.05

# Number of tests used in the Bonferroni correction.
n_tests <- 18

# Bonferroni-corrected alpha threshold.
alpha_bonf <- alpha_family / n_tests

# Name of the focal predictor.
focal_trait <- "dH"

# Set all nuisance fixed effects to zero by default.
beta_nuisance <- c(
  dH = 0.00,
  dE = 0.00,
  dX = 0.00,
  dA = 0.00,
  dC = 0.00,
  dO = 0.00
)

# Standard deviation of the random intercept for subject i.
sd_subj_i <- 1.00

# Standard deviation of the random intercept for subject j.
sd_subj_j <- 1.00

# Residual standard deviation.
sd_resid <- 1.00

# File path for the main checkpoint.
checkpoint_file <- "power_results_checkpoint.rds"

# File path for the backup checkpoint.
backup_file <- "power_results_backup.rds"

## --- Parallel settings for Mac ---

# Detect the number of available physical CPU cores.
available_cores <- parallel::detectCores(logical = FALSE)

# Fall back to 2 cores if detection fails.
if (is.na(available_cores)) available_cores <- 2L

# Number of cores to use.
n_cores <- 7L

## --- Faster lmer control ---

# Use a faster control configuration to reduce model-fitting overhead.
fit_control <- lmerControl(
  calc.derivs = FALSE,
  optimizer = "nloptwrap",
  optCtrl = list(algorithm = "NLOPT_LN_BOBYQA")
)

## --------------------------------------------------
## 2. Helper functions
## --------------------------------------------------

# Create all unordered pairs of targets.
make_pairs <- function(n_targets) {
  cmb <- utils::combn(n_targets, 2)
  list(
    subj_i = cmb[1, ],
    subj_j = cmb[2, ]
  )
}

# Build and store all pair-related objects for one sample size.
build_pair_object <- function(n_targets) {
  p <- make_pairs(n_targets)
  p$subj_i_f <- factor(p$subj_i, levels = seq_len(n_targets))
  p$subj_j_f <- factor(p$subj_j, levels = seq_len(n_targets))
  p$n_pairs <- length(p$subj_i)
  p
}

# Safely z-standardize a numeric vector.
# If the vector has zero or missing standard deviation, return zeros instead.
scale_safe <- function(x) {
  s <- sd(x)
  if (is.na(s) || s == 0) {
    rep(0, length(x))
  } else {
    (x - mean(x)) / s
  }
}

# Safely compute the mean of a vector.
safe_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

# Compute power-related summary statistics.
compute_power_stats <- function(n_detected,
                                usable,
                                mc_se_target = 0.01,
                                exact_ci = FALSE) {
  if (usable == 0) {
    return(list(
      n_detected = 0L,
      power = NA_real_,
      mc_se_power = NA_real_,
      ci_low_power = NA_real_,
      ci_high_power = NA_real_,
      precision_ok = NA
    ))
  }
  
  # Estimate power as the proportion of successful detections.
  power_hat <- n_detected / usable
  
  # Compute the Monte Carlo standard error of the power estimate.
  mc_se <- if (usable >= 2) {
    sqrt(power_hat * (1 - power_hat) / usable)
  } else {
    NA_real_
  }
  
  # Optionally compute an exact 95% binomial confidence interval.
  if (exact_ci) {
    bt <- binom.test(n_detected, usable, conf.level = 0.95)
    ci_low <- unname(bt$conf.int[1])
    ci_high <- unname(bt$conf.int[2])
  } else {
    ci_low <- NA_real_
    ci_high <- NA_real_
  }
  
  # Flag whether the Monte Carlo precision target has been reached.
  precision_ok <- !is.na(mc_se) && mc_se < mc_se_target
  
  list(
    n_detected = n_detected,
    power = power_hat,
    mc_se_power = mc_se,
    ci_low_power = ci_low,
    ci_high_power = ci_high,
    precision_ok = precision_ok
  )
}

# Build a one-line progress bar string with elapsed time and ETA.
progress_line_adaptive <- function(current,
                                   max_total,
                                   elapsed_sec,
                                   eta_sec,
                                   prefix = "",
                                   power = NA_real_,
                                   mc_se = NA_real_,
                                   width = 28) {
  frac <- if (max_total > 0) current / max_total else 0
  frac <- min(max(frac, 0), 1)
  done <- floor(width * frac)
  left <- width - done
  
  bar <- paste0(
    "[",
    paste(rep("=", done), collapse = ""),
    paste(rep("-", left), collapse = ""),
    "]"
  )
  
  sprintf(
    "\r%s %s %3d%% | sims: %d/%d | power: %s | MCSE: %s | elapsed: %.1f min | ETA: %.1f min",
    prefix,
    bar,
    round(frac * 100),
    current,
    max_total,
    ifelse(is.na(power), "NA", sprintf("%.3f", power)),
    ifelse(is.na(mc_se), "NA", sprintf("%.3f", mc_se)),
    elapsed_sec / 60,
    eta_sec / 60
  )
}

# Build standardized pairwise predictors for one replicate.
# Each target receives ordinal values from 1 to 7 on the six HEXACO traits.
# Pairwise predictors are absolute differences, then z-standardized.
build_standardized_predictors <- function(pair_obj, n_targets) {
  H <- sample.int(7, n_targets, replace = TRUE)
  E <- sample.int(7, n_targets, replace = TRUE)
  X <- sample.int(7, n_targets, replace = TRUE)
  A <- sample.int(7, n_targets, replace = TRUE)
  C <- sample.int(7, n_targets, replace = TRUE)
  O <- sample.int(7, n_targets, replace = TRUE)
  
  data.frame(
    dH = scale_safe(abs(H[pair_obj$subj_i] - H[pair_obj$subj_j])),
    dE = scale_safe(abs(E[pair_obj$subj_i] - E[pair_obj$subj_j])),
    dX = scale_safe(abs(X[pair_obj$subj_i] - X[pair_obj$subj_j])),
    dA = scale_safe(abs(A[pair_obj$subj_i] - A[pair_obj$subj_j])),
    dC = scale_safe(abs(C[pair_obj$subj_i] - C[pair_obj$subj_j])),
    dO = scale_safe(abs(O[pair_obj$subj_i] - O[pair_obj$subj_j]))
  )
}

# Estimate the expected variance components used for calibration.
# var_x stores the expected variance of the focal standardized predictor.
# var_other stores the expected variance of everything else in the outcome:
# nuisance fixed effects, random intercepts, and residual noise.
estimate_calibration_constants <- function(pair_obj,
                                           n_targets,
                                           focal_trait,
                                           beta_nuisance,
                                           sd_subj_i,
                                           sd_subj_j,
                                           sd_resid,
                                           n_calib = 200,
                                           seed = 999) {
  set.seed(seed + n_targets)
  
  vx <- numeric(n_calib)
  v_other <- numeric(n_calib)
  
  for (r in seq_len(n_calib)) {
    pred <- build_standardized_predictors(pair_obj, n_targets)
    
    # Set the focal coefficient to zero during calibration.
    beta_vec <- beta_nuisance
    beta_vec[focal_trait] <- 0
    
    # Compute the fixed part excluding the focal effect.
    fixed_other <- beta_vec["dH"] * pred$dH +
      beta_vec["dE"] * pred$dE +
      beta_vec["dX"] * pred$dX +
      beta_vec["dA"] * pred$dA +
      beta_vec["dC"] * pred$dC +
      beta_vec["dO"] * pred$dO
    
    # Draw random intercepts for the two members of each pair.
    b_i <- rnorm(n_targets, mean = 0, sd = sd_subj_i)
    b_j <- rnorm(n_targets, mean = 0, sd = sd_subj_j)
    
    # Map random intercepts to the pair rows.
    rand_part <- b_i[pair_obj$subj_i] + b_j[pair_obj$subj_j]
    
    # Draw residual noise.
    resid_part <- rnorm(pair_obj$n_pairs, mean = 0, sd = sd_resid)
    
    # Combine all non-focal components.
    other_part <- fixed_other + rand_part + resid_part
    
    # Store the focal predictor variance and the non-focal outcome variance.
    vx[r] <- var(pred[[focal_trait]])
    v_other[r] <- var(other_part)
  }
  
  list(
    var_x = mean(vx, na.rm = TRUE),
    var_other = mean(v_other, na.rm = TRUE),
    n_calib = n_calib
  )
}

# Convert a target standardized beta into the raw coefficient needed
# in the data-generating model.
#
# If:
# beta_std = b_raw * sd_x / sqrt(b_raw^2 * var_x + var_other)
#
# then solving for b_raw gives:
# b_raw = beta_std * sqrt(var_other) / sqrt(var_x * (1 - beta_std^2))
raw_from_target_beta <- function(beta_std, var_x, var_other) {
  if (is.na(beta_std) || is.na(var_x) || is.na(var_other)) return(NA_real_)
  if (beta_std <= 0 || beta_std >= 1) stop("beta_std must be in (0, 1).")
  beta_std * sqrt(var_other) / sqrt(var_x * (1 - beta_std^2))
}

## --------------------------------------------------
## 3. Build the design grid and precompute pair structures
## --------------------------------------------------

# Build the full grid of simulation conditions.
design_grid <- expand.grid(
  n_targets = n_grid,
  beta_true = beta_grid,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

# Extract and sort the unique n_targets values.
unique_n <- sort(unique(design_grid$n_targets))

# Precompute all pair structures for each n_targets value.
pair_cache <- setNames(
  lapply(unique_n, build_pair_object),
  as.character(unique_n)
)

# Precompute calibration constants for each n_targets value.
calibration_cache <- setNames(
  lapply(unique_n, function(nn) {
    estimate_calibration_constants(
      pair_obj = pair_cache[[as.character(nn)]],
      n_targets = nn,
      focal_trait = focal_trait,
      beta_nuisance = beta_nuisance,
      sd_subj_i = sd_subj_i,
      sd_subj_j = sd_subj_j,
      sd_resid = sd_resid,
      n_calib = 200,
      seed = 999
    )
  }),
  as.character(unique_n)
)

# Print a message confirming that the precomputation step has finished.
cat("Precomputation completed for n_targets values:",
    paste(unique_n, collapse = ", "), "\n")

# Print the number of CPU cores used.
cat("Parallel cores used:", n_cores, "\n")

## --------------------------------------------------
## 4. Simulate one dataset
## --------------------------------------------------

# Simulate one dataset for a single replication.
simulate_one_dataset_fast <- function(pair_obj,
                                      n_targets,
                                      beta_focal_std,
                                      focal_trait,
                                      beta_nuisance,
                                      sd_subj_i,
                                      sd_subj_j,
                                      sd_resid,
                                      calibration_obj) {
  
  # Generate the standardized pairwise predictors for this replicate.
  pred <- build_standardized_predictors(pair_obj, n_targets)
  
  # Start from the nuisance beta vector.
  beta_vec <- beta_nuisance
  
  # Calibrate the raw focal coefficient from the target standardized beta.
  beta_focal_raw <- raw_from_target_beta(
    beta_std = beta_focal_std,
    var_x = calibration_obj$var_x,
    var_other = calibration_obj$var_other
  )
  
  # Insert the calibrated focal coefficient into the beta vector.
  beta_vec[focal_trait] <- beta_focal_raw
  
  # Compute the fixed-effects linear predictor.
  eta <- beta_vec["dH"] * pred$dH +
    beta_vec["dE"] * pred$dE +
    beta_vec["dX"] * pred$dX +
    beta_vec["dA"] * pred$dA +
    beta_vec["dC"] * pred$dC +
    beta_vec["dO"] * pred$dO
  
  # Draw random intercepts for the first and second target in each pair.
  b_i <- rnorm(n_targets, mean = 0, sd = sd_subj_i)
  b_j <- rnorm(n_targets, mean = 0, sd = sd_subj_j)
  
  # Generate the raw outcome by combining fixed effects, random effects, and residual noise.
  y_raw <- eta + b_i[pair_obj$subj_i] + b_j[pair_obj$subj_j] +
    rnorm(pair_obj$n_pairs, mean = 0, sd = sd_resid)
  
  # Standardize the outcome within this replicate.
  y <- scale_safe(y_raw)
  
  # Return the simulated dataset.
  data.frame(
    y = y,
    dH = pred$dH,
    dE = pred$dE,
    dX = pred$dX,
    dA = pred$dA,
    dC = pred$dC,
    dO = pred$dO,
    subj_i = pair_obj$subj_i_f,
    subj_j = pair_obj$subj_j_f
  )
}

## --------------------------------------------------
## 5. Fit one model
## --------------------------------------------------

# Fit one mixed-effects model and extract the focal estimate and p-value.
fit_one_model_fast <- function(dat, focal_trait, fit_control) {
  
  # Track whether a convergence-related warning occurred.
  conv_warn <- FALSE
  
  # Fit the model while capturing warnings and converting errors to NULL.
  fit <- tryCatch(
    withCallingHandlers(
      lmer(
        y ~ dH + dE + dX + dA + dC + dO + (1 | subj_i) + (1 | subj_j),
        data = dat,
        REML = FALSE,
        control = fit_control
      ),
      warning = function(w) {
        msg <- conditionMessage(w)
        if (grepl("failed to converge|max\\|grad\\||unable to evaluate scaled gradient|degenerate",
                  msg, ignore.case = TRUE)) {
          conv_warn <<- TRUE
        }
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) NULL
  )
  
  # Return missing values if the model failed.
  if (is.null(fit)) {
    return(c(
      estimate = NA_real_,
      p_value = NA_real_,
      singular = NA_real_,
      conv_warning = NA_real_
    ))
  }
  
  # Extract the coefficient table.
  coefs <- coef(summary(fit))
  
  # Detect the p-value column name.
  p_col <- grep("^Pr\\(", colnames(coefs), value = TRUE)
  
  # Return missing values if the focal coefficient is unavailable.
  if (!(focal_trait %in% rownames(coefs)) || length(p_col) == 0) {
    return(c(
      estimate = NA_real_,
      p_value = NA_real_,
      singular = as.numeric(isSingular(fit, tol = 1e-4)),
      conv_warning = as.numeric(conv_warn)
    ))
  }
  
  # Return the focal estimate, p-value, and fitting diagnostics.
  c(
    estimate = unname(coefs[focal_trait, "Estimate"]),
    p_value = unname(coefs[focal_trait, p_col[1]]),
    singular = as.numeric(isSingular(fit, tol = 1e-4)),
    conv_warning = as.numeric(conv_warn)
  )
}

# Run one full simulation replicate: simulate data, fit the model, extract results.
simulate_one_rep <- function(pair_obj,
                             n_targets,
                             beta_focal_std,
                             focal_trait,
                             beta_nuisance,
                             sd_subj_i,
                             sd_subj_j,
                             sd_resid,
                             fit_control,
                             calibration_obj) {
  
  # Simulate one dataset.
  dat <- simulate_one_dataset_fast(
    pair_obj = pair_obj,
    n_targets = n_targets,
    beta_focal_std = beta_focal_std,
    focal_trait = focal_trait,
    beta_nuisance = beta_nuisance,
    sd_subj_i = sd_subj_i,
    sd_subj_j = sd_subj_j,
    sd_resid = sd_resid,
    calibration_obj = calibration_obj
  )
  
  # Fit the model and return the extracted results.
  fit_one_model_fast(dat, focal_trait = focal_trait, fit_control = fit_control)
}

# Run a batch of simulation replications in parallel.
run_batch_parallel <- function(n_batch,
                               pair_obj,
                               n_targets,
                               beta_focal_std,
                               focal_trait,
                               beta_nuisance,
                               sd_subj_i,
                               sd_subj_j,
                               sd_resid,
                               fit_control,
                               n_cores,
                               calibration_obj) {
  
  res_list <- parallel::mclapply(
    X = seq_len(n_batch),
    FUN = function(i) {
      simulate_one_rep(
        pair_obj = pair_obj,
        n_targets = n_targets,
        beta_focal_std = beta_focal_std,
        focal_trait = focal_trait,
        beta_nuisance = beta_nuisance,
        sd_subj_i = sd_subj_i,
        sd_subj_j = sd_subj_j,
        sd_resid = sd_resid,
        fit_control = fit_control,
        calibration_obj = calibration_obj
      )
    },
    mc.cores = n_cores,
    mc.set.seed = TRUE,
    mc.preschedule = FALSE
  )
  
  do.call(rbind, res_list)
}

## --------------------------------------------------
## 6. Simulate one design condition
## --------------------------------------------------

# Simulate all replications for one design condition with adaptive stopping.
simulate_condition_adaptive <- function(n_targets,
                                        beta_focal,
                                        focal_trait,
                                        beta_nuisance,
                                        sd_subj_i,
                                        sd_subj_j,
                                        sd_resid,
                                        alpha_bonf,
                                        cond_index,
                                        n_conditions_total,
                                        global_start_time,
                                        pair_cache,
                                        calibration_cache,
                                        batch_size,
                                        min_nsim,
                                        max_nsim,
                                        mc_se_target,
                                        fit_control,
                                        n_cores,
                                        show_progress = TRUE) {
  
  # Retrieve the precomputed pair object and calibration constants.
  pair_obj <- pair_cache[[as.character(n_targets)]]
  calibration_obj <- calibration_cache[[as.character(n_targets)]]
  
  # Preallocate the simulation-results matrix.
  sim_mat <- matrix(NA_real_, nrow = max_nsim, ncol = 4)
  colnames(sim_mat) <- c("estimate", "p_value", "singular", "conv_warning")
  
  # Initialize counters and metadata for this condition.
  cond_start_time <- Sys.time()
  stop_reason <- "max_nsim_reached"
  sims_run <- 0L
  usable <- 0L
  n_detected <- 0L
  
  # Keep simulating until a stopping rule is triggered.
  while (sims_run < max_nsim) {
    
    # Determine the size of the next batch.
    n_batch <- min(batch_size, max_nsim - sims_run)
    
    # Compute the row indices where this batch will be stored.
    start_idx <- sims_run + 1L
    end_idx <- sims_run + n_batch
    
    # Run the next batch of parallel simulations.
    batch_mat <- run_batch_parallel(
      n_batch = n_batch,
      pair_obj = pair_obj,
      n_targets = n_targets,
      beta_focal_std = beta_focal,
      focal_trait = focal_trait,
      beta_nuisance = beta_nuisance,
      sd_subj_i = sd_subj_i,
      sd_subj_j = sd_subj_j,
      sd_resid = sd_resid,
      fit_control = fit_control,
      n_cores = n_cores,
      calibration_obj = calibration_obj
    )
    
    # Store the batch results in the preallocated matrix.
    sim_mat[start_idx:end_idx, ] <- batch_mat
    
    # Update the total number of simulations run.
    sims_run <- end_idx
    
    # Extract p-values and estimates from the current batch.
    batch_p <- batch_mat[, "p_value"]
    batch_est <- batch_mat[, "estimate"]
    
    # Count usable fits and successful detections in the batch.
    usable_batch <- sum(!is.na(batch_p))
    detected_batch <- sum(!is.na(batch_p) & batch_p < alpha_bonf & batch_est > 0)
    
    # Update cumulative counts.
    usable <- usable + usable_batch
    n_detected <- n_detected + detected_batch
    
    # Recompute power statistics based on all completed usable runs.
    stats_now <- compute_power_stats(
      n_detected = n_detected,
      usable = usable,
      mc_se_target = mc_se_target,
      exact_ci = FALSE
    )
    
    # Update the progress bar if requested.
    if (show_progress) {
      elapsed_global_sec <- as.numeric(difftime(Sys.time(), global_start_time, units = "secs"))
      elapsed_cond_sec <- as.numeric(difftime(Sys.time(), cond_start_time, units = "secs"))
      avg_sec_per_sim <- elapsed_cond_sec / sims_run
      eta_cond_sec <- avg_sec_per_sim * (max_nsim - sims_run)
      
      cat(progress_line_adaptive(
        current = sims_run,
        max_total = max_nsim,
        elapsed_sec = elapsed_global_sec,
        eta_sec = eta_cond_sec,
        prefix = paste0("Cond ", cond_index, "/", n_conditions_total),
        power = stats_now$power,
        mc_se = stats_now$mc_se_power
      ))
      flush.console()
    }
    
    # Stop early if the minimum usable sample size is reached and the MCSE is small enough.
    if (usable >= min_nsim &&
        !is.na(stats_now$mc_se_power) &&
        stats_now$mc_se_power < mc_se_target) {
      stop_reason <- "mc_se_target_reached"
      break
    }
  }
  
  # Move to a new line after the progress bar.
  if (show_progress) cat("\n")
  
  # Keep only the rows corresponding to simulations that actually ran.
  sim_used <- sim_mat[seq_len(sims_run), , drop = FALSE]
  
  # Compute final usable-fit and detection counts.
  final_usable <- sum(!is.na(sim_used[, "p_value"]))
  final_detected <- sum(!is.na(sim_used[, "p_value"]) &
                          sim_used[, "p_value"] < alpha_bonf &
                          sim_used[, "estimate"] > 0)
  
  # Compute final power statistics, including the exact confidence interval.
  stats <- compute_power_stats(
    n_detected = final_detected,
    usable = final_usable,
    mc_se_target = mc_se_target,
    exact_ci = TRUE
  )
  
  # Extract the vector of focal estimates across usable runs.
  est_vec <- sim_used[, "estimate"]
  n_est_ok <- sum(!is.na(est_vec))
  
  # Compute the mean, standard deviation, and MCSE of the focal estimate.
  mean_est <- safe_mean(est_vec)
  sd_est <- if (n_est_ok >= 2) sd(est_vec, na.rm = TRUE) else NA_real_
  mc_se_est <- if (n_est_ok >= 2) sd_est / sqrt(n_est_ok) else NA_real_
  
  # Return one summary row for this design condition.
  data.frame(
    n_targets = n_targets,
    beta_true = beta_focal,
    sims_run = sims_run,
    usable_fits = final_usable,
    n_detected = stats$n_detected,
    power = stats$power,
    mc_se_power = stats$mc_se_power,
    ci_low_power = stats$ci_low_power,
    ci_high_power = stats$ci_high_power,
    mean_estimate = mean_est,
    bias = mean_est - beta_focal,
    sd_estimate = sd_est,
    mc_se_estimate = mc_se_est,
    singular_rate = safe_mean(sim_used[, "singular"]),
    conv_warn_rate = safe_mean(sim_used[, "conv_warning"]),
    precision_ok = stats$precision_ok,
    stop_reason = stop_reason
  )
}

## --------------------------------------------------
## 7. Start from scratch and remove old checkpoints
## --------------------------------------------------

# Delete the main checkpoint file if it already exists.
if (file.exists(checkpoint_file)) file.remove(checkpoint_file)

# Delete the backup checkpoint file if it already exists.
if (file.exists(backup_file)) file.remove(backup_file)

# Initialize an empty results data frame.
power_results <- data.frame()

# Create the vector of row indices to run.
remaining_idx <- seq_len(nrow(design_grid))

# Print a startup message.
cat("Starting from scratch. Checkpoints will be saved only during this run.\n")

## --------------------------------------------------
## 8. Outer loop with checkpoint saving
## --------------------------------------------------

# Record the global start time.
global_start_time <- Sys.time()

# Count the total number of design conditions.
n_total <- length(remaining_idx)

# Loop over all design conditions.
for (k in seq_along(remaining_idx)) {
  
  # Get the design-grid row index for the current condition.
  i <- remaining_idx[k]
  
  # Extract the current design condition.
  this_row <- design_grid[i, ]
  
  # Print a header for the current condition.
  cat(
    "\n====================================================\n",
    "Condition ", k, " / ", n_total,
    " | n = ", this_row$n_targets,
    " | target standardized beta = ", this_row$beta_true,
    "\n====================================================\n",
    sep = ""
  )
  
  # Run the adaptive simulation for the current condition.
  cond_res <- simulate_condition_adaptive(
    n_targets = this_row$n_targets,
    beta_focal = this_row$beta_true,
    focal_trait = focal_trait,
    beta_nuisance = beta_nuisance,
    sd_subj_i = sd_subj_i,
    sd_subj_j = sd_subj_j,
    sd_resid = sd_resid,
    alpha_bonf = alpha_bonf,
    cond_index = k,
    n_conditions_total = n_total,
    global_start_time = global_start_time,
    pair_cache = pair_cache,
    calibration_cache = calibration_cache,
    batch_size = batch_size,
    min_nsim = min_nsim,
    max_nsim = max_nsim,
    mc_se_target = mc_se_target,
    fit_control = fit_control,
    n_cores = n_cores,
    show_progress = TRUE
  )
  
  # Append the condition-level result to the cumulative results table.
  power_results <- rbind(power_results, cond_res)
  
  # Keep results ordered by sample size and target beta.
  power_results <- power_results[order(power_results$n_targets, power_results$beta_true), , drop = FALSE]
  
  # Save both the main and backup checkpoints.
  saveRDS(power_results, checkpoint_file)
  saveRDS(power_results, backup_file)
  
  # Compute total elapsed time in minutes.
  elapsed_total_min <- as.numeric(difftime(Sys.time(), global_start_time, units = "mins"))
  
  # Print a checkpoint confirmation message.
  cat(sprintf(
    "Saved checkpoint after condition %d/%d | total elapsed: %.1f min\n",
    k, n_total, elapsed_total_min
  ))
  
  # Flush console output so progress is updated immediately.
  flush.console()
}

## --------------------------------------------------
## 9. Final summary
## --------------------------------------------------

# Build a minimum detectable effect table for each sample size.
mde_table <- do.call(
  rbind,
  lapply(sort(unique(power_results$n_targets)), function(nn) {
    sub <- power_results[power_results$n_targets == nn, , drop = FALSE]
    ok <- !is.na(sub$power) & sub$power >= target_power
    
    data.frame(
      n_targets = nn,
      min_detectable_beta = if (any(ok)) min(sub$beta_true[ok], na.rm = TRUE) else NA_real_
    )
  })
)

# Copy the results table for formatted printing.
power_results_print <- power_results

# Columns to round to 3 decimal places.
num_cols_3 <- c("power", "mc_se_power", "ci_low_power", "ci_high_power",
                "singular_rate", "conv_warn_rate")

# Columns to round to 6 decimal places.
num_cols_6 <- c("mean_estimate", "bias", "sd_estimate", "mc_se_estimate")

# Round selected columns for easier reading.
power_results_print[num_cols_3] <- lapply(power_results_print[num_cols_3], round, 3)
power_results_print[num_cols_6] <- lapply(power_results_print[num_cols_6], round, 6)

# Print the minimum detectable effect table.
print(mde_table)

# Print the formatted results table.
print(power_results_print)

## --------------------------------------------------
## 10. Save final outputs
## --------------------------------------------------

# Uncomment to save the full condition-level results.
# saveRDS(power_results, "power_results_final.rds")

# Uncomment to save the minimum detectable effect table.
# saveRDS(mde_table, "mde_table_final.rds")