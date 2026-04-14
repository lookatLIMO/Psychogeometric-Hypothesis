# ============================================================
# MONTE CARLO POWER SIMULATION FOR A LLM IN THE ANCHOR DESIGN
# ------------------------------------------------------------
#
# Purpose:
# This script estimates statistical power for detecting a focal
# pairwise-difference predictor (e.g., dH) in a linear mixed-effects
# model with pairwise differences on a singular MDS dimension as outcome.
# LMM is fit with random intercepts for the two members of each pair.
#
# General idea:
# For each design condition, the script repeatedly simulates synthetic
# pairwise data, fits the same mixed model, extracts the focal estimate
# and p-value, and computes the proportion of simulations in which the
# focal effect is detected in the expected positive direction.
#
# Model structure:
# The simulated outcome y depends on six standardized pairwise-difference
# predictors (dH, dE, dX, dA, dC, dO), plus random intercepts for subj_i
# and subj_j, and residual error.
#
# Simulation design:
# The script evaluates one or more combinations of:
# - number of targets (n_targets)
# - true focal effect size (beta_true)
#
# Statistical decision rule:
# A result is counted as detected when the focal predictor has:
# - a non-missing model fit,
# - a p-value below the Bonferroni-corrected alpha threshold,
# - and an estimated coefficient in the expected positive direction.
#
# Adaptive stopping:
# For each condition, simulations continue until either:
# - the Monte Carlo standard error (MCSE) of the estimated power drops
#   below the target threshold after a minimum number of usable runs, or
# - the maximum number of simulations is reached.
#
# Outputs:
# The script returns:
# - condition-level power estimates,
# - Monte Carlo uncertainty measures,
# - confidence intervals for power,
# - estimate bias and variability,
# - convergence and singularity rates,
# - and a final minimum detectable effect (MDE) summary table.
#
# Notes:
# - Checkpoints are saved after each condition.
# - Pair structures are precomputed for efficiency.
# - Parallel execution can be enabled by increasing n_cores.
# ============================================================

## --------------------------------------------------
## 0. Packages
## --------------------------------------------------

# Install lme4 if it is not already available.
if (!requireNamespace("lme4", quietly = TRUE)) install.packages("lme4")

# Install lmerTest if it is not already available.
if (!requireNamespace("lmerTest", quietly = TRUE)) install.packages("lmerTest")

# Load lme4 for mixed-effects model estimation.
library(lme4)

# Load lmerTest so that p-values are available in model summaries.
library(lmerTest)

# Set a global random seed for reproducibility.
set.seed(123)

## --------------------------------------------------
## 1. Global settings
## --------------------------------------------------

# Number of simulations to run in each batch.
batch_size <- 1

# Minimum number of usable simulations required before early stopping is allowed.
min_nsim <- 100

# Maximum number of simulations allowed for each condition.
max_nsim <- 1000

# Monte Carlo standard error target for stopping.
mc_se_target <- 0.01

# Grid of sample sizes to evaluate.
n_grid <- c(1505)

# Grid of true focal beta values to evaluate.
beta_grid <- seq(0.003, 0.004, by = 0.001)

# Desired target power.
target_power <- 0.80

# Familywise alpha level.
alpha_family <- 0.05

# Number of tests used for Bonferroni correction.
n_tests <- 18 #assuming 6 predictors (traits) per dimension and 3 MDS dimensions

# Bonferroni-corrected alpha threshold.
alpha_bonf <- alpha_family / n_tests

# Name of the focal predictor whose power is being evaluated.
focal_trait <- "dH"

# Set all nuisance effects to zero by default.
beta_nuisance <- c(
  dH = 0.00,
  dE = 0.00,
  dX = 0.00,
  dA = 0.00,
  dC = 0.00,
  dO = 0.00
)

# Standard deviation of the random intercept for subject i.
sd_subj_i <- 0.40

# Standard deviation of the random intercept for subject j.
sd_subj_j <- 0.40

# Residual standard deviation.
sd_resid <- 1.00

# File used to save the main checkpoint during the run.
checkpoint_file <- "power_results_checkpoint.rds"

# File used to save a backup checkpoint during the run.
backup_file <- "power_results_backup.rds"

## --- Parallel settings for Mac ---

# Detect the number of available physical CPU cores.
available_cores <- parallel::detectCores(logical = FALSE)

# Fall back to 2 cores if detection fails.
if (is.na(available_cores)) available_cores <- 2L

# Number of cores to actually use.
n_cores <- 7L

## --- Fast lmer control ---

# Define a faster lmer control object to reduce fitting overhead.
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
  # Generate all 2-combinations of target indices.
  cmb <- utils::combn(n_targets, 2)
  
  # Return the two elements of each pair as separate vectors.
  list(
    subj_i = cmb[1, ],
    subj_j = cmb[2, ]
  )
}

# Build and store all pair-related objects for one sample size.
build_pair_object <- function(n_targets) {
  # Generate the pair indices.
  p <- make_pairs(n_targets)
  
  # Convert the first element of each pair to a factor for random effects.
  p$subj_i_f <- factor(p$subj_i, levels = seq_len(n_targets))
  
  # Convert the second element of each pair to a factor for random effects.
  p$subj_j_f <- factor(p$subj_j, levels = seq_len(n_targets))
  
  # Store the total number of pairs.
  p$n_pairs <- length(p$subj_i)
  
  # Return the completed pair object.
  p
}

# Safely standardize a numeric vector.
scale_safe <- function(x) {
  # Compute the standard deviation.
  s <- sd(x)
  
  # If the standard deviation is missing or zero, return zeros.
  if (is.na(s) || s == 0) {
    rep(0, length(x))
  } else {
    # Otherwise return the z-scored vector.
    (x - mean(x)) / s
  }
}

# Safely compute the mean of a vector.
safe_mean <- function(x) {
  # Return NA if all values are missing; otherwise return the mean.
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

# Compute power-related summary statistics.
compute_power_stats <- function(n_detected,
                                usable,
                                mc_se_target = 0.01,
                                exact_ci = FALSE) {
  
  # If there are no usable fits, return missing statistics.
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
  
  # Estimate power as the proportion of detected effects.
  power_hat <- n_detected / usable
  
  # Compute Monte Carlo standard error when at least two usable fits exist.
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
  
  # Check whether the Monte Carlo precision target has been reached.
  precision_ok <- !is.na(mc_se) && mc_se < mc_se_target
  
  # Return all power summary statistics.
  list(
    n_detected = n_detected,
    power = power_hat,
    mc_se_power = mc_se,
    ci_low_power = ci_low,
    ci_high_power = ci_high,
    precision_ok = precision_ok
  )
}

# Build a progress bar string with adaptive simulation information.
progress_line_adaptive <- function(current,
                                   max_total,
                                   elapsed_sec,
                                   eta_sec,
                                   prefix = "",
                                   power = NA_real_,
                                   mc_se = NA_real_,
                                   width = 28) {
  
  # Compute progress as a fraction of the maximum.
  frac <- if (max_total > 0) current / max_total else 0
  
  # Clamp the fraction to the [0, 1] interval.
  frac <- min(max(frac, 0), 1)
  
  # Compute the number of completed bar characters.
  done <- floor(width * frac)
  
  # Compute the number of remaining bar characters.
  left <- width - done
  
  # Build the text progress bar.
  bar <- paste0(
    "[",
    paste(rep("=", done), collapse = ""),
    paste(rep("-", left), collapse = ""),
    "]"
  )
  
  # Return a formatted one-line progress string.
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

## --------------------------------------------------
## 3. Build design grid and precompute pair structures
## --------------------------------------------------

# Build the full design grid of conditions to simulate.
design_grid <- expand.grid(
  n_targets = n_grid,
  beta_true = beta_grid,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

# Extract and sort the unique n_targets values.
unique_n <- sort(unique(design_grid$n_targets))

# Precompute pair objects for each unique n_targets value.
pair_cache <- setNames(
  lapply(unique_n, build_pair_object),
  as.character(unique_n)
)

# Print a message confirming that precomputation is complete.
cat("Precomputation completed for n_targets values:",
    paste(unique_n, collapse = ", "), "\n")

# Print the number of parallel cores in use.
cat("Parallel cores used:", n_cores, "\n")

## --------------------------------------------------
## 4. Simulate one dataset (FAST)
## --------------------------------------------------

# Simulate one dataset for a single replication.
simulate_one_dataset_fast <- function(pair_obj,
                                      n_targets,
                                      beta_focal,
                                      focal_trait,
                                      beta_nuisance,
                                      sd_subj_i,
                                      sd_subj_j,
                                      sd_resid) {
  
  # Draw ordinal trait values from 1 to 7 for H.
  H <- sample.int(7, n_targets, replace = TRUE)
  
  # Draw ordinal trait values from 1 to 7 for E.
  E <- sample.int(7, n_targets, replace = TRUE)
  
  # Draw ordinal trait values from 1 to 7 for X.
  X <- sample.int(7, n_targets, replace = TRUE)
  
  # Draw ordinal trait values from 1 to 7 for A.
  A <- sample.int(7, n_targets, replace = TRUE)
  
  # Draw ordinal trait values from 1 to 7 for C.
  C <- sample.int(7, n_targets, replace = TRUE)
  
  # Draw ordinal trait values from 1 to 7 for O.
  O <- sample.int(7, n_targets, replace = TRUE)
  
  # Compute standardized absolute pairwise differences for H.
  dH <- scale_safe(abs(H[pair_obj$subj_i] - H[pair_obj$subj_j]))
  
  # Compute standardized absolute pairwise differences for E.
  dE <- scale_safe(abs(E[pair_obj$subj_i] - E[pair_obj$subj_j]))
  
  # Compute standardized absolute pairwise differences for X.
  dX <- scale_safe(abs(X[pair_obj$subj_i] - X[pair_obj$subj_j]))
  
  # Compute standardized absolute pairwise differences for A.
  dA <- scale_safe(abs(A[pair_obj$subj_i] - A[pair_obj$subj_j]))
  
  # Compute standardized absolute pairwise differences for C.
  dC <- scale_safe(abs(C[pair_obj$subj_i] - C[pair_obj$subj_j]))
  
  # Compute standardized absolute pairwise differences for O.
  dO <- scale_safe(abs(O[pair_obj$subj_i] - O[pair_obj$subj_j]))
  
  # Start from the nuisance beta vector.
  beta_vec <- beta_nuisance
  
  # Replace the focal trait coefficient with the focal beta value.
  beta_vec[focal_trait] <- beta_focal
  
  # Compute the fixed-effects linear predictor.
  eta <- beta_vec["dH"] * dH +
    beta_vec["dE"] * dE +
    beta_vec["dX"] * dX +
    beta_vec["dA"] * dA +
    beta_vec["dC"] * dC +
    beta_vec["dO"] * dO
  
  # Draw random intercepts for the first member of each pair.
  b_i <- rnorm(n_targets, mean = 0, sd = sd_subj_i)
  
  # Draw random intercepts for the second member of each pair.
  b_j <- rnorm(n_targets, mean = 0, sd = sd_subj_j)
  
  # Add random intercept contributions to the linear predictor.
  eta <- eta + b_i[pair_obj$subj_i] + b_j[pair_obj$subj_j]
  
  # Add residual noise to generate the outcome.
  y <- eta + rnorm(pair_obj$n_pairs, mean = 0, sd = sd_resid)
  
  # Return the simulated dataset.
  data.frame(
    y = y,
    dH = dH,
    dE = dE,
    dX = dX,
    dA = dA,
    dC = dC,
    dO = dO,
    subj_i = pair_obj$subj_i_f,
    subj_j = pair_obj$subj_j_f
  )
}

## --------------------------------------------------
## 5. Fit one model
## --------------------------------------------------

# Fit one mixed-effects model and extract the focal estimate and p-value.
fit_one_model_fast <- function(dat, focal_trait, fit_control) {
  
  # Track whether a convergence-related warning occurs.
  conv_warn <- FALSE
  
  # Fit the model while catching warnings and errors.
  fit <- tryCatch(
    withCallingHandlers(
      lmer(
        y ~ dH + dE + dX + dA + dC + dO + (1 | subj_i) + (1 | subj_j),
        data = dat,
        REML = FALSE,
        control = fit_control
      ),
      warning = function(w) {
        # Extract the warning message text.
        msg <- conditionMessage(w)
        
        # Flag convergence-related warnings.
        if (grepl("failed to converge|max\\|grad\\||unable to evaluate scaled gradient|degenerate",
                  msg, ignore.case = TRUE)) {
          conv_warn <<- TRUE
        }
        
        # Suppress the warning from printing to the console.
        invokeRestart("muffleWarning")
      }
    ),
    # Return NULL if model fitting fails.
    error = function(e) NULL
  )
  
  # Return missing values if the model failed completely.
  if (is.null(fit)) {
    return(c(
      estimate = NA_real_,
      p_value = NA_real_,
      singular = NA_real_,
      conv_warning = NA_real_
    ))
  }
  
  # Extract the coefficient table from the fitted model summary.
  coefs <- coef(summary(fit))
  
  # Identify the p-value column name.
  p_col <- grep("^Pr\\(", colnames(coefs), value = TRUE)
  
  # Return missing values if the focal coefficient or p-value column is unavailable.
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
                             beta_focal,
                             focal_trait,
                             beta_nuisance,
                             sd_subj_i,
                             sd_subj_j,
                             sd_resid,
                             fit_control) {
  
  # Simulate one dataset.
  dat <- simulate_one_dataset_fast(
    pair_obj = pair_obj,
    n_targets = n_targets,
    beta_focal = beta_focal,
    focal_trait = focal_trait,
    beta_nuisance = beta_nuisance,
    sd_subj_i = sd_subj_i,
    sd_subj_j = sd_subj_j,
    sd_resid = sd_resid
  )
  
  # Fit the model and return the extracted results.
  fit_one_model_fast(dat, focal_trait = focal_trait, fit_control = fit_control)
}

# Run a batch of simulation replicates in parallel.
run_batch_parallel <- function(n_batch,
                               pair_obj,
                               n_targets,
                               beta_focal,
                               focal_trait,
                               beta_nuisance,
                               sd_subj_i,
                               sd_subj_j,
                               sd_resid,
                               fit_control,
                               n_cores) {
  
  # Run one simulation replicate for each batch index.
  res_list <- parallel::mclapply(
    X = seq_len(n_batch),
    FUN = function(i) {
      simulate_one_rep(
        pair_obj = pair_obj,
        n_targets = n_targets,
        beta_focal = beta_focal,
        focal_trait = focal_trait,
        beta_nuisance = beta_nuisance,
        sd_subj_i = sd_subj_i,
        sd_subj_j = sd_subj_j,
        sd_resid = sd_resid,
        fit_control = fit_control
      )
    },
    mc.cores = n_cores,
    mc.set.seed = TRUE,
    mc.preschedule = FALSE
  )
  
  # Combine the batch results into a matrix.
  do.call(rbind, res_list)
}

## --------------------------------------------------
## 6. Simulate one design condition (ADAPTIVE + ETA + PARALLEL)
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
                                        batch_size,
                                        min_nsim,
                                        max_nsim,
                                        mc_se_target,
                                        fit_control,
                                        n_cores,
                                        show_progress = TRUE) {
  
  # Retrieve the precomputed pair object for this sample size.
  pair_obj <- pair_cache[[as.character(n_targets)]]
  
  # Preallocate a matrix to store simulation results.
  sim_mat <- matrix(NA_real_, nrow = max_nsim, ncol = 4)
  
  # Name the columns of the result matrix.
  colnames(sim_mat) <- c("estimate", "p_value", "singular", "conv_warning")
  
  # Store the starting time of this condition.
  cond_start_time <- Sys.time()
  
  # Default stopping reason if the loop reaches max_nsim.
  stop_reason <- "max_nsim_reached"
  
  # Initialize the number of simulations run.
  sims_run <- 0L
  
  # Initialize the number of usable model fits.
  usable <- 0L
  
  # Initialize the number of detected significant effects.
  n_detected <- 0L
  
  # Continue until the maximum number of simulations is reached or stopping occurs.
  while (sims_run < max_nsim) {
    
    # Determine the size of the next batch.
    n_batch <- min(batch_size, max_nsim - sims_run)
    
    # Compute the starting row index for this batch.
    start_idx <- sims_run + 1L
    
    # Compute the ending row index for this batch.
    end_idx <- sims_run + n_batch
    
    # Run the next batch of simulations.
    batch_mat <- run_batch_parallel(
      n_batch = n_batch,
      pair_obj = pair_obj,
      n_targets = n_targets,
      beta_focal = beta_focal,
      focal_trait = focal_trait,
      beta_nuisance = beta_nuisance,
      sd_subj_i = sd_subj_i,
      sd_subj_j = sd_subj_j,
      sd_resid = sd_resid,
      fit_control = fit_control,
      n_cores = n_cores
    )
    
    # Store batch results in the preallocated result matrix.
    sim_mat[start_idx:end_idx, ] <- batch_mat
    
    # Update the total number of simulations run.
    sims_run <- end_idx
    
    # Extract p-values from the current batch.
    batch_p <- batch_mat[, "p_value"]
    
    # Extract estimates from the current batch.
    batch_est <- batch_mat[, "estimate"]
    
    # Count the number of usable fits in this batch.
    usable_batch <- sum(!is.na(batch_p))
    
    # Count the number of detected effects in this batch.
    detected_batch <- sum(!is.na(batch_p) & batch_p < alpha_bonf & batch_est > 0)
    
    # Update the cumulative number of usable fits.
    usable <- usable + usable_batch
    
    # Update the cumulative number of detected effects.
    n_detected <- n_detected + detected_batch
    
    # Recompute current power statistics.
    stats_now <- compute_power_stats(
      n_detected = n_detected,
      usable = usable,
      mc_se_target = mc_se_target,
      exact_ci = FALSE
    )
    
    # Print progress information if requested.
    if (show_progress) {
      # Compute total elapsed time since the entire run started.
      elapsed_global_sec <- as.numeric(difftime(Sys.time(), global_start_time, units = "secs"))
      
      # Compute elapsed time for the current condition.
      elapsed_cond_sec <- as.numeric(difftime(Sys.time(), cond_start_time, units = "secs"))
      
      # Compute the average time per simulation so far within this condition.
      avg_sec_per_sim <- elapsed_cond_sec / sims_run
      
      # Estimate remaining time for this condition.
      eta_cond_sec <- avg_sec_per_sim * (max_nsim - sims_run)
      
      # Print the progress line.
      cat(progress_line_adaptive(
        current = sims_run,
        max_total = max_nsim,
        elapsed_sec = elapsed_global_sec,
        eta_sec = eta_cond_sec,
        prefix = paste0("Cond ", cond_index, "/", n_conditions_total),
        power = stats_now$power,
        mc_se = stats_now$mc_se_power
      ))
      
      # Flush the console so the progress bar updates immediately.
      flush.console()
    }
    
    # Stop early if the minimum number of usable fits is reached and MCSE is small enough.
    if (usable >= min_nsim &&
        !is.na(stats_now$mc_se_power) &&
        stats_now$mc_se_power < mc_se_target) {
      stop_reason <- "mc_se_target_reached"
      break
    }
  }
  
  # Move to a new line after the progress bar.
  if (show_progress) cat("\n")
  
  # Keep only the rows that correspond to simulations that were actually run.
  sim_used <- sim_mat[seq_len(sims_run), , drop = FALSE]
  
  # Count the final number of usable fits.
  final_usable <- sum(!is.na(sim_used[, "p_value"]))
  
  # Count the final number of detected effects.
  final_detected <- sum(!is.na(sim_used[, "p_value"]) &
                          sim_used[, "p_value"] < alpha_bonf &
                          sim_used[, "estimate"] > 0)
  
  # Compute final power statistics, now including an exact confidence interval.
  stats <- compute_power_stats(
    n_detected = final_detected,
    usable = final_usable,
    mc_se_target = mc_se_target,
    exact_ci = TRUE
  )
  
  # Extract the vector of focal estimates.
  est_vec <- sim_used[, "estimate"]
  
  # Count non-missing estimates.
  n_est_ok <- sum(!is.na(est_vec))
  
  # Compute the mean estimate.
  mean_est <- safe_mean(est_vec)
  
  # Compute the standard deviation of estimates if possible.
  sd_est <- if (n_est_ok >= 2) sd(est_vec, na.rm = TRUE) else NA_real_
  
  # Compute the Monte Carlo standard error of the estimate mean if possible.
  mc_se_est <- if (n_est_ok >= 2) sd_est / sqrt(n_est_ok) else NA_real_
  
  # Return a one-row data frame summarizing this design condition.
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
## 7. Start from zero and remove old checkpoints
## --------------------------------------------------

# Delete the main checkpoint file if it already exists.
if (file.exists(checkpoint_file)) file.remove(checkpoint_file)

# Delete the backup checkpoint file if it already exists.
if (file.exists(backup_file)) file.remove(backup_file)

# Initialize an empty results data frame.
power_results <- data.frame()

# Create the index vector for all conditions to run.
remaining_idx <- seq_len(nrow(design_grid))

# Print a message indicating that the run starts from scratch.
cat("Starting from scratch. Checkpoints will be saved only during this run.\n")

## --------------------------------------------------
## 8. Outer loop with checkpoint saving
## --------------------------------------------------

# Store the global starting time.
global_start_time <- Sys.time()

# Count the total number of design conditions.
n_total <- length(remaining_idx)

# Loop over all design conditions.
for (k in seq_along(remaining_idx)) {
  
  # Get the row index for the current condition.
  i <- remaining_idx[k]
  
  # Extract the current condition from the design grid.
  this_row <- design_grid[i, ]
  
  # Print a header for the current condition.
  cat(
    "\n====================================================\n",
    "Condition ", k, " / ", n_total,
    " | n = ", this_row$n_targets,
    " | beta = ", this_row$beta_true,
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
    batch_size = batch_size,
    min_nsim = min_nsim,
    max_nsim = max_nsim,
    mc_se_target = mc_se_target,
    fit_control = fit_control,
    n_cores = n_cores,
    show_progress = TRUE
  )
  
  # Append the current condition result to the cumulative results.
  power_results <- rbind(power_results, cond_res)
  
  # Sort results by sample size and true beta.
  power_results <- power_results[order(power_results$n_targets, power_results$beta_true), , drop = FALSE]
  
  # Save the main checkpoint file.
  saveRDS(power_results, checkpoint_file)
  
  # Save the backup checkpoint file.
  saveRDS(power_results, backup_file)
  
  # Compute total elapsed time in minutes.
  elapsed_total_min <- as.numeric(difftime(Sys.time(), global_start_time, units = "mins"))
  
  # Print a checkpoint confirmation message.
  cat(sprintf(
    "Saved checkpoint after condition %d/%d | total elapsed: %.1f min\n",
    k, n_total, elapsed_total_min
  ))
  
  # Flush the console output.
  flush.console()
}

## --------------------------------------------------
## 9. Final summary
## --------------------------------------------------

# Build a minimum detectable effect table for each sample size.
mde_table <- do.call(
  rbind,
  lapply(sort(unique(power_results$n_targets)), function(nn) {
    # Subset results for the current sample size.
    sub <- power_results[power_results$n_targets == nn, , drop = FALSE]
    
    # Identify conditions that reached the target power.
    ok <- !is.na(sub$power) & sub$power >= target_power
    
    # Return the minimum detectable beta for this sample size.
    data.frame(
      n_targets = nn,
      min_detectable_beta = if (any(ok)) min(sub$beta_true[ok], na.rm = TRUE) else NA_real_
    )
  })
)

# Copy the results table for prettier printing.
power_results_print <- power_results

# Columns to round to 3 decimal places.
num_cols_3 <- c("power", "mc_se_power", "ci_low_power", "ci_high_power",
                "singular_rate", "conv_warn_rate")

# Columns to round to 6 decimal places.
num_cols_6 <- c("mean_estimate", "bias", "sd_estimate", "mc_se_estimate")

# Round the selected columns to 3 decimals.
power_results_print[num_cols_3] <- lapply(power_results_print[num_cols_3], round, 3)

# Round the selected columns to 6 decimals.
power_results_print[num_cols_6] <- lapply(power_results_print[num_cols_6], round, 6)

# Print the minimum detectable effect table.
print(mde_table)

# Print the rounded results table.
print(power_results_print)

## --------------------------------------------------
## 10. Save final outputs
## --------------------------------------------------

# Uncomment to save the final detailed power results.
# saveRDS(power_results, "power_results_final.rds")

# Uncomment to save the final minimum detectable effect table.
# saveRDS(mde_table, "mde_table_final.rds")