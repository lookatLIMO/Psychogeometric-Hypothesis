# ============================================================
# MONTE CARLO SIMULATION FOR GPA STABILITY WITH ANCHOR-BASED ALIGNMENT
# ------------------------------------------------------------
#
# Purpose:
# This script evaluates how well a common latent geometric structure can
# be recovered across multiple judges when each judge provides noisy
# dissimilarities, MDS is fit at a fixed dimensionality, and the resulting
# configurations are aligned using anchor-based generalized Procrustes
# analysis (GPA).
#
# General idea:
# For each simulation, the script generates a latent configuration in a
# chosen dimensionality (fit_dim), made of:
# - a set of common anchors shared across judges,
# - and a set of judge-specific unique targets.
#
# For each judge:
# 1. A latent configuration is generated.
# 2. A noisy dissimilarity matrix is constructed from that configuration.
# 3. An ordinal MDS solution is fit directly in fit_dim dimensions.
#
# After MDS estimation:
# The recovered judge-specific configurations are aligned using only the
# anchor points through a partial GPA procedure. The aligned solutions are
# then evaluated against the true latent structure.
#
# What is being assessed:
# The simulation measures whether the anchor-based alignment produces a
# stable and accurate reconstruction of the common latent space under
# different design conditions.
#
# Design conditions:
# The script evaluates combinations of:
# - number of judges (n_judges),
# - number of unique targets per judge (n_unique),
# - fitted dimensionality (fit_dim),
# - number of anchors (n_anchors).
#
# Only admissible conditions are simulated:
# A condition is included only if the number of anchors is at least as
# large as the fitted dimensionality (n_anchors >= fit_dim).
#
# Single-simulation evaluation criteria:
# Each simulated dataset is evaluated using:
# - correlation between true and recovered global pairwise distances,
# - RMSE between true and recovered global coordinates,
# - RMSE for the recovered anchors,
# - spread of aligned anchor locations across judges,
# - MDS stress and GPA convergence diagnostics.
#
# Stability rule:
# A single simulation is labeled as stable when all selected recovery
# criteria exceed the predefined thresholds for geometric accuracy and
# anchor consistency.
#
# Condition-level summary:
# After repeating the simulation many times for a condition, the script
# computes summary statistics such as mean recovery quality, proportion
# of stable simulations, GPA convergence rate, and average number of
# usable judges.
#
# Final classification:
# Each condition is classified as:
# - "feasible" if recovery is good and stability is high,
# - "borderline" if recovery is mixed or only moderately stable,
# - "unstable" otherwise.
#
# Outputs:
# The script produces:
# - raw simulation-level results for every condition,
# - a condition-level summary table,
# - and checkpoint files saved during the run.
#
# Notes:
# - The latent structure is generated directly in fit_dim.
# - The alignment uses anchors only, not all targets.
# - The final decision is based on empirical recovery and stability
#   rather than on inferential significance testing.
# ============================================================

# Install pacman if it is not already available.
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

# Load required packages.
pacman::p_load(smacof, dplyr, tibble, readr)

# Set a global random seed for reproducibility.
set.seed(123)

# ------------------------------------------------------------
# 0. Global settings
# ------------------------------------------------------------

# Number of simulations to run for each design condition.
nsim <- 1000

# Grid of numbers of judges to evaluate.
n_judges_grid <- c(150, 200)

# Grid of numbers of unique targets per judge to evaluate.
n_unique_grid <- c(15)

# Grid of possible numbers of anchors.
anchor_grid <- 1:6

# Grid of fitted MDS dimensionalities to evaluate.
fit_dim_grid <- 1:6

# Standard deviation of noise added to dissimilarities.
noise_sd <- 0.20

# Maximum number of iterations allowed in MDS fitting.
itmax_mds <- 1000

# Stability criteria at the single-simulation level.
thr_dist_cor <- 0.90
thr_anchor_rmse <- 0.25
thr_anchor_spread <- 0.20

# Empirical decision thresholds at the condition level.
thr_prop_stable_good <- 0.80
thr_prop_stable_mid  <- 0.50

# File name for the main checkpoint.
checkpoint_file <- "alignment_fixed_dim_checkpoint.rds"

# File name for the backup checkpoint.
backup_file <- "alignment_fixed_dim_backup.rds"

# ------------------------------------------------------------
# 1. Helpers
# ------------------------------------------------------------

# Convert all pairwise distances in a configuration into a vector.
vec_dist <- function(x) {
  as.vector(stats::dist(x))
}

# Center and normalize a configuration matrix.
normalize_config <- function(x) {
  # Convert input to matrix form.
  x <- as.matrix(x)
  
  # Center each column by subtracting its mean.
  xc <- sweep(x, 2, colMeans(x), "-")
  
  # Compute the total sum-of-squares norm.
  ss <- sqrt(sum(xc^2))
  
  # If the norm is essentially zero, return the centered matrix unchanged.
  if (ss < 1e-12) return(xc)
  
  # Otherwise scale the centered matrix to unit norm.
  xc / ss
}

# Build a noisy dissimilarity matrix from a coordinate matrix.
make_noisy_delta <- function(coords, noise_sd) {
  # Compute Euclidean distances between rows.
  d <- as.matrix(stats::dist(coords))
  
  # Generate a symmetric noise matrix.
  e <- matrix(rnorm(length(d), 0, noise_sd), nrow = nrow(d))
  e <- (e + t(e)) / 2
  diag(e) <- 0
  
  # Add noise to the distance matrix.
  delta <- d + e
  
  # Truncate negative dissimilarities to zero.
  delta[delta < 0] <- 0
  
  # Force zeros on the diagonal.
  diag(delta) <- 0
  
  # Return the noisy dissimilarity matrix.
  delta
}

# Safely fit an ordinal MDS solution.
fit_mds_safe <- function(diss_matrix, ndim, itmax = 1000) {
  tryCatch(
    smacof::mds(
      delta = diss_matrix,
      ndim = ndim,
      type = "ordinal",
      ties = "secondary",
      init = "torgerson",
      itmax = itmax
    ),
    error = function(e) NULL
  )
}

# ------------------------------------------------------------
# 2. Procrustes / partial GPA helpers
# ------------------------------------------------------------

# Fit a Procrustes transformation from source anchors to target anchors,
# and apply it to the full source configuration.
procrustes_fit <- function(source_anchor, target_anchor, source_full, allow_scaling = TRUE) {
  # Convert all inputs to matrices.
  X <- as.matrix(source_anchor)
  Y <- as.matrix(target_anchor)
  Z <- as.matrix(source_full)
  
  # Store the number of dimensions.
  p <- ncol(X)
  
  # Special case: with one anchor, only translation is identifiable.
  if (nrow(X) == 1) {
    # Compute the translation needed to align the single source anchor to the target anchor.
    shift <- as.numeric(Y[1, ] - X[1, ])
    
    # Apply the translation to the full configuration.
    Z_aligned <- sweep(Z, 2, shift, "+")
    
    # Apply the same translation to the anchor configuration.
    X_aligned <- sweep(X, 2, shift, "+")
    
    # Return the translation-only alignment object.
    return(list(
      coords = Z_aligned,
      anchor = X_aligned,
      mode = "translation",
      shift = shift,
      scale = 1,
      R = diag(p),
      center_source = rep(0, p),
      center_target = rep(0, p)
    ))
  }
  
  # Compute the centroid of the source anchors.
  mx <- colMeans(X)
  
  # Compute the centroid of the target anchors.
  my <- colMeans(Y)
  
  # Center the source anchors.
  Xc <- sweep(X, 2, mx, "-")
  
  # Center the target anchors.
  Yc <- sweep(Y, 2, my, "-")
  
  # Compute the singular value decomposition of the cross-product matrix.
  sv <- svd(t(Xc) %*% Yc)
  
  # Compute the optimal rotation matrix.
  R <- sv$u %*% t(sv$v)
  
  # Compute the optimal scaling factor, unless scaling is disabled.
  b <- if (allow_scaling) sum(sv$d) / sum(Xc^2) else 1
  
  # Apply centering, rotation, and scaling to the full source configuration.
  Z_aligned <- b * (sweep(Z, 2, mx, "-") %*% R)
  
  # Translate the aligned configuration to the target centroid.
  Z_aligned <- sweep(Z_aligned, 2, my, "+")
  
  # Apply the same transformation to the source anchors.
  X_aligned <- b * (sweep(X, 2, mx, "-") %*% R)
  X_aligned <- sweep(X_aligned, 2, my, "+")
  
  # Return the full Procrustes transformation object.
  list(
    coords = Z_aligned,
    anchor = X_aligned,
    mode = "full",
    shift = rep(0, p),
    scale = b,
    R = R,
    center_source = mx,
    center_target = my
  )
}

# Apply a previously computed Procrustes transformation to a full configuration.
apply_procrustes <- function(source_full, tr) {
  # Convert the input to matrix form.
  Z <- as.matrix(source_full)
  
  # If the transformation is translation-only, just shift the configuration.
  if (tr$mode == "translation") {
    return(sweep(Z, 2, tr$shift, "+"))
  }
  
  # Otherwise apply centering, rotation, scaling, and translation.
  Z_aligned <- tr$scale * (sweep(Z, 2, tr$center_source, "-") %*% tr$R)
  sweep(Z_aligned, 2, tr$center_target, "+")
}

# Perform a partial generalized Procrustes analysis using anchors only.
partial_gpa <- function(config_list, anchor_idx, max_iter = 100, tol = 1e-7) {
  # Number of judges/configurations.
  J <- length(config_list)
  
  # Number of anchors.
  n_anchors <- length(anchor_idx)
  
  # Special case: with one anchor, only translation can be aligned.
  if (n_anchors == 1) {
    # Use the first configuration's anchor as the reference.
    ref_anchor <- config_list[[1]][anchor_idx, , drop = FALSE]
    
    # Preallocate lists for aligned full configurations and aligned anchors.
    aligned_configs <- vector("list", J)
    aligned_anchors <- vector("list", J)
    
    # Align every configuration to the reference anchor using translation only.
    for (j in seq_len(J)) {
      tr <- procrustes_fit(
        source_anchor = config_list[[j]][anchor_idx, , drop = FALSE],
        target_anchor = ref_anchor,
        source_full = config_list[[j]],
        allow_scaling = FALSE
      )
      aligned_configs[[j]] <- tr$coords
      aligned_anchors[[j]] <- tr$anchor
    }
    
    # Compute the average aligned anchor position.
    consensus <- Reduce(`+`, aligned_anchors) / J
    
    # Return the translation-based GPA result.
    return(list(
      aligned_configs = aligned_configs,
      aligned_anchors = aligned_anchors,
      consensus = consensus,
      iterations = 1,
      converged = TRUE
    ))
  }
  
  # Initialize the reference as the normalized anchors from the first configuration.
  ref <- normalize_config(config_list[[1]][anchor_idx, , drop = FALSE])
  
  # Initialize convergence distance.
  delta <- Inf
  
  # Initialize iteration counter.
  iter <- 0
  
  # Iterate until convergence or until the maximum number of iterations is reached.
  while (iter < max_iter && delta > tol) {
    iter <- iter + 1
    
    # Fit each configuration to the current reference using anchor points only.
    fits <- lapply(config_list, function(cfg) {
      procrustes_fit(
        source_anchor = cfg[anchor_idx, , drop = FALSE],
        target_anchor = ref,
        source_full = cfg,
        allow_scaling = TRUE
      )
    })
    
    # Extract aligned full configurations.
    aligned_configs <- lapply(fits, `[[`, "coords")
    
    # Extract aligned anchors.
    aligned_anchors <- lapply(fits, `[[`, "anchor")
    
    # Compute the new mean anchor configuration across judges.
    new_ref <- Reduce(`+`, aligned_anchors) / J
    
    # Normalize the new reference.
    new_ref <- normalize_config(new_ref)
    
    # Measure change relative to the previous reference.
    delta <- sqrt(sum((new_ref - ref)^2))
    
    # Update the reference configuration.
    ref <- new_ref
  }
  
  # Return the GPA result.
  list(
    aligned_configs = aligned_configs,
    aligned_anchors = aligned_anchors,
    consensus = ref,
    iterations = iter,
    converged = delta <= tol
  )
}

# ------------------------------------------------------------
# 3. One simulation for one condition
# ------------------------------------------------------------

# Simulate one dataset for a single design condition.
simulate_one_dataset <- function(n_judges,
                                 n_unique,
                                 n_anchors,
                                 fit_dim,
                                 noise_sd,
                                 itmax_mds) {
  
  # Indices of the anchor rows within each judge-specific configuration.
  anchor_idx <- seq_len(n_anchors)
  
  # Generate the latent anchor coordinates directly in the fitted dimensionality.
  true_anchors <- matrix(rnorm(n_anchors * fit_dim), nrow = n_anchors, ncol = fit_dim)
  
  # Generate judge-specific latent coordinates for the unique targets.
  true_unique_list <- replicate(
    n_judges,
    matrix(rnorm(n_unique * fit_dim), nrow = n_unique, ncol = fit_dim),
    simplify = FALSE
  )
  
  # Combine common anchors and judge-specific unique targets for each judge.
  true_config_list <- lapply(true_unique_list, function(u) rbind(true_anchors, u))
  
  # Build one noisy dissimilarity matrix for each judge.
  diss_matrices <- lapply(true_config_list, make_noisy_delta, noise_sd = noise_sd)
  
  # Fit each judge's MDS solution directly in fit_dim dimensions.
  mds_fits <- lapply(diss_matrices, fit_mds_safe, ndim = fit_dim, itmax = itmax_mds)
  
  # Track which MDS fits succeeded.
  ok <- vapply(mds_fits, function(x) !is.null(x), logical(1))
  
  # If fewer than 3 fits succeed, return a failed simulation record.
  if (sum(ok) < 3) {
    return(tibble(
      fit_dim = fit_dim,
      n_anchors = n_anchors,
      cor_global_dist = NA_real_,
      rmse_global = NA_real_,
      anchor_rmse = NA_real_,
      anchor_spread = NA_real_,
      mean_stress_final = NA_real_,
      prop_final_stress_lt_10 = NA_real_,
      stable = FALSE,
      gpa_converged = FALSE,
      gpa_iterations = NA_real_,
      usable_judges = sum(ok)
    ))
  }
  
  # Extract estimated configurations from successful MDS fits only.
  est_configs <- lapply(mds_fits[ok], function(x) as.matrix(x$conf))
  
  # Keep the matching true unique configurations for successful judges only.
  true_unique_kept <- true_unique_list[ok]
  
  # Run GPA using only the anchor rows.
  gpa <- partial_gpa(
    config_list = est_configs,
    anchor_idx = anchor_idx
  )
  
  # Prepare objects for evaluation against the true latent space.
  true_anchors_eval <- true_anchors
  true_unique_eval <- true_unique_kept
  aligned_eval_input <- gpa$aligned_configs
  consensus_eval <- gpa$consensus
  
  # Align the GPA consensus anchor configuration to the true anchor configuration.
  eval_tr <- procrustes_fit(
    source_anchor = consensus_eval,
    target_anchor = true_anchors_eval,
    source_full = consensus_eval,
    allow_scaling = TRUE
  )
  
  # Apply that evaluation transformation to all aligned judge configurations.
  aligned_eval <- lapply(aligned_eval_input, apply_procrustes, tr = eval_tr)
  
  # Extract aligned anchor configurations for all retained judges.
  aligned_anchors_eval <- lapply(aligned_eval, function(x) x[anchor_idx, , drop = FALSE])
  
  # Compute the mean recovered anchor configuration across judges.
  anchor_mean <- Reduce(`+`, aligned_anchors_eval) / length(aligned_anchors_eval)
  
  # Stack aligned anchor configurations into a 3D array.
  anchor_array <- array(
    unlist(aligned_anchors_eval),
    dim = c(n_anchors, fit_dim, length(aligned_anchors_eval))
  )
  
  # Compute the standard deviation of each anchor coordinate across judges.
  anchor_sd <- apply(anchor_array, c(1, 2), stats::sd)
  
  # Summarize anchor variability as root mean square spread.
  anchor_spread <- sqrt(mean(anchor_sd^2))
  
  # Build the full true configuration across anchors and all retained unique targets.
  true_global <- rbind(true_anchors_eval, do.call(rbind, true_unique_eval))
  
  # Build the full recovered configuration using the mean anchors plus aligned unique targets.
  rec_global <- rbind(
    anchor_mean,
    do.call(rbind, lapply(aligned_eval, function(x) x[-anchor_idx, , drop = FALSE]))
  )
  
  # Compute the correlation between true and recovered pairwise distances.
  cor_global_dist <- suppressWarnings(stats::cor(vec_dist(true_global), vec_dist(rec_global)))
  
  # Compute overall RMSE between true and recovered global coordinates.
  rmse_global <- sqrt(mean((true_global - rec_global)^2))
  
  # Compute RMSE between the recovered mean anchors and the true anchors.
  anchor_rmse <- sqrt(mean((anchor_mean - true_anchors_eval)^2))
  
  # Compute the mean final stress across successful MDS fits.
  mean_stress_final <- mean(vapply(mds_fits[ok], function(x) x$stress, numeric(1)), na.rm = TRUE)
  
  # Compute the proportion of successful MDS fits with stress below .10.
  prop_final_stress_lt_10 <- mean(vapply(mds_fits[ok], function(x) x$stress < .10, logical(1)), na.rm = TRUE)
  
  # Determine whether this simulation is stable according to all criteria.
  stable <- isTRUE(cor_global_dist >= thr_dist_cor) &&
    isTRUE(anchor_rmse <= thr_anchor_rmse) &&
    isTRUE(anchor_spread <= thr_anchor_spread)
  
  # Return one-row simulation results.
  tibble(
    fit_dim = fit_dim,
    n_anchors = n_anchors,
    cor_global_dist = cor_global_dist,
    rmse_global = rmse_global,
    anchor_rmse = anchor_rmse,
    anchor_spread = anchor_spread,
    mean_stress_final = mean_stress_final,
    prop_final_stress_lt_10 = prop_final_stress_lt_10,
    stable = stable,
    gpa_converged = gpa$converged,
    gpa_iterations = gpa$iterations,
    usable_judges = length(est_configs)
  )
}

# ------------------------------------------------------------
# 4. Simulate one condition
# ------------------------------------------------------------

# Run all simulations for one design condition.
simulate_condition <- function(nsim,
                               n_judges,
                               n_unique,
                               n_anchors,
                               fit_dim,
                               noise_sd,
                               itmax_mds) {
  
  # Preallocate a list to store all simulation outputs.
  sim_list <- vector("list", nsim)
  
  # Store the start time of this design condition.
  start_time <- Sys.time()
  
  # Run all simulations one by one.
  for (s in seq_len(nsim)) {
    sim_list[[s]] <- simulate_one_dataset(
      n_judges = n_judges,
      n_unique = n_unique,
      n_anchors = n_anchors,
      fit_dim = fit_dim,
      noise_sd = noise_sd,
      itmax_mds = itmax_mds
    )
    
    # Compute elapsed time in minutes.
    elapsed_min <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    
    # Estimate remaining time in minutes.
    eta_min <- elapsed_min / s * (nsim - s)
    
    # Print one-line simulation progress.
    cat(
      sprintf(
        "\r  sim %d/%d | elapsed: %.1f min | ETA: %.1f min",
        s, nsim, elapsed_min, eta_min
      )
    )
    flush.console()
  }
  
  # Move to a new line after the progress bar.
  cat("\n")
  
  # Combine all simulation rows into one tibble.
  bind_rows(sim_list)
}

# ------------------------------------------------------------
# 5. Design grid
# ------------------------------------------------------------

# Build the full design grid and keep only admissible conditions.
design_grid <- expand.grid(
  n_judges = n_judges_grid,
  n_unique = n_unique_grid,
  fit_dim = fit_dim_grid,
  n_anchors = anchor_grid
) %>%
  as_tibble() %>%
  arrange(n_judges, n_unique, fit_dim, n_anchors) %>%
  mutate(
    n_targets_per_judge = n_unique + n_anchors,
    admissible = n_anchors >= fit_dim
  ) %>%
  filter(admissible) %>%
  select(-admissible)

# Print the design grid to the console.
print(design_grid)

# Print the total number of conditions to simulate.
cat(sprintf("\nTotal conditions to simulate: %d\n\n", nrow(design_grid)))

# ------------------------------------------------------------
# 6. Run simulation
# ------------------------------------------------------------

# Remove the main checkpoint file if it already exists.
if (file.exists(checkpoint_file)) file.remove(checkpoint_file)

# Remove the backup checkpoint file if it already exists.
if (file.exists(backup_file)) file.remove(backup_file)

# Initialize a list to collect raw condition-level results.
raw_results <- list()

# Store the global start time.
global_start_time <- Sys.time()

# Loop over all design conditions.
for (i in seq_len(nrow(design_grid))) {
  
  # Extract the current row of the design grid.
  row <- design_grid[i, ]
  
  # Print a header for the current condition.
  cat(
    "\n====================================================\n",
    "Condition ", i, "/", nrow(design_grid),
    " | judges = ", row$n_judges,
    " | unique = ", row$n_unique,
    " | fit_dim = ", row$fit_dim,
    " | anchors = ", row$n_anchors,
    "\n====================================================\n",
    sep = ""
  )
  
  # Run all simulations for the current condition.
  cond_res <- simulate_condition(
    nsim = nsim,
    n_judges = row$n_judges,
    n_unique = row$n_unique,
    n_anchors = row$n_anchors,
    fit_dim = row$fit_dim,
    noise_sd = noise_sd,
    itmax_mds = itmax_mds
  ) %>%
    mutate(
      n_judges = row$n_judges,
      n_unique = row$n_unique
    )
  
  # Store the current condition results.
  raw_results[[i]] <- cond_res
  
  # Bind all condition results accumulated so far.
  raw_df <- bind_rows(raw_results)
  
  # Save the current raw results to the main checkpoint file.
  saveRDS(raw_df, checkpoint_file)
  
  # Save the current raw results to the backup checkpoint file.
  saveRDS(raw_df, backup_file)
  
  # Compute total elapsed time in minutes.
  elapsed_total_min <- as.numeric(difftime(Sys.time(), global_start_time, units = "mins"))
  
  # Print a checkpoint confirmation message.
  cat(sprintf(
    "Saved checkpoint after condition %d/%d | total elapsed: %.1f min\n",
    i, nrow(design_grid), elapsed_total_min
  ))
}

# Bind all raw simulation results into a single data frame.
raw_df <- bind_rows(raw_results)

# ------------------------------------------------------------
# 7. Summaries
# ------------------------------------------------------------

# Summarize simulation results at the condition level.
summary_df <- raw_df %>%
  group_by(n_judges, n_unique, fit_dim, n_anchors) %>%
  summarise(
    mean_cor_global_dist = mean(cor_global_dist, na.rm = TRUE),
    sd_cor_global_dist = sd(cor_global_dist, na.rm = TRUE),
    mean_rmse_global = mean(rmse_global, na.rm = TRUE),
    mean_anchor_rmse = mean(anchor_rmse, na.rm = TRUE),
    mean_anchor_spread = mean(anchor_spread, na.rm = TRUE),
    mean_stress_final = mean(mean_stress_final, na.rm = TRUE),
    mean_prop_final_stress_lt_10 = mean(prop_final_stress_lt_10, na.rm = TRUE),
    prop_stable = mean(stable, na.rm = TRUE),
    prop_gpa_converged = mean(gpa_converged, na.rm = TRUE),
    mean_gpa_iterations = mean(gpa_iterations, na.rm = TRUE),
    mean_usable_judges = mean(usable_judges, na.rm = TRUE),
    usable_sims = sum(!is.na(cor_global_dist)),
    .groups = "drop"
  ) %>%
  mutate(
    mean_recovers_well = mean_cor_global_dist >= thr_dist_cor &
      mean_anchor_rmse <= thr_anchor_rmse &
      mean_anchor_spread <= thr_anchor_spread,
    final_status = case_when(
      prop_stable >= thr_prop_stable_good & mean_recovers_well ~ "feasible",
      prop_stable >= thr_prop_stable_mid  | mean_recovers_well ~ "borderline",
      TRUE ~ "unstable"
    )
  ) %>%
  arrange(n_judges, n_unique, fit_dim, n_anchors)

# Print the summary table.
print(summary_df)

# ------------------------------------------------------------
# 8. Save final outputs
# ------------------------------------------------------------

# Uncomment to save the raw results as an RDS file.
# saveRDS(raw_df, "alignment_fixed_dim_raw_results_multicond.rds")

# Uncomment to save the summary results as an RDS file.
# saveRDS(summary_df, "alignment_fixed_dim_summary_results_multicond.rds")

# Uncomment to save the raw results as a CSV file.
# readr::write_csv(raw_df, "alignment_fixed_dim_raw_results_multicond.csv")

# Uncomment to save the summary results as a CSV file.
# readr::write_csv(summary_df, "alignment_fixed_dim_summary_results_multicond.csv")