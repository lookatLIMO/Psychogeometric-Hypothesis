# ============================================================
# MONTE CARLO POWER SIMULATION FOR A LINEAR MIXED MODEL
# ON SIMULATED PAIRWISE SIMILARITY
# ------------------------------------------------------------
#
# Purpose:
# This script estimates statistical power for a linear mixed-effects
# model fitted to simulated pairwise similarity judgments.
#
# Design:
# - 150 judges
# - 120 observations per judge
# - predictors defined as pairwise differences between targets (demographic and HEXACO score)
# - outcome generated from a mixed model with a random intercept for judge
#
# Analyses:
# 1) omnibus test of the HEXACO block via reduced-vs-full model comparison - fcompare()
# 2) test of a single fixed effect (ZH) in the full model - fixed()
#
# Implementation details:
# - continuous predictors are standardized once in the base dataset
# - the simulated outcome is standardized within each replication
# - simulations are split into chunks and parallelized across CPU cores
# - power is estimated with exact binomial confidence intervals
# ============================================================

# -------------------------
# Packages
# -------------------------

# Load lme4 for fitting linear mixed-effects models
library(lme4)

# Load lmerTest for approximate p-values of fixed effects
library(lmerTest)

# Load ggplot2 for plotting power curves
library(ggplot2)

# Load pbapply for progress bars in iterative computations
library(pbapply)

# Load parallel for multi-core execution
library(parallel)

# -------------------------
# Reproducibility
# -------------------------

# Set a global random seed so the initial data-generating setup is reproducible
set.seed(1)

# Use a timer-style progress bar with estimated time remaining
pboptions(type = "timer")

# -------------------------
# Simulation design
# -------------------------

# Number of judges
n_judges <- 150

# Number of observations contributed by each judge
n_obs_per_judge <- 120

# Total number of observations in the simulated dataset
N <- n_judges * n_obs_per_judge

# Generate a raw matrix for the first target in each pair
X <- matrix(rnorm(N * 9), nrow = N, ncol = 9)

# Generate a raw matrix for the second target in each pair
Y <- matrix(rnorm(N * 9), nrow = N, ncol = 9)

# Assign variable names to both target matrices
colnames(X) <- colnames(Y) <- c(
  "gender", "age", "context",
  "H", "E", "X", "A", "C", "O"
)

# -------------------------
# Demographic variables
# -------------------------

# Convert the first target's gender into a binary factor with levels F and M
X[, "gender"] <- cut(
  X[, "gender"],
  breaks = c(-Inf, 0, Inf),
  labels = c("F", "M")
)

# Convert the second target's gender into a binary factor with levels F and M
Y[, "gender"] <- cut(
  Y[, "gender"],
  breaks = c(-Inf, 0, Inf),
  labels = c("F", "M")
)

# Convert the first target's age into 9 ordinal categories based on normal quantiles
X[, "age"] <- cut(
  X[, "age"],
  breaks = qnorm(seq(0, 1, length.out = 10)),
  include.lowest = TRUE,
  labels = FALSE
)

# Convert the second target's age into 9 ordinal categories based on normal quantiles
Y[, "age"] <- cut(
  Y[, "age"],
  breaks = qnorm(seq(0, 1, length.out = 10)),
  include.lowest = TRUE,
  labels = FALSE
)

# Convert the first target's context into 6 ordinal categories
X[, "context"] <- cut(
  X[, "context"],
  breaks = qnorm(seq(0, 1, length.out = 7)),
  include.lowest = TRUE,
  labels = FALSE
)

# Convert the second target's context into 6 ordinal categories
Y[, "context"] <- cut(
  Y[, "context"],
  breaks = qnorm(seq(0, 1, length.out = 7)),
  include.lowest = TRUE,
  labels = FALSE
)

# -------------------------
# Pairwise predictors
# -------------------------

# Allocate a matrix that will store all pairwise predictors
Z <- matrix(NA, nrow = N, ncol = 9)

# Name the pairwise predictors
colnames(Z) <- c(
  "Zgender", "Zage", "Zcontext",
  "ZH", "ZE", "ZX", "ZA", "ZC", "ZO"
)

# Code whether the two targets differ in gender: 1 = different, 0 = same
Z[, "Zgender"] <- as.numeric(X[, "gender"] != Y[, "gender"])

# Compute the absolute difference in ordinal age category
Z[, "Zage"] <- abs(as.numeric(X[, "age"]) - as.numeric(Y[, "age"]))

# Code whether the two targets share the same context: 1 = same, 0 = different
Z[, "Zcontext"] <- as.numeric(X[, "context"] == Y[, "context"])

# Compute absolute pairwise differences on the six HEXACO traits
Z[, "ZH"] <- abs(X[, "H"] - Y[, "H"])
Z[, "ZE"] <- abs(X[, "E"] - Y[, "E"])
Z[, "ZX"] <- abs(X[, "X"] - Y[, "X"])
Z[, "ZA"] <- abs(X[, "A"] - Y[, "A"])
Z[, "ZC"] <- abs(X[, "C"] - Y[, "C"])
Z[, "ZO"] <- abs(X[, "O"] - Y[, "O"])

# Create the judge identifier, repeating each judge for their number of observations
judges <- factor(rep(1:n_judges, each = n_obs_per_judge))

# Combine judge IDs and predictors into the working dataset
dat <- data.frame(judges, Z)

# -------------------------
# Standardize continuous predictors
# -------------------------

# List the predictors that should be z-standardized once in the base dataset
continuous_vars <- c("Zage", "ZH", "ZE", "ZX", "ZA", "ZC", "ZO")

# Standardize those predictors to mean 0 and standard deviation 1
dat[continuous_vars] <- lapply(dat[continuous_vars], function(v) as.numeric(scale(v)))

# -------------------------
# Fixed data-generating parameters
# -------------------------

# Set the intercept of the latent outcome model
beta_intercept <- 0.00

# Set the demographic effect of gender mismatch
beta_gender <- -0.02

# Set the demographic effect of age difference
beta_age <- -0.02

# Set the demographic effect of sharing the same context
beta_context <- 0.02

# Set the standard deviation of the random intercept for judge
sd_judge <- 1

# Set the residual standard deviation of the latent outcome
sd_resid <- 1

# -------------------------
# Helper functions
# -------------------------

# Split a total number of simulations into chunks of approximately equal size
split_nsim <- function(nsim_total, n_chunks) {
  base <- nsim_total %/% n_chunks
  rem  <- nsim_total %% n_chunks
  out <- rep(base, n_chunks)
  if (rem > 0) out[seq_len(rem)] <- out[seq_len(rem)] + 1L
  out[out > 0]
}

# Simulate one standardized dataset for a given common HEXACO effect size
simulate_dataset_std <- function(effect_size) {
  
  # Convert the positive effect-size input into a negative similarity effect:
  # larger differences imply lower similarity
  b_hex <- -effect_size
  
  # Draw one random intercept per judge
  u_judge <- rnorm(n_judges, mean = 0, sd = sd_judge)
  
  # Build the fixed-effects linear predictor
  eta <- beta_intercept +
    beta_gender  * dat$Zgender +
    beta_age     * dat$Zage +
    beta_context * dat$Zcontext +
    b_hex * dat$ZH +
    b_hex * dat$ZE +
    b_hex * dat$ZX +
    b_hex * dat$ZA +
    b_hex * dat$ZC +
    b_hex * dat$ZO
  
  # Add judge-specific random intercepts and residual error to obtain the raw outcome
  y_raw <- eta +
    u_judge[as.integer(dat$judges)] +
    rnorm(nrow(dat), mean = 0, sd = sd_resid)
  
  # Copy the predictor dataset
  dat_sim <- dat
  
  # Standardize the outcome within this replication
  dat_sim$y <- as.numeric(scale(y_raw))
  
  # Return the simulated dataset
  dat_sim
}

# Run a single simulation replicate and return whether the focal test is significant
run_one_rep <- function(effect_size, mode = c("block", "single")) {
  
  # Match the requested testing mode
  mode <- match.arg(mode)
  
  # Simulate one dataset under the requested effect size
  dat_sim <- simulate_dataset_std(effect_size)
  
  # Omnibus block test: compare reduced and full models
  if (mode == "block") {
    
    # Fit the reduced model with demographics only
    fit_red <- lme4::lmer(
      y ~ Zgender + Zage + Zcontext + (1 | judges),
      data = dat_sim,
      REML = FALSE
    )
    
    # Fit the full model with demographics plus all HEXACO difference predictors
    fit_full <- lme4::lmer(
      y ~ Zgender + Zage + Zcontext + ZH + ZE + ZX + ZA + ZC + ZO + (1 | judges),
      data = dat_sim,
      REML = FALSE
    )
    
    # Extract the likelihood-ratio-test p-value for the added HEXACO block
    p_val <- anova(fit_red, fit_full)$`Pr(>Chisq)`[2]
    
  } else {
    
    # Single-effect test: fit the full model and test the ZH coefficient
    fit_full <- lmerTest::lmer(
      y ~ Zgender + Zage + Zcontext + ZH + ZE + ZX + ZA + ZC + ZO + (1 | judges),
      data = dat_sim,
      REML = FALSE
    )
    
    # Extract the p-value associated with ZH
    p_val <- summary(fit_full)$coefficients["ZH", "Pr(>|t|)"]
  }
  
  # Return 1 if the test is finite, not missing, and significant at alpha = .05
  as.integer(is.finite(p_val) && !is.na(p_val) && p_val < 0.05)
}

# Run all simulations for one effect size, parallelizing over chunks of nsim
run_one_effect_parallel_nsim <- function(effect_size,
                                         nsim_total = 1000,
                                         mode = c("block", "single"),
                                         seed_base = 1000,
                                         chunk_multiplier = 3L) {
  
  # Match the requested testing mode
  mode <- match.arg(mode)
  
  # Store local copies of the inputs for safe export to workers
  effect_size_local <- unname(effect_size)
  mode_local <- mode
  seed_local <- seed_base
  
  # Define how many chunks to create, capped by the total number of simulations
  n_chunks <- min(nsim_total, max(n_cores, n_cores * chunk_multiplier))
  
  # Split the total number of simulations across chunks
  chunk_sizes <- split_nsim(nsim_total, n_chunks)
  
  # Run all chunks in parallel across the cluster
  chunk_res <- pblapply(seq_along(chunk_sizes), cl = cl, FUN = function(j) {
    
    # Set a chunk-specific seed for reproducibility
    set.seed(seed_local + j)
    
    # Number of simulations assigned to this chunk
    nsim_chunk <- chunk_sizes[j]
    
    # Counter for successful detections in this chunk
    x_success <- 0L
    
    # Counter for completed non-error simulations in this chunk
    n_done <- 0L
    
    # Run each simulation replicate inside the chunk
    for (k in seq_len(nsim_chunk)) {
      
      # Run one replication and convert failures to NA instead of stopping the whole job
      out <- tryCatch(
        suppressWarnings(run_one_rep(effect_size_local, mode_local)),
        error = function(e) NA_integer_
      )
      
      # Update counters only if the replication completed successfully
      if (!is.na(out)) {
        x_success <- x_success + out
        n_done <- n_done + 1L
      }
    }
    
    # Return chunk-level counts of significant results and usable replications
    data.frame(
      x = x_success,
      n = n_done
    )
  })
  
  # Combine all chunk-level results into one data frame
  chunk_res <- do.call(rbind, chunk_res)
  
  # Total number of significant detections across all chunks
  total_x <- sum(chunk_res$x)
  
  # Total number of usable simulation replications across all chunks
  total_n <- sum(chunk_res$n)
  
  # Compute an exact binomial confidence interval for power
  bt <- binom.test(total_x, total_n)
  
  # Return the power summary for this effect size
  data.frame(
    std_beta_hexaco = effect_size_local,
    signed_std_beta_hexaco = -effect_size_local,
    power = total_x / total_n,
    power_pct = 100 * total_x / total_n,
    ci_lower = unname(bt$conf.int[1]),
    ci_upper = unname(bt$conf.int[2]),
    nsim = total_n
  )
}

# Run a full grid of effect sizes sequentially, parallelizing nsim within each effect size
run_effect_grid_parallel_nsim <- function(effect_grid,
                                          nsim_total = 1000,
                                          mode = c("block", "single"),
                                          seed_offset = 1000,
                                          chunk_multiplier = 3L) {
  
  # Match the requested testing mode
  mode <- match.arg(mode)
  
  # Preallocate a list to store one result table per effect size
  out <- vector("list", length(effect_grid))
  
  # Loop over all effect sizes in the requested grid
  for (i in seq_along(effect_grid)) {
    
    # Add a blank line for readability in console output
    cat("\n")
    
    # Print the current effect size being processed
    cat("Running", mode, "analysis for standardized beta =", effect_grid[i], "\n")
    
    # Run the full Monte Carlo analysis for the current effect size
    out[[i]] <- run_one_effect_parallel_nsim(
      effect_size = effect_grid[i],
      nsim_total = nsim_total,
      mode = mode,
      seed_base = seed_offset + i * 10000L,
      chunk_multiplier = chunk_multiplier
    )
  }
  
  # Bind all effect-specific results into one power table
  power_table <- do.call(rbind, out)
  
  # Keep only effect sizes reaching at least 80% power
  min_b_80 <- power_table[power_table$power >= 0.80, , drop = FALSE]
  
  # Order them from smallest to largest effect size
  min_b_80 <- min_b_80[order(min_b_80$std_beta_hexaco), , drop = FALSE]
  
  # Extract the smallest detectable effect size with at least 80% power
  if (nrow(min_b_80) > 0) {
    min_detectable_80 <- min_b_80[1, , drop = FALSE]
  } else {
    min_detectable_80 <- NA
  }
  
  # Return the full power table and the minimum detectable effect
  list(
    power_table = power_table,
    min_detectable_80 = min_detectable_80
  )
}

# -------------------------
# Parallel setup
# -------------------------

# Detect the number of available CPU cores and leave one core free
n_cores <- max(1L, detectCores() - 1L)

# Create a PSOCK cluster using the selected number of cores
cl <- makeCluster(n_cores)

# Ensure the cluster is stopped automatically when the script exits
on.exit(stopCluster(cl), add = TRUE)

# Initialize parallel random-number streams across workers
clusterSetRNGStream(cl, 1)

# Load required packages on each worker
invisible(clusterEvalQ(cl, {
  library(lme4)
  library(lmerTest)
  NULL
}))

# Export all required objects and functions to the worker processes
clusterExport(
  cl,
  varlist = c(
    "dat", "n_judges",
    "beta_intercept", "beta_gender", "beta_age", "beta_context",
    "sd_judge", "sd_resid",
    "split_nsim", "simulate_dataset_std", "run_one_rep",
    "n_cores"
  ),
  envir = environment()
)

# =====================================================
# 1) Omnibus block test
# =====================================================

# Define the effect-size grid for the omnibus HEXACO block analysis
hexaco_grid_block <- seq(0.010, 0.015, by = 0.001)

# Total number of simulation replications per effect size
nsim_val_block <- 1000

# Run the power analysis for the omnibus block test
block_results <- run_effect_grid_parallel_nsim(
  effect_grid = hexaco_grid_block,
  nsim_total = nsim_val_block,
  mode = "block",
  seed_offset = 100,
  chunk_multiplier = 3L
)

# Extract the full power table for the omnibus analysis
power_table_block <- block_results$power_table

# Extract the minimum detectable standardized effect achieving 80% power
min_detectable_80_block <- block_results$min_detectable_80

# Print the omnibus power table
print(power_table_block)

# Print the minimum detectable omnibus effect size
print(min_detectable_80_block)

# Plot the omnibus power curve
p_block <- ggplot(power_table_block, aes(x = std_beta_hexaco, y = power)) +
  geom_line(color = "#2C7FB8", linewidth = 1) +
  geom_point(color = "#2C7FB8", size = 2) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.0005, alpha = 0.5) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "firebrick") +
  labs(
    x = "Standardized HEXACO beta (same beta for all six traits)",
    y = "Power",
    title = "Power for the omnibus HEXACO block effect"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal(base_size = 12)

# Print the omnibus plot
print(p_block)

# =====================================================
# 2) Single fixed effect in the full model
# =====================================================

# Define the effect-size grid for the single-coefficient analysis
hexaco_grid_single <- seq(0.015, 0.025, by = 0.001)

# Total number of simulation replications per effect size
nsim_val_single <- 1000

# Run the power analysis for the single fixed effect test on ZH
single_results <- run_effect_grid_parallel_nsim(
  effect_grid = hexaco_grid_single,
  nsim_total = nsim_val_single,
  mode = "single",
  seed_offset = 500,
  chunk_multiplier = 3L
)

# Extract the full power table for the single-effect analysis
power_table_single <- single_results$power_table

# Extract the minimum detectable standardized effect achieving 80% power
min_detectable_80_single <- single_results$min_detectable_80

# Print the single-effect power table
print(power_table_single)

# Print the minimum detectable single-effect size
print(min_detectable_80_single)

# Plot the single-effect power curve
p_single <- ggplot(power_table_single, aes(x = std_beta_hexaco, y = power)) +
  geom_line(color = "#1B9E77", linewidth = 1) +
  geom_point(color = "#1B9E77", size = 2) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.0005, alpha = 0.5) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "firebrick") +
  labs(
    x = "Standardized HEXACO beta (same beta for all six traits)",
    y = "Power",
    title = "Power to detect a single standardized HEXACO effect"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal(base_size = 12)

# Print the single-effect plot
print(p_single)