# ============================================================
# MONTE CARLO POWER SIMULATION FOR A LLM ON RAW SIMILARITIES
# ------------------------------------------------------------

## =====================================================
## Power analysis with parallelization over nsim
## 1) Omnibus HEXACO block effect: fcompare()
## 2) Minimum detectable fixed effect in the full model: fixed()
##
## This script performs a simulation-based power analysis for linear
## mixed-effects models using simr. It simulates pairwise similarity
## data for 150 judges with 120 observations each (pairing 16 targets), 
## defines predictors based on demographic and HEXACO score differences between targets,
## and estimates power under two scenarios:
## (1) an omnibus test of the incremental contribution of the HEXACO
##     block beyond demographics using fcompare();
## (2) a test of the minimum detectable individual fixed effect in the
##     full model using fixed().
## To reduce computation time, the script parallelizes simulations
## across nsim chunks, aggregates results across workers, computes
## confidence intervals for estimated power, and produces summary
## tables and plots for each analysis.
## =====================================================

# Load the package for fitting linear mixed-effects models
library(lme4)

# Load the package for simulation-based power analysis in mixed models
library(simr)

# Load the package for plotting results
library(ggplot2)

# Load the package for progress bars with apply-style functions
library(pbapply)

# Load the package for parallel computation
library(parallel)

# Set the random seed to ensure reproducibility
set.seed(1)

## Avoid nested progress bars from simr
# Disable simr's internal progress bar to prevent overlap with pbapply's progress bar
simrOptions(progress = FALSE)

## Progress bar with ETA
# Use a timer-style progress bar that shows elapsed and estimated remaining time
pboptions(type = "timer")

## -------------------------
## Design
## -------------------------

# Define the number of judges
n_judges <- 150

# Define the number of observations per judge
n_obs_per_judge <- 120

# Compute the total number of observations
N <- n_judges * n_obs_per_judge

# Generate a matrix of simulated characteristics for the first target in each pair
X <- matrix(rnorm(N * 9), nrow = N, ncol = 9)

# Generate a matrix of simulated characteristics for the second target in each pair
Y <- matrix(rnorm(N * 9), nrow = N, ncol = 9)

# Assign variable names to the columns of both matrices
colnames(X) <- colnames(Y) <- c(
  "gender", "age", "context",
  "H", "E", "X", "A", "C", "O"
)

## Demographic variables

# Convert the simulated gender variable in X into a binary categorical variable
X[, "gender"] <- cut(X[, "gender"], breaks = c(-Inf, 0, Inf), labels = c("F", "M"))

# Convert the simulated gender variable in Y into a binary categorical variable
Y[, "gender"] <- cut(Y[, "gender"], breaks = c(-Inf, 0, Inf), labels = c("F", "M"))

# Convert the simulated age variable in X into 9 ordered age categories
X[, "age"] <- cut(
  X[, "age"],
  breaks = qnorm(seq(0, 1, length.out = 10)),
  include.lowest = TRUE,
  labels = FALSE
)

# Convert the simulated age variable in Y into 9 ordered age categories
Y[, "age"] <- cut(
  Y[, "age"],
  breaks = qnorm(seq(0, 1, length.out = 10)),
  include.lowest = TRUE,
  labels = FALSE
)

# Convert the simulated context variable in X into 6 categorical groups
X[, "context"] <- cut(
  X[, "context"],
  breaks = qnorm(seq(0, 1, length.out = 7)),
  include.lowest = TRUE,
  labels = FALSE
)

# Convert the simulated context variable in Y into 6 categorical groups
Y[, "context"] <- cut(
  Y[, "context"],
  breaks = qnorm(seq(0, 1, length.out = 7)),
  include.lowest = TRUE,
  labels = FALSE
)

## Predictors

# Create an empty matrix that will contain pairwise difference/similarity predictors
Z <- matrix(NA, nrow = N, ncol = 9)

# Assign variable names to the columns of the predictor matrix
colnames(Z) <- c("Zgender", "Zage", "Zcontext", "ZH", "ZE", "ZX", "ZA", "ZC", "ZO")

# Code gender dissimilarity as 1 when the two targets differ in gender, 0 otherwise
Z[, "Zgender"]  <- as.numeric(X[, "gender"] != Y[, "gender"])   # 1 = different gender

# Compute absolute age difference between the two targets
Z[, "Zage"]     <- abs(as.numeric(X[, "age"]) - as.numeric(Y[, "age"]))

# Code contextual similarity as 1 when the two targets belong to the same context, 0 otherwise
Z[, "Zcontext"] <- as.numeric(X[, "context"] == Y[, "context"]) # 1 = same context

# Compute the absolute difference in Honesty-Humility between the two targets
Z[, "ZH"] <- abs(X[, "H"] - Y[, "H"])

# Compute the absolute difference in Emotionality between the two targets
Z[, "ZE"] <- abs(X[, "E"] - Y[, "E"])

# Compute the absolute difference in Extraversion between the two targets
Z[, "ZX"] <- abs(X[, "X"] - Y[, "X"])

# Compute the absolute difference in Agreeableness between the two targets
Z[, "ZA"] <- abs(X[, "A"] - Y[, "A"])

# Compute the absolute difference in Conscientiousness between the two targets
Z[, "ZC"] <- abs(X[, "C"] - Y[, "C"])

# Compute the absolute difference in Openness between the two targets
Z[, "ZO"] <- abs(X[, "O"] - Y[, "O"])

# Create a factor identifying the judge associated with each observation
judges <- factor(rep(1:n_judges, each = n_obs_per_judge))

# Combine the judge identifier and all predictors into a single data frame
dat <- data.frame(judges, Z)

## -------------------------
## Fixed parameters kept constant
## -------------------------

# Set the intercept of the simulated model
beta_intercept <-  0.00

# Set the effect of gender dissimilarity
beta_gender    <- -0.02

# Set the effect of age dissimilarity
beta_age       <- -0.02

# Set the effect of contextual similarity
beta_context   <-  0.02

# Set the random intercept variance
V1 <- 1

# Set the residual standard deviation
s  <- 1

## -------------------------
## Parallel setup
## -------------------------

# Detect the number of available CPU cores and leave one core free
n_cores <- max(1L, detectCores() - 1L)

# Create a parallel socket cluster using the selected number of cores
cl <- makeCluster(n_cores)

# Ensure that the cluster is automatically stopped when the script exits
on.exit(stopCluster(cl), add = TRUE)

# Set reproducible random number streams across cluster workers
clusterSetRNGStream(cl, 1)

# Load the required packages on each worker node
invisible(clusterEvalQ(cl, {
  library(lme4)
  library(simr)
  NULL
}))

# Export the necessary objects from the current environment to each worker node
clusterExport(
  cl,
  varlist = c(
    "dat",
    "V1", "s",
    "beta_intercept", "beta_gender", "beta_age", "beta_context"
  ),
  envir = environment()
)

## -------------------------
## Helper functions
## -------------------------

# Define a function to divide the total number of simulations into chunks
split_nsim <- function(nsim_total, n_chunks) {
  
  # Compute the base number of simulations per chunk
  base <- nsim_total %/% n_chunks
  
  # Compute the remainder to be distributed across the first chunks
  rem  <- nsim_total %% n_chunks
  
  # Create a vector with the base number of simulations for each chunk
  out <- rep(base, n_chunks)
  
  # Distribute the remainder across the first chunks
  if (rem > 0) out[seq_len(rem)] <- out[seq_len(rem)] + 1L
  
  # Return only positive chunk sizes
  out[out > 0]
}

# Define a function to run the power analysis for one effect size,
# parallelizing across simulation chunks
run_one_effect_parallel_nsim <- function(effect_size,
                                         nsim_total = 1000,
                                         mode = c("block", "single"),
                                         seed_base = 1000,
                                         chunk_multiplier = 3L) {
  
  # Match the requested analysis mode
  mode <- match.arg(mode)
  
  # Store the negative version of the tested effect size
  b_hex_local   <- -effect_size
  
  # Store the selected analysis mode locally
  mode_local    <- mode
  
  # Store the base seed locally
  seed_local    <- seed_base
  
  # Determine how many chunks to create
  n_chunks      <- min(nsim_total, max(n_cores, n_cores * chunk_multiplier))
  
  # Split the total number of simulations into chunks
  chunk_sizes   <- split_nsim(nsim_total, n_chunks)
  
  # Run the chunks in parallel using pbapply
  chunk_res <- pblapply(seq_along(chunk_sizes), cl = cl, FUN = function(j) {
    
    # Extract the number of simulations assigned to the current chunk
    nsim_chunk <- chunk_sizes[j]
    
    # Define the fixed-effect vector for the simulated model
    b <- c(
      beta_intercept,
      beta_gender,
      beta_age,
      beta_context,
      b_hex_local, b_hex_local, b_hex_local,
      b_hex_local, b_hex_local, b_hex_local
    )
    
    # Build the full mixed-effects model used for data simulation
    model_full <- makeLmer(
      y ~ Zgender + Zage + Zcontext + ZH + ZE + ZX + ZA + ZC + ZO + (1 | judges),
      fixef   = b,
      VarCorr = V1,
      sigma   = s,
      data    = dat
    )
    
    # If the selected mode is "block", test the joint contribution of HEXACO traits
    if (mode_local == "block") {
      ps <- powerSim(
        model_full,
        test = fcompare(~ Zgender + Zage + Zcontext, method = "lr"),
        nsim = nsim_chunk,
        seed = seed_local + j,
        progress = FALSE
      )
    } else {
      # Otherwise, test a single fixed effect within the full model
      ps <- powerSim(
        model_full,
        test = fixed("ZH", method = "t"),
        nsim = nsim_chunk,
        seed = seed_local + j,
        progress = FALSE
      )
    }
    
    # Return the number of successful detections and the number of simulations for this chunk
    data.frame(
      x = ps$x,
      n = ps$n
    )
  })
  
  # Combine the results from all chunks into a single data frame
  chunk_res <- do.call(rbind, chunk_res)
  
  # Sum the number of successful detections across chunks
  total_x <- sum(chunk_res$x)
  
  # Sum the total number of simulations across chunks
  total_n <- sum(chunk_res$n)
  
  # Compute the binomial confidence interval for the estimated power
  bt <- binom.test(total_x, total_n)
  
  # Return a summary data frame for the current effect size
  data.frame(
    abs_b_hexaco    = effect_size,
    signed_b_hexaco = -effect_size,
    power           = total_x / total_n,
    power_pct       = 100 * total_x / total_n,
    ci_lower        = unname(bt$conf.int[1]),
    ci_upper        = unname(bt$conf.int[2]),
    nsim            = total_n
  )
}

# Define a function to run the power analysis across a grid of effect sizes
run_effect_grid_parallel_nsim <- function(effect_grid,
                                          nsim_total = 1000,
                                          mode = c("block", "single"),
                                          seed_offset = 1000,
                                          chunk_multiplier = 3L) {
  
  # Match the requested analysis mode
  mode <- match.arg(mode)
  
  # Create an empty list to store the results for each effect size
  out <- vector("list", length(effect_grid))
  
  # Loop over all effect sizes in the grid
  for (i in seq_along(effect_grid)) {
    
    # Print a blank line for readability
    cat("\n")
    
    # Print a message indicating which effect size is currently being processed
    cat("Running", mode, "analysis for |b| =", effect_grid[i], "\n")
    
    # Run the power analysis for the current effect size
    out[[i]] <- run_one_effect_parallel_nsim(
      effect_size = effect_grid[i],
      nsim_total = nsim_total,
      mode = mode,
      seed_base = seed_offset + i * 10000L,
      chunk_multiplier = chunk_multiplier
    )
  }
  
  # Combine all effect-size-specific results into one data frame
  power_table <- do.call(rbind, out)
  
  # Select the effect sizes that reach at least 80% power
  min_b_80 <- power_table[power_table$power >= 0.80, , drop = FALSE]
  
  # Order those effect sizes from smallest to largest
  min_b_80 <- min_b_80[order(min_b_80$abs_b_hexaco), , drop = FALSE]
  
  # Extract the smallest effect size achieving at least 80% power
  if (nrow(min_b_80) > 0) {
    min_detectable_80 <- min_b_80[1, , drop = FALSE]
  } else {
    # Return NA if no tested effect size reaches 80% power
    min_detectable_80 <- NA
  }
  
  # Return both the full power table and the minimum detectable effect size
  list(
    power_table = power_table,
    min_detectable_80 = min_detectable_80
  )
}

## =====================================================
## 1) Omnibus block test: fcompare()
## Parallelized over nsim within each effect size
## =====================================================

# Define the effect-size grid for the omnibus block test
hexaco_grid_block <- seq(0.010, 0.015, by = 0.001)

# Define the total number of simulations for each effect size
nsim_val_block <- 1000

# Run the omnibus power analysis across the selected effect sizes
block_results <- run_effect_grid_parallel_nsim(
  effect_grid = hexaco_grid_block,
  nsim_total = nsim_val_block,
  mode = "block",
  seed_offset = 100,
  chunk_multiplier = 3L
)

# Extract the full power table for the omnibus block test
power_table_block <- block_results$power_table

# Extract the smallest effect size that reaches at least 80% power
min_detectable_80_block <- block_results$min_detectable_80

# Print the omnibus power table
print(power_table_block)

# Print the minimum detectable omnibus effect size
print(min_detectable_80_block)

# Build a plot of power as a function of effect size for the omnibus block test
p_block <- ggplot(power_table_block, aes(x = abs_b_hexaco, y = power)) +
  geom_line(color = "#2C7FB8", linewidth = 1) +
  geom_point(color = "#2C7FB8", size = 2) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.0003, alpha = 0.5) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "firebrick") +
  labs(
    x = "Absolute HEXACO effect size (same b for all six traits)",
    y = "Power",
    title = "Power for the omnibus HEXACO block effect"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal(base_size = 12)

# Print the omnibus power plot
print(p_block)

## =====================================================
## 2) Single fixed effect in the full model
## Parallelized over nsim within each effect size
## =====================================================

# Define the effect-size grid for the single-coefficient test
hexaco_grid_single <- seq(0.024, 0.025, by = 0.001)

# Define the total number of simulations for each effect size
nsim_val_single <- 1000

# Run the single fixed-effect power analysis across the selected effect sizes
single_results <- run_effect_grid_parallel_nsim(
  effect_grid = hexaco_grid_single,
  nsim_total = nsim_val_single,
  mode = "single",
  seed_offset = 500,
  chunk_multiplier = 3L
)

# Extract the full power table for the single fixed-effect test
power_table_single <- single_results$power_table

# Extract the smallest effect size that reaches at least 80% power
min_detectable_80_single <- single_results$min_detectable_80

# Print the single-effect power table
print(power_table_single)

# Print the minimum detectable single fixed effect
print(min_detectable_80_single)

# Build a plot of power as a function of effect size for the single fixed-effect test
p_single <- ggplot(power_table_single, aes(x = abs_b_hexaco, y = power)) +
  geom_line(color = "#1B9E77", linewidth = 1) +
  geom_point(color = "#1B9E77", size = 2) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.0003, alpha = 0.5) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "firebrick") +
  labs(
    x = "Absolute fixed effect size (same b for all six traits)",
    y = "Power",
    title = "Power to detect a HEXACO fixed effect in the full model"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal(base_size = 12)

# Print the single-effect power plot
print(p_single)