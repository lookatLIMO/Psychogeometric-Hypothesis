## ============================================================
## GPA through anchors + ENVFIT vectors
## ------------------------------------------------------------
## Purpose:
## This script simulates a common psychogeometric space across judges
## using anchor-based alignment of individual ordinal MDS solutions.
##
## General workflow:
## 1. Simulate latent personality profiles for anchors and judge-specific targets.
## 2. Map those profiles into a true 3D latent space using fixed beta weights.
## 3. Generate judge-specific noisy dissimilarity matrices from that space.
## 4. Fit a separate 3D ordinal MDS solution for each judge.
## 5. Align the judge-specific configurations through anchor-based GPA/Procrustes.
## 6. Orient the final common space post hoc so that the axes match the
##    intended composite trait dimensions.
## 7. Compute agreement diagnostics on anchors.
## 8. Draw envfit-style trait vectors directly from beta_true, with small
##    guided perturbations to make the arrows more realistic and readable.
##
## Important note:
## No final LMMs are fitted here. The trait vectors shown in the final plot
## are derived directly from beta_true after post-hoc orientation of the
## common MDS/GPA space.
## ============================================================


## ------------------------------------------------------------
## 0. Packages
## ------------------------------------------------------------

# Install pacman if it is not already available.
if (!require("pacman")) install.packages("pacman")

# Load the required packages.
pacman::p_load(
  smacof, MASS, dplyr, tidyr, psych, plotly
)

# Set a global random seed for reproducibility.
set.seed(123)


## ------------------------------------------------------------
## 1. General parameters
## ------------------------------------------------------------

# Number of judges.
n_judges <- 3

# Total number of targets shown to each judge, including anchors.
n_targets_per_judge <- 15

# Number of common anchors shared across judges.
n_anchors <- 5

# Require at least 3 anchors.
stopifnot(n_anchors >= 3)

# Require at least one judge-specific target beyond the anchors.
stopifnot(n_targets_per_judge > n_anchors)

# Number of unique judge-specific targets.
n_unique <- n_targets_per_judge - n_anchors

# Judge IDs.
judge_ids <- paste0("J", 1:n_judges)

# Anchor IDs.
anchor_ids <- paste0("A", 1:n_anchors)

# Dimensionality of the true latent space.
ndim_true <- 3

# Dimensionality used in MDS fitting.
ndim_mds <- 3

# Names of the latent dimensions.
dim_names <- paste0("D", 1:ndim_true)

# Trait names.
trait_names <- c("H","E","X","A","C","O")

# Number of traits.
n_traits <- length(trait_names)

# Noise added to simulated personality ratings.
rating_noise_sd <- 0.20

# Noise added to pairwise distances before ordinal conversion.
distance_noise_sd <- 0.01

# Simulation parameters.
within_judge_cov_scale <- 0.60
coord_noise <- 0.01
coord_spread <- 2.1
plot_spread <- 2.2


## ------------------------------------------------------------
## 2. Utility functions
## ------------------------------------------------------------

# Convert a numeric vector into 7 ordinal categories using empirical quantiles.
to_ordinal_7 <- function(x) {
  # Compute 8 cut points to create 7 ordinal bins.
  cuts <- quantile(x, probs = seq(0, 1, length.out = 8), na.rm = TRUE)
  
  # If quantile cut points are duplicated, jitter them slightly.
  if (any(duplicated(cuts))) {
    dup <- duplicated(cuts)
    cuts[dup] <- cuts[dup] + seq_len(sum(dup)) * 1e-8
  }
  
  # Return integer scores from 1 to 7.
  as.integer(cut(x, breaks = cuts, include.lowest = TRUE, labels = 1:7))
}

# Generate a random orthogonal rotation matrix.
rand_orth <- function(p) {
  # Compute an orthogonal matrix through QR decomposition.
  Q <- qr.Q(qr(matrix(rnorm(p * p), p, p)))
  
  # Force a proper rotation by fixing negative determinant when needed.
  if (det(Q) < 0) Q[, 1] <- -Q[, 1]
  
  # Return the rotation matrix.
  Q
}

# Convert a continuous dissimilarity matrix into a 9-point ordinal matrix.
to_ordinal_9_matrix <- function(dmat) {
  # Extract the upper-triangular distances.
  upper_vals <- dmat[upper.tri(dmat)]
  
  # Compute quantile cut points for 9 ordinal bins.
  breaks <- quantile(upper_vals, probs = seq(0, 1, length.out = 10), na.rm = TRUE)
  
  # Jitter duplicated cut points if needed.
  if (any(duplicated(breaks))) {
    dup <- duplicated(breaks)
    breaks[dup] <- breaks[dup] + seq_len(sum(dup)) * 1e-8
  }
  
  # Initialize the output matrix.
  out <- matrix(0, nrow(dmat), ncol(dmat))
  
  # Fill the upper triangle with ordinal values from 1 to 9.
  out[upper.tri(out)] <- cut(
    upper_vals,
    breaks = breaks,
    include.lowest = TRUE,
    labels = 1:9
  )
  
  # Mirror the upper triangle to obtain a symmetric matrix.
  out <- out + t(out)
  
  # Set diagonal values to zero.
  diag(out) <- 0
  
  # Store the matrix as numeric.
  storage.mode(out) <- "numeric"
  
  # Copy row and column names from the original matrix.
  rownames(out) <- rownames(dmat)
  colnames(out) <- colnames(dmat)
  
  # Return the ordinal dissimilarity matrix.
  out
}

# Fit a 3D ordinal MDS solution.
fit_mds_3d <- function(delta) {
  smacof::smacofSym(
    delta,
    ndim = 3,
    type = "ordinal",
    ties = "secondary",
    init = "torgerson",
    principal = TRUE,
    itmax = 1000
  )
}

# Select k points that are as spread out as possible from a candidate set.
select_spread_points <- function(candidates, k, first_index = NULL) {
  # Coerce candidates to a matrix.
  candidates <- as.matrix(candidates)
  
  # Number of candidate points.
  n <- nrow(candidates)
  
  # Stop if k is larger than the number of candidates.
  if (k > n) stop("k cannot be greater than the number of candidates.")
  
  # If no first point is supplied, start from the point farthest from the center.
  if (is.null(first_index)) {
    center <- colMeans(candidates)
    d0 <- apply(candidates, 1, function(x) sqrt(sum((x - center)^2)))
    selected <- which.max(d0)
  } else {
    # Otherwise start from the supplied first index.
    selected <- first_index[1]
  }
  
  # Iteratively add the point farthest from the currently selected set.
  while (length(selected) < k) {
    remaining <- setdiff(seq_len(n), selected)
    
    min_d_to_selected <- sapply(remaining, function(i) {
      min(sapply(selected, function(j) {
        sqrt(sum((candidates[i, ] - candidates[j, ])^2))
      }))
    })
    
    next_idx <- remaining[which.max(min_d_to_selected)]
    selected <- c(selected, next_idx)
  }
  
  # Return the selected points.
  candidates[selected, , drop = FALSE]
}

# Align one configuration to a reference using Procrustes transformation on anchors.
procrustes_align <- function(ref_anchor, src_anchor, src_all) {
  # Convert inputs to matrices.
  ref_anchor <- as.matrix(ref_anchor)
  src_anchor <- as.matrix(src_anchor)
  src_all <- as.matrix(src_all)
  
  # Compute anchor centroids.
  ref_mean <- colMeans(ref_anchor)
  src_mean <- colMeans(src_anchor)
  
  # Center anchor coordinates.
  ref_c <- sweep(ref_anchor, 2, ref_mean, "-")
  src_c <- sweep(src_anchor, 2, src_mean, "-")
  
  # Compute the SVD of the anchor cross-product.
  sv <- svd(t(src_c) %*% ref_c)
  
  # Compute the rotation matrix.
  R <- sv$u %*% t(sv$v)
  
  # Enforce a proper rotation if determinant is negative.
  if (det(R) < 0) {
    sv$u[, ncol(sv$u)] <- -sv$u[, ncol(sv$u)]
    R <- sv$u %*% t(sv$v)
  }
  
  # Compute the optimal scaling factor.
  scale_fac <- sum(sv$d) / sum(src_c^2)
  
  # Center the full source configuration.
  src_all_c <- sweep(src_all, 2, src_mean, "-")
  
  # Apply scale and rotation.
  aligned <- scale_fac * src_all_c %*% R
  
  # Translate back to the reference centroid.
  aligned <- sweep(aligned, 2, ref_mean, "+")
  
  # Restore row and column names.
  rownames(aligned) <- rownames(src_all)
  colnames(aligned) <- colnames(src_all)
  
  # Return the aligned full configuration.
  aligned
}

# Run iterative GPA using anchors as the alignment landmarks.
gpa_from_anchors <- function(configs, anchor_ids, max_iter = 1000, tol = 1e-9) {
  # Initialize the reference anchor configuration from the first judge.
  ref <- configs[[1]][anchor_ids, , drop = FALSE]
  
  # Preallocate a list for aligned configurations.
  aligned <- vector("list", length(configs))
  names(aligned) <- names(configs)
  
  # Initialize the convergence delta.
  delta <- NA_real_
  
  # Iterate GPA updates until convergence or max_iter.
  for (iter in 1:max_iter) {
    # Align every configuration to the current reference anchors.
    for (j in seq_along(configs)) {
      aligned[[j]] <- procrustes_align(
        ref_anchor = ref,
        src_anchor = configs[[j]][anchor_ids, , drop = FALSE],
        src_all = configs[[j]]
      )
    }
    
    # Compute the new reference as the average aligned anchor configuration.
    new_ref <- Reduce(
      `+`,
      lapply(aligned, function(m) m[anchor_ids, , drop = FALSE])
    ) / length(aligned)
    
    # Measure the maximum absolute change in the anchors.
    delta <- max(abs(new_ref - ref))
    
    # Stop if the reference has converged.
    if (delta < tol) break
    
    # Update the reference.
    ref <- new_ref
  }
  
  # Return the GPA result object.
  list(
    aligned = aligned,
    consensus_anchors = new_ref,
    iterations = iter,
    delta_final = delta,
    converged = isTRUE(delta < tol)
  )
}

# Simulate one judge's perceived dissimilarities from the true coordinates.
simulate_judge_dissimilarities <- function(target_ids, coords_true, noise_sd = 0.06) {
  # Draw a random rotation.
  R <- rand_orth(ndim_true)
  
  # Draw judge-specific axis weights.
  w <- runif(ndim_true, 0.85, 1.15)
  
  # Extract the true coordinates of this judge's targets.
  base <- as.matrix(coords_true[target_ids, dim_names, drop = FALSE])
  
  # Apply judge-specific rotation.
  perceived <- base %*% R
  
  # Apply judge-specific axis weighting.
  perceived <- sweep(perceived, 2, w, "*")
  
  # Compute continuous Euclidean distances in the perceived space.
  dtrue <- as.matrix(dist(perceived, method = "euclidean"))
  
  # Create a symmetric Gaussian noise matrix.
  noise_mat <- matrix(rnorm(length(dtrue), 0, noise_sd), nrow(dtrue), ncol(dtrue))
  noise_mat <- (noise_mat + t(noise_mat)) / 2
  diag(noise_mat) <- 0
  
  # Add noise to the continuous distances.
  delta_cont <- dtrue + noise_mat
  
  # Truncate negative distances to zero.
  delta_cont[delta_cont < 0] <- 0
  
  # Set diagonal values to zero.
  diag(delta_cont) <- 0
  
  # Add target names.
  rownames(delta_cont) <- target_ids
  colnames(delta_cont) <- target_ids
  
  # Convert continuous distances to 9-point ordinal dissimilarities.
  delta_ord <- to_ordinal_9_matrix(delta_cont)
  
  # Return all judge-specific dissimilarity objects.
  list(
    target_ids = target_ids,
    delta_cont = delta_cont,
    delta_ord = delta_ord,
    perceived = perceived,
    weights = w
  )
}

# Z-standardize a numeric vector and coerce it to plain numeric.
znum <- function(x) as.numeric(scale(x))

# Return all permutations of the three axis indices.
all_axis_perms <- function() {
  rbind(
    c(1, 2, 3),
    c(1, 3, 2),
    c(2, 1, 3),
    c(2, 3, 1),
    c(3, 1, 2),
    c(3, 2, 1)
  )
}


## ------------------------------------------------------------
## 3. Simulate latent personality profiles
## ------------------------------------------------------------

# Define the latent trait covariance matrix.
Sigma <- matrix(c(
  1.00, 0.05, 0.05, 0.45, 0.05, 0.05,
  0.05, 1.00, 0.05, 0.05, 0.45, 0.05,
  0.05, 0.05, 1.00, 0.05, 0.05, 0.45,
  0.45, 0.05, 0.05, 1.00, 0.05, 0.05,
  0.05, 0.45, 0.05, 0.05, 1.00, 0.05,
  0.05, 0.05, 0.45, 0.05, 0.05, 1.00
), 6, 6, byrow = TRUE)

# Name rows and columns of the covariance matrix.
colnames(Sigma) <- rownames(Sigma) <- trait_names

# Build judge-specific target IDs.
unique_ids <- unlist(lapply(judge_ids, function(j) paste0(j, "_T", 1:n_unique)))

# Concatenate anchor IDs and all unique target IDs.
all_ids <- c(anchor_ids, unique_ids)

## ------------------------------------------------------------
## 3A. Ad hoc anchors: 5 canonical person types
## ------------------------------------------------------------

# If there are exactly 5 anchors, use manually designed canonical anchor profiles.
if (n_anchors == 5) {
  
  anchor_profiles <- data.frame(
    id = anchor_ids,
    profile = c(
      "Prosocial_Pole",
      "Plasticity_Pole",
      "Stability_Pole",
      "Low_Prosociality",
      "Low_Stability_Low_Plasticity"
    ),
    H = c( 2.5, -0.5, -0.3, -2.3,  0.0),
    E = c(-0.6, -0.3,  2.4,  0.0, -2.1),
    X = c(-0.4,  2.5, -0.4,  0.2, -1.9),
    A = c( 2.3, -0.4,  0.0, -2.1,  0.1),
    C = c( 0.4, -0.3,  2.3, -0.2, -2.2),
    O = c(-0.3,  2.4, -0.4,  0.0, -1.8)
  )
  
  # Keep the trait columns as the latent anchor profiles.
  anchor_latent <- anchor_profiles[, trait_names]
  rownames(anchor_latent) <- anchor_profiles$id
  
} else {
  
  # Otherwise draw a large candidate pool from a multivariate normal distribution.
  anchor_candidate_pool <- MASS::mvrnorm(
    n = max(100, n_anchors * 20),
    mu = rep(0, n_traits),
    Sigma = Sigma * 1.4
  )
  
  # Select a spread-out subset of candidates as anchors.
  anchor_latent <- select_spread_points(anchor_candidate_pool, k = n_anchors)
  
  # Standardize and stretch the anchor profiles.
  anchor_latent <- as.data.frame(scale(anchor_latent) * 2.5)
  colnames(anchor_latent) <- trait_names
  rownames(anchor_latent) <- anchor_ids
  
  # Create generic anchor profile labels.
  anchor_profiles <- data.frame(
    id = anchor_ids,
    profile = paste0("Anchor_", seq_len(n_anchors))
  )
}

# Print anchor profile information.
cat("\nAnchor profiles:\n")
print(anchor_profiles)

## ------------------------------------------------------------
## 3B. More moderate judge centers
## ------------------------------------------------------------

# If there are 3 judges, use hand-tuned judge trait centers.
if (n_judges == 3) {
  judge_trait_centers <- rbind(
    J1 = c( 0.7, -0.1,  0.3,  0.5,  0.4,  0.1),
    J2 = c(-0.3,  0.4, -0.2,  0.1, -0.4,  0.5),
    J3 = c( 0.1,  0.0,  0.6, -0.3,  0.1, -0.2)
  )
} else {
  # Otherwise sample judge centers from a multivariate normal pool.
  judge_center_pool <- MASS::mvrnorm(
    n = max(50, n_judges * 10),
    mu = rep(0, n_traits),
    Sigma = Sigma
  )
  
  # Select spread-out judge centers.
  judge_trait_centers <- select_spread_points(judge_center_pool, k = n_judges)
  judge_trait_centers <- as.matrix(scale(judge_trait_centers))
}

# Label judge center rows and columns.
rownames(judge_trait_centers) <- judge_ids
colnames(judge_trait_centers) <- trait_names

# Generate judge-specific unique latent personality profiles.
unique_latent_list <- lapply(judge_ids, function(j) {
  tmp <- MASS::mvrnorm(
    n = n_unique,
    mu = judge_trait_centers[j, ],
    Sigma = Sigma * within_judge_cov_scale
  )
  
  tmp <- as.data.frame(tmp)
  colnames(tmp) <- trait_names
  rownames(tmp) <- paste0(j, "_T", 1:n_unique)
  tmp
})

# Combine all judge-specific unique profiles.
unique_latent <- do.call(rbind, unique_latent_list)

# Combine anchors and unique targets into one latent trait table.
traits_latent <- rbind(
  anchor_latent,
  unique_latent
)

# Convert latent traits to ordinal 1-7 ratings.
traits <- as.data.frame(lapply(traits_latent, to_ordinal_7))
traits$id <- rownames(traits_latent)
traits <- traits[, c("id", trait_names)]
rownames(traits) <- traits$id


## ------------------------------------------------------------
## 4. Fixed true betas
## ------------------------------------------------------------

# D1 = prosociality
# D2 = plasticity
# D3 = stability

# Define the true trait-to-dimension mapping matrix.
beta_true <- matrix(c(
  1.00,  0.10,  0.05,  0.90,  0.10,  0.05,   # D1: H, A
  0.05,  0.05,  0.95,  0.05,  0.05,  1.00,   # D2: X, O
  0.05,  0.95,  0.05,  0.05,  0.90,  0.10    # D3: E, C
), nrow = 3, byrow = TRUE)

# Label rows and columns of the beta matrix.
rownames(beta_true) <- dim_names
colnames(beta_true) <- trait_names

# Print the fixed beta matrix used to generate the latent space.
cat("\nFixed beta_true used for latent space generation:\n")
print(beta_true)


## ------------------------------------------------------------
## 5. Build the true 3D latent space
## ------------------------------------------------------------

# Standardize the ordinal trait matrix.
X_traits <- scale(as.matrix(traits[, trait_names]))

# Map traits into the true latent coordinates through beta_true and small noise.
coords_latent <- X_traits %*% t(beta_true) +
  matrix(rnorm(nrow(X_traits) * ndim_true, 0, coord_noise), nrow(X_traits), ndim_true)

# Preserve target IDs as row names.
rownames(coords_latent) <- rownames(traits)

# If there are 3 judges, use hand-tuned judge-specific spatial offsets.
if (n_judges == 3) {
  judge_coord_offsets <- rbind(
    J1 = c( 0.35, -0.20,  0.10),
    J2 = c(-0.25,  0.30, -0.10),
    J3 = c(-0.10, -0.10,  0.25)
  )
} else {
  # Otherwise draw judge-specific offsets from a multivariate normal distribution.
  judge_coord_offsets <- MASS::mvrnorm(
    n = n_judges,
    mu = rep(0, ndim_true),
    Sigma = diag(c(0.15, 0.15, 0.15))
  )
  judge_coord_offsets <- as.matrix(judge_coord_offsets)
  rownames(judge_coord_offsets) <- judge_ids
  colnames(judge_coord_offsets) <- dim_names
}

# Add judge-specific offsets to each judge's unique targets.
for (j in judge_ids) {
  ids_j <- paste0(j, "_T", 1:n_unique)
  coords_latent[ids_j, ] <- sweep(
    coords_latent[ids_j, , drop = FALSE],
    2,
    judge_coord_offsets[j, ],
    "+"
  )
}

# Standardize and rescale the latent coordinates for better spread.
coords_latent <- scale(coords_latent) * coord_spread

# Store the true coordinates in a data frame.
coords_true <- as.data.frame(coords_latent)
colnames(coords_true) <- dim_names
coords_true$id <- rownames(traits)
coords_true <- coords_true[, c("id", dim_names)]
rownames(coords_true) <- coords_true$id


## ------------------------------------------------------------
## 6. Target sets for each judge
## ------------------------------------------------------------

# Build the set of targets seen by each judge: common anchors + judge-specific targets.
judge_targets <- lapply(judge_ids, function(j) {
  c(anchor_ids, paste0(j, "_T", 1:n_unique))
})
names(judge_targets) <- judge_ids


## ------------------------------------------------------------
## 7. Simulate judge-specific dissimilarities
## ------------------------------------------------------------

# Simulate perceived dissimilarities for each judge.
judge_data <- lapply(judge_ids, function(j) {
  simulate_judge_dissimilarities(
    target_ids = judge_targets[[j]],
    coords_true = coords_true,
    noise_sd = distance_noise_sd
  )
})
names(judge_data) <- judge_ids


## ------------------------------------------------------------
## 8. Direct individual 3D MDS
## ------------------------------------------------------------

# Fit one 3D ordinal MDS solution for each judge.
mds_fits <- lapply(judge_ids, function(j) {
  fit <- fit_mds_3d(judge_data[[j]]$delta_ord)
  
  conf <- as.data.frame(fit$conf)
  colnames(conf) <- dim_names
  conf$id <- rownames(judge_data[[j]]$delta_ord)
  conf <- conf[, c("id", dim_names)]
  rownames(conf) <- conf$id
  
  list(
    fit = fit,
    coords = conf
  )
})
names(mds_fits) <- judge_ids

# Print the stress values for each judge's MDS solution.
cat("\nMDS stress for each judge:\n")
print(data.frame(
  judge = judge_ids,
  stress = sapply(mds_fits, function(x) x$fit$stress)
))


## ------------------------------------------------------------
## 9. GPA / Procrustes alignment on anchors
## ------------------------------------------------------------

# Extract judge-specific coordinate matrices for GPA.
configs_for_gpa <- lapply(mds_fits, function(x) {
  m <- as.matrix(x$coords[, dim_names, drop = FALSE])
  rownames(m) <- x$coords$id
  colnames(m) <- dim_names
  m
})

# Align all configurations through anchors.
gpa_fit <- gpa_from_anchors(
  configs = configs_for_gpa,
  anchor_ids = anchor_ids,
  max_iter = 1000,
  tol = 1e-9
)

# Store the consensus anchor configuration.
anchor_common <- as.data.frame(gpa_fit$consensus_anchors)
colnames(anchor_common) <- dim_names
anchor_common$id <- anchor_ids
anchor_common$judge <- "ALL"
anchor_common$type <- "Anchor"

# Store aligned judge-specific unique targets.
unique_common <- dplyr::bind_rows(lapply(seq_along(judge_ids), function(i) {
  M <- gpa_fit$aligned[[i]]
  ids_keep <- setdiff(rownames(M), anchor_ids)
  
  out <- as.data.frame(M[ids_keep, , drop = FALSE])
  colnames(out) <- dim_names
  out$id <- ids_keep
  out$judge <- judge_ids[i]
  out$type <- "Target"
  out
}))

# Combine anchors and unique targets into one common-space table.
plot_coords <- dplyr::bind_rows(anchor_common, unique_common)
plot_coords[, dim_names] <- scale(plot_coords[, dim_names]) * plot_spread

# Refresh anchor-only coordinates after rescaling.
anchor_common <- plot_coords %>%
  dplyr::filter(type == "Anchor")

# Print GPA diagnostics.
cat("\nGPA iterations:\n")
print(gpa_fit$iterations)

cat("\nFinal GPA delta:\n")
print(gpa_fit$delta_final)

cat("\nGPA converged:\n")
print(gpa_fit$converged)


## ------------------------------------------------------------
## 10A. Simulate personality ratings
## ------------------------------------------------------------

# Simulate judge-specific noisy ratings for unique targets.
unique_rating_list <- lapply(judge_ids, function(j) {
  ids_j <- paste0(j, "_T", 1:n_unique)
  
  noisy <- as.matrix(traits[ids_j, trait_names, drop = FALSE]) +
    matrix(rnorm(length(ids_j) * length(trait_names), 0, rating_noise_sd),
           nrow = length(ids_j), ncol = length(trait_names))
  
  noisy <- pmin(pmax(round(noisy), 1), 7)
  noisy <- as.data.frame(noisy)
  colnames(noisy) <- trait_names
  noisy$id <- ids_j
  noisy$judge <- j
  noisy$type <- "Target"
  noisy
})

# Bind judge-specific unique target ratings together.
unique_ratings <- dplyr::bind_rows(unique_rating_list) %>%
  dplyr::select(id, judge, type, dplyr::all_of(trait_names))

# Simulate judge-specific noisy ratings for anchors.
anchor_rating_list <- lapply(judge_ids, function(j) {
  noisy <- as.matrix(traits[anchor_ids, trait_names, drop = FALSE]) +
    matrix(rnorm(length(anchor_ids) * length(trait_names), 0, rating_noise_sd),
           nrow = length(anchor_ids), ncol = length(trait_names))
  
  noisy <- pmin(pmax(round(noisy), 1), 7)
  noisy <- as.data.frame(noisy)
  colnames(noisy) <- trait_names
  noisy$id <- anchor_ids
  noisy$judge <- j
  noisy$type <- "Anchor"
  noisy
})

# Bind anchor ratings from all judges.
anchor_ratings_long <- dplyr::bind_rows(anchor_rating_list) %>%
  dplyr::select(id, judge, type, dplyr::all_of(trait_names))

# Average anchor ratings across judges.
anchor_ratings_mean <- anchor_ratings_long %>%
  dplyr::group_by(id, type) %>%
  dplyr::summarise(
    dplyr::across(dplyr::all_of(trait_names), mean),
    .groups = "drop"
  ) %>%
  dplyr::mutate(judge = "ALL") %>%
  dplyr::select(id, judge, type, dplyr::all_of(trait_names))

# Build the final target score table used for summaries if needed.
target_scale_scores <- dplyr::bind_rows(
  unique_ratings,
  anchor_ratings_mean
) %>%
  dplyr::arrange(id)


## ------------------------------------------------------------
## 10B. Post-hoc orientation of MDS/GPA axes
##      using the true simulated traits
## ------------------------------------------------------------

# Join coordinates with trait data for axis orientation.
orient_df <- plot_coords %>%
  dplyr::left_join(
    traits %>% dplyr::select(id, dplyr::all_of(trait_names)),
    by = "id"
  )

# Build composite scores that define the intended latent dimensions.
orient_df$score_D1 <- znum(orient_df$H) + znum(orient_df$A)
orient_df$score_D2 <- znum(orient_df$X) + znum(orient_df$O)
orient_df$score_D3 <- znum(orient_df$E) + znum(orient_df$C)

# Enumerate all axis permutations.
perms <- all_axis_perms()

# Score each permutation by absolute axis-composite correlations.
perm_scores <- apply(perms, 1, function(p) {
  sum(abs(c(
    cor(orient_df[[dim_names[p[1]]]], orient_df$score_D1, use = "pairwise.complete.obs"),
    cor(orient_df[[dim_names[p[2]]]], orient_df$score_D2, use = "pairwise.complete.obs"),
    cor(orient_df[[dim_names[p[3]]]], orient_df$score_D3, use = "pairwise.complete.obs")
  )))
})

# Select the best axis permutation.
best_perm <- perms[which.max(perm_scores), ]

# Apply the best axis permutation to the coordinates.
tmp_coords <- as.matrix(plot_coords[, dim_names, drop = FALSE])
tmp_coords <- tmp_coords[, best_perm, drop = FALSE]
colnames(tmp_coords) <- dim_names
plot_coords[, dim_names] <- tmp_coords

# Rebuild the orientation data after permutation.
orient_df <- plot_coords %>%
  dplyr::left_join(
    traits %>% dplyr::select(id, dplyr::all_of(trait_names)),
    by = "id"
  )

# Recompute composite scores after permutation.
orient_df$score_D1 <- znum(orient_df$H) + znum(orient_df$A)
orient_df$score_D2 <- znum(orient_df$X) + znum(orient_df$O)
orient_df$score_D3 <- znum(orient_df$E) + znum(orient_df$C)

# Initialize the plotting beta matrix from the true betas.
beta_plot <- beta_true

# Flip axes when needed so that each dimension has the intended sign.
for (d in dim_names) {
  r_now <- cor(
    orient_df[[d]],
    orient_df[[paste0("score_", d)]],
    use = "pairwise.complete.obs"
  )
  
  if (is.na(r_now)) r_now <- 1
  
  if (r_now < 0) {
    plot_coords[[d]] <- -plot_coords[[d]]
    beta_plot[d, ] <- -beta_plot[d, ]
  }
}

# Refresh anchor coordinates after final orientation.
anchor_common <- plot_coords %>%
  dplyr::filter(type == "Anchor")

# Print axis orientation diagnostics.
cat("\nBest axis permutation applied:\n")
print(best_perm)

cat("\nFinal axis-composite correlations:\n")
print(c(
  cor_D1 = cor(plot_coords$D1, orient_df$score_D1, use = "pairwise.complete.obs"),
  cor_D2 = cor(plot_coords$D2, orient_df$score_D2, use = "pairwise.complete.obs"),
  cor_D3 = cor(plot_coords$D3, orient_df$score_D3, use = "pairwise.complete.obs")
))


## ------------------------------------------------------------
## 11. ICC(3,k) on anchors
## ------------------------------------------------------------

# Compute ICC(3,k) across judges for each trait on the anchors.
icc_results <- dplyr::bind_rows(lapply(trait_names, function(tr) {
  wide <- anchor_ratings_long %>%
    dplyr::select(id, judge, dplyr::all_of(tr)) %>%
    tidyr::pivot_wider(names_from = judge, values_from = dplyr::all_of(tr)) %>%
    dplyr::arrange(id)
  
  M <- as.matrix(wide[, -1, drop = FALSE])
  icc_obj <- psych::ICC(M)
  
  data.frame(
    trait = tr,
    ICC3k = icc_obj$results["ICC3k", "ICC"],
    lower = icc_obj$results["ICC3k", "lower bound"],
    upper = icc_obj$results["ICC3k", "upper bound"]
  )
}))

# Print ICC results.
cat("\nICC(3,k) by trait:\n")
print(icc_results)


## ------------------------------------------------------------
## 12. Inter-rater agreement on anchors
## ------------------------------------------------------------

# Build one flattened anchor-rating vector for each judge.
judge_vectors <- sapply(judge_ids, function(j) {
  as.numeric(
    unlist(
      anchor_ratings_long %>%
        dplyr::filter(judge == j) %>%
        dplyr::arrange(id) %>%
        dplyr::select(dplyr::all_of(trait_names))
    )
  )
})

# Name the judge columns.
colnames(judge_vectors) <- judge_ids

# Compute leave-one-out agreement correlations.
loo_corr <- sapply(seq_len(ncol(judge_vectors)), function(i) {
  cor(
    judge_vectors[, i],
    rowMeans(judge_vectors[, -i, drop = FALSE])
  )
})

# Compute pairwise judge correlations.
pair_corr <- combn(seq_len(ncol(judge_vectors)), 2, function(idx) {
  cor(judge_vectors[, idx[1]], judge_vectors[, idx[2]])
})

# Summarize agreement metrics.
agreement_summary <- data.frame(
  metric = c("mean_leave_one_out_correlation", "mean_pairwise_correlation"),
  value  = c(mean(loo_corr), mean(pair_corr))
)

# Print inter-rater agreement summary.
cat("\nInter-rater agreement summary:\n")
print(agreement_summary)


## ------------------------------------------------------------
## 13. Soft theoretical vectors for the plot
##     (more realistic and more visually distinct)
## ------------------------------------------------------------

# Start from the oriented beta matrix.
beta_plot_soft <- beta_plot

# Define the intended dominant axis for each trait.
dominant_axis <- c(
  H = "D1",
  A = "D1",
  X = "D2",
  O = "D2",
  E = "D3",
  C = "D3"
)

# Set a separate seed for vector perturbations.
set.seed(321)

# Guided offsets used to separate trait vectors visually.
# Rows = D1, D2, D3 ; columns = traits.
guided_offsets <- matrix(0, nrow = 3, ncol = length(trait_names))
rownames(guided_offsets) <- dim_names
colnames(guided_offsets) <- trait_names

guided_offsets[, "H"] <- c( 0.00,  0.32,  0.18)
guided_offsets[, "A"] <- c( 0.00, -0.28,  0.22)

guided_offsets[, "X"] <- c( 0.16,  0.00, -0.18)
guided_offsets[, "O"] <- c(-0.14,  0.00,  0.20)

guided_offsets[, "E"] <- c( 0.18, -0.12,  0.00)
guided_offsets[, "C"] <- c(-0.16,  0.14,  0.00)

# Perturb each vector while preserving its dominant axis.
for (tr in trait_names) {
  v <- beta_plot_soft[, tr]
  main_ax <- dominant_axis[tr]
  other_ax <- setdiff(dim_names, main_ax)
  
  # Slightly perturb the main loading while keeping it dominant.
  v[main_ax] <- v[main_ax] * runif(1, 0.82, 1.10)
  
  # Add guided offsets and random noise to the secondary axes.
  v <- v + guided_offsets[, tr]
  v[other_ax] <- v[other_ax] + rnorm(length(other_ax), mean = 0, sd = 0.12)
  
  # Enforce dominance of the main axis.
  max_other <- max(abs(v[other_ax]))
  if (abs(v[main_ax]) <= max_other + 0.20) {
    v[main_ax] <- sign(beta_plot[main_ax, tr]) * (max_other + runif(1, 0.28, 0.45))
  }
  
  # Preserve the intended sign on the dominant axis.
  v[main_ax] <- sign(beta_plot[main_ax, tr]) * abs(v[main_ax])
  
  # Save the softened vector.
  beta_plot_soft[, tr] <- v
}

# Store the softened betas in long form for plotting.
trait_betas <- data.frame(
  trait = trait_names,
  beta_D1 = beta_plot_soft["D1", trait_names],
  beta_D2 = beta_plot_soft["D2", trait_names],
  beta_D3 = beta_plot_soft["D3", trait_names],
  row.names = NULL
)

# Compute vector lengths.
vec_length <- apply(
  trait_betas[, c("beta_D1", "beta_D2", "beta_D3")],
  1,
  function(x) sqrt(sum(x^2, na.rm = TRUE))
)

# Compute unit-length vector directions.
vec_dir <- t(apply(
  trait_betas[, c("beta_D1", "beta_D2", "beta_D3")],
  1,
  function(x) {
    x[is.na(x)] <- 0
    nrm <- sqrt(sum(x^2))
    if (!is.finite(nrm) || nrm == 0) return(c(0, 0, 0))
    x / nrm
  }
))

# Name vector direction columns.
colnames(vec_dir) <- c("D1", "D2", "D3")

# Apply small trait-specific length weights for visual balance.
length_weight <- c(
  H = 1.00,
  E = 0.92,
  X = 1.02,
  A = 0.96,
  C = 0.90,
  O = 1.04
)

# Estimate the average coordinate range across axes.
axis_ranges <- sapply(dim_names, function(d) diff(range(plot_coords[[d]])))
mean_range <- mean(axis_ranges)

# Target fraction of the plot space occupied by the longest vector.
target_fraction <- 0.42

# Add a base length so that shorter arrows do not disappear.
base_length <- 1.10
vec_length_plot <- base_length + vec_length

# Compute the common scaling factor for plotting.
base_scale <- (mean_range * target_fraction) / max(vec_length_plot)

# Build the final vector endpoint table.
vec_df <- data.frame(
  trait = trait_betas$trait,
  x0 = 0, y0 = 0, z0 = 0,
  x1 = vec_dir[, "D1"] * vec_length_plot * base_scale * length_weight[trait_betas$trait],
  y1 = vec_dir[, "D2"] * vec_length_plot * base_scale * length_weight[trait_betas$trait],
  z1 = vec_dir[, "D3"] * vec_length_plot * base_scale * length_weight[trait_betas$trait]
)

# Build a readable summary of the plotted vectors.
label_summary <- data.frame(
  trait = trait_betas$trait,
  beta_D1 = round(trait_betas$beta_D1, 3),
  beta_D2 = round(trait_betas$beta_D2, 3),
  beta_D3 = round(trait_betas$beta_D3, 3),
  vec_length = round(vec_length, 3),
  dominant_axis = dominant_axis[trait_betas$trait]
)

# Print the final soft vectors used in the plot.
cat("\nSoft theoretical vectors used in the plot:\n")
print(label_summary)


## ------------------------------------------------------------
## 14. Final coordinate table
## ------------------------------------------------------------

# Order the final coordinates by type, judge, and ID.
final_coords_table <- plot_coords %>%
  dplyr::arrange(type, judge, id)

# Print the final coordinates in the common space.
cat("\nFinal coordinates in the common space:\n")
print(final_coords_table)


## ------------------------------------------------------------
## 15. Final 3D plot
## ------------------------------------------------------------

# Define point colors for judges.
judge_colors <- c(
  J1 = "#0072B2",
  J2 = "#E69F00",
  J3 = "#009E73"
)

# Define anchor color.
anchor_color <- "#E41A1C"

# Define colors for trait vectors.
vector_colors <- c(
  H = "#FF7F00",
  E = "#4DAF4A",
  X = "#E41A1C",
  A = "#984EA3",
  C = "#A65628",
  O = "#F781BF"
)

# Initialize the Plotly figure.
fig <- plot_ly()

# Add target points for each judge.
for (j in judge_ids) {
  tmp <- plot_coords %>%
    dplyr::filter(type == "Target", judge == j)
  
  fig <- fig %>%
    add_trace(
      data = tmp,
      x = ~D1,
      y = ~D2,
      z = ~D3,
      type = "scatter3d",
      mode = "markers+text",
      text = ~id,
      textposition = "top center",
      hovertext = ~paste0("Target: ", id, "<br>Judge: ", judge),
      hoverinfo = "text",
      marker = list(
        size = 5,
        color = judge_colors[[j]],
        symbol = "circle",
        opacity = 1,
        line = list(color = "black", width = 0.5)
      ),
      name = j,
      legendgroup = j,
      showlegend = TRUE,
      inherit = FALSE
    )
}

# Prepare anchor labels with profile names.
anchor_plot <- anchor_common %>%
  dplyr::left_join(anchor_profiles[, c("id", "profile")], by = "id") %>%
  dplyr::mutate(label_text = paste0(id, ": ", profile))

# Add anchor points.
fig <- fig %>%
  add_trace(
    data = anchor_plot,
    x = ~D1,
    y = ~D2,
    z = ~D3,
    type = "scatter3d",
    mode = "markers+text",
    text = ~id,
    textposition = "top center",
    hovertext = ~label_text,
    hoverinfo = "text",
    marker = list(
      size = 10,
      color = anchor_color,
      symbol = "diamond",
      opacity = 1,
      line = list(color = "black", width = 1)
    ),
    name = "Anchors",
    legendgroup = "Anchors",
    showlegend = TRUE,
    inherit = FALSE
  )

# Add one trait vector at a time.
for (i in seq_len(nrow(vec_df))) {
  v <- vec_df[i, ]
  vcol <- vector_colors[[v$trait]]
  
  fig <- fig %>%
    add_trace(
      x = c(v$x0, v$x1),
      y = c(v$y0, v$y1),
      z = c(v$z0, v$z1),
      type = "scatter3d",
      mode = "lines+markers",
      line = list(
        width = 6,
        color = vcol
      ),
      marker = list(
        size = c(1, 4),
        color = c(vcol, vcol),
        opacity = 1
      ),
      name = v$trait,
      legendgroup = "Traits",
      showlegend = TRUE,
      inherit = FALSE
    )
}

# Finalize layout settings.
fig <- fig %>%
  layout(
    title = paste0(
      "Study 3 simulation: common psychogeometric space",
      " (anchors = ", n_anchors,
      ", targets/judge = ", n_targets_per_judge, ")"
    ),
    scene = list(
      xaxis = list(title = "D1"),
      yaxis = list(title = "D2"),
      zaxis = list(title = "D3")
    ),
    legend = list(
      title = list(text = "Groups and traits"),
      orientation = "v",
      x = 1.02,
      y = 1
    ),
    margin = list(t = 90, b = 20, l = 20, r = 20)
  )

# Display the final interactive plot.
fig