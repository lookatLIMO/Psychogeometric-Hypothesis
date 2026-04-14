# ============================================================
# JACKKNIFE SIMULATION ON A MDS SOLUTION
# ------------------------------------------------------------
# This script simulates a symmetric dissimilarity matrix for
# 16 objects, fits a 3-dimensional ordinal MDS solution using
# the SMACOF algorithm, and applies a leave-one-out jackknife
# procedure to evaluate the stability of the configuration.
#
# It then extracts the full MDS configuration, the jackknife
# centroid configuration, and the leave-one-out estimates for
# each object, and visualizes them in an interactive 3D plot.
#
# In the final plot, each object is represented by:
# - multiple leave-one-out jackknife estimates,
# - a centroid summarizing its jackknife positions,
# - and line segments connecting each estimate to its centroid.
#
# The code therefore provides both a graphical and numerical
# inspection of the robustness, replicability, and dispersion
# of the MDS solution.
# =========================================================

# --- Required packages ---
if(!require("pacman")) install.packages("pacman")
pacman::p_load("smacof", "vegan", "plotly", "colorspace", "htmlwidgets")

# Set the random seed to ensure reproducibility
set.seed(123)

# =========================
# 1. Simulate dissimilarities and run MDS 
# =========================

# Number of objects/targets
n <- 16

# Number of unique object pairs
n_pairs <- choose(n, 2)

# Randomly generate ordinal dissimilarities from 1 to 9
dissimilarities <- sample(1:9, n_pairs, replace = TRUE)

# Number of MDS dimensions
ndim = 3

# Build a symmetric dissimilarity matrix
diss_matrix <- matrix(0, n, n)

# Fill the lower triangle with the simulated dissimilarities
diss_matrix[lower.tri(diss_matrix)] <- dissimilarities

# Mirror the lower triangle to the upper triangle
diss_matrix <- diss_matrix + t(diss_matrix)

# Fit a 3-dimensional ordinal MDS solution using SMACOF
fit_matrix_final <- smacofSym(
  diss_matrix,          # Input symmetric dissimilarity matrix
  type = "ordinal",     # Ordinal MDS
  ndim = ndim,          # Number of dimensions
  ties = "secondary",   # Secondary treatment of ties
  init = "torgerson",   # Torgerson starting configuration
  principal = TRUE,     # Apply principal axis transformation
  itmax = 2000,         # Maximum number of iterations
  eps = 1e-7            # Convergence tolerance
)

# =========================
# 2. SMACOF jackknife
# =========================

# Compute leave-one-out jackknife estimates from the fitted MDS solution
jackMatrix <- jackmds(fit_matrix_final, itmax = 1000)

# Plot the default jackknife output
plot(jackMatrix, legend = TRUE, cex = 1.2, cex.legend = 0.8, inset = c(-0.2, 0.1))

# =========================
# 3. Extract configurations
# =========================

# Retrieve the dimensionality of the fitted configuration
ndim <- ncol(fit_matrix_final$conf)

# Retrieve the number of objects
n <- nrow(fit_matrix_final$conf)

# Full SMACOF configuration
X_ref <- jackMatrix$smacof.conf

# Jackknife centroid configuration
Xbar_star <- jackMatrix$comparison.conf

# =========================
# 4. Prepare data for the 3D plot
# =========================

# Define a color palette with one color per object
palette <- rainbow(n)

# Initialize an empty data frame for leave-one-out coordinates
plot_data <- data.frame()

# Initialize an empty data frame for segments connecting estimates to centroids
line_data <- data.frame()

# Loop over all objects
for (obj in seq_len(n)) {
  
  # Extract the object label
  obj_name <- rownames(X_ref)[obj]
  
  # Extract jackknife coordinates for the current object
  coords <- jackMatrix$jackknife.conf[obj, , ]
  
  # Transpose so that rows correspond to replications
  coords <- t(coords)
  
  # Convert to data frame
  coords <- as.data.frame(coords)
  
  # Assign dimension names
  colnames(coords) <- paste0("D", 1:ndim)
  
  # Add object label
  coords$Object <- obj_name
  
  # Append to the full plot data
  plot_data <- rbind(plot_data, coords)
  
  # Loop over all jackknife replications for the current object
  for (k in seq_len(nrow(coords))) {
    
    # Add a segment from each leave-one-out estimate to its jackknife centroid
    line_data <- rbind(line_data, data.frame(
      x = c(coords[k, 1], Xbar_star[obj, 1]),
      y = c(coords[k, 2], Xbar_star[obj, 2]),
      z = if (ndim >= 3) c(coords[k, 3], Xbar_star[obj, 3]) else 0,
      Object = obj_name
    ))
  }
}

# Convert centroid coordinates into a data frame
mean_df <- as.data.frame(Xbar_star)

# Assign dimension names to centroid data
colnames(mean_df) <- paste0("D",1:ndim)

# Add object labels to centroid data
mean_df$Object <- rownames(X_ref)

# =========================
# 5. Adaptive 3D plot
# =========================

# Initialize an empty Plotly object
p <- plot_ly()

# Loop over objects to add lines, points, and centroids
for (obj in seq_len(n)) {
  
  # Extract the object label
  obj_name <- rownames(X_ref)[obj]
  
  # Subset line segments for the current object
  lines_obj  <- subset(line_data, Object == obj_name)
  
  # Subset leave-one-out points for the current object
  points_obj <- subset(plot_data, Object == obj_name)
  
  # Subset jackknife centroid for the current object
  mean_obj   <- subset(mean_df, Object == obj_name)
  
  # Add line segments from leave-one-out estimates to the centroid
  p <- add_trace(p, data=lines_obj, x=~x, y=~y, z=~z,
                 type="scatter3d", mode="lines",
                 line=list(width=1,color=palette[obj]),
                 name=obj_name, legendgroup=obj_name, showlegend=TRUE)
  
  # Add leave-one-out points
  p <- add_trace(p, data=points_obj, x=~D1, y=~D2, z=~D3,
                 type="scatter3d", mode="markers",
                 marker=list(size=3,color=palette[obj]),
                 legendgroup=obj_name, showlegend=FALSE)
  
  # Add the jackknife centroid point with label
  p <- add_trace(p, data=mean_obj, x=~D1, y=~D2, z=~D3,
                 type="scatter3d", mode="markers+text",
                 marker=list(size=6,color=palette[obj]),
                 text=obj_name, textposition="top center",
                 legendgroup=obj_name, showlegend=FALSE)
}

# Customize the 3D layout
p <- p %>% layout(scene=list(xaxis=list(title="D1"),
                             yaxis=list(title="D2"),
                             zaxis=list(title="D3"),
                             aspectmode="data"),
                  legend=list(x=1.05,y=0.95))

# Display the interactive 3D plot
p

# Print the full jackknife object
jackMatrix