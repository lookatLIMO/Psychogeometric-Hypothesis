# Psychogeometric-Hypothesis
This repository collects the code developed for my PhD application, including power simulations, modeling, and supporting scripts. It documents the methods, workflows, and computational tools used across the projects presented in my application.

Anchor-ENVFIT-sim_Commented.R
This script simulates a shared psychogeometric space across judges by generating latent personality profiles for common anchors and judge-specific targets, projecting them into a true 3D latent structure, and transforming continuous distances into noisy ordinal dissimilarity matrices.
It then fits separate 3D ordinal MDS solutions for each judge, aligns the recovered configurations through anchor-based GPA/Procrustes, and reorients the common space so that the final axes correspond to theoretically intended trait composites.
The script also evaluates anchor agreement and adds envfit-style trait vectors derived from the true beta structure, producing a visually interpretable representation of the latent psychogeometric space without fitting final mixed models

GPA-Stability-sim_Commented.R
This script performs a Monte Carlo simulation to assess how accurately a common latent geometric structure can be recovered across multiple judges when noisy dissimilarities are analyzed with ordinal MDS and aligned using anchor-based generalized Procrustes analysis.
Across design conditions varying the number of judges, unique targets, fitted dimensions, and anchors, it evaluates recovery quality using global distance correlations, coordinate RMSE, anchor RMSE, anchor spread, MDS stress, GPA convergence, and the proportion of stable simulations.
The final output classifies each condition as feasible, borderline, or unstable, making the script useful for identifying design settings in which anchor-based alignment is robust enough for empirical use.

Jackknife-3D-sim_Commented.R
This script simulates a symmetric ordinal dissimilarity matrix for 16 objects, fits a 3-dimensional ordinal MDS solution with SMACOF, and applies a leave-one-out jackknife procedure to evaluate the local stability of the recovered configuration.
It extracts the original MDS solution, the jackknife centroid configuration, and the leave-one-out estimates for each object, then displays them in an interactive 3D plot with centroids and connecting segments.
The result is a compact diagnostic tool for visually and numerically inspecting the robustness, dispersion, and replicability of a psychogeometric solution.

Pa-raw-similarities_Commented.R
This script conducts a simulation-based power analysis for linear mixed-effects models built on raw pairwise similarity data, using demographic predictors and HEXACO-based pairwise differences between targets.
It evaluates two inferential scenarios: an omnibus test of the incremental HEXACO block beyond demographic covariates and a single fixed-effect test within the full mixed model, while parallelizing simulations across chunks to reduce computation time.
For each effect-size condition, it aggregates simulated results, computes confidence intervals for estimated power, identifies the minimum effect needed to reach target power, and produces summary tables and plots.

PA-Anchor_Commented.R
This script estimates statistical power in the anchor design for detecting a focal pairwise-difference predictor in a linear mixed-effects model, where the outcome is a singular MDS dimension and the model includes random intercepts for the two members of each pair.
It repeatedly simulates synthetic pairwise datasets, fits the same mixed model under different combinations of target number and true focal effect size, and uses an adaptive stopping rule based on Monte Carlo standard error to balance precision and computation time.
In addition to power estimates, the script reports confidence intervals, bias, estimate variability, singularity and convergence rates, checkpointed results, and a final table of minimum detectable effects for the evaluated design conditions.
