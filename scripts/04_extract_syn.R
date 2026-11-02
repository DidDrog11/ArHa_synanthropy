# ==============================================================================
# Script: 04_extract_syn.R
# Project: ArHa Synanthropy
# Purpose: Extract continuous synanthropy metrics from the Bayesian occupancy
#          posteriors, using Rhat filtering. Primary index: occupancy-weighted
#          mean gHM over each species' sampled gHM range (decision log
#          2026-11-02). Sensitivity: gHM at maximum occupancy, high-minus-low
#          contrast, and the old clamped vertex.
# ==============================================================================

source("scripts/00_setup.R")

# 1. Locate Fitted Models ------------------------------------------------------
model_dir <- here("outputs", "models", "occupancy")
model_files <- list.files(model_dir, pattern = "_occ\\.rds$", full.names = TRUE)

# Grid settings: summaries are evaluated between each species' own 2.5% and
# 97.5% gHM quantiles, so the quadratic is never extrapolated
grid_n <- 101
range_probs <- c(0.025, 0.975)

# Posterior summary: median and 95% interval, named <prefix>_median etc.
summarise_draws <- function(draws, prefix) {
  out <- tibble(median = median(draws, na.rm = TRUE),
                lower = quantile(draws, 0.025, na.rm = TRUE, names = FALSE),
                upper = quantile(draws, 0.975, na.rm = TRUE, names = FALSE))
  names(out) <- paste(prefix, names(out), sep = "_")
  out
}

# Create an empty list to store traits
synanthropy_traits <- list()

cli_alert_info(paste("Starting extraction for", length(model_files), "models..."))

# 2. Execute Extraction Loop ---------------------------------------------------
for (i in seq_along(model_files)) {
  file_path <- model_files[i]
  spp_safe <- basename(file_path) |> str_remove("_occ\\.rds$")
  spp_name <- str_replace_all(spp_safe, "_", " ")

  # Safely load the model
  mod <- tryCatch(read_rds(file_path), error = function(e) NULL)

  # --- A. Check for Missing/Insufficient Data ---
  # Metric columns are left out here; bind_rows() fills them with NA
  if (is.null(mod) || is.character(mod)) {
    synanthropy_traits[[i]] <- tibble(species = spp_name, status = "Skipped/Insufficient Data", max_rhat = NA)
    next
  }

  # --- B. Diagnostic Check (Rhat Filter) ---
  worst_rhat <- max(mod$rhat$beta, na.rm = TRUE)

  if (is.na(worst_rhat) || worst_rhat >= 1.15) {
    synanthropy_traits[[i]] <- tibble(species = spp_name, status = "Failed Convergence", max_rhat = worst_rhat)
    next
  }

  # --- C. Predicted Occupancy Curve per Posterior Draw ---
  beta_mcmc <- mod$beta.samples

  # gHM grid across this species' sampled range (gHM values kept in the design matrix)
  spp_ghm <- mod$X[, "gHM"]
  ghm_range <- quantile(spp_ghm, range_probs, names = FALSE)
  ghm_grid <- seq(ghm_range[1], ghm_range[2], length.out = grid_n)

  # psi_draws: one row per posterior draw, one column per grid point
  logit_psi <- beta_mcmc[, "(Intercept)"] +
    outer(beta_mcmc[, "gHM"], ghm_grid) +
    outer(beta_mcmc[, "gHM_sq"], ghm_grid^2)
  psi_draws <- plogis(logit_psi)

  # --- D. Primary Index: Occupancy-Weighted Mean gHM and Tolerance (SD) ---
  psi_total <- rowSums(psi_draws)
  wmean_draws <- as.vector(psi_draws %*% ghm_grid) / psi_total
  wsd_draws <- sqrt(rowSums(psi_draws * (outer(wmean_draws, ghm_grid, "-"))^2) / psi_total)

  # --- E. Sensitivity Summaries ---
  # gHM at the highest predicted occupancy on the grid
  max_draws <- ghm_grid[apply(psi_draws, 1, which.max)]

  # Occupancy at the top minus the bottom of the sampled range
  contrast_draws <- psi_draws[, grid_n] - psi_draws[, 1]

  # Old index: parabola vertex (-b/2a), clamped to 0-1
  vertex_draws <- pmax(0, pmin(1, -beta_mcmc[, "gHM"] / (2 * beta_mcmc[, "gHM_sq"])))

  # Curve shape: the vertex is a peak only when gHM_sq < 0; when gHM_sq > 0
  # it is the trough of a U-shaped response
  prob_concave <- mean(beta_mcmc[, "gHM_sq"] < 0)

  synanthropy_traits[[i]] <- bind_cols(
    tibble(species = spp_name, status = "Pass", max_rhat = worst_rhat, prob_concave = prob_concave,
           ghm_range_lower = ghm_range[1], ghm_range_upper = ghm_range[2]),
    summarise_draws(wmean_draws, "wmean_gHM"),
    summarise_draws(wsd_draws, "wsd_gHM"),
    summarise_draws(max_draws, "max_gHM"),
    summarise_draws(contrast_draws, "contrast_psi"),
    summarise_draws(vertex_draws, "peak_gHM")
  )
}

# 3. Combine and Review Results ------------------------------------------------
final_traits <- bind_rows(synanthropy_traits)

# Print diagnostic breakdown to the console
print(table(final_traits$status))

# Isolate the pristine dataset for Machine Learning
clean_bart_data <- final_traits |>
  filter(status == "Pass") |>
  arrange(desc(wmean_gHM_median))

cli_alert_success(paste("Successfully extracted highly confident traits for", nrow(clean_bart_data), "species."))

# 4. Save Target Variables -----------------------------------------------------
write_rds(clean_bart_data, here("data", "processed", "BART_ready_synanthropy_traits.rds"))
