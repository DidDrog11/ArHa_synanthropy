# ==============================================================================
# Script: 04_extract_syn.R
# Project: ArHa Synanthropy
# Purpose: Extract continuous synanthropy metrics (optimal gHM peak) 
#          from Bayesian occupancy posteriors, using strict Rhat filtering.
# ==============================================================================

source("scripts/00_setup.R")

# 1. Locate Fitted Models ------------------------------------------------------
model_dir <- here("outputs", "models", "occupancy")
model_files <- list.files(model_dir, pattern = "_occ\\.rds$", full.names = TRUE)

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
  if (is.null(mod) || is.character(mod)) {
    
    synanthropy_traits[[i]] <- tibble(
      species = spp_name, status = "Skipped/Insufficient Data", 
      max_rhat = NA, peak_gHM_median = NA, peak_gHM_lower = NA, peak_gHM_upper = NA
    )
    
  } else {
    
    # --- B. Diagnostic Check (Rhat Filter) ---
    occ_rhats <- mod$rhat$beta
    worst_rhat <- max(occ_rhats, na.rm = TRUE)
    
    if (is.na(worst_rhat) || worst_rhat >= 1.15) {
      
      synanthropy_traits[[i]] <- tibble(
        species = spp_name, status = "Failed Convergence", 
        max_rhat = worst_rhat, peak_gHM_median = NA, peak_gHM_lower = NA, peak_gHM_upper = NA
      )
      
    } else {
      
      # --- C. Trait Extraction (Only runs if Rhat < 1.15) ---
      beta_mcmc <- mod$beta.samples
      
      # Calculate vertex (-b/2a) for all 3000 samples
      peak_posteriors <- -beta_mcmc[, "gHM"] / (2 * beta_mcmc[, "gHM_sq"])
      
      # Clamp to biological reality (gHM only exists from 0 to 1)
      peak_posteriors_clamped <- pmax(0, pmin(1, peak_posteriors))
      
      synanthropy_traits[[i]] <- tibble(
        species = spp_name,
        status = "Pass",
        max_rhat = worst_rhat,
        peak_gHM_median = median(peak_posteriors_clamped, na.rm = TRUE),
        peak_gHM_lower = quantile(peak_posteriors_clamped, 0.025, na.rm = TRUE),
        peak_gHM_upper = quantile(peak_posteriors_clamped, 0.975, na.rm = TRUE)
      )
    }
  }
}

# 3. Combine and Review Results ------------------------------------------------
final_traits <- bind_rows(synanthropy_traits)

# Print diagnostic breakdown to the console
print(table(final_traits$status))

# Isolate the pristine dataset for Machine Learning
clean_bart_data <- final_traits |> 
  filter(status == "Pass") |> 
  arrange(desc(peak_gHM_median))

cli_alert_success(paste("Successfully extracted highly confident traits for", nrow(clean_bart_data), "species."))

# 4. Save Target Variables -----------------------------------------------------
write_rds(clean_bart_data, here("data", "processed", "BART_ready_synanthropy_traits.rds"))
