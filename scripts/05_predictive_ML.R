# ==============================================================================
# Script: 05_pglmm_trait_model.R
# Project: ArHa Synanthropy
# Purpose: Train a Bayesian PGLMM using brms to predict continuous synanthropy 
#          across the rodent phylogeny, natively accounting for evolutionary history.
# ==============================================================================

source("scripts/00_setup.R")

# 1. Load Data -----------------------------------------------------------------
cli_h1("Loading Biological and Phylogenetic Data")

# Load Traits
traits_raw <- read_rds(here("data", "raw", "host_traits_final.rds")) |> 
  filter(order == "Rodentia") |> 
  mutate(species = str_replace_all(tip_label, "_", " ")) |> 
  select(species, tip_label, pace_of_life_pc1, repro_strategy_pc2, maternal_investment_pc3,
         geographic_range_km2, habitat_breadth, diet_breadth)

# Load Ground Truth Synanthropy (The 84 anchor species)
syn_targets <- read_rds(here("data", "processed", "BART_ready_synanthropy_traits.rds"))

# Load the Phylogenetic Tree
rodent_tree <- read.tree(here("data", "raw", "mammal_tree_matched.tre"))
valid_tips <- intersect(rodent_tree$tip.label, traits_raw$tip_label)
rodent_tree <- keep.tip(rodent_tree, valid_tips)

# Convert tree to a covariance matrix for the PGLMM
A_matrix <- vcv.phylo(rodent_tree)

# 2. Impute, Merge, and SCALE --------------------------------------------------
cli_h1("Imputing Predictors & Formatting Data")

# Impute missing ecological predictors
set.seed(42)
predictor_complete <- missRanger(traits_raw, formula = . ~ . - species - tip_label, 
                                 num.trees = 500, verbose = 0)

# Merge everything
master_data <- predictor_complete |> 
  left_join(syn_targets, by = "species") |> 
  select(-max_rhat, -status) |> 
  # THE FIX 1: Scale the massive raw numbers into Z-scores!
  mutate(
    range_scaled = scale(geographic_range_km2)[,1],
    habitat_scaled = scale(habitat_breadth)[,1],
    diet_scaled = scale(diet_breadth)[,1]
  )

# BETA REGRESSION SQUEEZE
master_data <- master_data |> 
  mutate(
    peak_gHM_squeezed = ifelse(
      !is.na(peak_gHM_median),
      (peak_gHM_median * (n() - 1) + 0.5) / n(), 
      NA
    )
  )

train_df <- master_data |> filter(!is.na(peak_gHM_squeezed))
shadow_df <- master_data |> filter(is.na(peak_gHM_squeezed))

cli_alert_info(paste("Training set:", nrow(train_df), "anchor species."))
cli_alert_info(paste("Shadow set:", nrow(shadow_df), "species to predict."))


# 3. Fit the Phylogenetic GLMM -------------------------------------------------
cli_h1("Fitting Bayesian PGLMM")

# The Formula: Using the new _scaled variables!
pglmm_formula <- bf(
  peak_gHM_squeezed ~ pace_of_life_pc1 + repro_strategy_pc2 + maternal_investment_pc3 + 
    range_scaled + habitat_scaled + diet_scaled + 
    (1 | gr(tip_label, cov = A_matrix))
)

# Run the model
pglmm_model <- brm(
  formula = pglmm_formula,
  data = train_df,
  data2 = list(A_matrix = A_matrix), 
  family = Beta(),                   
  chains = 4, cores = 4, iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.95), 
  init = "0",  # THE FIX 2: Forces Stan to start safely in the middle
  seed = 123,
  backend = "cmdstanr" 
)

print(summary(pglmm_model))

# 4. Predict the Shadow Species ------------------------------------------------
cli_h1("Predicting Global Synanthropy with Credible Intervals")

# Predict the shadow species. 
# allow_new_levels = TRUE lets brms use the phylogenetic covariance matrix 
# to estimate traits for species it didn't see in the training data!
predictions <- fitted(
  pglmm_model, 
  newdata = shadow_df, 
  data2 = list(A_matrix = A_matrix),
  allow_new_levels = TRUE, 
  sample_new_levels = "gaussian"
) |> as_tibble()

shadow_formatted <- shadow_df |> 
  select(species) |> 
  bind_cols(predictions) |> 
  select(
    species, 
    final_gHM_median = Estimate, 
    final_gHM_lower = Q2.5, 
    final_gHM_upper = Q97.5
  ) |> 
  mutate(data_source = "Predicted (PGLMM)")

train_formatted <- train_df |> 
  select(species, final_gHM_median = peak_gHM_median, 
         final_gHM_lower = peak_gHM_lower, final_gHM_upper = peak_gHM_upper) |> 
  mutate(data_source = "Empirical (Occupancy)")

final_synanthropy_db <- bind_rows(train_formatted, shadow_formatted)

# 5. Save Outputs --------------------------------------------------------------
write_rds(final_synanthropy_db, here("outputs", "models", "global_synanthropy_pglmm.rds"))
write_csv(final_synanthropy_db, here("outputs", "tables", "global_synanthropy_database.csv"))
write_rds(pglmm_model, here("outputs", "models", "pglmm_fitted_model.rds"))

cli_alert_success("Phase 5 complete. Robust Bayesian dataset saved!")