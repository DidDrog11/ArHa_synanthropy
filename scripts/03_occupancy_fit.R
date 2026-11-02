# ==============================================================================
# Script: 03_occupancy_fit.R
# Project: ArHa Synanthropy
# Purpose: Fit integrated Bayesian occupancy models (ArHa + GBIF) using 
#          phylogenetically stratified, temporally aligned, session-based data.
# ==============================================================================

source("scripts/00_setup.R")

# 1. Load Data & Filter --------------------------------------------------------
# 1A. Environmental Points (From Data-Driven MCPs in Script 02)
env_points <- read_rds(here("data", "processed", "model_ready_spatial_points.rds")) |> 
  mutate(longitude = round(longitude, 5), latitude = round(latitude, 5))

# 1B. GBIF Raw Data (Temporally Aligned & Rodents Only)
gbif_raw <- read_rds(here("data", "processed", "gbif_target_group_background.rds")) |> 
  drop_na(year) |>
  filter(year >= 2000) |> 
  inner_join(read_rds(here("data", "raw", "host_traits_final.rds")) |> 
               filter(order == "Rodentia") |> 
               mutate(species = str_replace_all(tip_label, "_", " ")) |> 
               select(species), by = "species") |> 
  select(species, gbif_id = query_gbif_id, year, longitude, latitude) |> 
  # Apply rounding
  mutate(longitude = round(longitude, 5), latitude = round(latitude, 5))

# 1C. ArHa Structured Database (Temporally Aligned & Rodents Only)
arha_db <- read_rds(here("data", "raw", "Project_ArHa_database_2026-10-23.rds"))
arha_host <- arha_db$host |> 
  drop_na(longitude, latitude) |>
  filter(coordinate_resolution_processed %in% c("site", "village", "town", "city", "study_area")) |>
  mutate(trap_year = year(start_date)) |> 
  filter(trap_year >= 2000) |>
  filter(host_order == "Rodentia") |> 
  # Apply rounding
  mutate(longitude = round(longitude, 5), latitude = round(latitude, 5))

# --- Define Target Species ---
# ArHa Core: Require >= 15 sampling events (unique site + session)
arha_counts <- arha_host |> 
  filter(number_of_hosts > 0) |> 
  drop_na(host_species) |>
  group_by(host_species) |> 
  summarise(n_arha = n_distinct(longitude, latitude, start_date)) |> 
  filter(n_arha >= 15)

arha_core <- arha_counts$host_species

# GBIF Infill: Phylogenetically Stratified
rodent_tips <- read_rds(here("data", "raw", "host_traits_final.rds")) |> 
  filter(order == "Rodentia") |> 
  pull(tip_label)
rodent_tree <- read.tree(here("data", "raw", "mammal_tree_matched.tre"))
valid_rodent_tips <- intersect(rodent_tree$tip.label, rodent_tips)
rodent_tree <- keep.tip(rodent_tree, valid_rodent_tips)

dist_matrix <- cophenetic.phylo(rodent_tree)
tree_clusters <- hclust(as.dist(dist_matrix))
clades <- cutree(tree_clusters, k = 150)

clade_df <- tibble(species = str_replace_all(names(clades), "_", " "), 
                   clade_id = unname(clades))

gbif_counts <- gbif_raw |> 
  filter(!species %in% arha_core) |> 
  group_by(species) |> 
  summarise(n_gbif = n_distinct(longitude, latitude, year)) |> 
  filter(n_gbif >= 50)

gbif_stratified <- gbif_counts |> 
  inner_join(clade_df, by = "species") |> 
  group_by(clade_id) |> 
  slice_max(order_by = n_gbif, n = 1, with_ties = FALSE) |> 
  ungroup() |> 
  pull(species)

target_species <- unique(c(arha_core, gbif_stratified))

cli_alert_success(paste("Total target species for occupancy modelling:", length(target_species)))

# Create Output Directory
model_out_dir <- here("outputs", "models", "occupancy")
dir.create(model_out_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Pre-Compute ArHa Site Effort ----------------------------------------------
arha_sites <- arha_host |> 
  group_by(longitude, latitude, start_date) |> 
  summarise(total_trapped = sum(number_of_hosts, na.rm = TRUE), .groups = "drop") |> 
  mutate(effort = log1p(total_trapped)) 

# 3. Define the Modelling Function ---------------------------------------------
fit_species_occupancy <- function(spp_name) {
  
  safe_name <- str_replace_all(spp_name, " ", "_")
  outfile <- file.path(model_out_dir, paste0(safe_name, "_occ.rds"))
  
  if (file.exists(outfile)) {
    cli_alert_success(paste("  ->", spp_name, "already complete. Skipping."))
    return(NULL)
  }
  
  # --- A. Environmental Filter ---
  spp_env <- env_points |> 
    filter(species_range == spp_name) |> 
    distinct(longitude, latitude, .keep_all = TRUE)
  
  if(nrow(spp_env) == 0) return(NULL)
  
  # --- B. ArHa Structured Data ---
  spp_arha_counts <- arha_host |> 
    filter(host_species == spp_name) |> 
    group_by(longitude, latitude, start_date) |> 
    summarise(spp_trapped = sum(number_of_hosts, na.rm = TRUE), .groups = "drop")
  
  arha_model_data <- arha_sites |> 
    left_join(spp_arha_counts, by = c("longitude", "latitude", "start_date")) |> 
    mutate(spp_trapped = replace_na(spp_trapped, 0),
           y_arha = ifelse(spp_trapped > 0, 1, 0)) |> 
    inner_join(spp_env, by = c("longitude", "latitude")) |> 
    distinct(longitude, latitude, start_date, .keep_all = TRUE)
  
  # --- C. GBIF Unstructured Data ---
  gbif_model_data <- gbif_raw |> 
    inner_join(spp_env, by = c("longitude", "latitude")) |> 
    mutate(y_gbif = ifelse(species == spp_name, 1, 0)) |> 
    distinct(longitude, latitude, year, .keep_all = TRUE)
  
  # --- D. Clean NAs and Check Sufficiency ---
  arha_clean <- arha_model_data |> drop_na(gHM)
  gbif_clean <- gbif_model_data |> drop_na(gHM)
  
  gbif_pres <- gbif_clean |> filter(y_gbif == 1)
  gbif_abs <- gbif_clean |> filter(y_gbif == 0)
  
  # Cap the background at 5,000 points, but keep ALL presences
  if(nrow(gbif_abs) > 5000) {
    gbif_abs <- gbif_abs |> slice_sample(n = 5000)
  }
  
  # Recombine
  gbif_clean <- bind_rows(gbif_pres, gbif_abs)
  
  n_arha_pres <- sum(arha_clean$y_arha, na.rm = TRUE)
  n_gbif_pres <- sum(gbif_clean$y_gbif, na.rm = TRUE)
  
  if (n_arha_pres < 5 && n_gbif_pres < 20) {
    write_rds("Insufficient Data", file.path(model_out_dir, paste0(safe_name, "_FAILED.rds")))
    return(NULL)
  }
  
  # --- E. Prepare Data List ---
  occ_covs <- bind_rows(select(arha_clean, gHM), select(gbif_clean, gHM)) |> 
    mutate(gHM_sq = gHM^2)
  
  n_arha <- nrow(arha_clean)
  n_gbif <- nrow(gbif_clean)
  
  data_list <- list(
    y = list(ArHa = arha_clean$y_arha, GBIF = gbif_clean$y_gbif),
    occ.covs = as.matrix(occ_covs),
    det.covs = list(
      ArHa = list(effort = arha_clean$effort),
      GBIF = list(int = rep(1, n_gbif))
    ),
    sites = list(ArHa = 1:n_arha, GBIF = (n_arha + 1):(n_arha + n_gbif)),
    coords = as.matrix(bind_rows(
      select(arha_clean, longitude, latitude),
      select(gbif_clean, longitude, latitude)
    ))
  )
  
  # --- F. Fit the Model ---
  occ_formula <- ~ gHM + gHM_sq
  det_formulas <- list(ArHa = ~ effort, GBIF = ~ 1)
  
  tryCatch({
    out <- suppressMessages(
      intPGOcc(
        occ.formula = occ_formula,
        det.formula = det_formulas,
        data = data_list,
        n.samples = 20000,   
        n.burn = 10000,      
        n.thin = 10,
        n.chains = 3,
        verbose = FALSE
      )
    )
    
    # STRIP THE MATRICES BEFORE SAVING
    out$z.samples <- NULL
    out$y.rep.samples <- NULL
    out$psi.samples <- NULL
    out$like.samples <- NULL
    out$p.samples <- NULL
    
    write_rds(out, outfile)
    
  }, error = function(e) {
    write_rds(e$message, file.path(model_out_dir, paste0(safe_name, "_FAILED.rds")))
  })
}

# 4. Execute Batch Fit ---------------------------------------------------------
cli_h1("Executing Production MCMC Model Fitting")

walk2(target_species, seq_along(target_species), function(spp_name, i) {
  cli_alert_info("[{i}/{length(target_species)}] Starting: {.emph {spp_name}}")
  fit_species_occupancy(spp_name)
})

cli_alert_success("All occupancy models complete!")
