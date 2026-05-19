# ==============================================================================
# Script: 01_gbif_pull.R
# Project: ArHa Synanthropy
# Purpose: Fetch high-quality GBIF occurrence data for all target species to 
#          generate the Target-Group Background (TGB) pseudo-absences.
# ==============================================================================

source("scripts/00_setup.R") # Ensures libraries, plan(multisession), and paths are loaded

if (!require("CoordinateCleaner")) install.packages("CoordinateCleaner")
library(CoordinateCleaner)

# 1. Load Target Species -------------------------------------------------------
traits <- read_rds(here("data", "raw", "host_traits_final.rds"))

# Extract valid GBIF IDs for our target species
# Filter out any missing gbif_ids
target_taxa <- traits |> 
  drop_na(gbif_id) |> 
  distinct(species, gbif_id)

# 2. Setup Cache Directory -----------------------------------------------------
cache_dir <- here("data", "processed", "gbif_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

# This creates a single server-side query for all species simultaneously
gbif_req <- occ_download(
  pred_in("taxonKey", target_taxa$gbif_id),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  pred_in("basisOfRecord", c("PRESERVED_SPECIMEN", "HUMAN_OBSERVATION", 
                             "MACHINE_OBSERVATION", "MATERIAL_SAMPLE")),
  format = "SIMPLE_CSV",
  user = Sys.getenv("GBIF_USER"),
  pwd = Sys.getenv("GBIF_PWD"),
  email = Sys.getenv("GBIF_EMAIL")
)

# Wait for the download to finish preparing on the GBIF servers
occ_download_wait(gbif_req)

# 3. Download and Import -------------------------------------------------------
# Download the zip file to our processed directory
download_path <- occ_download_get(gbif_req, path = here("data", "processed"), overwrite = TRUE)

# Import the CSV
gbif_raw <- occ_download_import(download_path) |> as_tibble()

# 4. Clean and Randomly Subset -------------------------------------------------
gbif_clean <- gbif_raw |> 
  select(
    query_gbif_id = speciesKey,
    species, 
    decimalLongitude, 
    decimalLatitude, 
    coordinateUncertaintyInMeters, 
    year
  ) |> 
  drop_na(decimalLongitude, decimalLatitude, species) |> 
  filter(is.na(coordinateUncertaintyInMeters) | coordinateUncertaintyInMeters <= 10000)

# CoordinateCleaner (Flags capitals, centroids, zeros)
gbif_scrubbed <- gbif_clean |> 
  clean_coordinates(
    lon = "decimalLongitude", 
    lat = "decimalLatitude",
    species = "species",
    tests = c("capitals", "centroids", "equal", "gbif", "institutions", "zeros"),
    value = "clean"
  )

# TRUE RANDOM SAMPLING: Prevent Mus/Rattus from dominating the background
# We take a random sample of up to 5000 records per species. 
# If a species has fewer than 5000, it keeps all of them.
set.seed(123)
gbif_final <- gbif_scrubbed |> 
  group_by(species) |> 
  slice_sample(n = 5000, replace = FALSE) |> 
  ungroup() |> 
  rename(longitude = decimalLongitude, latitude = decimalLatitude)

# 5. Save and Cleanup ----------------------------------------------------------
write_rds(gbif_final, here("data", "processed", "gbif_target_group_background.rds"))

# Optionally remove the massive raw zip file to save disk space
file.remove(download_path)