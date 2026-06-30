# ==============================================================================
# Script: 00_setup.R
# Project: ArHa Synanthropy - Continuous Spatial Framework
# Purpose: Initialise repository structure, load dependencies, and configure 
#          spatial processing environments.
# ==============================================================================

# 1. Package Management --------------------------------------------------------
if (!require("pacman")) install.packages("pacman")

pacman::p_load(
  # Data wrangling and workflow
  tidyverse, here, fs, cli,
  
  # Spatial data handling
  terra, tidyterra, sf,
  
  # Phylogenetics
  ape,
  
  # Modelling and APIs
  spOccupancy, brms, rgbif, missRanger, brms, tidybayes, ape,
  
  # Parallel processing
  future, furrr
)

# 2. Directory Scaffolding -----------------------------------------------------
# Define the core project structure
dirs_to_create <- c(
  here("data", "raw"),                  # Hand-off from macroecology paper
  here("data", "spatial_covariates"),   # Heavy 1km rasters (e.g., GHM, GHSL)
  here("data", "processed"),            # Extracted data, GBIF pulls
  here("data", "temp_terra"),           # Custom temp directory for spatial processing
  here("outputs", "models"),            # Saved .rds model objects
  here("outputs", "figures"),           # Plots and maps
  here("outputs", "tables")             # Summary statistics
)

# Create directories if they don't exist
walk(dirs_to_create, ~dir.create(.x, recursive = TRUE, showWarnings = FALSE))

# 3. Configure .gitignore ------------------------------------------------------
# Prevent massive spatial files and model posteriors from breaking GitHub
gitignore_path <- here(".gitignore")

ignore_lines <- c(
  ".Rproj.user",
  ".Rhistory",
  ".RData",
  ".Ruserdata",
  "data/spatial_covariates/*",
  "data/temp_terra/*",
  "outputs/models/*",
  "*.rds",
  "*.tif",
  "*.shp",
  "*.shx",
  "*.dbf",
  "*.prj"
)

# Write to .gitignore if it doesn't exist, or append if needed
if (!file.exists(gitignore_path)) {
  write_lines(ignore_lines, gitignore_path)
}

# 4. Spatial Processing & Parallelization Configuration ------------------------

# A. Set custom temp directory for terra
# 1km global rasters will generate massive temporary files. 
terraOptions(tempdir = here("data", "temp_terra"))

# Set memory limits for terra
terraOptions(memfrac = 0.8) 

# B. Set up parallel processing plan
plan(multisession, workers = 4)

# Increase the maximum allowable size of global variables exported to workers
options(future.globals.maxSize = 8000 * 1024^2)

# 5. Define Global Constants ---------------------------------------------------
# Mollweide equal-area projection string
mollweide_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

# Target taxonomic orders
target_orders <- c("Rodentia")
