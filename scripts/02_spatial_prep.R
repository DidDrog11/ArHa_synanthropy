# ==============================================================================
# Script: 02_extract_covars.R
# Project: ArHa Synanthropy
# Purpose: Pool all ArHa and GBIF points, extract 1km Global Human Modification,
#          and spatially intersect with data-driven MCPs to define background pools.
# ==============================================================================

source("scripts/00_setup.R")

# 1. Load the Global Raster ----------------------------------------------------
raster_path <- here("data", "spatial_covariates", "gHM_1km_global.tif")
r_modify <- rast(raster_path)

# 2. Load Raw Data (Temporally Filtered & Rodents Only) ------------------------
# We need the species names to build the individual polygons later
gbif_raw <- read_rds(here("data", "processed", "gbif_target_group_background.rds")) |> 
  drop_na(year) |>
  filter(year >= 2000) |> 
  inner_join(read_rds(here("data", "raw", "host_traits_final.rds")) |> 
               filter(order == "Rodentia") |> 
               mutate(species = str_replace_all(tip_label, "_", " ")) |> 
               select(species), by = "species") |> 
  select(species, longitude, latitude)

arha_db <- read_rds(here("data", "raw", "Project_ArHa_database_2026-10-23.rds"))
arha_host <- arha_db$host |> 
  drop_na(longitude, latitude) |> 
  mutate(trap_year = year(start_date)) |> 
  filter(trap_year >= 2000, host_order == "Rodentia") |> 
  select(host_species, longitude, latitude)

# 3. Master Coordinate Pool & Extraction ---------------------------------------

# Combine and get UNIQUE coordinates to minimize extraction time
all_coords_df <- bind_rows(gbif_raw |> select(longitude, latitude),
                           arha_host |> select(longitude, latitude)) |> 
  distinct(longitude, latitude)

all_coords_vect <- vect(all_coords_df, geom = c("longitude", "latitude"), crs = "EPSG:4326")

# Reproject points to match raster just for the extraction step
all_coords_proj <- project(all_coords_vect, crs(r_modify))
extracted_vals <- terra::extract(r_modify, all_coords_proj)

# Bind the extracted values back to the WGS84 vector
all_coords_vect$gHM <- extracted_vals[, 2]

# Drop points that fell in the ocean or outside raster bounds (NAs)
all_coords_clean <- all_coords_vect[!is.na(all_coords_vect$gHM), ]

# 4. Spatially Intersect with Data-Driven MCPs ---------------------------------
# We want backgrounds for all rodent species in our traits database
target_species <- read_rds(here("data", "raw", "host_traits_final.rds")) |> 
  filter(order == "Rodentia") |> 
  pull(tip_label) |> 
  str_replace_all("_", " ") |> 
  unique()

extract_spp_background <- function(spp_name) {
  
  # A. Gather coords for this specific species
  spp_coords <- bind_rows(arha_host |> filter(host_species == spp_name) |> select(longitude, latitude),
                          gbif_raw |> filter(species == spp_name) |> select(longitude, latitude)) |> 
    distinct()
  
  # Cannot draw a polygon with fewer than 3 points
  if(nrow(spp_coords) < 3) return(NULL) 
  
  # B. Generate terra MCP and Buffer
  spp_vect <- vect(spp_coords, geom = c("longitude", "latitude"), crs = "EPSG:4326")
  
  # terra handles WGS84 buffering geodesically (width is in meters)
  spp_hull <- convHull(spp_vect)
  spp_buffer <- buffer(spp_hull, width = 50000) # 50km buffer
  
  # C. Target-Group Intersection
  # Which points from the MASTER pool fall inside THIS species' buffer?
  intersected <- terra::intersect(all_coords_clean, spp_buffer)
  
  if(nrow(intersected) == 0) return(NULL)
  
  # D. Format output
  coords_df <- as_tibble(crds(intersected)) |> rename(longitude = x, latitude = y)
  
  result <- coords_df |> 
    mutate(species_range = spp_name,
           gHM = intersected$gHM) |> 
    select(species_range, longitude, latitude, gHM)
  
  return(result)
}

# Run the extraction for all species
points_to_ranges <- map_dfr(target_species, extract_spp_background, .progress = TRUE)

# 5. Save Model-Ready Spatial Database -----------------------------------------
write_rds(points_to_ranges, here("data", "processed", "model_ready_spatial_points.rds"))