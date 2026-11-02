# ==============================================================================
# Script: scoping_inat_counts.R
# Project: ArHa Synanthropy
# Purpose: Scope iNaturalist as external validation. For every rodent in the
#          trait table, count GBIF occurrences from 2000 onwards and the share
#          that are iNaturalist research-grade observations. Counts come from
#          the GBIF API (no records downloaded). Not part of the numbered
#          pipeline.
# ==============================================================================

library(tidyverse)
library(here)
library(jsonlite)

# 1. Species List --------------------------------------------------------------
# All rodents in the trait table, with whether they have an occupancy fit and
# their Ecke et al. category
traits <- read_rds(here("data", "raw", "host_traits_final.rds")) |>
  filter(order == "Rodentia", !is.na(gbif_id)) |>
  select(tip_label, gbif_id, synanthropy_status)

occ_files <- list.files(here("outputs", "models", "occupancy"))
modelled <- tibble(tip_label = str_remove(occ_files, "_(occ|FAILED)\\.rds$"),
                   occupancy_fit = str_detect(occ_files, "_occ\\.rds$"))

species <- traits |> left_join(modelled, by = "tip_label")

# 2. GBIF API Query ------------------------------------------------------------
# Same record filters as 01_gbif_pull.R, restricted to 2000 onwards as in 02/03.
# The iNaturalist research-grade dataset is checked by name before use.
inat_key <- "50c9509d-22c7-4a22-a47d-8c48425ef4a7"
inat_title <- fromJSON(paste0("https://api.gbif.org/v1/dataset/", inat_key))$title
stopifnot(str_detect(inat_title, "iNaturalist"))

base_query <- paste0("https://api.gbif.org/v1/occurrence/search?limit=0",
                     "&hasCoordinate=true&hasGeospatialIssue=false&year=2000,2026",
                     "&basisOfRecord=PRESERVED_SPECIMEN&basisOfRecord=HUMAN_OBSERVATION",
                     "&basisOfRecord=MACHINE_OBSERVATION&basisOfRecord=MATERIAL_SAMPLE")

# GBIF rate-limits rapid calls (HTTP 429): pause between calls, and on failure
# wait longer and retry, up to 5 attempts
gbif_count <- function(taxon_key, extra = "") {
  url <- paste0(base_query, "&taxonKey=", taxon_key, extra)
  for (attempt in 1:5) {
    Sys.sleep(0.3)
    result <- tryCatch(fromJSON(url)$count, error = function(e) NULL)
    if (!is.null(result)) return(result)
    Sys.sleep(10 * attempt)
  }
  stop("GBIF query failed after 5 attempts: ", url)
}

# One call for all records and one for iNaturalist only, per species
species$n_gbif <- NA_integer_
species$n_inat <- NA_integer_
for (i in seq_len(nrow(species))) {
  species$n_gbif[i] <- gbif_count(species$gbif_id[i])
  species$n_inat[i] <- gbif_count(species$gbif_id[i], paste0("&datasetKey=", inat_key))
  if (i %% 100 == 0) message(i, " / ", nrow(species))
}

# 3. Save ----------------------------------------------------------------------
species <- species |> mutate(inat_share = ifelse(n_gbif > 0, n_inat / n_gbif, NA))
write_csv(species, here("data", "processed", "inat_gbif_counts.csv"))
