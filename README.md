# ArHa Synanthropy: A Continuous Measure of Rodent Synanthropy

## Project Overview
This project builds on the ArHa macroecology dataset and on the synanthropy classification of Ecke et al. (2022), which classified rodent species as synanthropic or not and showed that synanthropy helps explain the risk of rodent-borne zoonoses. We extend that work in two directions: we estimate synanthropy as a continuous response to human modification of the landscape, and we predict it for species that have not yet been classified. The Ecke et al. categories are the reference against which the continuous index is checked.

By integrating global surveillance records with Target-Group Background (TGB) data from GBIF, we use multi-source Bayesian occupancy models to quantify species-specific tolerance to anthropogenic modification. We then use imputed life-history traits and phylogeny to predict a continuous synanthropy index for over 2,000 unclassified rodent species. The index is framed ecologically first: what it shows about how rodent communities change with human modification, and about pest activity. Zoonotic hazard is discussed as an implication, not mapped.

## Analysis stance
The analysis is not formally preregistered. To avoid HARKing, analytical choices are recorded in a dated decision log before the results they affect are seen, and the manuscript will separate analyses fixed in advance from those decided after results were seen. The occupancy models and early versions of the trait model were fitted before the log was started; the log records which of their results had already been seen.

## Data Handoff (Inputs)
The core analytical datasets are inherited from the `arha-macroecology` pipeline. Placed in `data/raw/`:
- `host_traits_final.rds`: Includes COMBINE life-history traits, PC scores, and the original categorical `synanthropy_status` for downstream validation.
- `Project_ArHa_database_[DATE].rds`: The harmonised global surveillance database (GPS coordinates, effort, assay results).
- `mammal_tree_matched.tre`: The pruned Upham consensus phylogeny.
- `ecke_synanthropy_clean.rds`: Ecke synanthropy categories, for the planned validation step.

Also required locally but not in the repository: the 1 km Global Human Modification raster (`data/spatial_covariates/gHM_1km_global.tif`), and IUCN range polygons (`data/raw/iucn/`, under IUCN terms of use; not currently used by the pipeline).

---

## Repository structure

```
ArHa_synanthropy/
│
├── data/
│   ├── raw/                 # Hand-off from arha-macroecology, IUCN polygons (read-only)
│   ├── spatial_covariates/  # Global rasters (gHM) - gitignored
│   └── processed/           # GBIF background, model-ready points, synanthropy indices
│
├── scripts/
│   ├── 00_setup.R           # Packages, directories, terra temp directory, parallel plan
│   ├── 01_gbif_pull.R       # GBIF download and cleaning for the target-group background
│   ├── 02_spatial_prep.R    # gHM extraction and per-species background pools
│   ├── 03_occupancy_fit.R   # Integrated ArHa + GBIF occupancy models, one per species
│   ├── 04_extract_syn.R     # Synanthropy indices from the occupancy posteriors
│   └── 05_predictive_ML.R   # Phylogenetic Beta model (brms) predicting unmodelled species
│
├── outputs/
│   ├── models/              # Fitted occupancy and trait models
│   ├── figures/
│   └── tables/
│
└── manuscript/              # Not yet started
```

Scripts are run in numbered order from the project root; each sources `00_setup.R`. `01_gbif_pull.R` needs the `GBIF_USER`, `GBIF_PWD` and `GBIF_EMAIL` environment variables.

## Pipeline as implemented

1. **GBIF background (01).** Records for the species in `host_traits_final.rds`, coordinate uncertainty ≤ 10 km, cleaned with CoordinateCleaner, at most 5,000 random records per species.
2. **Spatial preparation (02).** ArHa and GBIF rodent records from 2000 onwards are pooled and gHM is extracted at each unique location. Each species' background pool is the set of pooled points inside the convex hull of its own records plus a 50 km buffer.
3. **Occupancy models (03).** One integrated, non-spatial model per species (`spOccupancy::intPGOcc`): occupancy `~ gHM + gHM²`; detection `~ effort` for ArHa (log trapping effort per site-session) and intercept-only for GBIF. Species: ArHa species with ≥ 15 site-session events, plus one GBIF-recorded species from each of 150 phylogenetic clusters.
4. **Synanthropy index (04).** Species with max Rhat < 1.15 are kept. For each posterior draw, predicted occupancy is evaluated across the species' own 2.5–97.5% gHM range. The primary index is the occupancy-weighted mean gHM, with the weighted SD as tolerance. Sensitivity summaries: gHM at maximum occupancy, the high-minus-low occupancy contrast, and the quadratic vertex.
5. **Trait model (05).** Beta regression in brms with a phylogenetic random effect, using life-history PCs, range size, habitat breadth and diet breadth (missing values imputed with missRanger), to predict the index for unmodelled rodent species.

## Status

- [x] 01 GBIF background generated.
- [x] 02 Model-ready spatial points generated.
- [x] 03 Occupancy models fitted (150 species attempted; 132 fitted, 18 failed or had insufficient data).
- [ ] 04 Re-run with the new index definition (code updated, not yet run).
- [ ] 05 Re-run on the new index (code updated, not yet run). An earlier version used BART, and a PGLMM on the vertex index was fitted interactively.
- [ ] Validation against Ecke synanthropy categories.
- [ ] Occupancy model: add NNGP spatial random effects.
- [ ] Occupancy model: add night-time lights and accessibility to the detection sub-model.
- [ ] Background pools: use IUCN range polygons.
- [ ] Covariates: add GHSL and natural-habitat covariates (tree cover/NDVI).
- [ ] Extend to Eulipotyphla.
- [ ] Sensitivity refit of 03 with gHM entered as a binned (20 equal-width bins) random intercept.
- [ ] Manuscript.

## Dependencies
- **R packages** (loaded with `pacman` in `00_setup.R`): `tidyverse`, `here`, `terra`, `tidyterra`, `ape`, `spOccupancy`, `brms` (with `cmdstanr`), `rgbif`, `CoordinateCleaner`, `missRanger`, `tidybayes`, `future`, `furrr`.
