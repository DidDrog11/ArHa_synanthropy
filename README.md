# ArHa Synanthropy: A Continuous Spatial Framework for Predicting Rodent Zoonotic Hazard

## Project Overview
This project builds upon the ArHa macroecology dataset to explore synanthropy in small mammals. Instead of treating synanthropy as a rigid categorical trait (e.g., Totally, Occasionally, or Not Synanthropic), we model it as a continuous, spatial reaction process. 

By integrating global surveillance records with Target-Group Background (TGB) data from GBIF, we use multi-source Bayesian occupancy models to quantify species-specific tolerance to anthropogenic modification at a 1km resolution. We then leverage imputed life-history traits and phylogeny to predict a continuous "Synanthropy Index" for over 2,000 unclassified rodent species, assisting in the mapping of global zoonotic contact hazard.

## Data Handoff (Inputs)
The core analytical datasets are inherited from the `arha-macroecology` pipeline. Placed in `data/raw/`:
- `host_traits_final.rds`: Includes COMBINE life-history traits, PC scores, and the original categorical `synanthropy_status` for downstream validation.
- `Project_ArHa_database_[DATE].rds`: The harmonised global surveillance database (GPS coordinates, effort, assay results).
- `mammal_tree_matched.tre`: The pruned Upham consensus phylogeny.

---

## Analytical Pipeline & To-Do List

```
arha_synanthropy/
│
├── data/
│   ├── raw/                 # The cleaned dataset from Paper 1, IUCN polygons
│   ├── spatial_covariates/  # Downloaded global rasters (GHSL, VIIRS, gHM) - *Add to .gitignore!*
│   └── processed/           # Extracted raster values, target-group backgrounds, model-ready data
│
├── scripts/
│   ├── 00_setup.R           # Library loading, setting terra tmp directories, API keys (GBIF)
│   ├── 01_gbif_pull.R       # Querying GBIF for target-group background (TGB) absences
│   ├── 02_extract_covars.R  # Extracting environmental values at ArHa and GBIF points
│   ├── 03_occupancy_fit.R   # Running spoccupancy multi-source models (The Heavy Lifter)
│   ├── 04_extract_syn.R     # Extracting continuous synanthropy parameters from posteriors
│   ├── 05_predictive_ML.R   # Using brms/RF to predict synanthropy for the "shadow commensals"
│   └── 06_hazard_map.R      # Intersecting Competence (Paper 1) * Contact (Paper 2)
│
├── outputs/
│   ├── models/              # Saved .rds files of fitted spoccupancy and brms models
│   ├── figures/             # Density curves, hazard maps, tanglegrams
│   └── tables/              # Model summaries
│
├── manuscript/
│   ├── references.bib
│   └── synanthropy_ms.qmd   # Quarto or RMarkdown document
│
├── .gitignore               # Ignore large rasters, model .rds files, and sensitive keys
└── README.md                # Overview of the pipeline and link to Paper 1's repo
```

### Phase 1: Setup and Data Acquisition
- [ ] **`00_setup.R`**
  - [ ] Initialize repository structure (create `data/spatial_covariates`, `outputs/models`, etc.).
  - [ ] Load standard libraries (`terra`, `sf`, `tidyverse`, `spoccupancy`, `brms`).
  - [ ] Configure `terra` temp directories to handle heavy spatial processing.
- [ ] **`01_gbif_pull.R`**
  - [ ] Authenticate `rgbif` using API credentials.
  - [ ] Pull coordinate data for all target `Rodentia` and `Eulipotyphla` species.
  - [ ] Clean and filter GBIF data (remove fossil records, spatial outliers, high-uncertainty coordinates).
  - [ ] Generate the Target-Group Background (TGB) dataset to act as structured pseudo-absences.

### Phase 2: Spatial Processing
- [ ] **`02_extract_covars.R`**
  - [ ] Download and store global 1km anthropogenic rasters (e.g., Global Human Modification Index, GHSL) in `data/spatial_covariates/`.
  - [ ] Load and clean IUCN range polygons (`MAMMALS_TERRESTRIAL_ONLY.shp`), filtering to target orders.
  - [ ] Build a parallel `terra` workflow (via `furrr`) to loop through species, using their IUCN polygon to crop/mask the global 1km rasters.
  - [ ] Intersect ArHa surveillance points (presences/absences) and GBIF points with the cropped 1km rasters.
  - [ ] Extract point-level environmental covariates (Human Modification, baseline natural habitat like Tree Cover/NDVI).
  - [ ] Extract sampling effort covariates (e.g., distance to roads, night-time lights) to control for the "streetlight effect" in the observation sub-model.

### Phase 3: Occupancy Modelling
- [ ] **`03_occupancy_fit.R`**
  - [ ] Structure the data for multi-source joint likelihood modelling.
  - [ ] Fit `spoccupancy` models for data-rich species:
    - *Biological Process ($\psi$)*: Function of continuous human modification.
    - *Observation Process ($p$)*: Function of ArHa tested counts (offset), night-time lights, and accessibility.
  - [ ] Incorporate spatial random effects (NNGPs) to absorb unmeasured spatial autocorrelation.
- [ ] **`04_extract_syn.R`**
  - [ ] Extract posterior distributions of the $\beta$ coefficients for the human modification covariate.
  - [ ] Generate probability density functions across the modification gradient for each species.
  - [ ] Calculate the continuous **Synanthropy Index** (e.g., the optimum modification peak and tolerance variance) and save as a new species-level trait dataset.
  - [ ] *Validation step:* Plot these continuous curves grouped by the original Ecke `synanthropy_status` to confirm biological alignment.

### Phase 4: Prediction and Hazard Mapping
- [ ] **`05_predictive_ML.R`**
  - [ ] Merge the new continuous Synanthropy Index with `host_traits_final.rds`.
  - [ ] Fit a phylogenetic mixed model (`brms`) or phylogenetic Random Forest using `pace_of_life_pc1`, body mass, diet breadth, and phylogeny to predict the Synanthropy Index.
  - [ ] Project predictions onto the >2,000 unclassified species ("shadow commensals"), capturing posterior uncertainty.
- [ ] **`06_hazard_map.R`**
  - [ ] Decouple competence from contact: Multiply the *Competence Hazard* (from Paper 1) by the *Contact Hazard* (predicted synanthropy from Phase 4).
  - [ ] Generate global 20km bivariate hazard maps using `terra`.
  - [ ] Identify high-risk ecotones (e.g., peri-domestic interfaces and agricultural frontiers) where high competence meets high synanthropic plasticity.

## Dependencies
- **R Packages**: `terra`, `tidyterra`, `sf`, `spoccupancy`, `brms`, `rgbif`, `furrr`.
