
# ============================================================
# EAST AFRICAN SAND FLY ATLAS
# Demonstration Species Distribution Modelling Pipeline
#
# Purpose
# Demonstrate how the atlas can support predictive
# habitat suitability modelling for medically important
# sand fly species.
#
# Target species
# • Phlebotomus orientalis
# • Phlebotomus martini
#
# Environmental predictors
# • CHELSA bioclimatic variables
# • Elevation
# • Soil clay content
# • Continuous tree-cover proxy
#
# Background strategy
# • Kernel-density background
#
# Algorithms
# • Generalized Additive Model (GAM)
# • Random Forest (RF)
# • MaxEnt (maxnet)
#
# Model evaluation
# • Spatial block cross-validation (5 folds)
# • AUC
# • TSS
#
# Ensemble
# • TSS-weighted ensemble prediction
#
# Outputs
# • Habitat suitability maps
# • Prediction uncertainty maps
# • Ensemble variable importance
# • Ensemble response curves
# • Publication-quality figures
#
# ============================================================

rm(list = ls())
graphics.off()
cat("\014")

set.seed(12345)

# ============================================================
# 01. LOAD REQUIRED PACKAGES
# ============================================================

required_packages <- c(
  
  # Spatial
  "terra",
  "sf",
  
  # Data manipulation
  "dplyr",
  "tidyr",
  
  # Modelling
  "mgcv",
  "ranger",
  "maxnet",
  
  # Variable selection
  "usdm",
  "corrplot",
  
  # Model evaluation
  "pROC",
  "blockCV",
  
  # Mapping
  "rnaturalearth",
  "rnaturalearthdata",
  
  # Graphics
  "ggplot2",
  "patchwork",
  "viridis",
  "tidyterra"
)

missing_packages <- required_packages[
  !required_packages %in%
    installed.packages()[,"Package"]
]

if(length(missing_packages) > 0){
  install.packages(missing_packages)
}

invisible(
  lapply(
    required_packages,
    library,
    character.only = TRUE
  )
)

sf::sf_use_s2(FALSE)

# ============================================================
# 02. PROJECT DIRECTORIES AND INPUT DATA
# ============================================================

project_dir <- "/Users/kagboka/Desktop/ACTA-REVISION"

#------------------------------------------------------------
# Occurrence data
#------------------------------------------------------------

occ_path <- file.path(
  project_dir,
  "dashboard_1_data.csv"
)

#------------------------------------------------------------
# CHELSA bioclimatic variables
#------------------------------------------------------------

bio_dir <- file.path(
  project_dir,
  "CHELSA_bio11_1981-2010_V.2.1"
)

#------------------------------------------------------------
# Elevation
#------------------------------------------------------------

elevation_path <- file.path(
  project_dir,
  "elevation_1KMmd_GMTEDmd.tif"
)

#------------------------------------------------------------
# Soil clay content
#------------------------------------------------------------

clay_path <- file.path(
  project_dir,
  "Clay.tif"
)

#------------------------------------------------------------
# Continuous tree-cover proxy
# EarthEnv Consensus Land Cover:
# Class 4 = Mixed/Other Trees
#------------------------------------------------------------

treecover_path <- file.path(
  project_dir,
  "consensus_full_class_4.tif"
)

#------------------------------------------------------------
# Output directory
#------------------------------------------------------------

out_dir <- file.path(
  project_dir,
  "SDM_RESULTS"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

output_folders <- c(
  "00_QC",
  "01_Data",
  "02_Variable_Selection",
  "03_Background_Sampling",
  "04_Cross_Validation",
  "05_Fitted_Models",
  "06_Habitat_Maps",
  "07_Ensembles",
  "08_Response_Curves",
  "09_Variable_Importance",
  "10_Uncertainty",
  "11_Model_Evaluation",
  "12_Publication_Figures",
  "13_Supplementary"
)

for(folder in output_folders){
  
  dir.create(
    file.path(
      out_dir,
      folder
    ),
    recursive = TRUE,
    showWarnings = FALSE
  )
}

#------------------------------------------------------------
# Check that all required files exist
#------------------------------------------------------------

required_files <- c(
  occ_path,
  elevation_path,
  clay_path,
  treecover_path
)

missing_files <- required_files[
  !file.exists(required_files)
]

if(length(missing_files) > 0){
  
  stop(
    "The following required input files were not found:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
  
}

# ============================================================
# 03. TARGET SPECIES
# ============================================================

species_list <- c(
  "Phlebotomus orientalis",
  "Phlebotomus martini"
)

# ============================================================
# 04. GLOBAL SETTINGS
# ============================================================

seed_number <- 12345

set.seed(seed_number)

#------------------------------------------------------------
# Spatial block cross-validation
#------------------------------------------------------------

k_folds <- 5

#------------------------------------------------------------
# Kernel-density background
#------------------------------------------------------------

n_background <- 1000

background_method <- "Kernel"

#------------------------------------------------------------
# Algorithms
#------------------------------------------------------------

algorithm_list <- c(
  "GAM",
  "RF",
  "MAXNET"
)

#------------------------------------------------------------
# Ensemble
#------------------------------------------------------------

ensemble_method <- "TSS_weighted"

#------------------------------------------------------------
# Predictor selection
#------------------------------------------------------------

correlation_threshold <- 0.30

vif_threshold <- 5

#------------------------------------------------------------
# MaxEnt settings
#------------------------------------------------------------

maxent_features <- "lqph"

# ============================================================
# MODELLING STRATEGY
# ============================================================

# Candidate predictors:
#
# 1. CHELSA bioclimatic variables
# 2. Elevation
# 3. Soil clay content
# 4. Continuous tree-cover proxy
#
# All predictors are continuous.
#
# Predictor screening:
#
# 1. Correlation clustering
# 2. Variance Inflation Factor (VIF)
#
# The final SDMs are calibrated using:
#
# • Kernel-density background
# • GAM
# • Random Forest
# • MaxEnt
#
# Algorithm-specific predictions are combined using
# cross-validated TSS weights.
#
# No categorical land-cover variable is used.
#
# ============================================================
# END OF INITIAL SETTINGS




# ============================================================
# 05. OCCURRENCE DATA QUALITY CONTROL
# ============================================================

message("Loading occurrence data...")

#------------------------------------------------------------
# Read occurrence dataset
#------------------------------------------------------------

occ_raw <- read.csv(
  occ_path,
  stringsAsFactors = FALSE
)

names(occ_raw) <- tolower(names(occ_raw))

#------------------------------------------------------------
# Detect required columns
#------------------------------------------------------------

species_col <- grep(
  "species",
  names(occ_raw),
  value = TRUE
)[1]

lon_col <- grep(
  "lon|longitude|decimallongitude",
  names(occ_raw),
  value = TRUE
)[1]

lat_col <- grep(
  "lat|latitude|decimallatitude",
  names(occ_raw),
  value = TRUE
)[1]

occ <- occ_raw |>
  
  dplyr::rename(
    
    species = all_of(species_col),
    lon = all_of(lon_col),
    lat = all_of(lat_col)
    
  )

#------------------------------------------------------------
# Keep target species
#------------------------------------------------------------

occ <- occ |>
  
  dplyr::filter(
    
    grepl(
      "orientalis|martini",
      species,
      ignore.case = TRUE
    )
    
  )

#------------------------------------------------------------
# Standardize species names
#------------------------------------------------------------

occ$species <- dplyr::case_when(
  
  grepl(
    "orientalis",
    occ$species,
    ignore.case = TRUE
  ) ~ "Phlebotomus orientalis",
  
  grepl(
    "martini",
    occ$species,
    ignore.case = TRUE
  ) ~ "Phlebotomus martini",
  
  TRUE ~ occ$species
  
)

#------------------------------------------------------------
# Remove incomplete records
#------------------------------------------------------------

occ <- occ |>
  
  dplyr::filter(
    
    !is.na(lon),
    !is.na(lat)
    
  )

#------------------------------------------------------------
# Remove impossible coordinates
#------------------------------------------------------------

occ <- occ |>
  
  dplyr::filter(
    
    lon >= -180,
    lon <= 180,
    lat >= -90,
    lat <= 90
    
  )

#------------------------------------------------------------
# Remove duplicate occurrences
#------------------------------------------------------------

occ <- occ |>
  
  dplyr::distinct(
    
    species,
    lon,
    lat,
    .keep_all = TRUE
    
  )

#------------------------------------------------------------
# East African study area
#------------------------------------------------------------

east_africa <- rnaturalearth::ne_countries(
  
  country = c(
    "Sudan",
    "South Sudan",
    "Eritrea",
    "Djibouti",
    "Ethiopia",
    "Somalia",
    "Kenya",
    "Uganda",
    "Tanzania",
    "Rwanda",
    "Burundi"
  ),
  
  returnclass = "sf"
)

#------------------------------------------------------------
# Keep occurrences inside study area
#------------------------------------------------------------

occ_sf <- sf::st_as_sf(
  
  occ,
  
  coords = c("lon","lat"),
  
  crs = 4326,
  
  remove = FALSE
  
)

occ_sf <- sf::st_filter(
  
  occ_sf,
  
  east_africa,
  
  .predicate = sf::st_intersects
  
)

occ <- sf::st_drop_geometry(
  occ_sf
)

#------------------------------------------------------------
# Save cleaned dataset
#------------------------------------------------------------

write.csv(
  
  occ,
  
  file.path(
    out_dir,
    "01_Data",
    "Clean_occurrences.csv"
  ),
  
  row.names = FALSE
  
)

#------------------------------------------------------------
# Summary
#------------------------------------------------------------

cat("\nOccurrence summary\n")
print(table(occ$species))
cat("\nTotal records:", nrow(occ), "\n")

message("Section 05 completed successfully.")# ============================================================





# ============================================================
# 06. LOAD AND PREPARE ENVIRONMENTAL PREDICTORS
# ============================================================

message("Loading environmental predictors...")

#------------------------------------------------------------
# East African study area
#------------------------------------------------------------

east_vect <- terra::vect(east_africa)

#------------------------------------------------------------
# CHELSA bioclimatic variables
#------------------------------------------------------------

bio_files <- list.files(
  bio_dir,
  pattern = "\\.tif$",
  full.names = TRUE
)

bio_number <- as.numeric(
  gsub(
    ".*bio_?([0-9]+).*",
    "\\1",
    basename(bio_files)
  )
)

bio_files <- bio_files[
  order(bio_number)
]

stopifnot(length(bio_files) == 19)

bio <- terra::rast(bio_files)

names(bio) <- paste0("bio",1:19)

#------------------------------------------------------------
# Additional continuous predictors
#------------------------------------------------------------

elevation <- terra::rast(elevation_path)
names(elevation) <- "Elevation"

clay <- terra::rast(clay_path)
names(clay) <- "Clay"

treecover <- terra::rast(treecover_path)
names(treecover) <- "TreeCover"

#------------------------------------------------------------
# Crop CHELSA
#------------------------------------------------------------

bio <- terra::crop(
  bio,
  east_vect
)

bio <- terra::mask(
  bio,
  east_vect
)

#------------------------------------------------------------
# Match projections
#------------------------------------------------------------

if(!terra::same.crs(elevation,bio)){
  
  elevation <- terra::project(
    elevation,
    bio[[1]],
    method="bilinear"
  )
  
}

if(!terra::same.crs(clay,bio)){
  
  clay <- terra::project(
    clay,
    bio[[1]],
    method="bilinear"
  )
  
}

if(!terra::same.crs(treecover,bio)){
  
  treecover <- terra::project(
    treecover,
    bio[[1]],
    method="bilinear"
  )
  
}

#------------------------------------------------------------
# Crop and mask
#------------------------------------------------------------

predictor_list <- list(
  elevation,
  clay,
  treecover
)

predictor_list <- lapply(
  
  predictor_list,
  
  function(r){
    
    r <- terra::crop(
      r,
      east_vect
    )
    
    r <- terra::mask(
      r,
      east_vect
    )
    
    r
    
  }
)

elevation <- predictor_list[[1]]
clay <- predictor_list[[2]]
treecover <- predictor_list[[3]]

#------------------------------------------------------------
# Aggregate to ~5 km
#------------------------------------------------------------

bio <- terra::aggregate(
  bio,
  fact = 5,
  fun = mean,
  na.rm = TRUE
)

elevation <- terra::aggregate(
  elevation,
  fact = 5,
  fun = mean,
  na.rm = TRUE
)

clay <- terra::aggregate(
  clay,
  fact = 5,
  fun = mean,
  na.rm = TRUE
)

treecover <- terra::aggregate(
  treecover,
  fact = 5,
  fun = mean,
  na.rm = TRUE
)

#------------------------------------------------------------
# Align grids
#------------------------------------------------------------

elevation <- terra::resample(
  elevation,
  bio[[1]],
  method="bilinear"
)

clay <- terra::resample(
  clay,
  bio[[1]],
  method="bilinear"
)

treecover <- terra::resample(
  treecover,
  bio[[1]],
  method="bilinear"
)

#------------------------------------------------------------
# Build predictor stack
#------------------------------------------------------------

predictor_stack <- c(
  bio,
  elevation,
  clay,
  treecover
)

#------------------------------------------------------------
# Save prepared predictors
#------------------------------------------------------------

terra::writeRaster(
  predictor_stack,
  file.path(
    out_dir,
    "01_Data",
    "Environmental_predictors_5km.tif"
  ),
  overwrite = TRUE
)

cat("\nEnvironmental predictors loaded\n")
cat("---------------------------------\n")
cat("Bioclimatic variables :", terra::nlyr(bio), "\n")
cat("Additional predictors :", 3, "\n")
cat("Total predictors      :", terra::nlyr(predictor_stack), "\n")

message("Section 06 completed successfully.")





# ============================================================
# 07. VARIABLE SELECTION
# ============================================================

message("Selecting environmental predictors...")

#------------------------------------------------------------
# Random sample of predictor values
#------------------------------------------------------------

env_sample <- terra::spatSample(
  predictor_stack,
  size = 10000,
  method = "random",
  na.rm = TRUE,
  as.df = TRUE
)

env_sample <- na.omit(env_sample)

#------------------------------------------------------------
# Correlation matrix
#------------------------------------------------------------

cor_matrix <- cor(
  env_sample,
  method = "pearson"
)

write.csv(
  cor_matrix,
  file.path(
    out_dir,
    "02_Variable_Selection",
    "Correlation_matrix.csv"
  )
)

#------------------------------------------------------------
# Correlation heatmap
#------------------------------------------------------------

png(
  file.path(
    out_dir,
    "02_Variable_Selection",
    "Correlation_heatmap.png"
  ),
  width = 2200,
  height = 2000,
  res = 300
)

corrplot::corrplot(
  cor_matrix,
  method = "color",
  type = "upper",
  tl.cex = 0.7,
  number.cex = 0.5
)

dev.off()

#------------------------------------------------------------
# Correlation clustering
#------------------------------------------------------------

distance_matrix <- as.dist(
  1 - abs(cor_matrix)
)

hc <- hclust(
  distance_matrix,
  method = "average"
)

clusters <- cutree(
  hc,
  h = correlation_threshold
)

cluster_table <- data.frame(
  Variable = names(clusters),
  Cluster = clusters
)

selected_variables <- c()

for(cl in unique(cluster_table$Cluster)){
  
  vars <- cluster_table$Variable[
    cluster_table$Cluster == cl
  ]
  
  if(length(vars) == 1){
    
    selected_variables <- c(
      selected_variables,
      vars
    )
    
  } else {
    
    sub_cor <- abs(
      cor_matrix[vars, vars]
    )
    
    mean_cor <- rowMeans(sub_cor)
    
    selected_variables <- c(
      selected_variables,
      names(which.min(mean_cor))
    )
  }
}

cat("\nVariables after correlation filtering\n")
print(selected_variables)

#------------------------------------------------------------
# Variance Inflation Factor
#------------------------------------------------------------

vif_results <- usdm::vifstep(
  env_sample[, selected_variables],
  th = vif_threshold
)

final_variables <- vif_results@results$Variables

write.csv(
  vif_results@results,
  file.path(
    out_dir,
    "02_Variable_Selection",
    "VIF_results.csv"
  ),
  row.names = FALSE
)

cat("\nVariables retained after VIF\n")
print(final_variables)

#------------------------------------------------------------
# Final predictor stack
#------------------------------------------------------------

final_predictors <- predictor_stack[[final_variables]]

terra::writeRaster(
  final_predictors,
  file.path(
    out_dir,
    "02_Variable_Selection",
    "Final_predictor_stack.tif"
  ),
  overwrite = TRUE
)

write.csv(
  data.frame(
    Variable = names(final_predictors)
  ),
  file.path(
    out_dir,
    "02_Variable_Selection",
    "Selected_predictors.csv"
  ),
  row.names = FALSE
)

cat("\n===============================\n")
cat("FINAL PREDICTORS\n")
cat("===============================\n")

print(names(final_predictors))

message("Section 07 completed successfully.")






# ============================================================
# 08. KERNEL-DENSITY BACKGROUND SAMPLING
# ============================================================

message("Generating kernel-density background...")

library(MASS)

#------------------------------------------------------------
# Study area
#------------------------------------------------------------

study_area <- final_predictors[[1]]

#------------------------------------------------------------
# Occurrence coordinates
#------------------------------------------------------------

occ_xy <- occ[, c("lon","lat")]

#------------------------------------------------------------
# Kernel density estimation
#------------------------------------------------------------

ext <- terra::ext(study_area)

kernel <- MASS::kde2d(
  occ_xy$lon,
  occ_xy$lat,
  n = 300,
  lims = c(
    ext$xmin,
    ext$xmax,
    ext$ymin,
    ext$ymax
  )
)

kernel_raster <- terra::rast(
  nrows = nrow(kernel$z),
  ncols = ncol(kernel$z),
  xmin = min(kernel$x),
  xmax = max(kernel$x),
  ymin = min(kernel$y),
  ymax = max(kernel$y),
  crs = terra::crs(study_area)
)

values(kernel_raster) <-
  as.vector(
    t(kernel$z[nrow(kernel$z):1, ])
  )

#------------------------------------------------------------
# Match modelling grid
#------------------------------------------------------------

kernel_raster <- terra::resample(
  kernel_raster,
  study_area,
  method = "bilinear"
)

kernel_raster <-
  kernel_raster /
  terra::global(
    kernel_raster,
    "max",
    na.rm = TRUE
  )[1,1]

kernel_raster <- terra::mask(
  kernel_raster,
  study_area
)

#------------------------------------------------------------
# Sample background
#------------------------------------------------------------

bg_kernel <- terra::spatSample(
  kernel_raster,
  size = n_background,
  method = "weights",
  xy = TRUE,
  na.rm = TRUE
)

bg_kernel <- as.data.frame(bg_kernel)

names(bg_kernel)[1:2] <- c(
  "lon",
  "lat"
)

cat(
  "\nBackground points:",
  nrow(bg_kernel),
  "\n"
)

#------------------------------------------------------------
# Save background
#------------------------------------------------------------

write.csv(
  bg_kernel,
  file.path(
    out_dir,
    "03_Background_Sampling",
    "Kernel_background.csv"
  ),
  row.names = FALSE
)

#------------------------------------------------------------
# Publication figure
#------------------------------------------------------------

p <- ggplot() +
  
  geom_spatraster(
    data = kernel_raster
  ) +
  
  scale_fill_viridis_c(
    name = "Kernel density"
  ) +
  
  geom_point(
    data = bg_kernel,
    aes(lon, lat),
    colour = "grey40",
    alpha = 0.35,
    size = 0.25
  ) +
  
  geom_point(
    data = occ,
    aes(lon, lat),
    colour = "red",
    size = 0.8
  ) +
  
  theme_bw() +
  
  labs(
    title = "Kernel-density background sampling"
  )

ggsave(
  file.path(
    out_dir,
    "12_Publication_Figures",
    "Kernel_background_sampling.png"
  ),
  p,
  width = 7,
  height = 7,
  dpi = 600
)

message("Section 08 completed successfully.")






# ============================================================
# 09. SPATIAL BLOCK CROSS-VALIDATION
# ============================================================

message("Creating spatial block cross-validation...")

#------------------------------------------------------------
# Presence locations
#------------------------------------------------------------

presence_sf <- sf::st_as_sf(
  occ,
  coords = c("lon","lat"),
  crs = 4326,
  remove = FALSE
)

#------------------------------------------------------------
# Create spatial blocks
#------------------------------------------------------------

sb <- blockCV::cv_spatial(
  x = presence_sf,
  r = final_predictors,
  size = 100000,
  k = k_folds,
  selection = "random",
  iteration = 100,
  seed = seed_number,
  plot = FALSE,
  progress = TRUE
)

#------------------------------------------------------------
# Assign fold IDs
#------------------------------------------------------------

presence_sf$Fold <- sb$folds_ids

presence_fold_data <- sf::st_drop_geometry(
  presence_sf
)

#------------------------------------------------------------
# Save folds
#------------------------------------------------------------

write.csv(
  presence_fold_data,
  file.path(
    out_dir,
    "04_Cross_Validation",
    "Presence_folds.csv"
  ),
  row.names = FALSE
)

#------------------------------------------------------------
# Save polygons
#------------------------------------------------------------

sf::st_write(
  sb$blocks,
  file.path(
    out_dir,
    "04_Cross_Validation",
    "Spatial_blocks.gpkg"
  ),
  delete_dsn = TRUE,
  quiet = TRUE
)

#------------------------------------------------------------
# Fold summary
#------------------------------------------------------------

fold_summary <- presence_fold_data |>
  dplyr::count(
    Fold,
    name = "Occurrences"
  )

write.csv(
  fold_summary,
  file.path(
    out_dir,
    "04_Cross_Validation",
    "Fold_summary.csv"
  ),
  row.names = FALSE
)

print(fold_summary)

#------------------------------------------------------------
# Publication figure
#------------------------------------------------------------

p <- ggplot() +
  
  geom_sf(
    data = east_africa,
    fill = "grey95",
    colour = "grey40",
    linewidth = 0.3
  ) +
  
  geom_sf(
    data = sb$blocks,
    aes(fill = factor(folds)),
    alpha = 0.35,
    colour = "black",
    linewidth = 0.2
  ) +
  
  geom_sf(
    data = presence_sf,
    aes(colour = factor(Fold)),
    size = 1.1
  ) +
  
  scale_fill_viridis_d(
    name = "Fold"
  ) +
  
  scale_colour_viridis_d(
    guide = "none"
  ) +
  
  coord_sf() +
  
  theme_bw(base_size = 13) +
  
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  
  labs(
    title = "Spatial block cross-validation",
    subtitle = paste(
      k_folds,
      "folds (100-km blocks)"
    )
  )

ggsave(
  file.path(
    out_dir,
    "12_Publication_Figures",
    "Figure1_SpatialBlocks.png"
  ),
  p,
  width = 8,
  height = 7,
  dpi = 600
)

message("Section 09 completed successfully.")







# ============================================================
# 10. MODELLING FUNCTIONS
# ============================================================

message("Preparing modelling functions...")

#------------------------------------------------------------
# Extract environmental predictors
#------------------------------------------------------------

extract_environment <- function(points){
  
  pts <- terra::vect(
    points,
    geom = c("lon","lat"),
    crs = "EPSG:4326"
  )
  
  env <- terra::extract(
    final_predictors,
    pts
  )
  
  env <- env[, -1, drop = FALSE]
  
  dplyr::bind_cols(
    points,
    env
  )
  
}

#------------------------------------------------------------
# Prepare modelling dataset
#------------------------------------------------------------

prepare_model_frame <- function(dat){
  
  dat %>%
    dplyr::select(
      Presence,
      all_of(names(final_predictors))
    ) %>%
    na.omit()
  
}

#------------------------------------------------------------
# Generalized Additive Model
#------------------------------------------------------------

fit_gam <- function(dat){
  
  dat <- prepare_model_frame(dat)
  
  smooth_terms <- paste0(
    "s(",
    names(final_predictors),
    ",k=5)"
  )
  
  form <- as.formula(
    paste(
      "Presence ~",
      paste(
        smooth_terms,
        collapse = " + "
      )
    )
  )
  
  mgcv::gam(
    form,
    data = dat,
    family = binomial(),
    method = "REML"
  )
  
}

#------------------------------------------------------------
# Random Forest
#------------------------------------------------------------

fit_rf <- function(dat){
  
  dat <- prepare_model_frame(dat)
  
  dat$Presence <- factor(
    dat$Presence,
    levels = c(0,1)
  )
  
  ranger::ranger(
    Presence ~ .,
    data = dat,
    probability = TRUE,
    num.trees = 500,
    importance = "permutation",
    seed = seed_number
  )
  
}

#------------------------------------------------------------
# MaxEnt
#------------------------------------------------------------

fit_maxnet <- function(dat){
  
  dat <- prepare_model_frame(dat)
  
  y <- dat$Presence
  
  x <- dat |>
    dplyr::select(-Presence)
  
  model <- maxnet::maxnet(
    
    p = y,
    
    data = x,
    
    f = maxnet::maxnet.formula(
      p = y,
      data = x,
      classes = maxent_features
    )
    
  )
  
  model
  
}

#------------------------------------------------------------




#------------------------------------------------------------
# Master fitting function
#------------------------------------------------------------

fit_model <- function(algorithm, dat){
  
  switch(
    
    algorithm,
    
    GAM = fit_gam(dat),
    
    RF = fit_rf(dat),
    
    MAXNET = fit_maxnet(dat)
    
  )
  
}

#------------------------------------------------------------
# Prediction function
#------------------------------------------------------------

predict_model <- function(model, newdata, algorithm){
  
  x <- newdata %>%
    dplyr::select(
      all_of(names(final_predictors))
    )
  
  if(algorithm == "GAM"){
    
    return(
      as.numeric(
        predict(
          model,
          newdata = x,
          type = "response"
        )
      )
    )
    
  }
  
  if(algorithm == "RF"){
    
    pr <- predict(
      model,
      data = x
    )$predictions
    
    if(is.matrix(pr)){
      
      if("1" %in% colnames(pr)){
        
        return(as.numeric(pr[, "1"]))
        
      }else{
        
        return(as.numeric(pr[, ncol(pr)]))
        
      }
      
    }else{
      
      return(as.numeric(pr))
      
    }
    
  }
  
  if(algorithm == "MAXNET"){
    
    return(
      as.numeric(
        predict(
          model,
          x,
          type = "cloglog"
        )
      )
    )
    
  }
  
  stop("Unknown algorithm.")
  
}

message("Section 10 completed successfully.")




# ============================================================
# 11. SPATIAL CROSS-VALIDATION
# ============================================================

message("Running spatial block cross-validation...")

#------------------------------------------------------------
# Extract predictors
#------------------------------------------------------------

extract_environment <- function(points){
  
  pts <- terra::vect(
    points,
    geom = c("lon","lat"),
    crs = "EPSG:4326"
  )
  
  env <- terra::extract(
    final_predictors,
    pts
  )
  
  env <- env[,-1,drop=FALSE]
  
  dplyr::bind_cols(
    points,
    env
  )
  
}

#------------------------------------------------------------
# Assign kernel background to folds
#------------------------------------------------------------

assign_background_folds <- function(bg){
  
  bg_sf <- sf::st_as_sf(
    bg,
    coords = c("lon","lat"),
    crs = 4326,
    remove = FALSE
  )
  
  bg_sf <- sf::st_join(
    bg_sf,
    sb$blocks[, "folds", drop = FALSE],
    join = sf::st_intersects,
    left = TRUE
  )
  
  bg <- sf::st_drop_geometry(bg_sf)
  
  names(bg)[names(bg)=="folds"] <- "Fold"
  
  bg %>%
    dplyr::filter(!is.na(Fold))
  
}

background_fold <- assign_background_folds(
  bg_kernel
)

presence_fold_data <- sf::st_drop_geometry(
  presence_sf
)

#------------------------------------------------------------
# Evaluation function
#------------------------------------------------------------

evaluate_predictions <- function(obs,pred){
  
  keep <- complete.cases(obs,pred)
  
  obs <- obs[keep]
  pred <- pred[keep]
  
  if(length(unique(obs))<2)
    return(NULL)
  
  auc <- as.numeric(
    pROC::auc(obs,pred)
  )
  
  thresholds <- seq(
    0.01,
    0.99,
    by=0.01
  )
  
  tss <- sapply(
    thresholds,
    function(th){
      
      cls <- ifelse(pred>=th,1,0)
      
      TP <- sum(cls==1 & obs==1)
      TN <- sum(cls==0 & obs==0)
      FP <- sum(cls==1 & obs==0)
      FN <- sum(cls==0 & obs==1)
      
      sens <- TP/(TP+FN)
      
      spec <- TN/(TN+FP)
      
      sens+spec-1
      
    }
  )
  
  id <- which.max(tss)
  
  data.frame(
    
    AUC=auc,
    
    TSS=tss[id],
    
    Threshold=thresholds[id]
    
  )
  
}

#------------------------------------------------------------
# Results
#------------------------------------------------------------

cv_results <- list()

counter <- 1

#============================================================
# Species
#============================================================

for(sp in species_list){
  
  message("--------------------------------")
  message(sp)
  
  pres <- presence_fold_data %>%
    
    dplyr::filter(
      species==sp
    )
  
  pres <- extract_environment(
    pres
  )
  
  pres$Presence <- 1
  
  bg <- extract_environment(
    background_fold
  )
  
  bg$Presence <- 0
  
  model_data <- dplyr::bind_rows(
    pres,
    bg
  )
  
  model_data <- model_data %>%
    
    dplyr::select(
      Fold,
      Presence,
      all_of(names(final_predictors))
    ) %>%
    
    na.omit()
  
  #----------------------------------------------------------
  # Folds
  #----------------------------------------------------------
  
  for(fold in seq_len(k_folds)){
    
    message(
      paste(" Fold",fold)
    )
    
    train <- model_data %>%
      dplyr::filter(
        Fold!=fold
      )
    
    test <- model_data %>%
      dplyr::filter(
        Fold==fold
      )
    
    if(length(unique(train$Presence))<2) next
    if(length(unique(test$Presence))<2) next
    
    #--------------------------------------------------------
    # Algorithms
    #--------------------------------------------------------
    
    for(alg in algorithm_list){
      
      message(
        paste("   ",alg)
      )
      
      model <- tryCatch(
        
        fit_model(
          alg,
          train
        ),
        
        error=function(e){
          
          message(e$message)
          
          NULL
          
        }
        
      )
      
      if(is.null(model))
        next
      
      pred <- tryCatch(
        
        predict_model(
          model,
          test,
          alg
        ),
        
        error=function(e){
          
          message(e$message)
          
          NULL
          
        }
        
      )
      
      if(is.null(pred))
        next
      
      metrics <- evaluate_predictions(
        
        test$Presence,
        
        pred
        
      )
      
      if(is.null(metrics))
        next
      
      metrics$Species <- sp
      metrics$Algorithm <- alg
      metrics$Fold <- fold
      
      cv_results[[counter]] <- metrics
      
      counter <- counter+1
      
    }
    
  }
  
}

#------------------------------------------------------------
# Save results
#------------------------------------------------------------

cv_results <- dplyr::bind_rows(
  cv_results
)

write.csv(
  
  cv_results,
  
  file.path(
    
    out_dir,
    
    "11_Model_Evaluation",
    
    "CrossValidation_Results.csv"
    
  ),
  
  row.names=FALSE
  
)

#------------------------------------------------------------
# Mean performance
#------------------------------------------------------------

cv_summary <- cv_results %>%
  
  dplyr::group_by(
    
    Species,
    
    Algorithm
    
  ) %>%
  
  dplyr::summarise(
    
    Mean_AUC=mean(AUC),
    
    SD_AUC=sd(AUC),
    
    Mean_TSS=mean(TSS),
    
    SD_TSS=sd(TSS),
    
    .groups="drop"
    
  )

write.csv(
  
  cv_summary,
  
  file.path(
    
    out_dir,
    
    "11_Model_Evaluation",
    
    "CrossValidation_Summary.csv"
    
  ),
  
  row.names=FALSE
  
)

print(cv_summary)

#------------------------------------------------------------
# TSS ensemble weights
#------------------------------------------------------------

algorithm_weights <- cv_summary %>%
  
  dplyr::group_by(
    
    Species
    
  ) %>%
  
  dplyr::mutate(
    
    Weight = Mean_TSS /
      
      sum(Mean_TSS)
    
  ) %>%
  
  dplyr::ungroup()

write.csv(
  
  algorithm_weights,
  
  file.path(
    
    out_dir,
    
    "11_Model_Evaluation",
    
    "Algorithm_Weights.csv"
    
  ),
  
  row.names=FALSE
  
)

cat("\nTSS weights\n")

print(algorithm_weights)

message("Section 11 completed successfully.")


# ============================================================
# 12. FINAL MODELS AND HABITAT SUITABILITY PREDICTIONS
# ============================================================

message("Fitting final models...")

final_models <- list()

for(sp in species_list){
  
  message("--------------------------------")
  message(sp)
  
  #----------------------------------------------------------
  # Presence
  #----------------------------------------------------------
  
  presence <- occ %>%
    dplyr::filter(
      species == sp
    )
  
  presence <- extract_environment(
    presence
  )
  
  presence$Presence <- 1
  
  #----------------------------------------------------------
  # Background
  #----------------------------------------------------------
  
  background <- extract_environment(
    bg_kernel
  )
  
  background$Presence <- 0
  
  #----------------------------------------------------------
  # Build modelling dataset
  #----------------------------------------------------------
  
  model_data <- dplyr::bind_rows(
    presence,
    background
  )
  
  model_data <- prepare_model_frame(
    model_data
  )
  
  #----------------------------------------------------------
  # Fit algorithms
  #----------------------------------------------------------
  
  for(alg in algorithm_list){
    
    message(paste("  ",alg))
    
    model <- fit_model(
      alg,
      model_data
    )
    
    final_models[[paste(
      sp,
      alg,
      sep="_"
    )]] <- model
    
    saveRDS(
      
      model,
      
      file.path(
        
        out_dir,
        
        "05_Fitted_Models",
        
        paste0(
          
          gsub(" ","_",sp),
          
          "_",
          
          alg,
          
          ".rds"
          
        )
        
      )
      
    )
    
  }
  
}

message("Final models completed.")

#------------------------------------------------------------
# Raster prediction function
#------------------------------------------------------------

predict_raster_model <- function(model,
                                 algorithm){
  
  env <- terra::as.data.frame(
    final_predictors,
    na.rm = FALSE
  )
  
  pred <- rep(
    NA,
    nrow(env)
  )
  
  keep <- complete.cases(env)
  
  pred[keep] <-
    
    predict_model(
      
      model,
      
      env[keep,],
      
      algorithm
      
    )
  
  r <- final_predictors[[1]]
  
  terra::values(r) <- pred
  
  r
  
}

#------------------------------------------------------------
# Predict habitat suitability
#------------------------------------------------------------

prediction_stack <- list()

for(name in names(final_models)){
  
  message(name)
  
  model <- final_models[[name]]
  
  alg <- tail(
    strsplit(name,"_")[[1]],
    1
  )
  
  pred <- predict_raster_model(
    
    model,
    
    alg
    
  )
  
  prediction_stack[[name]] <- pred
  
  terra::writeRaster(
    
    pred,
    
    file.path(
      
      out_dir,
      
      "06_Habitat_Maps",
      
      paste0(name,".tif")
      
    ),
    
    overwrite = TRUE
    
  )
  
}

message("Habitat suitability maps completed.")




# ============================================================
# 13. TSS-WEIGHTED ENSEMBLE
# ============================================================

message("Generating TSS-weighted ensembles...")

ensemble_dir <- file.path(
  out_dir,
  "07_Ensembles"
)

dir.create(
  ensemble_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

#------------------------------------------------------------
# Read weights
#------------------------------------------------------------

weights <- read.csv(
  file.path(
    out_dir,
    "11_Model_Evaluation",
    "Algorithm_Weights.csv"
  )
)

#------------------------------------------------------------
# Species loop
#------------------------------------------------------------

for(sp in species_list){
  
  message("--------------------------------")
  message(sp)
  
  sp_weights <- weights %>%
    dplyr::filter(
      Species == sp
    )
  
  #----------------------------------------------------------
  # Load predictions
  #----------------------------------------------------------
  
  predictions <- list()
  

    for(alg in algorithm_list){
      
      file <- file.path(
        
        out_dir,
        
        "06_Habitat_Maps",
        
        paste0(
          sp,
          "_",
          alg,
          ".tif"
        )
        
       
      )
      
      print(file)
      
      if(!file.exists(file)){
        stop("Cannot find: ", file)
      }
      
      predictions[[alg]] <- terra::rast(file)
      
    }
  
  #----------------------------------------------------------
  # Weighted ensemble
  #----------------------------------------------------------
  
  ensemble <- predictions[[1]] * 0
  
  for(i in seq_along(algorithm_list)){
    
    alg <- algorithm_list[i]
    
    w <- sp_weights$Weight[
      sp_weights$Algorithm == alg
    ]
    
    ensemble <- ensemble +
      predictions[[alg]] * w
    
  }
  
  names(ensemble) <- "Suitability"
  
  #----------------------------------------------------------
  # Prediction uncertainty
  #----------------------------------------------------------
  
  pred_stack <- terra::rast(
    predictions
  )
  
  sd_map <- terra::app(
    pred_stack,
    sd,
    na.rm = TRUE
  )
  
  names(sd_map) <- "SD"
  
  #----------------------------------------------------------
  # Save rasters
  #----------------------------------------------------------
  
  prefix <- gsub(
    " ",
    "_",
    sp
  )
  
  terra::writeRaster(
    
    ensemble,
    
    file.path(
      
      ensemble_dir,
      
      paste0(
        prefix,
        "_Ensemble.tif"
      )
      
    ),
    
    overwrite = TRUE
    
  )
  
  terra::writeRaster(
    
    sd_map,
    
    file.path(
      
      ensemble_dir,
      
      paste0(
        prefix,
        "_Uncertainty.tif"
      )
      
    ),
    
    overwrite = TRUE
    
  )
  
}

message("Section 13 completed successfully.")




#============================================================
# 14. TSS-WEIGHTED ENSEMBLE PREDICTIONS
#============================================================

message("Building TSS-weighted ensembles...")

ensemble_dir <- file.path(
  out_dir,
  "07_Ensembles"
)

dir.create(
  ensemble_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

#------------------------------------------------------------
# Species loop
#------------------------------------------------------------

for(sp in species_list){
  
  message("--------------------------------")
  message(sp)
  
  sp_weights <- algorithm_weights |>
    dplyr::filter(
      Species == sp
    )
  
  predictions <- list()
  
  for(alg in algorithm_list){
    
    file <- file.path(
      
      out_dir,
      
      "06_Habitat_Maps",
      
      paste0(
        sp,
        "_",
        alg,
        ".tif"
      )
      
    )
    
    predictions[[alg]] <- terra::rast(file)
    
  }
  
  #----------------------------------------------------------
  # Weighted ensemble
  #----------------------------------------------------------
  
  ensemble <- predictions[[1]] * 0
  
  for(alg in algorithm_list){
    
    w <- sp_weights$Weight[
      sp_weights$Algorithm == alg
    ]
    
    ensemble <- ensemble +
      predictions[[alg]] * w
    
  }
  
  names(ensemble) <- "Suitability"
  
  #----------------------------------------------------------
  # Prediction uncertainty
  #----------------------------------------------------------
  
  pred_stack <- terra::rast(
    predictions
  )
  
  sd_map <- terra::app(
    pred_stack,
    sd,
    na.rm = TRUE
  )
  
  names(sd_map) <- "SD"
  
  #----------------------------------------------------------
  # Save rasters
  #----------------------------------------------------------
  
  terra::writeRaster(
    
    ensemble,
    
    file.path(
      
      ensemble_dir,
      
      paste0(
        sp,
        "_Ensemble.tif"
      )
      
    ),
    
    overwrite = TRUE
    
  )
  
  terra::writeRaster(
    
    sd_map,
    
    file.path(
      
      ensemble_dir,
      
      paste0(
        sp,
        "_SD.tif"
      )
      
    ),
    
    overwrite = TRUE
    
  )
  
}

message("Section 14 completed successfully.")





#============================================================
# 15. ENSEMBLE VARIABLE IMPORTANCE
#============================================================

message("Computing ensemble variable importance...")

importance_results <- list()

for(sp in species_list){
  
  message("--------------------------------")
  message(sp)
  
  #----------------------------------------------------------
  # Algorithm weights
  #----------------------------------------------------------
  
  sp_weights <- algorithm_weights |>
    dplyr::filter(
      Species == sp
    )
  
  #----------------------------------------------------------
  # Load fitted models
  #----------------------------------------------------------
  
  models <- list()
  
  for(alg in algorithm_list){
    
    models[[alg]] <- readRDS(
      
      file.path(
        
        out_dir,
        
        "05_Fitted_Models",
        
        paste0(
          gsub(" ","_",sp),
          "_",
          alg,
          ".rds"
        )
        
      )
      
    )
    
  }
  
  #----------------------------------------------------------
  # Environmental data
  #----------------------------------------------------------
  
  env <- terra::as.data.frame(
    final_predictors,
    na.rm = TRUE
  )
  
  #----------------------------------------------------------
  # Baseline ensemble prediction
  #----------------------------------------------------------
  
  base <- rep(0,nrow(env))
  
  for(alg in algorithm_list){
    
    pred <- predict_model(
      models[[alg]],
      env,
      alg
    )
    
    w <- sp_weights$Weight[
      sp_weights$Algorithm==alg
    ]
    
    base <- base + pred*w
    
  }
  
  #----------------------------------------------------------
  # Permutation importance
  #----------------------------------------------------------
  
  imp <- data.frame()
  
  for(v in names(final_predictors)){
    
    env_perm <- env
    
    env_perm[[v]] <- sample(
      env_perm[[v]]
    )
    
    perm <- rep(0,nrow(env))
    
    for(alg in algorithm_list){
      
      pred <- predict_model(
        models[[alg]],
        env_perm,
        alg
      )
      
      w <- sp_weights$Weight[
        sp_weights$Algorithm==alg
      ]
      
      perm <- perm + pred*w
      
    }
    
    importance <- mean(
      abs(base-perm),
      na.rm=TRUE
    )
    
    imp <- rbind(
      
      imp,
      
      data.frame(
        
        Species=sp,
        
        Variable=v,
        
        Importance=importance
        
      )
      
    )
    
  }
  
  importance_results[[sp]] <- imp
  
}

importance_results <- dplyr::bind_rows(
  importance_results
)

importance_results <- importance_results |>
  
  dplyr::group_by(Species) |>
  
  dplyr::mutate(
    
    Importance=
      Importance/sum(Importance)*100
    
  )

write.csv(
  
  importance_results,
  
  file.path(
    
    out_dir,
    
    "09_Variable_Importance",
    
    "Ensemble_Variable_Importance.csv"
    
  ),
  
  row.names=FALSE
  
)

print(importance_results)

message("Section 15 completed successfully.")









#============================================================
# 16. ENSEMBLE RESPONSE CURVES
#============================================================

message("Generating ensemble response curves...")

response_dir <- file.path(
  out_dir,
  "08_Response_Curves"
)

dir.create(
  response_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

env <- terra::as.data.frame(
  final_predictors,
  na.rm = TRUE
)

for(sp in species_list){
  
  message("--------------------------------")
  message(sp)
  
  #----------------------------------------------------------
  # Algorithm weights
  #----------------------------------------------------------
  
  sp_weights <- algorithm_weights |>
    dplyr::filter(
      Species == sp
    )
  
  #----------------------------------------------------------
  # Load fitted models
  #----------------------------------------------------------
  
  models <- list()
  
  for(alg in algorithm_list){
    
    models[[alg]] <- readRDS(
      
      file.path(
        
        out_dir,
        
        "05_Fitted_Models",
        
        paste0(
          gsub(" ","_",sp),
          "_",
          alg,
          ".rds"
        )
        
      )
      
    )
    
  }
  
  #----------------------------------------------------------
  # Generate one response curve per predictor
  #----------------------------------------------------------
  
  for(v in names(final_predictors)){
    
    base <- env[rep(1,100),]
    
    for(col in names(env)){
      
      if(col==v){
        
        base[[col]] <- seq(
          
          min(env[[col]],na.rm=TRUE),
          
          max(env[[col]],na.rm=TRUE),
          
          length.out=100
          
        )
        
      }else{
        
        base[[col]] <- median(
          
          env[[col]],
          
          na.rm=TRUE
          
        )
        
      }
      
    }
    
    ensemble_pred <- rep(0,nrow(base))
    
    for(alg in algorithm_list){
      
      pred <- predict_model(
        
        models[[alg]],
        
        base,
        
        alg
        
      )
      
      w <- sp_weights$Weight[
        sp_weights$Algorithm==alg
      ]
      
      ensemble_pred <- ensemble_pred +
        pred*w
      
    }
    
    df <- data.frame(
      
      Predictor=base[[v]],
      
      Suitability=ensemble_pred
      
    )
    
    p <- ggplot(
      
      df,
      
      aes(
        Predictor,
        Suitability
      )
      
    )+
      
      geom_line(
        linewidth=1.2,
        colour="steelblue"
      )+
      
      theme_bw(base_size=13)+
      
      labs(
        
        title=paste(sp,"-",v),
        
        y="Ensemble suitability"
        
      )
    
    ggsave(
      
      file.path(
        
        response_dir,
        
        paste0(
          
          gsub(" ","_",sp),
          
          "_",
          
          v,
          
          ".png"
          
        )
        
      ),
      
      p,
      
      width=5.5,
      
      height=4,
      
      dpi=500
      
    )
    
  }
  
}

message("Section 16 completed successfully.")




#============================================================
# 17. PUBLICATION FIGURE
# Ensemble suitability and prediction uncertainty
#============================================================

library(ggplot2)
library(tidyterra)
library(terra)
library(sf)
library(dplyr)
library(viridis)
library(patchwork)

#------------------------------------------------------------
# East African countries
#------------------------------------------------------------

countries <- rnaturalearth::ne_countries(

  country = c(

    "Sudan",
    "South Sudan",
    "Eritrea",
    "Djibouti",
    "Ethiopia",
    "Somalia",
    "Kenya",
    "Uganda",
    "Tanzania",
    "Rwanda",
    "Burundi"

  ),

  returnclass = "sf"

)

#------------------------------------------------------------
# Publication theme
#------------------------------------------------------------

theme_map <-

theme_void(base_size = 13)+

theme(

  plot.title = element_text(

    face = "italic",

    size = 15,

    hjust = 0.5

  ),

  legend.position = "right",

  legend.title = element_text(

    face = "bold",

    size = 11

  ),

  legend.text = element_text(size = 10),

  plot.margin = margin(3,3,3,3)

)

#------------------------------------------------------------
# Build one species panel
#------------------------------------------------------------

make_map <- function(species){

  occ_sf <-

    occ |>

    filter(species == !!species) |>

    st_as_sf(

      coords = c("lon","lat"),

      crs = 4326

    )

  ensemble <- rast(

    file.path(

      out_dir,

      "07_Ensembles",

      paste0(

        species,

        "_Ensemble.tif"

      )

    )

  )

  uncertainty <- rast(

    file.path(

      out_dir,

      "07_Ensembles",

      paste0(

        species,

        "_SD.tif"

      )

    )

  )

  title_expression <-

    if(species=="Phlebotomus orientalis"){

      expression(italic("Phlebotomus orientalis"))

    }else{

      expression(italic("Phlebotomus martini"))

    }

  #----------------------------------------------------------
  # Suitability
  #----------------------------------------------------------

  p1 <-

    ggplot() +

    geom_spatraster(

      data = ensemble

    ) +

    scale_fill_viridis_c(

      option = "inferno",

      limits = c(0,1),

      breaks = seq(0,1,0.2),

      name = "Suitability"

    ) +

    geom_sf(

      data = countries,

      fill = NA,

      colour = "grey65",

      linewidth = 0.25

    ) +

    geom_sf(

      data = occ_sf,

      shape = 21,

      fill = "black",

      colour = "white",

      stroke = 0.2,

      size = 1.0

    ) +

    labs(

      title = title_expression,

      tag = "A"

    ) +

    coord_sf(expand = FALSE) +

    theme_map

  #----------------------------------------------------------
  # Uncertainty
  #----------------------------------------------------------

  p2 <-

    ggplot() +

    geom_spatraster(

      data = uncertainty

    ) +

    scale_fill_viridis_c(

      option = "cividis",

      name = "SD"

    ) +

    geom_sf(

      data = countries,

      fill = NA,

      colour = "grey65",

      linewidth = 0.25

    ) +

    geom_sf(

      data = occ_sf,

      shape = 21,

      fill = "black",

      colour = "white",

      stroke = 0.2,

      size = 1.0

    ) +

    labs(

      tag = "B"

    ) +

    coord_sf(expand = FALSE) +

    theme_map

  p1 / p2

}

#------------------------------------------------------------
# Final figure
#------------------------------------------------------------

figure2 <-

  make_map("Phlebotomus orientalis") |

  make_map("Phlebotomus martini")

figure2 <-

figure2 +

plot_annotation(

  theme = theme(

    plot.margin = margin(5,5,5,5)

  )

)

ggsave(

  file.path(

    out_dir,

    "12_Publication_Figures",

    "Figure2_Ensemble_Maps.png"

  ),

  figure2,

  width = 14,

  height = 10,

  dpi = 800,

  bg = "white"

)

message("Figure 2 completed.")







#============================================================
# 18. FIGURE 3
# ENSEMBLE VARIABLE IMPORTANCE
#============================================================

message("Creating Figure 3...")

library(ggplot2)
library(dplyr)
library(patchwork)

importance <- read.csv(
  
  file.path(
    
    out_dir,
    
    "09_Variable_Importance",
    
    "Ensemble_Variable_Importance.csv"
    
  )
  
)

importance$Variable <- factor(
  
  importance$Variable,
  
  levels = rev(
    
    unique(
      
      importance$Variable
      
    )
    
  )
  
)

theme_bar <- theme_bw(base_size = 13) +
  
  theme(
    
    panel.grid.major.y = element_blank(),
    
    legend.position = "none",
    
    plot.title = element_text(face = "bold")
    
  )

#------------------------------------------------------------
# P. orientalis
#------------------------------------------------------------

p1 <- importance |>
  
  filter(
    
    Species == "Phlebotomus orientalis"
    
  ) |>
  
  ggplot(
    
    aes(
      
      Variable,
      
      Importance
      
    )
    
  ) +
  
  geom_col(
    
    fill = "#2C7FB8",
    
    width = 0.75
    
  ) +
  
  coord_flip() +
  
  labs(
    
    title = expression(
      
      italic("Phlebotomus orientalis")
      
    ),
    
    x = NULL,
    
    y = "Relative importance (%)"
    
  ) +
  
  theme_bar

#------------------------------------------------------------
# P. martini
#------------------------------------------------------------

p2 <- importance |>
  
  filter(
    
    Species == "Phlebotomus martini"
    
  ) |>
  
  ggplot(
    
    aes(
      
      Variable,
      
      Importance
      
    )
    
  ) +
  
  geom_col(
    
    fill = "#D95F0E",
    
    width = 0.75
    
  ) +
  
  coord_flip() +
  
  labs(
    
    title = expression(
      
      italic("Phlebotomus martini")
      
    ),
    
    x = NULL,
    
    y = "Relative importance (%)"
    
  ) +
  
  theme_bar

#------------------------------------------------------------
# Combine
#------------------------------------------------------------

figure3 <-
  
  p1 | p2 +
  
  plot_annotation(
    
    title = "Ensemble variable importance"
    
  )

ggsave(
  
  file.path(
    
    out_dir,
    
    "12_Publication_Figures",
    
    "Figure3_Ensemble_Variable_Importance.png"
    
  ),
  
  figure3,
  
  width = 11,
  
  height = 5,
  
  dpi = 600
  
)

message("Figure 3 completed.")







#============================================================
# 19. FIGURE 4
# ENSEMBLE RESPONSE CURVES
#============================================================

message("Creating Figure 4...")

library(ggplot2)
library(patchwork)
library(dplyr)

#------------------------------------------------------------
# Read importance results
#------------------------------------------------------------

importance <- read.csv(
  
  file.path(
    
    out_dir,
    
    "09_Variable_Importance",
    
    "Ensemble_Variable_Importance.csv"
    
  )
  
)

#------------------------------------------------------------
# Top four variables per species
#------------------------------------------------------------

top_vars <- importance |>
  
  group_by(Species) |>
  
  slice_max(
    
    Importance,
    
    n = 4
    
  ) |>
  
  pull(Variable) |>
  
  unique()

plots <- list()

counter <- 1

#------------------------------------------------------------
# Loop over species
#------------------------------------------------------------

for(sp in species_list){
  
  for(v in top_vars){
    
    file <- file.path(
      
      out_dir,
      
      "08_Response_Curves",
      
      paste0(
        
        gsub(" ","_",sp),
        
        "_",
        
        v,
        
        ".png"
        
      )
      
    )
    
    if(file.exists(file)){
      
      img <- png::readPNG(file)
      
      p <- ggplot() +
        
        annotation_raster(
          
          img,
          
          xmin=-Inf,
          
          xmax=Inf,
          
          ymin=-Inf,
          
          ymax=Inf
          
        ) +
        
        theme_void() +
        
        labs(
          
          title=paste(sp,"-",v)
          
        )
      
      plots[[counter]] <- p
      
      counter <- counter + 1
      
    }
    
  }
  
}

#------------------------------------------------------------
# Combine figure
#------------------------------------------------------------

figure4 <- wrap_plots(
  
  plots,
  
  ncol = 2
  
) +
  
  plot_annotation(
    
    title = "Ensemble response curves for the four most influential predictors"
    
  )

ggsave(
  
  file.path(
    
    out_dir,
    
    "12_Publication_Figures",
    
    "Figure4_Response_Curves.png"
    
  ),
  
  figure4,
  
  width = 10,
  
  height = 12,
  
  dpi = 600
  
)

message("Figure 4 completed.")








#============================================================
# TABLE 1. MODEL PERFORMANCE
#============================================================

message("Creating Table 1...")

performance_table <-
  
  algorithm_weights |>
  
  dplyr::select(
    
    Species,
    
    Algorithm,
    
    Mean_AUC,
    
    SD_AUC,
    
    Mean_TSS,
    
    SD_TSS,
    
    Weight
    
  ) |>
  
  dplyr::arrange(
    
    Species,
    
    desc(Mean_TSS)
    
  )

write.csv(
  
  performance_table,
  
  file.path(
    
    out_dir,
    
    "11_Model_Evaluation",
    
    "Table1_Model_Performance.csv"
    
  ),
  
  row.names = FALSE
  
)

print(performance_table)





#============================================================
# TABLE 2. FINAL ENSEMBLE SUMMARY
#============================================================

message("Creating Table 2...")

importance <- read.csv(
  
  file.path(
    
    out_dir,
    
    "09_Variable_Importance",
    
    "Ensemble_Variable_Importance.csv"
    
  )
  
)

top_predictor <-
  
  importance |>
  
  dplyr::group_by(
    
    Species
    
  ) |>
  
  dplyr::slice_max(
    
    Importance,
    
    n = 1
    
  ) |>
  
  dplyr::ungroup() |>
  
  dplyr::select(
    
    Species,
    
    TopPredictor = Variable,
    
    Importance
    
  )

best_algorithm <-
  
  algorithm_weights |>
  
  dplyr::group_by(
    
    Species
    
  ) |>
  
  dplyr::slice_max(
    
    Mean_TSS,
    
    n = 1
    
  ) |>
  
  dplyr::ungroup() |>
  
  dplyr::select(
    
    Species,
    
    BestAlgorithm = Algorithm,
    
    Mean_AUC,
    
    Mean_TSS
    
  )

ensemble_table <-
  
  dplyr::left_join(
    
    best_algorithm,
    
    top_predictor,
    
    by = "Species"
    
  ) |>
  
  dplyr::rename(
    
    Ensemble_AUC = Mean_AUC,
    
    Ensemble_TSS = Mean_TSS,
    
    TopPredictorImportance = Importance
    
  )

write.csv(
  
  ensemble_table,
  
  file.path(
    
    out_dir,
    
    "11_Model_Evaluation",
    
    "Table2_Ensemble_Summary.csv"
    
  ),
  
  row.names = FALSE
  
)

print(ensemble_table)
