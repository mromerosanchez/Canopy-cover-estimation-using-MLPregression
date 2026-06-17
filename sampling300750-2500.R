library(terra)

# Load Hansen canopy cover raster
cover <- rast("MexicoCity_ForestCover_Hansen.tif")

# Reclassify into bins
rcl <- matrix(c(
  0,20,1,
  20,40,2,
  40,60,3,
  60,80,4,
  80,100,5
), ncol=3, byrow=TRUE)

cover_class <- classify(cover, rcl)

# Stratified sampling
samples_1500 <- spatSample(
  cover_class,
  size = 1500,
  method = "stratified",
  as.points = TRUE,
  na.rm = TRUE
)

###Create layer stack

# =========================================================
# COMBINE SPECTRAL BANDS + VEGETATION INDICES INTO ONE STACK
# =========================================================

# Load library
library(terra)

# ---------------------------------------------------------
# 1. LOAD RASTER STACKS
# ---------------------------------------------------------

# Spectral reflectance bands
spectral <- rast("Landsat7_Reflectance_Only_2001.tif")

# Vegetation indices
indices <- rast("vegetation_indices_2001.tif")

# ---------------------------------------------------------
# 2. CHECK SPATIAL COMPATIBILITY
# ---------------------------------------------------------

geom_check <- compareGeom(spectral, indices)

# ---------------------------------------------------------
# 3. ALIGN STACKS IF NECESSARY
# ---------------------------------------------------------

if (!geom_check) {
  
  cat("Raster geometries do not match.\n")
  cat("Resampling vegetation indices to match spectral stack...\n")
  
  indices <- resample(
    indices,
    spectral,
    method = "bilinear"
  )
  
} else {
  
  cat("Raster geometries match.\n")
  
}

# ---------------------------------------------------------
# 4. COMBINE STACKS
# ---------------------------------------------------------

predictor_stack <- c(spectral, indices)

# ---------------------------------------------------------
# 5. CHECK RESULT
# ---------------------------------------------------------

print(predictor_stack)

cat("\nLayer names:\n")
print(names(predictor_stack))

cat("\nNumber of layers:\n")
print(nlyr(predictor_stack))

# ---------------------------------------------------------
# 6. SAVE COMBINED STACK
# ---------------------------------------------------------

writeRaster(
  predictor_stack,
  "predictor_stack_spectral_indices.tif",
  overwrite = TRUE
)

cat("\nCombined predictor stack saved successfully.\n")

# =========================================================
# OPTIONAL: EXTRACT VALUES FROM SAMPLE POINTS
# =========================================================

# Example:
# sample_points <- vect("training_points.shp")
# extracted_data <- extract(predictor_stack, sample_points)

# View extracted values
# head(extracted_data)

# =========================================================




###########################

library(terra)

# ---------------------------------------------------------
# 1. LOAD DATA
# ---------------------------------------------------------

predictors <- rast("predictor_stack_spectral_indices.tif")
hansen <- rast("MexicoCity_ForestCover_Hansen.tif")  # values 0–100

# Align Hansen to predictor stack if needed
if (!compareGeom(predictors, hansen, stopOnError = FALSE)) {
  hansen <- resample(hansen, predictors, method = "bilinear")
}

# ---------------------------------------------------------
# 2. CREATE CANOPY-COVER STRATA
# ---------------------------------------------------------

rcl <- matrix(c(
  0, 20, 1,
  20, 40, 2,
  40, 60, 3,
  60, 80, 4,
  80, 100, 5
), ncol = 3, byrow = TRUE)

hansen_class <- classify(hansen, rcl, include.lowest = TRUE)

names(hansen_class) <- "cover_class"

# ---------------------------------------------------------
# 3. FUNCTION TO CREATE STRATIFIED TRAINING DATA
# ---------------------------------------------------------

create_training_data <- function(n_samples, predictors, hansen, hansen_class) {
  
  set.seed(123)
  
  pts <- spatSample(
    hansen_class,
    size = n_samples,
    method = "stratified",
    as.points = TRUE,
    na.rm = TRUE
  )
  
  x_vals <- extract(predictors, pts)
  y_vals <- extract(hansen, pts)
  coords <- crds(pts)
  
  df <- data.frame(
    x = coords[,1],
    y = coords[,2],
    treecover = y_vals[,2],
    x_vals[, -1]
  )
  
  df <- na.omit(df)
  
  return(df)
}

# ---------------------------------------------------------
# 4. CREATE TRAINING DATASETS
# ---------------------------------------------------------

sample_sizes <- c(330, 750, 1500, 2500)

for (n in sample_sizes) {
  
  training_df <- create_training_data(
    n_samples = n,
    predictors = predictors,
    hansen = hansen,
    hansen_class = hansen_class
  )
  
  out_name <- paste0("training_data_", n, ".csv")
  write.csv(training_df, out_name, row.names = FALSE)
  
  cat("Saved:", out_name, "with", nrow(training_df), "valid samples\n")
}
