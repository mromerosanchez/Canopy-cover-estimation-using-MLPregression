# ============================================================
# REVISED CANOPY-COVER MODELLING WORKFLOW
# ============================================================
#
# Main improvements:
# 1. Exact stratified sampling: 300 samples per canopy class
# 2. Explicit predictor names and geometry checks
# 3. Complete-case sampling mask
# 4. Spatial five-fold cross-validation
# 5. Comparison of:
#       - Linear regression
#       - Random Forest
#       - Gradient boosting
#       - Multilayer Perceptron
# 6. Overall and class-specific error metrics
# 7. Regression-to-the-mean diagnostics
# 8. Combined uncertainty for two predicted maps
# 9. Change-threshold sensitivity analysis
#
# IMPORTANT:
# - Update all file paths and predictor names.
# - Check whether training imagery is from 2000 or 2001.
# - Do not report results until the complete workflow is rerun.
# ============================================================


# ============================================================
# 0. INSTALL AND LOAD PACKAGES
# ============================================================

required_packages <- c(
  "terra",
  "sf",
  "blockCV",
  "ranger",
  "xgboost",
  "torch"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before continuing: ",
    paste(missing_packages, collapse = ", ")
  )
}


install.packages("blockCV")
install.packages("sf")
install.packages("ranger")
install.packages("xgboost")

library(terra)
library(sf)
library(blockCV)
library(ranger)
library(xgboost)
library(torch)


# ============================================================
# 1. USER CONFIGURATION
# ============================================================

SEED <- 123L

# Equal stratified sample: 1000 points x 5 classes = 5000
N_PER_CLASS <- 1000L
N_CLASSES <- 5L
EXPECTED_TOTAL_N <- N_PER_CLASS * N_CLASSES

# Spatial block size in metres.
# Replace after evaluating spatial autocorrelation.
SPATIAL_BLOCK_SIZE_M <- 3000

# Cross-validation folds
N_FOLDS <- 5L

# MLP settings
MLP_EPOCHS <- 500L
MLP_LEARNING_RATE <- 0.001
MLP_WEIGHT_DECAY <- 1e-5
MLP_HIDDEN_1 <- 64L
MLP_HIDDEN_2 <- 32L

# Input rasters
SPECTRAL_FILE <- "Landsat7_Reflectance_Only_2001.tif"
INDICES_FILE  <- "vegetation_indices_2001_r.tif"
HANSEN_FILE   <- "MexicoCity_ForestCover_Hansen.tif"

# Output directory
OUTPUT_DIR <- "reviewer1_revised_analysis"

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

set.seed(SEED)
torch_manual_seed(SEED)


# ============================================================
# 2. LOAD DATA
# ============================================================

spectral <- rast(SPECTRAL_FILE)
indices  <- rast(INDICES_FILE)
hansen   <- rast(HANSEN_FILE)

if (nlyr(hansen) != 1L) {
  stop("The Hansen canopy-cover raster must contain one layer.")
}

cat("\nSpectral raster:\n")
print(spectral)

cat("\nVegetation-index raster:\n")
print(indices)

cat("\nHansen canopy-cover raster:\n")
print(hansen)


# ============================================================
# 3. ASSIGN EXPLICIT PREDICTOR NAMES
# ============================================================

# IMPORTANT:
# Replace these names if your layers are arranged differently.
# The manuscript must identify all predictors and their formulas.

spectral_names <- c(
  "blue",
  "green",
  "red",
  "nir",
  "swir1",
  "swir2"
)

# Replace index_4 and index_5 with the real predictor names.
index_names <- c(
  "ndvi",
  "savi",
  "evi",
  "Brightness",
  "Greenness"
)

if (nlyr(spectral) != length(spectral_names)) {
  stop(
    "Expected ", length(spectral_names),
    " spectral layers, but found ", nlyr(spectral), "."
  )
}

if (nlyr(indices) != length(index_names)) {
  stop(
    "Expected ", length(index_names),
    " vegetation-index layers, but found ", nlyr(indices), "."
  )
}

names(spectral) <- spectral_names
names(indices) <- index_names


# ============================================================
# 4. ALIGN VEGETATION INDICES TO SPECTRAL GRID
# ============================================================

indices_match <- terra::compareGeom(
  spectral,
  indices,
  crs = TRUE,
  ext = TRUE,
  rowcol = TRUE,
  res = TRUE,
  stopOnError = FALSE
)

if (!indices_match) {
  
  message("Vegetation-index geometry differs from spectral geometry.")
  
  if (!terra::same.crs(spectral, indices)) {
    
    message("Projecting vegetation indices to spectral CRS and grid.")
    
    indices <- terra::project(
      indices,
      spectral,
      method = "bilinear"
    )
    
  } else {
    
    message("Resampling vegetation indices to spectral grid.")
    
    indices <- terra::resample(
      indices,
      spectral,
      method = "bilinear"
    )
  }
  
} else {
  
  message("Vegetation-index geometry matches the spectral grid.")
}

# Verify alignment after project/resample
indices_match_after <- terra::compareGeom(
  spectral,
  indices,
  crs = TRUE,
  ext = TRUE,
  rowcol = TRUE,
  res = TRUE,
  stopOnError = FALSE
)

if (!indices_match_after) {
  stop(
    "Vegetation indices still do not match the spectral grid ",
    "after alignment."
  )
}

# Combine spectral bands and vegetation indices
predictors <- c(spectral, indices)

if (terra::nlyr(predictors) != 11L) {
  warning(
    "The combined stack contains ",
    terra::nlyr(predictors),
    " layers rather than the 11 predictors stated in the manuscript."
  )
}

if (anyDuplicated(names(predictors))) {
  stop("Predictor names must be unique.")
}

terra::writeRaster(
  predictors,
  filename = file.path(
    OUTPUT_DIR,
    "predictor_stack_training_year.tif"
  ),
  overwrite = TRUE
)

# ============================================================
# 5. ALIGN HANSEN CANOPY COVER
# ============================================================

hansen_match <- terra::compareGeom(
  predictors,
  hansen,
  crs = TRUE,
  ext = TRUE,
  rowcol = TRUE,
  res = TRUE,
  stopOnError = FALSE
)

# Check raster origins separately
origin_match <- isTRUE(
  all.equal(
    terra::origin(predictors),
    terra::origin(hansen),
    tolerance = 1e-08
  )
)

# Require both geometry and origin to match
full_hansen_match <- isTRUE(hansen_match) && origin_match

if (!full_hansen_match) {
  
  message("Hansen geometry differs from the predictor grid.")
  
  cat(
    "Same CRS:",
    terra::same.crs(predictors, hansen),
    "\n"
  )
  
  cat(
    "Predictor resolution:",
    terra::res(predictors),
    "\n"
  )
  
  cat(
    "Hansen resolution:",
    terra::res(hansen),
    "\n"
  )
  
  cat(
    "Predictor origin:",
    terra::origin(predictors),
    "\n"
  )
  
  cat(
    "Hansen origin:",
    terra::origin(hansen),
    "\n"
  )
  
  if (!terra::same.crs(predictors, hansen)) {
    
    message(
      "Projecting Hansen canopy cover to the predictor CRS and grid."
    )
    
    hansen_aligned <- terra::project(
      hansen,
      predictors,
      method = "bilinear"
    )
    
  } else {
    
    message(
      "Resampling Hansen canopy cover to the predictor grid."
    )
    
    hansen_aligned <- terra::resample(
      hansen,
      predictors,
      method = "bilinear"
    )
  }
  
} else {
  
  message("Hansen canopy cover already matches the predictor grid.")
  
  hansen_aligned <- hansen
}

# Assign a clear layer name
names(hansen_aligned) <- "canopy_cover"

# Restrict canopy-cover values to the valid range
hansen_aligned <- terra::clamp(
  hansen_aligned,
  lower = 0,
  upper = 100,
  values = TRUE
)

# Verify geometry after alignment
final_geom_match <- terra::compareGeom(
  predictors,
  hansen_aligned,
  crs = TRUE,
  ext = TRUE,
  rowcol = TRUE,
  res = TRUE,
  stopOnError = FALSE
)

final_origin_match <- isTRUE(
  all.equal(
    terra::origin(predictors),
    terra::origin(hansen_aligned),
    tolerance = 1e-08
  )
)

if (!isTRUE(final_geom_match) || !final_origin_match) {
  stop(
    "Hansen canopy cover still does not match the predictor grid ",
    "after alignment."
  )
}

message("Hansen canopy cover was aligned successfully.")

cat("\nAligned Hansen canopy-cover range:\n")

print(
  terra::global(
    hansen_aligned,
    fun = "range",
    na.rm = TRUE
  )
)


# ============================================================
# 6. CREATE FIVE CANOPY-COVER CLASSES
# ============================================================

# Classes:
# 1 = [0, 20)
# 2 = [20, 40)
# 3 = [40, 60)
# 4 = [60, 80)
# 5 = [80, 100]

rcl <- matrix(
  c(
    0,  20, 1,
    20,  40, 2,
    40,  60, 3,
    60,  80, 4,
    80, 101, 5
  ),
  ncol = 3,
  byrow = TRUE
)

hansen_class <- classify(
  hansen_aligned,
  rcl = rcl,
  right = FALSE,
  include.lowest = TRUE,
  others = NA
)

names(hansen_class) <- "cover_class"

cat("\nAvailable cells by canopy class:\n")
print(freq(hansen_class))


# ============================================================
# 7. CREATE COMPLETE-CASE SAMPLING MASK
# ============================================================

# Number of non-missing predictors in each cell
non_missing_count <- app(
  predictors,
  fun = function(x) sum(!is.na(x))
)

complete_predictor_mask <- ifel(
  non_missing_count == nlyr(predictors),
  1,
  NA
)

valid_sampling_mask <- ifel(
  !is.na(hansen_aligned) &
    !is.na(hansen_class) &
    !is.na(complete_predictor_mask),
  1,
  NA
)

sampling_classes <- mask(
  hansen_class,
  valid_sampling_mask
)

cat("\nValid cells available for sampling:\n")
print(freq(sampling_classes))


# ============================================================
# 8. SAMPLE EXACTLY 1000 POINTS PER CLASS
# ============================================================

sample_one_class <- function(
    class_raster,
    class_value,
    sample_n,
    seed
) {
  
  set.seed(seed + class_value)
  
  class_mask <- ifel(
    class_raster == class_value,
    class_value,
    NA
  )
  
  available_n <- global(
    !is.na(class_mask),
    sum,
    na.rm = TRUE
  )[1, 1]
  
  if (available_n < sample_n) {
    stop(
      "Class ", class_value,
      " contains only ", available_n,
      " valid cells; requested ", sample_n, "."
    )
  }
  
  points <- spatSample(
    class_mask,
    size = sample_n,
    method = "random",
    as.points = TRUE,
    values = TRUE,
    cells = TRUE,
    na.rm = TRUE,
    replace = FALSE
  )
  
  names(points)[1] <- "cover_class"
  
  points
}

point_list <- lapply(
  X = 1:N_CLASSES,
  FUN = function(k) {
    sample_one_class(
      class_raster = sampling_classes,
      class_value = k,
      sample_n = N_PER_CLASS,
      seed = SEED
    )
  }
)

sample_points <- do.call(rbind, point_list)

if (nrow(sample_points) != EXPECTED_TOTAL_N) {
  stop(
    "Expected ", EXPECTED_TOTAL_N,
    " points, but obtained ", nrow(sample_points), "."
  )
}

cat("\nSampled points by class:\n")
print(table(values(sample_points)$cover_class))



# ============================================================
# 9. EXTRACT CONTINUOUS RESPONSE, CLASS, AND PREDICTORS
# ============================================================

# Coordinates of sampled points
xy <- terra::crds(
  sample_points,
  df = FALSE
)

# Extract continuous canopy-cover response
response_values <- terra::extract(
  hansen_aligned,
  sample_points,
  ID = FALSE
)

# Extract canopy-cover class directly from the class raster
class_values <- terra::extract(
  hansen_class,
  sample_points,
  ID = FALSE
)

# Extract all predictor values
predictor_values <- terra::extract(
  predictors,
  sample_points,
  ID = FALSE
)

# Check dimensions before constructing the data frame
cat("Number of sample points:", nrow(sample_points), "\n")
cat("Coordinate rows:", nrow(xy), "\n")
cat("Response rows:", nrow(response_values), "\n")
cat("Class rows:", nrow(class_values), "\n")
cat("Predictor rows:", nrow(predictor_values), "\n")

# All extracted objects must have the same number of rows
expected_n <- nrow(sample_points)

row_counts <- c(
  points = expected_n,
  coordinates = nrow(xy),
  response = nrow(response_values),
  classes = nrow(class_values),
  predictors = nrow(predictor_values)
)

print(row_counts)

if (length(unique(row_counts)) != 1L) {
  stop(
    "The number of rows differs among sampled points and extracted values: ",
    paste(
      names(row_counts),
      row_counts,
      sep = "=",
      collapse = ", "
    )
  )
}

# Construct training data
training_df <- data.frame(
  point_id = seq_len(expected_n),
  
  cell_id = terra::cellFromXY(
    hansen_aligned,
    xy
  ),
  
  x = xy[, 1],
  y = xy[, 2],
  
  cover_class = class_values[, 1],
  
  canopy_cover = response_values[, 1],
  
  predictor_values,
  
  check.names = FALSE
)

# Inspect result
cat("\nTraining-data dimensions before cleaning:\n")
print(dim(training_df))

cat("\nCanopy classes before cleaning:\n")
print(table(
  training_df$cover_class,
  useNA = "ifany"
))

cat("\nCanopy-cover summary:\n")
print(summary(training_df$canopy_cover))

# ============================================================
# 9B. QUALITY CONTROL OF EXTRACTED TRAINING DATA
# ============================================================

rows_before <- nrow(training_df)

# Remove rows containing NA in any required variable
training_df <- training_df[
  complete.cases(training_df),
  ,
  drop = FALSE
]

rows_after_complete_cases <- nrow(training_df)

# Remove duplicate sampled raster cells
duplicate_cells <- duplicated(training_df$cell_id)

if (any(duplicate_cells)) {
  
  warning(
    sum(duplicate_cells),
    " duplicate raster cells were detected and removed."
  )
  
  training_df <- training_df[
    !duplicate_cells,
    ,
    drop = FALSE
  ]
}

cat("\nRows before cleaning:", rows_before, "\n")
cat(
  "Rows after removing incomplete cases:",
  rows_after_complete_cases,
  "\n"
)
cat(
  "Rows after removing duplicate cells:",
  nrow(training_df),
  "\n"
)

cat("\nFinal samples per canopy class:\n")
print(table(
  training_df$cover_class,
  useNA = "ifany"
))

write.csv(
  training_df,
  file.path(
    OUTPUT_DIR,
    "training_data_stratified_5000.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 10. DATA-QUALITY SUMMARY
# ============================================================

predictor_names <- names(predictors)

numeric_variables <- c(
  "canopy_cover",
  predictor_names
)

quality_summary <- data.frame(
  variable = numeric_variables,
  minimum = vapply(
    training_df[numeric_variables],
    min,
    numeric(1),
    na.rm = TRUE
  ),
  maximum = vapply(
    training_df[numeric_variables],
    max,
    numeric(1),
    na.rm = TRUE
  ),
  mean = vapply(
    training_df[numeric_variables],
    mean,
    numeric(1),
    na.rm = TRUE
  ),
  standard_deviation = vapply(
    training_df[numeric_variables],
    sd,
    numeric(1),
    na.rm = TRUE
  )
)

write.csv(
  quality_summary,
  file.path(
    OUTPUT_DIR,
    "training_data_quality_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 11. CREATE SPATIAL CROSS-VALIDATION FOLDS
# ============================================================

# ============================================================
# 11. CREATE SPATIAL CROSS-VALIDATION FOLDS
# ============================================================

library(sf)
library(blockCV)

# Create sf points using the predictor CRS
sample_sf <- sf::st_as_sf(
  training_df,
  coords = c("x", "y"),
  crs = terra::crs(predictors)
)

# Confirm that the source CRS exists
if (is.na(sf::st_crs(sample_sf))) {
  stop(
    "The sample points have no valid CRS. ",
    "Check terra::crs(predictors)."
  )
}

cat("Original CRS:\n")
print(sf::st_crs(sample_sf))

# Transform geographic coordinates to a metric CRS
if (isTRUE(sf::st_is_longlat(sample_sf))) {
  
  message(
    "Sample points use longitude/latitude coordinates. ",
    "Transforming them to WGS 84 / UTM zone 14N (EPSG:32614)."
  )
  
  sample_sf_metric <- sf::st_transform(
    sample_sf,
    crs = 32614
  )
  
} else {
  
  message("Sample points already use a projected CRS.")
  
  sample_sf_metric <- sample_sf
}

# Verify that the transformed CRS is projected
if (isTRUE(sf::st_is_longlat(sample_sf_metric))) {
  stop(
    "The points are still in geographic coordinates after transformation."
  )
}

cat("\nMetric CRS:\n")
print(sf::st_crs(sample_sf_metric))

cat("\nCoordinate range after transformation:\n")
print(sf::st_bbox(sample_sf_metric))


set.seed(SEED)

spatial_blocks <- cv_spatial(
  x = sample_sf,
  column = "canopy_cover",
  r = predictors[[1]],
  k = N_FOLDS,
  size = SPATIAL_BLOCK_SIZE_M,
  selection = "random",
  iteration = 100,
  seed = SEED,
  plot = FALSE,
  progress = FALSE
)

# blockCV output normally includes folds_ids
if (is.null(spatial_blocks$folds_ids)) {
  stop(
    "Could not find 'folds_ids' in the blockCV output. ",
    "Run names(spatial_blocks) and inspect your installed blockCV version."
  )
}

training_df$spatial_fold <- spatial_blocks$folds_ids

cat("\nSpatial fold counts:\n")
print(table(training_df$spatial_fold))

if (length(unique(training_df$spatial_fold)) != N_FOLDS) {
  stop("The expected number of spatial folds was not created.")
}

write.csv(
  training_df,
  file.path(
    OUTPUT_DIR,
    "training_data_with_spatial_folds.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. EVALUATION METRICS
# ============================================================

regression_metrics <- function(observed, predicted) {
  
  keep <- complete.cases(observed, predicted)
  
  observed <- observed[keep]
  predicted <- predicted[keep]
  
  if (length(observed) < 3L) {
    return(
      data.frame(
        n = length(observed),
        RMSE = NA_real_,
        MAE = NA_real_,
        Bias = NA_real_,
        Pearson_r = NA_real_,
        R2 = NA_real_,
        calibration_intercept = NA_real_,
        calibration_slope = NA_real_
      )
    )
  }
  
  residual_error <- predicted - observed
  calibration_fit <- lm(predicted ~ observed)
  
  r_value <- suppressWarnings(
    cor(
      observed,
      predicted,
      method = "pearson"
    )
  )
  
  data.frame(
    n = length(observed),
    RMSE = sqrt(mean(residual_error^2)),
    MAE = mean(abs(residual_error)),
    Bias = mean(residual_error),
    Pearson_r = r_value,
    R2 = r_value^2,
    calibration_intercept =
      unname(coef(calibration_fit)[1]),
    calibration_slope =
      unname(coef(calibration_fit)[2])
  )
}


clamp_predictions <- function(x) {
  pmin(100, pmax(0, x))
}


# ============================================================
# 13. MLP FUNCTIONS USING TORCH
# ============================================================

mlp_module <- nn_module(
  classname = "CanopyMLP",
  
  initialize = function(
    input_size,
    hidden_1 = 64,
    hidden_2 = 32
  ) {
    
    self$network <- nn_sequential(
      nn_linear(input_size, hidden_1),
      nn_relu(),
      nn_linear(hidden_1, hidden_2),
      nn_relu(),
      nn_linear(hidden_2, 1)
    )
  },
  
  forward = function(x) {
    self$network(x)
  }
)


fit_mlp <- function(
    train_x,
    train_y,
    hidden_1 = 64L,
    hidden_2 = 32L,
    epochs = 500L,
    learning_rate = 0.001,
    weight_decay = 1e-5,
    seed = 123L
) {
  
  torch_manual_seed(seed)
  
  train_x <- as.matrix(train_x)
  train_y <- as.numeric(train_y)
  
  predictor_means <- colMeans(train_x)
  predictor_sds <- apply(train_x, 2, sd)
  
  predictor_sds[
    is.na(predictor_sds) | predictor_sds == 0
  ] <- 1
  
  x_scaled <- sweep(
    train_x,
    2,
    predictor_means,
    FUN = "-"
  )
  
  x_scaled <- sweep(
    x_scaled,
    2,
    predictor_sds,
    FUN = "/"
  )
  
  response_mean <- mean(train_y)
  response_sd <- sd(train_y)
  
  if (is.na(response_sd) || response_sd == 0) {
    stop("Training response has zero variance.")
  }
  
  y_scaled <- (
    train_y - response_mean
  ) / response_sd
  
  x_tensor <- torch_tensor(
    x_scaled,
    dtype = torch_float()
  )
  
  y_tensor <- torch_tensor(
    matrix(y_scaled, ncol = 1),
    dtype = torch_float()
  )
  
  model <- mlp_module(
    input_size = ncol(train_x),
    hidden_1 = hidden_1,
    hidden_2 = hidden_2
  )
  
  optimizer <- optim_adam(
    model$parameters,
    lr = learning_rate,
    weight_decay = weight_decay
  )
  
  loss_function <- nn_mse_loss()
  
  loss_history <- numeric(epochs)
  
  model$train()
  
  for (epoch in seq_len(epochs)) {
    
    optimizer$zero_grad()
    
    predictions <- model(x_tensor)
    
    loss <- loss_function(
      predictions,
      y_tensor
    )
    
    loss$backward()
    optimizer$step()
    
    loss_history[epoch] <- loss$item()
  }
  
  list(
    model = model,
    predictor_means = predictor_means,
    predictor_sds = predictor_sds,
    response_mean = response_mean,
    response_sd = response_sd,
    loss_history = loss_history
  )
}


predict_mlp <- function(fitted_object, new_x) {
  
  new_x <- as.matrix(new_x)
  
  x_scaled <- sweep(
    new_x,
    2,
    fitted_object$predictor_means,
    FUN = "-"
  )
  
  x_scaled <- sweep(
    x_scaled,
    2,
    fitted_object$predictor_sds,
    FUN = "/"
  )
  
  x_tensor <- torch_tensor(
    x_scaled,
    dtype = torch_float()
  )
  
  fitted_object$model$eval()
  
  with_no_grad({
    
    scaled_predictions <-
      as.numeric(
        fitted_object$model(x_tensor)
      )
  })
  
  predictions <-
    scaled_predictions *
    fitted_object$response_sd +
    fitted_object$response_mean
  
  clamp_predictions(predictions)
}


# ============================================================
# 14. SPATIAL CROSS-VALIDATION OF ALL MODELS
# ============================================================

model_variables <- c(
  "canopy_cover",
  predictor_names
)

fold_ids <- sort(
  unique(training_df$spatial_fold)
)

prediction_records <- list()
record_index <- 1L

for (fold in fold_ids) {
  
  message(
    "Processing spatial fold ",
    fold,
    " of ",
    length(fold_ids)
  )
  
  train_data <- training_df[
    training_df$spatial_fold != fold,
  ]
  
  test_data <- training_df[
    training_df$spatial_fold == fold,
  ]
  
  train_x <- train_data[, predictor_names, drop = FALSE]
  test_x <- test_data[, predictor_names, drop = FALSE]
  
  train_y <- train_data$canopy_cover
  test_y <- test_data$canopy_cover
  
  
  # ----------------------------------------------------------
  # 14A. MULTIPLE LINEAR REGRESSION
  # ----------------------------------------------------------
  
  linear_formula <- as.formula(
    paste(
      "canopy_cover ~",
      paste(predictor_names, collapse = " + ")
    )
  )
  
  linear_model <- lm(
    linear_formula,
    data = train_data
  )
  
  linear_predictions <- predict(
    linear_model,
    newdata = test_data
  )
  
  prediction_records[[record_index]] <- data.frame(
    point_id = test_data$point_id,
    fold = fold,
    model = "Linear regression",
    observed = test_y,
    predicted = clamp_predictions(
      linear_predictions
    ),
    cover_class = test_data$cover_class
  )
  
  record_index <- record_index + 1L
  
  
  # ----------------------------------------------------------
  # 14B. RANDOM FOREST
  # ----------------------------------------------------------
  
  rf_model <- ranger(
    formula = linear_formula,
    data = train_data[, model_variables],
    num.trees = 500,
    mtry = max(
      1L,
      floor(sqrt(length(predictor_names)))
    ),
    min.node.size = 5,
    importance = "permutation",
    seed = SEED + fold
  )
  
  rf_predictions <- predict(
    rf_model,
    data = test_data[, predictor_names, drop = FALSE]
  )$predictions
  
  prediction_records[[record_index]] <- data.frame(
    point_id = test_data$point_id,
    fold = fold,
    model = "Random Forest",
    observed = test_y,
    predicted = clamp_predictions(
      rf_predictions
    ),
    cover_class = test_data$cover_class
  )
  
  record_index <- record_index + 1L
  
  
  # ----------------------------------------------------------
  # 14C. GRADIENT BOOSTING
  # ----------------------------------------------------------
  
  xgb_train <- xgb.DMatrix(
    data = as.matrix(train_x),
    label = train_y
  )
  
  xgb_test <- xgb.DMatrix(
    data = as.matrix(test_x)
  )
  
  xgb_model <- xgb.train(
    params = list(
      objective = "reg:squarederror",
      eval_metric = "rmse",
      eta = 0.03,
      max_depth = 4,
      min_child_weight = 5,
      subsample = 0.8,
      colsample_bytree = 0.8,
      seed = SEED + fold
    ),
    data = xgb_train,
    nrounds = 500,
    verbose = 0
  )
  
  xgb_predictions <- predict(
    xgb_model,
    xgb_test
  )
  
  prediction_records[[record_index]] <- data.frame(
    point_id = test_data$point_id,
    fold = fold,
    model = "Gradient boosting",
    observed = test_y,
    predicted = clamp_predictions(
      xgb_predictions
    ),
    cover_class = test_data$cover_class
  )
  
  record_index <- record_index + 1L
  
  
  # ----------------------------------------------------------
  # 14D. MULTILAYER PERCEPTRON
  # ----------------------------------------------------------
  
  mlp_fit <- fit_mlp(
    train_x = train_x,
    train_y = train_y,
    hidden_1 = MLP_HIDDEN_1,
    hidden_2 = MLP_HIDDEN_2,
    epochs = MLP_EPOCHS,
    learning_rate = MLP_LEARNING_RATE,
    weight_decay = MLP_WEIGHT_DECAY,
    seed = SEED + fold
  )
  
  mlp_predictions <- predict_mlp(
    fitted_object = mlp_fit,
    new_x = test_x
  )
  
  prediction_records[[record_index]] <- data.frame(
    point_id = test_data$point_id,
    fold = fold,
    model = "MLP",
    observed = test_y,
    predicted = mlp_predictions,
    cover_class = test_data$cover_class
  )
  
  record_index <- record_index + 1L
}

cv_predictions <- do.call(
  rbind,
  prediction_records
)

write.csv(
  cv_predictions,
  file.path(
    OUTPUT_DIR,
    "spatial_cv_all_predictions.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 15. OVERALL MODEL PERFORMANCE
# ============================================================

model_split <- split(
  cv_predictions,
  cv_predictions$model
)

overall_results <- do.call(
  rbind,
  lapply(
    model_split,
    function(d) {
      
      metrics <- regression_metrics(
        observed = d$observed,
        predicted = d$predicted
      )
      
      data.frame(
        model = unique(d$model),
        metrics,
        row.names = NULL
      )
    }
  )
)

rownames(overall_results) <- NULL

overall_results <- overall_results[
  order(overall_results$RMSE),
]

cat("\nOverall spatial cross-validation results:\n")
print(overall_results)

write.csv(
  overall_results,
  file.path(
    OUTPUT_DIR,
    "model_comparison_overall.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16. PERFORMANCE BY SPATIAL FOLD
# ============================================================

fold_model_groups <- split(
  cv_predictions,
  interaction(
    cv_predictions$model,
    cv_predictions$fold,
    drop = TRUE
  )
)

fold_results <- do.call(
  rbind,
  lapply(
    fold_model_groups,
    function(d) {
      
      metrics <- regression_metrics(
        observed = d$observed,
        predicted = d$predicted
      )
      
      data.frame(
        model = unique(d$model),
        fold = unique(d$fold),
        metrics,
        row.names = NULL
      )
    }
  )
)

rownames(fold_results) <- NULL

write.csv(
  fold_results,
  file.path(
    OUTPUT_DIR,
    "model_comparison_by_fold.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. CLASS-SPECIFIC PERFORMANCE
# ============================================================

class_model_groups <- split(
  cv_predictions,
  interaction(
    cv_predictions$model,
    cv_predictions$cover_class,
    drop = TRUE
  )
)

class_results <- do.call(
  rbind,
  lapply(
    class_model_groups,
    function(d) {
      
      metrics <- regression_metrics(
        observed = d$observed,
        predicted = d$predicted
      )
      
      data.frame(
        model = unique(d$model),
        cover_class = unique(d$cover_class),
        mean_observed = mean(d$observed),
        mean_predicted = mean(d$predicted),
        metrics,
        row.names = NULL
      )
    }
  )
)

rownames(class_results) <- NULL

class_results$class_mean_bias <-
  class_results$mean_predicted -
  class_results$mean_observed

class_results <- class_results[
  order(
    class_results$model,
    class_results$cover_class
  ),
]

cat("\nClass-specific results:\n")
print(class_results)

write.csv(
  class_results,
  file.path(
    OUTPUT_DIR,
    "model_comparison_by_canopy_class.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 18. GRADIENT-COMPRESSION DIAGNOSTICS
# ============================================================

gradient_results <- do.call(
  rbind,
  lapply(
    model_split,
    function(d) {
      
      fit <- lm(
        predicted ~ observed,
        data = d
      )
      
      data.frame(
        model = unique(d$model),
        intercept = unname(coef(fit)[1]),
        slope = unname(coef(fit)[2]),
        predicted_at_observed_0 =
          unname(coef(fit)[1]),
        predicted_at_observed_100 =
          unname(coef(fit)[1] +
                   100 * coef(fit)[2])
      )
    }
  )
)

rownames(gradient_results) <- NULL

cat("\nGradient-compression diagnostics:\n")
print(gradient_results)

write.csv(
  gradient_results,
  file.path(
    OUTPUT_DIR,
    "gradient_compression_diagnostics.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 19. BASIC DIAGNOSTIC PLOTS
# ============================================================

pdf(
  file.path(
    OUTPUT_DIR,
    "model_diagnostic_plots.pdf"
  ),
  width = 10,
  height = 8
)

models_to_plot <- unique(
  cv_predictions$model
)

par(mfrow = c(2, 2))

for (model_name in models_to_plot) {
  
  d <- cv_predictions[
    cv_predictions$model == model_name,
  ]
  
  plot(
    d$observed,
    d$predicted,
    pch = 16,
    cex = 0.55,
    col = rgb(0, 0, 0, 0.3),
    xlim = c(0, 100),
    ylim = c(0, 100),
    xlab = "Reference canopy cover (%)",
    ylab = "Predicted canopy cover (%)",
    main = model_name
  )
  
  abline(
    a = 0,
    b = 1,
    col = "red",
    lwd = 2,
    lty = 2
  )
  
  abline(
    lm(predicted ~ observed, data = d),
    col = "blue",
    lwd = 2
  )
}

dev.off()


# ============================================================
# 20. SELECT MODEL FOR FINAL MAPPING
# ============================================================

# This automatically selects the model with the lowest
# spatial-cross-validation RMSE.
#
# You may instead select a model based on multiple criteria,
# but the choice must be explained in the manuscript.

selected_model_name <- overall_results$model[
  which.min(overall_results$RMSE)
]

selected_model_rmse <- overall_results$RMSE[
  which.min(overall_results$RMSE)
]

cat(
  "\nSelected model:",
  selected_model_name,
  "\nSpatial-CV RMSE:",
  selected_model_rmse,
  "\n"
)


# ============================================================
# 21. FIT FINAL MODEL TO ALL TRAINING DATA
# ============================================================

full_x <- training_df[
  ,
  predictor_names,
  drop = FALSE
]

full_y <- training_df$canopy_cover

full_formula <- as.formula(
  paste(
    "canopy_cover ~",
    paste(predictor_names, collapse = " + ")
  )
)

if (selected_model_name == "Linear regression") {
  
  final_model <- lm(
    full_formula,
    data = training_df
  )
  
} else if (selected_model_name == "Random Forest") {
  
  final_model <- ranger(
    formula = full_formula,
    data = training_df[, model_variables],
    num.trees = 500,
    mtry = max(
      1L,
      floor(sqrt(length(predictor_names)))
    ),
    min.node.size = 5,
    importance = "permutation",
    seed = SEED
  )
  
} else if (selected_model_name == "Gradient boosting") {
  
  final_model <- xgb.train(
    params = list(
      objective = "reg:squarederror",
      eval_metric = "rmse",
      eta = 0.03,
      max_depth = 4,
      min_child_weight = 5,
      subsample = 0.8,
      colsample_bytree = 0.8,
      seed = SEED
    ),
    data = xgb.DMatrix(
      as.matrix(full_x),
      label = full_y
    ),
    nrounds = 500,
    verbose = 0
  )
  
} else if (selected_model_name == "MLP") {
  
  final_model <- fit_mlp(
    train_x = full_x,
    train_y = full_y,
    hidden_1 = MLP_HIDDEN_1,
    hidden_2 = MLP_HIDDEN_2,
    epochs = MLP_EPOCHS,
    learning_rate = MLP_LEARNING_RATE,
    weight_decay = MLP_WEIGHT_DECAY,
    seed = SEED
  )
  
} else {
  
  stop("The selected model name was not recognized.")
}

saveRDS(
  final_model,
  file.path(
    OUTPUT_DIR,
    "selected_final_model.rds"
  )
)


# ============================================================
# 22. FUNCTIONS FOR ANNUAL RASTER PREDICTION
# ============================================================

predict_final_model_matrix <- function(
    model_name,
    fitted_model,
    newdata
) {
  
  newdata <- as.data.frame(newdata)
  
  if (model_name == "Linear regression") {
    
    output <- predict(
      fitted_model,
      newdata = newdata
    )
    
  } else if (model_name == "Random Forest") {
    
    output <- predict(
      fitted_model,
      data = newdata
    )$predictions
    
  } else if (model_name == "Gradient boosting") {
    
    output <- predict(
      fitted_model,
      xgb.DMatrix(
        as.matrix(newdata)
      )
    )
    
  } else if (model_name == "MLP") {
    
    output <- predict_mlp(
      fitted_object = fitted_model,
      new_x = newdata
    )
    
  } else {
    
    stop("Unsupported model.")
  }
  
  clamp_predictions(output)
}


# Wrapper compatible with terra::app()
predict_raster_block <- function(
    values_matrix,
    model_name,
    fitted_model,
    predictor_names
) {
  
  values_matrix <- as.matrix(values_matrix)
  
  output <- rep(
    NA_real_,
    nrow(values_matrix)
  )
  
  complete_rows <- complete.cases(
    values_matrix
  )
  
  if (any(complete_rows)) {
    
    newdata <- as.data.frame(
      values_matrix[
        complete_rows,
        ,
        drop = FALSE
      ]
    )
    
    names(newdata) <- predictor_names
    
    output[complete_rows] <-
      predict_final_model_matrix(
        model_name = model_name,
        fitted_model = fitted_model,
        newdata = newdata
      )
  }
  
  output
}


# ============================================================
# 23. PREPARE ANNUAL PREDICTOR STACKS
# ============================================================

reflectance_files <- c(
  "1994" = "Landsat5_Reflectance_Only_1994.tif",
  "2003" = "Landsat7_Reflectance_Only_2003.tif",
  "2014" = "Landsat8_Reflectance_Only_2014.tif",
  "2024" = "Landsat8_Reflectance_Only_2024.tif"
)

indices_files <- c(
  "1994" = "vegetation_indices_1994_r.tif",
  "2003" = "vegetation_indices_2003_r.tif",
  "2014" = "vegetation_indices_2014_r.tif",
  "2024" = "vegetation_indices_2024_r.tif"
)

annual_predictions <- list()

for (year_name in names(reflectance_files)) {
  
  message("Processing ", year_name)
  
  spectral <- terra::rast(reflectance_files[[year_name]])
  indices  <- terra::rast(indices_files[[year_name]])
  
  #----------------------------------------------------------
  # Align vegetation indices to reflectance
  #----------------------------------------------------------
  
  geom_ok <- terra::compareGeom(
    spectral,
    indices,
    crs = TRUE,
    ext = TRUE,
    rowcol = TRUE,
    res = TRUE,
    stopOnError = FALSE
  )
  
  if (!isTRUE(geom_ok)) {
    
    if (!terra::same.crs(spectral, indices)) {
      
      indices <- terra::project(
        indices,
        spectral,
        method = "bilinear"
      )
      
    } else {
      
      indices <- terra::resample(
        indices,
        spectral,
        method = "bilinear"
      )
      
    }
    
  }
  
  #----------------------------------------------------------
  # Build predictor stack
  #----------------------------------------------------------
  
  annual_stack <- c(spectral, indices)
  
  ## Use exactly the same names as the training stack
  names(annual_stack) <- predictor_names
  
  #----------------------------------------------------------
  # Verify compatibility with training raster
  #----------------------------------------------------------
  
  geom_ok <- terra::compareGeom(
    predictors,
    annual_stack,
    crs = TRUE,
    ext = TRUE,
    rowcol = TRUE,
    res = TRUE,
    stopOnError = FALSE
  )
  
  if (!isTRUE(geom_ok)) {
    
    annual_stack <- terra::resample(
      annual_stack,
      predictors,
      method = "bilinear"
    )
    
  }
  
  #----------------------------------------------------------
  # Prediction
  #----------------------------------------------------------
  
  prediction_raster <- terra::app(
    
    annual_stack,
    
    fun = function(x){
      
      predict_raster_block(
        
        values_matrix = x,
        
        model_name = selected_model_name,
        
        fitted_model = final_model,
        
        predictor_names = predictor_names
        
      )
      
    },
    
    filename = file.path(
      OUTPUT_DIR,
      paste0(
        "predicted_canopy_",
        year_name,
        ".tif"
      )
    ),
    
    overwrite = TRUE
    
  )
  
  names(prediction_raster) <- paste0("canopy_", year_name)
  
  annual_predictions[[year_name]] <- prediction_raster
  
}


# ============================================================
# 24. CHANGE-UNCERTAINTY FUNCTIONS
# ============================================================

# The error of a difference between two predictions is
# approximated as:
#
# sqrt(RMSE_t1^2 + RMSE_t2^2)
#
# If both years use the same model RMSE:
#
# sqrt(2) * RMSE

combined_change_rmse <- function(
    rmse_earlier,
    rmse_later
) {
  
  sqrt(
    rmse_earlier^2 +
      rmse_later^2
  )
}


classify_canopy_change <- function(
    earlier_raster,
    later_raster,
    threshold,
    output_prefix = NULL
) {
  
  if (!compareGeom(
    earlier_raster,
    later_raster,
    stopOnError = FALSE
  )) {
    stop("The two annual prediction rasters are not aligned.")
  }
  
  difference_raster <-
    later_raster -
    earlier_raster
  
  names(difference_raster) <-
    "canopy_difference"
  
  # -1 = detectable loss
  #  0 = no detectable change / uncertain
  #  1 = detectable gain
  
  change_raster <- ifel(
    difference_raster < -threshold,
    -1,
    ifel(
      difference_raster > threshold,
      1,
      0
    )
  )
  
  names(change_raster) <- "change_class"
  
  if (!is.null(output_prefix)) {
    
    writeRaster(
      difference_raster,
      paste0(
        output_prefix,
        "_difference.tif"
      ),
      overwrite = TRUE
    )
    
    writeRaster(
      change_raster,
      paste0(
        output_prefix,
        "_classes.tif"
      ),
      overwrite = TRUE
    )
  }
  
  list(
    difference = difference_raster,
    classes = change_raster,
    threshold = threshold
  )
}


summarize_change_area <- function(
    change_raster,
    period_name,
    threshold_name,
    threshold_value
) {
  
  area_raster <- cellSize(
    change_raster,
    unit = "km"
  )
  
  categories <- c(-1, 0, 1)
  
  category_labels <- c(
    "Detectable loss",
    "No detectable change / uncertain",
    "Detectable gain"
  )
  
  areas <- vapply(
    categories,
    function(category_value) {
      
      category_area <- ifel(
        change_raster == category_value,
        area_raster,
        NA
      )
      
      global(
        category_area,
        sum,
        na.rm = TRUE
      )[1, 1]
    },
    numeric(1)
  )
  
  total_area <- sum(areas)
  
  data.frame(
    period = period_name,
    threshold_method = threshold_name,
    threshold_percentage_points = threshold_value,
    category_code = categories,
    category = category_labels,
    area_km2 = areas,
    percentage = 100 * areas / total_area
  )
}


## ============================================================
# 25. THRESHOLD-SENSITIVITY ANALYSIS
# ============================================================

if (length(annual_predictions) == 4L) {
  
  period_pairs <- list(
    "1994_2003" = c("1994", "2003"),
    "2003_2014" = c("2003", "2014"),
    "2014_2024" = c("2014", "2024"),
    "1994_2024" = c("1994", "2024")
  )
  
  threshold_values <- c(
    "Single_map_RMSE" =
      selected_model_rmse,
    
    "Combined_two_map_RMSE" =
      combined_change_rmse(
        selected_model_rmse,
        selected_model_rmse
      ),
    
    "Conservative_95_percent" =
      1.96 * combined_change_rmse(
        selected_model_rmse,
        selected_model_rmse
      )
  )
  
  sensitivity_outputs <- list()
  sensitivity_index <- 1L
  
  for (period_name in names(period_pairs)) {
    
    years <- period_pairs[[period_name]]
    
    earlier_map <- annual_predictions[[years[1]]]
    later_map   <- annual_predictions[[years[2]]]
    
    for (threshold_name in names(threshold_values)) {
      
      threshold_value <- threshold_values[[threshold_name]]
      
      output_prefix <- file.path(
        OUTPUT_DIR,
        paste0(
          "change_",
          period_name,
          "_",
          threshold_name
        )
      )
      
      change_result <- classify_canopy_change(
        earlier_raster = earlier_map,
        later_raster = later_map,
        threshold = threshold_value,
        output_prefix = output_prefix
      )
      
      sensitivity_outputs[[sensitivity_index]] <-
        summarize_change_area(
          change_raster = change_result$classes,
          period_name = period_name,
          threshold_name = threshold_name,
          threshold_value = threshold_value
        )
      
      sensitivity_index <- sensitivity_index + 1L
    }
  }
  
  sensitivity_table <- do.call(
    rbind,
    sensitivity_outputs
  )
  
  rownames(sensitivity_table) <- NULL
  
  print(sensitivity_table)
  
  write.csv(
    sensitivity_table,
    file.path(
      OUTPUT_DIR,
      "change_threshold_sensitivity_all_periods.csv"
    ),
    row.names = FALSE
  )
  
} else {
  
  warning(
    "Threshold-sensitivity analysis was not run because ",
    "annual_predictions contains ",
    length(annual_predictions),
    " maps instead of 4."
  )
}


# ============================================================
# 26. SAVE SESSION INFORMATION
# ============================================================

capture.output(
  sessionInfo(),
  file = file.path(
    OUTPUT_DIR,
    "session_information.txt"
  )
)

cat(
  "\nAnalysis completed.\n",
  "Review all warnings and output tables before updating the manuscript.\n",
  "Selected model: ", selected_model_name, "\n",
  "Spatial-CV RMSE: ", selected_model_rmse, "\n",
  sep = ""
)
