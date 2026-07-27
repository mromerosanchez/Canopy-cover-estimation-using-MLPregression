# ============================================================
# FINAL STATISTICS AND MANUSCRIPT TABLES
# Suelo de Conservacion canopy-cover analysis
#
# Outputs:
#   1. Masked annual canopy-cover rasters
#   2. Masked change-class rasters
#   3. Annual_canopy_statistics.csv
#   4. Table_change_by_period.csv
#   5. Table_threshold_sensitivity_all_periods.csv
#   6. Table_threshold_sensitivity_1994_2024.csv
#
# Required packages:
#   terra
#   dplyr
#   tidyr
#   readr
#   stringr
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(terra)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)


# ------------------------------------------------------------
# 2. USER SETTINGS
# ------------------------------------------------------------

# Folder containing the annual and change rasters
raster_dir <- "D:/Sabatico_geomatics/reviewer1_revised_analysis"

# Suelo de Conservacion shapefile
mask_file <- "D:/Sabatico_geomatics/SC_CDMX.shp"

# Folder for final outputs
output_dir <- file.path(
  raster_dir,
  "Suelo_de_Conservacion_outputs"
)

# Create output directories
dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

masked_annual_dir <- file.path(
  output_dir,
  "masked_annual_canopy"
)

masked_change_dir <- file.path(
  output_dir,
  "masked_change_classes"
)

table_dir <- file.path(
  output_dir,
  "manuscript_tables"
)

dir.create(
  masked_annual_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  masked_change_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 3. CLASS-CODE SETTINGS
# ------------------------------------------------------------

# Change these only if your change rasters use different codes.
#
# Expected coding:
#   -1 = Loss
#    0 = Stable / no detectable change
#    1 = Gain

LOSS_CODE   <- -1
STABLE_CODE <- 0
GAIN_CODE   <- 1


# ------------------------------------------------------------
# 4. EXPECTED ANNUAL RASTER FILES
# ------------------------------------------------------------

annual_files <- c(
  "1994" = "predicted_canopy_1994.tif",
  "2003" = "predicted_canopy_2003.tif",
  "2014" = "predicted_canopy_2014.tif",
  "2024" = "predicted_canopy_2024.tif"
)


# ------------------------------------------------------------
# 5. EXPECTED CHANGE PERIODS AND THRESHOLDS
# ------------------------------------------------------------

periods <- c(
  "1994_2003",
  "2003_2014",
  "2014_2024",
  "1994_2024"
)

threshold_types <- c(
  "Single_map_RMSE",
  "Combined_two_map_RMSE",
  "Conservative_95_percent"
)


# ------------------------------------------------------------
# 6. HELPER FUNCTIONS
# ------------------------------------------------------------

check_file_exists <- function(path) {
  
  if (!file.exists(path)) {
    stop(
      "File not found:\n",
      path,
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}


align_vector_to_raster <- function(vector_object, raster_object) {
  
  if (!same.crs(vector_object, raster_object)) {
    
    message(
      "Reprojecting mask to raster CRS."
    )
    
    vector_object <- project(
      vector_object,
      crs(raster_object)
    )
  }
  
  vector_object
}


mask_raster_to_polygon <- function(
    raster_file,
    mask_vector,
    output_file
) {
  
  check_file_exists(raster_file)
  
  r <- rast(raster_file)
  
  mask_aligned <- align_vector_to_raster(
    mask_vector,
    r
  )
  
  r_crop <- crop(
    r,
    mask_aligned,
    snap = "out"
  )
  
  r_mask <- mask(
    r_crop,
    mask_aligned
  )
  
  writeRaster(
    r_mask,
    output_file,
    overwrite = TRUE,
    wopt = list(
      datatype = datatype(r),
      gdal = c(
        "COMPRESS=LZW",
        "TILED=YES"
      )
    )
  )
  
  r_mask
}


safe_global_stat <- function(r, fun) {
  
  result <- global(
    r,
    fun = fun,
    na.rm = TRUE
  )
  
  as.numeric(result[1, 1])
}


detect_canopy_scale <- function(r) {
  
  r_max <- safe_global_stat(
    r,
    "max"
  )
  
  if (is.na(r_max)) {
    stop(
      "Raster contains no valid values.",
      call. = FALSE
    )
  }
  
  if (r_max <= 1.01) {
    return("proportion")
  }
  
  if (r_max <= 100.5) {
    return("percent")
  }
  
  warning(
    "Canopy raster contains values greater than 100. ",
    "Check the units before using the results."
  )
  
  "unknown"
}


convert_raster_to_percent <- function(r) {
  
  scale_type <- detect_canopy_scale(r)
  
  if (scale_type == "proportion") {
    
    message(
      "Converting canopy values from proportions to percentages."
    )
    
    r <- r * 100
  }
  
  r
}


calculate_polygon_area_ha <- function(mask_vector) {
  
  mask_area <- expanse(
    mask_vector,
    unit = "ha"
  )
  
  sum(
    mask_area,
    na.rm = TRUE
  )
}


calculate_valid_raster_area_ha <- function(r) {
  
  cell_areas <- cellSize(
    r,
    unit = "ha"
  )
  
  valid_area <- mask(
    cell_areas,
    r
  )
  
  safe_global_stat(
    valid_area,
    "sum"
  )
}


calculate_annual_statistics <- function(
    raster_object,
    year,
    study_area_ha
) {
  
  r_percent <- convert_raster_to_percent(
    raster_object
  )
  
  valid_values <- values(
    r_percent,
    mat = FALSE,
    na.rm = TRUE
  )
  
  if (length(valid_values) == 0) {
    stop(
      "No valid canopy values found for year ",
      year,
      ".",
      call. = FALSE
    )
  }
  
  valid_area_ha <- calculate_valid_raster_area_ha(
    r_percent
  )
  
  tibble(
    Year = as.integer(year),
    Valid_cells = length(valid_values),
    Study_area_ha = study_area_ha,
    Valid_raster_area_ha = valid_area_ha,
    Mean_canopy_percent = mean(valid_values),
    SD_canopy_percent = sd(valid_values),
    Median_canopy_percent = median(valid_values),
    Minimum_canopy_percent = min(valid_values),
    Maximum_canopy_percent = max(valid_values),
    Q25_canopy_percent = unname(
      quantile(valid_values, 0.25)
    ),
    Q75_canopy_percent = unname(
      quantile(valid_values, 0.75)
    )
  )
}


get_class_area_ha <- function(
    class_raster,
    class_code
) {
  
  class_mask <- class_raster == class_code
  
  area_raster <- cellSize(
    class_raster,
    unit = "ha"
  )
  
  selected_area <- mask(
    area_raster,
    class_mask,
    maskvalues = 0,
    updatevalue = NA
  )
  
  area_value <- global(
    selected_area,
    fun = "sum",
    na.rm = TRUE
  )
  
  if (
    nrow(area_value) == 0 ||
    is.na(area_value[1, 1])
  ) {
    return(0)
  }
  
  as.numeric(
    area_value[1, 1]
  )
}


calculate_change_statistics <- function(
    class_raster,
    period,
    threshold_name
) {
  
  raster_values <- unique(
    values(
      class_raster,
      mat = FALSE,
      na.rm = TRUE
    )
  )
  
  raster_values <- sort(
    raster_values
  )
  
  expected_values <- c(
    LOSS_CODE,
    STABLE_CODE,
    GAIN_CODE
  )
  
  unexpected_values <- setdiff(
    raster_values,
    expected_values
  )
  
  if (length(unexpected_values) > 0) {
    warning(
      "Unexpected class values in ",
      period,
      " / ",
      threshold_name,
      ": ",
      paste(
        unexpected_values,
        collapse = ", "
      )
    )
  }
  
  loss_ha <- get_class_area_ha(
    class_raster,
    LOSS_CODE
  )
  
  stable_ha <- get_class_area_ha(
    class_raster,
    STABLE_CODE
  )
  
  gain_ha <- get_class_area_ha(
    class_raster,
    GAIN_CODE
  )
  
  classified_area_ha <- sum(
    gain_ha,
    stable_ha,
    loss_ha
  )
  
  tibble(
    Period = str_replace_all(
      period,
      "_",
      "-"
    ),
    Threshold = threshold_name,
    Gain_ha = gain_ha,
    Stable_ha = stable_ha,
    Loss_ha = loss_ha,
    Classified_area_ha = classified_area_ha,
    Gain_percent = ifelse(
      classified_area_ha > 0,
      100 * gain_ha / classified_area_ha,
      NA_real_
    ),
    Stable_percent = ifelse(
      classified_area_ha > 0,
      100 * stable_ha / classified_area_ha,
      NA_real_
    ),
    Loss_percent = ifelse(
      classified_area_ha > 0,
      100 * loss_ha / classified_area_ha,
      NA_real_
    )
  )
}


format_area_percent <- function(
    area,
    percentage,
    area_digits = 1,
    percent_digits = 1
) {
  
  paste0(
    formatC(
      area,
      format = "f",
      digits = area_digits,
      big.mark = ","
    ),
    " (",
    formatC(
      percentage,
      format = "f",
      digits = percent_digits
    ),
    "%)"
  )
}


# ------------------------------------------------------------
# 7. READ SUelo DE CONSERVACION MASK
# ------------------------------------------------------------

check_file_exists(mask_file)

sc_mask <- vect(
  mask_file
)

if (is.lonlat(sc_mask)) {
  
  warning(
    "The Suelo de Conservacion shapefile uses geographic ",
    "coordinates. Raster cell areas will still be calculated ",
    "using terra::cellSize(), but a projected CRS is preferred."
  )
}

study_area_ha <- calculate_polygon_area_ha(
  sc_mask
)

message(
  "Suelo de Conservacion polygon area: ",
  round(study_area_ha, 2),
  " ha"
)


# ------------------------------------------------------------
# 8. MASK ANNUAL CANOPY-COVER RASTERS
# ------------------------------------------------------------

annual_results <- list()
annual_masked_rasters <- list()

for (year in names(annual_files)) {
  
  input_file <- file.path(
    raster_dir,
    annual_files[[year]]
  )
  
  output_file <- file.path(
    masked_annual_dir,
    paste0(
      "predicted_canopy_",
      year,
      "_Suelo_Conservacion.tif"
    )
  )
  
  message(
    "\nProcessing annual canopy raster: ",
    year
  )
  
  masked_raster <- mask_raster_to_polygon(
    raster_file = input_file,
    mask_vector = sc_mask,
    output_file = output_file
  )
  
  # Convert to percentage only for statistics.
  # Original values are preserved in the exported raster.
  annual_results[[year]] <- calculate_annual_statistics(
    raster_object = masked_raster,
    year = year,
    study_area_ha = study_area_ha
  )
  
  annual_masked_rasters[[year]] <- masked_raster
}


annual_statistics <- bind_rows(
  annual_results
) |>
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 3)
    )
  )


# Save annual statistics
write_csv(
  annual_statistics,
  file.path(
    table_dir,
    "Annual_canopy_statistics.csv"
  )
)

print(
  annual_statistics
)


# ------------------------------------------------------------
# 9. MASK ALL CHANGE-CLASS RASTERS
# ------------------------------------------------------------

change_results <- list()
result_counter <- 1L

for (period in periods) {
  
  for (threshold_type in threshold_types) {
    
    input_filename <- paste0(
      "change_",
      period,
      "_",
      threshold_type,
      "_classes.tif"
    )
    
    input_file <- file.path(
      raster_dir,
      input_filename
    )
    
    output_filename <- paste0(
      "change_",
      period,
      "_",
      threshold_type,
      "_classes_Suelo_Conservacion.tif"
    )
    
    output_file <- file.path(
      masked_change_dir,
      output_filename
    )
    
    message(
      "\nProcessing change raster: ",
      input_filename
    )
    
    if (!file.exists(input_file)) {
      
      warning(
        "Skipping missing file:\n",
        input_file
      )
      
      next
    }
    
    masked_change <- mask_raster_to_polygon(
      raster_file = input_file,
      mask_vector = sc_mask,
      output_file = output_file
    )
    
    change_results[[result_counter]] <-
      calculate_change_statistics(
        class_raster = masked_change,
        period = period,
        threshold_name = threshold_type
      )
    
    result_counter <- result_counter + 1L
  }
}


if (length(change_results) == 0) {
  stop(
    "No change-class rasters were processed.",
    call. = FALSE
  )
}


all_change_statistics <- bind_rows(
  change_results
) |>
  mutate(
    Threshold = recode(
      Threshold,
      "Single_map_RMSE" =
        "RMSE",
      "Combined_two_map_RMSE" =
        "RMSE_delta",
      "Conservative_95_percent" =
        "95% threshold"
    )
  ) |>
  arrange(
    factor(
      Period,
      levels = c(
        "1994-2003",
        "2003-2014",
        "2014-2024",
        "1994-2024"
      )
    ),
    factor(
      Threshold,
      levels = c(
        "RMSE",
        "RMSE_delta",
        "95% threshold"
      )
    )
  )


# Save full numerical results
write_csv(
  all_change_statistics,
  file.path(
    table_dir,
    "Change_statistics_all_periods_all_thresholds.csv"
  )
)


# ------------------------------------------------------------
# 10. MAIN MANUSCRIPT CHANGE TABLE
# Uses propagated two-map RMSE
# ------------------------------------------------------------

main_change_table_numeric <- all_change_statistics |>
  filter(
    Threshold == "RMSE_delta"
  ) |>
  select(
    Period,
    Gain_ha,
    Stable_ha,
    Loss_ha,
    Gain_percent,
    Stable_percent,
    Loss_percent
  ) |>
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 2)
    )
  )


write_csv(
  main_change_table_numeric,
  file.path(
    table_dir,
    "Table_change_by_period_numeric.csv"
  )
)


# Manuscript-ready version:
# area in hectares, percentage in parentheses

main_change_table_manuscript <-
  main_change_table_numeric |>
  transmute(
    Period,
    Gain = format_area_percent(
      Gain_ha,
      Gain_percent
    ),
    `No detectable change` =
      format_area_percent(
        Stable_ha,
        Stable_percent
      ),
    Loss = format_area_percent(
      Loss_ha,
      Loss_percent
    )
  )


write_csv(
  main_change_table_manuscript,
  file.path(
    table_dir,
    "Table_change_by_period_manuscript.csv"
  )
)

print(
  main_change_table_manuscript
)


# ------------------------------------------------------------
# 11. THRESHOLD-SENSITIVITY TABLE: ALL PERIODS
# ------------------------------------------------------------

threshold_table_all_periods_numeric <-
  all_change_statistics |>
  select(
    Period,
    Threshold,
    Gain_ha,
    Stable_ha,
    Loss_ha,
    Gain_percent,
    Stable_percent,
    Loss_percent
  ) |>
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 2)
    )
  )


write_csv(
  threshold_table_all_periods_numeric,
  file.path(
    table_dir,
    "Table_threshold_sensitivity_all_periods_numeric.csv"
  )
)


threshold_table_all_periods_manuscript <-
  threshold_table_all_periods_numeric |>
  transmute(
    Period,
    Threshold,
    Gain = format_area_percent(
      Gain_ha,
      Gain_percent
    ),
    Stable = format_area_percent(
      Stable_ha,
      Stable_percent
    ),
    Loss = format_area_percent(
      Loss_ha,
      Loss_percent
    )
  )


write_csv(
  threshold_table_all_periods_manuscript,
  file.path(
    table_dir,
    "Table_threshold_sensitivity_all_periods_manuscript.csv"
  )
)


# ------------------------------------------------------------
# 12. THRESHOLD-SENSITIVITY TABLE: 1994–2024 ONLY
# Matches the three-row manuscript table supplied
# ------------------------------------------------------------

threshold_table_1994_2024_numeric <-
  all_change_statistics |>
  filter(
    Period == "1994-2024"
  ) |>
  select(
    Threshold,
    Gain_ha,
    Stable_ha,
    Loss_ha,
    Gain_percent,
    Stable_percent,
    Loss_percent
  ) |>
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 2)
    )
  )


write_csv(
  threshold_table_1994_2024_numeric,
  file.path(
    table_dir,
    "Table_threshold_sensitivity_1994_2024_numeric.csv"
  )
)


threshold_table_1994_2024_manuscript <-
  threshold_table_1994_2024_numeric |>
  transmute(
    Threshold,
    Gain = format_area_percent(
      Gain_ha,
      Gain_percent
    ),
    Stable = format_area_percent(
      Stable_ha,
      Stable_percent
    ),
    Loss = format_area_percent(
      Loss_ha,
      Loss_percent
    )
  )


write_csv(
  threshold_table_1994_2024_manuscript,
  file.path(
    table_dir,
    "Table_threshold_sensitivity_1994_2024_manuscript.csv"
  )
)

print(
  threshold_table_1994_2024_manuscript
)


# ------------------------------------------------------------
# 13. OPTIONAL WIDE TABLE FOR ALL PERIODS AND THRESHOLDS
# ------------------------------------------------------------

threshold_table_wide <-
  threshold_table_all_periods_manuscript |>
  pivot_wider(
    names_from = Threshold,
    values_from = c(
      Gain,
      Stable,
      Loss
    )
  )


write_csv(
  threshold_table_wide,
  file.path(
    table_dir,
    "Table_threshold_sensitivity_wide.csv"
  )
)


# ------------------------------------------------------------
# 14. BASIC CONSISTENCY CHECKS
# ------------------------------------------------------------

message(
  "\nChecking whether class percentages sum to 100..."
)

percentage_checks <- all_change_statistics |>
  mutate(
    Percentage_sum =
      Gain_percent +
      Stable_percent +
      Loss_percent
  ) |>
  select(
    Period,
    Threshold,
    Percentage_sum
  )

print(
  percentage_checks
)


if (
  any(
    abs(
      percentage_checks$Percentage_sum - 100
    ) > 0.1,
    na.rm = TRUE
  )
) {
  
  warning(
    "Some gain/stable/loss percentages do not sum to ",
    "approximately 100%. Inspect missing or unexpected ",
    "class values."
  )
}


message(
  "\nAnalysis complete."
)

message(
  "\nOutputs saved to:\n",
  normalizePath(
    output_dir,
    winslash = "/",
    mustWork = FALSE
  )
)
