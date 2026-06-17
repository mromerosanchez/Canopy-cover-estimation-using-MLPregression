# Load required library
library(terra)


# ------------------------------
# Define the index calculation function
# ------------------------------
calculate_indices <- function(img_stack, satellite) {
  if (!satellite %in% c("L5", "L7", "L8")) {
    stop("`satellite` must be one of 'L5', 'L7', or 'L8'.")
  }
  
  # Band indices per satellite type
  bands <- switch(satellite,
                  "L5" = list(blue = 1, green = 2, red = 3, nir = 4, swir1 = 5, swir2 = 6),
                  "L7" = list(blue = 1, green = 2, red = 3, nir = 4, swir1 = 5, swir2 = 6),
                  "L8" = list(blue = 1, green = 2, red = 3, nir = 4, swir1 = 5, swir2 = 6)
  )
  
  b <- img_stack[[bands$blue]]
  g <- img_stack[[bands$green]]
  r <- img_stack[[bands$red]] 
  n <- img_stack[[bands$nir]]
  s1 <- img_stack[[bands$swir1]]
  s2 <- img_stack[[bands$swir2]]
  
  ndvi <- (n - r) / (n + r)
  kndvi <- tanh((ndvi^2) / (0.1^2))
  savi <- ((n - r) / (n + r + 0.5)) * 1.5
  ndwi <- (g - n) / (g + n)
  evi <- 2.5 * (n - r) / (n + 6 * r - 7.5 * b + 1)
  
  coeffs <- switch(satellite,
                   "L5" = list(
                     brightness = c(0.3037, 0.2793, 0.4743, 0.5585, 0.5082, 0.1863),
                     greenness  = c(-0.2848, -0.2435, -0.5436, 0.7243, 0.0840, -0.1800)
                   ),
                   "L7" = list(
                     brightness = c(0.3037, 0.2793, 0.4743, 0.5585, 0.5082, 0.1863),
                     greenness  = c(-0.2848, -0.2435, -0.5436, 0.7243, 0.0840, -0.1800)
                   ),
                   "L8" = list(
                     brightness = c(0.3029, 0.2786, 0.4733, 0.5599, 0.5080, 0.1872),
                     greenness  = c(-0.2941, -0.2430, -0.5424, 0.7276, 0.0713, -0.1608)
                   )
  )
  
  brightness <- coeffs$brightness[[1]] * r + coeffs$brightness[[2]] * g +
    coeffs$brightness[[3]] * b + coeffs$brightness[[4]] * n +
    coeffs$brightness[[5]] * s1 + coeffs$brightness[[6]] * s2
  
  greenness <- coeffs$greenness[[1]] * r + coeffs$greenness[[2]] * g +
    coeffs$greenness[[3]] * b + coeffs$greenness[[4]] * n +
    coeffs$greenness[[5]] * s1 + coeffs$greenness[[6]] * s2
  
  indices <- c(ndvi, kndvi, savi, ndwi, evi, greenness, brightness)
  names(indices) <- c("NDVI", "kNDVI", "SAVI", "NDWI", "EVI", "Greenness", "Brightness")
  return(indices)
}

# ------------------------------
# Define file list with metadata
# ------------------------------
landsat_files <- list(
  list(file = "Landsat5_Reflectance_Only_1994.tif", satellite = "L5", year = 1994),
  list(file = "Landsat7_Reflectance_Only_2003.tif", satellite = "L7", year = 2003),
  list(file = "Landsat8_Reflectance_Only_2014.tif", satellite = "L8", year = 2014),
  list(file = "Landsat8_Reflectance_Only_2024.tif", satellite = "L8", year = 2024)
)

# ------------------------------
# Process and save indices
# ------------------------------
for (info in landsat_files) {
  cat("Processing", info$file, "...\n")
  
  stack <- rast(info$file)
  result <- calculate_indices(stack, satellite = info$satellite)
  
  output_file <- paste0("vegetation_indices_", info$year, ".tif")
  writeRaster(result, output_file, overwrite = TRUE)
  
  cat("Saved:", output_file, "\n")
}
