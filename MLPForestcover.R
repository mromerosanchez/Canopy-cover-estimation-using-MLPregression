# forest_cover_mlp.R
# ------------------------------------------------------------
# Canopy cover estimation using Multilayer Perceptron regression
# ------------------------------------------------------------

# 1. LOAD LIBRARIES
library(torch)
library(terra)
library(readr)
library(dplyr)
library(purrr)

# Optional, only used for session information
# install.packages("ggplot2") if needed
library(ggplot2)

# ------------------------------------------------------------
# 2. GLOBAL SETTINGS
# ------------------------------------------------------------

set.seed(42)
torch_manual_seed(42)

epochs <- 100
learning_rate <- 0.001
dropout_rate <- 0.30

years <- c(1994, 2003, 2014, 2024)

# ------------------------------------------------------------
# 3. HELPER FUNCTIONS
# ------------------------------------------------------------

compute_metrics <- function(observed, predicted) {
  data.frame(
    RMSE = sqrt(mean((predicted - observed)^2, na.rm = TRUE)),
    MAE  = mean(abs(predicted - observed), na.rm = TRUE),
    Bias = mean(predicted - observed, na.rm = TRUE),
    r    = cor(predicted, observed, use = "complete.obs")
  )
}

compute_class_metrics <- function(validation_df) {
  validation_df %>%
    group_by(canopy_class) %>%
    summarise(
      n = n(),
      RMSE = sqrt(mean((predicted - observed)^2, na.rm = TRUE)),
      MAE = mean(abs(predicted - observed), na.rm = TRUE),
      Bias = mean(predicted - observed, na.rm = TRUE),
      r = cor(predicted, observed, use = "complete.obs"),
      .groups = "drop"
    )
}

# ------------------------------------------------------------
# 4. LOAD AND PREPARE TRAINING DATA
# ------------------------------------------------------------

df <- read_csv("training_data_1500b.csv")

y_raw <- df[[1]]       # forest/canopy cover percentage, 0–100
X_df  <- df[, -1]      # 11 spectral predictors

# Remove incomplete records
complete_rows <- complete.cases(df)
df <- df[complete_rows, ]
y_raw <- df[[1]]
X_df  <- df[, -1]

# Standardize predictors
scaler_means <- map_dbl(X_df, mean, na.rm = TRUE)
scaler_sds   <- map_dbl(X_df, sd, na.rm = TRUE)

# Avoid division by zero
scaler_sds[scaler_sds == 0] <- 1

X_scaled <- sweep(as.matrix(X_df), 2, scaler_means, "-")
X_scaled <- sweep(X_scaled, 2, scaler_sds, "/")

# Scale target to 0–1
y_scaled <- y_raw / 100

# Convert to torch tensors
X_tensor <- torch_tensor(X_scaled, dtype = torch_float())
y_tensor <- torch_tensor(matrix(y_scaled, ncol = 1), dtype = torch_float())

# ------------------------------------------------------------
# 5. STRATIFIED TRAIN/TEST SPLIT
# ------------------------------------------------------------

df_split <- df %>%
  mutate(
    canopy_class = cut(
      y_raw,
      breaks = c(0, 20, 40, 60, 80, 100),
      labels = c("0-20", "20-40", "40-60", "60-80", "80-100"),
      include.lowest = TRUE
    ),
    row_id = row_number()
  )

train_idx <- df_split %>%
  group_by(canopy_class) %>%
  sample_frac(0.8) %>%
  pull(row_id)

test_idx <- setdiff(seq_len(nrow(df)), train_idx)

X_train <- X_tensor[train_idx, ]
y_train <- y_tensor[train_idx, ]

X_test <- X_tensor[test_idx, ]
y_test <- y_tensor[test_idx, ]

sample_distribution <- df_split %>%
  count(canopy_class)

write_csv(sample_distribution, "sample_distribution_by_canopy_class.csv")

# ------------------------------------------------------------
# 6. DEFINE MLP MODEL
# ------------------------------------------------------------

ForestNet <- nn_module(
  "ForestNet",
  
  initialize = function(input_dim = 11) {
    self$fc1 <- nn_linear(input_dim, 64)
    self$drop <- nn_dropout(p = dropout_rate)
    self$fc2 <- nn_linear(64, 32)
    self$output <- nn_linear(32, 1)
  },
  
  forward = function(x) {
    x %>%
      self$fc1() %>%
      nnf_relu() %>%
      self$drop() %>%
      self$fc2() %>%
      nnf_relu() %>%
      self$output()
  }
)

model <- ForestNet(input_dim = ncol(X_scaled))

# ------------------------------------------------------------
# 7. TRAIN MODEL
# ------------------------------------------------------------

optimizer <- optim_adam(model$parameters, lr = learning_rate)
loss_fn <- nn_mse_loss()

train_losses <- numeric(epochs)

for (epoch in 1:epochs) {
  
  model$train()
  optimizer$zero_grad()
  
  preds <- model(X_train)
  loss <- loss_fn(preds, y_train)
  
  loss$backward()
  optimizer$step()
  
  train_losses[epoch] <- loss$item()
  
  if (epoch %% 10 == 0) {
    cat(sprintf("Epoch %3d | Training loss = %.5f\n", epoch, loss$item()))
  }
}

# ------------------------------------------------------------
# 8. TRAINING DIAGNOSTIC FIGURE
# ------------------------------------------------------------

png(
  "training_loss_curve.png",
  width = 1200,
  height = 800,
  res = 150
)

plot(
  1:epochs,
  train_losses,
  type = "l",
  col = "forestgreen",
  lwd = 3,
  xlab = "Epoch",
  ylab = "Training loss",
  main = "MLP Training Loss"
)

dev.off()

write_csv(
  data.frame(epoch = 1:epochs, training_loss = train_losses),
  "training_loss_values.csv"
)

# ------------------------------------------------------------
# 9. VALIDATE MODEL
# ------------------------------------------------------------

model$eval()

with_no_grad({
  pred_test_scaled <- model(X_test)
})

pred_test <- as.numeric(pred_test_scaled$cpu()) * 100
obs_test  <- as.numeric(y_test$cpu()) * 100

# Restrict predictions to valid canopy-cover range
pred_test <- pmin(pmax(pred_test, 0), 100)

validation_df <- data.frame(
  observed = obs_test,
  predicted = pred_test
) %>%
  mutate(
    canopy_class = cut(
      observed,
      breaks = c(0, 20, 40, 60, 80, 100),
      labels = c("0-20", "20-40", "40-60", "60-80", "80-100"),
      include.lowest = TRUE
    )
  )

metrics <- compute_metrics(validation_df$observed, validation_df$predicted)
class_metrics <- compute_class_metrics(validation_df)

print(metrics)
print(class_metrics)

write_csv(metrics, "model_validation_metrics.csv")
write_csv(class_metrics, "per_class_validation_metrics.csv")
write_csv(validation_df, "validation_predictions.csv")

rmse <- metrics$RMSE[1]
mae  <- metrics$MAE[1]
bias <- metrics$Bias[1]
rval <- metrics$r[1]

# ------------------------------------------------------------
# 10. VALIDATION FIGURE
# ------------------------------------------------------------

png(
  "validation_observed_vs_predicted.png",
  width = 1200,
  height = 1200,
  res = 150
)

plot(
  validation_df$observed,
  validation_df$predicted,
  pch = 16,
  col = rgb(0, 0.45, 0, 0.45),
  xlim = c(0, 100),
  ylim = c(0, 100),
  xlab = "Observed canopy cover (%)",
  ylab = "Predicted canopy cover (%)",
  main = "Observed vs Predicted Canopy Cover"
)

abline(0, 1, col = "red", lwd = 2, lty = 2)

legend(
  "topleft",
  legend = c(
    paste0("RMSE = ", round(rmse, 2)),
    paste0("MAE = ", round(mae, 2)),
    paste0("Bias = ", round(bias, 2)),
    paste0("r = ", round(rval, 2))
  ),
  bty = "n"
)

dev.off()

# ------------------------------------------------------------
# 11. SAVE MODEL AND METADATA
# ------------------------------------------------------------

torch_save(model, "forest_model.pt")

saveRDS(
  list(
    means = scaler_means,
    sds = scaler_sds,
    rmse = rmse,
    mae = mae,
    bias = bias,
    r = rval,
    epochs = epochs,
    learning_rate = learning_rate,
    dropout_rate = dropout_rate,
    predictor_names = names(X_df)
  ),
  "forest_model_metadata.rds"
)

# ------------------------------------------------------------
# 12. PREDICT CANOPY COVER FROM RASTER STACKS
# ------------------------------------------------------------

predict_forest_cover <- function(
    rast_path,
    out_path,
    model,
    scaler_means,
    scaler_sds
) {
  
  st <- rast(rast_path)
  
  if (nlyr(st) != length(scaler_means)) {
    stop(
      "Expected ",
      length(scaler_means),
      " bands, got ",
      nlyr(st)
    )
  }
  
  cat("Predicting canopy cover for:", rast_path, "\n")
  
  mat <- as.matrix(st)
  valid <- complete.cases(mat)
  
  cat("Valid pixels:", sum(valid), "\n")
  
  scaled <- sweep(mat[valid, ], 2, scaler_means, "-")
  scaled <- sweep(scaled, 2, scaler_sds, "/")
  
  input_t <- torch_tensor(scaled, dtype = torch_float())
  
  model$eval()
  
  with_no_grad({
    preds_t <- model(input_t)$squeeze()
  })
  
  preds <- as.numeric(preds_t$cpu()) * 100
  preds <- pmin(pmax(preds, 0), 100)
  
  out_r <- st[[1]]
  vals <- rep(NA_real_, ncell(out_r))
  vals[valid] <- preds
  
  out_r <- setValues(out_r, vals)
  names(out_r) <- "canopy_cover"
  
  writeRaster(out_r, out_path, overwrite = TRUE)
  
  cat("Wrote:", out_path, "\n")
  
  return(out_r)
}

for (yr in years) {
  
  in_file  <- sprintf("Landsat%d.tif", yr)
  out_file <- sprintf("forest_cover_%d.tif", yr)
  
  predict_forest_cover(
    rast_path = in_file,
    out_path = out_file,
    model = model,
    scaler_means = scaler_means,
    scaler_sds = scaler_sds
  )
}

# ------------------------------------------------------------
# 13. UNCERTAINTY-AWARE CHANGE DETECTION
# ------------------------------------------------------------

change_class <- function(
    earlier_path,
    later_path,
    out_change_path,
    threshold
) {
  
  earlier <- rast(earlier_path)
  later   <- rast(later_path)
  
  delta <- later - earlier
  names(delta) <- "canopy_change"
  
  change <- classify(
    delta,
    rcl = matrix(
      c(
        -Inf, -threshold, -1,
        -threshold, threshold, 0,
        threshold, Inf, 1
      ),
      ncol = 3,
      byrow = TRUE
    )
  )
  
  names(change) <- "change_class"
  
  delta_path <- sub(".tif$", "_delta.tif", out_change_path)
  
  writeRaster(delta, delta_path, overwrite = TRUE)
  writeRaster(change, out_change_path, overwrite = TRUE)
  
  cat("Wrote delta raster:", delta_path, "\n")
  cat("Wrote change-class raster:", out_change_path, "\n")
  
  return(change)
}

threshold <- rmse

cat("Using RMSE-based change threshold:", threshold, "% canopy cover\n")

chg_1994_2003 <- change_class(
  earlier_path = "forest_cover_1994.tif",
  later_path = "forest_cover_2003.tif",
  out_change_path = "change_1994_2003.tif",
  threshold = threshold
)

chg_2003_2014 <- change_class(
  earlier_path = "forest_cover_2003.tif",
  later_path = "forest_cover_2014.tif",
  out_change_path = "change_2003_2014.tif",
  threshold = threshold
)

chg_2014_2024 <- change_class(
  earlier_path = "forest_cover_2014.tif",
  later_path = "forest_cover_2024.tif",
  out_change_path = "change_2014_2024.tif",
  threshold = threshold
)

# Optional full-period change
chg_1994_2024 <- change_class(
  earlier_path = "forest_cover_1994.tif",
  later_path = "forest_cover_2024.tif",
  out_change_path = "change_1994_2024.tif",
  threshold = threshold
)

# ------------------------------------------------------------
# 14. CHANGE AREA STATISTICS
# ------------------------------------------------------------

area_stats <- function(change_raster, out_csv) {
  
  freq_df <- as.data.frame(terra::freq(change_raster, useNA = "no"))
  
  print(freq_df)
  print(names(freq_df))
  
  if (is.null(freq_df) || nrow(freq_df) == 0) {
    warning("No valid pixels found in raster.")
    return(NULL)
  }
  
  # terra::freq often returns columns: layer, value, count
  if ("value" %in% names(freq_df) && "count" %in% names(freq_df)) {
    freq_df <- freq_df %>%
      dplyr::select(value, count)
  } else if (ncol(freq_df) == 3) {
    freq_df <- freq_df[, c(2, 3)]
  } else if (ncol(freq_df) == 2) {
    freq_df <- freq_df[, c(1, 2)]
  } else {
    stop("Unexpected terra::freq() output structure.")
  }
  
  names(freq_df) <- c("class_value", "pixels")
  
  freq_df$class_value <- as.numeric(freq_df$class_value)
  freq_df$pixels <- as.numeric(freq_df$pixels)
  
  pixel_area_ha <- prod(terra::res(change_raster)) / 10000
  
  freq_df <- freq_df %>%
    dplyr::mutate(
      area_ha = pixels * pixel_area_ha,
      percent = 100 * area_ha / sum(area_ha, na.rm = TRUE),
      class_label = dplyr::case_when(
        class_value == -1 ~ "Loss",
        class_value == 0 ~ "Persistence/uncertain",
        class_value == 1 ~ "Gain",
        TRUE ~ "Unknown"
      )
    ) %>%
    dplyr::select(class_label, class_value, pixels, area_ha, percent)
  
  readr::write_csv(freq_df, out_csv)
  
  return(freq_df)
}
