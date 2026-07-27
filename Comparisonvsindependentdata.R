
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# ------------------------------------------------------------
# 1. Read the CSV
# ------------------------------------------------------------

comparison_data <- read_csv(
  "comparisonpoints.csv",
  show_col_types = FALSE
)

# Assign clear names according to the column order
names(comparison_data)[1:3] <- c(
  "INFyS",
  "GFW",
  "FC2024"
)

# Check the data
str(comparison_data)
summary(comparison_data)

#Standarize units

convert_to_percent <- function(x) {
  if (max(x, na.rm = TRUE) <= 1) {
    x * 100
  } else {
    x
  }
}

comparison_data <- comparison_data |>
  mutate(
    across(
      c(INFyS, GFW, FC2024),
      convert_to_percent
    )
  )

#METRIC FUNCTION

calculate_metrics <- function(estimate, reference) {
  
  valid <- complete.cases(estimate, reference)
  
  estimate <- estimate[valid]
  reference <- reference[valid]
  
  if (length(estimate) < 3) {
    stop("Fewer than three valid paired observations.")
  }
  
  error <- estimate - reference
  r_value <- cor(estimate, reference, method = "pearson")
  
  ccc_value <- (
    2 * cov(estimate, reference)
  ) / (
    var(estimate) +
      var(reference) +
      (mean(estimate) - mean(reference))^2
  )
  
  tibble(
    n = length(estimate),
    RMSE = sqrt(mean(error^2)),
    MAE = mean(abs(error)),
    Bias = mean(error),
    r = r_value,
    R2 = r_value^2,
    CCC = ccc_value
  )
}


#CALCULATE ALL COMPARISON
comparison_metrics <- bind_rows(
  
  calculate_metrics(
    comparison_data$FC2024,
    comparison_data$INFyS
  ) |>
    mutate(
      Estimate = "FC2024",
      Reference = "INFyS"
    ),
  
   calculate_metrics(
    comparison_data$FC2024,
    comparison_data$GFW
  ) |>
    mutate(
      Estimate = "FC2024",
      Reference = "GFW"
    )
) |>
  select(
    Estimate,
    Reference,
    n,
    RMSE,
    MAE,
    Bias,
    r,
    R2,
    CCC
  ) |>
  mutate(
    across(
      c(RMSE, MAE, Bias, r, R2, CCC),
      ~ round(.x, 3)
    )
  )

# Create labels for the two panels
plot_labels <- comparison_metrics |>
  mutate(
    label = paste0(
      "r = ", sprintf("%.3f", r),
      "\nRMSE = ", sprintf("%.1f", RMSE),
      "\nBias = ", sprintf("%+.1f", Bias),
      "\nCCC = ", sprintf("%.3f", CCC)
    )
  ) |>
  select(Reference, label)


print(comparison_metrics)

write_csv(
  comparison_metrics,
  "independent_product_comparison_metrics.csv"
)

comparison_plot <- ggplot(
  plot_data,
  aes(
    x = Reference_cover,
    y = Estimated_cover
  )
) +
  geom_point(
    alpha = 0.5,
    size = 1.6,
    colour = "#2166AC"
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.7,
    colour = "black"
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    colour = "#B2182B",
    linewidth = 0.8
  ) +
  
  # Add statistics to each panel
  geom_text(
    data = plot_labels,
    aes(
      x = 4,
      y = 96,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.5
  ) +
  
  facet_wrap(
    ~Reference,
    nrow = 1
  ) +
  coord_equal(
    xlim = c(0, 100),
    ylim = c(0, 100),
    expand = FALSE
  ) +
  labs(
    x = "Independent-product canopy cover (%)",
    y = "MLP-estimated canopy cover for 2024 (%)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(
      fill = "grey90",
      colour = "grey40"
    ),
    strip.text = element_text(face = "bold")
  )

print(comparison_plot)

ggsave(
  filename = "independent_product_comparison_plot_annotated.png",
  plot = comparison_plot,
  width = 8,
  height = 4.5,
  dpi = 300
)

label = paste0(
  "r = ", sprintf("%.3f", r),
  "\nRMSE = ", sprintf("%.1f", RMSE), "%",
  "\nBias = ", sprintf("%+.1f", Bias), "%",
  "\nCCC = ", sprintf("%.3f", CCC)
)
