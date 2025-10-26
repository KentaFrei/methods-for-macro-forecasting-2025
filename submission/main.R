# TODO: Is the cutoff at 0.99 correlation common practice? --> Now is set to 0.95!!!
# 0.99 --> 0.95: Log Likelihood: from -167.4 to -149.6, Max root: from 1.007 to 1.004
# TODO compare with their forecasts, reconstruction error there
# TODO review if this model setup is right --> tested with "summary(pca_X)" and "roots(var_model)"
# TODO fix the missing gap in the data
install.packages("vars") 
install.packages("ggplot2")
install.packages("tidyr")
library(ggplot2)
library(dplyr)
library(vars) 
library(tidyr)   

df <- utils::read.csv("data/data_quarterly.csv")
head(df)
# To build our model, we focus on data in the past 
df_train <- df %>%
  mutate(
    date = as.Date(paste0(date, "-01"))
  ) %>%
  filter(date < as.Date("2025-10-01"))
# Safety check to make sure all vals are numeric
df_train <- df_train %>%
  mutate(across(-date, as.numeric))
# Check to confirm no missing values
sapply(df_train, function(x) sum(is.na(x))) %>% sort(decreasing = TRUE)

metadata <- utils::read.csv(“data/metadata_quarterly_en.csv”)

# We need to transform the raw values to make the series stationary
unit_to_transform <- c(
  "real in millions of francs" = "logdiff",
  "Index" = "logdiff",
  "in millions of francs at current prices" = "logdiff",
  "Gross domestic product (adjusted for sporting events) per hour worked" = "logdiff",
  "Gross domestic product (adjusted for sporting events) per full-time equivalent" = "logdiff",
  "in thousands of persons" = "diff",
  "Price deflator" = "logdiff",
  "in US dollars" = "logdiff",
  "Real value of gross domestic product per full-time equivalent" = "logdiff",
  "in percent" = "diff",
  "real in millions of US dollars" = "logdiff",
  "in thousands of Swiss francs at current prices" = "logdiff",
  "real in thousands of Swiss francs" = "logdiff",
  "rate" = "diff"
)
metadata$transform <- unit_to_transform[metadata$unit]
# Now each column has its transformation stored with it in metadata

#Helper function
transform_variable <- function(x, transform_type) {
  if (transform_type == "logdiff") {
    return(c(NA, diff(log(x)) * 100))
  } else if (transform_type == "diff") {
    return(c(NA, diff(x)))
  } else {
    return(x)
  }
}

# Apply transformations to the df_clean dataset
df_train_clean <- df_train

for (i in seq_len(nrow(metadata))) {
  var <- metadata$variable[i]
  t_type <- metadata$transform[i]
  
  # Skip if variable not in df_train
  if (!var %in% names(df_train)) next
  
  df_train_clean[[var]] <- transform_variable(df_train[[var]], t_type)
}
# First row using diff makes na
df_train_clean <- df_train_clean[-1, ]
head(df_train_clean)

# Now we need to standardize everything
df_train_std <- df_train_clean

df_train_std[-1] <- scale(df_train_clean[-1]) # scale subtracts col mean and divides by SD
head(df_train_std)

# Saving the means and sd just in case we need later
scaler_params <- data.frame(
  variable = names(df_train_clean)[-1],
  mean = apply(df_train_clean[-1], 2, mean, na.rm = TRUE),
  sd = apply(df_train_clean[-1], 2, sd, na.rm = TRUE)
)


# Remove the variables we are trying to predice
vars_Y <- c("gdp", "cpi", "unempoff")
vars_X <- setdiff(names(df_train_std)[-1], vars_Y)

# Remove variables which are highly correlated (nominal gdp and real gdp, for instance) to avoid overfitting
cor_mat <- cor(df_train_std[, vars_X], use = "pairwise.complete.obs")
threshold <- 0.95
hi_pairs <- which(abs(cor_mat) > threshold & upper.tri(cor_mat), arr.ind = TRUE)
redundant <- unique(rownames(cor_mat)[hi_pairs[, 1]])
vars_X_filtered <- setdiff(vars_X, redundant)
cat("Removed", length(redundant), "variables with |corr| >", threshold, "\n")

# run pca with 5 component vectors
pca_X <- prcomp(df_train_std[, vars_X_filtered], center = FALSE, scale. = FALSE)
r <- 5
F_hat <- pca_X$x[, 1:r]
Lambda_hat <- pca_X$rotation[, 1:r]
F_df <- as.data.frame(F_hat)

# fit a simple var model with 4 quarter (1yr lag)
var_model <- VAR(F_df, ic = "AIC", lag.max = 4)
summary(var_model)


# Regressions on factors (IN-SAMPLE)
Y_df <- df_train_std[, vars_Y]

beta_hat <- matrix(NA, nrow = ncol(F_df), ncol = ncol(Y_df))
colnames(beta_hat) <- vars_Y
for (j in seq_along(vars_Y)) {
  fit <- lm(Y_df[[j]] ~ as.matrix(F_df))
  beta_hat[, j] <- coef(fit)[-1]
}

# In-sample predictions (standardised scale)
Y_pred_std <- as.matrix(F_df) %*% beta_hat
colnames(Y_pred_std) <- vars_Y

# R^2 in-sample
r2_values <- sapply(seq_len(ncol(Y_df)), function(i) {
  1 - mean((Y_df[, i] - Y_pred_std[, i])^2, na.rm = TRUE) /
    mean((Y_df[, i] - mean(Y_df[, i], na.rm = TRUE))^2, na.rm = TRUE)
})
r2_table <- data.frame(variable = vars_Y, R2 = round(r2_values, 3))
print(r2_table)

# # Destandardise in-sample predictions
scaler_subset <- scaler_params[match(vars_Y, scaler_params$variable), ]
Y_pred_real <- sweep(Y_pred_std, 2, scaler_subset$sd, "*")
Y_pred_real <- sweep(Y_pred_real, 2, scaler_subset$mean, "+")


# Forecast multi-CI (50%, 80%, 95%) and chart
# Future outlook and dates
h <- 8
last_date <- max(df_train_clean$date)
future_dates <- seq(last_date, by = "quarter", length.out = h + 1)[-1]

# Compute forecast bands
forecast_bands <- function(var_model, beta_hat, scaler_params, vars_Y, h, ci, future_dates) {
  F_fc <- predict(var_model, n.ahead = h, ci = ci)
  
  # h x r matrices of predicted factors (fcst/low/up)
  F_forecast_mat <- as.matrix(sapply(F_fc$fcst, function(x) x[, "fcst"]))
  F_lower_mat    <- as.matrix(sapply(F_fc$fcst, function(x) x[, "lower"]))
  F_upper_mat    <- as.matrix(sapply(F_fc$fcst, function(x) x[, "upper"]))
  
  # Projection on targets (standardised scale)
  Y_forecast_std <- F_forecast_mat %*% beta_hat
  Y_lower_std    <- F_lower_mat    %*% beta_hat
  Y_upper_std    <- F_upper_mat    %*% beta_hat
  
  colnames(Y_forecast_std) <- vars_Y
  colnames(Y_lower_std)    <- vars_Y
  colnames(Y_upper_std)    <- vars_Y
  
  # Destandardize
  scaler_subset <- scaler_params[match(vars_Y, scaler_params$variable), ]
  Y_forecast_real <- sweep(Y_forecast_std, 2, scaler_subset$sd, "*")
  Y_forecast_real <- sweep(Y_forecast_real, 2, scaler_subset$mean, "+")
  Y_lower_real    <- sweep(Y_lower_std,    2, scaler_subset$sd, "*")
  Y_lower_real    <- sweep(Y_lower_real,   2, scaler_subset$mean, "+")
  Y_upper_real    <- sweep(Y_upper_std,    2, scaler_subset$sd, "*")
  Y_upper_real    <- sweep(Y_upper_real,   2, scaler_subset$mean, "+")
  
  # Data frame "long" for ggplot
  df_mid <- data.frame(date = future_dates, Y_forecast_real)
  df_low <- data.frame(date = future_dates, Y_lower_real)
  df_up  <- data.frame(date = future_dates, Y_upper_real)
  
  long_mid <- tidyr::pivot_longer(df_mid, cols = dplyr::all_of(vars_Y),
                                  names_to = "variable", values_to = "value")
  long_low <- tidyr::pivot_longer(df_low, cols = dplyr::all_of(vars_Y),
                                  names_to = "variable", values_to = "lower")
  long_up  <- tidyr::pivot_longer(df_up,  cols = dplyr::all_of(vars_Y),
                                  names_to = "variable", values_to = "upper")
  
  bands <- long_mid %>%
    dplyr::left_join(long_low, by = c("date", "variable")) %>%
    dplyr::left_join(long_up,  by = c("date", "variable")) %>%
    dplyr::mutate(ci = paste0(round(ci * 100), "%"))
  
  list(
    mid   = long_mid %>% dplyr::mutate(ci = paste0(round(ci * 100), "%")),
    bands = bands
  )
}

# Calculate bands for multiple levels
cis <- c(0.50, 0.80, 0.95)
res_list <- lapply(cis, function(x) forecast_bands(var_model, beta_hat, scaler_params, vars_Y, h, x, future_dates))

# Central forecast (mid) and bands in ‘long’ format
forecast_long_all <- dplyr::bind_rows(lapply(res_list, `[[`, "mid"))
bands_long_all    <- dplyr::bind_rows(lapply(res_list, `[[`, "bands"))

# Observed/predicted/forecast data for lines 
actual_df <- df_train_clean[, c("date", vars_Y)]
pred_df   <- data.frame(date = df_train_clean$date, Y_pred_real)

long_actual <- tidyr::pivot_longer(actual_df, cols = vars_Y,
                                   names_to = "variable", values_to = "value") %>%
  dplyr::mutate(type = "Observed")

long_pred <- tidyr::pivot_longer(pred_df, cols = vars_Y,
                                 names_to = "variable", values_to = "value") %>%
  dplyr::mutate(type = "Predicted")

long_fore <- forecast_long_all %>%
  dplyr::mutate(type = "Forecast")

combined_long <- dplyr::bind_rows(long_actual, long_pred, long_fore)

# Single plot with multi-CI bands + lines 
bands_long_all$ci <- factor(bands_long_all$ci, levels = c("95%", "80%", "50%"))

ggplot() +
  # Bands (first, so they remain below the lines)
  geom_ribbon(data = bands_long_all,
              aes(x = date, ymin = lower, ymax = upper, fill = ci),
              alpha = 0.2) +
  # Original lines
  geom_line(data = combined_long,
            aes(x = date, y = value, color = type, linetype = type),
            linewidth = 0.5) +
  facet_wrap(~ variable, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c("Observed" = "black",
                                "Predicted" = "darkgreen",
                                "Forecast"  = "steelblue")) +
  scale_linetype_manual(values = c("Observed" = "solid",
                                   "Predicted" = "solid",
                                   "Forecast"  = "solid")) +
  # Legend for CI
  scale_fill_manual(name = "Conf. Interval",
                    values = c("95%" = "grey60", "80%" = "grey50", "50%" = "grey40")) +
  labs(
    title = "GDP, CPI, and Unemployment",
    subtitle = paste0("(R²: ",
                      paste(r2_table$variable, r2_table$R2, collapse = ", "), ")"),
    x = "Date", y = "Value", color = "", linetype = ""
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14)
  )

## Below is just demo code
# Create sample data
data <- data.frame(
  x = rnorm(100),
  y = rnorm(100),
  group = sample(c("A", "B", "C"), 100, replace = TRUE)
)

# Basic analysis
cat("Data summary:\n")
print(summary(data))

cat("\nGroup counts:\n")
print(table(data$group))

# # Create a plot
p <- ggplot(data, aes(x = x, y = y, color = group)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Sample Scatter Plot", 
       x = "X values", 
       y = "Y values")

# Save plot if output directory exists
if (dir.exists("/app/output")) {
  ggsave("/app/output/sample_plot.png", p, width = 8, height = 6)
  cat("Plot saved to output/sample_plot.png\n")
} else {
  # Just display plot info
  print(p)
  cat("Plot created (output directory not mounted)\n")
}

cat("Script completed successfully!\n")
