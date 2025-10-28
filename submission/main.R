# TODO: Is the cutoff at 0.99 correlation common practice? --> Now is set to 0.95!!!
# 0.99 --> 0.95: Log Likelihood: from -167.4 to -149.6, Max root: from 1.007 to 1.004
# TODO compare with their forecasts, reconstruction error there
# TODO review if this model setup is right --> tested with "summary(pca_X)" and "roots(var_model)"
# TODO fix the missing gap in the data --> tested with "tail(df_train_clean$date)" "head(future_dates)", "tail(Y_pred_real)" "head(forecast_long_all)"
# added Bai-Ng criterion
library(ggplot2)
library(dplyr)
library(vars) 
library(tidyr)   
library(urca)

df <- utils::read.csv("C:/Users/kfree/OneDrive/Desktop/MASTER 3/MACRO FORECAST/DATASET/data_quarterly.csv")
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

metadata <- utils::read.csv("C:/Users/kfree/OneDrive/Desktop/MASTER 3/MACRO FORECAST/DATASET/metadata_quarterly_en.csv")

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
metadata <- metadata[!is.na(metadata$transform), ]
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

# Test Stationarity after transformation
adf_results <- sapply(df_train_clean[-1], function(x) {
  test <- tryCatch(ur.df(na.omit(x), type = "drift", selectlags = "AIC"), error = function(e) NULL)
  if (is.null(test)) return(NA)
  pval <- test@testreg$coefficients[2,4]
  return(pval)
})

adf_results <- sort(adf_results)
cat("P-value ADF test (first 10 variables less stationary):\n")
print(head(adf_results, 10))

# Delete Varibales with Variance 0
var_zero <- sapply(df_train_clean[-1], function(x) var(x, na.rm = TRUE) == 0)
if (any(var_zero)) {
  cat("Deleted", sum(var_zero), "variabiles with 0 variance:\n")
  print(names(var_zero[var_zero]))
  df_train_clean <- df_train_clean[, c(TRUE, !var_zero)]  # keep column "Date"
}

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

# NEW: Bai–Ng selectors 
bn_select_r <- function(X, rmax = NULL, crit = c("ICp3","ICp2","ICp1","PCp1","PCp2","PCp3")){
  crit <- match.arg(crit)
  X <- as.matrix(X)
  Tn <- nrow(X); Nn <- ncol(X)
  
  # NEW: prudent upper bound for r
  if (is.null(rmax)) {
    rmax <- min(15, Nn - 1, floor(0.2 * min(Nn, Tn)))
  }
  
  # center columns 
  Xc <- scale(X, center = TRUE, scale = FALSE)
  
  # SVD once
  sv <- svd(Xc)                # Xc = U D V'
  d2 <- sv$d^2                 # singular values squared
  
  # Residual SSR for r = 0,1,... (tail sum of d^2)
  tail_ss <- rev(cumsum(rev(d2)))         # length = min(N,T)
  rgrid   <- 0:min(rmax, length(d2) - 1)
  SSR_vec <- tail_ss[rgrid + 1]
  V_r     <- SSR_vec / (Nn * Tn)          # average residual variance
  
  # Penalties 
  C_NT  <- (Nn + Tn) / (Nn * Tn)
  L1    <- log(min(Nn, Tn))
  L2    <- log(Nn * Tn / (Nn + Tn))
  L3    <- log(Nn * Tn)
  
  # Information criteria (ICp*) and panel criteria (PCp*)
  ICp1 <- log(V_r) + rgrid * C_NT * L1
  ICp2 <- log(V_r) + rgrid * C_NT * L2
  ICp3 <- log(V_r) + rgrid * C_NT * L3
  
  PCp1 <- V_r + rgrid * C_NT * L1
  PCp2 <- V_r + rgrid * C_NT * L2
  PCp3 <- V_r + rgrid * C_NT * L3
  
  all_IC <- switch(crit,
                   ICp1 = ICp1, ICp2 = ICp2, ICp3 = ICp3,
                   PCp1 = PCp1, PCp2 = PCp2, PCp3 = PCp3)
  
  r_sel <- rgrid[which.min(all_IC)]
  list(r = r_sel, rgrid = rgrid,
       IC = data.frame(r = rgrid, ICp1 = ICp1, ICp2 = ICp2, ICp3 = ICp3,
                       PCp1 = PCp1, PCp2 = PCp2, PCp3 = PCp3))
}

Xmat <- df_train_std[, vars_X_filtered, drop = FALSE]

# NEW: stronger penalty by default (ICp3). 
bn   <- bn_select_r(Xmat, crit = "ICp3")   # try "ICp2" / "PCp2" to compare
r    <- max(1, bn$r)
message("Bai–Ng (", "ICp3", ") selected r = ", r)

# run PCA once, keep first r
pca_X <- prcomp(Xmat, center = FALSE, scale. = FALSE)
F_hat <- pca_X$x[, 1:r, drop = FALSE]
Lambda_hat <- pca_X$rotation[, 1:r, drop = FALSE]
F_df <- as.data.frame(F_hat)

# NEW: quick diagnostics 
print(head(bn$IC, 12))
# plot(bn$IC$r, bn$IC$ICp3, type="b", main="ICp3 vs r")

# Test Stationarity for each factor
adf_factors <- sapply(F_df, function(x) {
  test <- tryCatch(ur.df(na.omit(x), type = "drift", selectlags = "AIC"), error = function(e) NULL)
  if (is.null(test)) return(NA)
  pval <- test@testreg$coefficients[2,4]
  return(pval)
})

adf_factors <- sort(adf_factors)
cat("P-value ADF test for PCA factors:\n")
print(adf_factors)


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

