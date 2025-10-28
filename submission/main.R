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
                   
# =========================
# OOS setup (80/20)
# =========================

# Working dataset for OOS: keep only needed columns and drop NAs
needed_cols <- c("date", vars_Y, vars_X_filtered)
needed_cols <- intersect(needed_cols, names(df_train_clean))  # safety
df_oos <- df_train_clean[, needed_cols, drop = FALSE]
df_oos <- df_oos[stats::complete.cases(df_oos), , drop = FALSE]

# Sort by date and reindex
df_oos <- df_oos[order(df_oos$date), , drop = FALSE]
row.names(df_oos) <- NULL

# Sanity checks
stopifnot(length(vars_Y) >= 1)
stopifnot(exists("vars_X_filtered"), length(vars_X_filtered) >= 2)

# Dimensions and 80/20 split
T_total <- nrow(df_oos)
stopifnot(T_total >= 40)  # minimal reasonable size for stability

t0_idx <- floor(0.8 * T_total)      # end index of initial training
h_long <- 4                         # 1-year ahead (quarterly)
last_origin_idx <- T_total - h_long # last valid origin for h=4

# If 80% falls too close to the end, adjust to guarantee room for h=4
if (t0_idx > last_origin_idx) {
  t0_idx <- max(
    floor(0.75 * T_total),                 # try 75%
    min(T_total - (h_long + 8),            # leave some buffer
        T_total - h_long - 1)
  )
}

# Final guardrails
t0_idx <- max(1, t0_idx)
stopifnot(last_origin_idx >= t0_idx)

# Expanding-window origins: t = t0, t0+1, ..., T-4
oos_origin_idx <- seq(t0_idx, last_origin_idx)
stopifnot(length(oos_origin_idx) > 0)

# Useful dates (origins and targets)
origin_dates    <- df_oos$date[oos_origin_idx]
target_h1_dates <- df_oos$date[oos_origin_idx + 1]
target_h4_dates <- df_oos$date[oos_origin_idx + h_long]

# Quick report
cat("\n=== OOS setup (expanding) ===\n")
cat("Usable observations (post-NA):  ", T_total, "\n")
cat("Initial training up to idx:     ", t0_idx, " => ", as.character(df_oos$date[t0_idx]), "\n")
cat("Last OOS origin (idx):          ", last_origin_idx, " => ", as.character(df_oos$date[last_origin_idx]), "\n")
cat("Number of OOS origins:          ", length(oos_origin_idx), "\n")
cat("First origin -> targets h1/h4:  ", as.character(origin_dates[1]), " -> ",
    as.character(target_h1_dates[1]), " / ", as.character(target_h4_dates[1]), "\n")
cat("Last origin -> targets h1/h4:   ", as.character(origin_dates[length(origin_dates)]), " -> ",
    as.character(target_h1_dates[length(target_h1_dates)]), " / ", as.character(target_h4_dates[length(target_h4_dates)]), "\n")

# Store objects for next steps 
OOS_SETUP <- list(
  df_oos = df_oos,
  T_total = T_total,
  t0_idx = t0_idx,
  h_long = h_long,
  oos_origin_idx = oos_origin_idx,
  origin_dates = origin_dates,
  target_h1_dates = target_h1_dates,
  target_h4_dates = target_h4_dates,
  vars_Y = vars_Y,
  vars_X_filtered = vars_X_filtered
)


# =========================
# STEP 2 — Expanding OOS: PCA+VAR on factors, OLS mapping Y|F
# =========================

# Unpack setup
df_oos           <- OOS_SETUP$df_oos
t0_idx           <- OOS_SETUP$t0_idx
h_long           <- OOS_SETUP$h_long
oos_origin_idx   <- OOS_SETUP$oos_origin_idx
vars_Y           <- OOS_SETUP$vars_Y
# vars_X_filtered  <- OOS_SETUP$vars_X_filtered   # <- replaced by vars_X_oos below

# ---------- (1) Freeze predictors by filtering collinearity ONLY on initial training (≤ t0_idx) ----------
all_X_names <- setdiff(colnames(df_oos), c("date", vars_Y))
X0_for_filter <- as.matrix(df_oos[1:t0_idx, all_X_names, drop = FALSE])

cor_mat0   <- cor(X0_for_filter, use = "pairwise.complete.obs")
threshold  <- 0.95
hi_pairs0  <- which(abs(cor_mat0) > threshold & upper.tri(cor_mat0), arr.ind = TRUE)
redundant0 <- if (nrow(hi_pairs0) > 0) unique(rownames(cor_mat0)[hi_pairs0[, 1]]) else character(0)
vars_X_oos <- setdiff(all_X_names, redundant0)

# Safety: ensure we have at least 2 predictors for PCA
stopifnot(length(vars_X_oos) >= 2)

# ---------- Fix r using Bai–Ng on the initial training (≤ t0_idx) ----------
# Build initial training matrices
X0 <- as.matrix(df_oos[1:t0_idx, vars_X_oos, drop = FALSE])
Y0 <- as.matrix(df_oos[1:t0_idx, vars_Y,     drop = FALSE])

# Standardize X (and Y if needed) on initial training to compute r (Bai–Ng)
X0_mean <- apply(X0, 2, mean, na.rm = TRUE)
X0_sd   <- apply(X0, 2, sd,   na.rm = TRUE); X0_sd[X0_sd == 0] <- 1
X0_std  <- scale(X0, center = X0_mean, scale = X0_sd)

# Choose r once on the initial training window
bn0 <- bn_select_r(X0_std, crit = "ICp3")  # same selector as before
r_fixed <- max(1, bn0$r)
message("Fixed number of factors r (ICp3) on initial training = ", r_fixed)

# ---------- Prepare storage objects ----------
nY   <- length(vars_Y)
nOrg <- length(oos_origin_idx)

# Standardized forecasts and reals (for metrics)
pred_h1_std <- matrix(NA, nrow = nOrg, ncol = nY, dimnames = list(NULL, vars_Y))
pred_h4_std <- matrix(NA, nrow = nOrg, ncol = nY, dimnames = list(NULL, vars_Y))
real_h1_std <- matrix(NA, nrow = nOrg, ncol = nY, dimnames = list(NULL, vars_Y))
real_h4_std <- matrix(NA, nrow = nOrg, ncol = nY, dimnames = list(NULL, vars_Y))

# Real-scale forecasts and reals (destandardized per-iteration, handy for reporting)
pred_h1_real <- matrix(NA, nrow = nOrg, ncol = nY, dimnames = list(NULL, vars_Y))
pred_h4_real <- matrix(NA, nrow = nOrg, ncol = nY, dimnames = list(NULL, vars_Y))
real_h1_real <- matrix(NA, nrow = nOrg, ncol = nY, dimnames = list(NULL, vars_Y))
real_h4_real <- matrix(NA, nrow = nOrg, ncol = nY, dimnames = list(NULL, vars_Y))

origin_dates <- df_oos$date[oos_origin_idx]
h1_dates     <- df_oos$date[oos_origin_idx + 1]
h4_dates     <- df_oos$date[oos_origin_idx + h_long]

# ---------- Expanding-window loop ----------
p_fix <- 2  # (2) fixed VAR lag; set to 1 or 2
k <- 0
for (tcut in oos_origin_idx) {
  k <- k + 1
  
  # Training data ≤ tcut
  X_train <- as.matrix(df_oos[1:tcut, vars_X_oos, drop = FALSE])
  Y_train <- as.matrix(df_oos[1:tcut, vars_Y,     drop = FALSE])
  
  # Standardize using training-only moments (real-time)
  X_mean <- apply(X_train, 2, mean, na.rm = TRUE)
  X_sd   <- apply(X_train, 2, sd,   na.rm = TRUE); X_sd[X_sd == 0] <- 1
  X_std  <- scale(X_train, center = X_mean, scale = X_sd)
  
  Y_mean <- apply(Y_train, 2, mean, na.rm = TRUE)
  Y_sd   <- apply(Y_train, 2, sd,   na.rm = TRUE); Y_sd[Y_sd == 0] <- 1
  Y_std  <- scale(Y_train, center = Y_mean, scale = Y_sd)
  
  # PCA on standardized X (re-estimated each origin), keep r_fixed
  pca_fit <- prcomp(X_std, center = FALSE, scale. = FALSE)
  F_train <- pca_fit$x[, 1:r_fixed, drop = FALSE]   # factors (T_train x r)
  colnames(F_train) <- paste0("F", seq_len(r_fixed))
  
  # ---------- (3) OLS mapping with factor lags: Y_std ~ [1, F_t, F_{t-1}] ----------
  L_ols <- 1
  Tn <- nrow(F_train); r <- ncol(F_train)
  stopifnot(Tn > L_ols)
  
  # Build regressor matrix aligned with Y
  F_curr  <- F_train[(1+L_ols):Tn, , drop = FALSE]     # F_t
  F_lag1  <- F_train[(1):(Tn-L_ols), , drop = FALSE]   # F_{t-1}
  X_ols   <- cbind(F_curr, F_lag1)                     # [F_t | F_{t-1}]
  colnames(X_ols) <- c(paste0("F",1:r), paste0("L1_F",1:r))
  
  Y_ols   <- Y_std[(1+L_ols):Tn, , drop = FALSE]       # align Y
  X_ols_const <- cbind(const = 1, X_ols)
  
  beta_hat <- array(NA, dim = c(ncol(X_ols_const), ncol(Y_ols)),
                    dimnames = list(colnames(X_ols_const), colnames(Y_ols)))
  for (j in seq_len(ncol(Y_ols))) {
    fit_j <- lm(Y_ols[, j] ~ X_ols)  # includes intercept
    beta_hat[, j] <- coef(fit_j)
  }
  
  # VAR on factors (levels), fixed lag p_fix
  var_fit <- VAR(as.data.frame(F_train), p = p_fix, type = "const")
  
  # Forecast factors iteratively to h=4
  fc <- predict(var_fit, n.ahead = h_long)
  fac_names <- colnames(F_train)
  F_h1 <- sapply(fac_names, function(nm) fc$fcst[[nm]][1,        "fcst"])  # F_{t+1}
  F_h3 <- sapply(fac_names, function(nm) fc$fcst[[nm]][h_long-1, "fcst"])  # F_{t+3}
  F_h4 <- sapply(fac_names, function(nm) fc$fcst[[nm]][h_long,   "fcst"])  # F_{t+4}
  
  # Map factor forecasts to Y (standardized scale) using the same lag structure
  # h = 1: regressors = [1, F_{t+1}, F_t]
  F_t <- as.numeric(F_train[nrow(F_train), , drop = FALSE])      # last observed factors at t
  x_h1_const <- c(1, c(F_h1, F_t))
  pred_h1_std[k, ] <- as.numeric(x_h1_const %*% beta_hat)
  
  # h = 4: regressors = [1, F_{t+4}, F_{t+3}]
  x_h4_const <- c(1, c(F_h4, F_h3))
  pred_h4_std[k, ] <- as.numeric(x_h4_const %*% beta_hat)
  
  # Collect real Y at t+1 and t+4 (in transformed units)
  Y_h1_real_scale <- as.numeric(df_oos[tcut + 1,       vars_Y, drop = TRUE])
  Y_h4_real_scale <- as.numeric(df_oos[tcut + h_long,  vars_Y, drop = TRUE])
  
  # Convert reals to standardized scale of this iteration (coherent metrics)
  real_h1_std[k, ] <- (Y_h1_real_scale - Y_mean) / Y_sd
  real_h4_std[k, ] <- (Y_h4_real_scale - Y_mean) / Y_sd
  
  # De-standardize predictions back to "transformed" units
  pred_h1_real[k, ] <- pred_h1_std[k, ] * Y_sd + Y_mean
  pred_h4_real[k, ] <- pred_h4_std[k, ] * Y_sd + Y_mean
  real_h1_real[k, ] <- Y_h1_real_scale
  real_h4_real[k, ] <- Y_h4_real_scale
}

# ---------- OOS metrics ----------
rmse <- function(e) sqrt(mean(e^2, na.rm = TRUE))
mae  <- function(e) mean(abs(e), na.rm = TRUE)

# Metrics on standardized scale
rmse_h1_std <- sapply(seq_len(ncol(real_h1_std)), function(j) rmse(real_h1_std[, j] - pred_h1_std[, j]))
rmse_h4_std <- sapply(seq_len(ncol(real_h4_std)), function(j) rmse(real_h4_std[, j] - pred_h4_std[, j]))
mae_h1_std  <- sapply(seq_len(ncol(real_h1_std)), function(j) mae(real_h1_std[, j] - pred_h1_std[, j]))
mae_h4_std  <- sapply(seq_len(ncol(real_h4_std)), function(j) mae(real_h4_std[, j] - pred_h4_std[, j]))

metrics_std <- data.frame(
  variable = vars_Y,
  RMSE_h1_std = rmse_h1_std,
  RMSE_h4_std = rmse_h4_std,
  MAE_h1_std  = mae_h1_std,
  MAE_h4_std  = mae_h4_std
)
print(metrics_std)

# Metrics on "transformed" real scale
rmse_h1_real <- sapply(seq_len(ncol(real_h1_real)), function(j) rmse(real_h1_real[, j] - pred_h1_real[, j]))
rmse_h4_real <- sapply(seq_len(ncol(real_h4_real)), function(j) rmse(real_h4_real[, j] - pred_h4_real[, j]))
mae_h1_real  <- sapply(seq_len(ncol(real_h1_real)), function(j) mae(real_h1_real[, j] - pred_h1_real[, j]))
mae_h4_real  <- sapply(seq_len(ncol(real_h4_real)), function(j) mae(real_h4_real[, j] - pred_h4_real[, j]))

metrics_real <- data.frame(
  variable = vars_Y,
  RMSE_h1_real = rmse_h1_real,
  RMSE_h4_real = rmse_h4_real,
  MAE_h1_real  = mae_h1_real,
  MAE_h4_real  = mae_h4_real
)
print(metrics_real)

# ---------- Collect a tidy OOS predictions table ----------
oos_results <- data.frame(
  origin_date = origin_dates,
  h1_date = h1_dates,
  h4_date = h4_dates
)
# attach predictions and actuals (real scale)
oos_results <- cbind(
  oos_results,
  setNames(as.data.frame(pred_h1_real), paste0(vars_Y, "_pred_h1")),
  setNames(as.data.frame(real_h1_real), paste0(vars_Y, "_real_h1")),
  setNames(as.data.frame(pred_h4_real), paste0(vars_Y, "_pred_h4")),
  setNames(as.data.frame(real_h4_real), paste0(vars_Y, "_real_h4"))
)

# Quick preview
head(oos_results, 5)

# =========================
# STEP 3 — Plot OOS results (h = 1 and h = 4)
# =========================

library(ggplot2)
library(dplyr)
library(tidyr)

# --- Prepare long-format data for ggplot ---

# Select only relevant columns (real + pred for both horizons)
plot_df <- oos_results %>%
  dplyr::select(origin_date,
                dplyr::ends_with("_pred_h1"),
                dplyr::ends_with("_real_h1"),
                dplyr::ends_with("_pred_h4"),
                dplyr::ends_with("_real_h4")) %>%
  dplyr::rename(date = origin_date)

# Reshape to long format for each horizon
long_h1 <- plot_df %>%
  dplyr::select(date, dplyr::ends_with("_pred_h1"), dplyr::ends_with("_real_h1")) %>%
  tidyr::pivot_longer(
    cols = -date,
    names_to = c("variable", "type"),
    names_pattern = "(.*)_(pred|real)_h1",
    values_to = "value"
  ) %>%
  dplyr::mutate(horizon = "h=1 (1Q ahead)")

long_h4 <- plot_df %>%
  dplyr::select(date, dplyr::ends_with("_pred_h4"), dplyr::ends_with("_real_h4")) %>%
  tidyr::pivot_longer(
    cols = -date,
    names_to = c("variable", "type"),
    names_pattern = "(.*)_(pred|real)_h4",
    values_to = "value"
  ) %>%
  dplyr::mutate(horizon = "h=4 (1Y ahead)")


long_all <- bind_rows(long_h1, long_h4)

# Clean variable names for facets
long_all$variable <- recode(long_all$variable,
                            gdp = "GDP",
                            cpi = "CPI",
                            unempoff = "Unemployment"
)

# --- Plot OOS forecasts vs actuals ---
ggplot(long_all, aes(x = date, y = value, color = type, linetype = type)) +
  geom_line(linewidth = 0.7) +
  facet_grid(variable ~ horizon, scales = "free_y") +
  scale_color_manual(values = c("real" = "black", "pred" = "steelblue")) +
  scale_linetype_manual(values = c("real" = "solid", "pred" = "dashed")) +
  labs(
    title = "Out-of-Sample Forecasts (Expanding Window)",
    subtitle = "Dynamic Factor Model (PCA + VAR)",
    x = "Date", y = "Value",
    color = "", linetype = ""
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14)
  )

# =========================
# "NOW" FORECASTS from latest data (h=1 and h=4)
# =========================

# Keep only observed data up to today (adjust the date if needed)
df_train_clean <- df_train_clean %>%
  dplyr::filter(date <= as.Date("2025-07-01"))  # last observed quarter

# Decide predictor set for NOW: prefer vars_X_oos if available (frozen on initial 80%)
if (exists("vars_X_oos")) {
  X_now_names <- vars_X_oos
  message("NOW: using frozen predictor set vars_X_oos (from initial 80%).")
} else {
  X_now_names <- vars_X_filtered
  message("NOW: vars_X_oos not found, falling back to vars_X_filtered.")
}

# Build working frame
needed_cols <- c("date", vars_Y, X_now_names)
needed_cols <- intersect(needed_cols, names(df_train_clean))
df_now <- df_train_clean[, needed_cols, drop = FALSE]
df_now <- df_now[stats::complete.cases(df_now), , drop = FALSE]
df_now <- df_now[order(df_now$date), , drop = FALSE]
row.names(df_now) <- NULL

# Split into X and Y (full history up to last observed quarter)
X_all <- as.matrix(df_now[, X_now_names, drop = FALSE])
Y_all <- as.matrix(df_now[, vars_Y,      drop = FALSE])

# Standardize using full-history moments (we are producing NOW forecasts)
X_mean <- apply(X_all, 2, mean, na.rm = TRUE)
X_sd   <- apply(X_all, 2, sd,   na.rm = TRUE); X_sd[X_sd == 0] <- 1
X_std  <- scale(X_all, center = X_mean, scale = X_sd)

Y_mean <- apply(Y_all, 2, mean, na.rm = TRUE)
Y_sd   <- apply(Y_all, 2, sd,   na.rm = TRUE); Y_sd[Y_sd == 0] <- 1
Y_std  <- scale(Y_all, center = Y_mean, scale = Y_sd)

# Select number of factors r (reuse r_fixed if available; else choose on initial 80%)
if (exists("r_fixed")) {
  r_now <- r_fixed
  message("NOW: using r_fixed = ", r_now)
} else {
  t0_now <- max(1L, floor(0.8 * nrow(X_std)))
  X0_now <- X_std[1:t0_now, , drop = FALSE]
  bn_now <- bn_select_r(X0_now, crit = "ICp3")
  r_now  <- max(1, bn_now$r)
  message("NOW: r selected by ICp3 on initial 80% = ", r_now)
}

# PCA on standardized X (full history)
pca_now <- prcomp(X_std, center = FALSE, scale. = FALSE)
F_all   <- pca_now$x[, 1:r_now, drop = FALSE]
colnames(F_all) <- paste0("F", seq_len(r_now))

# OLS with factor lags (L=1): Y_std ~ [1, F_t, F_{t-1}]
L_ols <- 1
Tn <- nrow(F_all); r <- ncol(F_all)
stopifnot(Tn > L_ols)

F_curr <- F_all[(1+L_ols):Tn, , drop = FALSE]        # F_t
F_lag1 <- F_all[(1):(Tn-L_ols), , drop = FALSE]      # F_{t-1}
X_ols  <- cbind(F_curr, F_lag1)                      # [F_t | F_{t-1}]
colnames(X_ols) <- c(paste0("F",1:r), paste0("L1_F",1:r))
Y_ols  <- Y_std[(1+L_ols):Tn, , drop = FALSE]        # align Y

X_ols_const <- cbind(const = 1, X_ols)
beta_hat_now <- array(NA, dim = c(ncol(X_ols_const), ncol(Y_ols)),
                      dimnames = list(colnames(X_ols_const), colnames(Y_ols)))
for (j in seq_len(ncol(Y_ols))) {
  fit_j <- lm(Y_ols[, j] ~ X_ols)  # includes intercept
  beta_hat_now[, j] <- coef(fit_j)
}

# VAR on factors (levels) with fixed lag p=2 (more stable)
p_fix <- 2
var_now <- VAR(as.data.frame(F_all), p = p_fix, type = "const")

# Forecast factors to h=1 and h=4
h_long <- 4
fc_now <- predict(var_now, n.ahead = h_long)
fac_names <- colnames(F_all)

F_h1 <- sapply(fac_names, function(nm) fc_now$fcst[[nm]][1, "fcst"])  # F_{t+1}
F_h3 <- sapply(fac_names, function(nm) fc_now$fcst[[nm]][3, "fcst"])  # F_{t+3}
F_h4 <- sapply(fac_names, function(nm) fc_now$fcst[[nm]][4, "fcst"])  # F_{t+4}
F_t  <- as.numeric(F_all[nrow(F_all), , drop = FALSE])                # F_t

# Map factor forecasts to Y (std) using same lag structure
# h=1: regressors = [1, F_{t+1}, F_t]
Y_h1_std <- as.numeric(c(1, c(F_h1, F_t)) %*% beta_hat_now)
# h=4: regressors = [1, F_{t+4}, F_{t+3}]
Y_h4_std <- as.numeric(c(1, c(F_h4, F_h3)) %*% beta_hat_now)

# De-standardize back to transformed units
Y_h1_real <- Y_h1_std * Y_sd + Y_mean
Y_h4_real <- Y_h4_std * Y_sd + Y_mean
names(Y_h1_real) <- vars_Y
names(Y_h4_real) <- vars_Y

# Target dates (next quarter and +4 quarters)
last_date  <- max(df_now$date)
h1_date    <- seq(last_date, by = "quarter", length.out = 2)[2]
h4_date    <- seq(last_date, by = "quarter", length.out = 5)[5]

now_table <- data.frame(
  variable = vars_Y,
  target_h1_date = h1_date,
  forecast_h1    = as.numeric(Y_h1_real),
  target_h4_date = h4_date,
  forecast_h4    = as.numeric(Y_h4_real)
)
print(now_table)

