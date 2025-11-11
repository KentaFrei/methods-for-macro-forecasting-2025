# ================================
# 0. Initialize renv
# ================================
source("renv/activate.R")

# ================================
# 1. Setup, Load Packages and Data
# ================================

library(ggplot2)
library(dplyr)
library(vars)
library(tidyr)
library(urca)
library(rlang)
library(KFAS)


kfilepath_data <- "C:/Users/kfree/OneDrive/Desktop/MASTER 3/MACRO FORECAST/DATASET/data_quarterly.csv"
kfilepath_metadata <- "C:/Users/kfree/OneDrive/Desktop/MASTER 3/MACRO FORECAST/DATASET/metadata_quarterly_en.csv"
nfilepath_data <- "~/Documents/school/methods-for-macro-forecasting-2025/submission/data/data_quarterly.csv"
nfilepath_metadata <- "~/Documents/school/methods-for-macro-forecasting-2025/submission/data/metadata_quarterly_en.csv"
# Adjust depending on environment
df <- utils::read.csv(nfilepath_data)
metadata <- utils::read.csv(nfilepath_metadata)
cutoff_pre_covid <- as.Date("2018-12-31")
# ================================
# 2. Data Preparation
# ================================
names(metadata) <- tolower(names(metadata))
if (!"ts_key" %in% names(metadata)) stop("Missing ts_key in metadata")

# Filter to data we know (not projections)
df_train <- df %>%
  mutate(date = as.Date(paste0(date, "-01"))) %>%        
  filter(date <= cutoff_pre_covid) %>%
  mutate(across(-date, as.numeric))

df_post <- df %>%
  mutate(date = as.Date(paste0(date, "-01"))) %>%
  filter(date > cutoff_pre_covid) %>%
  mutate(across(-date, as.numeric))

# Check to confirm no missing values
sapply(df_train, function(x) sum(is.na(x))) %>% sort(decreasing = TRUE)

# normalize caps
metadata$unit_low <- tolower(metadata$unit)

# We need to transform the raw values to make the series stationary
unit_to_transform <- c(
  # Real and nominal quantities (logs capture proportional growth)
  "real in millions of francs" = "logdiff",
  "in millions of francs at current prices" = "logdiff",
  "real in millions of us dollars" = "logdiff",
  "in us dollars" = "logdiff",
  "real in thousands of swiss francs" = "logdiff",
  "in thousands of swiss francs at current prices" = "logdiff",
  
  # Indexes and deflators (inflation-type measures)
  "index" = "logdiff",
  "price deflator" = "logdiff",
  
  # Productivity measures (real per-hour or per-FTE)
  "gross domestic product (adjusted for sporting events) per hour worked" = "logdiff",
  "gross domestic product (adjusted for sporting events) per full-time equivalent" = "logdiff",
  "real value of gross domestic product per full-time equivalent" = "logdiff",
  
  # Labor quantities (use diff if absolute, logdiff if proportional)
  "in thousands of persons" = "diff",
  
  # Rates and percentages (stationary in differences or demeaned)
  "in percent" = "diff",
  "rate" = "diff"
)
metadata$transform <- unit_to_transform[metadata$unit_low]

df_train_clean <- df_train

transform_variable <- function(x, transform_type) {
  if (transform_type == "logdiff") {
    y <- rep(NA_real_, length(x))
    pos <- x > 0
    y[pos] <- log(x[pos])             # niente warning
    return(c(NA_real_, diff(y) * 100))
  } else if (transform_type == "diff") {
    return(c(NA_real_, diff(x)))
  } else {
    return(x)
  }
}

# Apply transformations using ts_key 
for (i in seq_len(nrow(metadata))) {
  key <- metadata$ts_key[i]
  t_type <- metadata$transform[i]
  if (!is.na(t_type) && key %in% names(df_train_clean)) {
    df_train_clean[[key]] <- transform_variable(df_train_clean[[key]], t_type)
  }
}

# Remove first raw (diff) e drop problematic columns (all NA or 0 variance)
df_train_clean <- df_train_clean[-1, ]
is_all_na  <- vapply(df_train_clean, function(x) all(is.na(x)), logical(1))
df_train_clean <- df_train_clean[, !is_all_na, drop = FALSE]
keep_idx <- c(TRUE, sapply(df_train_clean[-1], function(x) sd(x, na.rm=TRUE) > 1e-8))
df_train_clean <- df_train_clean[, keep_idx, drop = FALSE]

# Standardization
df_train_std <- df_train_clean
df_train_std[-1] <- scale(df_train_clean[-1])

# Define targets (ts_key) than X
vars_Y <- c("rvgdp", "cpi", "wkfreuro")
growth_vars <- c("rvgdp")  # stay as Δlog×100 (growth)
level_vars  <- setdiff(vars_Y, growth_vars)
missing_y <- setdiff(vars_Y, names(df_train_std))
if (length(missing_y)) stop(paste("Missing targets:", paste(missing_y, collapse=", ")))
vars_X <- setdiff(names(df_train_std)[-1], vars_Y)

# Test Stationarity after transformation
adf_results <- sapply(df_train_clean[-1], function(x) {
  z <- na.omit(x)
  if (length(z) < 20) return(NA_real_)
  test <- tryCatch(ur.df(z, type = "drift", selectlags = "AIC"), error = function(e) NULL)
  if (is.null(test)) return(NA_real_)
  pval <- tryCatch(test@testreg$coefficients[2,4], error = function(e) NA_real_)
  return(pval)
})
adf_results <- sort(adf_results)

# Saving means and sd 
scaler_params <- data.frame(
  variable = names(df_train_clean)[-1],
  mean = apply(df_train_clean[-1], 2, mean, na.rm = TRUE),
  sd = apply(df_train_clean[-1], 2, sd, na.rm = TRUE)
)

# Sort by date (safety)
df_train_std <- df_train_std[order(df_train_std$date), , drop = FALSE]

# Dimensions and 80/20 split
T_total <- nrow(df_train_std)
t0_idx <- floor(0.8 * T_total)      # end index of initial training 
h_long <- 4                         # 1-year ahead (quarterly)
last_origin_idx <- T_total - h_long # last valid origin for h=4
# If 80% falls too close to the end, adjust to guarantee room for h=4
if (t0_idx > last_origin_idx) {
  t0_idx <- max(
    floor(0.75 * T_total),
    min(T_total - (h_long + 8),
        T_total - h_long - 1)
  )
}

# Drop variables with <30% of NA from training set for Bai-Ng criterion
na_frac <- colMeans(is.na(as.matrix(df_train_std[1:t0_idx, vars_X, drop = FALSE])))
vars_X <- setdiff(vars_X, names(na_frac[na_frac > 0.30]))

# Final guardrails
t0_idx <- max(1, t0_idx)
stopifnot(last_origin_idx >= t0_idx)

# Expanding-window origins: t = t0, t0+1, ..., T-4
oos_origin_idx <- seq(t0_idx, last_origin_idx)
stopifnot(length(oos_origin_idx) > 0)

# Useful dates (origins and targets)
origin_dates    <- df_train_std$date[oos_origin_idx]
target_h1_dates <- df_train_std$date[oos_origin_idx + 1]
target_h4_dates <- df_train_std$date[oos_origin_idx + h_long]

# Quick report
cat("\n=== OOS setup (expanding) ===\n")
cat("Usable observations (post-NA):  ", T_total, "\n")
cat("Initial training up to idx:     ", t0_idx, " => ", as.character(df_train_std$date[t0_idx]), "\n")
cat("Last OOS origin (idx):          ", last_origin_idx, " => ", as.character(df_train_std$date[last_origin_idx]), "\n")
cat("Number of OOS origins:          ", length(oos_origin_idx), "\n")
cat("First origin -> targets h1/h4:  ", as.character(origin_dates[1]), " -> ",
    as.character(target_h1_dates[1]), " / ", as.character(target_h4_dates[1]), "\n")
cat("Last origin -> targets h1/h4:   ", as.character(origin_dates[length(origin_dates)]), " -> ",
    as.character(target_h1_dates[length(target_h1_dates)]), " / ", as.character(target_h4_dates[length(target_h4_dates)]), "\n")
# --- De-duplicate highly collinear series (cluster on |corr|)
X_train_for_sel <- as.matrix(df_train_std[1:t0_idx, vars_X, drop = FALSE])

# Guard: need at least 2 cols and some non-NA variance
ok <- apply(X_train_for_sel, 2, function(x) sd(x, na.rm = TRUE) > 0)
X_train_for_sel <- X_train_for_sel[, ok, drop = FALSE]
vars_X <- colnames(X_train_for_sel)

C  <- cor(X_train_for_sel, use = "pairwise.complete.obs")
d  <- as.dist(1 - abs(C))
hc <- hclust(d, method = "average")

# Cut height: 0.1 ≈ keep 1 from any group with |ρ| > 0.9
groups <- cutree(hc, h = 0.1)

# Keep the column in each cluster with the fewest NAs (or highest variance—your choice)
keepers <- tapply(colnames(X_train_for_sel), groups, function(cols) {
  nafrac <- colMeans(is.na(X_train_for_sel[, cols, drop = FALSE]))
  cols[which.min(nafrac)]
})

vars_X <- unname(unlist(keepers))
# Prepare initial training matrices
X0 <- as.matrix(df_train_std[1:t0_idx, vars_X, drop = FALSE])
Y0 <- as.matrix(df_train_std[1:t0_idx, vars_Y,     drop = FALSE])

# TESTS
sum(names(df_train_clean)[-1] != names(df_train)[-1])

# How many columns are transformed withing the remaining ones?
common <- intersect(names(df_train_clean)[-1], names(df_train)[-1])
changed <- sum(sapply(common, function(nm) !identical(df_train_clean[[nm]], df_train[[nm]][-1])))
cat("Serie effettivamente trasformate (tra quelle rimaste): ", changed, "\n")

# Which columns are deleted (es. all NA after logdiff or sd≈0)?
dropped <- setdiff(names(df_train)[-1], names(df_train_clean)[-1])
cat("Colonne rimosse: ", paste(dropped, collapse = ", "), "\n")

# Dimension X/Y for safety
cat("Dim X0: ", paste(dim(X0), collapse=" x "), " | Dim Y0: ", paste(dim(Y0), collapse=" x "), "\n")

# Are there NA in X0?
X0 <- as.matrix(df_train_std[1:t0_idx, vars_X, drop = FALSE])
sort(colSums(is.na(X0)), decreasing = TRUE)[1:10]
anyNA(X0)

# ================================
# 3. Bai-Ng factors selection
# ================================

# --- Definition before applying ---
bn_select_r_icp2 <- function(X, rmax = 10) {
  X <- as.matrix(X)
  Tn <- nrow(X); Nn <- ncol(X)
  rmax <- min(rmax, Nn - 1, floor(0.2 * min(Nn, Tn)))
  Xc <- scale(X, center = TRUE, scale = FALSE)
  sv <- svd(Xc)
  d2 <- sv$d^2
  tail_ss <- rev(cumsum(rev(d2)))
  rgrid <- 1:min(rmax, length(d2) - 1)
  SSR_vec <- tail_ss[rgrid + 1]
  V_r <- SSR_vec / (Nn * Tn)
  C_NT <- (Nn + Tn) / (Nn * Tn)
  Lmin <- log(min(Nn, Tn))
  ICp2 <- log(V_r) + rgrid * C_NT * Lmin
  r_sel <- rgrid[which.min(ICp2)]
  list(r = r_sel, IC = data.frame(r = rgrid, ICp2 = ICp2))
}

# --- Simple imputation ONLY for r selection ---
X0_imp <- X0
X0_imp[!is.finite(X0_imp)] <- NA
for (j in seq_len(ncol(X0_imp))) {
  m <- mean(X0_imp[, j], na.rm = TRUE); if (!is.finite(m)) m <- 0
  X0_imp[is.na(X0_imp[, j]), j] <- m
}
bn <- bn_select_r_icp2(X0_imp, rmax = 10)
r_fixed <- bn$r
message("Bai–Ng ICp2 (cap=10) selected r = ", r_fixed)

pca_res <- prcomp(X0_imp, center = TRUE, scale. = TRUE)

# Extract variance explained
eigvals <- pca_res$sdev^2
var_explained <- eigvals / sum(eigvals)
cumvar_explained <- cumsum(var_explained)

df_pca <- data.frame(
  Component = seq_along(eigvals),
  Variance = var_explained,
  Cumulative = cumvar_explained
)

# --- Step 3: Bai–Ng selection (if not yet run) ---
bn <- bn_select_r_icp2(X0_imp, rmax = 50)
r_selected <- bn$r

# --- Step 4: Plot cumulative variance explained ---
library(ggplot2)

ggplot(df_pca, aes(x = Component, y = Cumulative)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 2) +
  geom_vline(xintercept = r_selected, linetype = "dashed", color = "red", linewidth = 0.8) +
  annotate(
    "text", 
    x = r_selected + 0.3, 
    y = cumvar_explained[r_selected], 
    label = paste0("Bai–Ng selected r = ", r_selected),
    color = "red", 
    hjust = 0, 
    vjust = -0.5, 
    size = 4
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Cumulative Variance Explained by Principal Components",
    subtitle = "Dashed red line = Bai–Ng selected number of factors",
    x = "Number of Factors (Principal Components)",
    y = "Cumulative Variance Explained (%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )


# ================================
# 4. Forecasting with PCA + VAR (Dynamic Factor Model)
# ================================

nY <- length(vars_Y)
nOrg <- length(oos_origin_idx)
h_long <- 4
p_fix <- 2
L_ols <- 1

pred_h1_real <- pred_h4_real <- real_h1_real <- real_h4_real <-
  matrix(NA, nrow = nOrg, ncol = nY, dimnames = list(NULL, vars_Y))

# new matrices for standardised version (for metrics)
real_h1_std <- real_h4_std <- pred_h1_std_mat <- pred_h4_std_mat <-
  matrix(NA, nrow = nOrg, ncol = nY, dimnames = list(NULL, vars_Y))

origin_dates <- df_train_std$date[oos_origin_idx]
h1_dates <- df_train_std$date[oos_origin_idx + 1]
h4_dates <- df_train_std$date[oos_origin_idx + h_long]

cat("\n=== Expanding-window forecasting ===\n")

k <- 0
for (tcut in oos_origin_idx) {
  k <- k + 1
  cat("Origin", k, "/", nOrg, "— fino a", as.character(df_train_std$date[tcut]), "\n")
  
  X_train <- as.matrix(df_train_std[1:tcut, vars_X, drop = FALSE])
  Y_train <- as.matrix(df_train_std[1:tcut, vars_Y, drop = FALSE])
  
  X_mean <- colMeans(X_train, na.rm = TRUE)
  X_sd <- apply(X_train, 2, sd, na.rm = TRUE); X_sd[X_sd == 0] <- 1
  X_std <- scale(X_train, center = X_mean, scale = X_sd)
  
  for (j in seq_len(ncol(X_std))) {
    m <- mean(X_std[, j], na.rm = TRUE)
    if (!is.finite(m)) m <- 0
    X_std[is.na(X_std[, j]), j] <- m
  }
  
  Y_mean <- colMeans(Y_train, na.rm = TRUE)
  Y_sd <- apply(Y_train, 2, sd, na.rm = TRUE); Y_sd[Y_sd == 0] <- 1
  Y_std <- scale(Y_train, center = Y_mean, scale = Y_sd)
  
  pca_fit <- prcomp(X_std, center = FALSE, scale. = FALSE)
  F_train <- pca_fit$x[, 1:r_fixed, drop = FALSE]
  colnames(F_train) <- paste0("F", seq_len(r_fixed))
  
  Tn <- nrow(F_train)
  if (Tn <= L_ols + 1) next
  
  F_curr <- F_train[(1+L_ols):Tn, , drop = FALSE]
  F_lag1 <- F_train[1:(Tn-L_ols), , drop = FALSE]
  X_ols <- cbind(F_curr, F_lag1)
  Y_ols <- Y_std[(1+L_ols):Tn, , drop = FALSE]
  
  beta_hat <- sapply(seq_len(ncol(Y_ols)), function(j) coef(lm(Y_ols[, j] ~ X_ols)))
  rownames(beta_hat) <- c("const", colnames(X_ols))
  colnames(beta_hat) <- colnames(Y_ols)
  
  var_fit <- VAR(as.data.frame(F_train), p = p_fix, type = "const")
  fc <- predict(var_fit, n.ahead = h_long)
  fac_names <- colnames(F_train)
  F_h1 <- sapply(fac_names, function(nm) fc$fcst[[nm]][1, "fcst"])
  F_h3 <- sapply(fac_names, function(nm) fc$fcst[[nm]][h_long-1, "fcst"])
  F_h4 <- sapply(fac_names, function(nm) fc$fcst[[nm]][h_long, "fcst"])
  
  F_t <- as.numeric(F_train[nrow(F_train), ])
  x_h1_const <- c(1, F_h1, F_t)
  x_h4_const <- c(1, F_h4, F_h3)
  
  # standardised scale forecasts (window)
  pred_h1_std <- as.numeric(x_h1_const %*% beta_hat)
  pred_h4_std <- as.numeric(x_h4_const %*% beta_hat)
  pred_h1_std_mat[k, ] <- pred_h1_std
  pred_h4_std_mat[k, ] <- pred_h4_std
  
  # return to real scala (diff/logdiff of window)
  pred_h1_real[k, ] <- pred_h1_std * Y_sd + Y_mean
  pred_h4_real[k, ] <- pred_h4_std * Y_sd + Y_mean
  
  # actuals at real scale (for pred_h*_real)
  real_h1_real[k, ] <- as.numeric(df_train_clean[tcut + 1, vars_Y])
  real_h4_real[k, ] <- as.numeric(df_train_clean[tcut + h_long, vars_Y])
  
  # actuals in standardized scale  (same mean/sd of window, for metrics)
  real_h1_std[k, ] <- (as.numeric(df_train_clean[tcut + 1, vars_Y]) - Y_mean) / Y_sd
  real_h4_std[k, ] <- (as.numeric(df_train_clean[tcut + h_long, vars_Y]) - Y_mean) / Y_sd
}

cat("\n=== Forecast loop completed ===\n")


# ================================
# 5) OOS metrics + Back-transformation to levels
# ================================

# ---- Safety checks (avoid cryptic errors) ----
req_objs <- c("vars_Y","metadata","df_train",
              "origin_dates","h1_dates","h4_dates","oos_origin_idx","h_long",
              "pred_h1_real","pred_h4_real","real_h1_real","real_h4_real")
miss <- req_objs[!vapply(req_objs, exists, logical(1))]
if (length(miss)) stop("Those objects are missing: ", paste(miss, collapse=", "))

# ---- Helper: RMSE / MAE ----
rmse <- function(e) sqrt(mean(e^2, na.rm = TRUE))
mae  <- function(e) mean(abs(e), na.rm = TRUE)

# ---- Metrics in transformed scale ----
rmse_h1_tr <- sapply(seq_len(ncol(real_h1_real)), function(j) rmse(real_h1_real[, j] - pred_h1_real[, j]))
rmse_h4_tr <- sapply(seq_len(ncol(real_h4_real)), function(j) rmse(real_h4_real[, j] - pred_h4_real[, j]))
mae_h1_tr  <- sapply(seq_len(ncol(real_h1_real)), function(j) mae(real_h1_real[, j] - pred_h1_real[, j]))
mae_h4_tr  <- sapply(seq_len(ncol(real_h4_real)), function(j) mae(real_h4_real[, j] - pred_h4_real[, j]))

metrics_transformed <- data.frame(
  variable = colnames(real_h1_real),
  RMSE_h1  = rmse_h1_tr,
  RMSE_h4  = rmse_h4_tr,
  MAE_h1   = mae_h1_tr,
  MAE_h4   = mae_h4_tr,
  row.names = NULL
)
cat("\n=== Forecast Evaluation (Transformed Units) ===\n")
print(metrics_transformed)

# ---- Back-transform to levels (h=1 e h=4) ----
get_transform <- function(var) {
  i <- match(var, metadata$ts_key)
  tr <- if (is.na(i)) NA_character_ else metadata$transform[i]
  ifelse(is.na(tr), "none", tolower(tr))
}

# historical levels (non transformed)
levels_df <- df_train[, c("date", vars_Y)]

# Pre-allocation of containers (lists -> data.frame)
pred_h1_level <- setNames(vector("list", length(vars_Y)), vars_Y)
pred_h4_level <- setNames(vector("list", length(vars_Y)), vars_Y)
real_h1_level <- setNames(vector("list", length(vars_Y)), vars_Y)
real_h4_level <- setNames(vector("list", length(vars_Y)), vars_Y)

for (v in vars_Y) {
  tr_v <- get_transform(v)
  base_h1 <- levels_df[[v]][oos_origin_idx]
  base_h4 <- levels_df[[v]][oos_origin_idx + (h_long - 1)]
  L_t1 <- levels_df[[v]][oos_origin_idx + 1]
  L_t4 <- levels_df[[v]][oos_origin_idx + h_long]
  yhat_h1 <- pred_h1_real[, v]
  yhat_h4 <- pred_h4_real[, v]
  
  if (v %in% growth_vars) {
    # Keep in Δlog×100 (QoQ %)
    pred_h1_level[[v]] <- yhat_h1
    pred_h4_level[[v]] <- yhat_h4
    real_h1_level[[v]] <- real_h1_real[, v]
    real_h4_level[[v]] <- real_h4_real[, v]
  } else if (tr_v == "logdiff") {
    # Back to levels (index / rate)
    pred_h1_level[[v]] <- base_h1 * exp(yhat_h1 / 100)
    pred_h4_level[[v]] <- base_h4 * exp(yhat_h4 / 100)
    real_h1_level[[v]] <- L_t1
    real_h4_level[[v]] <- L_t4
  } else if (tr_v == "diff") {
    pred_h1_level[[v]] <- base_h1 + yhat_h1
    pred_h4_level[[v]] <- base_h4 + yhat_h4
    real_h1_level[[v]] <- L_t1
    real_h4_level[[v]] <- L_t4
  } else {
    pred_h1_level[[v]] <- yhat_h1
    pred_h4_level[[v]] <- yhat_h4
    real_h1_level[[v]] <- L_t1
    real_h4_level[[v]] <- L_t4
  }
}

# for (v in vars_Y) {
#   tr_v <- get_transform(v)
#   # base for back-transform
#   base_h1 <- levels_df[[v]][oos_origin_idx]             # L_t (per h=1)
#   base_h4 <- levels_df[[v]][oos_origin_idx + (h_long-1)]# L_{t+3} (per h=4)
#   # reals (levels) a t+1 e t+4
#   L_t1 <- levels_df[[v]][oos_origin_idx + 1]
#   L_t4 <- levels_df[[v]][oos_origin_idx + h_long]
#   
#   # forecasts in transformed units
#   yhat_h1 <- pred_h1_real[, v]
#   yhat_h4 <- pred_h4_real[, v]
#   
#   # back-transform
#   if (tr_v == "logdiff") {
#     # y is 100 * Δlog -> multiplier exp(y/100)
#     pred_h1_level[[v]] <- base_h1 * exp(yhat_h1 / 100)
#     pred_h4_level[[v]] <- base_h4 * exp(yhat_h4 / 100)
#   } else if (tr_v == "diff") {
#     pred_h1_level[[v]] <- base_h1 + yhat_h1
#     pred_h4_level[[v]] <- base_h4 + yhat_h4
#   } else { # none
#     pred_h1_level[[v]] <- yhat_h1
#     pred_h4_level[[v]] <- yhat_h4
#   }
#   
#   # reals in levels
#   real_h1_level[[v]] <- L_t1
#   real_h4_level[[v]] <- L_t4
# }

# Lists -> data.frame
pred_h1_level <- as.data.frame(pred_h1_level, check.names = FALSE)
pred_h4_level <- as.data.frame(pred_h4_level, check.names = FALSE)
real_h1_level <- as.data.frame(real_h1_level, check.names = FALSE)
real_h4_level <- as.data.frame(real_h4_level, check.names = FALSE)

# ---- Metrics in levels ----
rmse_h1_lvl <- sapply(vars_Y, function(v) rmse(real_h1_level[[v]] - pred_h1_level[[v]]))
rmse_h4_lvl <- sapply(vars_Y, function(v) rmse(real_h4_level[[v]] - pred_h4_level[[v]]))
mae_h1_lvl  <- sapply(vars_Y, function(v) mae(real_h1_level[[v]] - pred_h1_level[[v]]))
mae_h4_lvl  <- sapply(vars_Y, function(v) mae(real_h4_level[[v]] - pred_h4_level[[v]]))

metrics_levels <- data.frame(
  variable = vars_Y,
  RMSE_h1  = as.numeric(rmse_h1_lvl),
  RMSE_h4  = as.numeric(rmse_h4_lvl),
  MAE_h1   = as.numeric(mae_h1_lvl),
  MAE_h4   = as.numeric(mae_h4_lvl),
  row.names = NULL
)
cat("\n=== Forecast Evaluation (Levels) ===\n")
print(metrics_levels)

# ---- Build oos_results (trasformed + levels) ----
oos_results <- data.frame(
  origin_date = origin_dates,
  h1_date = h1_dates,
  h4_date = h4_dates
)

# transformed 
oos_results <- cbind(
  oos_results,
  setNames(as.data.frame(pred_h1_real), paste0(vars_Y, "_pred_h1_tr")),
  setNames(as.data.frame(real_h1_real), paste0(vars_Y, "_real_h1_tr")),
  setNames(as.data.frame(pred_h4_real), paste0(vars_Y, "_pred_h4_tr")),
  setNames(as.data.frame(real_h4_real), paste0(vars_Y, "_real_h4_tr"))
)

# levels (back-transformed)
oos_results <- cbind(
  oos_results,
  setNames(pred_h1_level, paste0(vars_Y, "_pred_h1_lvl")),
  setNames(real_h1_level, paste0(vars_Y, "_real_h1_lvl")),
  setNames(pred_h4_level, paste0(vars_Y, "_pred_h4_lvl")),
  setNames(real_h4_level, paste0(vars_Y, "_real_h4_lvl"))
)

cat("\n=== OOS results (head) ===\n")
print(head(oos_results, 5))

# ================================
#  6. Plot — OOS transformed scale (model units)
# ================================
graphics.off()

plot_df_h1_tr <- oos_results %>%
  dplyr::select(h1_date, dplyr::ends_with("_pred_h1_tr"), dplyr::ends_with("_real_h1_tr")) %>%
  dplyr::rename(date = h1_date)

plot_df_h4_tr <- oos_results %>%
  dplyr::select(h4_date, dplyr::ends_with("_pred_h4_tr"), dplyr::ends_with("_real_h4_tr")) %>%
  dplyr::rename(date = h4_date)

long_h1_tr <- plot_df_h1_tr %>%
  pivot_longer(
    cols = -date,
    names_to = c("variable", "type"),
    names_pattern = "(.+)_(pred|real)_h1_tr",
    values_to = "value"
  ) %>% mutate(horizon = "h=1 (1Q ahead)")

long_h4_tr <- plot_df_h4_tr %>%
  pivot_longer(
    cols = -date,
    names_to = c("variable", "type"),
    names_pattern = "(.+)_(pred|real)_h4_tr",
    values_to = "value"
  ) %>% mutate(horizon = "h=4 (1Y ahead)")

long_all_tr <- bind_rows(long_h1_tr, long_h4_tr) %>%
  mutate(
    date = as.Date(date),
    type = factor(type, levels = c("real","pred")),
    horizon = factor(horizon, levels = c("h=1 (1Q ahead)", "h=4 (1Y ahead)"))
  )

long_all_tr$variable <- recode(long_all_tr$variable,
                               rvgdp    = "Real GDP Growth",
                               cpi  = "CPI",
                               wkfreuro = "Exchange Rate"
)

p_tr <- ggplot(long_all_tr, aes(x = date, y = value, color = type, linetype = type)) +
  geom_line(linewidth = 0.7) +
  facet_grid(variable ~ horizon, scales = "free_y") +
  scale_color_manual(values = c("real" = "black", "pred" = "steelblue")) +
  scale_linetype_manual(values = c("real" = "solid", "pred" = "dashed")) +
  labs(title = "OOS Forecasts — Transformed Units",
       x = "Date", y = "Transformed value", color = "", linetype = "") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top", strip.text = element_text(face="bold"))
print(p_tr)

# ================================
#  6.1 Plot — OOS levels (back-transformed)
# ================================
plot_df_h1_lvl <- oos_results %>%
  dplyr::select(h1_date, dplyr::ends_with("_pred_h1_lvl"), dplyr::ends_with("_real_h1_lvl")) %>%
  dplyr::rename(date = h1_date)

plot_df_h4_lvl <- oos_results %>%
  dplyr::select(h4_date, dplyr::ends_with("_pred_h4_lvl"), dplyr::ends_with("_real_h4_lvl")) %>%
  dplyr::rename(date = h4_date)

long_h1_lvl <- plot_df_h1_lvl %>%
  pivot_longer(
    cols = -date,
    names_to = c("variable", "type"),
    names_pattern = "(.+)_(pred|real)_h1_lvl",
    values_to = "value"
  ) %>% mutate(horizon = "h=1 (1Q ahead)")

long_h4_lvl <- plot_df_h4_lvl %>%
  pivot_longer(
    cols = -date,
    names_to = c("variable", "type"),
    names_pattern = "(.+)_(pred|real)_h4_lvl",
    values_to = "value"
  ) %>% mutate(horizon = "h=4 (1Y ahead)")

long_all_lvl <- bind_rows(long_h1_lvl, long_h4_lvl) %>%
  mutate(
    date = as.Date(date),
    type = factor(type, levels = c("real","pred")),
    horizon = factor(horizon, levels = c("h=1 (1Q ahead)", "h=4 (1Y ahead)"))
  )

long_all_lvl$variable <- recode(long_all_lvl$variable,
                                rvgdp    = "Real GDP Growth",
                                cpi  = "CPI",
                                wkfreuro = "Exchange Rate"
)

graphics.off()
p_lvl <- ggplot(long_all_lvl, aes(x = date, y = value, color = type, linetype = type)) +
  geom_line(linewidth = 0.7) +
  facet_grid(variable ~ horizon, scales = "free_y") +
  scale_color_manual(values = c("real" = "black", "pred" = "steelblue")) +
  scale_linetype_manual(values = c("real" = "solid", "pred" = "dashed")) +
  labs(title = "OOS Forecasts — Back-transformed Levels",
       x = "Date", y = "Level", color = "", linetype = "") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top", strip.text = element_text(face="bold"))
print(p_lvl)

summary(long_all_lvl$date)
unique(diff(sort(unique(long_all_lvl$date))))


# =========================
#   7. "NOW" FORECASTS from latest data (h=1 and h=4)
# =========================

# Keep only observed data up to today 
# Use full dataset up to today for the NOW forecast (post-training)
df_train_clean <- df %>%
  mutate(date = as.Date(paste0(date, "-01"))) %>%
  mutate(across(-date, as.numeric)) %>%
  filter(date <= as.Date("2025-07-01"))

for (i in seq_len(nrow(metadata))) {
  key <- metadata$ts_key[i]
  t_type <- metadata$transform[i]
  if (!is.na(t_type) && key %in% names(df_train_clean)) {
    df_train_clean[[key]] <- transform_variable(df_train_clean[[key]], t_type)
  }
}

# Remove first row after differencing
df_train_clean <- df_train_clean[-1, ]

# Standardize using the same approach as before
df_train_std <- df_train_clean
df_train_std[-1] <- scale(df_train_clean[-1])
# Decide predictor set for NOW (frozen on initial 80%)
X_now_names <- vars_X

# Build working frame
needed_cols <- c("date", vars_Y, X_now_names)
needed_cols <- intersect(needed_cols, names(df_train_clean))
df_now <- df_train_clean[, needed_cols, drop = FALSE]
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

# Only impute standardised X values: NA -> 0 (which on a standard scale = mean of the variable)
X_std_imp <- X_std
for (j in seq_len(ncol(X_std_imp))) {
  X_std_imp[is.na(X_std_imp[, j]), j] <- 0
}

# Select number of factors r (reuse r_fixed if available; else choose on initial 80%)
if (exists("r_fixed")) {
  r_now <- r_fixed
  message("NOW: using r_fixed = ", r_now)
} else {
  t0_now <- max(1L, floor(0.8 * nrow(X_std)))
  X0_now <- X_std_imp[1:t0_now, , drop = FALSE]
  bn_now <- bn_select_r(X0_now, crit = "ICp3")
  r_now  <- max(1, bn_now$r)
  message("NOW: r selected by ICp3 on initial 80% = ", r_now)
}

# PCA on standardized X (full history)
pca_now <- prcomp(X_std_imp, center = FALSE, scale. = FALSE)

# Dimensional protection of r 
r_use <- min(r_now, ncol(pca_now$x))
F_all <- pca_now$x[, 1:r_use, drop = FALSE]
colnames(F_all) <- paste0("F", seq_len(r_use))

# OLS with factor lags (L=1): Y_std ~ [1, F_t, F_{t-1}]
L_ols <- 1
Tn <- nrow(F_all); r <- ncol(F_all)
stopifnot(Tn > L_ols)

F_curr <- F_all[(1+L_ols):Tn, , drop = FALSE]        # F_t
F_lag1 <- F_all[1:(Tn-L_ols), , drop = FALSE]        # F_{t-1}

# Design matrix with stable names [F1..Fr | L1_F1..L1_Fr]
X_ols <- cbind(F_curr, F_lag1)
colnames(X_ols) <- c(paste0("F", 1:r), paste0("L1_F", 1:r))
Xdf <- as.data.frame(X_ols)

Y_ols <- Y_std[(1+L_ols):Tn, , drop = FALSE]         # align Y

# --- Robust OLS to coefficients names---
X_ols_const_names <- c("(Intercept)", colnames(X_ols))
beta_hat_now <- array(
  NA_real_,
  dim = c(length(X_ols_const_names), ncol(Y_ols)),
  dimnames = list(X_ols_const_names, colnames(Y_ols))
)

for (j in seq_len(ncol(Y_ols))) {
  fit_j <- lm(Y_ols[, j] ~ ., data = Xdf)  # include intercept
  cf <- coef(fit_j)
  # reordnung/cleaning (if something missing, put 0)
  beta_vec <- setNames(numeric(length(X_ols_const_names)), X_ols_const_names)
  beta_vec[names(cf)] <- cf
  beta_hat_now[, j] <- beta_vec
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

# Map factor forecasts to Y (std) using same ordnung of coefficients
x_h1 <- c("(Intercept)" = 1,
          setNames(F_h1, paste0("F", 1:r)),
          setNames(F_t,  paste0("L1_F", 1:r)))
x_h4 <- c("(Intercept)" = 1,
          setNames(F_h4, paste0("F", 1:r)),
          setNames(F_h3, paste0("L1_F", 1:r)))

Y_h1_std <- as.numeric(t(x_h1[X_ols_const_names]) %*% beta_hat_now)
Y_h4_std <- as.numeric(t(x_h4[X_ols_const_names]) %*% beta_hat_now)

# De-standardize back to transformed units
Y_h1_real <- Y_h1_std * Y_sd + Y_mean
Y_h4_real <- Y_h4_std * Y_sd + Y_mean
names(Y_h1_real) <- vars_Y
names(Y_h4_real) <- vars_Y

# Target dates (next quarter and +4 quarters)
last_date <- max(df_train_clean$date, na.rm = TRUE)
h1_date   <- seq(last_date, by = "quarter", length.out = 2)[2]
h4_date   <- seq(last_date, by = "quarter", length.out = 5)[5]

now_table <- data.frame(
  variable = vars_Y,
  target_h1_date = h1_date,
  forecast_h1    = as.numeric(Y_h1_real),
  target_h4_date = h4_date,
  forecast_h4    = as.numeric(Y_h4_real)
)

# =========================
#   7.1. NOW: Route h=1..4 e back-transform to levels (recursive)
# =========================

# Factors predicted F_{t+1..t+4} by the VAR already estimated (fc_now)
fac_names <- colnames(F_all)
get_fc <- function(step) sapply(fac_names, function(nm) fc_now$fcst[[nm]][step, "fcst"])
F_h1 <- get_fc(1)
F_h2 <- get_fc(2)
F_h3 <- get_fc(3)
F_h4 <- get_fc(4)
F_t  <- as.numeric(F_all[nrow(F_all), , drop = FALSE])

# Same order as regressors of OLS: x_h = [ (Intercept), F_{t+h}, F_{t+h-1} ]
make_xh <- function(F_h, F_hm1, r) {
  c("(Intercept)" = 1,
    setNames(F_h,   paste0("F", 1:r)),
    setNames(F_hm1, paste0("L1_F", 1:r)))
}
r <- ncol(F_all)
x1 <- make_xh(F_h1, F_t,  r)
x2 <- make_xh(F_h2, F_h1, r)
x3 <- make_xh(F_h3, F_h2, r)
x4 <- make_xh(F_h4, F_h3, r)

# Predictions in TRANSFORMED UNITS (std -> real of the model)
yh1_std <- as.numeric(t(x1[X_ols_const_names]) %*% beta_hat_now)
yh2_std <- as.numeric(t(x2[X_ols_const_names]) %*% beta_hat_now)
yh3_std <- as.numeric(t(x3[X_ols_const_names]) %*% beta_hat_now)
yh4_std <- as.numeric(t(x4[X_ols_const_names]) %*% beta_hat_now)

yh1 <- yh1_std * Y_sd + Y_mean
yh2 <- yh2_std * Y_sd + Y_mean
yh3 <- yh3_std * Y_sd + Y_mean
yh4 <- yh4_std * Y_sd + Y_mean
names(yh1) <- names(yh2) <- names(yh3) <- names(yh4) <- vars_Y

# Recursive back-transform at LEVELS
get_transform <- function(var) {
  i <- match(var, metadata$ts_key)
  tr <- if (is.na(i)) NA_character_ else metadata$transform[i]
  ifelse(is.na(tr), "none", tolower(tr))
}

levels_df <- df %>%
  mutate(date = as.Date(paste0(date, "-01"))) %>%
  dplyr::select(date, all_of(vars_Y))
last_date <- max(df_train_clean$date, na.rm = TRUE)

# last observed level L_t for each Y
L_t_vec <- setNames(numeric(length(vars_Y)), vars_Y)
for (v in vars_Y) {
  L_t_val <- levels_df[[v]][levels_df$date == last_date]
  if (length(L_t_val) != 1 || !is.finite(L_t_val)) {
    stop("NOW back-transform: don't find L_t for variabile ", v)
  }
  L_t_vec[v] <- L_t_val
}

# Propagation L_{t+h} starting from L_t using yh1..yh4
L_h1 <- L_h2 <- L_h3 <- L_h4 <- setNames(numeric(length(vars_Y)), vars_Y)
for (v in vars_Y) {
  tr_v <- get_transform(v)
  Lt   <- L_t_vec[v]
  if (v %in% growth_vars) {
    # GDP stays in % growth
    L1 <- yh1[v]; L2 <- yh2[v]; L3 <- yh3[v]; L4 <- yh4[v]
  } else if (tr_v == "logdiff") {
    L1 <- Lt * exp(yh1[v] / 100)
    L2 <- L1 * exp(yh2[v] / 100)
    L3 <- L2 * exp(yh3[v] / 100)
    L4 <- L3 * exp(yh4[v] / 100)
  } else if (tr_v == "diff") {
    L1 <- Lt + yh1[v]
    L2 <- L1 + yh2[v]
    L3 <- L2 + yh3[v]
    L4 <- L3 + yh4[v]
  } else {
    L1 <- yh1[v]; L2 <- yh2[v]; L3 <- yh3[v]; L4 <- yh4[v]
  }
  L_h1[v] <- L1; L_h2[v] <- L2; L_h3[v] <- L3; L_h4[v] <- L4
}

# for (v in vars_Y) {
#   tr_v <- get_transform(v)
#   Lt   <- L_t_vec[v]
#   if (tr_v == "logdiff") {
#     L1 <- Lt * exp(yh1[v] / 100)
#     L2 <- L1 * exp(yh2[v] / 100)
#     L3 <- L2 * exp(yh3[v] / 100)
#     L4 <- L3 * exp(yh4[v] / 100)
#   } else if (tr_v == "diff") {
#     L1 <- Lt + yh1[v]
#     L2 <- L1 + yh2[v]
#     L3 <- L2 + yh3[v]
#     L4 <- L3 + yh4[v]
#   } else { # none
#     L1 <- yh1[v]; L2 <- yh2[v]; L3 <- yh3[v]; L4 <- yh4[v]
#   }
#   L_h1[v] <- L1; L_h2[v] <- L2; L_h3[v] <- L3; L_h4[v] <- L4
# }

# FINAL TABLES
h1_date <- seq(last_date, by = "quarter", length.out = 2)[2]
h4_date <- seq(last_date, by = "quarter", length.out = 5)[5]

now_table_transformed <- data.frame(
  variable       = vars_Y,
  target_h1_date = h1_date,
  forecast_h1    = as.numeric(yh1[vars_Y]),
  target_h4_date = h4_date,
  forecast_h4    = as.numeric(yh4[vars_Y]),
  row.names = NULL
)

now_table_levels <- data.frame(
  variable            = vars_Y,
  target_h1_date      = h1_date,
  forecast_h1_level   = as.numeric(L_h1[vars_Y]),
  target_h4_date      = h4_date,
  forecast_h4_level   = as.numeric(L_h4[vars_Y]),
  row.names = NULL, check.names = FALSE
)

print(now_table_transformed)  # Transformed (Δlog*100 or diff)
print(now_table_levels) # Levels

# =========================
#   8. NOW: PLOT VS KOF forecasts
# =========================

last_date <- max(df_train_clean$date, na.rm = TRUE)
h1_date   <- now_table_levels$target_h1_date[1]
h4_date   <- now_table_levels$target_h4_date[1]

levels_recent <- df_train %>%
  dplyr::select(date, all_of(vars_Y)) %>%
  dplyr::filter(date >= as.Date("2022-01-01")) %>%
  tidyr::pivot_longer(-date, names_to = "variable", values_to = "value")

our_pts <- bind_rows(
  now_table_levels %>% transmute(variable, date = target_h1_date, value = forecast_h1_level, horizon = "h1", source = "Our"),
  now_table_levels %>% transmute(variable, date = target_h4_date, value = forecast_h4_level, horizon = "h4", source = "Our")
)

df_growth <- df %>%
  mutate(date = as.Date(paste0(date, "-01"))) %>%
  arrange(date) %>%
  mutate(rvgdp_growth = c(NA, diff(log(rvgdp)) * 100))  # QoQ %

# Replace GDP level with growth in KOF comparison
df_kof_growth <- df_growth %>%
  dplyr::select(date, rvgdp_growth, cpi, wkfreuro) %>%
  rename(rvgdp = rvgdp_growth)

df_kof_levels <- df_kof_growth

get_kof_value <- function(date_target, varname) {
  if (!varname %in% names(df_kof_levels)) return(NA_real_)
  out <- df_kof_levels %>% dplyr::filter(date == date_target) %>% dplyr::pull(!!sym(varname))
  if (length(out) == 0) NA_real_ else out[1]
}

kof_pts <- bind_rows(
  tibble(variable = vars_Y, date = h1_date,
         value = sapply(vars_Y, function(v) get_kof_value(h1_date, v)),
         horizon = "h1", source = "KOF"),
  tibble(variable = vars_Y, date = h4_date,
         value = sapply(vars_Y, function(v) get_kof_value(h4_date, v)),
         horizon = "h4", source = "KOF")
) %>% dplyr::filter(is.finite(value))

pal_colors <- c("Our" = "#1f77b4", "KOF" = "#ff7f0e")
pal_shapes <- c("Our" = 16,       "KOF" = 17)

# ensure date are Date
levels_recent$date <- as.Date(levels_recent$date)
our_pts$date       <- as.Date(our_pts$date)
kof_pts$date       <- as.Date(kof_pts$date)

# Chart names
facet_labels <- c(
  "rvgdp"    = "Real GDP Growth",
  "cpi"      = "Inflation (CPI)",
  "wkfreuro" = "CHF/EUR Exchange Rate"
)

p <- ggplot() +
  geom_line(data = levels_recent,
            aes(x = date, y = value),
            linewidth = 0.6, color = "grey40") +
  geom_vline(xintercept = as.numeric(last_date),
             linetype = "solid", color = "grey50") +
  geom_point(data = bind_rows(our_pts, kof_pts),
             aes(x = date, y = value, color = source, shape = source),
             size = 3.2) +
  facet_wrap(
    ~ variable,
    scales = "free_y",
    ncol   = 1,
    labeller = as_labeller(facet_labels)   # <-- NUOVO: titoli pannelli
  ) +
# legenda 
  scale_color_manual(
    values = pal_colors,
    breaks = c("Our", "KOF"),
    labels = c("Dynamic Factor Model", "KOF forecast"),
    name   = ""
  ) +
  scale_shape_manual(
    values = pal_shapes,
    breaks = c("Our", "KOF"),
    labels = c("Dynamic Factor Model", "KOF forecast"),
    name   = ""
  ) +
  labs(title = "Forecast Comparison: Dynamic Factor Model (PCA+VAR) vs KOF",
       subtitle = paste0("Historical until last observed value: ",
                         format(last_date, "%Y-%m-%d")),
       x = "Date", y = "Levels") +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold")
  )

print(p)


# ============================================
#   9. STATE SPACE: Dynamic Factor Model Initialization (PCA/SVD → F0, Lambda0, R0, A0, Q0)
# Compatible with KFAS: x_t = Λ f_t + e_t ; f_t = A f_{t-1} + u_t
# ============================================
# 
# 
# # ============================================
# # STATE SPACE: Dynamic Factor Model Initialization (PCA/SVD → F0, Lambda0, R0, A0, Q0)
# # Compatible with KFAS: x_t = Λ f_t + e_t ; f_t = A f_{t-1} + u_t
# # ============================================
# #In our workflow, we currently have roughly these main parts:
# 
# # Data loading + cleaning
# 
# # Transformations (logdiff, diff, etc.)
# 
# # Standardization (df_train_std, vars_X, vars_Y)
# 
# # Number of factors selection (Bai–Ng, ICp3) → gives r_fixed
# 
# # now we switch to STATE SPACE MODEL
# 
# 
# 
# dfm_init <- function(df_train_std, vars_X, r_fixed, ar_cap = 0.98) {
#   stopifnot("date" %in% names(df_train_std))
#   X_std <- as.matrix(df_train_std[, vars_X, drop = FALSE])
#   Tn <- nrow(X_std); Nn <- ncol(X_std)
#   if (Tn < 5 || Nn < 1) stop("X_std has insufficient dimensions.")
#   
#   # --- 1) Prepare matrix for SVD: replace NA with 0 (mean=0 after standardization) ---
#   X_svd <- X_std
#   X_svd[!is.finite(X_svd)] <- 0
#   
#   # --- 2) SVD on standardized X (already zero mean) ---
#   #     Decomposition: X ≈ U D V'
#   #     Factors and loadings (Stock–Watson convention):
#   #     F0 = sqrt(T) * U_r         ;  Lambda0 = V_r * (D_r / sqrt(T))
#   sv <- svd(X_svd)  # X_svd = U D V'
#   r <- max(1, min(r_fixed, length(sv$d), Nn, floor(0.5 * min(Tn, Nn))))
#   
#   U_r <- sv$u[, 1:r, drop = FALSE]
#   D_r <- sv$d[1:r]
#   V_r <- sv$v[, 1:r, drop = FALSE]
#   
#   F0       <- sqrt(Tn) * U_r                            # T x r
#   Lambda0  <- V_r %*% diag(D_r / sqrt(Tn), nrow = r)    # N x r
#   
#   # --- 3) Idiosyncratic variances (R0, diagonal) from residuals X - F0 %*% t(Lambda0) ---
#   E0 <- X_svd - F0 %*% t(Lambda0)
#   R_diag <- colMeans(E0^2)  # N x 1
#   R_diag[!is.finite(R_diag)] <- 1e-4
#   R_diag <- pmax(R_diag, 1e-4)
#   R0 <- diag(R_diag, nrow = Nn)
#   
#   # --- 4) Factor dynamics: independent AR(1) for each factor ---
#   phi <- rep(NA_real_, r)
#   sig_u <- rep(NA_real_, r)
#   for (j in seq_len(r)) {
#     y  <- F0[-1, j]
#     x  <- F0[-Tn, j]
#     # OLS without intercept: f_t = phi * f_{t-1} + u_t
#     phi_j <- sum(x * y) / sum(x * x)
#     if (!is.finite(phi_j)) phi_j <- 0
#     phi_j <- max(-ar_cap, min(ar_cap, phi_j))  # stability cap
#     u_j   <- y - phi_j * x
#     s2_j  <- mean(u_j^2)
#     if (!is.finite(s2_j) || s2_j <= 0) s2_j <- 1e-4
#     phi[j]  <- phi_j
#     sig_u[j] <- s2_j
#   }
#   A0 <- diag(phi, nrow = r)                 # r x r transition matrix
#   Q0 <- diag(pmax(sig_u, 1e-4), nrow = r)   # r x r state noise covariance
#   
#   # --- 5) Build the state-space model for KFAS ---
#   # Measurement eq: x_t = Λ f_t + ε_t,  ε_t ~ N(0, R0)
#   # Transition eq:  f_t = A f_{t-1} + u_t,  u_t ~ N(0, Q0)
#   Zt <- array(0, dim = c(Nn, r, 1)); Zt[, , 1] <- Lambda0
#   Tt <- array(A0, dim = c(r, r, 1))
#   Rt <- diag(r)
#   
#   model <- SSModel(
#     X_std ~ -1 + SSMcustom(
#       Z = Zt, T = Tt, R = Rt,
#       Q = Q0, a1 = rep(0, r), P1 = diag(10, r),
#       state_names = paste0("F", 1:r)
#     ),
#     H = R0
#   )
#   
#   list(
#     X_std   = X_std,          # standardized data used in measurement equation
#     r       = r,
#     F0      = F0,             # initial factors (T x r)
#     Lambda0 = Lambda0,        # initial loadings (N x r)
#     R0      = R0,             # idiosyncratic variance matrix
#     A0      = A0,             # AR(1) transition matrix
#     Q0      = Q0,             # factor innovation variance
#     model   = model           # SSModel ready for EM / Kalman
#   )
# }
# 
# # ====== USAGE EXAMPLE ======
# # Assuming df_train_std, vars_X, r_fixed already exist:
# init <- dfm_init(df_train_std, vars_X, r_fixed = r_fixed)
# 
# # Optionally: first smoothing pass with initial parameters
# ks0 <- KFS(init$model, smoothing = c("state", "mean"))
# F_smooth0 <- ks0$alphahat            # T x r smoothed factors
# X_fitted0 <- signal(init$model)$signal  # fitted measurements
# 
# # Optional: run EM to optimize H and Q (and optionally Z, T if parameterized)
# # fit <- fitSSM(init$model,
# #               inits = c(log(diag(init$model$H)), log(diag(init$model$Q))),
# #               method = "BFGS")
# # mod_hat <- fit$model
# # ks_hat  <- KFS(mod_hat, smoothing = c("state", "mean"))
# 
# # ============================================================
# # ML estimation (quasi-ML) of H (idiosyncratic), Q (factor shocks), A (AR(1))
# # Λ (loadings) kept fixed to PCA for identification.
# # ============================================================
# 
# library(KFAS)
# 
# ssm0   <- init$model
# Nn     <- nrow(ssm0$y)             # number of series (N)
# r      <- dim(ssm0$Z)[2]           # number of factors (r)
# phi_cap <- 0.98                     # stability cap for AR(1)
# 
# # --- Parameterization & constraints ---
# # par = [ log(diag(H))_(N) , log(diag(Q))_(r) , atanh(phi/phi_cap)_(r) ]
# par0 <- c(
#   log(diag(ssm0$H)),                        # H diagonal (positive)
#   log(diag(ssm0$Q)),                        # Q diagonal (positive)
#   atanh(diag(ssm0$T[,,1]) / phi_cap)        # AR(1) φ_j ∈ (-phi_cap, phi_cap)
# )
# 
# # Map parameter vector -> model (H, Q, A updated; Λ fixed)
# updatefn <- function(par, model) {
#   idxH <- 1:Nn
#   idxQ <- (Nn + 1):(Nn + r)
#   idxA <- (Nn + r + 1):(Nn + 2*r)
#   
#   Hdiag <- exp(par[idxH])                    # ensure > 0
#   Qdiag <- exp(par[idxQ])                    # ensure > 0
#   phi   <- phi_cap * tanh(par[idxA])         # ensure |phi| < phi_cap
#   
#   model$H[,]     <- diag(Hdiag, nrow = Nn)
#   model$Q[,]     <- diag(Qdiag, nrow = r)
#   model$T[,,1]   <- diag(phi, nrow = r)      # AR(1) diagonal transition
#   # Z (loadings) and R (identity) are kept as-in (from PCA init)
#   model
# }
# 
# # --- Fit by BFGS (quasi-ML) ---
# fit <- fitSSM(
#   inits    = par0,
#   model    = ssm0,
#   updatefn = updatefn,
#   method   = "BFGS",
#   control  = list(maxit = 300, reltol = 1e-8)
# )
# 
# mod_hat <- fit$model
# 
# # --- Smoothing with estimated parameters ---
# ks_hat <- KFS(mod_hat, smoothing = c("state", "mean"))
# F_smooth <- ks_hat$alphahat               # T x r smoothed factors
# X_fitted <- signal(mod_hat)$signal        # T x N fitted measurements (in standardized units)
# 
# # --- Basic diagnostics ---
# ll    <- as.numeric(logLik(mod_hat))
# k_par <- length(par0)
# n_obs <- sum(is.finite(mod_hat$y))        # total observed measurements
# AIC   <- -2*ll + 2*k_par
# BIC   <- -2*ll + k_par*log(n_obs)
# 
# cat(sprintf("\n=== Quasi-ML fit ===\nlogLik = %.2f | AIC = %.2f | BIC = %.2f | k = %d | n = %d\n",
#             ll, AIC, BIC, k_par, n_obs))
# 
# # Standardized recursive residuals (one-step-ahead)
# res_std <- residuals(mod_hat, type = "recursive", standardized = TRUE)
# # You can inspect series-wise residual distribution / autocorrelation as needed.
# 
# # --- Extract estimated components (useful for reporting) ---
# H_hat   <- mod_hat$H[,]                   # diagonal idiosyncratic variances (N x N)
# Q_hat   <- mod_hat$Q[,]                   # diagonal factor shock variances (r x r)
# A_hat   <- mod_hat$T[,,1]                 # AR(1) coefficients (r x r, diagonal)
# Lambda  <- mod_hat$Z[,,1]                 # loadings (N x r) — fixed from PCA init
# 
# # (Optional) Reconstruct idiosyncratic residuals in standardized units:
# # E_hat = y - Λ * f |t (use smoothed states for f)
# Y_std_hat <- t(Lambda %*% t(F_smooth))    # T x N
# E_hat <- mod_hat$y - t(Y_std_hat)         # N x T in KFAS, so transpose back if needed
# # NOTE: KFAS stores y as N x T; signal() returns T x N; be mindful of dims if you use these.
# 
# # ============================================================
# # Forecasting h steps ahead (in transformed/standardized units)
# # ============================================================
# 
# h_steps <- 4L
# 
# # Extend the model with h rows of NA to obtain prediction for the measurements
# mod_fc <- mod_hat
# Y_obs  <- mod_fc$y                         # N x T
# Y_aug  <- cbind(Y_obs, matrix(NA_real_, nrow = Nn, ncol = h_steps))
# mod_fc$y <- Y_aug
# 
# # Filter forward to get one-step-ahead predictions
# kf_fc <- KFS(mod_fc, filtering = "state", smoothing = "none")
# 
# # Predicted signals (measurement means) for the whole augmented sample:
# sig_aug <- signal(mod_fc)$signal           # (T + h) x N  in standardized units
# 
# # Grab only the last h rows (forecasts)
# pred_transformed <- sig_aug[(nrow(sig_aug) - h_steps + 1):nrow(sig_aug), , drop = FALSE]
# colnames(pred_transformed) <- colnames(init$X_std)  # same measurement names (vars_X order)
# 
# # ============================================================
# # Map forecasts to your target variables and back-transform to levels (optional)
# # ============================================================
# # pred_transformed is in the *measurement space* (standardized transformed units).
# # If your targets vars_Y are among the measurement variables, subset and de-standardize using
# # the moments you use in your pipeline (mean/sd from the training window you prefer).
# # Otherwise, build a Y from the factors via a separate mapping (like your OLS step).
# 
# # Example: if vars_Y ⊆ colnames(pred_transformed) and you want to de-standardize with
# # full-history moments (Y_mean, Y_sd) as you already do for nowcasts:
# 
# pred_Y_tr <- pred_transformed[, vars_Y, drop = FALSE]  # transformed (Δlog*100 or diff), std units
# # De-standardize back to transformed units:
# # (Assumes you have Y_mean and Y_sd consistent with how you standardized Y in your pipeline.)
# pred_Y_transformed <- sweep(pred_Y_tr, 2, Y_sd, `*`)
# pred_Y_transformed <- sweep(pred_Y_transformed, 2, Y_mean, `+`)
# 
# # If you want back to *levels*, reuse your back-transform logic per variable:
# #   - logdiff:  Level_{t+h} = Base_{t+h-1} * exp( y_{t+h} / 100 )
# #   - diff:     Level_{t+h} = Base_{t+h-1} + y_{t+h}
# # where Base depends on the horizon (for h=1 use L_t, for h=4 use L_{t+3}).
# # You already implemented this in your PCA+VAR part; plug the same helper here.
# 
# # (Optional) Assemble a tidy forecast table
# fc_dates <- seq(from = max(df_train_std$date), by = "quarter", length.out = h_steps)
# now_table_transformed_ssm <- data.frame(
#   date = fc_dates,
#   pred_Y_transformed,
#   row.names = NULL,
#   check.names = FALSE
# )
# 
# print(head(now_table_transformed_ssm, 4))

# ============================================================================
# Enhanced State-Space Dynamic Factor Model Implementation
# Following analyst's detailed specification for KFAS implementation
# ============================================================================

cat("\n=== Enhanced State-Space Dynamic Factor Model (KFAS) ===\n")

# ============================================================================
# Data Preparation for State-Space Model
# Use all series (targets + predictors) for factor extraction
# ============================================================================

# Create the full measurement matrix Y_mat (T x m) including both targets and predictors
Y_cols <- c(vars_Y, vars_X)  # targets first, then predictors
Y_mat <- as.matrix(df_train_std[, Y_cols, drop = FALSE])
m <- ncol(Y_mat)  # number of series
r <- r_fixed      # number of factors from Bai-Ng selection

cat(sprintf("Building state-space model with m=%d series and r=%d factors\n", m, r))

# ============================================================================
# State-Space Model Construction with KFAS
# Observation equation: y_t = Λ f_t + ε_t
# State equation: f_t = Φ f_{t-1} + u_t
# ============================================================================

# Initialize system matrices with appropriate dimensions
Zt <- matrix(NA, m, r)          # loading matrix (m x r) with NAs for estimation
Tt <- matrix(NA, nrow = r, ncol = r) # AR(2) factor transition (r x r) with NAs on diagonal
Rt <- diag(1, r)                # identity matrix (r x r)
Qt <- diag(1, r)                # identity matrix (r x r) - fixed for identification
Ht <- diag(NA, m)               # diagonal observation covariance (m x m) with NAs

# Initial state setup for diffuse initialization
a1 <- rep(0, r)                 # initial state means
P1 <- matrix(0, r, r)           # diffuse part (will use P1inf for actual diffuse)
P1inf <- diag(1, r)             # diffuse initial covariance flags

# Build the SSModel object
model_ss <- SSModel(
  Y_mat ~ -1 + SSMcustom(
    Z = Zt, 
    T = Tt, 
    R = Rt, 
    Q = Qt,
    a1 = a1, 
    P1 = P1, 
    P1inf = P1inf,
    state_names = paste0("Factor", 1:r),
    index = 1:m
  ),
  H = Ht
)

cat("State-space model structure created\n")

# ============================================================================
# Parameter Estimation Setup
# Parameters to estimate: loadings (Λ), observation variances (H), AR coefficients (φ)
# ============================================================================

# Define parameter vector lengths
nZ   <- m * r            # number of loadings
nH   <- m                # number of unique H variances (diagonal)
nPhi <- r                # number of AR coefficients

# Update function to map parameter vector to model matrices
updatefn_ss <- function(par, model) {
  stopifnot(length(par) == nZ + nH + nPhi)
  
  # 1. Factor loadings (Z matrix) - fill by row-major order
  model$Z[,,1] <- matrix(par[1:nZ], nrow = m, ncol = r)
  
  # 2. Observation variances (H matrix diagonal) - exponential transformation
  h_vals <- exp(par[(nZ+1):(nZ+nH)])
  model$H[,,1] <- diag(h_vals, m)
  
  # 3. AR(1) coefficients (T matrix diagonal) - tanh transformation for stationarity
  phi_pars <- par[(nZ+nH+1):(nZ+nH+nPhi)]
  phi_vals <- tanh(phi_pars)  # maps to (-1, 1)
  model$T[,,1] <- diag(phi_vals, r)
  
  return(model)
}

# ============================================================================
# Initialize parameters using PCA results
# ============================================================================

cat("Initializing parameters using PCA...\n")

# Use PCA for initial loading estimates
Y_pca <- Y_mat
Y_pca[!is.finite(Y_pca)] <- 0  # replace NAs with 0 for PCA

pca_result <- prcomp(Y_pca, center = FALSE, scale. = FALSE)  # already standardized
init_loadings_pca <- pca_result$rotation[, 1:r, drop = FALSE]  # m x r

# Scale loadings to match our Q=I convention
eigenvals <- pca_result$sdev[1:r]^2
init_loadings_scaled <- init_loadings_pca %*% diag(sqrt(eigenvals), r)

# Initial parameter vector
init_loadings <- as.vector(init_loadings_scaled)  # vectorize by row-major order
init_H <- rep(log(0.5), nH)                      # log-variance initial guess
init_phi <- rep(0, nPhi)                         # start with no persistence

init_params <- c(init_loadings, init_H, init_phi)

cat(sprintf("Parameter vector length: %d (loadings: %d, H: %d, phi: %d)\n", 
            length(init_params), nZ, nH, nPhi))

# ============================================================================
# Maximum Likelihood Estimation
# ============================================================================

cat("Starting maximum likelihood estimation...\n")

# Fit the model using BFGS optimization
fit_ss <- fitSSM(
  model = model_ss, 
  inits = init_params, 
  updatefn = updatefn_ss, 
  method = "BFGS",
  control = list(maxit = 500, reltol = 1e-6, trace = 1)
)

# Check convergence
if (fit_ss$optim.out$convergence != 0) {
  cat("Warning: Optimization did not converge properly. Code:", fit_ss$optim.out$convergence, "\n")
} else {
  cat("Optimization converged successfully\n")
}

model_fit_ss <- fit_ss$model
logLik_val_ss <- logLik(model_fit_ss)

cat(sprintf("Final log-likelihood: %.2f\n", as.numeric(logLik_val_ss)))

# Extract estimated parameters
Lambda_hat <- model_fit_ss$Z[,,1]     # estimated loadings (m x r)
H_hat_ss <- model_fit_ss$H[,,1]       # estimated observation variances (m x m)
phi_hat <- diag(model_fit_ss$T[,,1])  # estimated AR coefficients (r x 1)

cat("Estimated AR(1) coefficients for factors:\n")
print(round(phi_hat, 3))

# ============================================================================
# Kalman Filtering and Smoothing
# ============================================================================

cat("Running Kalman filter and smoother...\n")

# Run Kalman filter and smoother
kfs_result <- KFS(model_fit_ss, smoothing = c("state", "mean"))

# Extract smoothed factors and fitted values
factors_smooth <- kfs_result$alphahat     # smoothed factors (T x r)
Y_fitted_ss <- fitted(model_fit_ss)       # fitted observations (T x m)

cat(sprintf("Extracted %d factors over %d time periods\n", ncol(factors_smooth), nrow(factors_smooth)))

# ============================================================================
# Forecasting with State-Space Model
# Generate nowcast (h=1) and 4-quarters ahead (h=4) forecasts
# ============================================================================

cat("Generating forecasts...\n")

# For KFAS, we need to extend the data matrix with NAs and use KFS for forecasting
h_forecast <- 4
T_obs <- nrow(Y_mat)

# Create extended data matrix with NAs for forecast periods
Y_extended <- rbind(Y_mat, matrix(NA, nrow = h_forecast, ncol = m))

# Create extended model
model_extended <- SSModel(
  Y_extended ~ -1 + SSMcustom(
    Z = model_fit_ss$Z, 
    T = model_fit_ss$T, 
    R = model_fit_ss$R, 
    Q = model_fit_ss$Q,
    a1 = model_fit_ss$a1, 
    P1 = model_fit_ss$P1, 
    P1inf = model_fit_ss$P1inf,
    state_names = paste0("Factor", 1:r),
    index = 1:m
  ),
  H = model_fit_ss$H
)

# Run Kalman filter on extended model to get forecasts
kfs_forecast <- KFS(model_extended, filtering = "state", smoothing = "none")

# Extract forecasts from the filtered estimates
# The forecasts are the fitted values for the NA periods
Y_forecast <- fitted(model_extended)[(T_obs + 1):(T_obs + h_forecast), , drop = FALSE]

# Extract forecasts for target variables
target_idx_ss <- match(vars_Y, colnames(Y_mat))

# Extract specific horizons
fc_h1_ss <- Y_forecast[1, target_idx_ss]  # nowcast (h=1)
fc_h4_ss <- Y_forecast[4, target_idx_ss]  # 4-quarters ahead (h=4)

# Create forecast results table (transformed units)
nowcast_results_ss <- data.frame(
  variable = vars_Y,
  last_obs_date = max(df_train_clean$date),
  nowcast_h1_date = seq(max(df_train_clean$date), by = "quarter", length.out = 2)[2],
  forecast_h1_tr = as.numeric(fc_h1_ss),
  forecast_h4_date = seq(max(df_train_clean$date), by = "quarter", length.out = 5)[5],
  forecast_h4_tr = as.numeric(fc_h4_ss)
)

cat("\n=== State-Space Model Forecasts (Transformed Units) ===\n")
print(nowcast_results_ss)

# ============================================================================
# Back-transformation to Levels (Optional)
# Convert forecasts from transformed units back to original levels
# ============================================================================

cat("\nBack-transforming forecasts to levels...\n")

# Get last observed values in original levels
last_levels_ss <- df_train[df_train$date == max(df_train_clean$date), vars_Y]
L_t_ss <- as.numeric(last_levels_ss)
names(L_t_ss) <- vars_Y

# Get transformation types from metadata
transform_type_ss <- setNames(metadata$transform[match(vars_Y, metadata$ts_key)], vars_Y)

# Back-transform forecasts
L_h1_ss <- numeric(length(vars_Y))
L_h4_ss <- numeric(length(vars_Y))

for (j in seq_along(vars_Y)) {
  v <- vars_Y[j]
  tr <- tolower(transform_type_ss[j])
  fh1 <- fc_h1_ss[j]  # h=1 forecast in transformed scale
  fh4 <- fc_h4_ss[j]  # h=4 forecast in transformed scale
  
  if (tr == "logdiff") {
    # For log-difference transformation: Level = last_level * exp(forecast/100)
    L_h1_ss[j] <- L_t_ss[v] * exp(fh1 / 100)
    L_h4_ss[j] <- L_t_ss[v] * exp(fh4 / 100)
  } else if (tr == "diff") {
    # For simple difference: Level = last_level + forecast
    L_h1_ss[j] <- L_t_ss[v] + fh1
    L_h4_ss[j] <- L_t_ss[v] + fh4
  } else {
    # No transformation case
    L_h1_ss[j] <- fh1
    L_h4_ss[j] <- fh4
  }
}

# Create back-transformed results table
backtransformed_results_ss <- data.frame(
  variable = vars_Y,
  nowcast_h1_date = nowcast_results_ss$nowcast_h1_date,
  forecast_h1_level = L_h1_ss,
  forecast_h4_date = nowcast_results_ss$forecast_h4_date,
  forecast_h4_level = L_h4_ss
)

cat("\n=== State-Space Model Forecasts (Level Units) ===\n")
print(backtransformed_results_ss)

# ============================================================================
# Model Diagnostics and Summary
# ============================================================================

cat("\n=== State-Space Model Diagnostics ===\n")

# Calculate information criteria
k_params_ss <- nZ + nH + nPhi
n_obs_ss <- sum(is.finite(model_fit_ss$y))
AIC_ss <- -2 * as.numeric(logLik_val_ss) + 2 * k_params_ss
BIC_ss <- -2 * as.numeric(logLik_val_ss) + k_params_ss * log(n_obs_ss)

cat(sprintf("Log-likelihood: %.2f\n", as.numeric(logLik_val_ss)))
cat(sprintf("AIC: %.2f\n", AIC_ss))
cat(sprintf("BIC: %.2f\n", BIC_ss))
cat(sprintf("Number of parameters: %d\n", k_params_ss))
cat(sprintf("Number of observations: %d\n", n_obs_ss))

# Residual diagnostics
# residuals_ss <- residuals(model_fit_ss, type = "recursive", standardized = TRUE)
# residual_mean <- apply(residuals_ss, 2, function(x) mean(x, na.rm = TRUE))
# residual_sd <- apply(residuals_ss, 2, function(x) sd(x, na.rm = TRUE))
# 
# cat(sprintf("\nResidual diagnostics - Mean absolute residual: %.4f\n", 
#             mean(abs(residual_mean), na.rm = TRUE)))
# cat(sprintf("Residual diagnostics - Average std deviation: %.4f\n", 
#             mean(residual_sd, na.rm = TRUE)))
# 
# # Factor loadings summary for target variables
# cat("\n=== Factor Loadings for Target Variables ===\n")
# loadings_targets <- Lambda_hat[target_idx_ss, , drop = FALSE]
# rownames(loadings_targets) <- vars_Y
# colnames(loadings_targets) <- paste0("Factor", 1:r)
# print(round(loadings_targets, 3))

# ============================================================================
# Comparison Summary
# Create summary for comparison with other methods
# ============================================================================

# Store state-space results for potential plotting/comparison
ssm_forecast_summary <- data.frame(
  variable = rep(vars_Y, 2),
  date = c(rep(nowcast_results_ss$nowcast_h1_date[1], length(vars_Y)),
           rep(nowcast_results_ss$forecast_h4_date[1], length(vars_Y))),
  value = c(L_h1_ss, L_h4_ss),
  horizon = rep(c("h1", "h4"), each = length(vars_Y)),
  source = "StateSpaceDFM"
)

cat("\n=== State-Space Dynamic Factor Model Implementation Complete ===\n")
cat("Key outputs:\n")
cat("- nowcast_results_ss: Forecasts in transformed units\n")
cat("- backtransformed_results_ss: Forecasts in level units\n")
cat("- factors_smooth: Extracted common factors\n")
cat("- Lambda_hat: Estimated factor loadings\n")
cat("- ssm_forecast_summary: Results formatted for comparison\n")

# ============================================================================
# State-Space Model vs KOF Forecast Comparison
# ============================================================================

cat("\n=== Comparing State-Space Model with KOF Forecasts ===\n")

# Get forecast dates from the state space model
h1_date_ss <- backtransformed_results_ss$nowcast_h1_date[1]
h4_date_ss <- backtransformed_results_ss$forecast_h4_date[1]

# Create state space forecast points for comparison
ss_pts <- bind_rows(
  tibble(variable = vars_Y, date = h1_date_ss,
         value = backtransformed_results_ss$forecast_h1_level,
         horizon = "h1", source = "StateSpace"),
  tibble(variable = vars_Y, date = h4_date_ss,
         value = backtransformed_results_ss$forecast_h4_level,
         horizon = "h4", source = "StateSpace")
)

# Get corresponding KOF forecasts for same dates
kof_pts_ss <- bind_rows(
  tibble(variable = vars_Y, date = h1_date_ss,
         value = sapply(vars_Y, function(v) get_kof_value(h1_date_ss, v)),
         horizon = "h1", source = "KOF"),
  tibble(variable = vars_Y, date = h4_date_ss,
         value = sapply(vars_Y, function(v) get_kof_value(h4_date_ss, v)),
         horizon = "h4", source = "KOF")
) %>% dplyr::filter(is.finite(value))

# Calculate differences (State Space - KOF)
differences_ss <- merge(
  ss_pts %>% dplyr::select(variable, date, horizon, value) %>% dplyr::rename(ss_value = value),
  kof_pts_ss %>% dplyr::select(variable, date, horizon, value) %>% dplyr::rename(kof_value = value),
  by = c("variable", "date", "horizon"),
  all.x = TRUE
) %>%
  dplyr::mutate(
    difference = ss_value - kof_value,
    abs_difference = abs(difference),
    relative_diff_pct = 100 * difference / abs(kof_value)
  )

cat("\n=== State-Space vs KOF Forecast Differences ===\n")
print(differences_ss)

# Summary statistics of differences
cat("\n=== Summary of Forecast Differences ===\n")
summary_stats_ss <- differences_ss %>%
  dplyr::group_by(horizon) %>%
  dplyr::summarise(
    mean_abs_diff = mean(abs_difference, na.rm = TRUE),
    median_abs_diff = median(abs_difference, na.rm = TRUE),
    max_abs_diff = max(abs_difference, na.rm = TRUE),
    mean_relative_diff = mean(abs(relative_diff_pct), na.rm = TRUE),
    .groups = 'drop'
  )
print(summary_stats_ss)

# Detailed comparison by variable
cat("\n=== Detailed Comparison by Variable ===\n")
detailed_comparison_ss <- differences_ss %>%
  dplyr::select(variable, horizon, ss_value, kof_value, difference, relative_diff_pct) %>%
  dplyr::arrange(variable, horizon)
print(detailed_comparison_ss)

# ============================================================================
# Visualization: State-Space Model vs KOF Forecasts
# ============================================================================

# Color palette for state space comparison
pal_colors_ss <- c("StateSpace" = "#2ca02c", "KOF" = "#ff7f0e")
pal_shapes_ss <- c("StateSpace" = 18, "KOF" = 17)

# Ensure dates are properly formatted
ss_pts$date <- as.Date(ss_pts$date)
kof_pts_ss$date <- as.Date(kof_pts_ss$date)

# Create comparison plot
p_ss_comparison <- ggplot() +
  geom_line(data = levels_recent,
            aes(x = date, y = value),
            linewidth = 0.6, color = "grey40") +
  geom_vline(xintercept = as.numeric(last_date),
             linetype = "solid", color = "grey50") +
  geom_point(data = bind_rows(ss_pts, kof_pts_ss),
             aes(x = date, y = value, color = source, shape = source),
             size = 3.2) +
  facet_wrap(
    ~ variable,
    scales = "free_y",
    ncol = 1,
    labeller = as_labeller(facet_labels)
  ) +
  scale_color_manual(
    values = pal_colors_ss,
    breaks = c("StateSpace", "KOF"),
    labels = c("State-Space DFM", "KOF forecast"),
    name = ""
  ) +
  scale_shape_manual(
    values = pal_shapes_ss,
    breaks = c("StateSpace", "KOF"),
    labels = c("State-Space DFM", "KOF forecast"),
    name = ""
  ) +
  labs(title = "Forecast Comparison: State-Space Dynamic Factor Model vs KOF",
       subtitle = paste0("Historical until last observed value: ",
                         format(last_date, "%Y-%m-%d")),
       x = "Date", y = "Levels") +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold")
  )

print(p_ss_comparison)

# ============================================================================
# Comparison of Both Models (PCA+VAR vs State-Space) vs KOF
# ============================================================================

cat("\n=== Three-Way Comparison: PCA+VAR vs State-Space vs KOF ===\n")

# Combine all forecasts for comparison
all_forecasts <- bind_rows(
  our_pts %>% mutate(source = "PCA_VAR"),  # Rename for clarity
  ss_pts,
  kof_pts_ss
)

# Calculate all pairwise differences
pca_vs_kof <- merge(
  our_pts %>% dplyr::select(variable, horizon, value) %>% dplyr::rename(pca_value = value),
  kof_pts_ss %>% dplyr::select(variable, horizon, value) %>% dplyr::rename(kof_value = value),
  by = c("variable", "horizon"), all.x = TRUE
) %>% dplyr::mutate(comparison = "PCA_VAR_vs_KOF", difference = pca_value - kof_value)

ss_vs_kof <- differences_ss %>% 
  dplyr::select(variable, horizon, difference) %>% 
  dplyr::mutate(comparison = "StateSpace_vs_KOF")

pca_vs_ss <- merge(
  our_pts %>% dplyr::select(variable, horizon, value) %>% dplyr::rename(pca_value = value),
  ss_pts %>% dplyr::select(variable, horizon, value) %>% dplyr::rename(ss_value = value),
  by = c("variable", "horizon"), all.x = TRUE
) %>% dplyr::mutate(comparison = "PCA_VAR_vs_StateSpace", difference = pca_value - ss_value)

all_differences <- dplyr::bind_rows(
  pca_vs_kof %>% dplyr::select(variable, horizon, comparison, difference),
  ss_vs_kof %>% dplyr::select(variable, horizon, comparison, difference),
  pca_vs_ss %>% dplyr::select(variable, horizon, comparison, difference)
)

cat("\n=== All Pairwise Differences ===\n")
print(all_differences)

# Summary of model performance vs KOF
performance_summary <- all_differences %>%
  dplyr::filter(grepl("_vs_KOF", comparison)) %>%
  dplyr::group_by(comparison, horizon) %>%
  dplyr::summarise(
    mean_abs_diff = mean(abs(difference), na.rm = TRUE),
    median_abs_diff = median(abs(difference), na.rm = TRUE),
    max_abs_diff = max(abs(difference), na.rm = TRUE),
    .groups = 'drop'
  )

cat("\n=== Model Performance Summary (vs KOF) ===\n")
print(performance_summary)

