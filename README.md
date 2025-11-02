# Methods of Macroeconomic Forecasting - Autumn 2025
Course repository for the ETH autumn 2025 course on "Methods of Macroeconomic Forecasting"

## General Information

- Lectures: [Dr. Samad Sarferaz](https://sites.google.com/site/samadsarferaz/home)
- Labs: Merlin Scherer, Lena Will
- [Tech Setup Checklist](/resources/tech_setup_checklist.pdf)
- [Project Description](/resources/project_description.pdf)

## Time Schedule

| Date | Time | Schedule |
|------|------| ---|
|02 Oct 2025 |10am - 3pm |Methods Boot Camp|
|02 Oct 2025 |3pm - 4pm |Lab|
|03 Oct 2025 |10am - 4pm |Methods Boot Camp with integrated Lab|
|10 Oct 2025 | EOB | Hand in your project preferences |
|14 Oct 2025 | EOB | Topic allocation |
|19 Nov 2025 | 11pm | Hand in your presentation slides |
|20 Nov 2025 |10am - 2pm |Presentations|
|21 Nov 2025 |10am - 2pm |Presentations|
|12 Dec 2025 |10am - 12pm |Exam|

# Documentation

### 1. Load Data
Straight forward. Install dependencies and load data

### 2. Data Preparation

* Converts metadata names to lowercase, ensures `ts_key` exists, and normalizes `unit` values for mapping
* Parses dates, filters pre-Oct 2025 data, coerces numeric types, and maintains temporal order
* Maps measurement units to transformations (`logdiff` for proportional growth, `diff` for absolute or rate data)
* Applies `transform_variable()` to each series, removing first-row NAs and dropping all-NA or zero-variance columns
* Standardizes variables with z-scores and saves scaling parameters for later rescaling
* Defines targets (`rvgdp`, `cpi`, `wkfreuro`) and predictors, excluding variables with over 30 % missing values
* Runs Augmented Dickey–Fuller tests to confirm stationarity of transformed series
* Splits data chronologically into training (≈ 80 %) and expanding out-of-sample windows with 1- and 4-quarter horizons

### 3. Bai Ng

* Defines bn_select_r_icp2() to choose the optimal number of factors using the ICp2 criterion
* Demeans data, runs SVD, and computes ICp2 across candidate ranks
* Imputes missing values in X0 with column means before selection
* Selects and reports the rank r_fixed minimizing ICp2 (capped at 10)

_Currently selecting 10_

### 4. Forecasting with PCA + VAR (Dynamic Factor Model)

* Initializes matrices to store 1-step and 4-step forecasts (both standardized and real scale)
* Iterates over each origin in an expanding-window setup, updating the training sample each time
* Standardizes predictors (X_train) and targets (Y_train) within each window, imputing missing values
* Extracts factors via PCA (r_fixed factors) and fits a VAR model (p_fix lags) on the factors
* Uses the VAR forecasts to predict future factors (h=1 and h=4)
* Estimates target forecasts by regressing Y on current and lagged factors
* Transforms predictions back to the real scale using window means and standard deviations
* Stores both standardized and real forecasts, along with corresponding actual values for evaluation

### 5. OOS metrics + Back-transformation to levels

* Verifies required objects exist before running evaluation to prevent runtime errors
* Defines RMSE and MAE helper functions for forecast accuracy metrics
* Computes RMSE and MAE in transformed units for both 1-step and 4-step horizons
* Retrieves each variable’s transformation type from metadata and reconstructs forecasts back to level scale
* Uses exponential transformation for logdiff and addition for diff to restore original values
* Calculates RMSE and MAE again in levels to assess real-world accuracy
* Combines all forecast outputs, actuals, and dates into a final oos_results dataframe containing both transformed and level-based forecasts

| variable |    RMSE_h1    |   RMSE_h4    |   MAE_h1    |   MAE_h4    |
|:----------|--------------:|-------------:|-------------:|-------------:|
| rvgdp     | 2903.4894353  | 11067.96     | 2176.0926407 | 4261.3372818 |
| cpi       | 0.8835048     | 1.682424     | 0.6185843    | 0.7090984    |
| wkfreuro  | 0.5412938     | 0.8717933    | 0.4018382    | 0.3010463    |

### 6. Plot — OOS transformed scale (model units)
* Prepares separate datasets for 1-quarter (h=1) and 1-year (h=4) horizons in both transformed and back-transformed levels
* Creates two main plots: Transformed scale: “OOS Forecasts — Transformed Units”, Level scale: “OOS Forecasts — Back-transformed Levels”

![Transformed Units](../images/6_transformed.png)
![Back Transformed Units](../images/6_backtransformed.png)

### 7. NOW Forecasts

* Filters the cleaned dataset to include only observed data up to July 2025
* Builds a standardized design matrix (`X_std`) and target matrix (`Y_std`) using full-history means and standard deviations
* Imputes missing standardized predictor values with zero (the variable mean)
* Selects number of factors (`r_now`) using `r_fixed` or re-estimates via Bai–Ng IC if absent
* Performs PCA on standardized predictors to extract factors and runs a VAR(2) model to forecast them to horizons h = 1–4
* Fits OLS regressions of targets on current and lagged factors to estimate mapping from factors to outputs
* Generates standardized and re-scaled forecasts (`Y_h1_real`, `Y_h4_real`) for near-term and one-year horizons
* Recursively back-transforms forecasts to original levels using last observed values and transformation rules (`logdiff` → exponential growth, `diff` → cumulative sums)
* Produces final nowcasting tables in both transformed and level scales (`now_table_transformed`, `now_table_levels`) with forecasted values for `rvgdp`, `cpi`, and `wkfreuro`

Difference table
| variable | target_h1_date | forecast_h1  | target_h4_date | forecast_h4  |
| -------- | -------------- | ------------ | -------------- | ------------ |
| rvgdp    | 2025-10-01     | -0.053006777 | 2026-07-01     | 0.029964283  |
| cpi      | 2025-10-01     | 0.077652466  | 2026-07-01     | 0.130492495  |
| wkfreuro | 2025-10-01     | -0.006250239 | 2026-07-01     | -0.008142884 |

Forecast (levels) table
| variable | target_h1_date | forecast_h1_level | target_h4_date | forecast_h4_level |
| -------- | -------------- | ----------------- | -------------- | ----------------- |
| rvgdp    | 2025-10-01     | 196566.7          | 2026-07-01     | 199853.7          |
| cpi      | 2025-10-01     | 107.5820          | 2026-07-01     | 108.1014          |
| wkfreuro | 2025-10-01     | 0.9290944         | 2026-07-01     | 0.9135523         |

### 8. NOW: PLOT VS KOF forecasts
* Filters recent observed data (from 2022) and extracts last observation date
* Builds comparison points from our model forecasts (`now_table_levels`) and KOF values at h=1 and h=4 horizons
* Retrieves matching KOF forecasts from the main dataset for the same target dates
* Combines historical series, our forecasts, and KOF forecasts into unified plotting data
* Plots historical levels as grey lines, adds forecast points for both sources, and highlights last observation with a vertical line
* Facets plots by variable (Real GDP, CPI, Exchange Rate) with consistent color and shape legends
* Produces a comparison chart titled “Forecast Comparison: Dynamic Factor Model (PCA+VAR) vs KOF”

![Comparison Chart](../images/8_comparison.png)

### 9. STATE SPACE: Dynamic Factor Model Initialization (PCA/SVD → F0, Lambda0, R0, A0, Q0)

* Workflow continues from PCA+VAR: after cleaning, transforming, standardizing, and selecting factors (`r_fixed`), we move to a **state-space DFM**
* Runs SVD on standardized predictors to get initial factors and loadings, then estimates AR(1) dynamics and variances for each factor
* Builds a **KFAS state-space model** (`SSModel`) linking factors to observed data
* Optionally smooths factors and fits model parameters (H, Q, A) via quasi-ML using BFGS
* Extracts estimated components and diagnostics (logLik, AIC, BIC, residuals)
* Forecasts h steps ahead by extending the model and filtering forward
* Maps factor forecasts to target variables, de-standardizes, and back-transforms to levels
* Outputs a simple forecast table (`now_table_transformed_ssm`) for the next quarters

* 