renv::status()
# Example R script to test the Docker setup
library(ggplot2)
library(dplyr)

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

metadata <- utils::read.csv("data/metadata_quarterly_en.csv")
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
