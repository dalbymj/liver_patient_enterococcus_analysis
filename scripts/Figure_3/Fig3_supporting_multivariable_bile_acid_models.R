# Figure 3 supporting analysis: Multivariable bile acid models
#
# Purpose:
# This script tests whether Enterococcus relative abundance is associated with
# selected conjugated primary bile acids after adjustment for antibiotic exposure
# and total bilirubin.
#
# Models:
# - Glycochenodeoxycholic_acid ~ Enterococcus + Receiving_antibiotics + Total_Bilirubin
# - Taurochenodeoxycholic_acid ~ Enterococcus + Receiving_antibiotics + Total_Bilirubin
#
# Input files:
# - data/raw_genus_data_15.04.2024.csv
# - data/clinical_metadata_and_data_renamed_30.01.2025.csv
# - data/Corrected bile acid data 01.12.2020.csv
#
# Output files:
# - outputs/Figure_3/Fig3_supporting_multivariable_bile_acid_model_data.csv
# - outputs/Figure_3/Fig3_supporting_multivariable_bile_acid_model_coefficients.csv
# - outputs/Figure_3/Fig3_supporting_multivariable_bile_acid_model_statistics.csv
# - outputs/Figure_3/Fig3_supporting_multivariable_bile_acid_bestNormalize_transformation_log.csv
# - outputs/Figure_3/Fig3_supporting_multivariable_bile_acid_models_summary.txt
#
# Run this script from the root folder of the repository:
# liver_patient_enterococcus_analysis/

library(tidyverse)
library(bestNormalize)

options(scipen = 999)

set.seed(123)

# Create output folder if it does not already exist
dir.create("outputs/Figure_3", recursive = TRUE, showWarnings = FALSE)

# Define input files
genus_file <- "data/raw_genus_data_15.04.2024.csv"
metadata_file <- "data/clinical_metadata_and_data_renamed_30.01.2025.csv"
bile_file <- "data/Corrected bile acid data 01.12.2020.csv"

# Check that input files exist
if (!file.exists(genus_file)) {
  stop("Missing input file: ", genus_file)
}

if (!file.exists(metadata_file)) {
  stop("Missing input file: ", metadata_file)
}

if (!file.exists(bile_file)) {
  stop("Missing input file: ", bile_file)
}

# Load data
genus_data_raw <- read_csv(genus_file, show_col_types = FALSE)
metadata_raw <- read_csv(metadata_file, show_col_types = FALSE)
bile_data_raw <- read_csv(bile_file, show_col_types = FALSE)

# Remove unnamed index columns if present
genus_data_raw <- genus_data_raw %>%
  select(-any_of("...1"))

metadata_raw <- metadata_raw %>%
  select(-any_of("...1"))

bile_data_raw <- bile_data_raw %>%
  select(-any_of("...1"))

# Define genus columns.
# This assumes that the first genus column is called 'Mucispirillum'
# and that all subsequent columns are genus abundance columns.
genus_start_col <- match("Mucispirillum", names(genus_data_raw))

if (is.na(genus_start_col)) {
  stop("Could not find genus column 'Mucispirillum' in genus_data_raw.")
}

genus_cols <- names(genus_data_raw)[genus_start_col:ncol(genus_data_raw)]

if (!"Enterococcus" %in% genus_cols) {
  stop("Could not find genus column 'Enterococcus' in genus_data_raw.")
}

# Prepare Enterococcus relative abundance data
genus_data <- genus_data_raw %>%
  select(study_no, all_of(genus_cols)) %>%
  mutate(
    across(
      all_of(genus_cols),
      ~ as.numeric(.x)
    )
  ) %>%
  drop_na()

genus_data <- genus_data %>%
  mutate(
    row_total = rowSums(across(all_of(genus_cols)), na.rm = TRUE)
  )

if (any(genus_data$row_total == 0)) {
  stop("At least one sample has a total genus abundance of zero.")
}

enterococcus_data <- genus_data %>%
  mutate(
    Enterococcus = Enterococcus / row_total * 100
  ) %>%
  select(
    study_no,
    Enterococcus
  )

# Prepare clinical metadata
required_metadata_cols <- c(
  "study_no",
  "Receiving_antibiotics",
  "Total_Br_SI_units"
)

missing_metadata_cols <- setdiff(required_metadata_cols, names(metadata_raw))

if (length(missing_metadata_cols) > 0) {
  stop(
    "The following metadata columns are missing: ",
    paste(missing_metadata_cols, collapse = ", ")
  )
}

metadata_subsampled <- metadata_raw %>%
  select(
    study_no,
    Receiving_antibiotics,
    Total_Br_SI_units
  ) %>%
  rename(
    Total_Bilirubin = Total_Br_SI_units
  ) %>%
  mutate(
    Total_Bilirubin = as.numeric(Total_Bilirubin),
    Receiving_antibiotics = as.character(Receiving_antibiotics),
    Receiving_antibiotics = case_when(
      Receiving_antibiotics == "Healthy control" ~ "Healthy_control",
      Receiving_antibiotics %in% c("0", "No", "No_antibiotics") ~ "No_antibiotics",
      Receiving_antibiotics %in% c("1", "Yes", "Receiving_antibiotics", "Reciving_antibiotics") ~ "Receiving_antibiotics",
      TRUE ~ Receiving_antibiotics
    ),
    Receiving_antibiotics = factor(
      Receiving_antibiotics,
      levels = c(
        "No_antibiotics",
        "Receiving_antibiotics",
        "Healthy_control"
      )
    )
  )

# Prepare bile acid data
selected_bile_acids <- c(
  "Glycochenodeoxycholic_acid",
  "Taurochenodeoxycholic_acid"
)

missing_bile_cols <- setdiff(c("study_no", selected_bile_acids), names(bile_data_raw))

if (length(missing_bile_cols) > 0) {
  stop(
    "The following bile acid columns are missing: ",
    paste(missing_bile_cols, collapse = ", ")
  )
}

bile_data <- bile_data_raw %>%
  select(
    study_no,
    all_of(selected_bile_acids)
  ) %>%
  mutate(
    across(
      all_of(selected_bile_acids),
      ~ as.numeric(.x)
    )
  )

# Merge data
model_data_raw <- metadata_subsampled %>%
  inner_join(bile_data, by = "study_no") %>%
  inner_join(enterococcus_data, by = "study_no") %>%
  filter(
    complete.cases(
      across(
        c(
          Receiving_antibiotics,
          Total_Bilirubin,
          Enterococcus,
          all_of(selected_bile_acids)
        )
      )
    )
  ) %>%
  droplevels()

# Check that model data exist
if (nrow(model_data_raw) == 0) {
  stop("No complete-case rows are available for the multivariable models.")
}

# Function to apply bestNormalize to selected continuous variables
apply_best_normalize <- function(data, variables, data_type) {
  
  transformed_data <- data
  transformation_log <- tibble(
    data_type = character(),
    feature = character(),
    transformation = character()
  )
  
  for (feature_name in variables) {
    
    feature_values <- transformed_data[[feature_name]]
    
    if (all(is.na(feature_values)) || sd(feature_values, na.rm = TRUE) == 0) {
      warning("Skipping transformation for ", feature_name, " because values are missing or constant.")
      next
    }
    
    transformed_feature <- bestNormalize::bestNormalize(
      feature_values,
      quiet = TRUE
    )
    
    transformed_data[[feature_name]] <- transformed_feature$x.t
    
    transformation_log <- transformation_log %>%
      add_row(
        data_type = data_type,
        feature = feature_name,
        transformation = class(transformed_feature$chosen_transform)[1]
      )
  }
  
  list(
    data = transformed_data,
    log = transformation_log
  )
}

continuous_variables <- c(
  "Total_Bilirubin",
  "Enterococcus",
  selected_bile_acids
)

transformed_output <- apply_best_normalize(
  data = model_data_raw,
  variables = continuous_variables,
  data_type = "multivariable_model_variable"
)

model_data <- transformed_output$data
transformation_log <- transformed_output$log

write_csv(
  transformation_log,
  "outputs/Figure_3/Fig3_supporting_multivariable_bile_acid_bestNormalize_transformation_log.csv"
)

write_csv(
  model_data,
  "outputs/Figure_3/Fig3_supporting_multivariable_bile_acid_model_data.csv"
)

# Define models
model_specs <- tribble(
  ~model_name, ~outcome,
  "Glycochenodeoxycholic acid model", "Glycochenodeoxycholic_acid",
  "Taurochenodeoxycholic acid model", "Taurochenodeoxycholic_acid"
)

# Function to fit one model
fit_multivariable_model <- function(data, outcome_name) {
  
  model_formula <- as.formula(
    paste0(
      outcome_name,
      " ~ Enterococcus + Receiving_antibiotics + Total_Bilirubin"
    )
  )
  
  lm(
    formula = model_formula,
    data = data
  )
}

# Fit models
fitted_models <- model_specs %>%
  mutate(
    model = map(
      outcome,
      ~ fit_multivariable_model(
        data = model_data,
        outcome_name = .x
      )
    )
  )

# Extract coefficient tables
model_coefficients <- fitted_models %>%
  mutate(
    coefficients = map(
      model,
      ~ summary(.x)$coefficients %>%
        as.data.frame() %>%
        rownames_to_column("term") %>%
        as_tibble()
    )
  ) %>%
  select(model_name, outcome, coefficients) %>%
  unnest(coefficients) %>%
  rename(
    estimate = Estimate,
    standard_error = `Std. Error`,
    t_value = `t value`,
    p_value = `Pr(>|t|)`
  )

write_csv(
  model_coefficients,
  "outputs/Figure_3/Fig3_supporting_multivariable_bile_acid_model_coefficients.csv"
)

# Extract model-level statistics
model_statistics <- fitted_models %>%
  mutate(
    model_summary = map(model, summary),
    n = map_int(model, nobs),
    r_squared = map_dbl(model_summary, ~ .x$r.squared),
    adjusted_r_squared = map_dbl(model_summary, ~ .x$adj.r.squared),
    f_statistic = map_dbl(model_summary, ~ unname(.x$fstatistic["value"])),
    f_num_df = map_dbl(model_summary, ~ unname(.x$fstatistic["numdf"])),
    f_den_df = map_dbl(model_summary, ~ unname(.x$fstatistic["dendf"])),
    model_p_value = pf(
      f_statistic,
      f_num_df,
      f_den_df,
      lower.tail = FALSE
    ),
    residual_standard_error = map_dbl(model_summary, ~ .x$sigma),
    AIC = map_dbl(model, AIC),
    BIC = map_dbl(model, BIC)
  ) %>%
  select(
    model_name,
    outcome,
    n,
    r_squared,
    adjusted_r_squared,
    f_statistic,
    f_num_df,
    f_den_df,
    model_p_value,
    residual_standard_error,
    AIC,
    BIC
  )

write_csv(
  model_statistics,
  "outputs/Figure_3/Fig3_supporting_multivariable_bile_acid_model_statistics.csv"
)

# Save summary output
summary_output_file <- "outputs/Figure_3/Fig3_supporting_multivariable_bile_acid_models_summary.txt"

sink(summary_output_file)

cat("Figure 3 supporting analysis: Multivariable bile acid models\n")
cat("============================================================\n\n")

cat("Input files:\n")
cat(genus_file, "\n")
cat(metadata_file, "\n")
cat(bile_file, "\n\n")

cat("Model formulas:\n")
cat("Glycochenodeoxycholic_acid ~ Enterococcus + Receiving_antibiotics + Total_Bilirubin\n")
cat("Taurochenodeoxycholic_acid ~ Enterococcus + Receiving_antibiotics + Total_Bilirubin\n\n")

cat("Number of complete-case samples included:\n")
print(nrow(model_data))

cat("\nReceiving_antibiotics levels:\n")
print(levels(model_data$Receiving_antibiotics))

cat("\nReceiving_antibiotics counts:\n")
print(table(model_data$Receiving_antibiotics, useNA = "ifany"))

cat("\nContinuous variables transformed with bestNormalize:\n")
print(continuous_variables)

cat("\nTransformation log:\n")
print(transformation_log)

cat("\nModel-level statistics:\n")
print(model_statistics)

cat("\nModel coefficients:\n")
print(model_coefficients)

cat("\nFull model summaries:\n\n")

for (i in seq_len(nrow(fitted_models))) {
  
  cat("------------------------------------------------------------\n")
  cat(fitted_models$model_name[i], "\n")
  cat("Outcome:", fitted_models$outcome[i], "\n")
  cat("------------------------------------------------------------\n\n")
  
  print(summary(fitted_models$model[[i]]))
  
  cat("\n")
}

sink()

message("Model data saved to: outputs/Figure_3/Fig3_supporting_multivariable_bile_acid_model_data.csv")
message("Model coefficients saved to: outputs/Figure_3/Fig3_supporting_multivariable_bile_acid_model_coefficients.csv")
message("Model statistics saved to: outputs/Figure_3/Fig3_supporting_multivariable_bile_acid_model_statistics.csv")
message("Transformation log saved to: outputs/Figure_3/Fig3_supporting_multivariable_bile_acid_bestNormalize_transformation_log.csv")
message("Summary saved to: ", summary_output_file)