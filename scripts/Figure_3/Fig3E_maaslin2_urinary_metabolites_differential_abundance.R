# Figure 3E: MaAsLin2 urinary metabolite differential abundance analysis
#
# Purpose:
# This script loads urinary NMR metabolite data and clinical metadata, applies
# bestNormalize transformation to urinary metabolites, and runs MaAsLin2 to
# identify urinary metabolites associated with ACLD disease severity group.
#
# Input files:
# - data/Updated_urinary_NMR_GLA_RifSys_8.03.2021 with edited sample_ID.csv
# - data/clinical_metadata_and_data_renamed_30.01.2025.csv
#
# Output files:
# - outputs/Figure_3/Fig3E_maaslin2_urinary_metabolites_all_results.csv
# - outputs/Figure_3/Fig3E_maaslin2_urinary_metabolites_filtered_results.csv
# - outputs/Figure_3/Fig3E_urinary_metabolites_bestNormalize_transformation_log.csv
# - outputs/Figure_3/Fig3E_maaslin2_urinary_metabolites_summary.txt
# - outputs/Figure_3/MaAsLin2_urinary_metabolites_inclusion_criteria_output/
#
# Run this script from the root folder of the repository:
# liver_patient_enterococcus_analysis/

library(tidyverse)
library(bestNormalize)
library(Maaslin2)

options(scipen = 999)

set.seed(123)

# Create output folder if it does not already exist
dir.create("outputs/Figure_3", recursive = TRUE, showWarnings = FALSE)

# Define input files
urinary_file <- "data/Urinary_NMR_GLA_RifSys_8.03.2021.csv"
metadata_file <- "data/clinical_metadata_and_data_renamed_30.01.2025.csv"

# Define MaAsLin2 output directory
maaslin2_output_dir <- "outputs/Figure_3/MaAsLin2_urinary_metabolites_inclusion_criteria_output"

# Check that input files exist
if (!file.exists(urinary_file)) {
  stop("Missing input file: ", urinary_file)
}

if (!file.exists(metadata_file)) {
  stop("Missing input file: ", metadata_file)
}

# Load data
metabolite_data_raw <- read_csv(urinary_file, show_col_types = FALSE)
metadata <- read_csv(metadata_file, show_col_types = FALSE)

# Remove unnamed index columns if present
metabolite_data_raw <- metabolite_data_raw %>%
  select(-any_of("...1"))

metadata <- metadata %>%
  select(-any_of("...1"))

# Standardise sample ID column name
if ("Study no" %in% names(metabolite_data_raw)) {
  names(metabolite_data_raw)[names(metabolite_data_raw) == "Study no"] <- "study_no"
}

# Replace "/" in metabolite names with "_", matching the original script
names(metabolite_data_raw) <- gsub(
  x = names(metabolite_data_raw),
  pattern = "/",
  replacement = "_"
)

# Define urinary metabolite columns
metabolite_start_col <- match("For_TS", names(metabolite_data_raw))
metabolite_end_col <- match("Lac_TS", names(metabolite_data_raw))

if (is.na(metabolite_start_col)) {
  stop("Could not find urinary metabolite column 'For_TS' in metabolite_data_raw.")
}

if (is.na(metabolite_end_col)) {
  stop("Could not find urinary metabolite column 'Lac_TS' in metabolite_data_raw.")
}

metabolite_cols <- names(metabolite_data_raw)[metabolite_start_col:metabolite_end_col]

# Prepare urinary metabolite data
metabolites <- metabolite_data_raw %>%
  select(study_no, all_of(metabolite_cols)) %>%
  mutate(
    across(
      all_of(metabolite_cols),
      ~ as.numeric(.x)
    )
  ) %>%
  drop_na(all_of(metabolite_cols))

# Prepare metadata
metadata_subsampled <- metadata %>%
  select(
    study_no,
    Inclusion_criteria
  ) %>%
  mutate(
    Inclusion_criteria = factor(
      Inclusion_criteria,
      levels = c(
        "Healthy control",
        "Stable cirrhotic",
        "Decomp cirrhotic",
        "ACLF"
      )
    )
  )

# Merge metadata and urinary metabolite data
merged_meta_metabolites <- metadata_subsampled %>%
  inner_join(metabolites, by = "study_no") %>%
  filter(!is.na(Inclusion_criteria))

# Prepare urinary metabolite table
metabolite_table_raw <- merged_meta_metabolites %>%
  select(study_no, all_of(metabolite_cols)) %>%
  column_to_rownames("study_no")

# Function to apply bestNormalize to each metabolite
apply_best_normalize <- function(data_table, data_type) {
  
  transformed_table <- data_table
  transformation_log <- tibble(
    data_type = character(),
    feature = character(),
    transformation = character()
  )
  
  for (feature_name in names(data_table)) {
    
    feature_values <- data_table[[feature_name]]
    
    if (all(is.na(feature_values)) || sd(feature_values, na.rm = TRUE) == 0) {
      warning("Skipping transformation for ", feature_name, " because values are missing or constant.")
      next
    }
    
    transformed_feature <- bestNormalize::bestNormalize(
      feature_values,
      quiet = TRUE
    )
    
    transformed_table[[feature_name]] <- transformed_feature$x.t
    
    transformation_log <- transformation_log %>%
      add_row(
        data_type = data_type,
        feature = feature_name,
        transformation = class(transformed_feature$chosen_transform)[1]
      )
  }
  
  list(
    data = transformed_table,
    log = transformation_log
  )
}

# Transform urinary metabolites
metabolite_transformed <- apply_best_normalize(
  data_table = metabolite_table_raw,
  data_type = "urinary_metabolite"
)

metabolite_table <- metabolite_transformed$data
transformation_log <- metabolite_transformed$log

write_csv(
  transformation_log,
  "outputs/Figure_3/Fig3E_urinary_metabolites_bestNormalize_transformation_log.csv"
)

# Prepare metadata table for MaAsLin2
metadata_table <- merged_meta_metabolites %>%
  select(
    study_no,
    Inclusion_criteria
  ) %>%
  column_to_rownames("study_no")

# Check that sample order matches between metabolite and metadata tables
if (!identical(rownames(metabolite_table), rownames(metadata_table))) {
  stop("Sample order does not match between metabolite table and metadata table.")
}

# Run MaAsLin2
fit_data <- Maaslin2(
  input_data = metabolite_table,
  input_metadata = metadata_table,
  analysis_method = "LM",
  min_abundance = 0,
  min_prevalence = 0,
  normalization = "NONE",
  transform = "NONE",
  output = maaslin2_output_dir,
  fixed_effects = c("Inclusion_criteria"),
  reference = c("Inclusion_criteria,Healthy control")
)

# Save full MaAsLin2 results
maaslin2_all_results <- fit_data$results %>%
  as.data.frame() %>%
  arrange(coef)

write_csv(
  maaslin2_all_results,
  "outputs/Figure_3/Fig3E_maaslin2_urinary_metabolites_all_results.csv"
)

# Save filtered results for convenience
maaslin2_filtered_results <- maaslin2_all_results %>%
  filter(qval <= 0.25) %>%
  filter(pval <= 0.05) %>%
  arrange(coef)

write_csv(
  maaslin2_filtered_results,
  "outputs/Figure_3/Fig3E_maaslin2_urinary_metabolites_filtered_results.csv"
)

# Save summary output
summary_output_file <- "outputs/Figure_3/Fig3E_maaslin2_urinary_metabolites_summary.txt"

sink(summary_output_file)

cat("Figure 3E: MaAsLin2 urinary metabolite differential abundance analysis\n")
cat("======================================================================\n\n")

cat("Input files:\n")
cat(urinary_file, "\n")
cat(metadata_file, "\n\n")

cat("Sample numbers by group:\n")
print(table(metadata_table$Inclusion_criteria, useNA = "ifany"))

cat("\nNumber of urinary metabolites included:\n")
print(ncol(metabolite_table))

cat("\nMaAsLin2 settings:\n")
cat("analysis_method = LM\n")
cat("min_abundance = 0\n")
cat("min_prevalence = 0\n")
cat("normalization = NONE\n")
cat("transform = NONE\n")
cat("fixed_effects = Inclusion_criteria\n")
cat("reference = Inclusion_criteria,Healthy control\n\n")

cat("bestNormalize transformation log:\n")
print(transformation_log)

cat("\nNumber of full MaAsLin2 result rows:\n")
print(nrow(maaslin2_all_results))

cat("\nNumber of filtered result rows, qval <= 0.25 and pval <= 0.05:\n")
print(nrow(maaslin2_filtered_results))

cat("\nFiltered MaAsLin2 results:\n")
print(maaslin2_filtered_results)

cat("\nSelected urinary metabolites used in Figure 3F-G:\n")
print(
  maaslin2_all_results %>%
    filter(feature %in% c("TMAO_TS", "Hip7p85_TS")) %>%
    select(any_of(c("feature", "metadata", "value", "coef", "pval", "qval"))) %>%
    arrange(feature, value)
)

sink()

message("Full MaAsLin2 results saved to: outputs/Figure_3/Fig3E_maaslin2_urinary_metabolites_all_results.csv")
message("Filtered MaAsLin2 results saved to: outputs/Figure_3/Fig3E_maaslin2_urinary_metabolites_filtered_results.csv")
message("Transformation log saved to: outputs/Figure_3/Fig3E_urinary_metabolites_bestNormalize_transformation_log.csv")
message("MaAsLin2 output directory: ", maaslin2_output_dir)
message("Summary saved to: ", summary_output_file)