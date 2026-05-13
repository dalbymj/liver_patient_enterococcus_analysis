# Figure 3: MaAsLin2 bile acid differential abundance analysis
#
# Purpose:
# This script loads plasma bile acid data and clinical metadata, then runs
# MaAsLin2 to identify bile acids associated with ACLD disease severity group.
#
# Input files:
# - data/Corrected bile acid data 01.12.2020.csv
# - data/clinical_metadata_and_data_renamed_30.01.2025.csv
#
# Output files:
# - outputs/Figure_3/Fig3A_maaslin2_bile_acids_all_results.csv
# - outputs/Figure_3/Fig3A_maaslin2_bile_acids_filtered_results.csv
# - outputs/Figure_3/Fig3A_maaslin2_bile_acids_summary.txt
# - outputs/Figure_3/Fig3A_MaAsLin2_bile_acids_inclusion_criteria_output/
#
# Run this script from the root folder of the repository:
# liver_patient_enterococcus_analysis/

library(tidyverse)
library(Maaslin2)

options(scipen = 999)

# Create output folder if it does not already exist
dir.create("outputs/Figure_3", recursive = TRUE, showWarnings = FALSE)

# Define input files
bile_file <- "data/Corrected bile acid data 01.12.2020.csv"
metadata_file <- "data/clinical_metadata_and_data_renamed_30.01.2025.csv"

# Define MaAsLin2 output directory
maaslin2_output_dir <- "outputs/Figure_3/Fig3A_MaAsLin2_bile_acids_inclusion_criteria_output"

# Check that input files exist
if (!file.exists(bile_file)) {
  stop("Missing input file: ", bile_file)
}

if (!file.exists(metadata_file)) {
  stop("Missing input file: ", metadata_file)
}

# Load data
bile_data_raw <- read_csv(bile_file)
metadata <- read_csv(metadata_file)

# Remove unnamed index columns if present
bile_data_raw <- bile_data_raw %>%
  select(-any_of("...1"))

metadata <- metadata %>%
  select(-any_of("...1"))

# Define bile acid columns.
# This assumes that the first bile acid column is 'Deoxycholic_acid'
# and the final bile acid column is 'Taurolithocholic_acid'.
bile_start_col <- match("Deoxycholic_acid", names(bile_data_raw))
bile_end_col <- match("Taurolithocholic_acid", names(bile_data_raw))

if (is.na(bile_start_col)) {
  stop("Could not find bile acid column 'Deoxycholic_acid' in bile_data_raw.")
}

if (is.na(bile_end_col)) {
  stop("Could not find bile acid column 'Taurolithocholic_acid' in bile_data_raw.")
}

bile_cols <- names(bile_data_raw)[bile_start_col:bile_end_col]

# Prepare bile acid data
bile_data <- bile_data_raw %>%
  select(study_no, all_of(bile_cols)) %>%
  mutate(
    across(
      all_of(bile_cols),
      ~ as.numeric(.x)
    )
  ) %>%
  filter(complete.cases(across(all_of(bile_cols))))

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

# Merge metadata and bile acid data
merged_meta_bile_acids <- metadata_subsampled %>%
  inner_join(bile_data, by = "study_no") %>%
  filter(!is.na(Inclusion_criteria))

# Prepare bile acid table for MaAsLin2
bile_acids_table <- merged_meta_bile_acids %>%
  select(study_no, all_of(bile_cols)) %>%
  column_to_rownames("study_no")

# Prepare metadata table for MaAsLin2
metadata_table <- merged_meta_bile_acids %>%
  select(
    study_no,
    Inclusion_criteria
  ) %>%
  column_to_rownames("study_no")

# Check that sample order matches between bile acid and metadata tables
if (!identical(rownames(bile_acids_table), rownames(metadata_table))) {
  stop("Sample order does not match between bile acid table and metadata table.")
}

# Run MaAsLin2
fit_data <- Maaslin2(
  input_data = bile_acids_table,
  input_metadata = metadata_table,
  analysis_method = "LM",
  min_abundance = 0,
  min_prevalence = 0,
  normalization = "CLR",
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
  "outputs/Figure_3/Fig3A_maaslin2_bile_acids_all_results.csv"
)

# Apply the filtering used in the original analysis
maaslin2_filtered_results <- maaslin2_all_results %>%
  filter(qval <= 0.25) %>%
  filter(pval <= 0.05) %>%
  select(
    feature,
    value,
    coef,
    pval,
    qval
  ) %>%
  arrange(coef)

write_csv(
  maaslin2_filtered_results,
  "outputs/Figure_3/Fig3A_maaslin2_bile_acids_filtered_results.csv"
)

# Save summary output
summary_output_file <- "outputs/Figure_3/Fig3A_maaslin2_bile_acids_summary.txt"

sink(summary_output_file)

cat("Figure 3: MaAsLin2 bile acid differential abundance analysis\n")
cat("===========================================================\n\n")

cat("Input files:\n")
cat(bile_file, "\n")
cat(metadata_file, "\n\n")

cat("Sample numbers by group:\n")
print(table(metadata_table$Inclusion_criteria, useNA = "ifany"))

cat("\nNumber of bile acids included:\n")
print(ncol(bile_acids_table))

cat("\nMaAsLin2 settings:\n")
cat("analysis_method = LM\n")
cat("min_abundance = 0\n")
cat("min_prevalence = 0\n")
cat("normalization = CLR\n")
cat("transform = NONE\n")
cat("fixed_effects = Inclusion_criteria\n")
cat("reference = Inclusion_criteria,Healthy control\n\n")

cat("Number of full MaAsLin2 result rows:\n")
print(nrow(maaslin2_all_results))

cat("\nNumber of filtered result rows, qval <= 0.25 and pval <= 0.05:\n")
print(nrow(maaslin2_filtered_results))

cat("\nFiltered MaAsLin2 results:\n")
print(maaslin2_filtered_results)

cat("\nUnique filtered bile acid features:\n")
print(unique(maaslin2_filtered_results$feature))

sink()

message("Full MaAsLin2 results saved to: outputs/Figure_3/Fig3A_maaslin2_bile_acids_all_results.csv")
message("Filtered MaAsLin2 results saved to: outputs/Figure_3/Fig3A_maaslin2_bile_acids_filtered_results.csv")
message("MaAsLin2 output directory: ", maaslin2_output_dir)
message("Summary saved to: ", summary_output_file)