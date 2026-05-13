# Figure 2B: MaAsLin2 genus differential abundance analysis
#
# Purpose:
# This script loads genus-level 16S rRNA sequencing data and clinical metadata,
# converts genus abundances to relative abundance percentages, and runs MaAsLin2
# to identify bacterial genera associated with ACLD disease severity group.
#
# Input files:
# - data/normalised_genus_data_16.04.2024.csv
# - data/clinical_metadata_and_data_renamed_30.01.2025.csv
#
# Output files:
# - outputs/Figure_2/Fig2B_maaslin2_genus_differential_abundance_results.csv
# - outputs/Figure_2/Fig2B_maaslin2_genus_differential_abundance_summary.txt
# - outputs/Figure_2/MaAsLin2_genus_inclusion_criteria_output/
#
# Run this script from the root folder of the repository:
# liver_patient_enterococcus_analysis/

library(tidyverse)
library(Maaslin2)

options(scipen = 999)

# Create output folder if it does not already exist
dir.create("outputs/Figure_2", recursive = TRUE, showWarnings = FALSE)

# Define input files
seq_file <- "data/normalised_genus_data_16.04.2024.csv"
metadata_file <- "data/clinical_metadata_and_data_renamed_30.01.2025.csv"

# Define MaAsLin2 output directory
maaslin2_output_dir <- "outputs/Figure_2/MaAsLin2_genus_inclusion_criteria_output"

# Check that input files exist
if (!file.exists(seq_file)) {
  stop("Missing input file: ", seq_file)
}

if (!file.exists(metadata_file)) {
  stop("Missing input file: ", metadata_file)
}

# Load data
seq_data <- read_csv(seq_file)
metadata <- read_csv(metadata_file)

# Remove unnamed index columns if present
seq_data <- seq_data %>%
  select(-any_of("...1"))

metadata <- metadata %>%
  select(-any_of("...1"))

# Define genus columns.
# This assumes that the first genus column is called 'Mucispirillum'
# and that all subsequent columns are genus abundance columns.
genus_start_col <- match("Mucispirillum", names(seq_data))

if (is.na(genus_start_col)) {
  stop("Could not find genus column 'Mucispirillum' in seq_data.")
}

genus_cols <- names(seq_data)[genus_start_col:ncol(seq_data)]

# Prepare genus abundance data
seq_data_genus <- seq_data %>%
  select(study_no, all_of(genus_cols)) %>%
  mutate(
    across(
      all_of(genus_cols),
      ~ as.numeric(.x)
    )
  ) %>%
  drop_na()

# Check that all samples have non-zero genus abundance
seq_data_genus <- seq_data_genus %>%
  mutate(
    row_total = rowSums(across(all_of(genus_cols)), na.rm = TRUE)
  )

if (any(seq_data_genus$row_total == 0)) {
  stop("At least one sample has a total genus abundance of zero.")
}

# Convert genus abundances to relative abundance percentage
seq_data_genus_percent <- seq_data_genus %>%
  mutate(
    across(
      all_of(genus_cols),
      ~ .x / row_total * 100
    )
  ) %>%
  select(-row_total)

# Remove genera with zero abundance across all samples
non_zero_genera <- genus_cols[colSums(seq_data_genus_percent[, genus_cols] != 0) > 0]

seq_data_genus_percent <- seq_data_genus_percent %>%
  select(study_no, all_of(non_zero_genera))

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

# Merge metadata and genus abundance data
merged_meta_seqs <- metadata_subsampled %>%
  inner_join(seq_data_genus_percent, by = "study_no") %>%
  filter(!is.na(Inclusion_criteria))

# Prepare sequence table for MaAsLin2
sequence_table <- merged_meta_seqs %>%
  select(study_no, all_of(non_zero_genera)) %>%
  column_to_rownames("study_no")

# Prepare metadata table for MaAsLin2
metadata_table <- merged_meta_seqs %>%
  select(
    study_no,
    Inclusion_criteria
  ) %>%
  column_to_rownames("study_no")

# Check that sample order matches between sequence and metadata tables
if (!identical(rownames(sequence_table), rownames(metadata_table))) {
  stop("Sample order does not match between sequence table and metadata table.")
}

# Run MaAsLin2
fit_data <- Maaslin2(
  input_data = sequence_table,
  input_metadata = metadata_table,
  analysis_method = "LM",
  min_abundance = 0.01,
  min_prevalence = 0.01,
  normalization = "NONE",
  output = maaslin2_output_dir,
  fixed_effects = c("Inclusion_criteria"),
  reference = c("Inclusion_criteria,Healthy control")
)

# Extract and save MaAsLin2 results
maaslin2_results <- fit_data$results %>%
  as.data.frame() %>%
  arrange(coef)

write_csv(
  maaslin2_results,
  "outputs/Figure_2/Fig2B_maaslin2_genus_differential_abundance_results.csv"
)

# Save summary output
summary_output_file <- "outputs/Figure_2/Fig2B_maaslin2_genus_differential_abundance_summary.txt"

sink(summary_output_file)

cat("Figure 2B: MaAsLin2 genus differential abundance analysis\n")
cat("=========================================================\n\n")

cat("Input files:\n")
cat(seq_file, "\n")
cat(metadata_file, "\n\n")

cat("Sample numbers by group:\n")
print(table(metadata_table$Inclusion_criteria, useNA = "ifany"))

cat("\nNumber of genera included after zero-sum filtering:\n")
print(ncol(sequence_table))

cat("\nMaAsLin2 settings:\n")
cat("analysis_method = LM\n")
cat("min_abundance = 0.01\n")
cat("min_prevalence = 0.01\n")
cat("normalization = NONE\n")
cat("fixed_effects = Inclusion_criteria\n")
cat("reference = Inclusion_criteria,Healthy control\n\n")

cat("Number of MaAsLin2 result rows:\n")
print(nrow(maaslin2_results))

cat("\nTop negative coefficients:\n")
print(
  maaslin2_results %>%
    arrange(coef) %>%
    select(feature, value, coef, pval, qval) %>%
    head(20)
)

cat("\nTop positive coefficients:\n")
print(
  maaslin2_results %>%
    arrange(desc(coef)) %>%
    select(feature, value, coef, pval, qval) %>%
    head(20)
)

cat("\nSelected genera used in Figure 2C-H:\n")
print(
  maaslin2_results %>%
    filter(
      feature %in% c(
        "Ruminococcus",
        "Roseburia",
        "Lachnospira",
        "Streptococcus",
        "Staphylococcus",
        "Enterococcus"
      )
    ) %>%
    select(feature, value, coef, pval, qval) %>%
    arrange(feature, value)
)

sink()

message("MaAsLin2 results saved to: outputs/Figure_2/Fig2B_maaslin2_genus_differential_abundance_results.csv")
message("MaAsLin2 output directory: ", maaslin2_output_dir)
message("Summary saved to: ", summary_output_file)