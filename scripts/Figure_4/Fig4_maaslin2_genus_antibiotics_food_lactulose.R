# Figure 4: MaAsLin2 genus analysis adjusted for antibiotics, food intake, and lactulose
#
# Purpose:
# This script loads genus-level 16S rRNA sequencing data and clinical metadata,
# converts genus abundances to relative abundance percentages, and runs MaAsLin2
# to test associations between genus abundance and disease severity,
# antibiotic exposure, oral food intake, and lactulose exposure.
#
# Input files:
# - data/normalised_genus_data_16.04.2024.csv
# - data/clinical_metadata_and_data_renamed_30.01.2025.csv
#
# Output files:
# - outputs/Figure_4/Fig4_maaslin2_genus_antibiotics_food_lactulose_results.csv
# - outputs/Figure_4/Fig4_maaslin2_selected_genus_results.csv
# - outputs/Figure_4/Fig4_maaslin2_genus_antibiotics_food_lactulose_summary.txt
# - outputs/Figure_4/MaAsLin2_genus_antibiotics_food_lactulose_output/
#
# Run this script from the root folder of the repository:
# liver_patient_enterococcus_analysis/

library(tidyverse)
library(Maaslin2)

options(scipen = 999)

# Create output folder if it does not already exist
dir.create("outputs/Figure_4", recursive = TRUE, showWarnings = FALSE)

# Define input files
genus_file <- "data/normalised_genus_data_16.04.2024.csv"
metadata_file <- "data/clinical_metadata_and_data_renamed_30.01.2025.csv"

# Check that input files exist
if (!file.exists(genus_file)) {
  stop("Missing input file: ", genus_file)
}

if (!file.exists(metadata_file)) {
  stop("Missing input file: ", metadata_file)
}

# Load data
genus_data_raw <- read_csv(genus_file, show_col_types = FALSE)
metadata_raw <- read_csv(metadata_file, show_col_types = FALSE)

# Remove unnamed index columns if present
genus_data_raw <- genus_data_raw %>%
  select(-any_of("...1"))

metadata_raw <- metadata_raw %>%
  select(-any_of("...1"))

# Define genus columns.
# This assumes that the first genus column is called 'Mucispirillum'
# and that all subsequent columns are genus abundance columns.
genus_start_col <- match("Mucispirillum", names(genus_data_raw))

if (is.na(genus_start_col)) {
  stop("Could not find genus column 'Mucispirillum' in genus_data_raw.")
}

genus_cols <- names(genus_data_raw)[genus_start_col:ncol(genus_data_raw)]

# Prepare genus abundance data
genus_data <- genus_data_raw %>%
  select(study_no, all_of(genus_cols)) %>%
  mutate(
    across(
      all_of(genus_cols),
      ~ as.numeric(.x)
    )
  ) %>%
  drop_na(all_of(genus_cols))

# Convert genus abundances to relative abundance percentage
genus_data <- genus_data %>%
  mutate(
    row_total = rowSums(across(all_of(genus_cols)), na.rm = TRUE)
  )

if (any(genus_data$row_total == 0)) {
  stop("At least one sample has a total genus abundance of zero.")
}

genus_data_percent <- genus_data %>%
  mutate(
    across(
      all_of(genus_cols),
      ~ .x / row_total * 100
    )
  ) %>%
  select(-row_total)

# Prepare metadata
metadata_subsampled <- metadata_raw %>%
  select(
    study_no,
    Inclusion_criteria,
    Receiving_antibiotics,
    Oral_intake,
    Lactulose
  ) %>%
  mutate(
    Inclusion_criteria = as.character(Inclusion_criteria),
    Inclusion_criteria = dplyr::recode(
      Inclusion_criteria,
      "AoCLF" = "ACLF",
      .default = Inclusion_criteria
    ),
    Inclusion_criteria = factor(
      Inclusion_criteria,
      levels = c(
        "Healthy control",
        "Stable cirrhotic",
        "Decomp cirrhotic",
        "ACLF"
      )
    ),
    
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
    ),
    
    Oral_intake = as.character(Oral_intake),
    Oral_intake = case_when(
      Oral_intake == "Healthy control" ~ "Healthy_control",
      Oral_intake %in% c("0", "Yes", "Oral_food_intake") ~ "Yes",
      Oral_intake %in% c("1", "No", "No_oral_food_intake") ~ "No",
      TRUE ~ Oral_intake
    ),
    Oral_intake = factor(
      Oral_intake,
      levels = c(
        "Yes",
        "No",
        "Healthy_control"
      )
    ),
    
    Lactulose = as.character(Lactulose),
    Lactulose = case_when(
      Lactulose == "Healthy control" ~ "Healthy_control",
      Lactulose %in% c("0", "No", "No_lactulose") ~ "No",
      Lactulose %in% c("1", "Yes", "Lactulose") ~ "Yes",
      TRUE ~ Lactulose
    ),
    Lactulose = factor(
      Lactulose,
      levels = c(
        "No",
        "Yes",
        "Healthy_control"
      )
    )
  ) %>%
  drop_na(
    Inclusion_criteria,
    Receiving_antibiotics,
    Oral_intake,
    Lactulose
  )

# Merge metadata and genus data
merged_meta_seqs <- metadata_subsampled %>%
  inner_join(genus_data_percent, by = "study_no")

# Prepare MaAsLin2 input tables
sequence_table <- merged_meta_seqs %>%
  select(study_no, all_of(genus_cols)) %>%
  column_to_rownames("study_no")

# Remove genera with zero abundance across all samples
sequence_table <- sequence_table[, colSums(sequence_table != 0) > 0, drop = FALSE]

metadata_table <- merged_meta_seqs %>%
  select(
    study_no,
    Inclusion_criteria,
    Receiving_antibiotics,
    Oral_intake,
    Lactulose
  ) %>%
  rename(
    Antibiotics = Receiving_antibiotics,
    Food_intake = Oral_intake
  ) %>%
  column_to_rownames("study_no")

# Check that sample order matches
if (!identical(rownames(sequence_table), rownames(metadata_table))) {
  stop("Sample order does not match between sequence_table and metadata_table.")
}

# MaAsLin2 output directory
maaslin2_output_dir <- "outputs/Figure_4/MaAsLin2_genus_antibiotics_food_lactulose_output"

dir.create(
  maaslin2_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Run MaAsLin2
fit_data <- Maaslin2(
  input_data = sequence_table,
  input_metadata = metadata_table,
  analysis_method = "LM",
  min_abundance = 0.01,
  min_prevalence = 0.01,
  normalization = "NONE",
  output = maaslin2_output_dir,
  fixed_effects = c(
    "Inclusion_criteria",
    "Antibiotics",
    "Food_intake",
    "Lactulose"
  ),
  reference = c(
    "Inclusion_criteria,Healthy control",
    "Antibiotics,No_antibiotics",
    "Food_intake,Yes",
    "Lactulose,No"
  )
)

# Save MaAsLin2 results
maaslin2_results <- fit_data$results %>%
  as_tibble() %>%
  arrange(coef)

write_csv(
  maaslin2_results,
  "outputs/Figure_4/Fig4_maaslin2_genus_antibiotics_food_lactulose_results.csv"
)

# Save selected results relevant to Figure 4A-E panels
selected_genus_results <- maaslin2_results %>%
  filter(
    feature %in% c(
      "Enterococcus",
      "Ruminococcus",
      "Roseburia",
      "Bifidobacterium"
    )
  ) %>%
  select(
    feature,
    metadata,
    value,
    coef,
    stderr,
    pval,
    qval
  )

write_csv(
  selected_genus_results,
  "outputs/Figure_4/Fig4_maaslin2_selected_genus_results.csv"
)

# Save summary output
summary_output_file <- "outputs/Figure_4/Fig4_maaslin2_genus_antibiotics_food_lactulose_summary.txt"

sink(summary_output_file)

cat("Figure 4: MaAsLin2 genus analysis adjusted for antibiotics, food intake, and lactulose\n")
cat("===================================================================================\n\n")

cat("Input files:\n")
cat(genus_file, "\n")
cat(metadata_file, "\n\n")

cat("Number of samples included:\n")
print(nrow(metadata_table))

cat("\nNumber of genera included after zero-sum filtering:\n")
print(ncol(sequence_table))

cat("\nSample numbers by disease severity group:\n")
print(table(metadata_table$Inclusion_criteria, useNA = "ifany"))

cat("\nSample numbers by antibiotic exposure:\n")
print(table(metadata_table$Antibiotics, useNA = "ifany"))

cat("\nSample numbers by food intake:\n")
print(table(metadata_table$Food_intake, useNA = "ifany"))

cat("\nSample numbers by lactulose exposure:\n")
print(table(metadata_table$Lactulose, useNA = "ifany"))

cat("\nMaAsLin2 fixed effects:\n")
print(c("Inclusion_criteria", "Antibiotics", "Food_intake", "Lactulose"))

cat("\nMaAsLin2 reference levels:\n")
print(c(
  "Inclusion_criteria,Healthy control",
  "Antibiotics,No_antibiotics",
  "Food_intake,Yes",
  "Lactulose,No"
))

cat("\nSelected Figure 4A-E genus results:\n")
print(selected_genus_results)

sink()

message("MaAsLin2 results saved to: outputs/Figure_4/Fig4_maaslin2_genus_antibiotics_food_lactulose_results.csv")
message("Selected genus results saved to: outputs/Figure_4/Fig4_maaslin2_selected_genus_results.csv")
message("Summary saved to: ", summary_output_file)