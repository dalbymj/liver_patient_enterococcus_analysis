# Figure 2H: Post-hoc power analysis for Enterococcus abundance
#
# Purpose:
# This script loads genus-level 16S rRNA sequencing data and clinical metadata,
# calculates Enterococcus relative abundance by disease severity group, estimates
# the observed ANOVA effect size, and performs a post-hoc power analysis.
#
# Input files:
# - data/normalised_genus_data_16.04.2024.csv
# - data/clinical_metadata_and_data_renamed_30.01.2025.csv
#
# Output files:
# - outputs/Figure_2/Fig2H_enterococcus_posthoc_power_analysis_summary.csv
# - outputs/Figure_2/Fig2H_enterococcus_posthoc_power_analysis.txt
#
# Run this script from the root folder of the repository:
# liver_patient_enterococcus_analysis/

library(tidyverse)
library(pwr)

options(scipen = 999)

# Create output folder if it does not already exist
dir.create("outputs/Figure_2", recursive = TRUE, showWarnings = FALSE)

# Define input files
seq_file <- "data/normalised_genus_data_16.04.2024.csv"
metadata_file <- "data/clinical_metadata_and_data_renamed_30.01.2025.csv"

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

if (!"Enterococcus" %in% genus_cols) {
  stop("Could not find genus column 'Enterococcus' in seq_data.")
}

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

# Merge metadata and Enterococcus abundance data
power_analysis_data <- metadata_subsampled %>%
  inner_join(
    seq_data_genus_percent %>%
      select(study_no, Enterococcus),
    by = "study_no"
  ) %>%
  filter(!is.na(Inclusion_criteria))

# Summary by group
summary_table <- power_analysis_data %>%
  group_by(Inclusion_criteria) %>%
  summarise(
    n = n(),
    mean_enterococcus = mean(Enterococcus, na.rm = TRUE),
    sd_enterococcus = sd(Enterococcus, na.rm = TRUE),
    median_enterococcus = median(Enterococcus, na.rm = TRUE),
    IQR_enterococcus = IQR(Enterococcus, na.rm = TRUE),
    max_enterococcus = max(Enterococcus, na.rm = TRUE),
    .groups = "drop"
  )

# One-way ANOVA used to estimate observed effect size
res_aov <- aov(Enterococcus ~ Inclusion_criteria, data = power_analysis_data)
aov_summary <- summary(res_aov)

# Extract sums of squares
ss_between <- aov_summary[[1]][["Sum Sq"]][1]
ss_within <- aov_summary[[1]][["Sum Sq"]][2]

partial_eta_squared <- ss_between / (ss_between + ss_within)
cohens_f <- sqrt(partial_eta_squared / (1 - partial_eta_squared))

# Harmonic mean sample size for unequal groups
group_sample_sizes <- summary_table$n
number_of_groups <- length(group_sample_sizes)
harmonic_mean_n <- number_of_groups / sum(1 / group_sample_sizes)

# Observed power based on harmonic mean sample size
observed_power <- pwr.anova.test(
  k = number_of_groups,
  n = harmonic_mean_n,
  f = cohens_f,
  sig.level = 0.05
)$power

# Minimum detectable Cohen's f for 80% and 90% power
detectable_f_80 <- uniroot(
  function(effect_size) {
    pwr.anova.test(
      k = number_of_groups,
      n = harmonic_mean_n,
      f = effect_size,
      sig.level = 0.05
    )$power - 0.80
  },
  c(1e-6, 5)
)$root

detectable_f_90 <- uniroot(
  function(effect_size) {
    pwr.anova.test(
      k = number_of_groups,
      n = harmonic_mean_n,
      f = effect_size,
      sig.level = 0.05
    )$power - 0.90
  },
  c(1e-6, 5)
)$root

# Save summary table
power_summary <- tibble(
  number_of_groups = number_of_groups,
  group_sample_sizes = paste(group_sample_sizes, collapse = ", "),
  harmonic_mean_n = harmonic_mean_n,
  ss_between = ss_between,
  ss_within = ss_within,
  partial_eta_squared = partial_eta_squared,
  cohens_f = cohens_f,
  observed_power = observed_power,
  detectable_f_80_percent_power = detectable_f_80,
  detectable_f_90_percent_power = detectable_f_90
)

write_csv(
  summary_table,
  "outputs/Figure_2/Fig2H_enterococcus_posthoc_power_group_summary.csv"
)

write_csv(
  power_summary,
  "outputs/Figure_2/Fig2H_enterococcus_posthoc_power_analysis_summary.csv"
)

# Save detailed text output
power_output_file <- "outputs/Figure_2/Fig2H_enterococcus_posthoc_power_analysis.txt"

sink(power_output_file)

cat("Figure 2H: Post-hoc power analysis for Enterococcus abundance\n")
cat("=============================================================\n\n")

cat("Input files:\n")
cat(seq_file, "\n")
cat(metadata_file, "\n\n")

cat("Sample numbers and Enterococcus relative abundance by group:\n")
print(summary_table)

cat("\nANOVA summary:\n")
print(aov_summary)

cat("\nEffect size and power summary:\n")
print(power_summary)

cat("\nInterpretation:\n")
cat("Cohen's f was calculated from the ANOVA sums of squares as:\n")
cat("partial eta squared = SS_between / (SS_between + SS_within)\n")
cat("Cohen's f = sqrt(partial eta squared / (1 - partial eta squared))\n")
cat("Observed power was estimated using pwr.anova.test with the harmonic mean sample size.\n")

sink()

message("Group summary saved to: outputs/Figure_2/Fig2H_enterococcus_posthoc_power_group_summary.csv")
message("Power analysis summary saved to: outputs/Figure_2/Fig2H_enterococcus_posthoc_power_analysis_summary.csv")
message("Detailed output saved to: ", power_output_file)