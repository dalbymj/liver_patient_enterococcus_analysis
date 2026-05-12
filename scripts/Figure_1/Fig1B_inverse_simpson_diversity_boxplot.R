# Figure 1B: Inverse Simpson diversity across disease severity groups
#
# Purpose:
# This script loads genus-level 16S rRNA sequencing data and clinical metadata,
# calculates Inverse Simpson diversity for each participant, generates the
# Figure 1B boxplot, and runs statistical tests.
#
# Input files:
# - data/raw_genus_data_15.04.2024.csv
# - data/clinical_metadata_and_data_renamed_30.01.2025.csv
#
# Output files:
# - outputs/Figure_1/Fig1B_inverse_simpson_diversity_boxplot.pdf
# - outputs/Figure_1/Fig1B_inverse_simpson_diversity_statistics.txt
#
# Run this script from the root folder of the repository:
# liver_patient_enterococcus_analysis/

library(tidyverse)
library(vegan)
library(car)
library(FSA)

options(scipen = 999)

# Create output folder if it does not already exist
dir.create("outputs/Figure_1", recursive = TRUE, showWarnings = FALSE)

# Define input files
seq_file <- "data/raw_genus_data_15.04.2024.csv"
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

# Define genus columns.
# This assumes that the first genus column is called 'Mucispirillum'
# and that all subsequent columns are genus abundance columns.
genus_start_col <- match("Mucispirillum", names(seq_data))

if (is.na(genus_start_col)) {
  stop("Could not find genus column 'Mucispirillum' in seq_data.")
}

genus_cols <- names(seq_data)[genus_start_col:ncol(seq_data)]

# Prepare sequence data
seq_data_genus <- seq_data %>%
  select(study_no, all_of(genus_cols))

# Prepare metadata
metadata_subsampled <- metadata %>%
  select(
    study_no,
    sample_date,
    Inclusion_criteria
  )

# Merge metadata and genus abundance data
merged_meta_seqs <- metadata_subsampled %>%
  inner_join(seq_data_genus, by = "study_no")

# Prepare metadata for plotting
metadata_for_plot <- merged_meta_seqs %>%
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

# Prepare abundance matrix
abund <- merged_meta_seqs %>%
  select(study_no, all_of(genus_cols)) %>%
  column_to_rownames("study_no")

# Remove genera with zero abundance across all samples
abund <- abund[, colSums(abund != 0) > 0]

# Check that all samples have non-zero total abundance
if (any(rowSums(abund) == 0)) {
  stop("At least one sample has a total genus abundance of zero.")
}

# Convert to relative abundance percentage.
# Inverse Simpson diversity is calculated from relative abundances.
abund_percent <- abund / rowSums(abund) * 100

# Calculate diversity measures
diversity_measures <- tibble(
  study_no = rownames(abund_percent),
  inverse_simpson = diversity(abund_percent, index = "invsimpson")
)

# Merge diversity measures with metadata
meta_diversity_measures <- metadata_for_plot %>%
  left_join(diversity_measures, by = "study_no")

# Generate Figure 1B
fig1b <- ggplot(
  meta_diversity_measures,
  aes(
    x = Inclusion_criteria,
    y = inverse_simpson,
    fill = Inclusion_criteria,
    shape = Inclusion_criteria
  )
) +
  geom_boxplot(
    position = position_dodge(),
    outlier.shape = NA,
    aes(
      colour = Inclusion_criteria,
      fill = Inclusion_criteria
    ),
    linewidth = 0.35
  ) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.8),
    aes(
      colour = Inclusion_criteria,
      fill = Inclusion_criteria
    ),
    size = 2,
    stroke = 0.35
  ) +
  scale_fill_manual(
    values = c(
      "lightskyblue",
      "yellow",
      "orange",
      "firebrick1"
    )
  ) +
  scale_color_manual(
    values = c(
      "blue4",
      "yellow4",
      "orange4",
      "red4"
    )
  ) +
  scale_shape_manual(values = c(1, 1, 1, 1)) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 12)
  ) +
  scale_x_discrete(
    labels = c(
      "Healthy\ncontrol",
      "Stable\ncirrhotic",
      "Decomp\ncirrhotic",
      "ACLF"
    )
  ) +
  labs(
    x = "Group",
    y = "Inverse Simpson diversity"
  ) +
  theme(
    plot.title = element_text(size = 8),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.position = "none",
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "grey"),
    legend.key = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 8, colour = "black"),
    axis.text.x = element_text(size = 6, colour = "black", angle = 0, hjust = 0.5),
    axis.title.y = element_text(size = 8, colour = "black"),
    axis.text.y = element_text(size = 8, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    aspect.ratio = 1.3
  ) +
  annotate("text", x = 2.5, y = 11.2, label = "p = 0.00011", size = 2.5) +
  annotate("segment", x = 1, xend = 4, y = 10.8, yend = 10.8, linewidth = 0.35) +
  annotate("text", x = 3, y = 10.4, label = "p = 0.037", size = 2.5) +
  annotate("segment", x = 2, xend = 4, y = 10, yend = 10, linewidth = 0.35) +
  annotate("text", x = 3.5, y = 9.6, label = "p = 0.04", size = 2) +
  annotate("segment", x = 3, xend = 4, y = 9.2, yend = 9.2, linewidth = 0.35)

# Save Figure 1B
ggsave(
  filename = "outputs/Figure_1/Fig1B_inverse_simpson_diversity_boxplot.pdf",
  plot = fig1b,
  width = 5,
  height = 6,
  units = "cm",
  dpi = 300
)

# Statistical testing
stats_output_file <- "outputs/Figure_1/Fig1B_inverse_simpson_diversity_statistics.txt"

sink(stats_output_file)

cat("Figure 1B: Inverse Simpson diversity\n")
cat("====================================\n\n")

cat("Sample numbers by group:\n")
print(table(meta_diversity_measures$Inclusion_criteria, useNA = "ifany"))

cat("\nSummary statistics by group:\n")
print(
  meta_diversity_measures %>%
    group_by(Inclusion_criteria) %>%
    summarise(
      n = n(),
      median = median(inverse_simpson, na.rm = TRUE),
      IQR = IQR(inverse_simpson, na.rm = TRUE),
      mean = mean(inverse_simpson, na.rm = TRUE),
      sd = sd(inverse_simpson, na.rm = TRUE),
      .groups = "drop"
    )
)

cat("\nShapiro-Wilk tests by group:\n")
print(
  meta_diversity_measures %>%
    group_by(Inclusion_criteria) %>%
    summarise(
      shapiro_p = shapiro.test(inverse_simpson)$p.value,
      .groups = "drop"
    )
)

cat("\nANOVA:\n")
res_aov <- aov(inverse_simpson ~ Inclusion_criteria, data = meta_diversity_measures)
print(summary(res_aov))

cat("\nTukey post-hoc test:\n")
print(TukeyHSD(res_aov))

cat("\nLevene's test for homogeneity of variance:\n")
print(leveneTest(inverse_simpson ~ Inclusion_criteria, data = meta_diversity_measures))

cat("\nShapiro-Wilk test of ANOVA residuals:\n")
aov_residuals <- residuals(object = res_aov)
print(shapiro.test(aov_residuals))

cat("\nKruskal-Wallis test:\n")
print(kruskal.test(inverse_simpson ~ Inclusion_criteria, data = meta_diversity_measures))

cat("\nDunn test with Benjamini-Hochberg correction:\n")
print(dunnTest(inverse_simpson ~ Inclusion_criteria, data = meta_diversity_measures, method = "bh"))

sink()

message("Figure saved to: outputs/Figure_1/Fig1B_inverse_simpson_diversity_boxplot.pdf")
message("Statistics saved to: ", stats_output_file)