# Supplementary Figure 1B: PCoA plot of genus-level faecal microbiota composition
#
# Purpose:
# This script loads normalised genus-level 16S rRNA sequencing data and clinical
# metadata, calculates Bray-Curtis dissimilarities, runs principal coordinates
# analysis (PCoA), generates Supplementary Figure 1B, and saves associated
# statistical output.
#
# Input files:
# - data/normalised_genus_data_16.04.2024.csv
# - data/clinical_metadata_and_data_renamed_30.01.2025.csv
#
# Output files:
# - outputs/Supplementary_Figures/Supplementary_figure_1/SuppFig1B_genus_microbiota_PCoA_plot.pdf
# - outputs/Supplementary_Figures/Supplementary_figure_1/SuppFig1B_genus_microbiota_PCoA_plot.png
# - outputs/Supplementary_Figures/Supplementary_figure_1/SuppFig1B_genus_microbiota_PCoA_statistics.txt
#
# Run this script from the root folder of the repository:
# liver_patient_enterococcus_analysis/

library(tidyverse)
library(vegan)
library(ape)

options(scipen = 999)

# Create output folder if it does not already exist
dir.create(
  "outputs/Supplementary_Figures/Supplementary_figure_1",
  recursive = TRUE,
  showWarnings = FALSE
)

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

# Prepare genus abundance data
seq_data_genus <- seq_data %>%
  select(study_no, all_of(genus_cols)) %>%
  drop_na()

# Convert genus abundances to relative abundance percentage
seq_data_genus_percent <- seq_data_genus %>%
  mutate(
    row_total = rowSums(across(all_of(genus_cols)), na.rm = TRUE)
  )

if (any(seq_data_genus_percent$row_total == 0)) {
  stop("At least one sample has a total genus abundance of zero.")
}

seq_data_genus_percent <- seq_data_genus_percent %>%
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
    sample_date,
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

# Prepare metadata table
meta_table <- merged_meta_seqs %>%
  select(
    study_no,
    Inclusion_criteria
  ) %>%
  column_to_rownames("study_no")

# Prepare abundance table
abund_table <- merged_meta_seqs %>%
  select(study_no, all_of(genus_cols)) %>%
  column_to_rownames("study_no")

# Remove genera with zero abundance across all samples
abund_table <- abund_table[, colSums(abund_table != 0) > 0]

# Check that sample order matches between abundance and metadata tables
if (!identical(rownames(abund_table), rownames(meta_table))) {
  stop("Sample order does not match between abundance table and metadata table.")
}

# Calculate Bray-Curtis distance matrix
bray_dist <- vegdist(abund_table, method = "bray")

# Run PCoA
pcoa_result <- ape::pcoa(bray_dist)

# Percentage variance explained by PCoA axes
variance_explained <- 100 * pcoa_result$values$Relative_eig

axis_1_label <- paste0("PCoA1 (", round(variance_explained[1], 1), "%)")
axis_2_label <- paste0("PCoA2 (", round(variance_explained[2], 1), "%)")

# Prepare PCoA plot data
plot_data <- as.data.frame(pcoa_result$vectors[, 1:2]) %>%
  setNames(c("PCoA1", "PCoA2")) %>%
  rownames_to_column("study_no") %>%
  left_join(
    meta_table %>%
      rownames_to_column("study_no"),
    by = "study_no"
  ) %>%
  mutate(
    plot_group = recode(
      as.character(Inclusion_criteria),
      "Healthy control" = "Healthy",
      "Stable cirrhotic" = "Stable",
      "Decomp cirrhotic" = "Decomp",
      "ACLF" = "ACLF"
    ),
    plot_group = factor(
      plot_group,
      levels = c("Healthy", "Stable", "Decomp", "ACLF")
    )
  )

# PERMANOVA: microbiota composition by disease severity group
permanova_inclusion <- adonis2(
  bray_dist ~ Inclusion_criteria,
  data = meta_table,
  permutations = 9999
)

# Check homogeneity of multivariate dispersion
dispersion_model <- betadisper(bray_dist, meta_table$Inclusion_criteria)
dispersion_anova <- anova(dispersion_model)

# Generate Supplementary Figure 1B PCoA plot
suppfig1b <- ggplot(
  plot_data,
  aes(
    x = PCoA1,
    y = PCoA2
  )
) +
  geom_point(
    aes(
      colour = plot_group,
      shape = plot_group,
      fill = plot_group
    ),
    size = 2.5,
    alpha = 0.85
  ) +
  scale_fill_manual(
    values = c(
      "Healthy" = "lightskyblue",
      "Stable" = "yellow",
      "Decomp" = "orange",
      "ACLF" = "red"
    )
  ) +
  scale_color_manual(
    values = c(
      "Healthy" = "blue4",
      "Stable" = "yellow4",
      "Decomp" = "orange4",
      "ACLF" = "red4"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Healthy" = 24,
      "Stable" = 23,
      "Decomp" = 21,
      "ACLF" = 22
    )
  ) +
  labs(
    title = "Genus microbiota",
    x = axis_1_label,
    y = axis_2_label
  ) +
  theme(
    plot.title = element_text(size = 8),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.position = c(0.8, 0.13),
    legend.spacing.x = unit(0, "cm"),
    legend.spacing.y = unit(0, "cm"),
    legend.key.height = unit(0.5, "line"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.box.margin = margin(t = 0, r = 2, b = 0, l = 0),
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "grey"),
    legend.key = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 8, colour = "black"),
    axis.text.x = element_text(size = 8, colour = "black", angle = 0, hjust = 0.5),
    axis.title.y = element_text(size = 8, colour = "black"),
    axis.text.y = element_text(size = 8, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    aspect.ratio = 1
  )

# Save Supplementary Figure 1B
ggsave(
  filename = "outputs/Supplementary_Figures/Supplementary_figure_1/SuppFig1B_genus_microbiota_PCoA_plot.pdf",
  plot = suppfig1b,
  width = 6.25,
  height = 6.25,
  units = "cm",
  dpi = 300
)

# Save statistical output
stats_output_file <- "outputs/Supplementary_Figures/Supplementary_figure_1/SuppFig1B_genus_microbiota_PCoA_statistics.txt"

sink(stats_output_file)

cat("Supplementary Figure 1B: Genus-level microbiota PCoA\n")
cat("====================================================\n\n")

cat("Input files:\n")
cat(seq_file, "\n")
cat(metadata_file, "\n\n")

cat("Sample numbers by group:\n")
print(table(meta_table$Inclusion_criteria, useNA = "ifany"))

cat("\nNumber of genera included after zero-sum filtering:\n")
print(ncol(abund_table))

cat("\nPCoA variance explained:\n")
cat(axis_1_label, "\n")
cat(axis_2_label, "\n")

cat("\nPERMANOVA: Inclusion criteria\n")
print(permanova_inclusion)

cat("\nHomogeneity of multivariate dispersion: betadisper ANOVA\n")
print(dispersion_anova)

sink()

message("Figure saved to: y outputs/Figure_1/Supplementary_figure_1/SuppFig1B_genus_microbiota_PCoA_plot.pdf")
message("Statistics saved to: ", stats_output_file)

