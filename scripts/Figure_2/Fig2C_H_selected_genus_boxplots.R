# Figure 2C-H: Boxplots of selected bacterial genera across disease severity groups
#
# Purpose:
# This script loads genus-level 16S rRNA sequencing data and clinical metadata,
# converts genus abundances to relative abundance percentages, and generates
# Figure 2C-H boxplots for selected genera associated with ACLD severity.
#
# Input files:
# - data/normalised_genus_data_16.04.2024.csv
# - data/clinical_metadata_and_data_renamed_30.01.2025.csv
# - outputs/Figure_2/Fig2B_maaslin2_genus_differential_abundance_results.csv
#
# Output files:
# - outputs/Figure_2/Fig2C_Ruminococcus_boxplot.pdf
# - outputs/Figure_2/Fig2D_Roseburia_boxplot.pdf
# - outputs/Figure_2/Fig2E_Lachnospira_boxplot.pdf
# - outputs/Figure_2/Fig2F_Streptococcus_boxplot.pdf
# - outputs/Figure_2/Fig2G_Staphylococcus_boxplot.pdf
# - outputs/Figure_2/Fig2H_Enterococcus_boxplot.pdf
# - outputs/Figure_2/Fig2C_H_selected_genus_boxplot_statistics.txt
# - outputs/Figure_2/Fig2C_H_selected_genus_maaslin2_results.csv
#
# Run this script from the root folder of the repository:
# liver_patient_enterococcus_analysis/

library(tidyverse)

options(scipen = 999)

# Create output folder if it does not already exist
dir.create("outputs/Figure_2", recursive = TRUE, showWarnings = FALSE)

# Define input files
seq_file <- "data/normalised_genus_data_16.04.2024.csv"
metadata_file <- "data/clinical_metadata_and_data_renamed_30.01.2025.csv"
maaslin2_results_file <- "outputs/Figure_2/Fig2B_maaslin2_genus_differential_abundance_results.csv"

# Check that required input files exist
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

# Selected genera and output mapping
selected_genera <- tibble(
  panel = c("Fig2C", "Fig2D", "Fig2E", "Fig2F", "Fig2G", "Fig2H"),
  genus = c(
    "Ruminococcus",
    "Roseburia",
    "Lachnospira",
    "Streptococcus",
    "Staphylococcus",
    "Enterococcus"
  ),
  y_limit = c(100, 40, 60, 59, 52, 115),
  output_file = c(
    "outputs/Figure_2/Fig2C_Ruminococcus_boxplot.pdf",
    "outputs/Figure_2/Fig2D_Roseburia_boxplot.pdf",
    "outputs/Figure_2/Fig2E_Lachnospira_boxplot.pdf",
    "outputs/Figure_2/Fig2F_Streptococcus_boxplot.pdf",
    "outputs/Figure_2/Fig2G_Staphylococcus_boxplot.pdf",
    "outputs/Figure_2/Fig2H_Enterococcus_boxplot.pdf"
  )
)

# Check that all selected genera are present
missing_genera <- setdiff(selected_genera$genus, names(merged_meta_seqs))

if (length(missing_genera) > 0) {
  stop(
    "The following selected genera are missing from the input data: ",
    paste(missing_genera, collapse = ", ")
  )
}

# Manual annotation positions and labels used in the original figure panels
annotation_table <- tribble(
  ~genus, ~text_x, ~text_y, ~label, ~segment_x, ~segment_xend, ~segment_y,
  "Enterococcus", 2.5, 112, "p < 0.001", 1, 4, 107,
  "Enterococcus", 2.0, 102, "p = 0.007", 1, 3, 97,
  
  "Ruminococcus", 2.5, 95, "p < 0.001", 1, 4, 91,
  "Ruminococcus", 2.0, 87, "p < 0.001", 1, 3, 83,
  "Ruminococcus", 1.5, 79, "p = 0.02", 1, 2, 75,
  
  "Roseburia", 2.5, 37, "p < 0.001", 1, 4, 35,
  "Roseburia", 2.0, 33, "p = 0.004", 1, 3, 31,
  
  "Lachnospira", 2.5, 58, "p < 0.001", 1, 4, 55,
  "Lachnospira", 2.0, 52, "p = 0.006", 1, 3, 49,
  
  "Staphylococcus", 2.5, 50, "p < 0.001", 1, 4, 47,
  
  "Streptococcus", 2.0, 56, "p = 0.039", 1, 3, 53,
  "Streptococcus", 1.5, 50, "p = 0.045", 1, 2, 47
)

# Shared plotting function
make_genus_boxplot <- function(data, genus_name, y_limit, annotations) {
  
  plot_object <- ggplot(
    data,
    aes(
      x = Inclusion_criteria,
      y = .data[[genus_name]],
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
      position = position_jitterdodge(jitter.width = 1),
      aes(
        colour = Inclusion_criteria,
        fill = Inclusion_criteria
      ),
      size = 1.5,
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
      limits = c(0, y_limit)
    ) +
    scale_x_discrete(
      labels = c("Healthy", "Stable", "Decomp", "ACLF")
    ) +
    labs(
      x = "Group",
      y = paste0(genus_name, " (%)")
    ) +
    theme(
      plot.title = element_text(size = 8),
      plot.margin = margin(t = 2, r = 2, b = 2, l = 2),
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
      axis.text.x = element_text(size = 8, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
      axis.title.y = element_text(size = 8, colour = "black"),
      axis.text.y = element_text(size = 8, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      aspect.ratio = 1.5
    )
  
  if (nrow(annotations) > 0) {
    plot_object <- plot_object +
      geom_text(
        data = annotations,
        aes(
          x = text_x,
          y = text_y,
          label = label
        ),
        inherit.aes = FALSE,
        size = 2.5
      ) +
      geom_segment(
        data = annotations,
        aes(
          x = segment_x,
          xend = segment_xend,
          y = segment_y,
          yend = segment_y
        ),
        inherit.aes = FALSE,
        linewidth = 0.35
      )
  }
  
  plot_object
}

# Generate and save each panel
for (i in seq_len(nrow(selected_genera))) {
  
  genus_name <- selected_genera$genus[i]
  output_file <- selected_genera$output_file[i]
  y_limit <- selected_genera$y_limit[i]
  
  annotations <- annotation_table %>%
    filter(genus == genus_name)
  
  genus_plot <- make_genus_boxplot(
    data = merged_meta_seqs,
    genus_name = genus_name,
    y_limit = y_limit,
    annotations = annotations
  )
  
  ggsave(
    filename = output_file,
    plot = genus_plot,
    width = 4,
    height = 6,
    units = "cm",
    dpi = 300
  )
}

# Load MaAsLin2 results if available
if (file.exists(maaslin2_results_file)) {
  
  maaslin2_results <- read_csv(maaslin2_results_file, show_col_types = FALSE)
  
  selected_maaslin2_results <- maaslin2_results %>%
    filter(feature %in% selected_genera$genus) %>%
    select(feature, value, coef, pval, qval) %>%
    arrange(feature, value)
  
  write_csv(
    selected_maaslin2_results,
    "outputs/Figure_2/Fig2C_H_selected_genus_maaslin2_results.csv"
  )
  
} else {
  
  selected_maaslin2_results <- tibble()
  
  warning(
    "MaAsLin2 results file was not found: ",
    maaslin2_results_file,
    "\nRun Fig2B_maaslin2_genus_differential_abundance.R first if MaAsLin2 results are required."
  )
}

# Save summary output
summary_output_file <- "outputs/Figure_2/Fig2C_H_selected_genus_boxplot_statistics.txt"

sink(summary_output_file)

cat("Figure 2C-H: Selected genus boxplots\n")
cat("====================================\n\n")

cat("Input files:\n")
cat(seq_file, "\n")
cat(metadata_file, "\n\n")

cat("Sample numbers by group:\n")
print(table(merged_meta_seqs$Inclusion_criteria, useNA = "ifany"))

cat("\nSelected genera:\n")
print(selected_genera)

cat("\nSummary statistics by group and genus:\n")
summary_statistics <- merged_meta_seqs %>%
  select(Inclusion_criteria, all_of(selected_genera$genus)) %>%
  pivot_longer(
    cols = all_of(selected_genera$genus),
    names_to = "Genus",
    values_to = "relative_abundance"
  ) %>%
  group_by(Genus, Inclusion_criteria) %>%
  summarise(
    n = n(),
    median = median(relative_abundance, na.rm = TRUE),
    IQR = IQR(relative_abundance, na.rm = TRUE),
    mean = mean(relative_abundance, na.rm = TRUE),
    sd = sd(relative_abundance, na.rm = TRUE),
    max = max(relative_abundance, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_statistics)

cat("\nMaAsLin2 results for selected genera:\n")

if (nrow(selected_maaslin2_results) > 0) {
  print(selected_maaslin2_results)
} else {
  cat("MaAsLin2 results file not available.\n")
}

sink()

message("Figure 2C-H boxplots saved to: outputs/Figure_2/")
message("Summary saved to: ", summary_output_file)