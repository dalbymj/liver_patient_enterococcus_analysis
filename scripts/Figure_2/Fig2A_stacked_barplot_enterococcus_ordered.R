# Figure 2A: Stacked bar plot of dominant faecal bacterial genera
#
# Purpose:
# This script loads genus-level 16S rRNA sequencing data and clinical metadata,
# converts genus counts/abundances to relative abundance percentages, selects
# the dominant genera shown in Figure 2A, groups remaining genera as "other",
# orders samples within disease groups by Enterococcus relative abundance, and
# generates the Figure 2A stacked bar plot.
#
# Input files:
# - data/normalised_genus_data_16.04.2024.csv
# - data/clinical_metadata_and_data_renamed_30.01.2025.csv
#
# Output files:
# - outputs/Figure_2/Fig2A_stacked_barplot_enterococcus_ordered.pdf
# - outputs/Figure_2/Fig2A_stacked_barplot_enterococcus_ordered_data.csv
# - outputs/Figure_2/Fig2A_stacked_barplot_enterococcus_ordered_summary.txt
#
# Run this script from the root folder of the repository:
# liver_patient_enterococcus_analysis/

library(tidyverse)

options(scipen = 999)

# Create output folder if it does not already exist
dir.create("outputs/Figure_2", recursive = TRUE, showWarnings = FALSE)

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

# Genera displayed individually in Figure 2A
selected_genera <- c(
  "Enterococcus",
  "Bacteroides",
  "Blautia",
  "Bifidobacterium",
  "Lactobacillus",
  "Streptococcus",
  "Ruminococcus",
  "Roseburia",
  "Citrobacter",
  "Prevotella"
)

# Check that all selected genera are present
missing_genera <- setdiff(selected_genera, names(merged_meta_seqs))

if (length(missing_genera) > 0) {
  stop(
    "The following selected genera are missing from the input data: ",
    paste(missing_genera, collapse = ", ")
  )
}

# Identify genera to group as "other"
other_genera <- setdiff(genus_cols, selected_genera)

# Prepare wide-format plotting data
plot_data_wide <- merged_meta_seqs %>%
  mutate(
    other = rowSums(across(all_of(other_genera)), na.rm = TRUE)
  ) %>%
  select(
    study_no,
    Inclusion_criteria,
    all_of(selected_genera),
    other
  ) %>%
  arrange(
    Inclusion_criteria,
    Enterococcus
  ) %>%
  mutate(
    study_no = factor(
      study_no,
      levels = unique(study_no)
    )
  )

# Check that relative abundance sums to approximately 100%
plot_data_wide <- plot_data_wide %>%
  mutate(
    relative_abundance_sum = rowSums(
      across(c(all_of(selected_genera), other)),
      na.rm = TRUE
    )
  )

if (any(abs(plot_data_wide$relative_abundance_sum - 100) > 0.01)) {
  warning("At least one sample does not sum to approximately 100%. Check input data.")
}

# Convert to long format for plotting
plot_data_long <- plot_data_wide %>%
  select(-relative_abundance_sum) %>%
  pivot_longer(
    cols = c(all_of(selected_genera), other),
    names_to = "Genus",
    values_to = "relative_abundance"
  ) %>%
  mutate(
    plot_group = dplyr::recode(
      as.character(Inclusion_criteria),
      "Healthy control" = "Healthy",
      "Stable cirrhotic" = "Stable",
      "Decomp cirrhotic" = "Decomp",
      "ACLF" = "ACLF"
    ),
    plot_group = factor(
      plot_group,
      levels = c("Healthy", "Stable", "Decomp", "ACLF")
    ),
    Genus = factor(
      Genus,
      levels = c(selected_genera, "other")
    )
  )

# Colour palette for stacked bar plot
genus_palette <- c(
  "Enterococcus" = "red",
  "Bacteroides" = "steelblue1",
  "Blautia" = "black",
  "Bifidobacterium" = "cadetblue1",
  "Lactobacillus" = "darkolivegreen1",
  "Streptococcus" = "pink",
  "Ruminococcus" = "purple",
  "Roseburia" = "darkgreen",
  "Citrobacter" = "red4",
  "Prevotella" = "orange",
  "other" = "grey"
)

# Generate Figure 2A
fig2a <- ggplot(
  plot_data_long,
  aes(
    x = study_no,
    y = relative_abundance,
    fill = Genus
  )
) +
  geom_bar(
    stat = "identity",
    width = 0.6
  ) +
  facet_grid(
    . ~ plot_group,
    drop = TRUE,
    scales = "free_x",
    space = "free"
  ) +
  scale_fill_manual(
    values = genus_palette
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 100.1)
  ) +
  labs(
    title = "",
    x = "Samples",
    y = "Relative abundance (%)"
  ) +
  theme(
    legend.key.size = unit(0.3, "cm"),
    panel.spacing = unit(0.1, "cm"),
    legend.position = "bottom",
    plot.title = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 8, face = "italic"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.spacing.x = unit(0.1, "cm"),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_text(size = 8),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8, colour = "black"),
    axis.line = element_line(linewidth = 0.35, linetype = "solid", colour = "black"),
    strip.text.x = element_text(size = 6, colour = "black", angle = 0),
    strip.background = element_rect(colour = "grey70", fill = "grey90")
  ) +
  guides(
    fill = guide_legend(nrow = 2)
  )

# Save Figure 2A
ggsave(
  filename = "outputs/Figure_2/Fig2A_stacked_barplot_enterococcus_ordered.pdf",
  plot = fig2a,
  width = 16.8,
  height = 6,
  units = "cm",
  dpi = 300
)

# Save processed plotting data
write_csv(
  plot_data_long,
  "outputs/Figure_2/Fig2A_stacked_barplot_enterococcus_ordered_data.csv"
)

# Save summary output
summary_output_file <- "outputs/Figure_2/Fig2A_stacked_barplot_enterococcus_ordered_summary.txt"

sink(summary_output_file)

cat("Figure 2A: Stacked bar plot of dominant faecal bacterial genera\n")
cat("================================================================\n\n")

cat("Input files:\n")
cat(seq_file, "\n")
cat(metadata_file, "\n\n")

cat("Sample numbers by group:\n")
print(table(plot_data_wide$Inclusion_criteria, useNA = "ifany"))

cat("\nGenera displayed individually:\n")
print(selected_genera)

cat("\nNumber of genera grouped as 'other':\n")
print(length(other_genera))

cat("\nRelative abundance sum check:\n")
print(summary(plot_data_wide$relative_abundance_sum))

cat("\nEnterococcus relative abundance by group:\n")
print(
  plot_data_wide %>%
    group_by(Inclusion_criteria) %>%
    summarise(
      n = n(),
      median = median(Enterococcus, na.rm = TRUE),
      IQR = IQR(Enterococcus, na.rm = TRUE),
      mean = mean(Enterococcus, na.rm = TRUE),
      sd = sd(Enterococcus, na.rm = TRUE),
      max = max(Enterococcus, na.rm = TRUE),
      .groups = "drop"
    )
)

sink()

message("Figure saved to: outputs/Figure_2/Fig2A_stacked_barplot_enterococcus_ordered.pdf")
message("Processed plotting data saved to: outputs/Figure_2/Fig2A_stacked_barplot_enterococcus_ordered_data.csv")
message("Summary saved to: ", summary_output_file)