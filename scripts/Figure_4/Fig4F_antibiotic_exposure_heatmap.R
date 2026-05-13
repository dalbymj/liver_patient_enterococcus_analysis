# Figure 4F: Antibiotic exposure heatmap for liver disease samples
#
# Purpose:
# This script loads genus-level 16S rRNA sequencing data and clinical metadata,
# converts genus abundances to relative abundance percentages, extracts antibiotic
# exposure variables for liver disease samples, clusters samples and antibiotic
# variables using binary Jaccard distance, and generates the Figure 4F
# ComplexHeatmap with Enterococcus abundance and clinical annotations.
#
# Input files:
# - data/normalised_genus_data_16.04.2024.csv
# - data/clinical_metadata_and_data_renamed_30.01.2025.csv
#
# Output files:
# - outputs/Figure_4/Fig4F_antibiotic_exposure_heatmap.pdf
# - outputs/Figure_4/Fig4F_antibiotic_exposure_matrix.csv
# - outputs/Figure_4/Fig4F_antibiotic_sample_annotations.csv
# - outputs/Figure_4/Fig4F_antibiotic_totals.csv
# - outputs/Figure_4/Fig4F_antibiotic_exposure_heatmap_summary.txt
#
# Run this script from the root folder of the repository:
# liver_patient_enterococcus_analysis/

library(tidyverse)
library(vegan)
library(ComplexHeatmap)
library(grid)

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

# Define genus columns
genus_start_col <- match("Mucispirillum", names(genus_data_raw))

if (is.na(genus_start_col)) {
  stop("Could not find genus column 'Mucispirillum' in genus_data_raw.")
}

genus_cols <- names(genus_data_raw)[genus_start_col:ncol(genus_data_raw)]

if (!"Enterococcus" %in% genus_cols) {
  stop("Could not find genus column 'Enterococcus' in genus_data_raw.")
}

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

# Define antibiotic columns
antibiotic_start_col <- match("Tazocin", names(metadata_raw))
antibiotic_end_col <- match("Trimethoprim_primazole", names(metadata_raw))

if (is.na(antibiotic_start_col)) {
  stop("Could not find antibiotic column 'Tazocin' in metadata_raw.")
}

if (is.na(antibiotic_end_col)) {
  stop("Could not find antibiotic column 'Trimethoprim_primazole' in metadata_raw.")
}

antibiotic_cols <- names(metadata_raw)[antibiotic_start_col:antibiotic_end_col]

# Prepare metadata
metadata_subsampled <- metadata_raw %>%
  select(
    study_no,
    sample_date,
    Inclusion_criteria,
    Oral_intake,
    Receiving_antibiotics,
    all_of(antibiotic_cols)
  ) %>%
  rename(
    Primazole = Trimethoprim_primazole
  ) %>%
  mutate(
    Inclusion_criteria = as.character(Inclusion_criteria),
    Inclusion_criteria = dplyr::recode(
      Inclusion_criteria,
      "AoCLF" = "ACLF",
      .default = Inclusion_criteria
    ),
    
    Receiving_antibiotics_raw = as.character(Receiving_antibiotics),
    
    Receiving_antibiotics = case_when(
      Receiving_antibiotics_raw == "Healthy control" ~ "Healthy control",
      Receiving_antibiotics_raw %in% c("0", "No", "No_antibiotics") ~ "No",
      Receiving_antibiotics_raw %in% c("1", "Yes", "Receiving_antibiotics", "Reciving_antibiotics") ~ "Yes",
      TRUE ~ Receiving_antibiotics_raw
    ),
    
    No_antibiotics = case_when(
      Receiving_antibiotics == "No" ~ 1,
      Receiving_antibiotics == "Yes" ~ 0,
      Receiving_antibiotics == "Healthy control" ~ NA_real_,
      TRUE ~ NA_real_
    ),
    
    Antibiotics = Receiving_antibiotics,
    
    Oral_intake = as.character(Oral_intake),
    Food_intake = case_when(
      Oral_intake == "Healthy control" ~ "Healthy control",
      Oral_intake %in% c("0", "Yes", "Oral_food_intake") ~ "Yes",
      Oral_intake %in% c("1", "No", "No_oral_food_intake") ~ "No",
      TRUE ~ Oral_intake
    )
  )

# Update antibiotic column names after renaming Trimethoprim_primazole
antibiotic_cols_renamed <- antibiotic_cols
antibiotic_cols_renamed[antibiotic_cols_renamed == "Trimethoprim_primazole"] <- "Primazole"

# Convert individual antibiotic columns to numeric
metadata_subsampled <- metadata_subsampled %>%
  mutate(
    across(
      all_of(antibiotic_cols_renamed),
      ~ as.numeric(.x)
    )
  )

# Merge metadata and genus abundance data
merged_meta_seqs <- metadata_subsampled %>%
  inner_join(genus_data_percent, by = "study_no")

# Keep liver disease samples only
selected_data <- merged_meta_seqs %>%
  filter(!Inclusion_criteria %in% c("Healthy control")) %>%
  select(
    study_no,
    sample_date,
    Inclusion_criteria,
    Food_intake,
    No_antibiotics,
    Receiving_antibiotics,
    Antibiotics,
    all_of(antibiotic_cols_renamed),
    Enterococcus,
    Roseburia
  )

# Prepare antibiotic exposure matrix
antibiotic_matrix_data <- selected_data %>%
  select(
    study_no,
    No_antibiotics,
    all_of(antibiotic_cols_renamed)
  )

# Remove antibiotic columns with zero exposure across all retained samples
nonzero_antibiotic_cols <- antibiotic_matrix_data %>%
  select(-study_no) %>%
  select(where(~ sum(.x != 0, na.rm = TRUE) > 0)) %>%
  names()

if (length(nonzero_antibiotic_cols) == 0) {
  stop("No non-zero antibiotic exposure columns remain after filtering.")
}

antibiotic_matrix_data <- antibiotic_matrix_data %>%
  select(study_no, all_of(nonzero_antibiotic_cols))

antibiotic_matrix <- antibiotic_matrix_data %>%
  column_to_rownames("study_no") %>%
  mutate(
    across(
      everything(),
      ~ as.numeric(.x)
    )
  ) %>%
  as.matrix()

# Replace any remaining missing values in the antibiotic matrix with 0
antibiotic_matrix[is.na(antibiotic_matrix)] <- 0

# Save matrix used for the heatmap
write_csv(
  antibiotic_matrix %>%
    as.data.frame() %>%
    rownames_to_column("study_no"),
  "outputs/Figure_4/Fig4F_antibiotic_exposure_matrix.csv"
)

# Prepare sample annotation data
sample_annotations <- selected_data %>%
  select(
    study_no,
    Inclusion_criteria,
    Antibiotics,
    Food_intake,
    Enterococcus
  ) %>%
  filter(study_no %in% rownames(antibiotic_matrix)) %>%
  arrange(match(study_no, rownames(antibiotic_matrix))) %>%
  mutate(
    Group = dplyr::recode(
      Inclusion_criteria,
      "Stable cirrhotic" = "Stable",
      "Decomp cirrhotic" = "Decomp",
      "ACLF" = "ACLF",
      "AoCLF" = "ACLF",
      .default = Inclusion_criteria
    ),
    Group = factor(
      Group,
      levels = c("Stable", "Decomp", "ACLF")
    ),
    Antibiotics = factor(
      Antibiotics,
      levels = c("No", "Yes")
    ),
    Food_intake = factor(
      Food_intake,
      levels = c("Yes", "No")
    )
  )

if (!identical(sample_annotations$study_no, rownames(antibiotic_matrix))) {
  stop("Sample order does not match between antibiotic matrix and sample annotations.")
}

write_csv(
  sample_annotations,
  "outputs/Figure_4/Fig4F_antibiotic_sample_annotations.csv"
)

# Save antibiotic totals
antibiotic_totals <- colSums(antibiotic_matrix, na.rm = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("antibiotic") %>%
  rename(total_samples = ".")

write_csv(
  antibiotic_totals,
  "outputs/Figure_4/Fig4F_antibiotic_totals.csv"
)

# Calculate binary Jaccard distance and hierarchical clustering
row_distance <- vegdist(
  antibiotic_matrix,
  method = "jaccard",
  binary = TRUE
)

row_clustering <- hclust(
  row_distance,
  method = "average"
)

column_distance <- vegdist(
  t(antibiotic_matrix),
  method = "jaccard",
  binary = TRUE
)

column_clustering <- hclust(
  column_distance,
  method = "average"
)

# ComplexHeatmap annotations
group_annotation <- data.frame(
  Group = sample_annotations$Group
)

antibiotics_annotation <- data.frame(
  Antibiotics = sample_annotations$Antibiotics
)

food_intake_annotation <- data.frame(
  Food_intake = sample_annotations$Food_intake
)

enterococcus_annotation <- data.frame(
  Enterococcus = sample_annotations$Enterococcus
) %>%
  as.matrix()

group_colours <- list(
  Group = c(
    "Stable" = "yellow",
    "Decomp" = "orange",
    "ACLF" = "red3"
  )
)

antibiotic_colours <- list(
  Antibiotics = c(
    "No" = "grey80",
    "Yes" = "grey20"
  )
)

food_intake_colours <- list(
  Food_intake = c(
    "Yes" = "lightskyblue",
    "No" = "red"
  )
)

# Define heatmap colour palette
heatmap_colour_palette <- colorRampPalette(
  colors = c("grey90", "black")
)(100)

# Row annotations
group_sidebar_annotation <- rowAnnotation(
  df = group_annotation,
  col = group_colours,
  simple_anno_size = unit(0.25, "cm"),
  annotation_legend_param = list(
    title = "Group",
    title_gp = gpar(fontsize = 7),
    labels_gp = gpar(fontsize = 7)
  ),
  show_annotation_name = c("Group" = TRUE),
  annotation_name_gp = gpar(fontsize = 7)
)

antibiotics_sidebar_annotation <- rowAnnotation(
  df = antibiotics_annotation,
  col = antibiotic_colours,
  simple_anno_size = unit(0.25, "cm"),
  annotation_legend_param = list(
    title = "Antibiotics",
    title_gp = gpar(fontsize = 7),
    labels_gp = gpar(fontsize = 7)
  ),
  show_annotation_name = c("Antibiotics" = TRUE),
  annotation_name_gp = gpar(fontsize = 7)
)

food_intake_sidebar_annotation <- rowAnnotation(
  df = food_intake_annotation,
  col = food_intake_colours,
  simple_anno_size = unit(0.25, "cm"),
  annotation_legend_param = list(
    title = "Food intake",
    title_gp = gpar(fontsize = 7),
    labels_gp = gpar(fontsize = 7)
  ),
  show_annotation_name = c("Food_intake" = TRUE),
  annotation_name_gp = gpar(fontsize = 7)
)

enterococcus_barplot_annotation <- rowAnnotation(
  "Enterococcus" = anno_barplot(
    enterococcus_annotation,
    width = unit(0.8, "cm"),
    which = "row",
    axis = TRUE,
    bar_width = 1,
    show_annotation_name = TRUE,
    axis_param = list(
      gp = gpar(fontsize = 7),
      at = c(0, 50, 100),
      labels_rot = 90
    )
  )
)

# Create heatmap
antibiotic_heatmap <- Heatmap(
  as.matrix(antibiotic_matrix),
  rect_gp = gpar(col = "white", lwd = 1),
  name = "Exposure",
  col = heatmap_colour_palette,
  cluster_rows = row_clustering,
  cluster_columns = column_clustering,
  show_row_names = FALSE,
  row_names_gp = gpar(fontsize = 3),
  column_names_gp = gpar(fontsize = 7),
  column_title_gp = gpar(fontsize = 7),
  row_title = "Samples",
  row_title_gp = gpar(fontsize = 7),
  column_title = "Antibiotic",
  column_title_side = "bottom",
  show_heatmap_legend = FALSE
)

# Combine heatmap and annotations
fig4f <- antibiotic_heatmap +
  enterococcus_barplot_annotation +
  group_sidebar_annotation +
  antibiotics_sidebar_annotation

# Add food_intake_sidebar_annotation here if needed:
# fig4f <- antibiotic_heatmap +
#   enterococcus_barplot_annotation +
#   group_sidebar_annotation +
#   food_intake_sidebar_annotation +
#   antibiotics_sidebar_annotation

# Save heatmap
pdf(
  file = "outputs/Figure_4/Fig4F_antibiotic_exposure_heatmap.pdf",
  width = 4,
  height = 4
)

print(fig4f)

dev.off()

# Save summary output
summary_output_file <- "outputs/Figure_4/Fig4F_antibiotic_exposure_heatmap_summary.txt"

sink(summary_output_file)

cat("Figure 4F: Antibiotic exposure heatmap for liver disease samples\n")
cat("================================================================\n\n")

cat("Input files:\n")
cat(genus_file, "\n")
cat(metadata_file, "\n\n")

cat("Number of liver disease samples included:\n")
print(nrow(antibiotic_matrix))

cat("\nAntibiotic exposure columns included after zero-sum filtering:\n")
print(colnames(antibiotic_matrix))

cat("\nAntibiotic totals:\n")
print(antibiotic_totals)

cat("\nSample numbers by disease group:\n")
print(table(sample_annotations$Group, useNA = "ifany"))

cat("\nSample numbers by antibiotic exposure:\n")
print(table(sample_annotations$Antibiotics, useNA = "ifany"))

cat("\nSample numbers by food intake:\n")
print(table(sample_annotations$Food_intake, useNA = "ifany"))

cat("\nEnterococcus relative abundance summary:\n")
print(summary(sample_annotations$Enterococcus))

cat("\nClustering methods:\n")
cat("Rows: binary Jaccard distance with average-linkage hierarchical clustering.\n")
cat("Columns: binary Jaccard distance with average-linkage hierarchical clustering.\n")

sink()

message("Figure saved to: outputs/Figure_4/Fig4F_antibiotic_exposure_heatmap.pdf")
message("Antibiotic matrix saved to: outputs/Figure_4/Fig4F_antibiotic_exposure_matrix.csv")
message("Sample annotations saved to: outputs/Figure_4/Fig4F_antibiotic_sample_annotations.csv")
message("Antibiotic totals saved to: outputs/Figure_4/Fig4F_antibiotic_totals.csv")
message("Summary saved to: ", summary_output_file)