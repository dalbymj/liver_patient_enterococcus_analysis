# Figure 3B-C: Bile acid boxplots across disease severity groups
#
# Purpose:
# This script loads plasma bile acid data and clinical metadata, then generates
# Figure 3B-C bile acid boxplots for selected bile acids associated with ACLD
# disease severity.
#
# Input files:
# - data/Corrected bile acid data 01.12.2020.csv
# - data/clinical_metadata_and_data_renamed_30.01.2025.csv
# - outputs/Figure_3/Fig3A_maaslin2_bile_acids_filtered_results.csv
#
# Output files:
# - outputs/Figure_3/Fig3B_Taurochenodeoxycholic_acid_boxplot.pdf
# - outputs/Figure_3/Fig3C_Glycochenodeoxycholic_acid_boxplot.pdf
# - outputs/Figure_3/Fig3B_C_bile_acid_boxplot_statistics.txt
# - outputs/Figure_3/Fig3B_C_selected_bile_acid_maaslin2_results.csv
#
# Run this script from the root folder of the repository:
# liver_patient_enterococcus_analysis/

library(tidyverse)

options(scipen = 999)

# Create output folder if it does not already exist
dir.create("outputs/Figure_3", recursive = TRUE, showWarnings = FALSE)

# Define input files
bile_file <- "data/Corrected bile acid data 01.12.2020.csv"
metadata_file <- "data/clinical_metadata_and_data_renamed_30.01.2025.csv"
maaslin2_results_file <- "outputs/Figure_3/Fig3A_maaslin2_bile_acids_filtered_results.csv"

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

# Define bile acid columns
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

# Bile acids shown in main Figure 3B-C
selected_bile_acids <- tribble(
  ~panel, ~bile_acid, ~display_name, ~y_label, ~y_limit, ~width_cm, ~height_cm, ~output_file,
  "Fig3B", "Taurochenodeoxycholic_acid", "Taurochenodeoxycholic acid", "Taurochenodeoxycholic acid (mmol/L)", 215, 8, 8, "outputs/Figure_3/Fig3B_Taurochenodeoxycholic_acid_boxplot.pdf",
  "Fig3C", "Glycochenodeoxycholic_acid", "Glycochenodeoxycholic acid", "Glycochenodeoxycholic acid (mmol/L)", 135, 8, 8, "outputs/Figure_3/Fig3C_Glycochenodeoxycholic_acid_boxplot.pdf"
)

# Check that all selected bile acids are present
missing_bile_acids <- setdiff(selected_bile_acids$bile_acid, names(merged_meta_bile_acids))

if (length(missing_bile_acids) > 0) {
  stop(
    "The following selected bile acids are missing from the input data: ",
    paste(missing_bile_acids, collapse = ", ")
  )
}

# Manual annotation positions and labels used in Figure 3B-C
annotation_table <- tribble(
  ~bile_acid, ~text_x, ~text_y, ~label, ~segment_x, ~segment_xend, ~segment_y,
  
  "Taurochenodeoxycholic_acid", 2.5, 209, "p < 0.001", 1, 4, 201,
  "Taurochenodeoxycholic_acid", 2.0, 193, "p = 0.02", 1, 3, 185,
  "Taurochenodeoxycholic_acid", 1.5, 176, "p = 0.01", 1, 2, 168,
  
  "Glycochenodeoxycholic_acid", 2.5, 128, "p < 0.001", 1, 4, 122,
  "Glycochenodeoxycholic_acid", 2.0, 116, "p = 0.02", 1, 3, 110,
  "Glycochenodeoxycholic_acid", 1.5, 104, "p = 0.01", 1, 2, 98
)

# Shared plotting function
make_bile_acid_boxplot <- function(data, bile_acid_name, y_label, y_limit, annotations) {
  
  plot_object <- ggplot(
    data,
    aes(
      x = Inclusion_criteria,
      y = .data[[bile_acid_name]],
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
      labels = c(
        "Healthy\ncontrol",
        "Stable\ncirrhotic",
        "Decomp\ncirrhotic",
        "ACLF"
      )
    ) +
    labs(
      x = "Group",
      y = y_label
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
      axis.text.x = element_text(size = 7, colour = "black", angle = 90, hjust = 0, vjust = 0.5),
      axis.title.y = element_text(size = 8, colour = "black"),
      axis.text.y = element_text(size = 7, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      aspect.ratio = 2
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

# Generate and save each bile acid panel
for (i in seq_len(nrow(selected_bile_acids))) {
  
  bile_acid_name <- selected_bile_acids$bile_acid[i]
  y_label <- selected_bile_acids$y_label[i]
  y_limit <- selected_bile_acids$y_limit[i]
  output_file <- selected_bile_acids$output_file[i]
  width_cm <- selected_bile_acids$width_cm[i]
  height_cm <- selected_bile_acids$height_cm[i]
  
  annotations <- annotation_table %>%
    filter(bile_acid == bile_acid_name)
  
  bile_acid_plot <- make_bile_acid_boxplot(
    data = merged_meta_bile_acids,
    bile_acid_name = bile_acid_name,
    y_label = y_label,
    y_limit = y_limit,
    annotations = annotations
  )
  
  ggsave(
    filename = output_file,
    plot = bile_acid_plot,
    width = width_cm,
    height = height_cm,
    units = "cm",
    dpi = 300
  )
}

# Load MaAsLin2 filtered results if available
if (file.exists(maaslin2_results_file)) {
  
  maaslin2_results <- read_csv(maaslin2_results_file, show_col_types = FALSE)
  
  selected_maaslin2_results <- maaslin2_results %>%
    filter(feature %in% selected_bile_acids$bile_acid) %>%
    select(feature, value, coef, pval, qval) %>%
    arrange(feature, value)
  
  write_csv(
    selected_maaslin2_results,
    "outputs/Figure_3/Fig3B_C_selected_bile_acid_maaslin2_results.csv"
  )
  
} else {
  
  selected_maaslin2_results <- tibble()
  
  warning(
    "MaAsLin2 results file was not found: ",
    maaslin2_results_file,
    "\nRun Fig3A_maaslin2_bile_acids_differential_abundance.R first if MaAsLin2 results are required."
  )
}

# Save summary output
summary_output_file <- "outputs/Figure_3/Fig3B_C_bile_acid_boxplot_statistics.txt"

sink(summary_output_file)

cat("Figure 3B-C: Bile acid boxplots\n")
cat("================================\n\n")

cat("Input files:\n")
cat(bile_file, "\n")
cat(metadata_file, "\n\n")

cat("Sample numbers by group:\n")
print(table(merged_meta_bile_acids$Inclusion_criteria, useNA = "ifany"))

cat("\nSelected bile acids:\n")
print(selected_bile_acids)

cat("\nSummary statistics by group and bile acid:\n")
summary_statistics <- merged_meta_bile_acids %>%
  select(Inclusion_criteria, all_of(selected_bile_acids$bile_acid)) %>%
  pivot_longer(
    cols = all_of(selected_bile_acids$bile_acid),
    names_to = "bile_acid",
    values_to = "concentration"
  ) %>%
  group_by(bile_acid, Inclusion_criteria) %>%
  summarise(
    n = n(),
    median = median(concentration, na.rm = TRUE),
    IQR = IQR(concentration, na.rm = TRUE),
    mean = mean(concentration, na.rm = TRUE),
    sd = sd(concentration, na.rm = TRUE),
    max = max(concentration, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_statistics)

cat("\nMaAsLin2 results for selected bile acids:\n")

if (nrow(selected_maaslin2_results) > 0) {
  print(selected_maaslin2_results)
} else {
  cat("MaAsLin2 results file not available.\n")
}

sink()

message("Figure 3B-C bile acid boxplots saved to: outputs/Figure_3/")
message("Summary saved to: ", summary_output_file)