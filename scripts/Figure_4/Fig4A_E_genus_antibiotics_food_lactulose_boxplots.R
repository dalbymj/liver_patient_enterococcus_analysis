# Figure 4A-E: Genus abundance boxplots by antibiotics, food intake, and lactulose
#
# Purpose:
# This script loads genus-level 16S rRNA sequencing data and clinical metadata,
# converts genus abundances to relative abundance percentages, and generates
# Figure 4A-E boxplots showing selected genera stratified by antibiotic exposure,
# oral food intake, and lactulose exposure.
#
# Input files:
# - data/normalised_genus_data_16.04.2024.csv
# - data/clinical_metadata_and_data_renamed_30.01.2025.csv
# - outputs/Figure_4/Fig4_maaslin2_genus_antibiotics_food_lactulose_results.csv
#
# Output files:
# - outputs/Figure_4/Fig4A_Enterococcus_antibiotics_boxplot.pdf
# - outputs/Figure_4/Fig4B_Ruminococcus_antibiotics_boxplot.pdf
# - outputs/Figure_4/Fig4C_Roseburia_food_intake_boxplot.pdf
# - outputs/Figure_4/Fig4D_Enterococcus_food_intake_boxplot.pdf
# - outputs/Figure_4/Fig4E_Bifidobacterium_lactulose_boxplot.pdf
# - outputs/Figure_4/Fig4A_E_boxplot_data.csv
# - outputs/Figure_4/Fig4A_E_boxplot_summary_statistics.csv
# - outputs/Figure_4/Fig4A_E_selected_maaslin2_results.csv
# - outputs/Figure_4/Fig4A_E_genus_boxplot_summary.txt
#
# Run this script from the root folder of the repository:
# liver_patient_enterococcus_analysis/

library(tidyverse)

options(scipen = 999)

# Create output folder if it does not already exist
dir.create("outputs/Figure_4", recursive = TRUE, showWarnings = FALSE)

# Define input files
genus_file <- "data/normalised_genus_data_16.04.2024.csv"
metadata_file <- "data/clinical_metadata_and_data_renamed_30.01.2025.csv"
maaslin2_results_file <- "outputs/Figure_4/Fig4_maaslin2_genus_antibiotics_food_lactulose_results.csv"

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

# Genera plotted in Figure 4A-E
selected_genera <- c(
  "Enterococcus",
  "Ruminococcus",
  "Roseburia",
  "Bifidobacterium"
)

missing_selected_genera <- setdiff(selected_genera, genus_cols)

if (length(missing_selected_genera) > 0) {
  stop(
    "The following selected genera are missing from the input data: ",
    paste(missing_selected_genera, collapse = ", ")
  )
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
        "Healthy_control",
        "No_antibiotics",
        "Receiving_antibiotics"
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
        "Healthy_control",
        "Yes",
        "No"
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
        "Healthy_control",
        "No",
        "Yes"
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
plot_data_wide <- metadata_subsampled %>%
  inner_join(genus_data_percent, by = "study_no")

# Save plot data used by Figure 4A-E
write_csv(
  plot_data_wide,
  "outputs/Figure_4/Fig4A_E_boxplot_data.csv"
)

# Shared visual settings
group_fill_values <- c(
  "Healthy_control" = "lightskyblue",
  "No_antibiotics" = "yellow",
  "Receiving_antibiotics" = "firebrick1",
  "Yes" = "yellow",
  "No" = "firebrick1"
)

group_colour_values <- c(
  "Healthy_control" = "blue4",
  "No_antibiotics" = "yellow4",
  "Receiving_antibiotics" = "red4",
  "Yes" = "yellow4",
  "No" = "red4"
)

group_shape_values <- c(
  "Healthy_control" = 1,
  "No_antibiotics" = 1,
  "Receiving_antibiotics" = 1,
  "Yes" = 1,
  "No" = 1
)

# Plot specifications
plot_specs <- tibble(
  panel = c(
    "Fig4A",
    "Fig4B",
    "Fig4C",
    "Fig4D",
    "Fig4E"
  ),
  feature = c(
    "Enterococcus",
    "Ruminococcus",
    "Roseburia",
    "Enterococcus",
    "Bifidobacterium"
  ),
  grouping_variable = c(
    "Receiving_antibiotics",
    "Receiving_antibiotics",
    "Oral_intake",
    "Oral_intake",
    "Lactulose"
  ),
  y_label = c(
    "Enterococcus (%)",
    "Ruminococcus (%)",
    "Roseburia (%)",
    "Enterococcus (%)",
    "Bifidobacterium (%)"
  ),
  y_limit = c(
    114,
    90,
    32,
    114,
    65
  ),
  width_cm = c(
    3.4,
    3.4,
    3.4,
    3.4,
    3.4
  ),
  height_cm = c(
    6,
    6,
    5.8,
    5.8,
    5.8
  ),
  output_file = c(
    "outputs/Figure_4/Fig4A_Enterococcus_antibiotics_boxplot.pdf",
    "outputs/Figure_4/Fig4B_Ruminococcus_antibiotics_boxplot.pdf",
    "outputs/Figure_4/Fig4C_Roseburia_food_intake_boxplot.pdf",
    "outputs/Figure_4/Fig4D_Enterococcus_food_intake_boxplot.pdf",
    "outputs/Figure_4/Fig4E_Bifidobacterium_lactulose_boxplot.pdf"
  ),
  x_labels = list(
    c(
      "Healthy_control" = "Healthy\ncontrol",
      "No_antibiotics" = "No antibiotics",
      "Receiving_antibiotics" = "Receiving\nantibiotics"
    ),
    c(
      "Healthy_control" = "Healthy\ncontrol",
      "No_antibiotics" = "No antibiotics",
      "Receiving_antibiotics" = "Receiving\nantibiotics"
    ),
    c(
      "Healthy_control" = "Healthy\ncontrol",
      "Yes" = "Oral food\nintake",
      "No" = "No oral\nfood intake"
    ),
    c(
      "Healthy_control" = "Healthy\ncontrol",
      "Yes" = "Oral food\nintake",
      "No" = "No oral\nfood intake"
    ),
    c(
      "Healthy_control" = "Healthy\ncontrol",
      "No" = "No lactulose",
      "Yes" = "Lactulose"
    )
  )
)

# Manual annotations used in the original Figure 4A-E panels
annotation_table <- tribble(
  ~panel, ~text_x, ~text_y, ~label, ~segment_x, ~segment_xend, ~segment_y,
  "Fig4A", 2.5, 108, "p = 0.02", 2, 3, 103,
  "Fig4B", 2.5, 84, "p = 0.02", 2, 3, 80,
  "Fig4C", 2.5, 28.5, "p = 0.03", 2, 3, 27,
  "Fig4D", 2.5, 108, "p = 0.048", 2, 3, 103,
  "Fig4E", 2.5, 61, "p = 0.03", 2, 3, 58
)

# Shared plotting function
make_genus_boxplot <- function(data, spec_row, annotations) {
  
  feature_name <- spec_row$feature
  grouping_variable <- spec_row$grouping_variable
  y_label <- spec_row$y_label
  y_limit <- spec_row$y_limit
  x_labels <- spec_row$x_labels[[1]]
  group_levels <- names(x_labels)
  
  plot_data <- data %>%
    mutate(
      plot_group = factor(
        .data[[grouping_variable]],
        levels = group_levels
      )
    ) %>%
    filter(!is.na(plot_group))
  
  plot_object <- ggplot(
    plot_data,
    aes(
      x = plot_group,
      y = .data[[feature_name]],
      fill = plot_group,
      shape = plot_group
    )
  ) +
    geom_boxplot(
      position = position_dodge(),
      outlier.shape = NA,
      aes(
        colour = plot_group,
        fill = plot_group
      ),
      linewidth = 0.35
    ) +
    geom_point(
      position = position_jitterdodge(jitter.width = 1),
      aes(
        colour = plot_group,
        fill = plot_group
      ),
      size = 1.5,
      stroke = 0.35
    ) +
    scale_fill_manual(
      values = group_fill_values[group_levels],
      drop = FALSE
    ) +
    scale_color_manual(
      values = group_colour_values[group_levels],
      drop = FALSE
    ) +
    scale_shape_manual(
      values = group_shape_values[group_levels],
      drop = FALSE
    ) +
    scale_y_continuous(
      expand = c(0, 0)
    ) +
    coord_cartesian(
      ylim = c(0, y_limit),
      clip = "off"
    ) +
    scale_x_discrete(
      labels = x_labels,
      drop = FALSE
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

# Generate and save each Figure 4A-E panel
for (i in seq_len(nrow(plot_specs))) {
  
  spec_row <- plot_specs[i, ]
  
  panel_annotations <- annotation_table %>%
    filter(panel == spec_row$panel)
  
  panel_plot <- make_genus_boxplot(
    data = plot_data_wide,
    spec_row = spec_row,
    annotations = panel_annotations
  )
  
  ggsave(
    filename = spec_row$output_file,
    plot = panel_plot,
    width = spec_row$width_cm,
    height = spec_row$height_cm,
    units = "cm",
    dpi = 300
  )
}

# Save summary statistics used for checking outputs
summary_statistics <- plot_specs %>%
  select(panel, feature, grouping_variable) %>%
  pmap_dfr(function(panel, feature, grouping_variable) {
    
    plot_data_wide %>%
      select(
        group = all_of(grouping_variable),
        value = all_of(feature)
      ) %>%
      filter(!is.na(group)) %>%
      group_by(group) %>%
      summarise(
        panel = panel,
        feature = feature,
        grouping_variable = grouping_variable,
        n = n(),
        median = median(value, na.rm = TRUE),
        IQR = IQR(value, na.rm = TRUE),
        mean = mean(value, na.rm = TRUE),
        sd = sd(value, na.rm = TRUE),
        max = max(value, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      select(
        panel,
        feature,
        grouping_variable,
        group,
        n,
        median,
        IQR,
        mean,
        sd,
        max
      )
  })

write_csv(
  summary_statistics,
  "outputs/Figure_4/Fig4A_E_boxplot_summary_statistics.csv"
)

# Load MaAsLin2 results if available and save selected results for Figure 4A-E
if (file.exists(maaslin2_results_file)) {
  
  maaslin2_results <- read_csv(maaslin2_results_file, show_col_types = FALSE)
  
  selected_maaslin2_results <- maaslin2_results %>%
    filter(
      feature %in% unique(plot_specs$feature)
    ) %>%
    select(
      feature,
      metadata,
      value,
      coef,
      stderr,
      pval,
      qval
    ) %>%
    arrange(feature, metadata, value)
  
  write_csv(
    selected_maaslin2_results,
    "outputs/Figure_4/Fig4A_E_selected_maaslin2_results.csv"
  )
  
} else {
  
  selected_maaslin2_results <- tibble()
  
  warning(
    "MaAsLin2 results file was not found: ",
    maaslin2_results_file,
    "\nRun Fig4_maaslin2_genus_antibiotics_food_lactulose.R first if MaAsLin2 results are required."
  )
}

# Save summary output
summary_output_file <- "outputs/Figure_4/Fig4A_E_genus_boxplot_summary.txt"

sink(summary_output_file)

cat("Figure 4A-E: Genus abundance boxplots by antibiotics, food intake, and lactulose\n")
cat("================================================================================\n\n")

cat("Input files:\n")
cat(genus_file, "\n")
cat(metadata_file, "\n\n")

cat("Number of samples included:\n")
print(nrow(plot_data_wide))

cat("\nFigure 4A-E plot specifications:\n")
print(plot_specs)

cat("\nSample numbers by antibiotic exposure:\n")
print(table(plot_data_wide$Receiving_antibiotics, useNA = "ifany"))

cat("\nSample numbers by oral food intake:\n")
print(table(plot_data_wide$Oral_intake, useNA = "ifany"))

cat("\nSample numbers by lactulose exposure:\n")
print(table(plot_data_wide$Lactulose, useNA = "ifany"))

cat("\nSummary statistics by panel:\n")
print(summary_statistics)

cat("\nSelected MaAsLin2 results:\n")

if (nrow(selected_maaslin2_results) > 0) {
  print(selected_maaslin2_results)
} else {
  cat("MaAsLin2 results file not available.\n")
}

sink()

message("Figure 4A-E boxplots saved to: outputs/Figure_4/")
message("Summary statistics saved to: outputs/Figure_4/Fig4A_E_boxplot_summary_statistics.csv")
message("Summary saved to: ", summary_output_file)