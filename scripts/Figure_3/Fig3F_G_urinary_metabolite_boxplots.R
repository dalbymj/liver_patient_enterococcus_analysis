# Figure 3F-G: Urinary metabolite boxplots across disease severity groups
#
# Purpose:
# This script loads urinary NMR metabolite data and clinical metadata, then
# generates Figure 3F-G boxplots for selected urinary metabolites associated
# with ACLD disease severity.
#
# Input files:
# - data/Updated_urinary_NMR_GLA_RifSys_8.03.2021 with edited sample_ID.csv
# - data/clinical_metadata_and_data_renamed_30.01.2025.csv
# - outputs/Figure_3/Fig3E_maaslin2_urinary_metabolites_all_results.csv
#
# Output files:
# - outputs/Figure_3/Fig3F_TMAO_boxplot.pdf
# - outputs/Figure_3/Fig3G_Hippurate_boxplot.pdf
# - outputs/Figure_3/Fig3F_G_urinary_metabolite_boxplot_statistics.txt
# - outputs/Figure_3/Fig3F_G_selected_urinary_metabolite_maaslin2_results.csv
#
# Run this script from the root folder of the repository:
# liver_patient_enterococcus_analysis/

library(tidyverse)

options(scipen = 999)

# Create output folder if it does not already exist
dir.create("outputs/Figure_3", recursive = TRUE, showWarnings = FALSE)

# Define input files
urinary_file <- "data/Urinary_NMR_GLA_RifSys_8.03.2021.csv"
metadata_file <- "data/clinical_metadata_and_data_renamed_30.01.2025.csv"
maaslin2_results_file <- "outputs/Figure_3/Fig3E_maaslin2_urinary_metabolites_all_results.csv"

# Check that input files exist
if (!file.exists(urinary_file)) {
  stop("Missing input file: ", urinary_file)
}

if (!file.exists(metadata_file)) {
  stop("Missing input file: ", metadata_file)
}

# Load data
metabolite_data_raw <- read_csv(urinary_file, show_col_types = FALSE)
metadata <- read_csv(metadata_file, show_col_types = FALSE)

# Remove unnamed index columns if present
metabolite_data_raw <- metabolite_data_raw %>%
  select(-any_of("...1"))

metadata <- metadata %>%
  select(-any_of("...1"))

# Standardise sample ID column name
if ("Study no" %in% names(metabolite_data_raw)) {
  names(metabolite_data_raw)[names(metabolite_data_raw) == "Study no"] <- "study_no"
}

# Replace "/" in metabolite names with "_", matching the original script
names(metabolite_data_raw) <- gsub(
  x = names(metabolite_data_raw),
  pattern = "/",
  replacement = "_"
)

# Define urinary metabolite columns
metabolite_start_col <- match("For_TS", names(metabolite_data_raw))
metabolite_end_col <- match("Lac_TS", names(metabolite_data_raw))

if (is.na(metabolite_start_col)) {
  stop("Could not find urinary metabolite column 'For_TS' in metabolite_data_raw.")
}

if (is.na(metabolite_end_col)) {
  stop("Could not find urinary metabolite column 'Lac_TS' in metabolite_data_raw.")
}

metabolite_cols <- names(metabolite_data_raw)[metabolite_start_col:metabolite_end_col]

# Prepare urinary metabolite data
metabolites <- metabolite_data_raw %>%
  select(study_no, all_of(metabolite_cols)) %>%
  mutate(
    across(
      all_of(metabolite_cols),
      ~ as.numeric(.x)
    )
  ) %>%
  drop_na(all_of(metabolite_cols))

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

# Merge metadata and urinary metabolite data
merged_meta_metabolites <- metadata_subsampled %>%
  inner_join(metabolites, by = "study_no") %>%
  filter(!is.na(Inclusion_criteria))

# Selected urinary metabolites shown in Figure 3F-G
selected_metabolites <- tribble(
  ~panel, ~metabolite, ~display_name, ~y_label, ~y_limit, ~width_cm, ~height_cm, ~output_file,
  "Fig3F", "TMAO_TS", "TMAO", "TMAO (a.u.)", 0.145, 3.5, 7, "outputs/Figure_3/Fig3F_TMAO_boxplot.pdf",
  "Fig3G", "Hip7p85_TS", "Hippurate", "Hippurate (a.u.)", 0.038, 3.5, 7, "outputs/Figure_3/Fig3G_Hippurate_boxplot.pdf"
)

# Check that all selected metabolites are present
missing_metabolites <- setdiff(selected_metabolites$metabolite, names(merged_meta_metabolites))

if (length(missing_metabolites) > 0) {
  stop(
    "The following selected urinary metabolites are missing from the input data: ",
    paste(missing_metabolites, collapse = ", ")
  )
}

# Manual annotation positions and labels used in Figure 3F-G
annotation_table <- tribble(
  ~metabolite, ~text_x, ~text_y, ~label, ~segment_x, ~segment_xend, ~segment_y,
  
  "TMAO_TS", 2.5, 0.140, "p < 0.001", 1, 4, 0.135,
  "TMAO_TS", 2.0, 0.130, "p < 0.001", 1, 3, 0.125,
  
  "Hip7p85_TS", 2.5, 0.0370, "p < 0.001", 1, 4, 0.0355,
  "Hip7p85_TS", 2.0, 0.0340, "p < 0.001", 1, 3, 0.0325
)

# Shared plotting function
make_metabolite_boxplot <- function(data, metabolite_name, y_label, y_limit, annotations) {
  
  plot_object <- ggplot(
    data,
    aes(
      x = Inclusion_criteria,
      y = .data[[metabolite_name]],
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
      expand = c(0, 0)
    ) +
    coord_cartesian(
      ylim = c(0, y_limit)
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
        size = 2
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

# Generate and save each urinary metabolite panel
for (i in seq_len(nrow(selected_metabolites))) {
  
  metabolite_name <- selected_metabolites$metabolite[i]
  y_label <- selected_metabolites$y_label[i]
  y_limit <- selected_metabolites$y_limit[i]
  output_file <- selected_metabolites$output_file[i]
  width_cm <- selected_metabolites$width_cm[i]
  height_cm <- selected_metabolites$height_cm[i]
  
  annotations <- annotation_table %>%
    filter(metabolite == metabolite_name)
  
  metabolite_plot <- make_metabolite_boxplot(
    data = merged_meta_metabolites,
    metabolite_name = metabolite_name,
    y_label = y_label,
    y_limit = y_limit,
    annotations = annotations
  )
  
  ggsave(
    filename = output_file,
    plot = metabolite_plot,
    width = width_cm,
    height = height_cm,
    units = "cm",
    dpi = 300
  )
}

# Load MaAsLin2 results if available
if (file.exists(maaslin2_results_file)) {
  
  maaslin2_results <- read_csv(maaslin2_results_file, show_col_types = FALSE)
  
  selected_maaslin2_results <- maaslin2_results %>%
    filter(feature %in% selected_metabolites$metabolite) %>%
    select(any_of(c("feature", "metadata", "value", "coef", "pval", "qval"))) %>%
    arrange(feature, value)
  
  write_csv(
    selected_maaslin2_results,
    "outputs/Figure_3/Fig3F_G_selected_urinary_metabolite_maaslin2_results.csv"
  )
  
} else {
  
  selected_maaslin2_results <- tibble()
  
  warning(
    "MaAsLin2 results file was not found: ",
    maaslin2_results_file,
    "\nRun Fig3E_maaslin2_urinary_metabolites_differential_abundance.R first if MaAsLin2 results are required."
  )
}

# Save summary output
summary_output_file <- "outputs/Figure_3/Fig3F_G_urinary_metabolite_boxplot_statistics.txt"

sink(summary_output_file)

cat("Figure 3F-G: Urinary metabolite boxplots\n")
cat("========================================\n\n")

cat("Input files:\n")
cat(urinary_file, "\n")
cat(metadata_file, "\n\n")

cat("Sample numbers by group:\n")
print(table(merged_meta_metabolites$Inclusion_criteria, useNA = "ifany"))

cat("\nSelected urinary metabolites:\n")
print(selected_metabolites)

cat("\nSummary statistics by group and urinary metabolite:\n")
summary_statistics <- merged_meta_metabolites %>%
  select(Inclusion_criteria, all_of(selected_metabolites$metabolite)) %>%
  pivot_longer(
    cols = all_of(selected_metabolites$metabolite),
    names_to = "metabolite",
    values_to = "abundance"
  ) %>%
  group_by(metabolite, Inclusion_criteria) %>%
  summarise(
    n = n(),
    median = median(abundance, na.rm = TRUE),
    IQR = IQR(abundance, na.rm = TRUE),
    mean = mean(abundance, na.rm = TRUE),
    sd = sd(abundance, na.rm = TRUE),
    max = max(abundance, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_statistics)

cat("\nMaAsLin2 results for selected urinary metabolites:\n")

if (nrow(selected_maaslin2_results) > 0) {
  print(selected_maaslin2_results)
} else {
  cat("MaAsLin2 results file not available.\n")
}

sink()

message("Figure 3F-G urinary metabolite boxplots saved to: outputs/Figure_3/")
message("Summary saved to: ", summary_output_file)