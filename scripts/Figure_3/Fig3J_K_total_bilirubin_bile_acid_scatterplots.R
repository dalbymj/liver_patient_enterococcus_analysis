# Figure 3J-K: Total bilirubin and bile acid scatterplots
#
# Purpose:
# This script loads clinical metadata and plasma bile acid data, applies
# bestNormalize transformation to total bilirubin and selected bile acids,
# calculates Pearson correlations, and generates Figure 3J-K scatterplots.
#
# Input files:
# - data/clinical_metadata_and_data_renamed_30.01.2025.csv
# - data/Corrected bile acid data 01.12.2020.csv
#
# Output files:
# - outputs/Figure_3/Fig3J_total_bilirubin_glycochenodeoxycholic_acid_scatterplot.pdf
# - outputs/Figure_3/Fig3K_total_bilirubin_taurochenodeoxycholic_acid_scatterplot.pdf
# - outputs/Figure_3/Fig3J_K_total_bilirubin_bile_acid_scatterplot_data.csv
# - outputs/Figure_3/Fig3J_K_total_bilirubin_bile_acid_scatterplot_statistics.csv
# - outputs/Figure_3/Fig3J_K_bestNormalize_transformation_log.csv
# - outputs/Figure_3/Fig3J_K_total_bilirubin_bile_acid_scatterplot_summary.txt
#
# Run this script from the root folder of the repository:
# liver_patient_enterococcus_analysis/

library(tidyverse)
library(bestNormalize)

options(scipen = 999)

set.seed(123)

# Create output folder if it does not already exist
dir.create("outputs/Figure_3", recursive = TRUE, showWarnings = FALSE)

# Define input files
metadata_file <- "data/clinical_metadata_and_data_renamed_30.01.2025.csv"
bile_file <- "data/Corrected bile acid data 01.12.2020.csv"

# Check that input files exist
if (!file.exists(metadata_file)) {
  stop("Missing input file: ", metadata_file)
}

if (!file.exists(bile_file)) {
  stop("Missing input file: ", bile_file)
}

# Load data
metadata_raw <- read_csv(metadata_file, show_col_types = FALSE)
bile_data_raw <- read_csv(bile_file, show_col_types = FALSE)

# Remove unnamed index columns if present
metadata_raw <- metadata_raw %>%
  select(-any_of("...1"))

bile_data_raw <- bile_data_raw %>%
  select(-any_of("...1"))

# Bile acids plotted in Figure 3J-K
selected_bile_acids <- c(
  "Glycochenodeoxycholic_acid",
  "Taurochenodeoxycholic_acid"
)

# Check required columns are present
missing_metadata_cols <- setdiff(c("study_no", "Total_Br_SI_units"), names(metadata_raw))
missing_bile_cols <- setdiff(c("study_no", selected_bile_acids), names(bile_data_raw))

if (length(missing_metadata_cols) > 0) {
  stop(
    "The following metadata columns are missing: ",
    paste(missing_metadata_cols, collapse = ", ")
  )
}

if (length(missing_bile_cols) > 0) {
  stop(
    "The following bile acid columns are missing: ",
    paste(missing_bile_cols, collapse = ", ")
  )
}

# Prepare clinical data
clinical_data <- metadata_raw %>%
  select(
    study_no,
    Total_Br_SI_units
  ) %>%
  rename(
    Total_Bilirubin = Total_Br_SI_units
  ) %>%
  mutate(
    Total_Bilirubin = as.numeric(Total_Bilirubin)
  )

# Prepare bile acid data
bile_data <- bile_data_raw %>%
  select(study_no, all_of(selected_bile_acids)) %>%
  mutate(
    across(
      all_of(selected_bile_acids),
      ~ as.numeric(.x)
    )
  )

# Merge data and retain complete cases for plotted variables
merged_data <- clinical_data %>%
  inner_join(bile_data, by = "study_no") %>%
  filter(
    complete.cases(
      across(c(Total_Bilirubin, all_of(selected_bile_acids)))
    )
  )

# Prepare analysis table
analysis_table_raw <- merged_data %>%
  select(
    study_no,
    Total_Bilirubin,
    all_of(selected_bile_acids)
  ) %>%
  column_to_rownames("study_no")

# Function to apply bestNormalize to each column
apply_best_normalize <- function(data_table, data_type) {
  
  transformed_table <- data_table
  transformation_log <- tibble(
    data_type = character(),
    feature = character(),
    transformation = character()
  )
  
  for (feature_name in names(data_table)) {
    
    feature_values <- data_table[[feature_name]]
    
    if (all(is.na(feature_values)) || sd(feature_values, na.rm = TRUE) == 0) {
      warning("Skipping transformation for ", feature_name, " because values are missing or constant.")
      next
    }
    
    transformed_feature <- bestNormalize::bestNormalize(
      feature_values,
      quiet = TRUE
    )
    
    transformed_table[[feature_name]] <- transformed_feature$x.t
    
    transformation_log <- transformation_log %>%
      add_row(
        data_type = data_type,
        feature = feature_name,
        transformation = class(transformed_feature$chosen_transform)[1]
      )
  }
  
  list(
    data = transformed_table,
    log = transformation_log
  )
}

# Transform variables
transformed_output <- apply_best_normalize(
  data_table = analysis_table_raw,
  data_type = "scatterplot_variable"
)

plot_data <- transformed_output$data %>%
  rownames_to_column("study_no")

transformation_log <- transformed_output$log

write_csv(
  transformation_log,
  "outputs/Figure_3/Fig3J_K_bestNormalize_transformation_log.csv"
)

write_csv(
  plot_data,
  "outputs/Figure_3/Fig3J_K_total_bilirubin_bile_acid_scatterplot_data.csv"
)

# Plot specifications
scatterplot_specs <- tribble(
  ~panel, ~bile_acid, ~x_label, ~y_label, ~output_file,
  "Fig3J", "Glycochenodeoxycholic_acid", "Normalised\nTotal Bilirubin", "Normalised\nGlycochenodeoxycholic acid", "outputs/Figure_3/Fig3J_total_bilirubin_glycochenodeoxycholic_acid_scatterplot.pdf",
  "Fig3K", "Taurochenodeoxycholic_acid", "Normalised\nTotal Bilirubin", "Normalised\nTaurochenodeoxycholic acid", "outputs/Figure_3/Fig3K_total_bilirubin_taurochenodeoxycholic_acid_scatterplot.pdf"
)

# Function to calculate statistics for one scatterplot
calculate_scatterplot_stats <- function(data, bile_acid_name) {
  
  complete_data <- data %>%
    select(Total_Bilirubin, all_of(bile_acid_name)) %>%
    drop_na()
  
  correlation_test <- cor.test(
    complete_data$Total_Bilirubin,
    complete_data[[bile_acid_name]],
    method = "pearson"
  )
  
  fit <- lm(
    complete_data[[bile_acid_name]] ~ complete_data$Total_Bilirubin
  )
  
  tibble(
    bile_acid = bile_acid_name,
    n = nrow(complete_data),
    r = as.numeric(correlation_test$estimate),
    p_value = correlation_test$p.value,
    r_squared = summary(fit)$r.squared
  )
}

# Calculate statistics
scatterplot_statistics <- map_dfr(
  scatterplot_specs$bile_acid,
  ~ calculate_scatterplot_stats(plot_data, .x)
)

write_csv(
  scatterplot_statistics,
  "outputs/Figure_3/Fig3J_K_total_bilirubin_bile_acid_scatterplot_statistics.csv"
)

# Function to generate one scatterplot
make_scatterplot <- function(data, bile_acid_name, x_label, y_label, statistics_table) {
  
  stats_row <- statistics_table %>%
    filter(bile_acid == bile_acid_name)
  
  label_text <- paste0(
    "r = ", round(stats_row$r, 2),
    "\np = ", format.pval(stats_row$p_value, digits = 2),
    "\nn = ", stats_row$n
  )
  
  ggplot(
    data,
    aes(
      x = Total_Bilirubin,
      y = .data[[bile_acid_name]]
    )
  ) +
    geom_point(
      shape = 21,
      size = 1
    ) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      colour = "black",
      fill = "grey80",
      linewidth = 0.4
    ) +
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = label_text,
      hjust = -0.1,
      vjust = 1.1,
      size = 2
    ) +
    coord_cartesian(
      clip = "off"
    ) +
    labs(
      x = x_label,
      y = y_label
    ) +
    theme(
      plot.title = element_text(size = 8),
      legend.position = "none",
      legend.key.size = unit(0.3, "cm"),
      legend.direction = "vertical",
      legend.spacing.x = unit(0.05, "cm"),
      plot.margin = margin(t = 5, r = 0, b = 0, l = 0),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
      legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
      legend.title = element_blank(),
      legend.text = element_text(size = 7),
      legend.background = element_blank(),
      legend.key = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.1, "cm"),
      strip.text.x = element_text(size = 8, face = "bold"),
      axis.title.x = element_text(size = 8),
      axis.text.x = element_text(size = 8, colour = "black", angle = 0, hjust = 0.5, vjust = 0.5),
      axis.title.y = element_text(size = 8, colour = "black"),
      axis.text.y = element_text(size = 8, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      aspect.ratio = 1
    )
}

# Generate and save scatterplots
for (i in seq_len(nrow(scatterplot_specs))) {
  
  bile_acid_name <- scatterplot_specs$bile_acid[i]
  x_label <- scatterplot_specs$x_label[i]
  y_label <- scatterplot_specs$y_label[i]
  output_file <- scatterplot_specs$output_file[i]
  
  scatterplot <- make_scatterplot(
    data = plot_data,
    bile_acid_name = bile_acid_name,
    x_label = x_label,
    y_label = y_label,
    statistics_table = scatterplot_statistics
  )
  
  ggsave(
    filename = output_file,
    plot = scatterplot,
    width = 4.5,
    height = 4.5,
    units = "cm",
    dpi = 300
  )
}

# Save summary output
summary_output_file <- "outputs/Figure_3/Fig3J_K_total_bilirubin_bile_acid_scatterplot_summary.txt"

sink(summary_output_file)

cat("Figure 3J-K: Total bilirubin and bile acid scatterplots\n")
cat("=======================================================\n\n")

cat("Input files:\n")
cat(metadata_file, "\n")
cat(bile_file, "\n\n")

cat("Number of matched complete-case samples:\n")
print(nrow(merged_data))

cat("\nVariables transformed and plotted:\n")
print(names(analysis_table_raw))

cat("\nTransformation log:\n")
print(transformation_log)

cat("\nScatterplot statistics:\n")
print(scatterplot_statistics)

sink()

message("Figure 3J saved to: outputs/Figure_3/Fig3J_total_bilirubin_glycochenodeoxycholic_acid_scatterplot.pdf")
message("Figure 3K saved to: outputs/Figure_3/Fig3K_total_bilirubin_taurochenodeoxycholic_acid_scatterplot.pdf")
message("Scatterplot statistics saved to: outputs/Figure_3/Fig3J_K_total_bilirubin_bile_acid_scatterplot_statistics.csv")
message("Summary saved to: ", summary_output_file)