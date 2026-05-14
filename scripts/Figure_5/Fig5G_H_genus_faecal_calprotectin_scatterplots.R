# Figure 5G-H: Genus and faecal calprotectin scatter plots
#
# Purpose:
# This script loads genus-level 16S rRNA sequencing data and faecal
# calprotectin data, applies bestNormalize transformation, and generates
# scatter plots for selected genus-versus-FCAL associations.
#
# Input files:
# - data/normalised_genus_data_16.04.2024.csv
# - data/data/FCAL_GLA_RifSys_combined_dataset_20.01.2022.csv
#
# Output files:
# - outputs/Figure_5/Fig5G_Lachnospira_faecal_calprotectin_scatterplot.pdf
# - outputs/Figure_5/Fig5H_Enterococcus_faecal_calprotectin_scatterplot.pdf
# - outputs/Figure_5/Fig5G_H_genus_faecal_calprotectin_scatterplot_data.csv
# - outputs/Figure_5/Fig5G_H_genus_faecal_calprotectin_scatterplot_statistics.csv
# - outputs/Figure_5/Fig5G_H_bestNormalize_transformation_log.csv
# - outputs/Figure_5/Fig5G_H_genus_faecal_calprotectin_scatterplot_summary.txt
#
# Run this script from the root folder of the repository:
# liver_patient_enterococcus_analysis/

library(tidyverse)
library(bestNormalize)
library(grid)

options(scipen = 999)

set.seed(123)

# Create output folder if it does not already exist
dir.create("outputs/Figure_5", recursive = TRUE, showWarnings = FALSE)

# Define input files
genus_file <- "data/normalised_genus_data_16.04.2024.csv"
fcal_file <- "data/FCAL_GLA_RifSys_combined_dataset_20.01.2022.csv"

# Check that input files exist
if (!file.exists(genus_file)) {
  stop("Missing input file: ", genus_file)
}

if (!file.exists(fcal_file)) {
  stop("Missing input file: ", fcal_file)
}

# Load data
genus_data_raw <- read_csv(genus_file, show_col_types = FALSE)
fcal_raw <- read_csv(fcal_file, show_col_types = FALSE)

# Remove unnamed index columns if present
genus_data_raw <- genus_data_raw %>%
  select(-any_of("...1"))

fcal_raw <- fcal_raw %>%
  select(-any_of("...1"))

# Rename FCAL columns if needed
if ("Study no" %in% names(fcal_raw)) {
  names(fcal_raw)[names(fcal_raw) == "Study no"] <- "study_no"
}

# Check required FCAL columns
required_fcal_cols <- c("study_no", "FCAL")
missing_fcal_cols <- setdiff(required_fcal_cols, names(fcal_raw))

if (length(missing_fcal_cols) > 0) {
  stop(
    "The following FCAL columns are missing: ",
    paste(missing_fcal_cols, collapse = ", ")
  )
}

# Define genus columns
genus_start_col <- match("Mucispirillum", names(genus_data_raw))

if (is.na(genus_start_col)) {
  stop("Could not find genus column 'Mucispirillum' in genus_data_raw.")
}

genus_cols <- names(genus_data_raw)[genus_start_col:ncol(genus_data_raw)]

# Genera used for Figure 5G-H scatter plots
scatter_genera <- c(
  "Lachnospira",
  "Enterococcus"
)

missing_scatter_genera <- setdiff(scatter_genera, genus_cols)

if (length(missing_scatter_genera) > 0) {
  stop(
    "The following scatter plot genera are missing from the input data: ",
    paste(missing_scatter_genera, collapse = ", ")
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

# Prepare FCAL data
fcal_data <- fcal_raw %>%
  select(study_no, FCAL) %>%
  mutate(
    FCAL = as.numeric(FCAL)
  ) %>%
  filter(complete.cases(.))

# Merge data
merged_data <- fcal_data %>%
  inner_join(genus_data_percent, by = "study_no")

if (nrow(merged_data) == 0) {
  stop("No matched samples were found after merging genus and FCAL data.")
}

# Select variables for scatter plots
scatter_data_raw <- merged_data %>%
  select(
    study_no,
    FCAL,
    all_of(scatter_genera)
  ) %>%
  drop_na()

# Function to apply bestNormalize to selected columns.
# If bestNormalize fails, fall back to simple z-score scaling.
apply_best_normalize <- function(data_table, variables, data_type) {
  
  transformed_table <- data_table
  
  transformation_log <- tibble(
    data_type = character(),
    feature = character(),
    transformation = character(),
    note = character()
  )
  
  for (feature_name in variables) {
    
    message("Transforming: ", feature_name)
    
    feature_values <- transformed_table[[feature_name]]
    non_missing_values <- feature_values[!is.na(feature_values)]
    unique_values <- unique(non_missing_values)
    
    if (
      length(non_missing_values) < 3 ||
      length(unique_values) < 2 ||
      sd(feature_values, na.rm = TRUE) == 0
    ) {
      
      warning(
        "Skipping transformation for ",
        feature_name,
        " because it has too few usable or distinct values."
      )
      
      transformation_log <- transformation_log %>%
        add_row(
          data_type = data_type,
          feature = feature_name,
          transformation = "not_transformed",
          note = "too_few_usable_or_distinct_values"
        )
      
      next
    }
    
    transformed_feature <- tryCatch(
      {
        bestNormalize::bestNormalize(
          feature_values,
          quiet = TRUE
        )
      },
      error = function(e) {
        warning(
          "bestNormalize failed for ",
          feature_name,
          ". Falling back to z-score scaling. Error was: ",
          conditionMessage(e)
        )
        NULL
      }
    )
    
    if (!is.null(transformed_feature)) {
      
      transformed_table[[feature_name]] <- transformed_feature$x.t
      
      transformation_log <- transformation_log %>%
        add_row(
          data_type = data_type,
          feature = feature_name,
          transformation = class(transformed_feature$chosen_transform)[1],
          note = "bestNormalize_success"
        )
      
    } else {
      
      transformed_table[[feature_name]] <- as.numeric(scale(feature_values))
      
      transformation_log <- transformation_log %>%
        add_row(
          data_type = data_type,
          feature = feature_name,
          transformation = "z_score",
          note = "bestNormalize_failed_z_score_used"
        )
    }
  }
  
  list(
    data = transformed_table,
    log = transformation_log
  )
}

# Transform variables
variables_to_transform <- c(
  "FCAL",
  scatter_genera
)

transformed_output <- apply_best_normalize(
  data_table = scatter_data_raw,
  variables = variables_to_transform,
  data_type = "scatterplot_variable"
)

scatter_data <- transformed_output$data
transformation_log <- transformed_output$log

write_csv(
  scatter_data,
  "outputs/Figure_5/Fig5G_H_genus_faecal_calprotectin_scatterplot_data.csv"
)

write_csv(
  transformation_log,
  "outputs/Figure_5/Fig5G_H_bestNormalize_transformation_log.csv"
)

# Plot specifications
scatter_specs <- tribble(
  ~panel, ~feature, ~x_label, ~y_label, ~output_file, ~label_position,
  "Fig5G", "Lachnospira", "Normalised\nLachnospira (%)", "Normalised FCAL (\u00B5g/g)", "outputs/Figure_5/Fig5G_Lachnospira_faecal_calprotectin_scatterplot.pdf", "top_right",
  "Fig5H", "Enterococcus", "Normalised\nEnterococcus (%)", "Normalised FCAL (\u00B5g/g)", "outputs/Figure_5/Fig5H_Enterococcus_faecal_calprotectin_scatterplot.pdf", "top_left"
)

# Calculate scatter plot statistics
calculate_scatter_statistics <- function(data, feature_column) {
  
  complete_data <- tibble(
    x_value = data[[feature_column]],
    y_value = data[["FCAL"]]
  ) %>%
    drop_na()
  
  correlation_test <- cor.test(
    complete_data$x_value,
    complete_data$y_value,
    method = "pearson"
  )
  
  linear_model <- lm(
    y_value ~ x_value,
    data = complete_data
  )
  
  tibble(
    n = nrow(complete_data),
    r = as.numeric(correlation_test$estimate),
    p_value = correlation_test$p.value,
    r_squared = summary(linear_model)$r.squared,
    adjusted_r_squared = summary(linear_model)$adj.r.squared
  )
}

scatter_statistics <- map_dfr(
  seq_len(nrow(scatter_specs)),
  function(i) {
    
    spec_row <- scatter_specs[i, ]
    feature_column <- spec_row$feature[[1]]
    
    calculate_scatter_statistics(
      data = scatter_data,
      feature_column = feature_column
    ) %>%
      mutate(
        panel = spec_row$panel[[1]],
        feature = feature_column,
        .before = 1
      )
  }
)

write_csv(
  scatter_statistics,
  "outputs/Figure_5/Fig5G_H_genus_faecal_calprotectin_scatterplot_statistics.csv"
)

# Shared scatter plotting function
make_scatterplot <- function(data, spec_row, statistics_row) {
  
  feature_column <- spec_row$feature[[1]]
  
  label_text <- paste0(
    "r = ", round(statistics_row$r[[1]], 2),
    "\np = ", format.pval(statistics_row$p_value[[1]], digits = 2),
    "\nn = ", statistics_row$n[[1]]
  )
  
  if (spec_row$label_position[[1]] == "top_right") {
    label_x <- Inf
    label_hjust <- 1.1
  } else {
    label_x <- -Inf
    label_hjust <- -0.1
  }
  
  # Important:
  # Create explicit columns before ggplot.
  # This avoids lazy evaluation errors involving feature_column or feature_name.
  plot_data <- tibble(
    x_value = data[[feature_column]],
    y_value = data[["FCAL"]]
  ) %>%
    drop_na()
  
  ggplot(
    plot_data,
    aes(
      x = x_value,
      y = y_value
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
      x = label_x,
      y = Inf,
      label = label_text,
      hjust = label_hjust,
      vjust = 1.1,
      size = 2
    ) +
    coord_cartesian(
      clip = "off"
    ) +
    labs(
      x = spec_row$x_label[[1]],
      y = spec_row$y_label[[1]]
    ) +
    theme(
      plot.title = element_text(size = 8),
      legend.position = "none",
      legend.key.size = unit(0.3, "cm"),
      legend.direction = "vertical",
      legend.spacing.x = unit(0.05, "cm"),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
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

# Generate and save scatter plots
for (i in seq_len(nrow(scatter_specs))) {
  
  spec_row <- scatter_specs[i, ]
  
  statistics_row <- scatter_statistics %>%
    filter(panel == spec_row$panel[[1]])
  
  scatter_plot <- make_scatterplot(
    data = scatter_data,
    spec_row = spec_row,
    statistics_row = statistics_row
  )
  
  ggsave(
    filename = spec_row$output_file[[1]],
    plot = scatter_plot,
    width = 4,
    height = 4,
    units = "cm",
    dpi = 300
  )
}

# Save summary output
summary_output_file <- "outputs/Figure_5/Fig5G_H_genus_faecal_calprotectin_scatterplot_summary.txt"

sink(summary_output_file)

cat("Figure 5G-H: Genus and faecal calprotectin scatter plots\n")
cat("=========================================================\n\n")

cat("Input files:\n")
cat(genus_file, "\n")
cat(fcal_file, "\n\n")

cat("Number of matched complete-case samples:\n")
print(nrow(scatter_data))

cat("\nVariables transformed with bestNormalize:\n")
print(variables_to_transform)

cat("\nTransformation log:\n")
print(transformation_log)

cat("\nScatter plot statistics:\n")
print(scatter_statistics)

sink()

message("Scatter plots saved to: outputs/Figure_5/")
message("Scatter plot data saved to: outputs/Figure_5/Fig5G_H_genus_faecal_calprotectin_scatterplot_data.csv")
message("Statistics saved to: outputs/Figure_5/Fig5G_H_genus_faecal_calprotectin_scatterplot_statistics.csv")
message("Summary saved to: ", summary_output_file)