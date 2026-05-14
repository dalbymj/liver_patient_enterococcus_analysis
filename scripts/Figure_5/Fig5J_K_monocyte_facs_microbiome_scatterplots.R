# Figure 5J-K: Enterococcus and monocyte TNF-alpha response scatter plots
#
# Purpose:
# This script loads monocyte FACS data and genus-level 16S rRNA sequencing data,
# applies bestNormalize transformation, and generates scatter plots for
# Enterococcus abundance against LPS- and PAM3-induced TNF-alpha responses.
#
# Input files:
# - data/Copy of FACS data_GLA TLR249 stimulation_07.09.2022.csv
# - data/normalised_genus_data_12.08.2020.csv
#
# Output files:
# - outputs/Figure_5/Fig5J_Enterococcus_LPS_TNFa_scatterplot.pdf
# - outputs/Figure_5/Fig5K_Enterococcus_PAM3_TNFa_scatterplot.pdf
# - outputs/Figure_5/Fig5J_K_monocyte_facs_microbiome_scatterplot_data.csv
# - outputs/Figure_5/Fig5J_K_monocyte_facs_microbiome_scatterplot_statistics.csv
# - outputs/Figure_5/Fig5J_K_bestNormalize_transformation_log.csv
# - outputs/Figure_5/Fig5J_K_monocyte_facs_microbiome_scatterplot_summary.txt
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
facs_file <- "data/FACS_data_GLA TLR249_stimulation_07.09.2022.csv"
genus_file <- "data/normalised_genus_data_16.04.2024.csv"

# Check that input files exist
input_files <- c(facs_file, genus_file)

missing_input_files <- input_files[!file.exists(input_files)]

if (length(missing_input_files) > 0) {
  stop(
    "Missing input file(s): ",
    paste(missing_input_files, collapse = ", ")
  )
}

# Helper function to clean FACS study IDs
clean_facs_study_no <- function(study_no) {
  
  study_no <- as.character(study_no)
  
  study_no <- gsub("GLA", "GLA_", study_no)
  study_no <- gsub("GLA_HC", "GLA HC00", study_no)
  
  study_no <- gsub("GLA_30", "GLA_030", study_no)
  study_no <- gsub("GLA_35", "GLA_035", study_no)
  study_no <- gsub("GLA_41", "GLA_041", study_no)
  study_no <- gsub("GLA_42", "GLA_042", study_no)
  study_no <- gsub("GLA_43", "GLA_043", study_no)
  study_no <- gsub("GLA_51", "GLA_051", study_no)
  study_no <- gsub("GLA_52", "GLA_052", study_no)
  study_no <- gsub("GLA_53", "GLA_053", study_no)
  study_no <- gsub("GLA_55", "GLA_055", study_no)
  study_no <- gsub("GLA_58", "GLA_058", study_no)
  study_no <- gsub("GLA_62", "GLA_062", study_no)
  study_no <- gsub("GLA_65", "GLA_065", study_no)
  study_no <- gsub("GLA_71", "GLA_071", study_no)
  study_no <- gsub("GLA_77", "GLA_077", study_no)
  study_no <- gsub("GLA_82", "GLA_082", study_no)
  study_no <- gsub("GLA_88", "GLA_088", study_no)
  study_no <- gsub("GLA_89", "GLA_089", study_no)
  study_no <- gsub("GLA_90", "GLA_090", study_no)
  study_no <- gsub("GLA_92", "GLA_092", study_no)
  study_no <- gsub("GLA_93", "GLA_093", study_no)
  
  study_no
}

# Helper function to apply bestNormalize safely
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
    non_missing_index <- !is.na(feature_values)
    non_missing_values <- feature_values[non_missing_index]
    unique_values <- unique(non_missing_values)
    
    if (
      length(non_missing_values) < 3 ||
      length(unique_values) < 2 ||
      sd(non_missing_values, na.rm = TRUE) == 0
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
          non_missing_values,
          quiet = TRUE,
          out_of_sample = FALSE
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
      
      transformed_values <- rep(NA_real_, length(feature_values))
      transformed_values[non_missing_index] <- transformed_feature$x.t
      transformed_table[[feature_name]] <- transformed_values
      
      transformation_log <- transformation_log %>%
        add_row(
          data_type = data_type,
          feature = feature_name,
          transformation = class(transformed_feature$chosen_transform)[1],
          note = "bestNormalize_success"
        )
      
    } else {
      
      transformed_values <- rep(NA_real_, length(feature_values))
      transformed_values[non_missing_index] <- as.numeric(scale(non_missing_values))
      transformed_table[[feature_name]] <- transformed_values
      
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

# Load FACS data
facs_raw <- read_csv(facs_file, show_col_types = FALSE) %>%
  select(-any_of("...1"))

facs_variables <- c(
  "lps_tnfa",
  "pam3_tnfa"
)

missing_facs_variables <- setdiff(c("donor", facs_variables), names(facs_raw))

if (length(missing_facs_variables) > 0) {
  stop(
    "The following FACS columns are missing: ",
    paste(missing_facs_variables, collapse = ", ")
  )
}

facs_data <- facs_raw %>%
  rename(study_no = donor) %>%
  select(
    study_no,
    all_of(facs_variables)
  ) %>%
  mutate(
    study_no = clean_facs_study_no(study_no),
    across(
      all_of(facs_variables),
      ~ as.numeric(.x)
    )
  ) %>%
  filter(complete.cases(across(all_of(facs_variables))))

# Load genus data
genus_raw <- read_csv(genus_file, show_col_types = FALSE) %>%
  select(-any_of("...1"))

genus_start_col <- match("Mucispirillum", names(genus_raw))
genus_end_col <- match("Methanosphaera", names(genus_raw))

if (is.na(genus_start_col)) {
  stop("Could not find genus column 'Mucispirillum' in genus data.")
}

if (is.na(genus_end_col)) {
  stop("Could not find genus column 'Methanosphaera' in genus data.")
}

genus_cols <- names(genus_raw)[genus_start_col:genus_end_col]

if (!"Enterococcus" %in% genus_cols) {
  stop("Could not find genus column 'Enterococcus' in genus data.")
}

genus_data <- genus_raw %>%
  select(study_no, all_of(genus_cols)) %>%
  mutate(
    across(
      all_of(genus_cols),
      ~ as.numeric(.x)
    )
  ) %>%
  filter(complete.cases(across(all_of(genus_cols)))) %>%
  mutate(
    row_total = rowSums(across(all_of(genus_cols)), na.rm = TRUE)
  )

if (any(genus_data$row_total == 0)) {
  stop("At least one sample has a total genus abundance of zero.")
}

genus_percent <- genus_data %>%
  mutate(
    across(
      all_of(genus_cols),
      ~ .x / row_total * 100
    )
  ) %>%
  select(
    study_no,
    Enterococcus
  )

# Merge data
merged_data <- facs_data %>%
  inner_join(
    genus_percent,
    by = "study_no"
  )

if (nrow(merged_data) == 0) {
  stop("No matched samples were found after merging FACS and genus data.")
}

# Transform variables
variables_to_transform <- c(
  "Enterococcus",
  facs_variables
)

transformed_output <- apply_best_normalize(
  data_table = merged_data,
  variables = variables_to_transform,
  data_type = "scatterplot_variable"
)

scatter_data <- transformed_output$data
transformation_log <- transformed_output$log

write_csv(
  scatter_data,
  "outputs/Figure_5/Fig5J_K_monocyte_facs_microbiome_scatterplot_data.csv"
)

write_csv(
  transformation_log,
  "outputs/Figure_5/Fig5J_K_bestNormalize_transformation_log.csv"
)

# Plot specifications
scatter_specs <- tribble(
  ~panel, ~x_feature, ~y_feature, ~plot_title, ~x_label, ~y_label, ~output_file, ~width_cm, ~height_cm, ~label_position,
  "Fig5J", "Enterococcus", "lps_tnfa", "Stimuli: LPS", "Normalised\nEnterococcus (%)", "Normalised % TNF-alpha+\nof CD14+ monocytes", "outputs/Figure_5/Fig5J_Enterococcus_LPS_TNFa_scatterplot.pdf", 4, 3.65, "top_left",
  "Fig5K", "Enterococcus", "pam3_tnfa", "Stimuli: PAM3", "Normalised\nEnterococcus (%)", "Normalised % TNF-alpha+\nof CD14+ monocytes", "outputs/Figure_5/Fig5K_Enterococcus_PAM3_TNFa_scatterplot.pdf", 4, 4, "top_left"
)

# Calculate scatter plot statistics
calculate_scatter_statistics <- function(data, x_feature, y_feature) {
  
  complete_data <- tibble(
    x_value = data[[x_feature]],
    y_value = data[[y_feature]]
  ) %>%
    drop_na()
  
  if (
    nrow(complete_data) < 3 ||
    sd(complete_data$x_value) == 0 ||
    sd(complete_data$y_value) == 0
  ) {
    return(
      tibble(
        n = nrow(complete_data),
        r = NA_real_,
        p_value = NA_real_,
        r_squared = NA_real_,
        adjusted_r_squared = NA_real_
      )
    )
  }
  
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
    
    calculate_scatter_statistics(
      data = scatter_data,
      x_feature = spec_row$x_feature[[1]],
      y_feature = spec_row$y_feature[[1]]
    ) %>%
      mutate(
        panel = spec_row$panel[[1]],
        x_feature = spec_row$x_feature[[1]],
        y_feature = spec_row$y_feature[[1]],
        .before = 1
      )
  }
)

write_csv(
  scatter_statistics,
  "outputs/Figure_5/Fig5J_K_monocyte_facs_microbiome_scatterplot_statistics.csv"
)

# Shared scatter plotting function
make_scatterplot <- function(data, spec_row, statistics_row) {
  
  x_feature <- spec_row$x_feature[[1]]
  y_feature <- spec_row$y_feature[[1]]
  
  if (is.na(statistics_row$r[[1]])) {
    label_text <- paste0("n = ", statistics_row$n[[1]])
  } else {
    label_text <- paste0(
      "r = ", round(statistics_row$r[[1]], 2),
      "\np = ", format.pval(statistics_row$p_value[[1]], digits = 2),
      "\nn = ", statistics_row$n[[1]]
    )
  }
  
  if (spec_row$label_position[[1]] == "top_right") {
    label_x <- Inf
    label_hjust <- 1.1
  } else {
    label_x <- -Inf
    label_hjust <- -0.1
  }
  
  plot_data <- tibble(
    x_value = data[[x_feature]],
    y_value = data[[y_feature]]
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
      title = spec_row$plot_title[[1]],
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
    width = spec_row$width_cm[[1]],
    height = spec_row$height_cm[[1]],
    units = "cm",
    dpi = 300
  )
}

# Save summary output
summary_output_file <- "outputs/Figure_5/Fig5J_K_monocyte_facs_microbiome_scatterplot_summary.txt"

sink(summary_output_file)

cat("Figure 5J-K: Enterococcus and monocyte TNF-alpha response scatter plots\n")
cat("=======================================================================\n\n")

cat("Input files:\n")
cat(facs_file, "\n")
cat(genus_file, "\n\n")

cat("Number of matched samples:\n")
print(nrow(scatter_data))

cat("\nVariables transformed with bestNormalize:\n")
print(variables_to_transform)

cat("\nTransformation log:\n")
print(transformation_log)

cat("\nScatter plot statistics:\n")
print(scatter_statistics)

sink()

message("Scatter plots saved to: outputs/Figure_5/")
message("Scatter plot data saved to: outputs/Figure_5/Fig5J_K_monocyte_facs_microbiome_scatterplot_data.csv")
message("Statistics saved to: outputs/Figure_5/Fig5J_K_monocyte_facs_microbiome_scatterplot_statistics.csv")
message("Summary saved to: ", summary_output_file)