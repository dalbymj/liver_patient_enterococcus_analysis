# Figure 3I: Clinical cholestasis marker and bile acid correlation heatmap
#
# Purpose:
# This script loads clinical metadata and plasma bile acid data, applies
# bestNormalize transformation to selected clinical cholestasis markers and bile
# acids, calculates Pearson correlations, and generates the Figure 3I heatmap.
#
# Input files:
# - data/clinical_metadata_and_data_renamed_30.01.2025.csv
# - data/Corrected bile acid data 01.12.2020.csv
#
# Output files:
# - outputs/Figure_3/Fig3I_cholestasis_bile_acid_correlation_heatmap.pdf
# - outputs/Figure_3/Fig3I_cholestasis_bile_acid_correlation_table.csv
# - outputs/Figure_3/Fig3I_cholestasis_bile_acid_correlation_plot_data.csv
# - outputs/Figure_3/Fig3I_bestNormalize_transformation_log.csv
# - outputs/Figure_3/Fig3I_cholestasis_bile_acid_correlation_summary.txt
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

# Clinical cholestasis markers used in Figure 3I
clinical_marker_cols <- c(
  "ALP",
  "Total_Br_SI_units"
)

clinical_marker_labels <- c(
  "ALP" = "ALP",
  "Total_Bilirubin" = "Total bilirubin"
)

# Bile acids used in Figure 3I
selected_bile_acids <- c(
  "Glycocholic_acid",
  "Taurocholic_acid",
  "Glycochenodeoxycholic_acid",
  "Taurochenodeoxycholic_acid"
)

bile_acid_labels <- c(
  "Glycocholic_acid" = "Glycocholic",
  "Taurocholic_acid" = "Taurocholic",
  "Glycochenodeoxycholic_acid" = "Glycochenodeoxycholic",
  "Taurochenodeoxycholic_acid" = "Taurochenodeoxycholic"
)

# Check required columns are present
missing_metadata_cols <- setdiff(c("study_no", clinical_marker_cols), names(metadata_raw))
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

# Prepare clinical marker data
clinical_data <- metadata_raw %>%
  select(
    study_no,
    ALP,
    Total_Br_SI_units
  ) %>%
  rename(
    Total_Bilirubin = Total_Br_SI_units
  ) %>%
  mutate(
    across(
      c(ALP, Total_Bilirubin),
      ~ as.numeric(.x)
    )
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

# Merge data and retain complete cases for the selected analysis variables
merged_data <- clinical_data %>%
  inner_join(bile_data, by = "study_no") %>%
  filter(
    complete.cases(
      across(c(ALP, Total_Bilirubin, all_of(selected_bile_acids)))
    )
  )

# Prepare tables
clinical_table <- merged_data %>%
  select(study_no, ALP, Total_Bilirubin) %>%
  column_to_rownames("study_no")

bile_table <- merged_data %>%
  select(study_no, all_of(selected_bile_acids)) %>%
  column_to_rownames("study_no")

# Check that sample order matches
if (!identical(rownames(clinical_table), rownames(bile_table))) {
  stop("Sample order does not match between clinical marker and bile acid tables.")
}

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

# Transform clinical markers and bile acids
clinical_transformed <- apply_best_normalize(
  data_table = clinical_table,
  data_type = "clinical_cholestasis_marker"
)

bile_transformed <- apply_best_normalize(
  data_table = bile_table,
  data_type = "bile_acid"
)

clinical_for_correlation <- clinical_transformed$data
bile_for_correlation <- bile_transformed$data

transformation_log <- bind_rows(
  clinical_transformed$log,
  bile_transformed$log
)

write_csv(
  transformation_log,
  "outputs/Figure_3/Fig3I_bestNormalize_transformation_log.csv"
)

# Correlation method
correlation_method <- "pearson"

# Calculate correlations
correlation_table <- expand_grid(
  Marker = names(clinical_for_correlation),
  Bile_acid = names(bile_for_correlation)
) %>%
  mutate(
    correlation_result = map2(
      Marker,
      Bile_acid,
      function(marker_name, bile_acid_name) {
        
        correlation_data <- tibble(
          marker_value = clinical_for_correlation[[marker_name]],
          bile_acid_value = bile_for_correlation[[bile_acid_name]]
        ) %>%
          drop_na()
        
        if (
          nrow(correlation_data) < 3 ||
          sd(correlation_data$marker_value) == 0 ||
          sd(correlation_data$bile_acid_value) == 0
        ) {
          return(
            tibble(
              Correlation = NA_real_,
              Pvalue = NA_real_,
              n = nrow(correlation_data)
            )
          )
        }
        
        correlation_test <- cor.test(
          correlation_data$marker_value,
          correlation_data$bile_acid_value,
          method = correlation_method
        )
        
        tibble(
          Correlation = as.numeric(correlation_test$estimate),
          Pvalue = correlation_test$p.value,
          n = nrow(correlation_data)
        )
      }
    )
  ) %>%
  unnest(correlation_result) %>%
  filter(!is.na(Correlation), !is.na(Pvalue)) %>%
  mutate(
    # Original script used no multiple-testing correction for this heatmap.
    AdjPvalue = Pvalue,
    Significance = cut(
      AdjPvalue,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("***", "**", "*", "")
    ),
    Significance = as.character(Significance),
    Significance = replace_na(Significance, "")
  )

# Save full correlation table
write_csv(
  correlation_table,
  "outputs/Figure_3/Fig3I_cholestasis_bile_acid_correlation_table.csv"
)

# Prepare plot data
plot_data <- correlation_table %>%
  mutate(
    Marker_label = dplyr::recode(
      Marker,
      !!!clinical_marker_labels
    ),
    Marker_label = factor(
      Marker_label,
      levels = unname(clinical_marker_labels[c("ALP", "Total_Bilirubin")])
    ),
    Bile_acid_label = dplyr::recode(
      Bile_acid,
      !!!bile_acid_labels
    ),
    Bile_acid_label = factor(
      Bile_acid_label,
      levels = unname(bile_acid_labels[selected_bile_acids])
    )
  )

write_csv(
  plot_data,
  "outputs/Figure_3/Fig3I_cholestasis_bile_acid_correlation_plot_data.csv"
)

# Generate Figure 3I heatmap
fig3i <- ggplot(
  data = plot_data,
  aes(
    x = Bile_acid_label,
    y = Marker_label,
    fill = Correlation
  )
) +
  geom_tile() +
  scale_fill_gradient2(
    low = "slateblue",
    mid = "white",
    high = "red"
  ) +
  geom_text(
    aes(label = Significance),
    color = "black",
    size = 4,
    vjust = 0.75,
    hjust = 0.5
  ) +
  labs(
    y = NULL,
    x = NULL,
    fill = correlation_method
  ) +
  facet_grid(
    . ~ Bile_acid_label,
    drop = TRUE,
    scales = "free",
    space = "free_x"
  ) +
  theme(
    plot.title = element_text(size = 5),
    legend.key.size = unit(0.3, "cm"),
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
    panel.spacing = unit(0.05, "lines"),
    legend.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(size = 7, colour = "black", angle = 90, vjust = 0.5, hjust = 0),
    axis.title.x = element_text(size = 8, colour = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 8, colour = "black"),
    axis.text.y = element_text(size = 7, colour = "black")
  )

# Save Figure 3I
ggsave(
  filename = "outputs/Figure_3/Fig3I_cholestasis_bile_acid_correlation_heatmap.pdf",
  plot = fig3i,
  width = 7,
  height = 5,
  units = "cm",
  dpi = 300
)

# Save summary output
summary_output_file <- "outputs/Figure_3/Fig3I_cholestasis_bile_acid_correlation_summary.txt"

sink(summary_output_file)

cat("Figure 3I: Clinical cholestasis marker and bile acid correlation heatmap\n")
cat("========================================================================\n\n")

cat("Input files:\n")
cat(metadata_file, "\n")
cat(bile_file, "\n\n")

cat("Number of matched complete-case samples:\n")
print(nrow(merged_data))

cat("\nClinical markers:\n")
print(names(clinical_table))

cat("\nBile acids:\n")
print(selected_bile_acids)

cat("\nCorrelation method:\n")
print(correlation_method)

cat("\nMultiple testing correction:\n")
cat("None. AdjPvalue is set equal to Pvalue to match the original script.\n\n")

cat("Transformation log:\n")
print(transformation_log)

cat("\nCorrelation table:\n")
print(correlation_table)

sink()

message("Figure saved to: outputs/Figure_3/Fig3I_cholestasis_bile_acid_correlation_heatmap.pdf")
message("Correlation table saved to: outputs/Figure_3/Fig3I_cholestasis_bile_acid_correlation_table.csv")
message("Plot data saved to: outputs/Figure_3/Fig3I_cholestasis_bile_acid_correlation_plot_data.csv")
message("Summary saved to: ", summary_output_file)