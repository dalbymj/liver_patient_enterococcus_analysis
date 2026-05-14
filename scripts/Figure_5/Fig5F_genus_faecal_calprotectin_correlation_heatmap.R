# Figure 5F: Genus and faecal calprotectin correlation heatmap
#
# Purpose:
# This script loads genus-level 16S rRNA sequencing data and faecal
# calprotectin data, selects the top 20 genera by total relative abundance,
# calculates Pearson correlations after bestNormalize transformation, and
# generates the smaller Figure 5F correlation heatmap.
#
# Input files:
# - data/normalised_genus_data_16.04.2024.csv
# - data/data/FCAL_GLA_RifSys_combined_dataset_20.01.2022.csv
#
# Output files:
# - outputs/Figure_5/Fig5F_genus_faecal_calprotectin_correlation_heatmap.pdf
# - outputs/Figure_5/Fig5F_genus_faecal_calprotectin_correlation_table.csv
# - outputs/Figure_5/Fig5F_genus_faecal_calprotectin_correlation_plot_data.csv
# - outputs/Figure_5/Fig5F_bestNormalize_transformation_log.csv
# - outputs/Figure_5/Fig5F_genus_faecal_calprotectin_correlation_summary.txt
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

# Merge genus and FCAL data
merged_data <- fcal_data %>%
  inner_join(genus_data_percent, by = "study_no")

if (nrow(merged_data) == 0) {
  stop("No matched samples were found after merging genus and FCAL data.")
}

# Prepare genus table
genus_table <- merged_data %>%
  select(study_no, all_of(genus_cols)) %>%
  column_to_rownames("study_no")

# Prepare FCAL table
fcal_table <- merged_data %>%
  select(study_no, FCAL) %>%
  column_to_rownames("study_no")

# Check sample order
if (!identical(rownames(genus_table), rownames(fcal_table))) {
  stop("Sample order does not match between genus table and FCAL table.")
}

# Select top 20 genera by total relative abundance
top_genera <- genus_table %>%
  select(where(is.numeric)) %>%
  select(where(~ sum(.x, na.rm = TRUE) > 0)) %>%
  colSums(na.rm = TRUE) %>%
  sort(decreasing = TRUE) %>%
  names() %>%
  head(20)

# Smaller Figure 5F heatmap genus set, in manuscript order
figure_5f_genera_order <- c(
  "Lachnospira",
  "Roseburia",
  "Faecalibacterium",
  "Enterococcus"
)

missing_figure_5f_genera <- setdiff(figure_5f_genera_order, names(genus_table))

if (length(missing_figure_5f_genera) > 0) {
  warning(
    "The following Figure 5F genera are missing from the input data and will be omitted: ",
    paste(missing_figure_5f_genera, collapse = ", ")
  )
}

figure_5f_genera_not_in_top20 <- setdiff(figure_5f_genera_order, top_genera)

if (length(figure_5f_genera_not_in_top20) > 0) {
  warning(
    "The following Figure 5F genera are not in the top 20 genera by total relative abundance and will be omitted: ",
    paste(figure_5f_genera_not_in_top20, collapse = ", ")
  )
}

figure_5f_genera <- figure_5f_genera_order[
  figure_5f_genera_order %in% top_genera &
    figure_5f_genera_order %in% names(genus_table)
]

if (length(figure_5f_genera) == 0) {
  stop("No Figure 5F genera remain after filtering to top 20 genera present in the input data.")
}

# Subset to top 20 genera for the correlation analysis
genus_selected <- genus_table %>%
  select(all_of(top_genera))

fcal_selected <- fcal_table %>%
  select(FCAL)

# Function to apply bestNormalize to each column.
# If bestNormalize fails, fall back to simple z-score scaling.
apply_best_normalize <- function(data_table, data_type) {
  
  transformed_table <- data_table
  
  transformation_log <- tibble(
    data_type = character(),
    feature = character(),
    transformation = character(),
    note = character()
  )
  
  for (feature_name in names(data_table)) {
    
    message("Transforming: ", feature_name)
    
    feature_values <- data_table[[feature_name]]
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

# Transform data
genus_transformed <- apply_best_normalize(
  data_table = genus_selected,
  data_type = "genus"
)

fcal_transformed <- apply_best_normalize(
  data_table = fcal_selected,
  data_type = "faecal_calprotectin"
)

genus_for_correlation <- genus_transformed$data
fcal_for_correlation <- fcal_transformed$data

transformation_log <- bind_rows(
  genus_transformed$log,
  fcal_transformed$log
)

write_csv(
  transformation_log,
  "outputs/Figure_5/Fig5F_bestNormalize_transformation_log.csv"
)

# Correlation method
correlation_method <- "pearson"

# Calculate correlations
correlation_table <- expand_grid(
  Taxa = names(genus_for_correlation),
  Metabolite = names(fcal_for_correlation)
) %>%
  mutate(
    correlation_result = map2(
      Taxa,
      Metabolite,
      function(taxon_name, metabolite_name) {
        
        correlation_data <- tibble(
          genus_value = genus_for_correlation[[taxon_name]],
          fcal_value = fcal_for_correlation[[metabolite_name]]
        ) %>%
          drop_na()
        
        if (
          nrow(correlation_data) < 3 ||
          sd(correlation_data$genus_value) == 0 ||
          sd(correlation_data$fcal_value) == 0
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
          correlation_data$genus_value,
          correlation_data$fcal_value,
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
  group_by(Taxa) %>%
  mutate(
    AdjPvalue = p.adjust(Pvalue, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    Significance = cut(
      AdjPvalue,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("***", "**", "*", "")
    ),
    Significance = as.character(Significance),
    Significance = replace_na(Significance, "")
  )

# Save full correlation table for top 20 genera
write_csv(
  correlation_table,
  "outputs/Figure_5/Fig5F_genus_faecal_calprotectin_correlation_table.csv"
)

# Prepare smaller Figure 5F plot data
plot_data <- correlation_table %>%
  filter(
    Taxa %in% figure_5f_genera,
    Metabolite == "FCAL"
  ) %>%
  mutate(
    Metabolite_label = "FCAL",
    Metabolite_label = factor(
      Metabolite_label,
      levels = "FCAL"
    ),
    Taxa = factor(
      Taxa,
      levels = figure_5f_genera_order
    )
  ) %>%
  filter(!is.na(Taxa))

if (nrow(plot_data) == 0) {
  stop("No plot data were generated. Check genus and FCAL column names.")
}

write_csv(
  plot_data,
  "outputs/Figure_5/Fig5F_genus_faecal_calprotectin_correlation_plot_data.csv"
)

# Generate Figure 5F heatmap
fig5f <- ggplot(
  data = plot_data,
  aes(
    x = Metabolite_label,
    y = Taxa,
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
    size = 3
  ) +
  labs(
    y = NULL,
    x = NULL,
    fill = correlation_method
  ) +
  facet_grid(
    . ~ Metabolite_label,
    drop = TRUE
  ) +
  theme(
    legend.position = c(-1.5, 1.3),
    legend.direction = "horizontal",
    plot.title = element_text(size = 8),
    legend.key.size = unit(0.3, "cm"),
    legend.title = element_blank(),
    legend.text = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 0.5),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
    panel.spacing = unit(0.1, "lines"),
    legend.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(size = 8, colour = "black", angle = 90, vjust = 0.5, hjust = 0.5),
    axis.title.x = element_text(size = 8, colour = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 8, colour = "black"),
    axis.text.y = element_text(size = 8, colour = "black")
  )

# Save Figure 5F
ggsave(
  filename = "outputs/Figure_5/Fig5F_genus_faecal_calprotectin_correlation_heatmap.pdf",
  plot = fig5f,
  height = 4,
  width = 3.75,
  units = "cm",
  dpi = 300
)

# Save summary output
summary_output_file <- "outputs/Figure_5/Fig5F_genus_faecal_calprotectin_correlation_summary.txt"

sink(summary_output_file)

cat("Figure 5F: Genus and faecal calprotectin correlation heatmap\n")
cat("=============================================================\n\n")

cat("Input files:\n")
cat(genus_file, "\n")
cat(fcal_file, "\n\n")

cat("Number of matched samples:\n")
print(nrow(merged_data))

cat("\nCorrelation method:\n")
print(correlation_method)

cat("\nMultiple testing correction:\n")
cat("Benjamini-Hochberg correction applied separately within each genus.\n\n")

cat("Top 20 genera used in correlation analysis:\n")
print(top_genera)

cat("\nGenera shown in Figure 5F:\n")
print(figure_5f_genera)

cat("\nTransformation log:\n")
print(transformation_log)

cat("\nCorrelations shown in Figure 5F:\n")
print(
  plot_data %>%
    select(Taxa, Metabolite, Correlation, Pvalue, AdjPvalue, Significance, n)
)

sink()

message("Figure saved to: outputs/Figure_5/Fig5F_genus_faecal_calprotectin_correlation_heatmap.pdf")
message("Correlation table saved to: outputs/Figure_5/Fig5F_genus_faecal_calprotectin_correlation_table.csv")
message("Plot data saved to: outputs/Figure_5/Fig5F_genus_faecal_calprotectin_correlation_plot_data.csv")
message("Summary saved to: ", summary_output_file)