# Supplementary Figure 3F: All selected genus and urinary metabolite correlation heatmap
#
# Purpose:
# This script loads genus-level 16S rRNA sequencing data and urinary NMR
# metabolite data, selects the top 20 genera by total relative abundance,
# calculates Pearson correlations between those genera and selected urinary
# metabolites after bestNormalize transformation, applies Benjamini-Hochberg
# correction by genus, orders genera by their correlation with hippurate, and
# generates the Supplementary Figure 3F heatmap.
#
# Input files:
# - data/normalised_genus_data_16.04.2024.csv
# - data/Updated_urinary_NMR_GLA_RifSys_8.03.2021 with edited sample_ID.csv
#
# Output files:
# - outputs/Supplementary_Figures/Supplementary_figure_3/SuppFig3F_all_genus_urinary_metabolite_correlation_heatmap.pdf
# - outputs/Supplementary_Figures/Supplementary_figure_3/SuppFig3F_all_genus_urinary_metabolite_correlation_table.csv
# - outputs/Supplementary_Figures/Supplementary_figure_3/SuppFig3F_all_genus_urinary_metabolite_correlation_plot_data.csv
# - outputs/Supplementary_Figures/Supplementary_figure_3/SuppFig3F_bestNormalize_transformation_log.csv
# - outputs/Supplementary_Figures/Supplementary_figure_3/SuppFig3F_all_genus_urinary_metabolite_correlation_summary.txt
#
# Run this script from the root folder of the repository:
# liver_patient_enterococcus_analysis/

library(tidyverse)
library(bestNormalize)

options(scipen = 999)

set.seed(123)

# Create output folder
dir.create(
  "outputs/Supplementary_Figures/Supplementary_figure_3",
  recursive = TRUE,
  showWarnings = FALSE
)

# Define input files
genus_file <- "data/normalised_genus_data_16.04.2024.csv"
urinary_file <- "data/Urinary_NMR_GLA_RifSys_8.03.2021.csv"

# Check that input files exist
if (!file.exists(genus_file)) {
  stop("Missing input file: ", genus_file)
}

if (!file.exists(urinary_file)) {
  stop("Missing input file: ", urinary_file)
}

# Load data
genus_data_raw <- read_csv(genus_file, show_col_types = FALSE)
metabolite_data_raw <- read_csv(urinary_file, show_col_types = FALSE)

# Remove unnamed index columns if present
genus_data_raw <- genus_data_raw %>%
  select(-any_of("...1"))

metabolite_data_raw <- metabolite_data_raw %>%
  select(-any_of("...1"))

# Standardise urinary sample ID column name
if ("Study no" %in% names(metabolite_data_raw)) {
  names(metabolite_data_raw)[names(metabolite_data_raw) == "Study no"] <- "study_no"
}

# Replace "/" in metabolite names with "_"
names(metabolite_data_raw) <- gsub(
  x = names(metabolite_data_raw),
  pattern = "/",
  replacement = "_"
)

# Define genus columns
genus_start_col <- match("Mucispirillum", names(genus_data_raw))

if (is.na(genus_start_col)) {
  stop("Could not find genus column 'Mucispirillum' in genus_data_raw.")
}

genus_cols <- names(genus_data_raw)[genus_start_col:ncol(genus_data_raw)]

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

# Prepare genus data
genus_data <- genus_data_raw %>%
  select(study_no, all_of(genus_cols)) %>%
  mutate(
    across(
      all_of(genus_cols),
      ~ as.numeric(.x)
    )
  ) %>%
  drop_na()

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

# Prepare urinary metabolite data
metabolite_data <- metabolite_data_raw %>%
  select(study_no, all_of(metabolite_cols)) %>%
  mutate(
    across(
      all_of(metabolite_cols),
      ~ as.numeric(.x)
    )
  ) %>%
  drop_na(all_of(metabolite_cols))

# Merge genus and urinary metabolite data
merged_data <- genus_data_percent %>%
  inner_join(metabolite_data, by = "study_no")

# Prepare genus table
genus_table <- merged_data %>%
  select(study_no, all_of(genus_cols)) %>%
  column_to_rownames("study_no")

# Prepare urinary metabolite table
metabolite_table <- merged_data %>%
  select(study_no, all_of(metabolite_cols)) %>%
  column_to_rownames("study_no")

# Check sample order
if (!identical(rownames(genus_table), rownames(metabolite_table))) {
  stop("Sample order does not match between genus table and urinary metabolite table.")
}

# Select top 20 genera by total relative abundance
selected_genera <- genus_table %>%
  select(where(is.numeric)) %>%
  select(where(~ sum(.x, na.rm = TRUE) > 0)) %>%
  colSums(na.rm = TRUE) %>%
  sort(decreasing = TRUE) %>%
  names() %>%
  head(20)

# Urinary metabolites shown in Supplementary Figure 3F
selected_metabolites <- c(
  "Hip7p85_TS",
  "TMAO_TS",
  "Crn3p05_TS",
  "Cit_TS",
  "DMA_TS",
  "For_TS",
  "Man_TS"
)

metabolite_labels <- c(
  "Hip7p85_TS" = "Hippurate",
  "TMAO_TS" = "TMAO",
  "Crn3p05_TS" = "Creatinine",
  "Cit_TS" = "Citrate",
  "DMA_TS" = "Dimethylamine",
  "For_TS" = "Formate",
  "Man_TS" = "Mannitol"
)

# Check selected metabolites are present
missing_metabolites <- setdiff(selected_metabolites, names(metabolite_table))

if (length(missing_metabolites) > 0) {
  stop(
    "The following selected urinary metabolites are missing: ",
    paste(missing_metabolites, collapse = ", ")
  )
}

genera_for_correlation <- selected_genera

# Subset tables
genus_selected <- genus_table %>%
  select(all_of(genera_for_correlation))

metabolite_selected <- metabolite_table %>%
  select(all_of(selected_metabolites))

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

# Transform genus and urinary metabolite features
genus_transformed <- apply_best_normalize(
  data_table = genus_selected,
  data_type = "genus"
)

metabolite_transformed <- apply_best_normalize(
  data_table = metabolite_selected,
  data_type = "urinary_metabolite"
)

genus_for_correlation <- genus_transformed$data
metabolite_for_correlation <- metabolite_transformed$data

transformation_log <- bind_rows(
  genus_transformed$log,
  metabolite_transformed$log
)

write_csv(
  transformation_log,
  "outputs/Supplementary_Figures/Supplementary_figure_3/SuppFig3F_bestNormalize_transformation_log.csv"
)

# Correlation method
correlation_method <- "pearson"

# Calculate correlations
correlation_table <- expand_grid(
  Taxa = names(genus_for_correlation),
  Metabolite = names(metabolite_for_correlation)
) %>%
  mutate(
    correlation_result = map2(
      Taxa,
      Metabolite,
      function(taxon_name, metabolite_name) {
        
        correlation_data <- tibble(
          genus_value = genus_for_correlation[[taxon_name]],
          metabolite_value = metabolite_for_correlation[[metabolite_name]]
        ) %>%
          drop_na()
        
        if (
          nrow(correlation_data) < 3 ||
          sd(correlation_data$genus_value) == 0 ||
          sd(correlation_data$metabolite_value) == 0
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
          correlation_data$metabolite_value,
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

write_csv(
  correlation_table,
  "outputs/Supplementary_Figures/Supplementary_figure_3/SuppFig3F_all_genus_urinary_metabolite_correlation_table.csv"
)

# Order genera by correlation with hippurate
reference_metabolite <- "Hip7p85_TS"

genus_order_by_reference_metabolite <- correlation_table %>%
  filter(
    Metabolite == reference_metabolite,
    Taxa %in% genera_for_correlation
  ) %>%
  arrange(Correlation) %>%
  pull(Taxa)

if (length(genus_order_by_reference_metabolite) == 0) {
  stop("No genus order could be created. Check that Hip7p85_TS is present.")
}

# Prepare plot data
plot_data <- correlation_table %>%
  filter(
    Taxa %in% genera_for_correlation,
    Metabolite %in% selected_metabolites
  ) %>%
  mutate(
    Metabolite_label = dplyr::recode(Metabolite, !!!metabolite_labels),
    Metabolite_label = factor(
      Metabolite_label,
      levels = unname(metabolite_labels[selected_metabolites])
    ),
    Taxa = factor(
      Taxa,
      levels = genus_order_by_reference_metabolite
    )
  )

write_csv(
  plot_data,
  "outputs/Supplementary_Figures/Supplementary_figure_3/SuppFig3F_all_genus_urinary_metabolite_correlation_plot_data.csv"
)

# Generate Supplementary Figure 3F heatmap
suppfig3f <- ggplot(
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
    size = 2.5,
    vjust = 0.75,
    hjust = 0.5
  ) +
  labs(
    y = NULL,
    x = NULL,
    fill = correlation_method
  ) +
  facet_grid(
    . ~ Metabolite_label,
    drop = TRUE,
    scales = "free",
    space = "free_x"
  ) +
  theme(
    legend.position = c(-0.15, 1.25),
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

ggsave(
  filename = "outputs/Supplementary_Figures/Supplementary_figure_3/SuppFig3F_all_genus_urinary_metabolite_correlation_heatmap.pdf",
  plot = suppfig3f,
  width = 5.5,
  height = 6,
  units = "cm",
  dpi = 300
)

# Save summary
summary_output_file <- "outputs/Supplementary_Figures/Supplementary_figure_3/SuppFig3F_all_genus_urinary_metabolite_correlation_summary.txt"

sink(summary_output_file)

cat("Supplementary Figure 3F: All selected genus and urinary metabolite correlation heatmap\n")
cat("====================================================================================\n\n")

cat("Input files:\n")
cat(genus_file, "\n")
cat(urinary_file, "\n\n")

cat("Number of matched samples:\n")
print(nrow(merged_data))

cat("\nCorrelation method:\n")
print(correlation_method)

cat("\nMultiple testing correction:\n")
cat("Benjamini-Hochberg correction applied separately within each genus.\n\n")

cat("Top 20 genera by total relative abundance:\n")
print(selected_genera)

cat("\nReference metabolite used for genus ordering:\n")
print(reference_metabolite)

cat("\nGenus order used in Supplementary Figure 3F heatmap:\n")
print(genus_order_by_reference_metabolite)

cat("\nUrinary metabolites shown in Supplementary Figure 3F:\n")
print(selected_metabolites)

cat("\nTransformation log:\n")
print(transformation_log)

cat("\nSignificant correlations shown in Supplementary Figure 3F:\n")
print(
  plot_data %>%
    filter(Significance != "") %>%
    select(Taxa, Metabolite, Metabolite_label, Correlation, Pvalue, AdjPvalue, Significance) %>%
    arrange(Taxa, Metabolite)
)

sink()

message("Figure saved to: outputs/Supplementary_Figures/Supplementary_figure_3/SuppFig3F_all_genus_urinary_metabolite_correlation_heatmap.pdf")
message("Correlation table saved to: outputs/Supplementary_Figures/Supplementary_figure_3/SuppFig3F_all_genus_urinary_metabolite_correlation_table.csv")
message("Plot data saved to: outputs/Supplementary_Figures/Supplementary_figure_3/SuppFig3F_all_genus_urinary_metabolite_correlation_plot_data.csv")
message("Summary saved to: ", summary_output_file)