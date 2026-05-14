# Figure 5I: Monocyte FACS and microbiome correlation heatmap
#
# Purpose:
# This script loads monocyte FACS data, genus-level 16S rRNA sequencing data,
# and whole blood bacterial 16S rDNA qPCR data. It applies bestNormalize
# transformation, calculates Pearson correlations between monocyte response
# variables and selected microbiome/translocation variables, applies
# Benjamini-Hochberg correction within each FACS variable, and generates
# the Figure 5I correlation heatmap.
#
# Input files:
# - data/FACS_data_GLA TLR249_stimulation_07.09.2022.csv
# - data/normalised_genus_data_16.04.2024.csv
# - data/GLA_RifSys_16S_rRNA_QPCR_07.09.2020.csv
#
# Output files:
# - outputs/Figure_5/Fig5I_monocyte_facs_microbiome_correlation_heatmap.pdf
# - outputs/Figure_5/Fig5I_monocyte_facs_microbiome_correlation_table.csv
# - outputs/Figure_5/Fig5I_monocyte_facs_microbiome_correlation_plot_data.csv
# - outputs/Figure_5/Fig5I_bestNormalize_transformation_log.csv
# - outputs/Figure_5/Fig5I_monocyte_facs_microbiome_correlation_summary.txt
#
# Run this script from the root folder of the repository:
# liver_patient_enterococcus_analysis/

library(tidyverse)
library(vegan)
library(bestNormalize)
library(grid)

options(scipen = 999)

set.seed(123)

# Create output folder if it does not already exist
dir.create("outputs/Figure_5", recursive = TRUE, showWarnings = FALSE)

# Define input files
facs_file <- "data/FACS_data_GLA TLR249_stimulation_07.09.2022.csv"
genus_file <- "data/normalised_genus_data_16.04.2024.csv"
blood_16s_file <- "data/GLA_RifSys_16S_rRNA_QPCR_07.09.2020.csv"

# Check that input files exist
input_files <- c(facs_file, genus_file, blood_16s_file)

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
  "us_tnfa",
  "us_hladr",
  "us_il6",
  "lps_tnfa",
  "lps_hladr",
  "lps_il6",
  "pam3_tnfa",
  "pam3_hladr",
  "pam3_il6"
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

required_genus_cols <- c("Enterococcus", "Lachnospira")
missing_required_genera <- setdiff(required_genus_cols, genus_cols)

if (length(missing_required_genera) > 0) {
  stop(
    "The following required genera are missing: ",
    paste(missing_required_genera, collapse = ", ")
  )
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
  select(-row_total)

# Calculate number of detected genera from untransformed relative abundance data
genus_number_table <- genus_percent %>%
  column_to_rownames("study_no")

genus_number <- vegan::specnumber(
  genus_number_table,
  MARGIN = 1
)

genus_number <- tibble(
  study_no = names(genus_number),
  genus_number = as.numeric(genus_number)
)

selected_genus_data <- genus_percent %>%
  select(
    study_no,
    Enterococcus,
    Lachnospira
  ) %>%
  inner_join(
    genus_number,
    by = "study_no"
  )

# Load blood bacterial 16S rDNA qPCR data
blood_16s_raw <- read_csv(blood_16s_file, show_col_types = FALSE) %>%
  select(-any_of("...1"))

if ("Study no" %in% names(blood_16s_raw)) {
  names(blood_16s_raw)[names(blood_16s_raw) == "Study no"] <- "study_no"
}

if ("16S Bac DNA" %in% names(blood_16s_raw)) {
  names(blood_16s_raw)[names(blood_16s_raw) == "16S Bac DNA"] <- "Bac_16S_DNA"
}

required_blood_16s_cols <- c("study_no", "Bac_16S_DNA")
missing_blood_16s_cols <- setdiff(required_blood_16s_cols, names(blood_16s_raw))

if (length(missing_blood_16s_cols) > 0) {
  stop(
    "The following blood 16S columns are missing: ",
    paste(missing_blood_16s_cols, collapse = ", ")
  )
}

blood_16s_data <- blood_16s_raw %>%
  select(
    study_no,
    Bac_16S_DNA
  ) %>%
  mutate(
    Bac_16S_DNA = as.numeric(Bac_16S_DNA)
  ) %>%
  filter(complete.cases(.))

# Merge data
microbiome_data <- selected_genus_data %>%
  left_join(
    blood_16s_data,
    by = "study_no"
  )

merged_data <- facs_data %>%
  left_join(
    microbiome_data,
    by = "study_no"
  )

if (nrow(merged_data) == 0) {
  stop("No matched samples were found after merging FACS and microbiome data.")
}

microbiome_variables <- c(
  "Bac_16S_DNA",
  "genus_number",
  "Enterococcus",
  "Lachnospira"
)

# Transform FACS and microbiome variables
facs_transformed_output <- apply_best_normalize(
  data_table = merged_data,
  variables = facs_variables,
  data_type = "FACS"
)

microbiome_transformed_output <- apply_best_normalize(
  data_table = facs_transformed_output$data,
  variables = microbiome_variables,
  data_type = "microbiome"
)

analysis_data <- microbiome_transformed_output$data

transformation_log <- bind_rows(
  facs_transformed_output$log,
  microbiome_transformed_output$log
)

write_csv(
  transformation_log,
  "outputs/Figure_5/Fig5I_bestNormalize_transformation_log.csv"
)

# Calculate Pearson correlations
correlation_method <- "pearson"

correlation_table <- expand_grid(
  FACS_marker = facs_variables,
  Microbiome_variable = microbiome_variables
) %>%
  mutate(
    correlation_result = map2(
      FACS_marker,
      Microbiome_variable,
      function(facs_marker, microbiome_variable) {
        
        correlation_data <- tibble(
          facs_value = analysis_data[[facs_marker]],
          microbiome_value = analysis_data[[microbiome_variable]]
        ) %>%
          drop_na()
        
        if (
          nrow(correlation_data) < 3 ||
          sd(correlation_data$facs_value) == 0 ||
          sd(correlation_data$microbiome_value) == 0
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
          correlation_data$facs_value,
          correlation_data$microbiome_value,
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
  group_by(FACS_marker) %>%
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
  "outputs/Figure_5/Fig5I_monocyte_facs_microbiome_correlation_table.csv"
)

# Prepare plot labels and ordering
microbiome_label_levels <- c(
  "Bacterial\n16S rDNA",
  "No. of\ngenera",
  "Enterococcus",
  "Lachnospira"
)

facs_label_levels <- c(
  "PAM3:               HLA-DR",
  "PAM3:               IL-6",
  "PAM3:               TNF-alpha",
  "LPS:                  HLA-DR",
  "LPS:                  IL-6",
  "LPS:                  TNF-alpha",
  "Unstimulated:  HLA-DR",
  "Unstimulated:  IL-6",
  "Unstimulated:  TNF-alpha"
)

plot_data <- correlation_table %>%
  mutate(
    Microbiome_label = dplyr::recode(
      Microbiome_variable,
      "Bac_16S_DNA" = "Bacterial\n16S rDNA",
      "genus_number" = "No. of\ngenera",
      "Enterococcus" = "Enterococcus",
      "Lachnospira" = "Lachnospira"
    ),
    Microbiome_label = factor(
      Microbiome_label,
      levels = microbiome_label_levels
    ),
    FACS_label = dplyr::recode(
      FACS_marker,
      "pam3_hladr" = "PAM3:               HLA-DR",
      "pam3_il6" = "PAM3:               IL-6",
      "pam3_tnfa" = "PAM3:               TNF-alpha",
      "lps_hladr" = "LPS:                  HLA-DR",
      "lps_il6" = "LPS:                  IL-6",
      "lps_tnfa" = "LPS:                  TNF-alpha",
      "us_hladr" = "Unstimulated:  HLA-DR",
      "us_il6" = "Unstimulated:  IL-6",
      "us_tnfa" = "Unstimulated:  TNF-alpha"
    ),
    FACS_label = factor(
      FACS_label,
      levels = facs_label_levels
    ),
    panel_x = ""
  )

write_csv(
  plot_data,
  "outputs/Figure_5/Fig5I_monocyte_facs_microbiome_correlation_plot_data.csv"
)

# Generate Figure 5I heatmap
fig5i <- ggplot(
  data = plot_data,
  aes(
    x = panel_x,
    y = FACS_label,
    fill = Correlation
  )
) +
  geom_tile() +
  scale_fill_gradient2(
    "Pearson",
    low = "slateblue",
    mid = "white",
    high = "red"
  ) +
  geom_text(
    aes(label = Significance),
    color = "black",
    size = 5,
    hjust = 0.5,
    vjust = 0.8
  ) +
  labs(
    y = "",
    x = "",
    fill = correlation_method
  ) +
  facet_grid(
    . ~ Microbiome_label,
    drop = TRUE,
    scales = "free",
    space = "free_x"
  ) +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"),
    legend.margin = margin(-10, 0, 0, 0),
    plot.title = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    panel.spacing = unit(0.25, "lines"),
    legend.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(size = 7, colour = "black"),
    axis.title.x = element_text(size = 8, colour = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 8, colour = "black"),
    axis.text.y = element_text(size = 8, colour = "black", hjust = 0)
  )

ggsave(
  filename = "outputs/Figure_5/Fig5I_monocyte_facs_microbiome_correlation_heatmap.pdf",
  plot = fig5i,
  height = 7,
  width = 12.5,
  units = "cm",
  dpi = 300
)

# Save summary output
summary_output_file <- "outputs/Figure_5/Fig5I_monocyte_facs_microbiome_correlation_summary.txt"

sink(summary_output_file)

cat("Figure 5I: Monocyte FACS and microbiome correlation heatmap\n")
cat("===========================================================\n\n")

cat("Input files:\n")
cat(facs_file, "\n")
cat(genus_file, "\n")
cat(blood_16s_file, "\n\n")

cat("Number of FACS samples:\n")
print(nrow(facs_data))

cat("\nNumber of merged samples:\n")
print(nrow(merged_data))

cat("\nCorrelation method:\n")
print(correlation_method)

cat("\nMultiple testing correction:\n")
cat("Benjamini-Hochberg correction applied separately within each FACS variable.\n\n")

cat("FACS variables:\n")
print(facs_variables)

cat("\nMicrobiome variables:\n")
print(microbiome_variables)

cat("\nTransformation log:\n")
print(transformation_log)

cat("\nSignificant correlations:\n")
print(
  plot_data %>%
    filter(Significance != "") %>%
    select(
      FACS_marker,
      Microbiome_variable,
      Correlation,
      Pvalue,
      AdjPvalue,
      Significance,
      n
    ) %>%
    arrange(FACS_marker, Microbiome_variable)
)

cat("\nNote on bestNormalize:\n")
cat("bestNormalize was run with out_of_sample = FALSE to avoid k-fold cross-validation errors in variables with small matched sample numbers.\n")

sink()

message("Figure saved to: outputs/Figure_5/Fig5I_monocyte_facs_microbiome_correlation_heatmap.pdf")
message("Correlation table saved to: outputs/Figure_5/Fig5I_monocyte_facs_microbiome_correlation_table.csv")
message("Plot data saved to: outputs/Figure_5/Fig5I_monocyte_facs_microbiome_correlation_plot_data.csv")
message("Summary saved to: ", summary_output_file)