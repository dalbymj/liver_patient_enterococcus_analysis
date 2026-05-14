# Figure 5E: Faecal calprotectin by disease severity group
#
# Purpose:
# This script loads faecal calprotectin data and clinical metadata, generates
# the Figure 5E boxplot, and saves statistical test outputs.
#
# Input files:
# - data/clinical_metadata_and_data_renamed_30.01.2025.csv
# - data/FCAL_GLA_RifSys_combined_dataset_20.01.2022.csv
#
# Output files:
# - outputs/Figure_5/Fig5E_faecal_calprotectin_boxplot.pdf
# - outputs/Figure_5/Fig5E_faecal_calprotectin_plot_data.csv
# - outputs/Figure_5/Fig5E_faecal_calprotectin_summary_statistics.csv
# - outputs/Figure_5/Fig5E_faecal_calprotectin_statistics.txt
# - outputs/Figure_5/Fig5E_faecal_calprotectin_dunn_test.csv
#
# Run this script from the root folder of the repository:
# liver_patient_enterococcus_analysis/

library(tidyverse)

options(scipen = 999)

# Create output folder if it does not already exist
dir.create("outputs/Figure_5", recursive = TRUE, showWarnings = FALSE)

# Define input files
metadata_file <- "data/clinical_metadata_and_data_renamed_30.01.2025.csv"
fcal_file <- "data/FCAL_GLA_RifSys_combined_dataset_20.01.2022.csv"

# Check that input files exist
if (!file.exists(metadata_file)) {
  stop("Missing input file: ", metadata_file)
}

if (!file.exists(fcal_file)) {
  stop("Missing input file: ", fcal_file)
}

# Load data
metadata_raw <- read_csv(metadata_file, show_col_types = FALSE)
fcal_raw <- read_csv(fcal_file, show_col_types = FALSE)

# Remove unnamed index columns if present
metadata_raw <- metadata_raw %>%
  select(-any_of("...1"))

fcal_raw <- fcal_raw %>%
  select(-any_of("...1"))

# Rename columns if needed
if ("Study no" %in% names(fcal_raw)) {
  names(fcal_raw)[names(fcal_raw) == "Study no"] <- "study_no"
}

# Check required columns
required_fcal_cols <- c("study_no", "FCAL")
missing_fcal_cols <- setdiff(required_fcal_cols, names(fcal_raw))

if (length(missing_fcal_cols) > 0) {
  stop(
    "The following FCAL columns are missing: ",
    paste(missing_fcal_cols, collapse = ", ")
  )
}

# Prepare FCAL data
fcal_data <- fcal_raw %>%
  select(
    study_no,
    FCAL
  ) %>%
  mutate(
    FCAL = as.numeric(FCAL)
  ) %>%
  filter(
    complete.cases(.)
  )

# Prepare metadata
metadata_subsampled <- metadata_raw %>%
  select(
    study_no,
    Inclusion_criteria
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
    )
  )

# Merge metadata and FCAL data
plot_data <- metadata_subsampled %>%
  inner_join(fcal_data, by = "study_no") %>%
  filter(!is.na(Inclusion_criteria))

# Check that data are available
if (nrow(plot_data) == 0) {
  stop("No matched samples were found after merging metadata and FCAL data.")
}

# Save plot data
write_csv(
  plot_data,
  "outputs/Figure_5/Fig5E_faecal_calprotectin_plot_data.csv"
)

# Summary statistics
summary_statistics <- plot_data %>%
  group_by(Inclusion_criteria) %>%
  summarise(
    n = n(),
    median = median(FCAL, na.rm = TRUE),
    IQR = IQR(FCAL, na.rm = TRUE),
    mean = mean(FCAL, na.rm = TRUE),
    sd = sd(FCAL, na.rm = TRUE),
    min = min(FCAL, na.rm = TRUE),
    max = max(FCAL, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  summary_statistics,
  "outputs/Figure_5/Fig5E_faecal_calprotectin_summary_statistics.csv"
)

# Generate Figure 5E
fig5e <- ggplot(
  plot_data,
  aes(
    x = Inclusion_criteria,
    y = FCAL,
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
    position = position_jitter(width = 0.25),
    aes(
      colour = Inclusion_criteria,
      fill = Inclusion_criteria
    ),
    size = 1.5,
    stroke = 0.35
  ) +
  scale_fill_manual(
    values = c(
      "Healthy control" = "lightskyblue",
      "Stable cirrhotic" = "yellow",
      "Decomp cirrhotic" = "orange",
      "ACLF" = "firebrick1"
    )
  ) +
  scale_color_manual(
    values = c(
      "Healthy control" = "blue4",
      "Stable cirrhotic" = "yellow4",
      "Decomp cirrhotic" = "orange4",
      "ACLF" = "red4"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Healthy control" = 1,
      "Stable cirrhotic" = 1,
      "Decomp cirrhotic" = 1,
      "ACLF" = 1
    )
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 3000)
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
    y = "Faecal calprotectin (FCAL) (\u00B5g/g)"
  ) +
  theme(
    plot.title = element_text(size = 8),
    plot.margin = margin(t = 4, r = 4, b = 4, l = 4),
    legend.position = "none",
    legend.spacing.x = unit(0.2, "cm"),
    legend.spacing.y = unit(0.2, "cm"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.box.margin = margin(t = -8, r = -8, b = 0, l = -8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.background = element_blank(),
    legend.key = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 8, colour = "black"),
    axis.text.x = element_text(size = 7, colour = "black", angle = 0, hjust = 0.5),
    axis.title.y = element_text(size = 8, colour = "black"),
    axis.text.y = element_text(size = 8, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    aspect.ratio = 1
  ) +
  annotate("text", x = 2.5, y = 2900, label = "p < 0.001", size = 2.5) +
  annotate("segment", x = 1, xend = 4, y = 2750, yend = 2750, linewidth = 0.35) +
  annotate("text", x = 2, y = 2600, label = "p < 0.001", size = 2.5) +
  annotate("segment", x = 1, xend = 3, y = 2450, yend = 2450, linewidth = 0.35) +
  annotate("text", x = 1.5, y = 2300, label = "p = 0.01", size = 2.5) +
  annotate("segment", x = 1, xend = 2, y = 2150, yend = 2150, linewidth = 0.35)

# Save Figure 5E
ggsave(
  filename = "outputs/Figure_5/Fig5E_faecal_calprotectin_boxplot.pdf",
  plot = fig5e,
  width = 5.5,
  height = 5.5,
  units = "cm",
  dpi = 300
)

# Helper function for safe Shapiro-Wilk tests
safe_shapiro_p <- function(x) {
  
  x <- x[!is.na(x)]
  
  if (length(x) < 3 || length(unique(x)) < 2) {
    return(NA_real_)
  }
  
  shapiro.test(x)$p.value
}

# Statistical testing
anova_model <- aov(
  FCAL ~ Inclusion_criteria,
  data = plot_data
)

anova_summary <- summary(anova_model)
tukey_results <- TukeyHSD(anova_model)
anova_residuals <- residuals(anova_model)
residual_shapiro <- shapiro.test(anova_residuals)

kruskal_results <- kruskal.test(
  FCAL ~ Inclusion_criteria,
  data = plot_data
)

pairwise_wilcox_results <- pairwise.wilcox.test(
  plot_data$FCAL,
  plot_data$Inclusion_criteria,
  p.adjust.method = "fdr"
)

group_shapiro_results <- plot_data %>%
  group_by(Inclusion_criteria) %>%
  summarise(
    n = n(),
    shapiro_p_value = safe_shapiro_p(FCAL),
    .groups = "drop"
  )

# Levene test, if the car package is installed
levene_results <- NULL

if (requireNamespace("car", quietly = TRUE)) {
  levene_results <- car::leveneTest(
    FCAL ~ Inclusion_criteria,
    data = plot_data
  )
} else {
  warning("Package 'car' is not installed, so Levene's test was not run.")
}

# Dunn test, if the FSA package is installed
dunn_results <- NULL

if (requireNamespace("FSA", quietly = TRUE)) {
  
  dunn_results <- FSA::dunnTest(
    FCAL ~ Inclusion_criteria,
    data = plot_data,
    method = "bonferroni"
  )
  
  write_csv(
    dunn_results$res,
    "outputs/Figure_5/Fig5E_faecal_calprotectin_dunn_test.csv"
  )
  
} else {
  warning("Package 'FSA' is not installed, so Dunn's test was not run.")
}

# Save statistical output
statistics_output_file <- "outputs/Figure_5/Fig5E_faecal_calprotectin_statistics.txt"

sink(statistics_output_file)

cat("Figure 5E: Faecal calprotectin by disease severity group\n")
cat("========================================================\n\n")

cat("Input files:\n")
cat(metadata_file, "\n")
cat(fcal_file, "\n\n")

cat("Sample numbers by group:\n")
print(table(plot_data$Inclusion_criteria, useNA = "ifany"))

cat("\nSummary statistics by group:\n")
print(summary_statistics)

cat("\nGroup-level Shapiro-Wilk tests:\n")
print(group_shapiro_results)

cat("\nANOVA:\n")
print(anova_summary)

cat("\nTukey HSD:\n")
print(tukey_results)

cat("\nShapiro-Wilk test of ANOVA residuals:\n")
print(residual_shapiro)

cat("\nLevene test for homogeneity of variance:\n")
if (!is.null(levene_results)) {
  print(levene_results)
} else {
  cat("Levene test not run because package 'car' is not installed.\n")
}

cat("\nKruskal-Wallis test:\n")
print(kruskal_results)

cat("\nPairwise Wilcoxon tests with FDR correction:\n")
print(pairwise_wilcox_results)

cat("\nDunn test with Bonferroni correction:\n")
if (!is.null(dunn_results)) {
  print(dunn_results)
} else {
  cat("Dunn test not run because package 'FSA' is not installed.\n")
}

sink()

message("Figure saved to: outputs/Figure_5/Fig5E_faecal_calprotectin_boxplot.pdf")
message("Plot data saved to: outputs/Figure_5/Fig5E_faecal_calprotectin_plot_data.csv")
message("Summary statistics saved to: outputs/Figure_5/Fig5E_faecal_calprotectin_summary_statistics.csv")
message("Statistics saved to: ", statistics_output_file)