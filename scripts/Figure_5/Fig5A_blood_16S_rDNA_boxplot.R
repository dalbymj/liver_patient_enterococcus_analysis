# Figure 5A: Whole blood bacterial 16S rDNA by disease severity group
#
# Purpose:
# This script loads whole blood bacterial 16S rDNA qPCR data and clinical
# metadata, generates the Figure 5A boxplot, and saves statistical test outputs.
#
# Input files:
# - data/clinical_metadata_and_data_renamed_30.01.2025.csv
# - data/GLA_RifSys_16S_rRNA_QPCR_07.09.2020.csv
#
# Output files:
# - outputs/Figure_5/Fig5A_blood_16S_rDNA_boxplot.pdf
# - outputs/Figure_5/Fig5A_blood_16S_rDNA_plot_data.csv
# - outputs/Figure_5/Fig5A_blood_16S_rDNA_summary_statistics.csv
# - outputs/Figure_5/Fig5A_blood_16S_rDNA_statistics.txt
#
# Run this script from the root folder of the repository:
# liver_patient_enterococcus_analysis/

library(tidyverse)

options(scipen = 999)

# Create output folder if it does not already exist
dir.create("outputs/Figure_5", recursive = TRUE, showWarnings = FALSE)

# Define input files
metadata_file <- "data/clinical_metadata_and_data_renamed_30.01.2025.csv"
blood_16s_file <- "data/GLA_RifSys_16S_rRNA_QPCR_07.09.2020.csv"

# Check that input files exist
if (!file.exists(metadata_file)) {
  stop("Missing input file: ", metadata_file)
}

if (!file.exists(blood_16s_file)) {
  stop("Missing input file: ", blood_16s_file)
}

# Load data
metadata_raw <- read_csv(metadata_file, show_col_types = FALSE)
blood_16s_raw <- read_csv(blood_16s_file, show_col_types = FALSE)

# Remove unnamed index columns if present
metadata_raw <- metadata_raw %>%
  select(-any_of("...1"))

blood_16s_raw <- blood_16s_raw %>%
  select(-any_of("...1"))

# Rename qPCR columns
blood_16s_data <- blood_16s_raw %>%
  rename(
    study_no = `Study no`,
    Bac_16S_DNA = `16S Bac DNA`
  ) %>%
  select(
    study_no,
    Bac_16S_DNA
  ) %>%
  mutate(
    Bac_16S_DNA = as.numeric(Bac_16S_DNA)
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

# Merge metadata and blood 16S rDNA data
plot_data <- metadata_subsampled %>%
  inner_join(blood_16s_data, by = "study_no") %>%
  filter(!is.na(Inclusion_criteria))

# Check that data are available
if (nrow(plot_data) == 0) {
  stop("No matched samples were found after merging metadata and blood 16S rDNA data.")
}

# Save plot data
write_csv(
  plot_data,
  "outputs/Figure_5/Fig5A_blood_16S_rDNA_plot_data.csv"
)

# Summary statistics
summary_statistics <- plot_data %>%
  group_by(Inclusion_criteria) %>%
  summarise(
    n = n(),
    median = median(Bac_16S_DNA, na.rm = TRUE),
    IQR = IQR(Bac_16S_DNA, na.rm = TRUE),
    mean = mean(Bac_16S_DNA, na.rm = TRUE),
    sd = sd(Bac_16S_DNA, na.rm = TRUE),
    min = min(Bac_16S_DNA, na.rm = TRUE),
    max = max(Bac_16S_DNA, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  summary_statistics,
  "outputs/Figure_5/Fig5A_blood_16S_rDNA_summary_statistics.csv"
)

# Generate Figure 5A
fig5a <- ggplot(
  plot_data,
  aes(
    x = Inclusion_criteria,
    y = Bac_16S_DNA,
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
    position = position_jitter(width = 0.2),
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
    limits = c(0, 32)
  ) +
  scale_x_discrete(
    labels = c(
      "Healthy\ncontrols",
      "Stable\ncirrhotic",
      "Decomp\ncirrhotic",
      "ACLF"
    )
  ) +
  labs(
    x = "Group",
    y = "(Whole) blood 16S rDNA\n(16S copies/\u00B5g/L)"
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
  annotate("text", x = 2.5, y = 30.5, label = "p < 0.001", size = 2.5) +
  annotate("segment", x = 1, xend = 4, y = 29, yend = 29, linewidth = 0.35) +
  annotate("text", x = 3, y = 27.5, label = "p < 0.001", size = 2.5) +
  annotate("segment", x = 2, xend = 4, y = 26, yend = 26, linewidth = 0.35) +
  annotate("text", x = 3.5, y = 24.5, label = "p < 0.001", size = 2.5) +
  annotate("segment", x = 3, xend = 4, y = 23, yend = 23, linewidth = 0.35) +
  annotate("text", x = 2, y = 21.5, label = "p = 0.0029", size = 2.5) +
  annotate("segment", x = 1, xend = 3, y = 20, yend = 20, linewidth = 0.35)

# Save Figure 5A
ggsave(
  filename = "outputs/Figure_5/Fig5A_blood_16S_rDNA_boxplot.pdf",
  plot = fig5a,
  width = 5.5,
  height = 5.5,
  units = "cm",
  dpi = 300
)

# Statistical testing
anova_model <- aov(
  Bac_16S_DNA ~ Inclusion_criteria,
  data = plot_data
)

anova_summary <- summary(anova_model)
tukey_results <- TukeyHSD(anova_model)
anova_residuals <- residuals(anova_model)
residual_shapiro <- shapiro.test(anova_residuals)
kruskal_results <- kruskal.test(
  Bac_16S_DNA ~ Inclusion_criteria,
  data = plot_data
)

group_shapiro_results <- plot_data %>%
  group_by(Inclusion_criteria) %>%
  summarise(
    n = n(),
    shapiro_p_value = if_else(
      n >= 3,
      shapiro.test(Bac_16S_DNA)$p.value,
      NA_real_
    ),
    .groups = "drop"
  )

# Levene test, if the car package is installed
levene_results <- NULL

if (requireNamespace("car", quietly = TRUE)) {
  levene_results <- car::leveneTest(
    Bac_16S_DNA ~ Inclusion_criteria,
    data = plot_data
  )
} else {
  warning("Package 'car' is not installed, so Levene's test was not run.")
}

# Save statistical output
statistics_output_file <- "outputs/Figure_5/Fig5A_blood_16S_rDNA_statistics.txt"

sink(statistics_output_file)

cat("Figure 5A: Whole blood bacterial 16S rDNA by disease severity group\n")
cat("===================================================================\n\n")

cat("Input files:\n")
cat(metadata_file, "\n")
cat(blood_16s_file, "\n\n")

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

sink()

message("Figure saved to: outputs/Figure_5/Fig5A_blood_16S_rDNA_boxplot.pdf")
message("Plot data saved to: outputs/Figure_5/Fig5A_blood_16S_rDNA_plot_data.csv")
message("Summary statistics saved to: outputs/Figure_5/Fig5A_blood_16S_rDNA_summary_statistics.csv")
message("Statistics saved to: ", statistics_output_file)