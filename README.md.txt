# Liver patient Enterococcus analysis

This repository contains R scripts used for analysis and figure generation for the manuscript:

**Pathological expansion of gut microbiota-associated Enterococcus in advanced chronic liver disease corresponds with multilevel perturbations of the gut-liver-immune axis**

The scripts are organised by manuscript figure. Each script is intended to be run from the root folder of this repository.

## Data

The data required to run these scripts are provided as Supplementary Data files with the manuscript and/or are available from public repositories described in the manuscript.

Raw 16S rRNA amplicon sequencing data are available from NCBI SRA under BioProject PRJNA1223169.

Individual-level clinical metadata are not included directly in this repository because they derive from human participants and are subject to ethical and governance restrictions.

## Required R packages

The current scripts require:

```r
install.packages(c("tidyverse", "vegan", "car", "FSA", "ape"))

## Figure-to-script map

| Manuscript figure | Analysis | Script | Main output |
|---|---|---|---|
| Figure 1A | Observed number of bacterial genera | `scripts/Figure_1/Fig1A_number_of_genera_boxplot.R` | `outputs/Figure_1/Fig1A_number_of_genera_boxplot.pdf` |
| Figure 1B | Inverse Simpson diversity | `scripts/Figure_1/Fig1B_inverse_simpson_diversity_boxplot.R` | `outputs/Figure_1/Fig1B_inverse_simpson_diversity_boxplot.pdf` |
| Figure 1C | Genus-level microbiota NMDS | `scripts/Figure_1/Fig1C_genus_microbiota_NMDS_plot.R` | `outputs/Figure_1/Fig1C_genus_microbiota_NMDS_plot.pdf` |
| Figure 1D | Plasma bile acid NMDS | `scripts/Figure_1/Fig1D_plasma_bile_acids_NMDS_plot.R` | `outputs/Figure_1/Fig1D_plasma_bile_acids_NMDS_plot.pdf` |
| Figure 1E | Urinary metabolite NMDS | `scripts/Figure_1/Fig1E_urinary_metabolites_NMDS_plot.R` | `outputs/Figure_1/Fig1E_urinary_metabolites_NMDS_plot.pdf` |
| Figure 1F | Serum metabolite NMDS | `scripts/Figure_1/Fig1F_serum_metabolites_NMDS_plot.R` | `outputs/Figure_1/Fig1F_serum_metabolites_NMDS_plot.pdf` |
| Supplementary Figure 1A | Shannon diversity | `scripts/Supplementary_Figures/Supplementary_figure_1/SuppFig1A_shannon_diversity_boxplot.R` | `outputs/Supplementary_Figures/Supplementary_figure_1/SuppFig1A_shannon_diversity_boxplot.pdf` |
| Supplementary Figure 1B | Genus-level microbiota PCoA | `scripts/Supplementary_Figures/Supplementary_figure_1/SuppFig1B_genus_microbiota_PCoA_plot.R` | `outputs/Supplementary_Figures/Supplementary_figure_1/SuppFig1B_genus_microbiota_PCoA_plot.pdf` |



| Manuscript figure | Script | Run command |
|---|---|---|
| Figure 1A | `scripts/Figure_1/Fig1A_number_of_genera_boxplot.R` | `source("scripts/Figure_1/Fig1A_number_of_genera_boxplot.R")` |
| Figure 1B | `scripts/Figure_1/Fig1B_inverse_simpson_diversity_boxplot.R` | `source("scripts/Figure_1/Fig1B_inverse_simpson_diversity_boxplot.R")` |
| Figure 1C | `scripts/Figure_1/Fig1C_genus_microbiota_NMDS_plot.R` | `source("scripts/Figure_1/Fig1C_genus_microbiota_NMDS_plot.R")` |
| Figure 1D | `scripts/Figure_1/Fig1D_plasma_bile_acids_NMDS_plot.R` | `source("scripts/Figure_1/Fig1D_plasma_bile_acids_NMDS_plot.R")` |
| Figure 1E | `scripts/Figure_1/Fig1E_urinary_metabolites_NMDS_plot.R` | `source("scripts/Figure_1/Fig1E_urinary_metabolites_NMDS_plot.R")` |
| Figure 1F | `scripts/Figure_1/Fig1F_serum_metabolites_NMDS_plot.R` | `source("scripts/Figure_1/Fig1F_serum_metabolites_NMDS_plot.R")` |
| Supplementary Figure 1A | `scripts/Supplementary_Figures/Supplementary_figure_1/SuppFig1A_shannon_diversity_boxplot.R` | `source("scripts/Supplementary_Figures/Supplementary_figure_1/SuppFig1A_shannon_diversity_boxplot.R")` |
| Supplementary Figure 1B | `scripts/Supplementary_Figures/Supplementary_figure_1/SuppFig1B_genus_microbiota_PCoA_plot.R` | `source("scripts/Supplementary_Figures/Supplementary_figure_1/SuppFig1B_genus_microbiota_PCoA_plot.R")` |
| Figure 2B | `scripts/Figure_2/Fig2B_maaslin2_genus_differential_abundance.R` | `source("scripts/Figure_2/Fig2B_maaslin2_genus_differential_abundance.R")` |
| Figure 2C-H | `scripts/Figure_2/Fig2C_H_selected_genus_boxplots.R` | `source("scripts/Figure_2/Fig2C_H_selected_genus_boxplots.R")` |
| Figure 2H post-hoc power analysis | `scripts/Figure_2/Fig2H_enterococcus_posthoc_power_analysis.R` | `source("scripts/Figure_2/Fig2H_enterococcus_posthoc_power_analysis.R")` |
| Figure 3 bile acid MaAsLin2 analysis | `scripts/Figure_3/Fig3_maaslin2_bile_acids_differential_abundance.R` | `source("scripts/Figure_3/Fig3_maaslin2_bile_acids_differential_abundance.R")` |
| Figure 3 bile acid boxplots | `scripts/Figure_3/Fig3_bile_acid_boxplots.R` | `source("scripts/Figure_3/Fig3_bile_acid_boxplots.R")` |







