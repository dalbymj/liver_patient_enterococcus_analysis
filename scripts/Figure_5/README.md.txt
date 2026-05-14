# Figure 5 scripts

Scripts used to generate Figure 5 panels. Run all scripts from the root folder of the repository:

setwd("C:/Users/dalby/Documents/Hall Lab/GitHub/liver_patient_enterococcus_analysis")


## Figure-to-script map

| Manuscript figure | Script | Run command |
|---|---|---|
| Figure 5A | scripts/Figure_5/Fig5A_blood_16S_rDNA_boxplot.R | source("scripts/Figure_5/Fig5A_blood_16S_rDNA_boxplot.R") |
| Figure 5B | scripts/Figure_5/Fig5B_genus_blood_16S_rDNA_correlation_heatmap.R | source("scripts/Figure_5/Fig5B_genus_blood_16S_rDNA_correlation_heatmap.R") |
| Figure 5C-D | scripts/Figure_5/Fig5C_D_genus_blood_16S_rDNA_scatterplots.R | source("scripts/Figure_5/Fig5C_D_genus_blood_16S_rDNA_scatterplots.R") |
| Figure 5E | scripts/Figure_5/Fig5E_faecal_calprotectin_boxplot.R | source("scripts/Figure_5/Fig5E_faecal_calprotectin_boxplot.R") |
| Figure 5F | scripts/Figure_5/Fig5F_genus_faecal_calprotectin_correlation_heatmap.R | source("scripts/Figure_5/Fig5F_genus_faecal_calprotectin_correlation_heatmap.R") |
| Figure 5G-H | scripts/Figure_5/Fig5G_H_genus_faecal_calprotectin_scatterplots.R | source("scripts/Figure_5/Fig5G_H_genus_faecal_calprotectin_scatterplots.R") |
| Figure 5I | scripts/Figure_5/Fig5I_monocyte_facs_microbiome_correlation_heatmap.R | source("scripts/Figure_5/Fig5I_monocyte_facs_microbiome_correlation_heatmap.R") |
| Figure 5J-K | scripts/Figure_5/Fig5J_K_monocyte_facs_microbiome_scatterplots.R | source("scripts/Figure_5/Fig5J_K_monocyte_facs_microbiome_scatterplots.R") |


## Script descriptions

Fig5A_blood_16S_rDNA_boxplot.R

Generates the Figure 5A boxplot showing whole blood bacterial 16S rDNA across disease severity groups.

Main outputs:
- outputs/Figure_5/Fig5A_blood_16S_rDNA_boxplot.pdf
- outputs/Figure_5/Fig5A_blood_16S_rDNA_plot_data.csv
- outputs/Figure_5/Fig5A_blood_16S_rDNA_summary_statistics.csv
- outputs/Figure_5/Fig5A_blood_16S_rDNA_statistics.txt


Fig5B_genus_blood_16S_rDNA_correlation_heatmap.R

Generates the Figure 5B heatmap showing correlations between selected genera and whole blood bacterial 16S rDNA. Genus abundances and 16S rDNA values are transformed before Pearson correlation analysis.

Main outputs:
- outputs/Figure_5/Fig5B_genus_blood_16S_rDNA_correlation_heatmap.pdf
- outputs/Figure_5/Fig5B_genus_blood_16S_rDNA_correlation_table.csv
- outputs/Figure_5/Fig5B_genus_blood_16S_rDNA_correlation_plot_data.csv
- outputs/Figure_5/Fig5B_bestNormalize_transformation_log.csv
- outputs/Figure_5/Fig5B_genus_blood_16S_rDNA_correlation_summary.txt


Fig5C_D_genus_blood_16S_rDNA_scatterplots.R

Generates the Figure 5C-D scatter plots showing selected genus associations with whole blood bacterial 16S rDNA.

Main outputs:
- outputs/Figure_5/Fig5C_Lachnospira_blood_16S_rDNA_scatterplot.pdf
- outputs/Figure_5/Fig5D_Enterococcus_blood_16S_rDNA_scatterplot.pdf
- outputs/Figure_5/Fig5C_D_genus_blood_16S_rDNA_scatterplot_data.csv
- outputs/Figure_5/Fig5C_D_genus_blood_16S_rDNA_scatterplot_statistics.csv
- outputs/Figure_5/Fig5C_D_bestNormalize_transformation_log.csv
- outputs/Figure_5/Fig5C_D_genus_blood_16S_rDNA_scatterplot_summary.txt


Fig5E_faecal_calprotectin_boxplot.R

Generates the Figure 5E boxplot showing faecal calprotectin across disease severity groups.

Main outputs:
- outputs/Figure_5/Fig5E_faecal_calprotectin_boxplot.pdf
- outputs/Figure_5/Fig5E_faecal_calprotectin_plot_data.csv
- outputs/Figure_5/Fig5E_faecal_calprotectin_summary_statistics.csv
- outputs/Figure_5/Fig5E_faecal_calprotectin_statistics.txt


Fig5F_genus_faecal_calprotectin_correlation_heatmap.R

Generates the Figure 5F heatmap showing correlations between selected genera and faecal calprotectin. The analysis selects the top 20 genera by total relative abundance, transforms variables, and applies Pearson correlation testing with Benjamini-Hochberg correction within each genus.

Main outputs:
- outputs/Figure_5/Fig5F_genus_faecal_calprotectin_correlation_heatmap.pdf
- outputs/Figure_5/Fig5F_genus_faecal_calprotectin_correlation_table.csv
- outputs/Figure_5/Fig5F_genus_faecal_calprotectin_correlation_plot_data.csv
- outputs/Figure_5/Fig5F_bestNormalize_transformation_log.csv
- outputs/Figure_5/Fig5F_genus_faecal_calprotectin_correlation_summary.txt


Fig5G_H_genus_faecal_calprotectin_scatterplots.R

Generates the Figure 5G-H scatter plots showing selected genus associations with faecal calprotectin.

Main outputs:
- outputs/Figure_5/Fig5G_Lachnospira_faecal_calprotectin_scatterplot.pdf
- outputs/Figure_5/Fig5H_Enterococcus_faecal_calprotectin_scatterplot.pdf
- outputs/Figure_5/Fig5G_H_genus_faecal_calprotectin_scatterplot_data.csv
- outputs/Figure_5/Fig5G_H_genus_faecal_calprotectin_scatterplot_statistics.csv
- outputs/Figure_5/Fig5G_H_bestNormalize_transformation_log.csv
- outputs/Figure_5/Fig5G_H_genus_faecal_calprotectin_scatterplot_summary.txt


Fig5I_monocyte_facs_microbiome_correlation_heatmap.R

Generates the Figure 5I heatmap showing correlations between monocyte FACS response variables and selected microbiome/translocation variables: whole blood bacterial 16S rDNA, number of detected genera, Enterococcus, and Lachnospira.

Main outputs:
- outputs/Figure_5/Fig5I_monocyte_facs_microbiome_correlation_heatmap.pdf
- outputs/Figure_5/Fig5I_monocyte_facs_microbiome_correlation_table.csv
- outputs/Figure_5/Fig5I_monocyte_facs_microbiome_correlation_plot_data.csv
- outputs/Figure_5/Fig5I_bestNormalize_transformation_log.csv
- outputs/Figure_5/Fig5I_monocyte_facs_microbiome_correlation_summary.txt


Fig5J_K_monocyte_facs_microbiome_scatterplots.R

Generates the Figure 5J-K scatter plots showing Enterococcus associations with LPS- and PAM3-induced monocyte TNF-alpha responses.

Main outputs:
- outputs/Figure_5/Fig5J_Enterococcus_LPS_TNFa_scatterplot.pdf
- outputs/Figure_5/Fig5K_Enterococcus_PAM3_TNFa_scatterplot.pdf
- outputs/Figure_5/Fig5J_K_monocyte_facs_microbiome_scatterplot_data.csv
- outputs/Figure_5/Fig5J_K_monocyte_facs_microbiome_scatterplot_statistics.csv
- outputs/Figure_5/Fig5J_K_bestNormalize_transformation_log.csv
- outputs/Figure_5/Fig5J_K_monocyte_facs_microbiome_scatterplot_summary.txt


## Notes

The Figure 5 scripts expect the required input data files to be present in the repository data/ folder.

Several scripts use bestNormalize transformation before Pearson correlation analysis. Where matched sample numbers are small, bestNormalize() may require out_of_sample = FALSE to avoid k-fold cross-validation errors.

All generated figure panels and associated statistics are saved to:

outputs/Figure_5/