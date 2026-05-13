## Figure 1 scripts

Scripts used to generate Figure 1 panels and related Supplementary Figure 1 panels.

Run scripts from the root folder of the repository:

```r
getwd()
# Should return the path ending in:
# liver_patient_enterococcus_analysis
```

### Figure 1 script map

| Manuscript figure | Script | Run command |
|---|---|---|
| Figure 1A | `scripts/Figure_1/Fig1A_number_of_genera_boxplot.R` | `source("scripts/Figure_1/Fig1A_number_of_genera_boxplot.R")` |
| Figure 1B | `scripts/Figure_1/Fig1B_inverse_simpson_diversity_boxplot.R` | `source("scripts/Figure_1/Fig1B_inverse_simpson_diversity_boxplot.R")` |
| Figure 1C | `scripts/Figure_1/Fig1C_genus_microbiota_NMDS_plot.R` | `source("scripts/Figure_1/Fig1C_genus_microbiota_NMDS_plot.R")` |

### Related supplementary figure scripts

| Supplementary figure | Script | Run command |
|---|---|---|
| Supplementary Figure 1A | `scripts/Supplementary_Figures/Supplementary_figure_1/SuppFig1A_shannon_diversity_boxplot.R` | `source("scripts/Supplementary_Figures/Supplementary_figure_1/SuppFig1A_shannon_diversity_boxplot.R")` |
| Supplementary Figure 1B | `scripts/Supplementary_Figures/Supplementary_figure_1/SuppFig1B_genus_microbiota_PCoA_plot.R` | `source("scripts/Supplementary_Figures/Supplementary_figure_1/SuppFig1B_genus_microbiota_PCoA_plot.R")` |

### Notes

- Figure 1A calculates and plots observed genus richness.
- Figure 1B calculates and plots inverse Simpson diversity.
- Supplementary Figure 1A calculates and plots Shannon diversity.
- Figure 1C generates the genus-level NMDS plot and saves PERMANOVA outputs.
- Supplementary Figure 1B generates the genus-level PCoA plot.
- The NMDS and PCoA scripts use Bray-Curtis dissimilarity on genus-level faecal microbiota profiles.
- Figure 1 scripts are intended to be run independently from the repository root folder.