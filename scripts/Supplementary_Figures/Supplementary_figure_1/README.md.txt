## Supplementary Figure 1 scripts

Scripts used to generate Supplementary Figure 1 panels related to alpha diversity and genus-level microbiota composition.

Run scripts from the root folder of the repository:

```r
getwd()
# Should return the path ending in:
# liver_patient_enterococcus_analysis
```

### Supplementary Figure 1 script map

| Supplementary figure | Script | Run command |
|---|---|---|
| Supplementary Figure 1A | `scripts/Supplementary_Figures/Supplementary_figure_1/SuppFig1A_shannon_diversity_boxplot.R` | `source("scripts/Supplementary_Figures/Supplementary_figure_1/SuppFig1A_shannon_diversity_boxplot.R")` |
| Supplementary Figure 1B | `scripts/Supplementary_Figures/Supplementary_figure_1/SuppFig1B_genus_microbiota_PCoA_plot.R` | `source("scripts/Supplementary_Figures/Supplementary_figure_1/SuppFig1B_genus_microbiota_PCoA_plot.R")` |

### Notes

- Supplementary Figure 1A calculates Shannon diversity and generates the Shannon diversity boxplot.
- Supplementary Figure 1B generates the genus-level PCoA plot.
- The PCoA script uses Bray-Curtis dissimilarity on genus-level faecal microbiota profiles.
- Supplementary Figure 1 scripts are intended to be run independently from the repository root folder.