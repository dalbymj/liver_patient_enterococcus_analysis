## Supplementary Figure 2 scripts

Scripts used to generate Supplementary Figure 2 panels related to bile acid analyses.

Run scripts from the root folder of the repository:

```r
getwd()
# Should return the path ending in:
# liver_patient_enterococcus_analysis
```

### Supplementary Figure 2 script map

| Supplementary figure | Script | Run command |
|---|---|---|
| Supplementary Figure 2A-B | `scripts/Supplementary_Figures/Supplementary_figure_2/SuppFig2A_B_bile_acid_boxplots.R` | `source("scripts/Supplementary_Figures/Supplementary_figure_2/SuppFig2A_B_bile_acid_boxplots.R")` |
| Supplementary Figure 2C | `scripts/Supplementary_Figures/Supplementary_figure_2/SuppFig2C_all_genus_bile_acid_correlation_heatmap.R` | `source("scripts/Supplementary_Figures/Supplementary_figure_2/SuppFig2C_all_genus_bile_acid_correlation_heatmap.R")` |

### Notes

- Supplementary Figure 2A-B generates bile acid abundance boxplots for selected bile acids not shown in the main Figure 3 panels.
- Supplementary Figure 2C generates the full genus-versus-bile-acid correlation heatmap.
- Supplementary Figure 2C selects the top 20 genera by total relative abundance.
- Supplementary Figure 2C orders genus rows by their correlation with taurochenodeoxycholic acid.
- Correlations are calculated after `bestNormalize` transformation.
- Pearson correlations are used.
- Benjamini-Hochberg correction is applied separately within each genus.
- Supplementary Figure 2 scripts are intended to be run independently from the repository root folder.