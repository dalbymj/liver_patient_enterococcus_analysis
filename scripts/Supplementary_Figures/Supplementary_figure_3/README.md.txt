## Supplementary Figure 3 scripts

Scripts used to generate Supplementary Figure 3 panels related to urinary metabolite analyses.

Run scripts from the root folder of the repository:

```r
getwd()
# Should return the path ending in:
# liver_patient_enterococcus_analysis
```

### Supplementary Figure 3 script map

| Supplementary figure | Script | Run command |
|---|---|---|
| Supplementary Figure 3F | `scripts/Supplementary_Figures/Supplementary_figure_3/SuppFig3F_all_genus_urinary_metabolite_correlation_heatmap.R` | `source("scripts/Supplementary_Figures/Supplementary_figure_3/SuppFig3F_all_genus_urinary_metabolite_correlation_heatmap.R")` |

### Notes

- Supplementary Figure 3F generates the full genus-versus-urinary-metabolite correlation heatmap.
- The script selects the top 20 genera by total relative abundance.
- The script includes all urinary metabolites from `For_TS` to `Lac_TS`.
- Genus rows are ordered by their correlation with hippurate.
- Correlations are calculated after `bestNormalize` transformation.
- Pearson correlations are used.
- Benjamini-Hochberg correction is applied separately within each genus.
- Supplementary Figure 3 scripts are intended to be run independently from the repository root folder.