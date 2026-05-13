## Figure 2 scripts

Scripts used to generate Figure 2 panels and related supporting analyses.

Run scripts from the root folder of the repository:

```r
getwd()
# Should return the path ending in:
# liver_patient_enterococcus_analysis
```

### Recommended run order

Run the MaAsLin2 differential abundance script before the selected genus boxplot script, because the boxplot script can read and summarise the MaAsLin2 results.

| Manuscript figure | Script | Run command |
|---|---|---|
| Figure 2A | `scripts/Figure_2/Fig2A_genus_stacked_barplot.R` | `source("scripts/Figure_2/Fig2A_genus_stacked_barplot.R")` |
| Figure 2B | `scripts/Figure_2/Fig2B_maaslin2_genus_differential_abundance.R` | `source("scripts/Figure_2/Fig2B_maaslin2_genus_differential_abundance.R")` |
| Figure 2C-H | `scripts/Figure_2/Fig2C_H_selected_genus_boxplots.R` | `source("scripts/Figure_2/Fig2C_H_selected_genus_boxplots.R")` |
| Figure 2 supporting analysis | `scripts/Figure_2/Fig2H_enterococcus_posthoc_power_analysis.R` | `source("scripts/Figure_2/Fig2H_enterococcus_posthoc_power_analysis.R")` |

### Notes

- Figure 2A generates the stacked bar plot of the most abundant genera across samples.
- Figure 2B runs MaAsLin2 genus-level differential abundance testing by disease severity group.
- Figure 2C-H generates selected genus abundance boxplots, including Enterococcus and depleted commensal genera.
- The Enterococcus post-hoc power analysis script records the supporting power calculation and does not need to be run before the main plotting scripts.
- Figure 2 scripts use genus-level faecal microbiota relative abundance data.
- Figure 2 scripts are intended to be run from the repository root folder.