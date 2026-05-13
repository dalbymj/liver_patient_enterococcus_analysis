## Figure 4 scripts

Scripts used to generate Figure 4 panels and related supporting analyses.

Run scripts from the root folder of the repository:

```r
getwd()
# Should return the path ending in:
# liver_patient_enterococcus_analysis
```

### Recommended run order

Run the MaAsLin2 adjusted genus analysis before the Figure 4A-E boxplot script, because the boxplot script can read and summarise the MaAsLin2 results.

| Manuscript figure | Script | Run command |
|---|---|---|
| Figure 4A-E MaAsLin2 analysis | `scripts/Figure_4/Fig4_maaslin2_genus_antibiotics_food_lactulose.R` | `source("scripts/Figure_4/Fig4_maaslin2_genus_antibiotics_food_lactulose.R")` |
| Figure 4A-E | `scripts/Figure_4/Fig4A_E_genus_antibiotics_food_lactulose_boxplots.R` | `source("scripts/Figure_4/Fig4A_E_genus_antibiotics_food_lactulose_boxplots.R")` |
| Figure 4F | `scripts/Figure_4/Fig4F_antibiotic_exposure_heatmap.R` | `source("scripts/Figure_4/Fig4F_antibiotic_exposure_heatmap.R")` |

### Notes

- The Figure 4 MaAsLin2 script tests genus associations with disease severity, antibiotic exposure, oral food intake, and lactulose exposure.
- Figure 4A-E generates selected genus abundance boxplots stratified by antibiotic exposure, oral food intake, and lactulose exposure.
- Figure 4A-E includes Enterococcus, Ruminococcus, Roseburia, and Bifidobacterium panels.
- Figure 4F generates the antibiotic exposure heatmap for liver disease samples only.
- Figure 4F uses binary Jaccard distance and average-linkage hierarchical clustering for both samples and antibiotic variables.
- Figure 4F includes an Enterococcus relative abundance row annotation, disease group annotation, and antibiotic exposure annotation.
- Figure 4 scripts use genus-level faecal microbiota relative abundance data.
- Figure 4 scripts are intended to be run independently from the repository root folder.