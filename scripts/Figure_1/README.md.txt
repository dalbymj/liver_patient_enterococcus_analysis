# Figure 1 scripts

Scripts used to generate Figure 1 panels.

## Fig1A_number_of_genera_boxplot.R

Generates the Figure 1A observed genus richness boxplot and associated statistical output.

Input files:
- data/raw_genus_data_15.04.2024.csv
- data/clinical_metadata_and_data_renamed_30.01.2025.csv

Output files:
- outputs/Figure_1/Fig1A_number_of_genera_boxplot.pdf
- outputs/Figure_1/Fig1A_number_of_genera_statistics.txt

To run from the repository root:

```r
source("scripts/Figure_1/Fig1A_number_of_genera_boxplot.R")