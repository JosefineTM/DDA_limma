# LIMMA FOR DIFFERENTIAL ABUNDANCE ANALYSIS
## Limma
The objective of limma analysis is to detect statistically significant abundant proteins in two different condition groups. In this case, it will estimate Tumor and Normal group mean abundance for each protein and test whether the group means are significantly different. To this end, limma implements eBayes moderated t-statistics which 1) estimate the variance for each gene; 2) the common variance trend across all genes; 3) Shrinks each gene's variance estimate toward the common trend; 4) compute the statistical significance of the estimated differences.
The method takes as input the protein abundance matrix and returns as output a data table containing the Differential Abundant Proteins, each associated with log2FoldChange, p-values (adjusted) and other score metrics. 
## Requirments
To run limma locally, a conda env with R has been created
```bash
conda create --name myenv r-base=4.3
```
```bash
conda activate myenv r-base=4.3
``` 
The following packages should be installed to run the script:
- "BiocManager"
- "limma"
- "dplyr"
- "tibble"
- "optparse"

## Running the analysis
The script takes as inputs:
- The protein abundance matrix (nrows = protein; ncols = samples)
- A dataframe containing for each sample ID, its assigned group label: Tumor (T) or Normal (N)
- A dataframe containing for each gene name (hugo name), the corresponding protein database ID (ENSEMBL Protein ID)
The script can be run as follows:

```R
Rscript limma_method.r -d <path_to_dataset> -l <path_to_group_label_df> -m <path_to_gene_proteinIDs_df> 
``` 
The output is a limma data table containing for each protein, its corresponding gene name, the effect size value and the p.value computed by the limma analysis. The data table is automatically saved as .csv in this folder.

