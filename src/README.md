- `multiplet_exploration_241115.R` - Original script from Amy.

- `multiplet_exploration_241115.Rmd` - MD modified script for processing `VCU-PC-062_gem_classification.csv`.
  - Input: `VCU-PC-062_gem_classification.csv`
  - Output: knitted PDF with plots

- `01_multiplet_exploration_all.Rmd` - Plots for all `gem_classification.csv` files.

- `02_multiplet_raw.Rmd` - using raw data. Also includes fixed_threshold filtering.

- `cellranger_10k_hgmm_gemx.sh` - downloads and runs cellranger on the 5p and 3p 10k 1:1 Mixture of Human HEK293T and Mouse NIH3T3 Cells, Chromium GEM-X Single Cell datasets from [10x](https://www.10xgenomics.com/datasets?configure%5BhitsPerPage%5D=50&configure%5BmaxValuesPerFacet%5D=1000&refinementList%5Bplatform%5D%5B0%5D=Chromium%20Single%20Cell&refinementList%5Bspecies%5D%5B0%5D=Human%2C%20Mouse)
  - Input: none, script will download the two datasets
  - Output: standard output from cellranger for each sample within the two datasets

- `cellranger_10k_hgmm_gemx.sh` - downloads and runs cellranger on SRR8890633 from [GSE129578](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129578), [McGinnis, C.S., Patterson, D.M., Winkler, J. et al.](https://doi.org/10.1038/s41592-019-0433-8)
  - Input: none, script will download the dataset
  - Output: standard output from cellranger for the sample
