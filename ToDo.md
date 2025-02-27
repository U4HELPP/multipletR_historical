# `multiplet_exploration_241115.Rmd`

+ MD: For all samples, generate plots p1 (line 117) and p2 (line 130)
+ MD: Test distribution type for raw total number of reads per cell
+ MD: Try "gamma" distribution
+ MD: Try adding minimum thresholds to catch high human/mouse percent cells that 
cannot be captured by a distribution

# General

- Process test datasets to get 10x classification, https://www.10xgenomics.com/datasets?configure%5BhitsPerPage%5D=50&configure%5BmaxValuesPerFacet%5D=1000&refinementList%5Bplatform%5D%5B0%5D=Chromium%20Single%20Cell&refinementList%5Bspecies%5D%5B0%5D=Human%2C%20Mouse

- Classify test datasets with XenoCell, https://gitlab.com/XenoCell/XenoCell

# Biology

- How multiplet cells cluster? Detected by 10x and our methods.