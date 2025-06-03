- Run `MultipletR_Usage_Documentation.Rmd` and `example_walkthrough.R` in the `main` branch

# `multiplet_exploration_241115.Rmd`

+ MD: For all samples, generate plots p1 (line 117) and p2 (line 130)
+ MD: Test distribution type for raw total number of reads per cell
+ MD: Try "gamma" distribution
+ MD: Try adding minimum thresholds to catch high human/mouse percent cells that 
cannot be captured by a distribution

# General

- Process test datasets to get 10x classification, https://www.10xgenomics.com/datasets?configure%5BhitsPerPage%5D=50&configure%5BmaxValuesPerFacet%5D=1000&refinementList%5Bplatform%5D%5B0%5D=Chromium%20Single%20Cell&refinementList%5Bspecies%5D%5B0%5D=Human%2C%20Mouse

Universal 5' Gene Expression v3, Universal 3' Gene Expression v4

```
-rw-rw---- 1 jlmcclay pgxseq  7111683205 May  1  2024 CPF7_S3_R1_001.fastq.gz
-rw-rw---- 1 jlmcclay pgxseq 16111759831 May  1  2024 CPF7_S3_R2_001.fastq.gz

cellranger count --id CPF7 --transcriptome /lustre/home/mdozmorov/data/ExtData/10x/refdata-Rnor6-ensembl/Rnor_genome --fastqs /lustre/home/mdozmorov/data/WorkData/McClay/2024-05.scRNA-seq/00_raw --sample CPF7
```

- Classify test datasets with XenoCell, https://gitlab.com/XenoCell/XenoCell

# Biology

- How multiplet cells cluster? Detected by 10x and our methods.
