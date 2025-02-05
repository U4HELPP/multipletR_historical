## Example Workflow for multipletR package

library("dplyr")
library("ggplot2")
library("reshape2")
library("gridExtra")
library("grid")
library("caret")
library("e1071")



## Load dataset
gems <- read.delim("/lustre/home/harrell_lab/scRNASeq/exploratory_analyses/multiplet_exploration/VCUBC024lungmet_GEM.txt", header=FALSE, sep=",")
gem_data <- read.delim(paste0(gems[1,2], "/analysis/gem_classification.csv"), header = TRUE, row.names = 1, sep=",")

## Preprocess gem df
gem_df <- prep_gem_counts(gem_data)

## Generate summary plot before
summary_plot <- gem_classification_summary(gem_df)

## Find Semilog Threshold
gem_df <- find_semilog_threshold(gem_df, graft_max=.25, graft_min=.10)

summary_plot2 <- gem_classification_summary(gem_df, method="SemiLog", assigned_species_col="AssignedSpecies_SemiLog")

## Plot Semilog scatter
my_plot <- plot_semilog_threshold(gem_df, graft_max=.25, graft_min=.10, color="AssignedSpecies_SemiLog")

## Evaluate the classifications compared to 10X
metrics_df <- evaluate_classification(gem_df, "AssignedSpecies_10X")
print(metrics_df)



