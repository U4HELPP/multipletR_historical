## Example Workflow for multipletR package

library("dplyr")
library("ggplot2")
library("reshape2")
library("gridExtra")
library("grid")
library("caret")
library("e1071")

## Source R functions until package is importable
source("R/utils.R")
source("R/viz_utils.R")
source("R/method_semilog_linear.R")

## Load example dataset
## Two types of example data sets are provided, LowMouse and StandardMouse.
## Each indicates the relative amount of mouse content in the sample.
## Standard Mouse samples did not use a mouse depletion kit prior to sampling and LowMouse samples did, so they have very few mouse cells called.
## Available files:
## VCU-BC-024_110857-61_lungmet_StandardMouse_gem_classification.csv
## VCU-BC-043_StandardMouse_gem_classification.csv
## VCU-CO-098_4245_LowMouse_gem_classification.csv
## VCU-PC-081_4220_LowMouse_gem_classification.csv

gem_data <- read.delim("example_data/VCU-BC-043_StandardMouse_gem_classification.csv", header = TRUE, row.names = 1, sep=",")

## Preprocess gem df
gem_df <- prep_gem_counts(result$df1, graft_col="GRCh38", host_col="mm10", call_col="call", 
                          graft_label="GRCh38", host_label="mm10")

## Generate summary plot before
summary_plot <- gem_classification_summary(gem_df, graft_col="GRCh38", host_col="mm10")

## Find Semilog Threshold
gem_df <- find_semilog_threshold(gem_df, graft_max=.25, graft_min=.10)

summary_plot3 <- gem_classification_summary(gem_df, method="SemiLog", graft_col="GRCh38", host_col="mm10", assigned_species_col="AssignedSpecies_SemiLog")

## Plot Semilog scatter
my_plot <- plot_semilog_threshold(gem_df, graft_max=.25, graft_min=.10, color="AssignedSpecies_SemiLog")


## Evaluate the classifications compared to 10X
metrics_df <- evaluate_classification(gem_df, "AssignedSpecies_10X")
print(metrics_df)



