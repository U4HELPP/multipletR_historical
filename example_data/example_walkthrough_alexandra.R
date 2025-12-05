## =============================================================================
## Example Walkthrough for multipletR Three-Step Classification Package
##
## This script demonstrates the complete workflow for identifying and removing
## low-quality cells and multiplets using statistical outlier detection:
## Step 1: Low-read classification (abnormally low total reads)
## Step 2: High-read classification (potential doublets with high reads)
## Step 3: Multiplet classification (species mixing analysis using flipping approach)
## =============================================================================

## Load required libraries
library("dplyr")
library("fitdistrplus")

## Source R functions until package is importable
source("FUNCTIONS_ALEXANDRA.R")

## =============================================================================
## Load example dataset
## =============================================================================

## Two types of example datasets are available in the 10X_gem_classifications folder:
## 1. StandardMouse samples: Natural human-mouse xenograft samples with mixed species content
## 2. LowMouse samples: Samples processed with mouse depletion, resulting in fewer mouse cells
##
## Available example files:
## VCU-CO-063_gem_classification.csv (StandardMouse)
## VCU-BC-043_StandardMouse_gem_classification.csv
## VCU-CO-098_4245_LowMouse_gem_classification.csv
## VCU-PC-081_4220_LowMouse_gem_classification.csv

gem_data <- read.csv("10X_gem_classifications/VCU-CO-063_gem_classification.csv")

## =============================================================================
## Step 1: Initial Data Exploration
## =============================================================================

## Visualize raw data before any classification
initial_plot <- plot_initial_data(gem_data,
                                  title = "Raw Data Overview",
                                  show_stats = TRUE)

## =============================================================================
## Step 2: Low-Read Classification (Quality Control Step 1)
## =============================================================================

## Test multiple sigma values to find optimal threshold
low_sigma_comparison <- plot_low_reads_multi_sigma(gem_data,
                                                   sigma_levels = c(0.5, 0.75, 1.0, 1.25, 1.5, 2.0),
                                                   distribution = "auto")

## Apply low-read classification with chosen sigma (σ = 1.0 recommended)
gem_data <- classify_low_reads(gem_data,
                               n_sd = 1.0,
                               distribution = "auto",
                               verbose = TRUE)

## Generate detailed visualization for low-read classification
low_reads_plot <- plot_low_reads(gem_data, n_sd = 1.0, plot_type = "standard")

## =============================================================================
## Step 3: High-Read Classification (Quality Control Step 2)
## =============================================================================

## Test multiple sigma values for high-read detection
high_sigma_comparison <- plot_high_reads_multi_sigma(gem_data,
                                                     sigma_levels = c(1.0, 1.25, 1.5, 1.75, 2.0, 2.5),
                                                     distribution = "auto")

## Apply high-read classification with chosen sigma (σ = 1.5 recommended)
gem_data <- classify_high_reads(gem_data,
                                n_sd = 1.5,
                                distribution = "auto",
                                verbose = TRUE)

## Generate detailed visualization for high-read classification
high_reads_plot <- plot_high_reads(gem_data, n_sd = 1.5, plot_type = "standard")

## =============================================================================
## Step 4: Multiplet Classification (Species Separation)
## =============================================================================

## Test multiple sigma values for multiplet detection
multiplet_sigma_comparison <- plot_multiplets_multi_sigma(gem_data,
                                                          sigma_levels = c(2.0, 2.5, 3.0, 3.5, 4.0, 4.5),
                                                          distribution = "auto")

## Apply multiplet classification with chosen sigma (σ = 3.0 recommended)
gem_data <- classify_multiplets(gem_data,
                                n_sd = 3.0,
                                distribution = "auto",
                                verbose = TRUE)

## Generate detailed visualization for multiplet classification
multiplet_plot <- plot_multiplets(gem_data, n_sd = 3.0, plot_type = "standard")

## Visualize the flipping approach used for multiplet detection
flipping_analysis <- plot_flipped_space(gem_data,
                                        sigma_levels = c(2.0, 2.5, 3.0, 3.5, 4.0, 4.5),
                                        selected_sigma = 3.0)

## =============================================================================
## Step 5: Combined Analysis and Results
## =============================================================================

## Generate comprehensive combined classification plot
combined_plot <- plot_combined_classification(gem_data,
                                              low_sd = 1.0,
                                              high_sd = 1.5,
                                              multiplet_sd = 3.0,
                                              show_stats = TRUE,
                                              show_legend = TRUE,
                                              title = "Complete Three-Step Classification Results")

## =============================================================================
## Step 6: Alternative Visualization Options
## =============================================================================

## Create single plots suitable for presentations
single_low_plot <- plot_low_reads(gem_data, n_sd = 1.0, plot_type = "single")
single_high_plot <- plot_high_reads(gem_data, n_sd = 1.5, plot_type = "single")
single_multiplet_plot <- plot_multiplets(gem_data, n_sd = 3.0, plot_type = "single")

