## =============================================================================
## Complete Walkthrough: Alexandra's Median-Based Multiplet Detection Pipeline
##
## This walkthrough demonstrates ALL functions in FUNCTIONS_ALEXANDRA_MEDIAN.R
## for identifying multiplets and low-quality cells in human-mouse xenograft
## scRNA-seq data using median-based statistical approach with flipping transformation.
##
## =============================================================================

## Load required libraries
library(fitdistrplus)
library(ggplot2)
library(gridExtra)
library(dplyr)

## Source the functions
source("FUNCTIONS_ALEXANDRA_MEDIAN.R")

## =============================================================================
## SECTION 1: DATA LOADING AND INITIAL EXPLORATION
## =============================================================================

## Load the 10X Genomics cell classification data
data <- read.csv("10X_gem_classifications/VCU-CO-063_gem_classification.csv")

## Function 1: plot_initial_data - Visualize raw data before classification
initial_summary <- plot_initial_data(data,
                                     title = "VCU-CO-063: Raw Data Before Classification",
                                     show_stats = TRUE)

## =============================================================================
## SECTION 2: LOW-READ CLASSIFICATION
## =============================================================================

## Function 2: plot_low_reads_multi_sigma - Test multiple sigma values
low_sigma_results <- plot_low_reads_multi_sigma(data,
                                                sigma_levels = c(0.5, 0.75, 1.0, 1.25, 1.5, 2.0),
                                                distribution = "auto",
                                                verbose = FALSE)

## Function 3: classify_low_reads - Apply low-read classification
data <- classify_low_reads(data,
                           n_sd = 1.0,
                           distribution = "auto",
                           verbose = TRUE)

## Function 4: plot_low_reads - Visualize low-read classification
## Both plot types demonstrated
plot_low_reads(data, n_sd = 1.0, plot_type = "single")
plot_low_reads(data, n_sd = 1.0, plot_type = "standard")

## =============================================================================
## SECTION 3: HIGH-READ CLASSIFICATION
## =============================================================================

## Function 5: plot_high_reads_multi_sigma - Test multiple sigma values
high_sigma_results <- plot_high_reads_multi_sigma(data,
                                                  sigma_levels = c(0.5, 0.75, 1.0, 1.25, 1.5, 2.0),
                                                  distribution = "auto",
                                                  verbose = FALSE)

## Function 6: classify_high_reads - Apply high-read classification
data <- classify_high_reads(data,
                            n_sd = 1.5,
                            distribution = "auto",
                            verbose = TRUE)

## Function 7: plot_high_reads - Visualize high-read classification
## Both plot types demonstrated
plot_high_reads(data, n_sd = 1.5, plot_type = "single")
plot_high_reads(data, n_sd = 1.5, plot_type = "standard")

## =============================================================================
## SECTION 4: MULTIPLET CLASSIFICATION
## =============================================================================

## Function 8: plot_multiplets_multi_sigma - Test multiple sigma values
multiplet_sigma_results <- plot_multiplets_multi_sigma(data,
                                                       sigma_levels = c(1.0, 1.5, 2.0, 2.5, 3.0, 4.0),
                                                       distribution = "auto",
                                                       verbose = FALSE)

## Function 9: classify_multiplets - Apply multiplet classification
data <- classify_multiplets(data,
                            n_sd = 3.0,
                            distribution = "auto",
                            verbose = TRUE)

## Function 10: plot_multiplets - Visualize multiplet classification
## Both plot types demonstrated
plot_multiplets(data, n_sd = 3.0, plot_type = "single")
plot_multiplets(data, n_sd = 3.0, plot_type = "standard")

## Function 11: plot_flipped_space - Show mathematical basis
flipping_params <- plot_flipped_space(data,
                                      sigma_levels = c(1.0, 1.5, 2.0, 2.5, 3.0, 4.0),
                                      selected_sigma = 3.0)

## =============================================================================
## SECTION 5: MULTIPLET VALIDATION
## =============================================================================

## Function 12: validate_multiplet_distributions - Basic validation
validation_results <- validate_multiplet_distributions(data,
                                                       sigma_levels = c(2.0, 2.5, 3.0, 3.5),
                                                       use_log = TRUE,
                                                       distribution = "auto",
                                                       verbose = TRUE)

## Function 13: validate_multiplet_distributions_pub - Publication version
pub_validation <- validate_multiplet_distributions_pub(data,
                                                       sigma_levels = c(2.0, 2.5, 3.0, 3.5),
                                                       use_log = TRUE,
                                                       distribution = "auto",
                                                       color_scheme = "default",
                                                       verbose = TRUE)

## Function 14: validate_10X_multiplets - Validate 10X method
tenx_validation <- validate_10X_multiplets(data,
                                           use_log = TRUE,
                                           verbose = TRUE)

## =============================================================================
## SECTION 6: COMPREHENSIVE VALIDATION TABLE
## =============================================================================

## Function 15: create_multiplet_validation_table - Generate validation table
validation_table <- create_multiplet_validation_table(data,
                                                      sigma_levels = c(1.0, 1.5, 2.0, 2.5, 3.0, 3.5),
                                                      use_log = TRUE,
                                                      distribution = "auto",
                                                      verbose = TRUE)

## Function 16: export_validation_table - Export table in multiple formats
export_validation_table(validation_table, filename_base = "VCU-CO-063_validation")

## =============================================================================
## SECTION 7: COMBINED VISUALIZATION
## =============================================================================

## Function 17: plot_combined_classification - All three steps together
combined_results <- plot_combined_classification(data,
                                                 low_sd = 1.0,
                                                 high_sd = 1.5,
                                                 multiplet_sd = 3.0,
                                                 xlim_max = NULL,
                                                 show_stats = TRUE,
                                                 show_legend = TRUE,
                                                 point_size = 0.5,
                                                 point_alpha = 0.6,
                                                 title = "Complete Three-Step Classification",
                                                 run_missing_classifications = FALSE,
                                                 verbose = TRUE)

## =============================================================================
## SECTION 8: HELPER FUNCTIONS DEMONSTRATION
## =============================================================================

## These helper functions are used internally by other functions:
## - standardize_column_names() - Called by all classification functions
## - get_fitted_density() - Used by plotting functions
## - calculate_distribution_overlap() - Used in validation
## - calculate_low_threshold() - Used in low-read classification
## - calculate_high_threshold() - Used in high-read classification
## - calculate_multiplet_threshold() - Used in multiplet classification
## - get_threshold_value() - Used for threshold extraction
## - get_multiplet_thresholds() - Used for multiplet threshold extraction

## Note: The _viz versions of threshold calculation functions are also used internally

## =============================================================================
## SECTION 9: SAVE RESULTS
## =============================================================================

## Add combined classification column
data$final_classification <- paste(data$low_classification,
                                   data$high_classification,
                                   data$multiplet_classification,
                                   sep = "_")

## Save the fully classified dataset
write.csv(data, "VCU-CO-063_fully_classified_median_method.csv", row.names = FALSE)

## Save validation metrics
write.csv(validation_table, "VCU-CO-063_validation_metrics.csv", row.names = FALSE)

## =============================================================================
## ANALYSIS COMPLETE
## =============================================================================
