#' @title Prepare GEM Counts for Analysis and Plotting
#' @description This function takes a data frame of read counts from a GEM Classification file and adds additional columns for analysis and plotting.
#' It allows for dynamic column names and classification labels, generalized for PDX (Graft vs Host).
#'
#' @param gem_df A data frame containing read counts.
#' @param graft_col Character. The name of the column containing Graft read counts. Default is "GRCh38".
#' @param host_col Character. The name of the column containing Host read counts. Default is "GRCm39".
#' @param call_col Character. The name of the column containing the classification call. Default is "call".
#' @param graft_label Character. The value in `call_col` that identifies a Graft cell. Default is "GRCh38".
#' @param host_label Character. The value in `call_col` that identifies a Host cell. Default is "GRCm39".
#' @param multiplet_label Character. The value in `call_col` that identifies a Multiplet. Default is "Multiplet".
#' @param suffix Character. The suffix to append to the "AssignedSpecies_" column name. Default is "10X".
#'
#' @return A data frame with the following additional columns:
#' \itemize{
#'   \item \code{GraftDiff}: Difference between graft and host read counts.
#'   \item \code{HostDiff}: Difference between host and graft read counts.
#'   \item \code{percentHost}: Proportion of reads that map to host.
#'   \item \code{totalReads}: Total number of reads (sum of graft and host counts).
#'   \item \code{totalReadsLog}: Logarithm of the total number of reads.
#'   \item \code{AssignedSpecies_<suffix>}: Normalized classification ("Graft", "Host", "Multiplet", or NA).
#' }
#' @examples
#' \dontrun{
#' # Example with default names
#' gem_df <- data.frame(barcode = c("bc1", "bc2"), GRCh38 = c(100, 10), GRCm39 = c(10, 200), call = c("GRCh38", "GRCm39"))
#' result <- prep_gem_counts(gem_df)
#'
#' # Example with custom names
#' gem_df_custom <- data.frame(barcode = c("bc1"), g_reads = c(100), h_reads = c(50), assign = c("Graft"))
#' result <- prep_gem_counts(gem_df_custom, graft_col="g_reads", host_col="h_reads", call_col="assign", graft_label="Graft")
#' }
#' @export
prep_gem_counts <- function(gem_df, 
                            graft_col = "GRCh38", 
                            host_col = "GRCm39", 
                            call_col = "call",
                            graft_label = "GRCh38", 
                            host_label = "GRCm39", 
                            multiplet_label = "Multiplet",
                            suffix = "10X"){
  
  # Ensure the specified columns exist in the dataframe
  required_cols <- c(graft_col, host_col, call_col)
  missing_cols <- required_cols[!required_cols %in% names(gem_df)]
  if (length(missing_cols) > 0) {
    stop(paste("The following required columns are missing from gem_df:", paste(missing_cols, collapse = ", ")))
  }
  
  # Extract vectors for calculation to make formula cleaner
  g_counts <- gem_df[[graft_col]]
  h_counts <- gem_df[[host_col]]
  calls <- gem_df[[call_col]]
  
  # Perform Calculations
  gem_df$GraftDiff <- g_counts - h_counts
  gem_df$HostDiff <- h_counts - g_counts
  
  # Calculate total reads
  gem_df$totalReads <- g_counts + h_counts
  
  # Calculate percentage (handle potential division by zero if totalReads is 0)
  gem_df$percentHost <- ifelse(gem_df$totalReads > 0, h_counts / gem_df$totalReads, 0)
  
  # Log transform (log(0) is -Inf, usually acceptable in R, or add pseudocount if needed)
  gem_df$totalReadsLog <- log(gem_df$totalReads)
  
  # Assign standardized species labels based on dynamic input values
  assign_col_name <- paste0("AssignedSpecies_", suffix)
  gem_df[[assign_col_name]] <- ifelse(calls == graft_label, "Graft",
                                      ifelse(calls == host_label, "Host",
                                             ifelse(calls == multiplet_label, "Multiplet", NA)))
  
  return(gem_df)
}





#' @title Calculate Classification Metrics
#' @description This function calculates precision, recall, F1-score, and other classification metrics for each class and macro values.
#'
#' @param df A data frame containing the data to be processed.
#' @param gold_standard_col A string specifying the name of the gold standard column.
#' @param test_col A string specifying the name of the test column.
#'
#' @return A list containing precision, recall, F1-score, accuracy, kappa, and macro metrics.
#' @examples
#' \dontrun{
#' df <- data.frame(AssignedSpecies1 = sample(c("Human", "Mouse", "Multiplet"), 100, replace = TRUE),
#'                  AssignedSpecies2 = sample(c("Human", "Mouse", "Multiplet"), 100, replace = TRUE))
#' metrics <- calculate_metrics(df, "AssignedSpecies1", "AssignedSpecies2")
#' print(metrics)
#' }
#' @author Microsoft Copilot prompted by Amy Olex
#' @date 2025-02-05
#' @version 1.0
#' @keywords classification, metrics, precision, recall, F1, accuracy, kappa
#' @export
calculate_metrics <- function(df, gold_standard_col, test_col) {
  # Extract the gold standard and test columns
  gold_standard <- df[[gold_standard_col]]
  test <- df[[test_col]]

  # Create a confusion matrix
  cm <- confusionMatrix(as.factor(test), as.factor(gold_standard))

  # Calculate metrics for each class
  class_metrics <- data.frame(
    Class = rownames(cm$byClass),
    Precision = cm$byClass[, "Pos Pred Value"],
    Recall = cm$byClass[, "Sensitivity"],
    F1 = cm$byClass[, "F1"]
  )

  # Calculate macro metrics
  macro_precision <- mean(class_metrics$Precision, na.rm = TRUE)
  macro_recall <- mean(class_metrics$Recall, na.rm = TRUE)
  macro_f1 <- mean(class_metrics$F1, na.rm = TRUE)
  accuracy <- cm$overall["Accuracy"]
  kappa <- cm$overall["Kappa"]

  # Return metrics as a list
  return(list(
    class_metrics = class_metrics,
    macro_precision = macro_precision,
    macro_recall = macro_recall,
    macro_f1 = macro_f1,
    accuracy = accuracy,
    kappa = kappa
  ))
}



#' @title Evaluate Classification
#' @description This function evaluates classification performance by calculating metrics for each test column against a gold standard column.
#'
#' @param df A data frame containing the data to be processed.
#' @param gold_standard_col A string specifying the name of the gold standard column.
#'
#' @return A data frame containing precision, recall, F1-score, accuracy, kappa, and macro metrics for each test column.
#' @examples
#' \dontrun{
#' df <- data.frame(AssignedSpecies1 = sample(c("Human", "Mouse", "Multiplet"), 100, replace = TRUE),
#'                  AssignedSpecies2 = sample(c("Human", "Mouse", "Multiplet"), 100, replace = TRUE),
#'                  AssignedSpecies3 = sample(c("Human", "Mouse", "Multiplet"), 100, replace = TRUE))
#' metrics_df <- evaluate_classification(df, "AssignedSpecies1")
#' print(metrics_df)
#' }
#' @author Microsoft Copilot prompted by Amy Olex
#' @date 2025-02-05
#' @version 1.0
#' @keywords classification, evaluation, metrics, precision, recall, F1, accuracy, kappa
#' @export
evaluate_classification <- function(df, gold_standard_col) {
  # Extract columns with "AssignedSpecies" as the first part of the name
  assigned_species_cols <- grep("^AssignedSpecies", names(df), value = TRUE)

  # Check if the gold standard column is in the extracted columns
  if (!(gold_standard_col %in% assigned_species_cols)) {
    stop("Specified gold standard column is not in the extracted 'AssignedSpecies' columns.")
  }

  # Initialize a list to store results
  results <- list()

  # Loop through the test columns
  for (test_col in assigned_species_cols) {
    if (test_col != gold_standard_col) {
      metrics <- calculate_metrics(df, gold_standard_col, test_col)
      metrics$test_col <- test_col
      results[[test_col]] <- metrics
    }
  }

  # Convert results to a dataframe
  results_df <- do.call(rbind, lapply(results, function(x) {
    data.frame(
      Test_Column = x$test_col,
      Class = x$class_metrics$Class,
      Precision = x$class_metrics$Precision,
      Recall = x$class_metrics$Recall,
      F1 = x$class_metrics$F1,
      Macro_Precision = x$macro_precision,
      Macro_Recall = x$macro_recall,
      Macro_F1 = x$macro_f1,
      Accuracy = x$accuracy,
      Kappa = x$kappa
    )
  }))

  return(results_df)
}






#' @title Intersect and Subset Classification Data
#' @description Finds the intersection of barcodes between two data frames based on their row names and subsets them to include only shared cells.
#' @param df1 First data frame (barcodes must be in row.names).
#' @param df2 Second data frame (barcodes must be in row.names).
#' @return A list containing the two subsetted data frames (named df1 and df2).
#' @examples
#' \dontrun{
#' # Ensure data is read with row.names, e.g., read.csv(..., row.names = 1)
#' data1 <- read.csv("xenocell_classification_sampleA.csv", row.names = 1)
#' data2 <- read.csv("other_classification_sampleA.csv", row.names = 1)
#' result <- intersect_and_subset(data1, data2)
#' data1_clean <- result$df1
#' data2_clean <- result$df2
#' }
#' @export
intersect_and_subset <- function(df1, df2) {
  
  # Find common barcodes using row names
  common_barcodes <- intersect(rownames(df1), rownames(df2))
  
  # Report the intersection size
  message(paste("Found", length(common_barcodes), "common barcodes between the two datasets."))
  
  if (length(common_barcodes) == 0) {
    warning("No common barcodes found. Returning empty data frames.")
  }
  
  # Subset both data frames to keep only common barcodes
  # We subset by the character vector of row names directly
  df1_subset <- df1[common_barcodes, , drop = FALSE]
  df2_subset <- df2[common_barcodes, , drop = FALSE]
  
  # Return a list containing both processed data frames
  return(list(df1 = df1_subset, df2 = df2_subset))
}
