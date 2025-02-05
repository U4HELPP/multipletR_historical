#' @title Prepare GEM Counts for Analysis and Plotting
#' @description This function takes a data frame of read counts from a 10X GEM file and adds additional columns for analysis and plotting.
#' @param gem_df A data frame containing read counts with columns `barcode`, `GRCh38`, `mm10` and `call`.
#' @return A data frame with the following additional columns:
#' \itemize{
#'   \item \code{HumanDiff}: Difference between human (GRCh38) and mouse (mm10) read counts.
#'   \item \code{MouseDiff}: Difference between mouse (mm10) and human (GRCh38) read counts.
#'   \item \code{percentMouse}: Proportion of reads that map to mouse (mm10).
#'   \item \code{totalReads}: Total number of reads (sum of GRCh38 and mm10).
#'   \item \code{totalReadsLog}: Logarithm of the total number of reads.
#' }
#' @examples
#' \dontrun{
#' gem_df <- data.frame(barcode = c("AAACCAAAGCCATGCG-1", "AAACCCGCAATACTCT-1", "AAACGAATCAATGTGT-1"), GRCh38 = c(100, 200, 300), mm10 = c(50, 10, 250), call = c("GRCh38", "GRCh38", "Multiplet"))
#' result <- prep_gem_counts(gem_df)
#' print(result)
#' }
#' @export
prep_gem_counts <- function(gem_df){
  gem_df$HumanDiff <- gem_df$GRCh38 - gem_df$mm10
  gem_df$MouseDiff <- gem_df$mm10 - gem_df$GRCh38
  gem_df$percentMouse <- gem_df$mm10/(gem_df$mm10 + gem_df$GRCh38)
  gem_df$totalReads <- gem_df$mm10 + gem_df$GRCh38
  gem_df$totalReadsLog <- log(gem_df$mm10 + gem_df$GRCh38)
  gem_df$AssignedSpecies_10X <- ifelse(gem_df$call == "GRCh38", "Human",
                                   ifelse(gem_df$call == "mm10", "Mouse",
                                          ifelse(gem_df$call == "Multiplet", "Multiplet", NA)))

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
