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
  gem_df$Assigned10XSpecies <- ifelse(gem_df$call == "GRCh38", "Human",
                                   ifelse(gem_df$call == "mm10", "Mouse",
                                          ifelse(gem_df$call == "Multiplet", "Multiplet", NA)))

  return(gem_df)
}
