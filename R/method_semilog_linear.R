#' @title Semi-Log Linear Threshold Method Functions for GEM Analysis
#' @description This script contains functions for calculating and visualizing linear semi-logarithmic thresholds for GEM classification.
#' @details The functions in this script are designed to help users analyze and visualize GEM data using semi-logarithmic methods. Functions include:
#'          - Calculating log cutoffs
#'          - Determining custom slopes
#'          - Assigning new species classifications
#'          - Finding semi-log thresholds
#'          - Creating semi-log scatter plots with custom threshold lines
#' @author Amy Olex
#' @date 2025-02-05
#' @version 1.0
#' @keywords semi-log, GEM analysis, threshold, scatter plot, species classification


#' @title Find Log Cutoff
#' @description This function calculates the predicted y-value on a semi-logarithmic scale based on the given x-value and two reference points.
#'
#' @param x A numeric value representing the x-coordinate (in logarithmic scale) for which the y-value is to be predicted.
#' @param y1 A numeric value representing the y-coordinate of the first reference point.
#' @param y2 A numeric value representing the y-coordinate of the second reference point.
#' @param x_max A numeric value representing the maximum x-coordinate (in logarithmic scale).
#'
#' @return A numeric value representing the predicted y-value.
#' @examples
#' \dontrun{
#' y_pred <- find_log_cutoff(x = log(750), y1 = 0.25, y2 = 0.10, x_max = log(1000))
#' print(y_pred)
#' }
#' @export
find_log_cutoff <- function(x, y1, y2, x_max){
  x1 = log(500)
  x2 = x_max
  m = (y1-y2)/(log(500)-x_max)
  b = y1 - m*x1
  y_pred = (m*x) + b
  return(y_pred)
}


#' @title Get Custom Slope
#' @description This function calculates the slope and y-intercept for a line passing through two points on a semi-logarithmic scale.
#'
#' @param y1 A numeric value representing the y-coordinate of the first point.
#' @param y2 A numeric value representing the y-coordinate of the second point.
#' @param x_max A numeric value representing the maximum x-coordinate (in logarithmic scale).
#'
#' @return A named vector with the slope (`"slope"`) and y-intercept (`"y-intercept"`) of the line.
#' @examples
#' \dontrun{
#' slope_intercept <- get_custom_slope(y1 = 0.25, y2 = 0.10, x_max = log(1000))
#' print(slope_intercept)
#' }
#' @export
get_custom_slope <- function(y1, y2, x_max){
  x1 = log(500)
  x2 = x_max
  m = (y1-y2)/(log(500)-x_max)
  b = y1 - m*x1
  return(c("slope"=m, "y-intercept"=b))
}

#' @title Assign New Species Call
#' @description This function assigns new species classifications based on mouse percent, graft cutoff, and host cutoff values.
#'
#' @param df A data frame containing three columns: mouse percent, graft cutoff, and host cutoff. The columns must be in this order.
#'
#' @return A vector of assigned classifications (`"GRCh38"`, `"mm10"`, or `"Multiplet"`) based on the input thresholds.
#' @examples
#' \dontrun{
#' df <- data.frame(mouse_percent = runif(100), graft_cut = runif(100), host_cut = runif(100))
#' calls <- assign_new_call(df)
#' }
#' @export
assign_new_call <- function(df){
  names(df) <- c("host_percent", "graft_cut", "host_cut")
  df$call <- NA
  for(i in 1:dim(df)[1]){
    #print(df[i,"host_percent"])
    if(df[i,"host_percent"] <= df[i,"graft_cut"])
    {df[i,"call"] = "GRCh38"}
    else if(df[i,"host_percent"] >= df[i,"host_cut"])
    {df[i,"call"] = "mm10"}
    else
    {df[i,"call"] = "Multiplet"}
  }
  return(df$call)
}


#' @title Find Semi-Log Threshold
#' @description This function calculates semi-logarithmic threshold values for graft and host populations and assigns new classifications based on these thresholds.
#'
#' @param gem_df A data frame containing the data to be processed, with a column `totalReadsLog` representing the logarithm of total reads.
#' @param graft_max A numeric value indicating the maximum graft threshold. Default is 0.25. If greater than 1, it will be divided by 100.
#' @param graft_min A numeric value indicating the minimum graft threshold. Default is 0.10. If greater than 1, it will be divided by 100.
#'
#' @return A data frame with additional columns `graft_cutpoint_log`, `host_cutpoint_log`, `SemiLog_Call`, and `SemiLogSpecies` representing the calculated thresholds and new classifications.
#' @examples
#' \dontrun{
#' gem_df <- data.frame(totalReadsLog = rnorm(100), percentMouse = runif(100))
#' gem_df <- find_semilog_threshold(gem_df, graft_max = 0.25, graft_min = 0.10)
#' }
#' @importFrom stats sapply
#' @export
find_semilog_threshold <- function(gem_df, graft_max=.25, graft_min=.10){

  if(graft_max>1){
    graft_max <- graft_max/100
  }
  if(graft_min>1){
    graft_min <- graft_min/100
  }

  log_max = max(gem_df$totalReadsLog)
  gem_df$graft_cutpoint_log <- sapply(gem_df$totalReadsLog, find_log_cutoff, y1=graft_max, y2=graft_min, x_max=log_max)
  gem_df$host_cutpoint_log <- sapply(gem_df$totalReadsLog, find_log_cutoff, y1=1-graft_max, y2=1-graft_min, x_max=log_max)
  gem_df$SemiLog_Call <- assign_new_call(gem_df[,c("percentMouse", "graft_cutpoint_log", "host_cutpoint_log")])
  gem_df$SemiLogSpecies <- ifelse(gem_df$SemiLog_Call == "GRCh38", "Human",
                                   ifelse(gem_df$SemiLog_Call == "mm10", "Mouse",
                                          ifelse(gem_df$SemiLog_Call == "Multiplet", "Multiplet", NA)))

  return(gem_df)
}



#' @title Percent Scatter Plot with Semi-Log Threshold
#' @description This function creates a semi-logarithmic scatter plot with custom threshold lines for graft and host populations.
#'
#' @param gem_df A data frame containing the data to be plotted.
#' @param graft_max A numeric value indicating the maximum graft threshold. Default is 0.25.
#' @param graft_min A numeric value indicating the minimum graft threshold. Default is 0.10.
#' @param ... Additional arguments passed to the `semilog_percent_scatter` function.
#'
#' @return A ggplot object representing the semi-logarithmic scatter plot with custom threshold lines.
#' @importFrom ggplot2 geom_abline
#' @importFrom scales alpha
#' @export
#'
#' @examples
#' \dontrun{
#' gem_df <- data.frame(totalReadsLog = rnorm(100), SemiLogSpecies = sample(c("Human", "Mouse", "Multiplet"), 100, replace = TRUE))
#' plot_semilog_threshold(gem_df, graft_max = 0.25, graft_min = 0.10)
#' }
plot_semilog_threshold <- function(gem_df, graft_max=.25, graft_min=.10, ...){

  if(graft_max>1){
    graft_max <- graft_max/100
  }
  if(graft_min>1){
    graft_min <- graft_min/100
  }

  # Capture the additional arguments
  args <- list(...)

  # Check if 'w' is missing and set a default value if it is
  if (!"colormapping" %in% names(args)) {
    args$colormapping <- c("Human" = alpha("royalblue", 0.5), "Mouse" = alpha("red3", 0.5), "Multiplet" = alpha("purple", 0.5))
  }

  if(!"color" %in% names(args)){
    args$color = "SemiLogSpecies"
  }

  # Create the initial semilog scatter plot.
  p <- do.call(semilog_percent_scatter, c(list(gem_df), args))

  # Calculate the slope of the threshold lines in semilog space.
  log_max = max(gem_df$totalReadsLog)
  graft_line = get_custom_slope(y1=graft_max, y2=graft_min, x_max=log_max)
  host_line = get_custom_slope(y1=1-graft_max, y2=1-graft_min, x_max=log_max)

  # Add the additional threshold lines to the plot.
  p <- p + geom_abline(slope = graft_line["slope"], intercept = graft_line["y-intercept"], color = "black", linetype = "dashed", size = .75) +  # Add sloped line
    geom_abline(slope = host_line["slope"], intercept = host_line["y-intercept"], color = "black", linetype = "dashed", size = .75)  # Add sloped line

  return(p)

}



## Also want to add in methods to plot the semilog threshold using the viz_utils plots.



gem_df$barcode <- row.names(gem_df)
for_gem_file <- gem_df[,c("barcode", "GRCh38_count", "mm10_count", "MESSY_Call")]
names(for_gem_file) <- c("barcode", "GRCh38", "mm10", "call")

write.csv(gem_df, file=paste0("MultipletCalls_Matrix_", title, ".csv"), quote = FALSE, row.names = FALSE)
write.csv(for_gem_file, file=paste0(save_dir, "gem_classification.csv"), quote = FALSE, row.names = FALSE)

# write out Loupe Annotations
write.csv(gem_df[,c("barcode","x10X_Call")], file=paste0(title, "_10x_gem_class_Annotations.csv"), quote = FALSE, row.names = FALSE)
write.csv(gem_df[,c("barcode","MESSY_Call")], file=paste0(title, "_MESSY_class_Annotations.csv"), quote = FALSE, row.names = FALSE)

