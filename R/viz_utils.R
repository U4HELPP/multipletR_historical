#' @title Visualization Functions for MultipletR Analysis
#' @description This script contains functions for generating plots and visualizations
#'              to display single cell data before and after removing multiplet GEMs.
#' @details The functions in this script are designed to help users visualize the impact
#'          of multiplet removal on their single cell datasets. Functions include:
#'          - Plotting raw data
#'          - Visualizing cleaned data
#'          - Comparing pre- and post-cleaning visualizations
#' @author Amy Olex
#' @date 2025-01-30
#' @version 1.0
#' @keywords visualization, plotting, single cell, multiplet GEMs
#'
#'

library("dplyr")
library("ggplot2")
library("reshape2")
library("gridExtra")


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

  return(gem_df)
}




#' @title Linear Percent Scatter Plot
#' @description This function creates a scatter plot of the percent of total reads that map to a species versus the total library size (UMI count) for each single cell. The axes and colors can be customized.
#' @param A data frame that has been processed by `prep_gem_counts()`, containing read counts with columns `barcode`, `GRCh38`, `mm10`, `call`, and additional columns added by `prep_gem_counts()`.
#' @param Xdata A string specifying the column name for the x-axis data. Default is "totalReads".
#' @param Ydata A string specifying the column name for the y-axis data. Default is "percentMouse".
#' @param color A string specifying the column name for the color grouping. Default is "call".
#' @param title A string specifying the title of the plot. Default is "Linear Percent Scatter Plot".
#' @param umi_cutoff A numeric value specifying the UMI cutoff for a vertical dashed line. Default is NA (no line).
#' @param species_cutoff A numeric value specifying the species percent cutoff for horizontal dashed lines. Default is NA (no lines). If greater than 1, it will be divided by 100.
#' @param Xaxislab A string specifying the label for the x-axis. Default is "Total UMI".
#' @param Yaxislab A string specifying the label for the y-axis. Default is "Percent Mouse".
#' @param colormapping A named vector specifying the colors for the different groups. Default is NA, which uses predefined colors.
#' @return A ggplot object representing the scatter plot.
#' @examples
#' \dontrun{
#' gem_df <- data.frame(barcode = c("AAACCAAAGCCATGCG-1", "AAACCCGCAATACTCT-1", "AAACGAATCAATGTGT-1"), GRCh38 = c(100, 200, 300), mm10 = c(50, 10, 250), call = c("GRCh38", "GRCh38", "Multiplet"))
#' gem_df <- prep_gem_counts(gem_df)
#' plot <- linear_percent_scatter(gem_df)
#' print(plot)
#' }
#' @export
linear_percent_scatter <- function(gem_df, Xdata = "totalReads", Ydata = "percentMouse", color = "call",
                                   title = "Linear Percent Scatter Plot", umi_cutoff = NA, species_cutoff = NA,
                                   Xaxislab = "Total UMI", Yaxislab = "Percent Mouse", colormapping = NA){
  if(all(is.na(colormapping))){
    colormapping <- c("GRCh38" = alpha("royalblue", 0.5), "mm10" = alpha("forestgreen", 0.3), "Multiplet" = alpha("gold", 1))
  }
  p1 <- ggplot(gem_df, aes_string(x = Xdata, y = Ydata, color = color)) +
    geom_point() +  # Add points
    scale_color_manual(values = colormapping) +  # Specify colors
    labs(title = title,
         x = Xaxislab,
         y = Yaxislab) +
    theme_minimal()  # Optional: use a minimal theme

  if(!is.na(umi_cutoff)){
    p1 <- p1 + geom_vline(xintercept = umi_cutoff, linetype = "dashed", color = "grey", size = .75)
  }
  if(!is.na(species_cutoff)){
    if(species_cutoff>1){
      species_cutoff <- species_cutoff/100
    }
    p1 <- p1 + geom_hline(yintercept = species_cutoff, linetype = "dashed", color = "darkgrey", size = .75) +
                geom_hline(yintercept = 1-species_cutoff, linetype = "dashed", color = "darkgrey", size = .75)
  }

  return(p1)
}



#' @title Semi-Log Percent Scatter Plot
#' @description This function creates a scatter plot of the percent of total reads that map to a species versus the log-transformed total library size (UMI count) for each single cell. The axes and colors can be customized.
#' @param gem_df A data frame that has been processed by `prep_gem_counts()`, containing read counts with columns `barcode`, `GRCh38`, `mm10`, `call`, and additional columns added by `prep_gem_counts()`.
#' @param Xdata A string specifying the column name for the x-axis data. Default is "totalReadsLog".
#' @param Ydata A string specifying the column name for the y-axis data. Default is "percentMouse".
#' @param color A string specifying the column name for the color grouping. Default is "call".
#' @param title A string specifying the title of the plot. Default is "Semi-Log Percent Scatter Plot".
#' @param umi_cutoff A numeric value specifying the UMI cutoff for a vertical dashed line. Default is NA (no line).
#' @param species_cutoff A numeric value specifying the species percent cutoff for horizontal dashed lines. Default is NA (no lines). If greater than 1, it will be divided by 100.
#' @param Xaxislab A string specifying the label for the x-axis. Default is "Total UMI".
#' @param Yaxislab A string specifying the label for the y-axis. Default is "Percent Mouse".
#' @return A ggplot object representing the scatter plot.
#' @examples
#' \dontrun{
#' gem_df <- data.frame(barcode = c("AAACCAAAGCCATGCG-1", "AAACCCGCAATACTCT-1", "AAACGAATCAATGTGT-1"), GRCh38 = c(100, 200, 300), mm10 = c(50, 10, 250), call = c("GRCh38", "GRCh38", "Multiplet"))
#' gem_df <- prep_gem_counts(gem_df)
#' plot <- semilog_percent_scatter(gem_df)
#' print(plot)
#' }
#' @export
semilog_percent_scatter <- function(gem_df, Xdata = "totalReadsLog", Ydata = "percentMouse", color = "call",
                                   title = "Semi-Log Percent Scatter Plot", umi_cutoff = NA, species_cutoff = NA,
                                   Xaxislab = "Total UMI", Yaxislab = "Percent Mouse"){

  p1 <- ggplot(gem_df, aes_string(x = Xdata, y = Ydata, color = color)) +
    geom_point() +  # Add points
    scale_color_manual(values = c("GRCh38" = alpha("royalblue", 0.5), "mm10" = alpha("forestgreen", 0.3), "Multiplet" = alpha("gold", 1))) +  # Specify colors
    labs(title = title,
         x = Xaxislab,
         y = Yaxislab) +
    theme_minimal()  # Optional: use a minimal theme

  if(!is.na(umi_cutoff)){
    p1 <- p1 + geom_vline(xintercept = umi_cutoff, linetype = "dashed", color = "grey", size = .75)
  }
  if(!is.na(species_cutoff)){
    if(species_cutoff>1){
      species_cutoff <- species_cutoff/100
    }
    p1 <- p1 + geom_hline(yintercept = species_cutoff, linetype = "dashed", color = "darkgrey", size = .75) +
      geom_hline(yintercept = 1-species_cutoff, linetype = "dashed", color = "darkgrey", size = .75)
  }

  return(p1)
}




############# The following function is still on its way to being converted to individual plotting functions as above.
multiplet_plots <- function(gem_df, title, save_dir){

  gem_df <- prep_gem_counts(gem_df)

  log_max = max(gem_df$totalReadsLog)
  gem_df$graft_cutpoint_log <- sapply(gem_df$totalReadsLog, find_log_cutoff, y1=.25, y2=.10, x_max=log_max)
  gem_df$host_cutpoint_log <- sapply(gem_df$totalReadsLog, find_log_cutoff, y1=.75, y2=.90, x_max=log_max)
  gem_df$MESSY_Call <- assign_new_call(gem_df[,c("percentMouse", "graft_cutpoint_log", "host_cutpoint_log")])

  graft_line = get_custom_slope(y1=.25, y2=.10, x_max=log_max)
  host_line = get_custom_slope(y1=.75, y2=.90, x_max=log_max)

  p1 <- ggplot(gem_df, aes(x = totalReads, y = percentMouse, color = x10X_Call)) +
    geom_point() +  # Add points
    scale_color_manual(values = c("GRCh38" = alpha("royalblue", 0.5), "mm10" = alpha("forestgreen", 0.3), "Multiplet" = alpha("gold", 1))) +  # Specify colors
    geom_vline(xintercept = 500, linetype = "dashed", color = "grey", size = .75) +
    geom_hline(yintercept = .10, linetype = "dashed", color = "darkgrey", size = .75) +
    geom_hline(yintercept = .90, linetype = "dashed", color = "darkgrey", size = .75) +
    #geom_hline(yintercept = .25, linetype = "dashed", color = "black", size = .75) +
    #geom_hline(yintercept = .75, linetype = "dashed", color = "black", size = .75) +
    #geom_abline(slope = graft_line["slope"], intercept = graft_line["y-intercept"], color = "red", linetype = "dashed", size = .75) +  # Add sloped line
    #geom_abline(slope = host_line["slope"], intercept = host_line["y-intercept"], color = "red", linetype = "dashed", size = .75) +  # Add sloped line

    labs(title = paste0("10X GEM Multiplet Calls: ", title),
         x = "Total Number of Reads in Cell",
         y = "Percent Mouse") +
    theme_minimal()  # Optional: use a minimal theme

  p2 <- ggplot(gem_df, aes(x = totalReadsLog, y = percentMouse, color = x10X_Call)) +
    geom_point() +  # Add points
    scale_color_manual(values = c("GRCh38" = alpha("royalblue", 0.5), "mm10" = alpha("forestgreen", 0.3), "Multiplet" = alpha("gold", 1))) +  # Specify colors
    geom_vline(xintercept = log(500), linetype = "dashed", color = "grey", size = .75) +
    geom_hline(yintercept = .10, linetype = "dashed", color = "darkgrey", size = .75) +
    geom_hline(yintercept = .90, linetype = "dashed", color = "darkgrey", size = .75) +
    #geom_hline(yintercept = .25, linetype = "dashed", color = "black", size = .75) +
    #geom_hline(yintercept = .75, linetype = "dashed", color = "black", size = .75) +

    labs(title = paste0("10X GEM Multiplet Calls: ", title),
         x = "Log of Total Number of Reads in Cell",
         y = "Percent Mouse") +
    theme_minimal()  # Optional: use a minimal theme

  p11 <- ggplot(gem_df, aes(x = totalReads, y = percentMouse, color = MESSY_Call)) +
    geom_point() +  # Add points
    scale_color_manual(values = c("GRCh38" = alpha("royalblue", 0.5), "mm10" = alpha("red3", 0.5), "Multiplet" = alpha("purple", 0.5))) +  # Specify colors
    geom_vline(xintercept = 500, linetype = "dashed", color = "grey", size = .75) +
    geom_hline(yintercept = .10, linetype = "dashed", color = "darkgrey", size = .75) +
    geom_hline(yintercept = .90, linetype = "dashed", color = "darkgrey", size = .75) +
    #geom_hline(yintercept = .25, linetype = "dashed", color = "black", size = .75) +
    #geom_hline(yintercept = .75, linetype = "dashed", color = "black", size = .75) +
    labs(title = paste0("Olex Multiplet Calls: ", title),
         x = "Total Number of Reads in Cell",
         y = "Percent Mouse") +
    theme_minimal()  # Optional: use a minimal theme

  p22 <- ggplot(gem_df, aes(x = totalReadsLog, y = percentMouse, color = MESSY_Call)) +
    geom_point() +  # Add points
    scale_color_manual(values = c("GRCh38" = alpha("royalblue", 0.5), "mm10" = alpha("red3", 0.5), "Multiplet" = alpha("purple", 0.5))) +  # Specify colors
    geom_vline(xintercept = log(500), linetype = "dashed", color = "grey", size = .75) +
    geom_hline(yintercept = .10, linetype = "dashed", color = "darkgrey", size = .75) +
    geom_hline(yintercept = .90, linetype = "dashed", color = "darkgrey", size = .75) +
    #geom_hline(yintercept = .25, linetype = "dashed", color = "black", size = .75) +
    #geom_hline(yintercept = .75, linetype = "dashed", color = "black", size = .75) +
    geom_abline(slope = graft_line["slope"], intercept = graft_line["y-intercept"], color = "black", linetype = "dashed", size = .75) +  # Add sloped line
    geom_abline(slope = host_line["slope"], intercept = host_line["y-intercept"], color = "black", linetype = "dashed", size = .75) +  # Add sloped line

    labs(title = paste0("Olex Multiplet Calls: ", title),
         x = "Log of Total Number of Reads in Cell",
         y = "Percent Mouse") +
    theme_minimal()  # Optional: use a minimal theme

  gem_df_melt <- as.data.frame(reshape2::melt(gem_df[,c("GRCh38_count", "mm10_count","x10X_Call")]))
  names(gem_df_melt) <- c("CellCall", "Aligned2Ref", "ReadCount")

  # Visualize the number UMIs/transcripts per cell
  p3 <- gem_df_melt[gem_df_melt$CellCall=="GRCh38",] %>%
    ggplot(aes(color=Aligned2Ref, x=log(ReadCount), fill= Aligned2Ref, color= Aligned2Ref)) +
    geom_density(alpha = 0.2) +
    scale_color_manual(values = c("GRCh38_count" = "royalblue", "mm10_count" = "forestgreen")) +  # Specify colors
    scale_fill_manual(values = c("GRCh38_count" = "royalblue", "mm10_count" = "forestgreen")) +  # Specify colors
    #geom_histogram()
    theme_classic() +
    ylab("GRCh38 Cell Count Density") +
    xlab("Log Read Count") +
    geom_vline(xintercept = 25)


  p4 <- gem_df_melt[gem_df_melt$CellCall=="mm10",] %>%
    ggplot(aes(color=Aligned2Ref, x=log(ReadCount), fill= Aligned2Ref, color= Aligned2Ref)) +
    geom_density(alpha = 0.2) +
    scale_color_manual(values = c("GRCh38_count" = "royalblue", "mm10_count" = "forestgreen")) +  # Specify colors
    scale_fill_manual(values = c("GRCh38_count" = "royalblue", "mm10_count" = "forestgreen")) +  # Specify colors
    #geom_histogram()
    theme_classic() +
    ylab("mm10 Cell Count Density") +
    xlab("Log Read Count") +
    geom_vline(xintercept = 25)

  p5 <- gem_df_melt[gem_df_melt$CellCall=="Multiplet",] %>%
    ggplot(aes(color=Aligned2Ref, x=log(ReadCount), fill= Aligned2Ref, color= Aligned2Ref)) +
    geom_density(alpha = 0.2) +
    scale_color_manual(values = c("GRCh38_count" = "royalblue", "mm10_count" = "forestgreen")) +  # Specify colors
    scale_fill_manual(values = c("GRCh38_count" = "royalblue", "mm10_count" = "forestgreen")) +  # Specify colors
    #geom_histogram()
    theme_classic() +
    ylab("Multiplet Cell Count Density") +
    xlab("Log Read Count") +
    geom_vline(xintercept = 25)


  gem_df_melt2 <- as.data.frame(reshape2::melt(gem_df[,c("GRCh38_count","mm10_count","MESSY_Call")]))
  names(gem_df_melt2) <- c("OlexCellCall", "Aligned2Ref", "ReadCount")

  # Visualize the number UMIs/transcripts per cell
  p33 <- gem_df_melt2[gem_df_melt2$OlexCellCall=="GRCh38",] %>%
    ggplot(aes(color=Aligned2Ref, x=log(ReadCount), fill= Aligned2Ref, color= Aligned2Ref)) +
    geom_density(alpha = 0.2) +
    scale_color_manual(values = c("GRCh38_count" = "royalblue", "mm10_count" = "red3")) +  # Specify colors
    scale_fill_manual(values = c("GRCh38_count" = "royalblue", "mm10_count" = "red3")) +  # Specify colors
    #geom_histogram()
    theme_classic() +
    ylab("Human Cell Count Density") +
    xlab("Log Read Count") +
    geom_vline(xintercept = 25)


  p44 <- gem_df_melt2[gem_df_melt2$OlexCellCall=="mm10",] %>%
    ggplot(aes(color=Aligned2Ref, x=log(ReadCount), fill= Aligned2Ref, color= Aligned2Ref)) +
    geom_density(alpha = 0.2) +
    scale_color_manual(values = c("GRCh38_count" = "royalblue", "mm10_count" = "red3")) +  # Specify colors
    scale_fill_manual(values = c("GRCh38_count" = "royalblue", "mm10_count" = "red3")) +  # Specify colors
    #geom_histogram()
    theme_classic() +
    ylab("Mouse Cell Count Density") +
    xlab("Log Read Count") +
    geom_vline(xintercept = 25)

  p55 <- gem_df_melt2[gem_df_melt2$OlexCellCall=="Multiplet",] %>%
    ggplot(aes(color=Aligned2Ref, x=log(ReadCount), fill= Aligned2Ref, color= Aligned2Ref)) +
    geom_density(alpha = 0.2) +
    scale_color_manual(values = c("GRCh38_count" = "royalblue", "mm10_count" = "red3")) +  # Specify colors
    scale_fill_manual(values = c("GRCh38_count" = "royalblue", "mm10_count" = "red3")) +  # Specify colors
    #geom_histogram()
    theme_classic() +
    ylab("Multiplet Cell Count Density") +
    xlab("Log Read Count") +
    geom_vline(xintercept = 25)

  png(filename = paste0("MultipletCalls_", title, ".png"), width = 1000, height = 1300)
  grid.arrange(p1, p11, p2, p22, p3, p33, p4, p44, p5, p55, ncol = 2)  # Arrange plots in 1 row and 2 columns
  dev.off()

  gem_df$barcode <- row.names(gem_df)
  for_gem_file <- gem_df[,c("barcode", "GRCh38_count", "mm10_count", "MESSY_Call")]
  names(for_gem_file) <- c("barcode", "GRCh38", "mm10", "call")

  write.csv(gem_df, file=paste0("MultipletCalls_Matrix_", title, ".csv"), quote = FALSE, row.names = FALSE)
  write.csv(for_gem_file, file=paste0(save_dir, "gem_classification.csv"), quote = FALSE, row.names = FALSE)

  # write out Loupe Annotations
  write.csv(gem_df[,c("barcode","x10X_Call")], file=paste0(title, "_10x_gem_class_Annotations.csv"), quote = FALSE, row.names = FALSE)
  write.csv(gem_df[,c("barcode","MESSY_Call")], file=paste0(title, "_MESSY_class_Annotations.csv"), quote = FALSE, row.names = FALSE)




  #return(gem_df)

}
