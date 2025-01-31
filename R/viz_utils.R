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
#' @importFrom ggplot2 ggplot aes_string geom_point scale_color_manual labs theme_minimal geom_vline geom_hline
#' @importFrom scales alpha
#' @export
linear_percent_scatter <- function(gem_df, Xdata = "totalReads", Ydata = "percentMouse", color = "AssignedSpecies",
                                   title = "Linear Percent Scatter Plot", umi_cutoff = NA, species_cutoff = NA,
                                   Xaxislab = "Total UMI", Yaxislab = "Percent Mouse", colormapping = NA){
  if(all(is.na(colormapping))){
    colormapping <- c("Human" = alpha("royalblue", 0.5), "Mouse" = alpha("forestgreen", 0.3), "Multiplet" = alpha("gold", 1))
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
#' @param colormapping A named vector specifying the colors for the different groups. Default is NA, which uses predefined colors.
#' @return A ggplot object representing the scatter plot.
#' @examples
#' \dontrun{
#' gem_df <- data.frame(barcode = c("AAACCAAAGCCATGCG-1", "AAACCCGCAATACTCT-1", "AAACGAATCAATGTGT-1"), GRCh38 = c(100, 200, 300), mm10 = c(50, 10, 250), call = c("GRCh38", "GRCh38", "Multiplet"))
#' gem_df <- prep_gem_counts(gem_df)
#' plot <- semilog_percent_scatter(gem_df)
#' print(plot)
#' }
#' @importFrom ggplot2 ggplot aes_string geom_point scale_color_manual labs theme_minimal geom_vline geom_hline
#' @importFrom scales alpha
#' @export
semilog_percent_scatter <- function(gem_df, Xdata = "totalReadsLog", Ydata = "percentMouse", color = "AssignedSpecies",
                                   title = "Semi-Log Percent Scatter Plot", umi_cutoff = NA, species_cutoff = NA,
                                   Xaxislab = "Total UMI (Natural Log)", Yaxislab = "Percent Mouse", colormapping = NA){

  if(all(is.na(colormapping))){
    colormapping <- c("Human" = alpha("royalblue", 0.5), "Mouse" = alpha("forestgreen", 0.3), "Multiplet" = alpha("gold", 1))
  }

  p2 <- ggplot(gem_df, aes_string(x = Xdata, y = Ydata, color = color)) +
    geom_point() +  # Add points
    scale_color_manual(values = colormapping) +  # Specify colors
    labs(title = title,
         x = Xaxislab,
         y = Yaxislab) +
    theme_minimal()  # Optional: use a minimal theme

  if(!is.na(umi_cutoff)){
    p2 <- p2 + geom_vline(xintercept = log(umi_cutoff), linetype = "dashed", color = "grey", size = .75)
  }
  if(!is.na(species_cutoff)){
    if(species_cutoff>1){
      species_cutoff <- species_cutoff/100
    }
    p2 <- p2 + geom_hline(yintercept = species_cutoff, linetype = "dashed", color = "darkgrey", size = .75) +
      geom_hline(yintercept = 1-species_cutoff, linetype = "dashed", color = "darkgrey", size = .75)
  }

  return(p2)
}



#' @title Species Density Plot
#' @description This function creates a density plot of aligned read counts for a specified assigned species. The data frame must be processed by `prep_gem_counts()` first.
#' @param gem_df A data frame that has been processed by `prep_gem_counts()`, containing read counts with columns `barcode`, `GRCh38`, `mm10`, `call`, and additional columns added by `prep_gem_counts()`.
#' @param assigned_species A string specifying the assigned species to plot. Default is "Human".
#' @param refGenome_linecol A named vector specifying the line colors for the reference genomes. Default is NA, which uses predefined colors.
#' @param refGenome_fillcol A named vector specifying the fill colors for the reference genomes. Default is NA, which uses predefined colors.
#' @param title A string specifying the title of the plot. Default is "Assigned Species: Human".
#' @param Xaxislab A string specifying the label for the x-axis. Default is "Aligned Read Counts (Natural Log)".
#' @param Yaxislab A string specifying the label for the y-axis. Default is "Cell Density".
#' @return A ggplot object representing the density plot.
#' @examples
#' \dontrun{
#' gem_df <- data.frame(barcode = c("AAACCAAAGCCATGCG-1", "AAACCCGCAATACTCT-1", "AAACGAATCAATGTGT-1"), GRCh38 = c(100, 200, 300), mm10 = c(50, 10, 250), call = c("GRCh38", "GRCh38", "Multiplet"))
#' gem_df <- prep_gem_counts(gem_df)
#' plot <- species_density_plot(gem_df, assigned_species = "Human")
#' print(plot)
#' }
#' @importFrom ggplot2 ggplot aes geom_density scale_color_manual scale_fill_manual theme_classic labs geom_vline
#' @importFrom reshape2 melt
#' @importFrom dplyr %>%
#' @importFrom scales alpha
#' @export
species_density_plot <- function(gem_df, assigned_species = "Human",
                                  refGenome_linecol = NA, refGenome_fillcol = NA,
                                  title = paste("Assigned Species:", assigned_species),
                                  Xaxislab = "Aligned Read Counts (Natural Log)",
                                  Yaxislab = "Cell Density"){

  gem_df_melt <- as.data.frame(reshape2::melt(gem_df[,c("GRCh38", "mm10","AssignedSpecies")]))
  names(gem_df_melt) <- c("AssignedSpecies", "RefGenome", "AlignedCount")

  if(all(is.na(refGenome_linecol))){
    refGenome_linecol <- c("GRCh38" = "royalblue", "mm10" = "forestgreen")
  }
  if(all(is.na(refGenome_fillcol))){
    refGenome_fillcol <- c("GRCh38" = "royalblue", "mm10" = "forestgreen")
  }

  # Visualize the number UMIs/transcripts per cell
  p3 <- gem_df_melt[gem_df_melt$AssignedSpecies==assigned_species,] %>%
    ggplot(aes(color=RefGenome, x=log(AlignedCount), fill= RefGenome, color= RefGenome)) +
    geom_density(alpha = 0.2) +
    scale_color_manual(values = refGenome_linecol) +  # Specify colors
    scale_fill_manual(values = refGenome_fillcol) +  # Specify colors
    theme_classic() +
    labs(title = title,
         x = Xaxislab,
         y = Yaxislab) +
    geom_vline(xintercept = 25)

  return(p3)
}


#' @title GEM Classification Summary
#' @description This function creates a summary plot of GEM classification, including linear and semi-log scatter plots and density plots for human, mouse, and multiplet assigned species. The data frame must be processed by `prep_gem_counts()` first.
#' @param gem_df A data frame that has been processed by `prep_gem_counts()`, containing read counts with columns `barcode`, `GRCh38`, `mm10`, `call`, and additional columns added by `prep_gem_counts()`.
#' @param title A string specifying the title of the summary plot. Default is "GEM_Classification_Summary".
#' @param umi_cutoff A numeric value specifying the UMI cutoff for vertical dashed lines in the scatter plots. Default is 500 UMI.
#' @param species_cutoff A numeric value specifying the species percent cutoff for horizontal dashed lines in the scatter plots. Default is .10. If greater than 1, it will be divided by 100.
#' @return A grid of ggplot objects representing the summary plot.
#' @examples
#' \dontrun{
#' gem_df <- data.frame(barcode = c("AAACCAAAGCCATGCG-1", "AAACCCGCAATACTCT-1", "AAACGAATCAATGTGT-1"), GRCh38 = c(100, 200, 300), mm10 = c(50, 10, 250), call = c("GRCh38", "GRCh38", "Multiplet"))
#' gem_df <- prep_gem_counts(gem_df)
#' summary_plot <- gem_classification_summary(gem_df)
#' print(summary_plot)
#' }
#' @importFrom ggplot2 ggplot aes_string geom_point scale_color_manual labs theme_minimal geom_vline geom_hline
#' @importFrom scales alpha
#' @importFrom gridExtra grid.arrange
#' @export
gem_classification_summary <- function(gem_df, title="GEM_Classification_Summary", umi_cutoff=500, species_cutoff=.10){

  p1 <- linear_percent_scatter(gem_df, umi_cutoff = umi_cutoff, species_cutoff = species_cutoff)
  p2 <- semilog_percent_scatter(gem_df, umi_cutoff = umi_cutoff, species_cutoff = species_cutoff)
  p3_human <- species_density_plot(gem_df, assigned_species = "Human")
  p3_mouse <- species_density_plot(gem_df, assigned_species = "Mouse")
  p3_multiplet <- species_density_plot(gem_df, assigned_species = "Multiplet")
  summary_plot <- grid.arrange(p1, p2, p3_human, p3_mouse, p3_multiplet, ncol = 2)

  return(summary_plot)
}




