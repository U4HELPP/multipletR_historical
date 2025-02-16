## Amy Olex
## 10.23.2024
## Exploring the doublet and mutiplet detection from CellRanger as we are getting odd results for the Pancrease data.

library("dplyr")
library("ggplot2")
library("reshape2")
library("gridExtra")

setwd("/lustre/home/harrell_lab/scRNASeq/exploratory_analyses/multiplet_exploration")

find_log_cutoff <- function(x, y1, y2, x_max){
  x1 = log(500)
  x2 = x_max
  m = (y1-y2)/(log(500)-x_max)
  b = y1 - m*x1
  y_pred = (m*x) + b
  return(y_pred)
}

get_custom_slope <- function(y1, y2, x_max){
  x1 = log(500)
  x2 = x_max
  m = (y1-y2)/(log(500)-x_max)
  b = y1 - m*x1
  return(c("slope"=m, "y-intercept"=b))
}

## df must be mouse percent, graft cutoff, host cutoff columns only
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

setwd("10X_gem_classifications/")
pc62 <- read.delim("VCU-PC-062_gem_classification.csv", header = TRUE, row.names = 1, sep=",")
dir.create("Analyses/02_VCU-PC-062_1803_CRv8_241016/outs/analysis_olexClass/", recursive = TRUE)
gem_df = pc62; title = "PC-062"; save_dir = "Analyses/02_VCU-PC-062_1803_CRv8_241016/outs/analysis_olexClass/"

multiplet_plots <- function(gem_df, title, save_dir){
  names(gem_df) <- c("GRCh38_count", "mm10_count", "x10X_Call")
  gem_df$HumanDiff <- gem_df$GRCh38 - gem_df$mm10 # Human counts, without mouse
  gem_df$MouseDiff <- gem_df$mm10 - gem_df$GRCh38 # Mouse counts, without human
  gem_df$percentMouse <- gem_df$mm10/(gem_df$mm10 + gem_df$GRCh38) # Percent mouse
  gem_df$totalReads <- gem_df$mm10 + gem_df$GRCh38 # Total counts
  gem_df$totalReadsLog <- log10(gem_df$mm10 + gem_df$GRCh38) # Log total counts
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
  
  p1
  p2
  
  # Mikhail and ChatGPT code
  summary(gem_df$totalReadsLog)
  quantile(gem_df$totalReadsLog, probs = c(0, 0.25, 0.5, 0.75, 1))
  
  # Evaluating the effect of the number of regions
  library(dplyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
  library(MASS)
  
  # Select data
  data <- gem_df[, c("percentMouse", "totalReadsLog", "x10X_Call")]
  
  fit_distribution <- function(values, lower_bound = TRUE, dist = "lnorm") {
    if (lower_bound) {
      values <- values[values <= 0.5]
    } else {
      values <- 1 - values[values > 0.5]
    }
    
    if (length(values) > 1) {
      if (dist == "gamma") {
        fit <- fitdistr(values, "gamma")
      } else if (dist == "lnorm") {
        fit <- fitdistr(values, "lognormal")
      } else if (dist == "norm") {
        fit <- fitdistr(values, "normal")
      } else {
        stop("Unsupported distribution")
      }
      return(fit)
    }
    return(NULL)
  }
  
  classify_point <- function(value, fit_0, fit_1, dist = "lnorm") {
    if (value <= 0.5) {
      if (dist == "gamma") {
        p_fit_0 <- dgamma(value, shape = fit_0$estimate["shape"], rate = fit_0$estimate["rate"])
      } else if (dist == "lnorm") {
        p_fit_0 <- dlnorm(value, meanlog = fit_0$estimate["meanlog"], sdlog = fit_0$estimate["sdlog"])
      } else if (dist == "norm") {
        p_fit_0 <- dnorm(value, mean = fit_0$estimate["mean"], sd = fit_0$estimate["sd"])
      }
      return(ifelse(p_fit_0 > 0.05, "human", "mixed"))
    } else {
      if (dist == "gamma") {
        p_fit_1 <- dgamma(1 - value, shape = fit_1$estimate["shape"], rate = fit_1$estimate["rate"])
      } else if (dist == "lnorm") {
        p_fit_1 <- dlnorm(1 - value, meanlog = fit_1$estimate["meanlog"], sdlog = fit_1$estimate["sdlog"])
      } else if (dist == "norm") {
        p_fit_1 <- dnorm(1 - value, mean = fit_1$estimate["mean"], sd = fit_1$estimate["sd"])
      }
      return(ifelse(p_fit_1 > 0.05, "mouse", "mixed"))
    }
  }
  
  results_list <- tibble(n_regions = integer(), human = integer(), mouse = integer(), mixed = integer())
  
  for (num_regions in 1:20) {
    ranges <- tibble(
      Start = quantile(data$totalReadsLog, probs = seq(0, 1, length.out = num_regions + 1)[-length(seq(0, 1, length.out = num_regions + 1))]),
      End = quantile(data$totalReadsLog, probs = seq(0, 1, length.out = num_regions + 1)[-1])
    )
    
    results <- data %>%
      group_by(region = cut(totalReadsLog, breaks = c(ranges$Start, max(ranges$End)), include.lowest = TRUE)) %>%
      group_modify(~ {
        fit_0 <- fit_distribution(.x$percentMouse, lower_bound = TRUE, dist = "lnorm")
        fit_1 <- fit_distribution(.x$percentMouse, lower_bound = FALSE, dist = "lnorm")
        
        .x %>%
          rowwise() %>%
          mutate(
            GMM = if (!is.null(fit_0) && !is.null(fit_1)) {
              classify_point(percentMouse, fit_0, fit_1, dist = "lnorm")
            } else {
              "mixed"
            }
          )
      }) %>%
      ungroup()
    
    summary_counts <- results %>%
      count(GMM) %>%
      tidyr::pivot_wider(names_from = GMM, values_from = n, values_fill = list(n = 0)) %>%
      mutate(n_regions = num_regions) %>%
      dplyr::select(n_regions, human, mouse, mixed)
    
    results_list <- bind_rows(results_list, summary_counts)
  }
  
  p1 <- ggplot(results, aes(x = totalReadsLog, y = percentMouse, color = x10X_Call)) +
    geom_point(size = 3, alpha = 0.7) +
    theme_minimal() +
    labs(title = "Scatterplot of percentMouse by x10X_Call",
         x = "totalReadsLog",
         y = "percentMouse")
  
  p2 <- ggplot(results, aes(x = totalReadsLog, y = percentMouse, color = GMM)) +
    geom_point(size = 3, alpha = 0.7) +
    theme_minimal() +
    labs(title = "Scatterplot of percentMouse by GMM",
         x = "totalReadsLog",
         y = "percentMouse")
  
  p3_hm <- ggplot(results_list, aes(x = n_regions)) +
    geom_line(aes(y = human, color = "human"), size = 1) +
    geom_line(aes(y = mouse, color = "mouse"), size = 1) +
    geom_hline(yintercept = 3301, linetype = "dashed", color = "blue") +
    geom_hline(yintercept = 3429, linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(title = "Classification of Human and Mouse by Number of Regions",
         x = "Number of Regions",
         y = "Count",
         color = "Classification")
  
  p3_mixed <- ggplot(results_list, aes(x = n_regions)) +
    geom_line(aes(y = mixed, color = "mixed"), size = 1) +
    geom_hline(yintercept = 459, linetype = "dashed", color = "green") +
    theme_minimal() +
    labs(title = "Classification of Mixed by Number of Regions",
         x = "Number of Regions",
         y = "Count",
         color = "Classification")
  
  print(p1)
  print(p2)
  print(p3_hm)
  print(p3_mixed)
  
  print(results_list)
  
  # Evaluating distribution best fit
  library(fitdistrplus)
  library(MASS)
  
  # Separate data into groups
  grch38_data <- data[data$x10X_Call == "GRCh38", ]$percentMouse
  mm10_data <- data[data$x10X_Call == "mm10", ]$percentMouse
  
  # Invert mm10 values (left-tailed, starts near 1)
  mm10_inverted <- 1 - mm10_data
  
  # Plot distributions
  par(mfrow = c(2, 2))
  hist(grch38_data, breaks = 30, main = "GRCh38 percentMouse", xlab = "percentMouse")
  hist(mm10_inverted, breaks = 30, main = "mm10 (Inverted) percentMouse", xlab = "1 - percentMouse")
  
  # Fit distributions to GRCh38
  grch38_fits <- list(
    gamma = fitdist(grch38_data, "gamma"),
    exp = fitdist(grch38_data, "exp"),
    norm = fitdist(grch38_data, "norm"),
    lnorm = fitdist(grch38_data, "lnorm"),
    weibull = fitdist(grch38_data, "weibull"),
    beta = fitdist(grch38_data, "beta"),
    cauchy = fitdist(grch38_data, "cauchy")
  )
  
  # Fit distributions to mm10 (inverted)
  mm10_fits <- list(
    gamma = fitdist(mm10_inverted, "gamma"),
    exp = fitdist(mm10_inverted, "exp"),
    norm = fitdist(mm10_inverted, "norm"),
    lnorm = fitdist(mm10_inverted, "lnorm"),
    weibull = fitdist(mm10_inverted, "weibull"),
    beta = fitdist(mm10_inverted, "beta"),
    cauchy = fitdist(mm10_inverted, "cauchy")
  )
  
  # Assess goodness of fit
  rank_fits <- function(fit_list) {
    aic_values <- sapply(fit_list, function(fit) fit$aic)
    aic_values[order(aic_values)]
  }
  
  cat("GRCh38 AIC Rankings:\n")
  print(rank_fits(grch38_fits))
  
  cat("mm10 (Inverted) AIC Rankings:\n")
  print(rank_fits(mm10_fits))
  
  # Visual inspection
  par(mfrow = c(2, 2))
  denscomp(grch38_fits, legendtext = names(grch38_fits), main = "GRCh38 Fits")
  denscomp(mm10_fits, legendtext = names(mm10_fits), main = "mm10 Inverted Fits")
  
  # Summarize the best fits
  cat("Best Fit for GRCh38:\n")
  summary(grch38_fits[[names(rank_fits(grch38_fits))[1]]])
  
  cat("Best Fit for mm10 (Inverted):\n")
  summary(mm10_fits[[names(rank_fits(mm10_fits))[1]]])
  
  
  
  
  
  
  
    
  p11 <- ggplot(gem_df, aes(x = totalReads, y = percentMouse, color = MESSY_Call)) +
    geom_point() +  # Add points
    scale_color_manual(values = c("Human" = alpha("royalblue", 0.5), "Mouse" = alpha("red3", 0.5), "Multiplet" = alpha("purple", 0.5))) +  # Specify colors
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
    scale_color_manual(values = c("Human" = alpha("royalblue", 0.5), "Mouse" = alpha("red3", 0.5), "Multiplet" = alpha("purple", 0.5))) +  # Specify colors
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
    xlab("Log10 Read Count") +
    geom_vline(xintercept = 25)
  
  
  p4 <- gem_df_melt[gem_df_melt$CellCall=="mm10",] %>% 
    ggplot(aes(color=Aligned2Ref, x=log(ReadCount), fill= Aligned2Ref, color= Aligned2Ref)) + 
    geom_density(alpha = 0.2) + 
    scale_color_manual(values = c("GRCh38_count" = "royalblue", "mm10_count" = "forestgreen")) +  # Specify colors
    scale_fill_manual(values = c("GRCh38_count" = "royalblue", "mm10_count" = "forestgreen")) +  # Specify colors
    #geom_histogram()
    theme_classic() +
    ylab("mm10 Cell Count Density") +
    xlab("Log10 Read Count") +
    geom_vline(xintercept = 25)
  
  p5 <- gem_df_melt[gem_df_melt$CellCall=="Multiplet",] %>% 
    ggplot(aes(color=Aligned2Ref, x=log(ReadCount), fill= Aligned2Ref, color= Aligned2Ref)) + 
    geom_density(alpha = 0.2) + 
    scale_color_manual(values = c("GRCh38_count" = "royalblue", "mm10_count" = "forestgreen")) +  # Specify colors
    scale_fill_manual(values = c("GRCh38_count" = "royalblue", "mm10_count" = "forestgreen")) +  # Specify colors
    #geom_histogram()
    theme_classic() +
    ylab("Multiplet Cell Count Density") +
    xlab("Log10 Read Count") +
    geom_vline(xintercept = 25)
  
  
  gem_df_melt2 <- as.data.frame(reshape2::melt(gem_df[,c("GRCh38_count","mm10_count","MESSY_Call")]))
  names(gem_df_melt2) <- c("OlexCellCall", "Aligned2Ref", "ReadCount")
  
  # Visualize the number UMIs/transcripts per cell
  p33 <- gem_df_melt2[gem_df_melt2$OlexCellCall=="Human",] %>% 
    ggplot(aes(color=Aligned2Ref, x=log(ReadCount), fill= Aligned2Ref, color= Aligned2Ref)) + 
    geom_density(alpha = 0.2) + 
    scale_color_manual(values = c("GRCh38_count" = "royalblue", "mm10_count" = "red3")) +  # Specify colors
    scale_fill_manual(values = c("GRCh38_count" = "royalblue", "mm10_count" = "red3")) +  # Specify colors
    #geom_histogram()
    theme_classic() +
    ylab("Human Cell Count Density") +
    xlab("Log10 Read Count") +
    geom_vline(xintercept = 25)
  
  
  p44 <- gem_df_melt2[gem_df_melt2$OlexCellCall=="Mouse",] %>% 
    ggplot(aes(color=Aligned2Ref, x=log(ReadCount), fill= Aligned2Ref, color= Aligned2Ref)) + 
    geom_density(alpha = 0.2) + 
    scale_color_manual(values = c("GRCh38_count" = "royalblue", "mm10_count" = "red3")) +  # Specify colors
    scale_fill_manual(values = c("GRCh38_count" = "royalblue", "mm10_count" = "red3")) +  # Specify colors
    #geom_histogram()
    theme_classic() +
    ylab("Mouse Cell Count Density") +
    xlab("Log10 Read Count") +
    geom_vline(xintercept = 25)
  
  p55 <- gem_df_melt2[gem_df_melt2$OlexCellCall=="Multiplet",] %>% 
    ggplot(aes(color=Aligned2Ref, x=log(ReadCount), fill= Aligned2Ref, color= Aligned2Ref)) + 
    geom_density(alpha = 0.2) + 
    scale_color_manual(values = c("GRCh38_count" = "royalblue", "mm10_count" = "red3")) +  # Specify colors
    scale_fill_manual(values = c("GRCh38_count" = "royalblue", "mm10_count" = "red3")) +  # Specify colors
    #geom_histogram()
    theme_classic() +
    ylab("Multiplet Cell Count Density") +
    xlab("Log10 Read Count") +
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

setwd("10X_gem_classifications/")
## read ing GEM file for 2 samples
pc62 <- read.delim("VCU-PC-062_gem_classification.csv", header = TRUE, row.names = 1, sep=",")
pc65 <- read.delim("VCU-PC-065_gem_classification.csv", header = TRUE, row.names = 1, sep=",")
pc67 <- read.delim("VCU-PC-067_gem_classification.csv", header = TRUE, row.names = 1, sep=",")
pc75 <- read.delim("VCU-PC-075_gem_classification.csv", header = TRUE, row.names = 1, sep=",")

pc72 <- read.delim("VCU-PC-072_gem_classification.csv", header = TRUE, row.names = 1, sep=",")
pc68 <- read.delim("VCU-PC-068_gem_classification.csv", header = TRUE, row.names = 1, sep=",")
pc66 <- read.delim("VCU-PC-066_gem_classification.csv", header = TRUE, row.names = 1, sep=",")
pc61 <- read.delim("VCU-PC-061_gem_classification.csv", header = TRUE, row.names = 1, sep=",")
co63 <- read.delim("VCU-CO-063_gem_classification.csv", header = TRUE, row.names = 1, sep=",")
bc74 <- read.delim("VCU-BC-074_gem_classification.csv", header = TRUE, row.names = 1, sep=",")
bc43 <- read.delim("VCU-BC-043_gem_classification.csv", header = TRUE, row.names = 1, sep=",")
bc37 <- read.delim("VCU-BC-037_gem_classification.csv", header = TRUE, row.names = 1, sep=",")


dir.create("Analyses/02_VCU-PC-062_1803_CRv8_241016/outs/analysis_olexClass/", recursive = TRUE)
gem_df = pc62; title = "PC-062"; save_dir = "Analyses/02_VCU-PC-062_1803_CRv8_241016/outs/analysis_olexClass/"
multiplet_plots(gem_df = pc62, title = "PC-062", save_dir = "Analyses/02_VCU-PC-062_1803_CRv8_241016/outs/analysis_olexClass/")
multiplet_plots(pc65, "PC-065", "Analyses/02_VCU-PC-065_1815_CRv8_241016/outs/analysis_olexClass/")
multiplet_plots(pc67, "PC-067", "Analyses/02_VCU-PC-067_1823_CRv8_241016/outs/analysis_olexClass/")
multiplet_plots(pc75, "PC-075", "Analyses/02_VCU-PC-075_1836_CRv8_241016/outs/analysis_olexClass/")

multiplet_plots(pc72, "PC-072", "./")
multiplet_plots(pc68, "PC-068", "./")
multiplet_plots(pc66, "PC-066", "./")
multiplet_plots(pc61, "PC-061", "./")
multiplet_plots(co63, "CO-063", "./")
multiplet_plots(bc74, "BC-074", "./")
multiplet_plots(bc43, "BC-043", "./")
multiplet_plots(bc37, "BC-037", "./")


## Dead cell file paths:
pc62_dead <- read.delim("/lustre/home/harrell_lab/scRNASeq/04_cellrangerv8count_HumanOnly_grch38/04_VCU-PC-062_1803_CRv8_241016/outs/analysis_deadcells/FindDeadCells_mito13human_CRv8_NextSeq241009_241016_VCU-PC-062_1803_deadcells.csv", header = TRUE)
pc65_dead <- read.delim("/lustre/home/harrell_lab/scRNASeq/04_cellrangerv8count_HumanOnly_grch38/04_VCU-PC-065_1815_CRv8_241016/outs/analysis_deadcells/FindDeadCells_mito13human_CRv8_NextSeq241009_241016_VCU-PC-065_1815_deadcells.csv", header = TRUE)
pc67_dead <- read.delim("/lustre/home/harrell_lab/scRNASeq/04_cellrangerv8count_HumanOnly_grch38/04_VCU-PC-067_1823_CRv8_241016/outs/analysis_deadcells/FindDeadCells_mito13human_CRv8_NextSeq241009_241016_VCU-PC-067_1823_deadcells.csv", header = TRUE)
pc75_dead <- read.delim("/lustre/home/harrell_lab/scRNASeq/04_cellrangerv8count_HumanOnly_grch38/04_VCU-PC-075_1836_CRv8_241016/outs/analysis_deadcells/FindDeadCells_mito13human_CRv8_NextSeq241009_241016_VCU-PC-075_1836_deadcells.csv", header = TRUE)


pc62_dead$CellState <- "dead"
pc65_dead$CellState <- "dead"
pc67_dead$CellState <- "dead"
pc75_dead$CellState <- "dead"

write.csv(pc62_dead, file="VCU-PC-062_1803_deadcell_annotation.csv", quote = FALSE, row.names = FALSE)
write.csv(pc65_dead, file="VCU-PC-065_1815_deadcell_annotation.csv", quote = FALSE, row.names = FALSE)
write.csv(pc67_dead, file="VCU-PC-067_1823_deadcell_annotation.csv", quote = FALSE, row.names = FALSE)
write.csv(pc75_dead, file="VCU-PC-075_1836_deadcell_annotation.csv", quote = FALSE, row.names = FALSE)







## Get all barcodes where human > mouse
pc75_human <- pc75[pc75$GRCh38 > pc75$mm10,]
pc75_human <- pc75_human[order(pc75_human$MouseDiff, decreasing = TRUE),]
pc75_mouse <- pc75[pc75$mm10 > pc75$GRCh38,]
pc75_mouse <- pc75_mouse[order(pc75_mouse$HumanDiff, decreasing = TRUE),]

pc75_mouse_10p <- quantile(pc75_mouse$mm10, c(.10)) ## 1188.8
pc75_human_10p <- quantile(pc75_human$GRCh38, c(.10)) ## 1446


pc65_human <- pc65[pc65$GRCh38 > pc65$mm10,]
pc65_human <- pc65_human[order(pc65_human$MouseDiff, decreasing = TRUE),]
pc65_mouse <- pc65[pc65$mm10 > pc65$GRCh38,]
pc65_mouse <- pc65_mouse[order(pc65_mouse$HumanDiff, decreasing = TRUE),]

pc65_mouse_10p <- quantile(pc65_mouse$mm10, c(.10)) ## 1279.8
pc65_human_10p <- quantile(pc65_human$GRCh38, c(.10)) ## 866


pc62$HumanDiff <- pc62$GRCh38 - pc62$mm10
pc62$MouseDiff <- pc62$mm10 - pc62$GRCh38

pc62_human <- pc62[pc62$GRCh38 > pc62$mm10,]
pc62_human <- pc62_human[order(pc62_human$MouseDiff, decreasing = TRUE),]
pc62_mouse <- pc62[pc62$mm10 > pc62$GRCh38,]
pc62_mouse <- pc62_mouse[order(pc62_mouse$HumanDiff, decreasing = TRUE),]

pc62_mouse_10p <- quantile(pc62_mouse$mm10, c(.10)) ## 1946
pc62_human_10p <- quantile(pc62_human$GRCh38, c(.10)) ## 1002









