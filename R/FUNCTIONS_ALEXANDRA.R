# =============================================================================
# INITIAL DATA VISUALIZATION FUNCTION
# This function visualizes raw data BEFORE any classification
# Uses ground truth 'call' column if available for coloring
# =============================================================================

plot_initial_data <- function(data,
                              xlim_max = NULL,
                              ylim = c(0, 1),
                              point_size = 0.5,
                              alpha = 0.6,
                              title = "Initial Data Overview (Before Classification)",
                              show_stats = TRUE) {

  # =============================================================================
  # DATA PREPARATION
  # =============================================================================

  # Handle different column naming conventions
  # Check for mm10/hg19 and convert to GRCh38/GRCm39
  if ("hg19" %in% colnames(data)) {
    colnames(data)[colnames(data) == "hg19"] <- "GRCh38"
  }
  if ("mm10" %in% colnames(data)) {
    colnames(data)[colnames(data) == "mm10"] <- "GRCm39"
  }

  # Also check for human/mouse naming
  if ("human" %in% tolower(colnames(data))) {
    idx <- which(tolower(colnames(data)) == "human")
    colnames(data)[idx] <- "GRCh38"
  }
  if ("mouse" %in% tolower(colnames(data))) {
    idx <- which(tolower(colnames(data)) == "mouse")
    colnames(data)[idx] <- "GRCm39"
  }

  # Validate required columns exist
  required_cols <- c("GRCh38", "GRCm39")
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", "),
               "\nPlease ensure your data has GRCh38 and GRCm39 columns (or hg19/mm10)"))
  }

  # Calculate total reads and percent mouse if not present
  if (!"totalReads" %in% colnames(data)) {
    data$totalReads <- data$GRCh38 + data$GRCm39
  }
  if (!"percentMouse" %in% colnames(data)) {
    data$percentMouse <- data$GRCm39 / data$totalReads
    # Handle division by zero
    data$percentMouse[!is.finite(data$percentMouse)] <- 0
  }

  # Set xlim if not provided
  if (is.null(xlim_max)) {
    xlim_max <- quantile(data$totalReads, 0.99, na.rm = TRUE)
  }

  # =============================================================================
  # DETERMINE COLORING SCHEME
  # =============================================================================

  # Check if ground truth 'call' column exists
  has_ground_truth <- "call" %in% colnames(data)

  if (has_ground_truth) {
    # Use ground truth for coloring
    cat("Ground truth 'call' column found. Using it for visualization.\n")

    # Standardize call values (in case of different naming)
    data$call_standard <- data$call
    data$call_standard[data$call %in% c("GRCh38", "hg19", "human", "Human")] <- "Human"
    data$call_standard[data$call %in% c("GRCm39", "mm10", "mouse", "Mouse")] <- "Mouse"
    data$call_standard[data$call %in% c("Multiplet", "multiplet", "doublet", "Doublet")] <- "Multiplet"

    # Define colors for each category
    color_map <- c("Human" = rgb(0, 0, 1, alpha),      # Blue
                   "Mouse" = rgb(1, 0, 0, alpha),      # Red
                   "Multiplet" = rgb(0.5, 0, 0.5, alpha)) # Purple

    # Get colors for each point
    point_colors <- color_map[data$call_standard]
    point_colors[is.na(point_colors)] <- rgb(0.5, 0.5, 0.5, alpha) # Gray for unknown

  } else {
    # No ground truth - use gradient coloring based on percentMouse
    cat("No ground truth 'call' column found. Using gradient coloring based on percent mouse.\n")

    # Create gradient from blue (human) to red (mouse)
    gradient <- colorRampPalette(c("blue", "purple", "red"))(100)
    color_indices <- cut(data$percentMouse, breaks = 100, labels = FALSE)
    point_colors <- gradient[color_indices]
    point_colors[is.na(point_colors)] <- "gray"
  }

  # =============================================================================
  # CREATE PLOT
  # =============================================================================

  # Save current par settings
  old_par <- par(no.readonly = TRUE)
  par(mar = c(5, 4, 4, 2))

  # Create base plot
  plot(data$totalReads, data$percentMouse,
       xlim = c(0, xlim_max),
       ylim = ylim,
       type = "n",
       main = title,
       xlab = "Total Reads (GRCh38 + GRCm39)",
       ylab = "Percent Mouse (GRCm39 / Total)",
       cex.main = 1.3,
       cex.lab = 1.2,
       cex.axis = 1.1)

  # Add grid for better readability
  grid(col = "lightgray", lty = "dotted")

  # Add reference lines
  abline(h = 0.5, col = "gray40", lwd = 1.5, lty = 2)  # 50% line
  abline(h = c(0.1, 0.9), col = "gray60", lwd = 1, lty = 3)  # 10% and 90% lines

  # Plot points
  points(data$totalReads, data$percentMouse,
         col = point_colors,
         pch = 16,
         cex = point_size)

  # =============================================================================
  # ADD ANNOTATIONS AND LEGEND
  # =============================================================================

  # Add legend based on coloring scheme
  if (has_ground_truth) {
    # Count cells in each category
    human_count <- sum(data$call_standard == "Human", na.rm = TRUE)
    mouse_count <- sum(data$call_standard == "Mouse", na.rm = TRUE)
    multiplet_count <- sum(data$call_standard == "Multiplet", na.rm = TRUE)
    unknown_count <- sum(is.na(data$call_standard))

    legend_labels <- c(
      paste0("Human (n=", human_count, ")"),
      paste0("Mouse (n=", mouse_count, ")"),
      paste0("Multiplet (n=", multiplet_count, ")")
    )
    legend_colors <- c(rgb(0, 0, 1, 0.8), rgb(1, 0, 0, 0.8), rgb(0.5, 0, 0.5, 0.8))

    if (unknown_count > 0) {
      legend_labels <- c(legend_labels, paste0("Unknown (n=", unknown_count, ")"))
      legend_colors <- c(legend_colors, rgb(0.5, 0.5, 0.5, 0.8))
    }

    legend("topright",
           legend = legend_labels,
           col = legend_colors,
           pch = 16,
           pt.cex = 1.2,
           cex = 0.9,
           title = "Ground Truth Classification",
           bg = "white",
           box.lwd = 1)

  } else {
    # Add gradient legend for percent mouse
    legend("topright",
           legend = c("0% Mouse (Human)", "50% Mouse", "100% Mouse"),
           col = c("blue", "purple", "red"),
           pch = 16,
           pt.cex = 1.2,
           cex = 0.9,
           title = "Percent Mouse Gradient",
           bg = "white",
           box.lwd = 1)
  }

  # Add statistics box if requested
  if (show_stats) {
    total_cells <- nrow(data)
    median_reads <- median(data$totalReads, na.rm = TRUE)
    mean_reads <- mean(data$totalReads, na.rm = TRUE)

    stats_text <- paste0(
      "Total cells: ", format(total_cells, big.mark = ","), "\n",
      "Median reads: ", format(round(median_reads), big.mark = ","), "\n",
      "Mean reads: ", format(round(mean_reads), big.mark = ","), "\n"
    )

    # Add ground truth stats if available
    if (has_ground_truth) {
      stats_text <- paste0(
        stats_text,
        "Human: ", round(human_count/total_cells*100, 1), "%\n",
        "Mouse: ", round(mouse_count/total_cells*100, 1), "%\n",
        "Multiplet: ", round(multiplet_count/total_cells*100, 1), "%"
      )
    }

    # Position stats box based on data distribution
    x_pos <- xlim_max * 0.02
    y_pos <- 0.3

    legend(x_pos, y_pos,
           legend = stats_text,
           bty = "n",
           cex = 0.8,
           title = "Data Summary")
  }

  # Add title annotations
  if (!has_ground_truth) {
    mtext("Note: No ground truth available - using gradient coloring",
          side = 3, line = 0.3, cex = 0.8, col = "gray40")
  }

  # Restore par settings
  par(old_par)

  # Return data summary invisibly
  invisible(list(
    total_cells = nrow(data),
    has_ground_truth = has_ground_truth,
    median_reads = median(data$totalReads, na.rm = TRUE),
    mean_reads = mean(data$totalReads, na.rm = TRUE),
    percent_below_1000_reads = sum(data$totalReads < 1000, na.rm = TRUE) / nrow(data) * 100,
    percent_above_50000_reads = sum(data$totalReads > 50000, na.rm = TRUE) / nrow(data) * 100
  ))
}


library(fitdistrplus)  # For fitdist() function

# =============================================================================
# HELPER FUNCTION: STANDARDIZE COLUMN NAMES
# =============================================================================# =============================================================================
# SHARED HELPER FUNCTIONS (used by multiple plotting functions)
# =============================================================================

# Helper function for column name standardization
standardize_column_names <- function(data) {
  # Handle different column naming conventions
  if ("hg19" %in% colnames(data)) {
    colnames(data)[colnames(data) == "hg19"] <- "GRCh38"
  }
  if ("mm10" %in% colnames(data)) {
    colnames(data)[colnames(data) == "mm10"] <- "GRCm39"
  }
  col_lower <- tolower(colnames(data))
  if ("human" %in% col_lower) {
    idx <- which(col_lower == "human")
    colnames(data)[idx] <- "GRCh38"
  }
  if ("mouse" %in% col_lower) {
    idx <- which(col_lower == "mouse")
    colnames(data)[idx] <- "GRCm39"
  }
  return(data)
}

# Helper function for fitted density (used by plot_low_reads, plot_high_reads, plot_multiplets)
get_fitted_density <- function(x_seq, fit, distribution) {
  if (distribution == "weibull") {
    return(dweibull(x_seq, shape = fit$estimate["shape"], scale = fit$estimate["scale"]))
  } else if (distribution == "lnorm") {
    return(dlnorm(x_seq, meanlog = fit$estimate["meanlog"], sdlog = fit$estimate["sdlog"]))
  } else if (distribution == "gamma") {
    return(dgamma(x_seq, shape = fit$estimate["shape"], rate = fit$estimate["rate"]))
  } else if (distribution == "norm") {
    return(dnorm(x_seq, mean = fit$estimate["mean"], sd = fit$estimate["sd"]))
  } else if (distribution == "exp") {
    return(dexp(x_seq, rate = fit$estimate["rate"]))
  } else {
    stop(paste("Unsupported distribution:", distribution))
  }
}





# =============================================================================
# FUNCTION 1: CLASSIFY LOW READS
# =============================================================================

classify_low_reads <- function(data, n_sd = 1.0, distribution = "auto", verbose = FALSE) {
  # Classify cells with abnormally low total reads
  # Parameters:
  # - data: dataframe with GRCh38/GRCm39 (or hg19/mm10, or human/mouse) columns
  # - n_sd: number of standard deviations below mean for threshold
  # - distribution: "auto", "weibull", "lnorm", "gamma", or "norm"
  # - verbose: whether to print summary information
  # Returns:
  # - Same dataframe with added 'low_classification' column

  # STANDARDIZE COLUMN NAMES FIRST
  data <- standardize_column_names(data)

  # Validate that required columns now exist
  required_cols <- c("GRCh38", "GRCm39")
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", "),
               "\nPlease ensure your data has human and mouse read count columns",
               "\nAccepted names: GRCh38/GRCm39, hg19/mm10, human/mouse"))
  }

  # Calculate total reads if not present
  if (!"totalReads" %in% colnames(data)) {
    data$totalReads <- data$GRCh38 + data$GRCm39
  }

  # Validate input data
  if (nrow(data) == 0) {
    stop("Input data frame is empty")
  }

  # Remove any invalid values for fitting
  valid_reads <- data$totalReads[is.finite(data$totalReads) & data$totalReads > 0]

  if (length(valid_reads) == 0) {
    stop("No valid read counts found")
  }

  # Fit distribution
  best_fit <- NULL
  best_dist <- distribution

  if (distribution == "auto") {
    # Test different distributions and select best by AIC
    dist_names <- c("weibull", "lnorm", "gamma", "norm")
    best_aic <- Inf

    for (dist in dist_names) {
      fit <- tryCatch({
        fitdist(valid_reads, dist)
      }, error = function(e) NULL)

      if (!is.null(fit) && fit$aic < best_aic) {
        best_aic <- fit$aic
        best_fit <- fit
        best_dist <- dist
      }
    }
  } else {
    # Use specified distribution
    best_fit <- tryCatch({
      fitdist(valid_reads, distribution)
    }, error = function(e) {
      stop(paste("Failed to fit", distribution, "distribution"))
    })
  }

  if (is.null(best_fit)) {
    stop("Failed to fit any distribution to the data")
  }

  # Calculate threshold based on distribution type
  threshold <- calculate_low_threshold(best_fit, best_dist, n_sd)

  # Ensure threshold is positive
  threshold <- max(threshold, 0)

  # Check if low_classification column already exists
  if ("low_classification" %in% colnames(data)) {
    if (verbose) {
      message("Warning: 'low_classification' column already exists. Overwriting...")
    }
  }

  # Classify cells
  data$low_classification <- ifelse(data$totalReads < threshold, "Low", "Normal")

  # Calculate results
  n_low <- sum(data$low_classification == "Low", na.rm = TRUE)
  n_normal <- sum(data$low_classification == "Normal", na.rm = TRUE)
  low_percentage <- round(n_low / nrow(data) * 100, 1)

  # Print summary if requested
  if (verbose) {
    cat("\nLow-read classification summary:\n")
    cat("• Distribution used:", best_dist, "\n")
    cat("• Sigma (σ):", n_sd, "\n")
    cat("• Threshold:", format(round(threshold, 0), big.mark = ","), "reads\n")
    cat("• Low cells:", n_low, "/", nrow(data), "(", low_percentage, "%)\n")
    cat("• Normal cells:", n_normal, "/", nrow(data), "(", round(100 - low_percentage, 1), "%)\n")
  }

  return(data)
}

# Helper function for threshold calculation
calculate_low_threshold <- function(fit, distribution, n_sd) {

  if (distribution == "weibull") {
    shape <- fit$estimate["shape"]
    scale <- fit$estimate["scale"]

    # Calculate Weibull mean and SD
    weibull_mean <- scale * gamma(1 + 1/shape)
    weibull_var <- scale^2 * (gamma(1 + 2/shape) - gamma(1 + 1/shape)^2)
    weibull_sd <- sqrt(weibull_var)

    threshold <- weibull_mean - n_sd * weibull_sd

  } else if (distribution == "lnorm") {
    meanlog <- fit$estimate["meanlog"]
    sdlog <- fit$estimate["sdlog"]

    # For log-normal, work in log space
    threshold <- exp(meanlog - n_sd * sdlog)

  } else if (distribution == "gamma") {
    shape <- fit$estimate["shape"]
    rate <- fit$estimate["rate"]

    gamma_mean <- shape / rate
    gamma_sd <- sqrt(shape) / rate

    threshold <- gamma_mean - n_sd * gamma_sd

  } else if (distribution == "norm") {
    norm_mean <- fit$estimate["mean"]
    norm_sd <- fit$estimate["sd"]

    threshold <- norm_mean - n_sd * norm_sd

  } else {
    stop(paste("Unsupported distribution:", distribution))
  }

  return(threshold)
}

# =============================================================================
# FUNCTION 2: CLASSIFY HIGH READS
# =============================================================================

classify_high_reads <- function(data, n_sd = 1.5, distribution = "auto", verbose = FALSE) {
  # Classify cells with abnormally high total reads (potential doublets)
  # Parameters:
  # - data: dataframe with GRCh38/GRCm39 (or hg19/mm10, or human/mouse) columns
  # - n_sd: number of standard deviations above mean for threshold
  # - distribution: "auto", "weibull", "lnorm", "gamma", or "norm"
  # - verbose: whether to print summary information
  # Returns:
  # - Same dataframe with added 'high_classification' column

  # STANDARDIZE COLUMN NAMES FIRST
  data <- standardize_column_names(data)

  # Validate that required columns now exist
  required_cols <- c("GRCh38", "GRCm39")
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", "),
               "\nPlease ensure your data has human and mouse read count columns",
               "\nAccepted names: GRCh38/GRCm39, hg19/mm10, human/mouse"))
  }

  # Calculate total reads if not present
  if (!"totalReads" %in% colnames(data)) {
    data$totalReads <- data$GRCh38 + data$GRCm39
  }

  # Validate input data
  if (nrow(data) == 0) {
    stop("Input data frame is empty")
  }

  # Remove any invalid values for fitting
  valid_reads <- data$totalReads[is.finite(data$totalReads) & data$totalReads > 0]

  if (length(valid_reads) == 0) {
    stop("No valid read counts found")
  }

  # Fit distribution
  best_fit <- NULL
  best_dist <- distribution

  if (distribution == "auto") {
    # Test different distributions and select best by AIC
    dist_names <- c("weibull", "lnorm", "gamma", "norm")
    best_aic <- Inf

    for (dist in dist_names) {
      fit <- tryCatch({
        fitdist(valid_reads, dist)
      }, error = function(e) NULL)

      if (!is.null(fit) && fit$aic < best_aic) {
        best_aic <- fit$aic
        best_fit <- fit
        best_dist <- dist
      }
    }
  } else {
    # Use specified distribution
    best_fit <- tryCatch({
      fitdist(valid_reads, distribution)
    }, error = function(e) {
      stop(paste("Failed to fit", distribution, "distribution"))
    })
  }

  if (is.null(best_fit)) {
    stop("Failed to fit any distribution to the data")
  }

  # Calculate threshold based on distribution type
  threshold <- calculate_high_threshold(best_fit, best_dist, n_sd)

  # Check if high_classification column already exists
  if ("high_classification" %in% colnames(data)) {
    if (verbose) {
      message("Warning: 'high_classification' column already exists. Overwriting...")
    }
  }

  # Classify cells (key difference: > instead of <)
  data$high_classification <- ifelse(data$totalReads > threshold, "High", "Normal")

  # Calculate results
  n_high <- sum(data$high_classification == "High", na.rm = TRUE)
  n_normal <- sum(data$high_classification == "Normal", na.rm = TRUE)
  high_percentage <- round(n_high / nrow(data) * 100, 1)

  # Print summary if requested
  if (verbose) {
    cat("\nHigh-read classification summary:\n")
    cat("• Distribution used:", best_dist, "\n")
    cat("• Sigma (σ):", n_sd, "\n")
    cat("• Threshold:", format(round(threshold, 0), big.mark = ","), "reads\n")
    cat("• High cells:", n_high, "/", nrow(data), "(", high_percentage, "%)\n")
    cat("• Normal cells:", n_normal, "/", nrow(data), "(", round(100 - high_percentage, 1), "%)\n")
  }

  return(data)
}

# Helper function for high threshold calculation
calculate_high_threshold <- function(fit, distribution, n_sd) {

  if (distribution == "weibull") {
    shape <- fit$estimate["shape"]
    scale <- fit$estimate["scale"]

    # Calculate Weibull mean and SD
    weibull_mean <- scale * gamma(1 + 1/shape)
    weibull_var <- scale^2 * (gamma(1 + 2/shape) - gamma(1 + 1/shape)^2)
    weibull_sd <- sqrt(weibull_var)

    # For HIGH reads, we ADD standard deviations
    threshold <- weibull_mean + n_sd * weibull_sd

  } else if (distribution == "lnorm") {
    meanlog <- fit$estimate["meanlog"]
    sdlog <- fit$estimate["sdlog"]

    # For log-normal, work in log space - ADD for high threshold
    threshold <- exp(meanlog + n_sd * sdlog)

  } else if (distribution == "gamma") {
    shape <- fit$estimate["shape"]
    rate <- fit$estimate["rate"]

    gamma_mean <- shape / rate
    gamma_sd <- sqrt(shape) / rate

    # For HIGH reads, we ADD standard deviations
    threshold <- gamma_mean + n_sd * gamma_sd

  } else if (distribution == "norm") {
    norm_mean <- fit$estimate["mean"]
    norm_sd <- fit$estimate["sd"]

    # For HIGH reads, we ADD standard deviations
    threshold <- norm_mean + n_sd * norm_sd

  } else {
    stop(paste("Unsupported distribution:", distribution))
  }

  return(threshold)
}

# =============================================================================
# FUNCTION 3: CLASSIFY MULTIPLETS
# =============================================================================

classify_multiplets <- function(data, n_sd = 3.0, distribution = "auto", verbose = FALSE) {
  # Classify cells as Human, Mouse, or Multiplet based on percent mouse
  # Uses flipping approach for symmetric distribution
  # Parameters:
  # - data: dataframe with GRCh38/GRCm39 (or hg19/mm10, or human/mouse) columns
  # - n_sd: number of standard deviations for outlier detection
  # - distribution: "auto", "lnorm", "gamma", "norm", or "exp"
  # - verbose: whether to print summary information
  # Returns:
  # - Same dataframe with added 'multiplet_classification' column

  # STANDARDIZE COLUMN NAMES FIRST
  data <- standardize_column_names(data)

  # Validate that required columns now exist
  required_cols <- c("GRCh38", "GRCm39")
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", "),
               "\nPlease ensure your data has human and mouse read count columns",
               "\nAccepted names: GRCh38/GRCm39, hg19/mm10, human/mouse"))
  }

  # Calculate percent mouse if not present
  if (!"percentMouse" %in% colnames(data)) {
    data$percentMouse <- data$GRCm39 / (data$GRCh38 + data$GRCm39)
  }

  # Validate input data
  if (nrow(data) == 0) {
    stop("Input data frame is empty")
  }

  # Remove any NaN or infinite values
  valid_indices <- is.finite(data$percentMouse)
  n_invalid <- sum(!valid_indices)

  if (sum(valid_indices) == 0) {
    stop("No valid percentMouse values found")
  }

  # Work with valid data for fitting
  valid_data <- data[valid_indices, ]

  # =============================================================================
  # APPLY FLIPPING LOGIC
  # =============================================================================

  # Separate human and mouse cells
  human_cells <- valid_data$percentMouse[valid_data$percentMouse <= 0.5]
  mouse_cells_original <- valid_data$percentMouse[valid_data$percentMouse > 0.5]

  # Flip mouse cells to create symmetric distribution
  mouse_cells_flipped <- 1 - mouse_cells_original

  # Combine flipped data
  combined_data <- c(human_cells, mouse_cells_flipped)

  # Remove exact 0s and 1s for distribution fitting
  combined_data_clean <- combined_data[combined_data > 0 & combined_data < 1]

  if (length(combined_data_clean) == 0) {
    stop("No valid data points for distribution fitting after flipping")
  }

  # =============================================================================
  # FIT DISTRIBUTION TO COMBINED FLIPPED DATA
  # =============================================================================

  best_fit <- NULL
  best_dist <- distribution

  if (distribution == "auto") {
    # Test different distributions
    dist_names <- c("lnorm", "gamma", "norm", "exp")
    best_aic <- Inf

    for (dist in dist_names) {
      fit <- tryCatch({
        fitdist(combined_data_clean, dist)
      }, error = function(e) NULL)

      if (!is.null(fit) && fit$aic < best_aic) {
        best_aic <- fit$aic
        best_fit <- fit
        best_dist <- dist
      }
    }
  } else {
    # Use specified distribution
    best_fit <- tryCatch({
      fitdist(combined_data_clean, distribution)
    }, error = function(e) {
      stop(paste("Failed to fit", distribution, "distribution"))
    })
  }

  if (is.null(best_fit)) {
    stop("Failed to fit any distribution to the flipped data")
  }

  # =============================================================================
  # CALCULATE THRESHOLD BASED ON DISTRIBUTION
  # =============================================================================

  threshold_flipped <- calculate_multiplet_threshold(best_fit, best_dist, n_sd)

  # Ensure threshold is within bounds
  threshold_flipped <- min(max(threshold_flipped, 0), 1)

  # Convert back to original space
  human_threshold <- threshold_flipped
  mouse_threshold <- 1 - threshold_flipped

  # Check if multiplet_classification column already exists
  if ("multiplet_classification" %in% colnames(data)) {
    if (verbose) {
      message("Warning: 'multiplet_classification' column already exists. Overwriting...")
    }
  }

  # =============================================================================
  # CLASSIFY CELLS
  # =============================================================================

  # Initialize classification column
  data$multiplet_classification <- NA

  # Classify valid cells
  data$multiplet_classification[valid_indices] <- ifelse(
    data$percentMouse[valid_indices] < human_threshold, "Human",
    ifelse(data$percentMouse[valid_indices] > mouse_threshold, "Mouse", "Multiplet")
  )

  # Mark invalid cells
  if (n_invalid > 0) {
    data$multiplet_classification[!valid_indices] <- "Invalid"
  }

  # Calculate results
  class_counts <- table(data$multiplet_classification, useNA = "ifany")

  # Print summary if requested
  if (verbose) {
    cat("\nMultiplet classification summary:\n")
    cat("• Distribution used:", best_dist, "\n")
    cat("• Sigma (σ):", n_sd, "\n")
    cat("• Human threshold: <", round(human_threshold * 100, 1), "% mouse\n")
    cat("• Mouse threshold: >", round(mouse_threshold * 100, 1), "% mouse\n")
    cat("• Multiplet zone:", round(human_threshold * 100, 1), "% -", round(mouse_threshold * 100, 1), "%\n")

    for (class in names(class_counts)) {
      cat("• ", class, " cells: ", class_counts[class], " / ", nrow(data),
          " (", round(class_counts[class]/nrow(data)*100, 1), "%)\n", sep="")
    }
  }

  return(data)
}

# Helper function for multiplet threshold calculation
calculate_multiplet_threshold <- function(fit, distribution, n_sd) {

  if (distribution == "lnorm") {
    meanlog <- fit$estimate["meanlog"]
    sdlog <- fit$estimate["sdlog"]

    # Calculate threshold in flipped space
    threshold_flipped <- exp(meanlog + n_sd * sdlog)

  } else if (distribution == "gamma") {
    shape <- fit$estimate["shape"]
    rate <- fit$estimate["rate"]

    gamma_mean <- shape / rate
    gamma_sd <- sqrt(shape) / rate

    threshold_flipped <- gamma_mean + n_sd * gamma_sd

  } else if (distribution == "exp") {
    rate <- fit$estimate["rate"]

    exp_mean <- 1 / rate
    exp_sd <- 1 / rate

    threshold_flipped <- exp_mean + n_sd * exp_sd

  } else if (distribution == "norm") {
    norm_mean <- fit$estimate["mean"]
    norm_sd <- fit$estimate["sd"]

    threshold_flipped <- norm_mean + n_sd * norm_sd

  } else {
    stop(paste("Unsupported distribution:", distribution))
  }

  return(threshold_flipped)
}

# =============================================================================
# FUNCTION: PLOT MULTI-SIGMA COMPARISON FOR LOW READS
# =============================================================================

plot_low_reads_multi_sigma <- function(data, sigma_levels = c(0.5, 0.75, 1.0, 1.25, 1.5, 2.0),
                                       distribution = "auto", verbose = FALSE) {
  # Create multi-sigma comparison plots for low reads classification
  # Parameters:
  # - data: dataframe with GRCh38 and GRCm39 columns
  # - sigma_levels: vector of sigma values to test
  # - distribution: distribution type to use
  # - verbose: whether to print progress
  # Returns:
  # - results_summary: dataframe with results for each sigma level

  # STANDARDIZE COLUMN NAMES FIRST
  data <- standardize_column_names(data)
  # Validate input
  if (!"totalReads" %in% colnames(data)) {
    data$totalReads <- data$GRCh38 + data$GRCm39
  }

  # Create results summary table
  results_summary <- data.frame(
    Sigma_Level = numeric(),
    Threshold_Reads = numeric(),
    Low_Count = numeric(),
    Low_Percentage = numeric(),
    Normal_Count = numeric(),
    Normal_Percentage = numeric(),
    Human_Low = numeric(),
    Mouse_Low = numeric(),
    Multiplet_Low = numeric()
  )

  # Get valid data and fit distribution once
  valid_reads <- data$totalReads[is.finite(data$totalReads) & data$totalReads > 0]

  if (distribution == "auto") {
    best_fit <- fitdist(valid_reads, "weibull")  # Default to Weibull for consistency
    best_dist <- "weibull"
  } else {
    best_fit <- fitdist(valid_reads, distribution)
    best_dist <- distribution
  }

  # Calculate distribution parameters (assuming Weibull for multi-sigma)
  shape_param <- best_fit$estimate["shape"]
  scale_param <- best_fit$estimate["scale"]
  weibull_mean <- scale_param * gamma(1 + 1/shape_param)
  weibull_var <- scale_param^2 * (gamma(1 + 2/shape_param) - gamma(1 + 1/shape_param)^2)
  weibull_sd <- sqrt(weibull_var)

  # Set up plotting layout
  old_par <- par(no.readonly = TRUE)
  par(mfrow = c(2, 3),
      mar = c(4, 4, 3, 2),
      cex.main = 1.2,
      cex.lab = 1.1,
      cex.axis = 1.0)

  # Test each sigma level
  for (i in 1:length(sigma_levels)) {
    sigma <- sigma_levels[i]

    if (verbose) {
      cat("--- Testing σ =", sigma, "---\n")
    }

    # Apply classification
    data_temp <- classify_low_reads(data, n_sd = sigma, distribution = distribution, verbose = FALSE)

    # Calculate threshold
    threshold <- max(weibull_mean - sigma * weibull_sd, 0)

    # Count classifications
    low_count <- sum(data_temp$low_classification == "Low")
    normal_count <- sum(data_temp$low_classification == "Normal")
    low_percentage <- round(low_count / nrow(data) * 100, 1)
    normal_percentage <- round(normal_count / nrow(data) * 100, 1)

    # Create enhanced histogram
    hist(data$totalReads, breaks = 50, col = "lightblue",
         main = paste0("Step 1: Low-Read Classification (σ = ", sigma, ")"),
         xlab = "Total Reads", ylab = "Frequency",
         xlim = c(0, 50000))

    # Store histogram data for highlighting
    hist_data <- hist(data$totalReads, breaks = 50, plot = FALSE)

    # Create colored bars: pink for outliers, lightblue for normal
    bar_colors <- ifelse(hist_data$mids <= (threshold + 1),
                         rgb(1, 0.7, 0.7, 0.8),  # Light pink for outliers
                         "lightblue")             # Light blue for normal

    # Redraw with colored bars
    hist(data$totalReads, breaks = 50, col = bar_colors, border = "white",
         main = paste0("Step 1: Low-Read Classification (σ = ", sigma, ")"),
         xlab = "Total Reads", ylab = "Frequency",
         xlim = c(0, 50000), add = TRUE)

    # Add lines and annotations
    abline(v = weibull_mean, col = "darkgreen", lwd = 2.5, lty = 1)
    abline(v = threshold, col = "red", lwd = 2.5, lty = 2)

    # Add distance arrow
    if (threshold > 0 && threshold < weibull_mean) {
      arrow_y <- par("usr")[4] * 0.75
      arrows(x0 = threshold, y0 = arrow_y, x1 = weibull_mean, y1 = arrow_y,
             col = "purple", lwd = 2, code = 3, length = 0.05, angle = 90)
      text(x = (threshold + weibull_mean) / 2, y = arrow_y * 1.05,
           labels = paste0(sigma, "σ"), col = "purple", cex = 1, font = 2)
    }

    # Add annotations
    text(weibull_mean + 1000, par("usr")[4] * 0.85,
         paste0("Mean\n", format(round(weibull_mean), big.mark = ",")),
         col = "darkgreen", cex = 0.8, font = 2)

    if (threshold > 0) {
      text(threshold - 1000, par("usr")[4] * 0.5,
           paste0("Threshold\n", format(round(threshold), big.mark = ",")),
           col = "red", cex = 0.8, font = 2)
    }

    # Add subtitle
    mtext(paste0("Threshold: ", format(round(threshold), big.mark = ","),
                 " reads | ", low_count, " cells classified as 'low' (", low_percentage, "%)"),
          side = 3, line = 0.3, cex = 0.9)

    # Add legends
    legend("topright",
           legend = c("Mean", "Threshold", paste0(sigma, "σ distance"), "Outlier region"),
           col = c("darkgreen", "red", "purple", rgb(1, 0.6, 0.6, 0.6)),
           lty = c(1, 2, 1, NA), pch = c(NA, NA, NA, 15),
           lwd = c(2.5, 2.5, 2, NA), cex = 1.0, bg = "white", box.lty = 1,
           title = "Legend")

    legend(x = 35000, y = par("usr")[4] * 0.6,
           legend = c(paste("Low cells:", low_count, "(", low_percentage, "%)"),
                      paste("Normal cells:", normal_count, "(", normal_percentage, "%)"),
                      "Distribution: Weibull"),
           cex = 1.0, bg = "white", box.lty = 1, title = "Statistics")

    # Analyze ground truth if available
    human_low <- mouse_low <- multiplet_low <- 0
    if ("call" %in% colnames(data)) {
      low_cells <- data_temp[data_temp$low_classification == "Low", ]
      if (nrow(low_cells) > 0) {
        ground_truth_table <- table(low_cells$call)
        human_low <- ifelse("GRCh38" %in% names(ground_truth_table), ground_truth_table["GRCh38"], 0)
        mouse_low <- ifelse("GRCm39" %in% names(ground_truth_table), ground_truth_table["GRCm39"], 0)
        multiplet_low <- ifelse("Multiplet" %in% names(ground_truth_table), ground_truth_table["Multiplet"], 0)
      }
    }

    # Add to results summary
    results_summary <- rbind(results_summary, data.frame(
      Sigma_Level = sigma,
      Threshold_Reads = round(threshold),
      Low_Count = low_count,
      Low_Percentage = low_percentage,
      Normal_Count = normal_count,
      Normal_Percentage = normal_percentage,
      Human_Low = human_low,
      Mouse_Low = mouse_low,
      Multiplet_Low = multiplet_low
    ))

    if (verbose) {
      cat("Threshold:", format(round(threshold), big.mark = ","), "reads\n")
      cat("Low cells:", low_count, "/", nrow(data), "(", low_percentage, "%)\n")
      cat("Human low:", human_low, "| Mouse low:", mouse_low, "| Multiplet low:", multiplet_low, "\n\n")
    }
  }

  # Restore par settings
  par(old_par)

  return(results_summary)
}

# =============================================================================
# FUNCTION: PLOT MULTI-SIGMA COMPARISON FOR HIGH READS
# =============================================================================

plot_high_reads_multi_sigma <- function(data, sigma_levels = c(0.5, 0.75, 1.0, 1.25, 1.5, 2.0),
                                        distribution = "auto", verbose = FALSE) {
  # Create multi-sigma comparison plots for high reads classification
  # Similar structure to low reads version but for high thresholds

  # STANDARDIZE COLUMN NAMES FIRST
  data <- standardize_column_names(data)
  # Validate input
  if (!"totalReads" %in% colnames(data)) {
    data$totalReads <- data$GRCh38 + data$GRCm39
  }

  # Create results summary table
  results_summary <- data.frame(
    Sigma_Level = numeric(),
    Threshold_Reads = numeric(),
    High_Count = numeric(),
    High_Percentage = numeric(),
    Normal_Count = numeric(),
    Normal_Percentage = numeric(),
    Human_High = numeric(),
    Mouse_High = numeric(),
    Multiplet_High = numeric()
  )

  # Get valid data and fit distribution once
  valid_reads <- data$totalReads[is.finite(data$totalReads) & data$totalReads > 0]

  if (distribution == "auto") {
    best_fit <- fitdist(valid_reads, "weibull")
    best_dist <- "weibull"
  } else {
    best_fit <- fitdist(valid_reads, distribution)
    best_dist <- distribution
  }

  # Calculate Weibull parameters
  shape_param <- best_fit$estimate["shape"]
  scale_param <- best_fit$estimate["scale"]
  weibull_mean <- scale_param * gamma(1 + 1/shape_param)
  weibull_var <- scale_param^2 * (gamma(1 + 2/shape_param) - gamma(1 + 1/shape_param)^2)
  weibull_sd <- sqrt(weibull_var)

  # Set up plotting layout
  old_par <- par(no.readonly = TRUE)
  par(mfrow = c(2, 3),
      mar = c(4, 4, 3, 2),
      cex.main = 1.2,
      cex.lab = 1.1,
      cex.axis = 1.0)

  # Test each sigma level
  for (i in 1:length(sigma_levels)) {
    sigma <- sigma_levels[i]

    if (verbose) {
      cat("--- Testing σ =", sigma, "---\n")
    }

    # Apply classification
    data_temp <- classify_high_reads(data, n_sd = sigma, distribution = distribution, verbose = FALSE)

    # Calculate threshold (PLUS for high reads)
    threshold <- weibull_mean + sigma * weibull_sd

    # Count classifications
    high_count <- sum(data_temp$high_classification == "High")
    normal_count <- sum(data_temp$high_classification == "Normal")
    high_percentage <- round(high_count / nrow(data) * 100, 1)
    normal_percentage <- round(normal_count / nrow(data) * 100, 1)

    # Create enhanced histogram
    hist(data$totalReads, breaks = 50, col = "lightblue",
         main = paste0("Step 2: High-Read Classification (σ = ", sigma, ")"),
         xlab = "Total Reads", ylab = "Frequency",
         xlim = c(0, 50000))

    # Store histogram data for highlighting
    hist_data <- hist(data$totalReads, breaks = 50, plot = FALSE)

    # Create colored bars: orange for high read outliers, lightblue for normal
    bar_colors <- ifelse(hist_data$mids >= threshold,
                         rgb(1, 0.6, 0.2, 0.8),  # Orange for high read outliers
                         "lightblue")             # Light blue for normal

    # Redraw with colored bars
    hist(data$totalReads, breaks = 50, col = bar_colors, border = "white",
         main = paste0("Step 2: High-Read Classification (σ = ", sigma, ")"),
         xlab = "Total Reads", ylab = "Frequency",
         xlim = c(0, 50000), add = TRUE)

    # Add lines and annotations
    abline(v = weibull_mean, col = "darkgreen", lwd = 2.5, lty = 1)
    abline(v = threshold, col = "red", lwd = 2.5, lty = 2)

    # Add distance arrow (to the RIGHT for high reads)
    if (threshold > weibull_mean && threshold < 50000) {
      arrow_y <- par("usr")[4] * 0.75
      arrows(x0 = weibull_mean, y0 = arrow_y, x1 = threshold, y1 = arrow_y,
             col = "purple", lwd = 2, code = 3, length = 0.05, angle = 90)
      text(x = (threshold + weibull_mean) / 2, y = arrow_y * 1.05,
           labels = paste0(sigma, "σ"), col = "purple", cex = 1, font = 2)
    }

    # Add annotations
    text(weibull_mean, par("usr")[4] * 0.9,
         paste0("Mean\n", format(round(weibull_mean), big.mark = ",")),
         pos = 2, col = "darkgreen", cex = 0.8, font = 2)

    if (threshold < 50000) {
      text(threshold, par("usr")[4] * 0.5,
           paste0("Threshold\n", format(round(threshold), big.mark = ",")),
           pos = 4, col = "red", cex = 0.8, font = 2)
    }

    # Add subtitle
    mtext(paste0("Threshold: ", format(round(threshold), big.mark = ","),
                 " reads | ", high_count, " cells classified as 'high' (", high_percentage, "%)"),
          side = 3, line = 0.3, cex = 0.9)

    # Add legends
    legend("topright",
           legend = c("Mean", "Threshold", paste0(sigma, "σ distance"), "Outlier region"),
           col = c("darkgreen", "red", "purple", rgb(1, 0.6, 0.2, 0.8)),
           lty = c(1, 2, 1, NA), pch = c(NA, NA, NA, 15),
           lwd = c(2.5, 2.5, 2, NA), cex = 1.0, bg = "white", box.lty = 1,
           title = "Legend")

    legend(x = 25000, y = par("usr")[4] * 0.6,
           legend = c(paste("High cells:", high_count, "(", high_percentage, "%)"),
                      paste("Normal cells:", normal_count, "(", normal_percentage, "%)"),
                      "Distribution: Weibull"),
           cex = 1.0, bg = "white", box.lty = 1, title = "Statistics")

    # Analyze ground truth if available
    human_high <- mouse_high <- multiplet_high <- 0
    if ("call" %in% colnames(data)) {
      high_cells <- data_temp[data_temp$high_classification == "High", ]
      if (nrow(high_cells) > 0) {
        ground_truth_table <- table(high_cells$call)
        human_high <- ifelse("GRCh38" %in% names(ground_truth_table), ground_truth_table["GRCh38"], 0)
        mouse_high <- ifelse("GRCm39" %in% names(ground_truth_table), ground_truth_table["GRCm39"], 0)
        multiplet_high <- ifelse("Multiplet" %in% names(ground_truth_table), ground_truth_table["Multiplet"], 0)
      }
    }

    # Add to results summary
    results_summary <- rbind(results_summary, data.frame(
      Sigma_Level = sigma,
      Threshold_Reads = round(threshold),
      High_Count = high_count,
      High_Percentage = high_percentage,
      Normal_Count = normal_count,
      Normal_Percentage = normal_percentage,
      Human_High = human_high,
      Mouse_High = mouse_high,
      Multiplet_High = multiplet_high
    ))

    if (verbose) {
      cat("Threshold:", format(round(threshold), big.mark = ","), "reads\n")
      cat("High cells:", high_count, "/", nrow(data), "(", high_percentage, "%)\n")
      cat("Human high:", human_high, "| Mouse high:", mouse_high, "| Multiplet high:", multiplet_high, "\n\n")
    }
  }

  # Restore par settings
  par(old_par)

  return(results_summary)
}

# =============================================================================
# FUNCTION: PLOT MULTI-SIGMA COMPARISON FOR MULTIPLETS
# =============================================================================

plot_multiplets_multi_sigma <- function(data, sigma_levels = c(1.0, 1.5, 2.0, 2.5, 3.0, 4.0),
                                        distribution = "auto", verbose = FALSE) {
  # Create multi-sigma comparison plots for multiplet classification

  # STANDARDIZE COLUMN NAMES FIRST
  data <- standardize_column_names(data)
  # Validate input
  if (!"percentMouse" %in% colnames(data)) {
    data$percentMouse <- data$GRCm39 / (data$GRCh38 + data$GRCm39)
  }

  # Create results summary table
  results_summary <- data.frame(
    Sigma_Level = numeric(),
    Threshold_Flipped = numeric(),
    Human_Threshold = numeric(),
    Mouse_Threshold = numeric(),
    Multiplet_Zone_Width = numeric(),
    Human_Count = numeric(),
    Mouse_Count = numeric(),
    Multiplet_Count = numeric(),
    Multiplet_Percentage = numeric(),
    True_Multiplets_Captured = numeric(),
    Sensitivity = numeric(),
    Precision = numeric()
  )

  # Get valid data and apply flipping logic once
  valid_indices <- is.finite(data$percentMouse)
  valid_data <- data[valid_indices, ]
  valid_pm <- data$percentMouse[valid_indices]

  human_pm <- valid_pm[valid_pm <= 0.5]
  mouse_pm_original <- valid_pm[valid_pm > 0.5]
  mouse_pm_flipped <- 1 - mouse_pm_original
  combined_flipped <- c(human_pm, mouse_pm_flipped)
  combined_clean <- combined_flipped[combined_flipped > 0 & combined_flipped < 1]

  # Fit distribution once
  best_fit <- fitdist(combined_clean, "lnorm")
  meanlog <- best_fit$estimate["meanlog"]
  sdlog <- best_fit$estimate["sdlog"]

  if (verbose) {
    cat("Using fitted log-normal distribution:\n")
    cat("• meanlog =", round(meanlog, 4), "\n")
    cat("• sdlog =", round(sdlog, 4), "\n\n")
  }

  # Get total true multiplets if ground truth available
  total_true_multiplets <- 0
  if ("call" %in% colnames(data)) {
    total_true_multiplets <- sum(data$call == "Multiplet")
  }

  # Set up plotting layout
  old_par <- par(no.readonly = TRUE)
  par(mfrow = c(2, 3),
      mar = c(4, 4, 3, 2),
      cex.main = 1.2,
      cex.lab = 1.1,
      cex.axis = 1.0)

  # Test each sigma level
  for (i in 1:length(sigma_levels)) {
    sigma <- sigma_levels[i]

    if (verbose) {
      cat("--- Testing σ =", sigma, "---\n")
    }

    # Apply classification
    data_temp <- classify_multiplets(data, n_sd = sigma, distribution = distribution, verbose = FALSE)

    # Calculate threshold using pre-fitted distribution
    threshold_flipped <- exp(meanlog + sigma * sdlog)
    threshold_flipped <- min(max(threshold_flipped, 0.001), 0.999)

    # Convert back to original space
    human_threshold <- threshold_flipped
    mouse_threshold <- 1 - threshold_flipped

    # Create histogram
    hist(valid_pm, breaks = 100, col = "lightgreen",
         main = paste0("Step 3: Multiplet Detection (σ = ", sigma, ")"),
         xlab = "Percent Mouse", ylab = "Frequency", xlim = c(0, 1))

    # Add reference lines
    abline(v = 0.5, col = "black", lwd = 2, lty = 3)  # 50% line
    abline(v = human_threshold, col = "blue", lwd = 2.5, lty = 2)  # Human threshold
    abline(v = mouse_threshold, col = "red", lwd = 2.5, lty = 2)  # Mouse threshold

    # Shade multiplet zone
    rect(human_threshold, 0, mouse_threshold, par("usr")[4],
         col = rgb(0.5, 0.5, 0.5, 0.2), border = NA)

    # Add annotations
    text(0.15, par("usr")[4] * 0.9, "Human", col = "blue", cex = 1.2, font = 2)
    text(0.85, par("usr")[4] * 0.9, "Mouse", col = "red", cex = 1.2, font = 2)
    if (mouse_threshold - human_threshold > 0.1) {
      text(0.5, par("usr")[4] * 0.7, "Multiplet\nZone", col = "purple", cex = 1, font = 2)
    }

    # Count classifications
    human_count <- sum(data_temp$multiplet_classification == "Human", na.rm = TRUE)
    mouse_count <- sum(data_temp$multiplet_classification == "Mouse", na.rm = TRUE)
    multiplet_count <- sum(data_temp$multiplet_classification == "Multiplet", na.rm = TRUE)
    multiplet_percentage <- round(multiplet_count / nrow(data) * 100, 1)

    # Calculate performance metrics
    true_multiplets_captured <- 0
    sensitivity <- 0
    precision <- 0

    if ("call" %in% colnames(data)) {
      true_multiplets_captured <- sum(data$call == "Multiplet" &
                                        data_temp$multiplet_classification == "Multiplet", na.rm = TRUE)
      sensitivity <- if (total_true_multiplets > 0) {
        round(true_multiplets_captured / total_true_multiplets * 100, 1)
      } else { 0 }

      precision <- if (multiplet_count > 0) {
        round(true_multiplets_captured / multiplet_count * 100, 1)
      } else { 0 }
    }

    # Add subtitle with percentages
    mtext(paste0("H: ", human_count, " (", round(human_count/nrow(data)*100, 1), "%) | ",
                 "M: ", mouse_count, " (", round(mouse_count/nrow(data)*100, 1), "%) | ",
                 "Mult: ", multiplet_count, " (", multiplet_percentage, "%)"),
          side = 3, line = 0.3, cex = 0.9)

    # Add to results summary
    results_summary <- rbind(results_summary, data.frame(
      Sigma_Level = sigma,
      Threshold_Flipped = round(threshold_flipped, 4),
      Human_Threshold = round(human_threshold, 4),
      Mouse_Threshold = round(mouse_threshold, 4),
      Multiplet_Zone_Width = round(mouse_threshold - human_threshold, 4),
      Human_Count = human_count,
      Mouse_Count = mouse_count,
      Multiplet_Count = multiplet_count,
      Multiplet_Percentage = multiplet_percentage,
      True_Multiplets_Captured = true_multiplets_captured,
      Sensitivity = sensitivity,
      Precision = precision
    ))

    if (verbose) {
      cat("Thresholds: Human <", round(human_threshold, 4), "| Mouse >", round(mouse_threshold, 4), "\n")
      cat("Multiplet zone width:", round(mouse_threshold - human_threshold, 4), "\n")
      cat("Multiplets detected:", multiplet_count, "/", nrow(data), "(", multiplet_percentage, "%)\n")
      cat("True multiplets captured:", true_multiplets_captured, "/", total_true_multiplets, "(", sensitivity, "%)\n")
      cat("Precision:", precision, "%\n\n")
    }
  }

  # Restore par settings
  par(old_par)

  return(results_summary)
}

# =============================================================================
# FUNCTION: PLOT FLIPPED SPACE VISUALIZATION
# =============================================================================

plot_flipped_space <- function(data, sigma_levels = c(1.0, 1.5, 2.0, 2.5, 3.0, 4.0),
                               selected_sigma = 3.0) {
  # Create visualization showing flipped space analysis for different sigma levels
  # ADD THIS LINE:
  data <- standardize_column_names(data)
  # Validate input
  if (!"percentMouse" %in% colnames(data)) {
    data$percentMouse <- data$GRCm39 / (data$GRCh38 + data$GRCm39)
  }

  # Get valid data and apply flipping logic
  valid_indices <- is.finite(data$percentMouse)
  valid_pm <- data$percentMouse[valid_indices]

  human_pm <- valid_pm[valid_pm <= 0.5]
  mouse_pm_original <- valid_pm[valid_pm > 0.5]
  mouse_pm_flipped <- 1 - mouse_pm_original
  combined_flipped <- c(human_pm, mouse_pm_flipped)
  combined_clean <- combined_flipped[combined_flipped > 0 & combined_flipped < 1]

  # Fit distribution
  best_fit <- fitdist(combined_clean, "lnorm")
  meanlog <- best_fit$estimate["meanlog"]
  sdlog <- best_fit$estimate["sdlog"]
  geometric_mean <- exp(meanlog)

  # Save current par and set layout
  old_par <- par(no.readonly = TRUE)
  par(mfrow = c(2, 3),
      mar = c(4, 4, 3, 2),
      cex.main = 1.2,
      cex.lab = 1.1,
      cex.axis = 1.0)

  # Test different sigma values
  for (i in 1:length(sigma_levels)) {
    sigma <- sigma_levels[i]

    # Create histogram of flipped data
    hist(combined_clean, breaks = 100, col = "lightcyan", probability = TRUE,
         main = paste0("Flipped Space (σ = ", sigma, ")"),
         xlab = "Flipped Values", ylab = "Density",
         xlim = c(0, 1))

    # Add fitted log-normal curve
    x_seq <- seq(0.001, 1, length.out = 1000)
    y_fitted <- dlnorm(x_seq, meanlog = meanlog, sdlog = sdlog)
    lines(x_seq, y_fitted, col = "red", lwd = 2)

    # Calculate and show threshold
    threshold_flipped <- exp(meanlog + sigma * sdlog)
    threshold_flipped <- min(threshold_flipped, 0.999)

    # Add threshold line
    abline(v = threshold_flipped, col = "purple", lwd = 3, lty = 2)

    # Add geometric mean line
    abline(v = geometric_mean, col = "blue", lwd = 2, lty = 3)

    # Shade the outlier region
    x_shade <- seq(threshold_flipped, 1, length.out = 100)
    y_shade <- dlnorm(x_shade, meanlog = meanlog, sdlog = sdlog)
    polygon(c(threshold_flipped, x_shade, 1),
            c(0, y_shade, 0),
            col = rgb(1, 0, 0, 0.2), border = NA)

    # Highlight selected sigma
    if (sigma == selected_sigma) {
      box(col = "darkgreen", lwd = 3)
      mtext("SELECTED", side = 3, line = -1.5, col = "darkgreen", font = 2, cex = 0.8)
    }

    # Get current plot limits for text positioning
    usr <- par("usr")
    y_max <- usr[4]

    # Add annotations with dynamic positioning
    text(geometric_mean, y_max * 0.85, "μ", col = "blue", cex = 1.5, font = 2)

    if (threshold_flipped < 0.8) {
      text(threshold_flipped, y_max * 0.75, paste0("μ+", sigma, "σ"),
           col = "purple", cex = 1.2, font = 2, pos = 4)
    }

    # Add threshold value
    text(0.5, y_max * 0.6, paste("Threshold =", round(threshold_flipped, 3)),
         col = "purple", cex = 1.2)

    # Calculate percentage of outliers
    outlier_pct <- sum(combined_clean > threshold_flipped) / length(combined_clean) * 100
    text(0.5, y_max * 0.5, paste0("Outliers: ", round(outlier_pct, 1), "%"),
         col = "red", cex = 1.2, font = 2)
  }

  # Restore par settings
  par(old_par)

  invisible(list(meanlog = meanlog, sdlog = sdlog, geometric_mean = geometric_mean))
}

# =============================================================================
# FUNCTION: PLOT LOW READS CLASSIFICATION (UPDATED)
# =============================================================================

plot_low_reads <- function(data,
                           n_sd = 1.0,
                           distribution = "auto",
                           xlim_max = NULL,
                           plot_type = "standard",
                           run_classification = FALSE,
                           verbose = FALSE) {

  # Create visualizations for low reads classification
  # Parameters:
  # - data: dataframe with GRCh38/GRCm39 (or alternative naming) columns
  # - n_sd: sigma value to use for classification
  # - distribution: distribution type to use ("auto" or specific)
  # - xlim_max: maximum x-axis limit (auto if NULL)
  # - plot_type: "standard" for 2x2 layout, "single" for histogram only
  # - run_classification: if TRUE, run classification even if already done
  # - verbose: whether to print messages

  # Standardize column names first (for flexibility)
  data <- standardize_column_names(data)

  # Validate required columns exist
  required_cols <- c("GRCh38", "GRCm39")
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", "),
               "\nPlease ensure your data has human and mouse read count columns"))
  }

  # Calculate total reads if not present
  if (!"totalReads" %in% colnames(data)) {
    data$totalReads <- data$GRCh38 + data$GRCm39
    if (verbose) {
      cat("Calculated totalReads column\n")
    }
  }

  # Check if classification needs to be run
  needs_classification <- !"low_classification" %in% colnames(data)

  if (needs_classification && !run_classification) {
    # Ask user what to do
    message("Note: Data has not been classified yet.")
    message("Running classify_low_reads() with n_sd = ", n_sd)
    run_classification <- TRUE
  }

  # Run classification if needed or requested
  if (run_classification || needs_classification) {
    if (verbose) {
      cat("Running low reads classification with σ =", n_sd, "\n")
    }
    data <- classify_low_reads(data,
                               n_sd = n_sd,
                               distribution = distribution,
                               verbose = FALSE)
  } else if ("low_classification" %in% colnames(data)) {
    # Classification exists - check if parameters match
    if (verbose) {
      cat("Using existing low_classification column\n")
      cat("Note: If you want to re-run with different parameters, set run_classification = TRUE\n")
    }
  }

  # Get valid data for fitting
  valid_reads <- data$totalReads[is.finite(data$totalReads) & data$totalReads > 0]

  if (length(valid_reads) == 0) {
    stop("No valid read counts found for plotting")
  }

  # Fit distribution (same logic as classification function)
  if (distribution == "auto") {
    dist_names <- c("weibull", "lnorm", "gamma", "norm")
    best_aic <- Inf
    best_fit <- NULL
    best_dist <- NULL

    for (dist in dist_names) {
      fit <- tryCatch({
        fitdist(valid_reads, dist)
      }, error = function(e) NULL)

      if (!is.null(fit) && fit$aic < best_aic) {
        best_aic <- fit$aic
        best_fit <- fit
        best_dist <- dist
      }
    }
  } else {
    best_fit <- tryCatch({
      fitdist(valid_reads, distribution)
    }, error = function(e) {
      stop(paste("Failed to fit", distribution, "distribution"))
    })
    best_dist <- distribution
  }

  if (is.null(best_fit)) {
    stop("Failed to fit any distribution to the data")
  }

  # Calculate threshold
  threshold <- calculate_low_threshold_viz(best_fit, best_dist, n_sd)
  threshold <- max(threshold, 0)

  # Set xlim if not provided
  if (is.null(xlim_max)) {
    xlim_max <- quantile(data$totalReads, 0.99, na.rm = TRUE)
  }

  # Calculate statistics for display
  n_low <- sum(data$low_classification == "Low", na.rm = TRUE)
  n_normal <- sum(data$low_classification == "Normal", na.rm = TRUE)
  low_percentage <- round(n_low / nrow(data) * 100, 1)
  normal_percentage <- round(100 - low_percentage, 1)

  # Save current par settings
  old_par <- par(no.readonly = TRUE)

  if (plot_type == "standard") {
    # Standard 2x2 layout
    par(mfrow = c(2, 2),
        mar = c(4, 4, 3, 2),
        cex.main = 1.3,
        cex.lab = 1.1,
        cex.axis = 1.0)

    # Plot 1: Histogram with threshold and colored bars
    hist_data <- hist(data$totalReads, breaks = 50, plot = FALSE)

    # Create colored bars: pink for low reads, lightblue for normal
    bar_colors <- ifelse(hist_data$mids <= threshold,
                         rgb(1, 0.7, 0.7, 0.8),  # Light pink for low reads
                         "lightblue")             # Light blue for normal

    hist(data$totalReads, breaks = 50, col = bar_colors, border = "white",
         main = paste0("Step 1: Low-Read Classification (σ = ", n_sd, ")"),
         xlab = "Total Reads", ylab = "Frequency",
         xlim = c(0, xlim_max))

    abline(v = threshold, col = "red", lwd = 2.5, lty = 2)

    # Add text annotations
    text(threshold, par("usr")[4] * 0.9,
         paste("Threshold\n", format(round(threshold, 0), big.mark=",")),
         col = "red", pos = 4, cex = 1.1, font = 2)

    # Add statistics to plot
    legend("topright",
           legend = c(paste("Low:", n_low, "(", low_percentage, "%)"),
                      paste("Normal:", n_normal, "(", normal_percentage, "%)"),
                      paste("Distribution:", best_dist),
                      paste("Sigma:", n_sd)),
           bty = "n", cex = 0.9)

    # Plot 2: Density with fitted distribution
    plot(density(valid_reads), main = "Fitted Distribution vs Data",
         xlab = "Total Reads", ylab = "Density",
         xlim = c(0, quantile(valid_reads, 0.99)))

    # Add fitted distribution curve
    x_seq <- seq(min(valid_reads), quantile(valid_reads, 0.99), length.out = 1000)
    y_fitted <- get_fitted_density(x_seq, best_fit, best_dist)
    lines(x_seq, y_fitted, col = "red", lwd = 2)
    abline(v = threshold, col = "blue", lwd = 2, lty = 2)

    legend("topright",
           legend = c("Empirical density", paste("Fitted", best_dist), "Threshold"),
           col = c("black", "red", "blue"),
           lty = c(1, 1, 2),
           lwd = c(1, 2, 2),
           cex = 0.9)

    # Plot 3: Classification results bar chart
    class_table <- table(data$low_classification)
    barplot(class_table,
            col = c(rgb(1, 0.7, 0.7, 0.8), "lightblue"),
            main = "Classification Results",
            ylab = "Number of Cells",
            ylim = c(0, max(class_table) * 1.1))

    # Add count labels on bars
    text(x = seq(0.7, by = 1.2, length.out = length(class_table)),
         y = class_table + max(class_table) * 0.02,
         labels = paste0(class_table, "\n(",
                         round(class_table/sum(class_table)*100, 1), "%)"),
         cex = 0.9)

    # Plot 4: Box plot comparison (log scale)
    boxplot(totalReads ~ low_classification, data = data,
            col = c(rgb(1, 0.7, 0.7, 0.8), "lightblue"),
            main = "Total Reads by Classification",
            ylab = "Total Reads (log scale)",
            xlab = "Classification",
            log = "y",
            outline = FALSE)  # Don't show outliers as points

    # Add sample size to x-axis labels
    axis(1, at = 1:2,
         labels = paste0(c("Low", "Normal"), "\n(n=", class_table, ")"),
         tick = FALSE, line = 1)

  } else if (plot_type == "single") {
    # Single histogram plot
    par(mar = c(5, 4, 4, 2),
        cex.main = 1.4,
        cex.lab = 1.2,
        cex.axis = 1.1)

    # Create histogram with colored bars
    hist_data <- hist(data$totalReads, breaks = 50, plot = FALSE)

    bar_colors <- ifelse(hist_data$mids <= threshold,
                         rgb(1, 0.7, 0.7, 0.8),  # Light pink for low reads
                         "lightblue")             # Light blue for normal

    hist(data$totalReads, breaks = 50, col = bar_colors, border = "white",
         main = paste0("Step 1: Low-Read Classification (σ = ", n_sd, ")"),
         xlab = "Total Reads", ylab = "Frequency",
         xlim = c(0, xlim_max))

    abline(v = threshold, col = "red", lwd = 3, lty = 2)

    # Add threshold annotation
    text(threshold, par("usr")[4] * 0.9,
         paste("Threshold\n", format(round(threshold, 0), big.mark=",")),
         col = "red", pos = 4, cex = 1.2, font = 2)

    # Add subtitle with statistics
    mtext(paste0("Low cells: ", n_low, " (", low_percentage, "%) | ",
                 "Normal cells: ", n_normal, " (", normal_percentage, "%) | ",
                 "Distribution: ", best_dist),
          side = 3, line = 0.3, cex = 0.9)

    # Add legend
    legend("topright",
           legend = c("Low reads", "Normal reads", "Threshold"),
           fill = c(rgb(1, 0.7, 0.7, 0.8), "lightblue", NA),
           border = c("black", "black", NA),
           lty = c(NA, NA, 2),
           lwd = c(NA, NA, 3),
           col = c(NA, NA, "red"),
           cex = 1.0)
  }

  # Restore par settings
  par(old_par)

  # Print summary if verbose
  if (verbose) {
    cat("\n=== Low Reads Classification Plot Summary ===\n")
    cat("Distribution used:", best_dist, "\n")
    cat("Sigma (σ):", n_sd, "\n")
    cat("Threshold:", format(round(threshold, 0), big.mark = ","), "reads\n")
    cat("Low cells:", n_low, "/", nrow(data), "(", low_percentage, "%)\n")
    cat("Normal cells:", n_normal, "/", nrow(data), "(", normal_percentage, "%)\n")
  }

  # Return useful information invisibly
  invisible(list(
    data = data,  # Return the potentially classified data
    threshold = threshold,
    distribution = best_dist,
    n_sd = n_sd,
    low_count = n_low,
    normal_count = n_normal,
    low_percentage = low_percentage
  ))
}

# =============================================================================
# HELPER FUNCTION FOR THRESHOLD CALCULATION (needed by plot_low_reads)
# =============================================================================

calculate_low_threshold_viz <- function(fit, distribution, n_sd) {
  if (distribution == "weibull") {
    shape <- fit$estimate["shape"]
    scale <- fit$estimate["scale"]
    weibull_mean <- scale * gamma(1 + 1/shape)
    weibull_var <- scale^2 * (gamma(1 + 2/shape) - gamma(1 + 1/shape)^2)
    weibull_sd <- sqrt(weibull_var)
    threshold <- weibull_mean - n_sd * weibull_sd
  } else if (distribution == "lnorm") {
    meanlog <- fit$estimate["meanlog"]
    sdlog <- fit$estimate["sdlog"]
    threshold <- exp(meanlog - n_sd * sdlog)
  } else if (distribution == "gamma") {
    shape <- fit$estimate["shape"]
    rate <- fit$estimate["rate"]
    gamma_mean <- shape / rate
    gamma_sd <- sqrt(shape) / rate
    threshold <- gamma_mean - n_sd * gamma_sd
  } else if (distribution == "norm") {
    norm_mean <- fit$estimate["mean"]
    norm_sd <- fit$estimate["sd"]
    threshold <- norm_mean - n_sd * norm_sd
  } else {
    stop(paste("Unsupported distribution:", distribution))
  }
  return(threshold)
}


# =============================================================================
# FUNCTION: PLOT HIGH READS CLASSIFICATION (UPDATED)
# =============================================================================

plot_high_reads <- function(data,
                            n_sd = 1.5,
                            distribution = "auto",
                            xlim_max = NULL,
                            plot_type = "standard",
                            run_classification = FALSE,
                            verbose = FALSE) {

  # Create visualizations for high reads classification
  # Parameters:
  # - data: dataframe with GRCh38/GRCm39 (or alternative naming) columns
  # - n_sd: sigma value to use for classification
  # - distribution: distribution type to use ("auto" or specific)
  # - xlim_max: maximum x-axis limit (auto if NULL)
  # - plot_type: "standard" for 2x2 layout, "single" for histogram only
  # - run_classification: if TRUE, run classification even if already done
  # - verbose: whether to print messages

  # Standardize column names first (for flexibility)
  data <- standardize_column_names(data)

  # Validate required columns exist
  required_cols <- c("GRCh38", "GRCm39")
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", "),
               "\nPlease ensure your data has human and mouse read count columns"))
  }

  # Calculate total reads if not present
  if (!"totalReads" %in% colnames(data)) {
    data$totalReads <- data$GRCh38 + data$GRCm39
    if (verbose) {
      cat("Calculated totalReads column\n")
    }
  }

  # Check if classification needs to be run
  needs_classification <- !"high_classification" %in% colnames(data)

  if (needs_classification && !run_classification) {
    # Ask user what to do
    message("Note: Data has not been classified for high reads yet.")
    message("Running classify_high_reads() with n_sd = ", n_sd)
    run_classification <- TRUE
  }

  # Run classification if needed or requested
  if (run_classification || needs_classification) {
    if (verbose) {
      cat("Running high reads classification with σ =", n_sd, "\n")
    }
    data <- classify_high_reads(data,
                                n_sd = n_sd,
                                distribution = distribution,
                                verbose = FALSE)
  } else if ("high_classification" %in% colnames(data)) {
    # Classification exists - check if parameters match
    if (verbose) {
      cat("Using existing high_classification column\n")
      cat("Note: If you want to re-run with different parameters, set run_classification = TRUE\n")
    }
  }

  # Get valid data for fitting
  valid_reads <- data$totalReads[is.finite(data$totalReads) & data$totalReads > 0]

  if (length(valid_reads) == 0) {
    stop("No valid read counts found for plotting")
  }

  # Fit distribution (same logic as classification function)
  if (distribution == "auto") {
    dist_names <- c("weibull", "lnorm", "gamma", "norm")
    best_aic <- Inf
    best_fit <- NULL
    best_dist <- NULL

    for (dist in dist_names) {
      fit <- tryCatch({
        fitdist(valid_reads, dist)
      }, error = function(e) NULL)

      if (!is.null(fit) && fit$aic < best_aic) {
        best_aic <- fit$aic
        best_fit <- fit
        best_dist <- dist
      }
    }
  } else {
    best_fit <- tryCatch({
      fitdist(valid_reads, distribution)
    }, error = function(e) {
      stop(paste("Failed to fit", distribution, "distribution"))
    })
    best_dist <- distribution
  }

  if (is.null(best_fit)) {
    stop("Failed to fit any distribution to the data")
  }

  # Calculate threshold (HIGH version - ADD standard deviations)
  threshold <- calculate_high_threshold_viz(best_fit, best_dist, n_sd)

  # Set xlim if not provided
  if (is.null(xlim_max)) {
    xlim_max <- quantile(data$totalReads, 0.99, na.rm = TRUE)
  }

  # Calculate statistics for display
  n_high <- sum(data$high_classification == "High", na.rm = TRUE)
  n_normal <- sum(data$high_classification == "Normal", na.rm = TRUE)
  high_percentage <- round(n_high / nrow(data) * 100, 1)
  normal_percentage <- round(100 - high_percentage, 1)

  # Save current par settings
  old_par <- par(no.readonly = TRUE)

  if (plot_type == "standard") {
    # Standard 2x2 layout
    par(mfrow = c(2, 2),
        mar = c(4, 4, 3, 2),
        cex.main = 1.3,
        cex.lab = 1.1,
        cex.axis = 1.0)

    # Plot 1: Histogram with threshold and colored bars
    hist_data <- hist(data$totalReads, breaks = 50, plot = FALSE)

    # Create colored bars: orange for high reads, lightblue for normal
    bar_colors <- ifelse(hist_data$mids >= threshold,
                         rgb(1, 0.6, 0.2, 0.8),  # Orange for high reads
                         "lightblue")             # Light blue for normal

    hist(data$totalReads, breaks = 50, col = bar_colors, border = "white",
         main = paste0("Step 2: High-Read Classification (σ = ", n_sd, ")"),
         xlab = "Total Reads", ylab = "Frequency",
         xlim = c(0, xlim_max))

    abline(v = threshold, col = "red", lwd = 2.5, lty = 2)

    # Add text annotations
    text(threshold, par("usr")[4] * 0.9,
         paste("Threshold\n", format(round(threshold, 0), big.mark=",")),
         col = "red", pos = 2, cex = 1.1, font = 2)

    # Add statistics to plot
    legend("topright",
           legend = c(paste("High:", n_high, "(", high_percentage, "%)"),
                      paste("Normal:", n_normal, "(", normal_percentage, "%)"),
                      paste("Distribution:", best_dist),
                      paste("Sigma:", n_sd)),
           bty = "n", cex = 0.9)

    # Plot 2: Density with fitted distribution
    plot(density(valid_reads), main = "Fitted Distribution vs Data",
         xlab = "Total Reads", ylab = "Density",
         xlim = c(0, quantile(valid_reads, 0.99)))

    # Add fitted distribution curve
    x_seq <- seq(min(valid_reads), quantile(valid_reads, 0.99), length.out = 1000)
    y_fitted <- get_fitted_density(x_seq, best_fit, best_dist)
    lines(x_seq, y_fitted, col = "red", lwd = 2)
    abline(v = threshold, col = "blue", lwd = 2, lty = 2)

    legend("topright",
           legend = c("Empirical density", paste("Fitted", best_dist), "Threshold"),
           col = c("black", "red", "blue"),
           lty = c(1, 1, 2),
           lwd = c(1, 2, 2),
           cex = 0.9)

    # Plot 3: Classification results bar chart
    class_table <- table(data$high_classification)

    # Reorder to ensure "High" comes first if it exists
    if ("High" %in% names(class_table) && "Normal" %in% names(class_table)) {
      class_table <- class_table[c("High", "Normal")]
    }

    barplot(class_table,
            col = c(rgb(1, 0.6, 0.2, 0.8), "lightblue"),
            main = "Classification Results",
            ylab = "Number of Cells",
            ylim = c(0, max(class_table) * 1.1))

    # Add count labels on bars
    text(x = seq(0.7, by = 1.2, length.out = length(class_table)),
         y = class_table + max(class_table) * 0.02,
         labels = paste0(class_table, "\n(",
                         round(class_table/sum(class_table)*100, 1), "%)"),
         cex = 0.9)

    # Plot 4: Box plot comparison (log scale)
    # Reorder factor levels for consistent display
    data$high_classification <- factor(data$high_classification,
                                       levels = c("High", "Normal"))

    boxplot(totalReads ~ high_classification, data = data,
            col = c(rgb(1, 0.6, 0.2, 0.8), "lightblue"),
            main = "Total Reads by Classification",
            ylab = "Total Reads (log scale)",
            xlab = "Classification",
            log = "y",
            outline = FALSE)  # Don't show outliers as points

    # Add sample size to x-axis labels
    axis(1, at = 1:2,
         labels = paste0(c("High", "Normal"), "\n(n=",
                         c(n_high, n_normal), ")"),
         tick = FALSE, line = 1)

  } else if (plot_type == "single") {
    # Single histogram plot
    par(mar = c(5, 4, 4, 2),
        cex.main = 1.4,
        cex.lab = 1.2,
        cex.axis = 1.1)

    # Create histogram with colored bars
    hist_data <- hist(data$totalReads, breaks = 50, plot = FALSE)

    bar_colors <- ifelse(hist_data$mids >= threshold,
                         rgb(1, 0.6, 0.2, 0.8),  # Orange for high reads
                         "lightblue")             # Light blue for normal

    hist(data$totalReads, breaks = 50, col = bar_colors, border = "white",
         main = paste0("Step 2: High-Read Classification (σ = ", n_sd, ")"),
         xlab = "Total Reads", ylab = "Frequency",
         xlim = c(0, xlim_max))

    abline(v = threshold, col = "red", lwd = 3, lty = 2)

    # Add threshold annotation
    text(threshold, par("usr")[4] * 0.9,
         paste("Threshold\n", format(round(threshold, 0), big.mark=",")),
         col = "red", pos = 2, cex = 1.2, font = 2)

    # Add subtitle with statistics
    mtext(paste0("High cells: ", n_high, " (", high_percentage, "%) | ",
                 "Normal cells: ", n_normal, " (", normal_percentage, "%) | ",
                 "Distribution: ", best_dist),
          side = 3, line = 0.3, cex = 0.9)

    # Add legend
    legend("topleft",  # Changed position to avoid overlap with high threshold
           legend = c("High reads", "Normal reads", "Threshold"),
           fill = c(rgb(1, 0.6, 0.2, 0.8), "lightblue", NA),
           border = c("black", "black", NA),
           lty = c(NA, NA, 2),
           lwd = c(NA, NA, 3),
           col = c(NA, NA, "red"),
           cex = 1.0)
  }

  # Restore par settings
  par(old_par)

  # Print summary if verbose
  if (verbose) {
    cat("\n=== High Reads Classification Plot Summary ===\n")
    cat("Distribution used:", best_dist, "\n")
    cat("Sigma (σ):", n_sd, "\n")
    cat("Threshold:", format(round(threshold, 0), big.mark = ","), "reads\n")
    cat("High cells:", n_high, "/", nrow(data), "(", high_percentage, "%)\n")
    cat("Normal cells:", n_normal, "/", nrow(data), "(", normal_percentage, "%)\n")
  }

  # Return useful information invisibly
  invisible(list(
    data = data,  # Return the potentially classified data
    threshold = threshold,
    distribution = best_dist,
    n_sd = n_sd,
    high_count = n_high,
    normal_count = n_normal,
    high_percentage = high_percentage
  ))
}

# =============================================================================
# HELPER FUNCTION FOR HIGH THRESHOLD CALCULATION (needed by plot_high_reads)
# =============================================================================

calculate_high_threshold_viz <- function(fit, distribution, n_sd) {
  if (distribution == "weibull") {
    shape <- fit$estimate["shape"]
    scale <- fit$estimate["scale"]
    weibull_mean <- scale * gamma(1 + 1/shape)
    weibull_var <- scale^2 * (gamma(1 + 2/shape) - gamma(1 + 1/shape)^2)
    weibull_sd <- sqrt(weibull_var)
    # For HIGH reads, we ADD standard deviations
    threshold <- weibull_mean + n_sd * weibull_sd
  } else if (distribution == "lnorm") {
    meanlog <- fit$estimate["meanlog"]
    sdlog <- fit$estimate["sdlog"]
    # For log-normal, work in log space - ADD for high threshold
    threshold <- exp(meanlog + n_sd * sdlog)
  } else if (distribution == "gamma") {
    shape <- fit$estimate["shape"]
    rate <- fit$estimate["rate"]
    gamma_mean <- shape / rate
    gamma_sd <- sqrt(shape) / rate
    # For HIGH reads, we ADD standard deviations
    threshold <- gamma_mean + n_sd * gamma_sd
  } else if (distribution == "norm") {
    norm_mean <- fit$estimate["mean"]
    norm_sd <- fit$estimate["sd"]
    # For HIGH reads, we ADD standard deviations
    threshold <- norm_mean + n_sd * norm_sd
  } else {
    stop(paste("Unsupported distribution:", distribution))
  }
  return(threshold)
}

# Note: get_fitted_density() helper function should be defined elsewhere
# (same as in plot_low_reads) as it's shared between functions

# =============================================================================
# FUNCTION: PLOT MULTIPLET CLASSIFICATION (UPDATED)
# =============================================================================

plot_multiplets <- function(data,
                            n_sd = 3.0,
                            distribution = "auto",
                            plot_type = "standard",
                            run_classification = FALSE,
                            verbose = FALSE) {

  # Create visualizations for multiplet classification
  # Parameters:
  # - data: dataframe with GRCh38/GRCm39 (or alternative naming) columns
  # - n_sd: sigma value to use for classification
  # - distribution: distribution type to use ("auto", "lnorm", "gamma", "norm", "exp")
  # - plot_type: "standard" for 2x2 layout, "single" for histogram only
  # - run_classification: if TRUE, run classification even if already done
  # - verbose: whether to print messages

  # Standardize column names first (for flexibility)
  data <- standardize_column_names(data)

  # Validate required columns exist
  required_cols <- c("GRCh38", "GRCm39")
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", "),
               "\nPlease ensure your data has human and mouse read count columns"))
  }

  # Calculate percent mouse if not present
  if (!"percentMouse" %in% colnames(data)) {
    data$percentMouse <- data$GRCm39 / (data$GRCh38 + data$GRCm39)
    # Handle division by zero
    data$percentMouse[!is.finite(data$percentMouse)] <- 0
    if (verbose) {
      cat("Calculated percentMouse column\n")
    }
  }

  # Check if classification needs to be run
  needs_classification <- !"multiplet_classification" %in% colnames(data)

  if (needs_classification && !run_classification) {
    # Ask user what to do
    message("Note: Data has not been classified for multiplets yet.")
    message("Running classify_multiplets() with n_sd = ", n_sd)
    run_classification <- TRUE
  }

  # Run classification if needed or requested
  if (run_classification || needs_classification) {
    if (verbose) {
      cat("Running multiplet classification with σ =", n_sd, "\n")
    }
    data <- classify_multiplets(data,
                                n_sd = n_sd,
                                distribution = distribution,
                                verbose = FALSE)
  } else if ("multiplet_classification" %in% colnames(data)) {
    # Classification exists
    if (verbose) {
      cat("Using existing multiplet_classification column\n")
      cat("Note: If you want to re-run with different parameters, set run_classification = TRUE\n")
    }
  }

  # Get valid data for flipping approach
  valid_indices <- is.finite(data$percentMouse)
  valid_data <- data[valid_indices, ]

  if (nrow(valid_data) == 0) {
    stop("No valid percentMouse values found for plotting")
  }

  # Apply flipping logic
  human_cells <- valid_data$percentMouse[valid_data$percentMouse <= 0.5]
  mouse_cells_original <- valid_data$percentMouse[valid_data$percentMouse > 0.5]
  mouse_cells_flipped <- 1 - mouse_cells_original
  combined_data <- c(human_cells, mouse_cells_flipped)
  combined_data_clean <- combined_data[combined_data > 0 & combined_data < 1]

  if (length(combined_data_clean) == 0) {
    stop("No valid data points for distribution fitting after flipping")
  }

  # Fit distribution
  if (distribution == "auto") {
    dist_names <- c("lnorm", "gamma", "norm", "exp")
    best_aic <- Inf
    best_fit <- NULL
    best_dist <- NULL

    for (dist in dist_names) {
      fit <- tryCatch({
        fitdist(combined_data_clean, dist)
      }, error = function(e) NULL)

      if (!is.null(fit) && fit$aic < best_aic) {
        best_aic <- fit$aic
        best_fit <- fit
        best_dist <- dist
      }
    }
  } else {
    best_fit <- tryCatch({
      fitdist(combined_data_clean, distribution)
    }, error = function(e) {
      stop(paste("Failed to fit", distribution, "distribution"))
    })
    best_dist <- distribution
  }

  if (is.null(best_fit)) {
    stop("Failed to fit any distribution to the flipped data")
  }

  # Calculate threshold
  threshold_flipped <- calculate_multiplet_threshold_viz(best_fit, best_dist, n_sd)
  threshold_flipped <- min(max(threshold_flipped, 0.001), 0.999)
  human_threshold <- threshold_flipped
  mouse_threshold <- 1 - threshold_flipped

  # Calculate statistics for display
  class_counts <- table(data$multiplet_classification, useNA = "ifany")
  n_human <- ifelse("Human" %in% names(class_counts), class_counts["Human"], 0)
  n_mouse <- ifelse("Mouse" %in% names(class_counts), class_counts["Mouse"], 0)
  n_multiplet <- ifelse("Multiplet" %in% names(class_counts), class_counts["Multiplet"], 0)
  n_invalid <- ifelse("Invalid" %in% names(class_counts), class_counts["Invalid"], 0)

  human_percentage <- round(n_human / nrow(data) * 100, 1)
  mouse_percentage <- round(n_mouse / nrow(data) * 100, 1)
  multiplet_percentage <- round(n_multiplet / nrow(data) * 100, 1)

  # Save current par settings
  old_par <- par(no.readonly = TRUE)

  if (plot_type == "standard") {
    # Standard 2x2 layout
    par(mfrow = c(2, 2),
        mar = c(4, 4, 3, 2),
        cex.main = 1.3,
        cex.lab = 1.1,
        cex.axis = 1.0)

    # Plot 1: Original data with thresholds and colored histogram
    hist_data <- hist(data$percentMouse[valid_indices], breaks = 100, plot = FALSE)

    # Create colored bars based on thresholds
    bar_colors <- ifelse(hist_data$mids < human_threshold,
                         rgb(0, 0, 1, 0.6),  # Blue for human
                         ifelse(hist_data$mids > mouse_threshold,
                                rgb(1, 0, 0, 0.6),  # Red for mouse
                                rgb(0.5, 0, 0.5, 0.6)))  # Purple for multiplet zone

    hist(data$percentMouse[valid_indices], breaks = 100, col = bar_colors, border = "white",
         main = paste0("Step 3: Multiplet Classification (σ = ", n_sd, ")"),
         xlab = "Percent Mouse", ylab = "Frequency", xlim = c(0, 1))

    # Add lines
    abline(v = 0.5, col = "black", lwd = 2, lty = 3)
    abline(v = human_threshold, col = "darkblue", lwd = 2.5, lty = 2)
    abline(v = mouse_threshold, col = "darkred", lwd = 2.5, lty = 2)

    # Add zone labels
    text(human_threshold/2, par("usr")[4] * 0.9, "Human", col = "darkblue", cex = 1.1, font = 2)
    text((1 + mouse_threshold)/2, par("usr")[4] * 0.9, "Mouse", col = "darkred", cex = 1.1, font = 2)
    if (mouse_threshold - human_threshold > 0.05) {
      text(0.5, par("usr")[4] * 0.9, "Multiplet", col = "purple", cex = 1.1, font = 2)
    }

    # Add legend with counts
    legend("top",
           legend = c(paste("Human:", n_human, "(", human_percentage, "%)"),
                      paste("Mouse:", n_mouse, "(", mouse_percentage, "%)"),
                      paste("Multiplet:", n_multiplet, "(", multiplet_percentage, "%)")),
           fill = c(rgb(0, 0, 1, 0.6), rgb(1, 0, 0, 0.6), rgb(0.5, 0, 0.5, 0.6)),
           cex = 0.9, bty = "n")

    # Plot 2: Flipped data distribution with fit
    hist(combined_data_clean, breaks = 50, col = "lightcyan", probability = TRUE,
         main = "Flipped Data Distribution with Fit",
         xlab = "Flipped Values", ylab = "Density")

    # Add fitted distribution curve
    x_seq <- seq(min(combined_data_clean), max(combined_data_clean), length.out = 1000)
    y_fitted <- get_fitted_density(x_seq, best_fit, best_dist)
    lines(x_seq, y_fitted, col = "red", lwd = 2)

    abline(v = threshold_flipped, col = "purple", lwd = 2.5, lty = 2)

    # Add annotations
    text(threshold_flipped, par("usr")[4] * 0.9,
         paste("Threshold\n", round(threshold_flipped, 3)),
         col = "purple", pos = 4, cex = 1.1, font = 2)

    legend("topright",
           legend = c("Data", paste("Fitted", best_dist), "Threshold"),
           col = c("lightcyan", "red", "purple"),
           lty = c(NA, 1, 2),
           lwd = c(NA, 2, 2.5),
           pch = c(15, NA, NA),
           cex = 0.9)

    # Plot 3: Classification results bar chart
    # Order the bars logically
    ordered_counts <- class_counts[c("Human", "Mouse", "Multiplet", "Invalid")]
    ordered_counts <- ordered_counts[!is.na(ordered_counts)]

    bar_colors_plot <- c("Human" = rgb(0, 0, 1, 0.8),
                         "Mouse" = rgb(1, 0, 0, 0.8),
                         "Multiplet" = rgb(0.5, 0, 0.5, 0.8),
                         "Invalid" = rgb(0.5, 0.5, 0.5, 0.8))

    barplot(ordered_counts,
            col = bar_colors_plot[names(ordered_counts)],
            main = "Classification Results",
            ylab = "Number of Cells",
            ylim = c(0, max(ordered_counts) * 1.1))

    # Add count labels on bars
    text(x = seq(0.7, by = 1.2, length.out = length(ordered_counts)),
         y = ordered_counts + max(ordered_counts) * 0.02,
         labels = paste0(ordered_counts, "\n(",
                         round(ordered_counts/sum(ordered_counts)*100, 1), "%)"),
         cex = 0.9)

    # Plot 4: Scatter plot colored by classification
    # Calculate total reads if not present
    if (!"totalReads" %in% colnames(data)) {
      data$totalReads <- data$GRCh38 + data$GRCm39
    }

    # Define colors for each classification
    colors_scatter <- c("Human" = rgb(0, 0, 1, 0.5),
                        "Mouse" = rgb(1, 0, 0, 0.5),
                        "Multiplet" = rgb(0.5, 0, 0.5, 0.7),
                        "Invalid" = rgb(0.5, 0.5, 0.5, 0.3))

    plot(data$totalReads[valid_indices],
         data$percentMouse[valid_indices],
         col = colors_scatter[data$multiplet_classification[valid_indices]],
         pch = 16, cex = 0.5,
         main = "Classification by Total Reads vs Percent Mouse",
         xlab = "Total Reads", ylab = "Percent Mouse",
         xlim = c(0, quantile(data$totalReads, 0.99, na.rm = TRUE)),
         ylim = c(0, 1))

    # Add threshold lines
    abline(h = c(human_threshold, mouse_threshold), col = "gray40", lty = 2, lwd = 1.5)
    abline(h = 0.5, col = "gray60", lty = 3, lwd = 1)

    # Add shaded multiplet zone
    rect(0, human_threshold, par("usr")[2], mouse_threshold,
         col = rgb(0.5, 0.5, 0.5, 0.1), border = NA)

    legend("right",
           legend = names(colors_scatter)[names(colors_scatter) %in% unique(data$multiplet_classification)],
           col = colors_scatter[names(colors_scatter) %in% unique(data$multiplet_classification)],
           pch = 16,
           pt.cex = 1.2,
           cex = 0.9,
           title = "Classification")

  } else if (plot_type == "single") {
    # Single histogram plot
    par(mar = c(5, 4, 4, 2),
        cex.main = 1.4,
        cex.lab = 1.2,
        cex.axis = 1.1)

    # Create colored histogram
    hist_data <- hist(data$percentMouse[valid_indices], breaks = 100, plot = FALSE)

    bar_colors <- ifelse(hist_data$mids < human_threshold,
                         rgb(0, 0, 1, 0.6),  # Blue for human
                         ifelse(hist_data$mids > mouse_threshold,
                                rgb(1, 0, 0, 0.6),  # Red for mouse
                                rgb(0.5, 0, 0.5, 0.6)))  # Purple for multiplet

    hist(data$percentMouse[valid_indices], breaks = 100, col = bar_colors, border = "white",
         main = paste0("Step 3: Multiplet Detection (σ = ", n_sd, ")"),
         xlab = "Percent Mouse", ylab = "Frequency", xlim = c(0, 1))

    # Add lines
    abline(v = 0.5, col = "black", lwd = 2, lty = 3)
    abline(v = human_threshold, col = "darkblue", lwd = 3, lty = 2)
    abline(v = mouse_threshold, col = "darkred", lwd = 3, lty = 2)

    # Add annotations
    text(human_threshold/2, par("usr")[4] * 0.9, "Human", col = "darkblue", cex = 1.2, font = 2)
    text((1 + mouse_threshold)/2, par("usr")[4] * 0.9, "Mouse", col = "darkred", cex = 1.2, font = 2)
    if (mouse_threshold - human_threshold > 0.1) {
      text(0.5, par("usr")[4] * 0.7, "Multiplet\nZone", col = "purple", cex = 1.2, font = 2)
    }

    # Add subtitle with statistics
    mtext(paste0("Human: ", n_human, " (", human_percentage, "%) | ",
                 "Mouse: ", n_mouse, " (", mouse_percentage, "%) | ",
                 "Multiplet: ", n_multiplet, " (", multiplet_percentage, "%)"),
          side = 3, line = 0.3, cex = 0.9)

    # Add thresholds text
    text(human_threshold, par("usr")[3],
         paste0(round(human_threshold*100, 1), "%"),
         col = "darkblue", pos = 3, cex = 0.9, font = 2)
    text(mouse_threshold, par("usr")[3],
         paste0(round(mouse_threshold*100, 1), "%"),
         col = "darkred", pos = 3, cex = 0.9, font = 2)
  }

  # Restore par settings
  par(old_par)

  # Print summary if verbose
  if (verbose) {
    cat("\n=== Multiplet Classification Plot Summary ===\n")
    cat("Distribution used:", best_dist, "\n")
    cat("Sigma (σ):", n_sd, "\n")
    cat("Human threshold: <", round(human_threshold * 100, 1), "% mouse\n")
    cat("Mouse threshold: >", round(mouse_threshold * 100, 1), "% mouse\n")
    cat("Multiplet zone:", round(human_threshold * 100, 1), "% -",
        round(mouse_threshold * 100, 1), "%\n")
    cat("\nClassification counts:\n")
    cat("Human cells:", n_human, "/", nrow(data), "(", human_percentage, "%)\n")
    cat("Mouse cells:", n_mouse, "/", nrow(data), "(", mouse_percentage, "%)\n")
    cat("Multiplet cells:", n_multiplet, "/", nrow(data), "(", multiplet_percentage, "%)\n")
    if (n_invalid > 0) {
      cat("Invalid cells:", n_invalid, "\n")
    }
  }

  # Return useful information invisibly
  invisible(list(
    data = data,  # Return the potentially classified data
    human_threshold = human_threshold,
    mouse_threshold = mouse_threshold,
    threshold_flipped = threshold_flipped,
    distribution = best_dist,
    n_sd = n_sd,
    human_count = as.numeric(n_human),
    mouse_count = as.numeric(n_mouse),
    multiplet_count = as.numeric(n_multiplet),
    multiplet_percentage = multiplet_percentage
  ))
}

# =============================================================================
# HELPER FUNCTION FOR MULTIPLET THRESHOLD CALCULATION
# =============================================================================

calculate_multiplet_threshold_viz <- function(fit, distribution, n_sd) {
  if (distribution == "lnorm") {
    meanlog <- fit$estimate["meanlog"]
    sdlog <- fit$estimate["sdlog"]
    # Calculate threshold in flipped space
    threshold_flipped <- exp(meanlog + n_sd * sdlog)
  } else if (distribution == "gamma") {
    shape <- fit$estimate["shape"]
    rate <- fit$estimate["rate"]
    gamma_mean <- shape / rate
    gamma_sd <- sqrt(shape) / rate
    threshold_flipped <- gamma_mean + n_sd * gamma_sd
  } else if (distribution == "exp") {
    rate <- fit$estimate["rate"]
    exp_mean <- 1 / rate
    exp_sd <- 1 / rate
    threshold_flipped <- exp_mean + n_sd * exp_sd
  } else if (distribution == "norm") {
    norm_mean <- fit$estimate["mean"]
    norm_sd <- fit$estimate["sd"]
    threshold_flipped <- norm_mean + n_sd * norm_sd
  } else {
    stop(paste("Unsupported distribution:", distribution))
  }
  return(threshold_flipped)
}

# =============================================================================
# FUNCTION: PLOT COMBINED CLASSIFICATION (UPDATED)
# =============================================================================

plot_combined_classification <- function(data,
                                         low_sd = 1.0,
                                         high_sd = 1.5,
                                         multiplet_sd = 3.0,
                                         xlim_max = NULL,
                                         show_stats = TRUE,
                                         show_legend = TRUE,
                                         point_size = 0.5,
                                         point_alpha = 0.6,
                                         title = NULL,
                                         run_missing_classifications = FALSE,
                                         verbose = FALSE) {

  # Create combined classification plot showing all available classification steps
  # Parameters:
  # - data: dataframe with at least one classification (low/high/multiplet)
  # - low_sd: sigma value used/to use for low reads classification
  # - high_sd: sigma value used/to use for high reads classification
  # - multiplet_sd: sigma value used/to use for multiplet classification
  # - xlim_max: maximum x-axis limit for total reads (auto if NULL)
  # - show_stats: whether to show statistics box
  # - show_legend: whether to show classification legend
  # - point_size: size of points (single value or vector)
  # - point_alpha: transparency of points (single value or vector)
  # - title: main title for the plot (auto-generated if NULL)
  # - run_missing_classifications: if TRUE, run any missing classifications
  # - verbose: whether to print messages

  # Standardize column names first
  data <- standardize_column_names(data)

  # Validate required base columns exist
  required_cols <- c("GRCh38", "GRCm39")
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", "),
               "\nPlease ensure your data has human and mouse read count columns"))
  }

  # Calculate derived columns if not present
  if (!"totalReads" %in% colnames(data)) {
    data$totalReads <- data$GRCh38 + data$GRCm39
    if (verbose) cat("Calculated totalReads column\n")
  }
  if (!"percentMouse" %in% colnames(data)) {
    data$percentMouse <- data$GRCm39 / data$totalReads
    data$percentMouse[!is.finite(data$percentMouse)] <- 0
    if (verbose) cat("Calculated percentMouse column\n")
  }

  # Check which classifications are available
  has_low <- "low_classification" %in% colnames(data)
  has_high <- "high_classification" %in% colnames(data)
  has_multiplet <- "multiplet_classification" %in% colnames(data)

  # Count available classifications
  n_classifications <- sum(c(has_low, has_high, has_multiplet))

  if (n_classifications == 0) {
    if (!run_missing_classifications) {
      stop("No classification columns found. Set run_missing_classifications = TRUE to run all classifications,\n",
           "or run at least one classification function first.")
    } else {
      message("No classifications found. Running all three classification steps...")
    }
  } else {
    if (verbose) {
      cat("\nAvailable classifications:\n")
      if (has_low) cat("✓ Low reads classification\n")
      if (has_high) cat("✓ High reads classification\n")
      if (has_multiplet) cat("✓ Multiplet classification\n")
    }
  }

  # Run missing classifications if requested
  if (run_missing_classifications) {
    if (!has_low) {
      if (verbose) cat("Running low reads classification with σ =", low_sd, "\n")
      data <- classify_low_reads(data, n_sd = low_sd, verbose = FALSE)
      has_low <- TRUE
    }
    if (!has_high) {
      if (verbose) cat("Running high reads classification with σ =", high_sd, "\n")
      data <- classify_high_reads(data, n_sd = high_sd, verbose = FALSE)
      has_high <- TRUE
    }
    if (!has_multiplet) {
      if (verbose) cat("Running multiplet classification with σ =", multiplet_sd, "\n")
      data <- classify_multiplets(data, n_sd = multiplet_sd, verbose = FALSE)
      has_multiplet <- TRUE
    }
    n_classifications <- 3  # All classifications now available
  }

  # Generate title if not provided
  if (is.null(title)) {
    if (n_classifications == 3) {
      title <- "Complete Three-Step Classification Analysis"
    } else if (n_classifications == 2) {
      steps <- c()
      if (has_low) steps <- c(steps, "Low")
      if (has_high) steps <- c(steps, "High")
      if (has_multiplet) steps <- c(steps, "Multiplet")
      title <- paste("Classification Analysis:", paste(steps, collapse = " & "))
    } else if (n_classifications == 1) {
      if (has_low) title <- "Low Reads Classification Only"
      else if (has_high) title <- "High Reads Classification Only"
      else if (has_multiplet) title <- "Multiplet Classification Only"
    }
  }

  # Calculate thresholds only for available classifications
  thresholds <- list()

  if (has_low) {
    thresholds$low <- get_threshold_value(data, "low", low_sd)
  }
  if (has_high) {
    thresholds$high <- get_threshold_value(data, "high", high_sd)
  }
  if (has_multiplet) {
    mult_thresh <- get_multiplet_thresholds(data, multiplet_sd)
    thresholds$human <- mult_thresh$human
    thresholds$mouse <- mult_thresh$mouse
  }

  # Set xlim if not provided
  if (is.null(xlim_max)) {
    xlim_max <- quantile(data$totalReads, 0.99, na.rm = TRUE)
  }

  # Create combined classification status column
  data$combined_status <- "Normal"  # Default

  # Apply classifications in order of priority
  if (has_low) {
    data$combined_status[data$low_classification == "Low"] <- "Low"
  }
  if (has_high) {
    data$combined_status[data$high_classification == "High"] <- "High"
  }
  if (has_multiplet) {
    data$combined_status[data$multiplet_classification == "Multiplet"] <- "Multiplet"
    # Handle overlap cases
    if (has_low) {
      data$combined_status[data$low_classification == "Low" &
                             data$multiplet_classification == "Multiplet"] <- "Low+Multiplet"
    }
    if (has_high) {
      data$combined_status[data$high_classification == "High" &
                             data$multiplet_classification == "Multiplet"] <- "High+Multiplet"
    }
  }

  # Define colors for each status
  status_colors <- list(
    "Normal" = rgb(0.5, 0.5, 0.5, point_alpha * 0.5),
    "Low" = rgb(0, 0, 1, point_alpha),
    "High" = rgb(0, 0.7, 0, point_alpha),
    "Multiplet" = rgb(0.5, 0, 0.5, point_alpha),
    "Low+Multiplet" = rgb(1, 0.5, 0, point_alpha),
    "High+Multiplet" = rgb(1, 0, 0, point_alpha)
  )

  # Save current par settings
  old_par <- par(no.readonly = TRUE)
  par(mar = c(5, 4, 4, 2))

  # Create base plot
  plot(data$totalReads, data$percentMouse,
       type = "n",
       main = title,
       xlab = "Total Reads", ylab = "Percent Mouse",
       xlim = c(0, xlim_max), ylim = c(0, 1),
       cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.1)

  # Add grid for better readability
  grid(col = "lightgray", lty = "dotted", lwd = 0.5)

  # Add threshold lines based on what's available
  if (has_low && !is.null(thresholds$low)) {
    abline(v = thresholds$low, col = "blue", lwd = 2, lty = 3)
    text(thresholds$low, 0.02, "Low", col = "blue", cex = 0.9, pos = 4, srt = 90)
  }

  if (has_high && !is.null(thresholds$high)) {
    abline(v = thresholds$high, col = "darkgreen", lwd = 2, lty = 3)
    text(thresholds$high, 0.02, "High", col = "darkgreen", cex = 0.9, pos = 4, srt = 90)
  }

  if (has_multiplet) {
    if (!is.null(thresholds$human)) {
      abline(h = thresholds$human, col = "darkblue", lwd = 2, lty = 2)
      text(xlim_max * 0.02, thresholds$human - 0.02, "Human", col = "darkblue", cex = 0.9, pos = 3)
    }
    if (!is.null(thresholds$mouse)) {
      abline(h = thresholds$mouse, col = "darkred", lwd = 2, lty = 2)
      text(xlim_max * 0.02, thresholds$mouse + 0.02, "Mouse", col = "darkred", cex = 0.9, pos = 1)
    }
    # Shade multiplet zone
    if (!is.null(thresholds$human) && !is.null(thresholds$mouse)) {
      rect(0, thresholds$human, xlim_max, thresholds$mouse,
           col = rgb(0.5, 0, 0.5, 0.05), border = NA)
    }
  }

  # Always add 50% reference line
  abline(h = 0.5, col = "gray40", lwd = 1, lty = 3)

  # Plot points colored by combined status
  for (status in unique(data$combined_status)) {
    idx <- which(data$combined_status == status)
    if (length(idx) > 0) {
      # Adjust point size based on status
      pt_size <- if (status == "Normal") point_size * 0.8 else point_size
      points(data$totalReads[idx], data$percentMouse[idx],
             pch = 16, cex = pt_size,
             col = status_colors[[status]])
    }
  }

  # Add legend if requested
  if (show_legend) {
    # Get actual statuses present in data
    present_statuses <- unique(data$combined_status)
    present_statuses <- present_statuses[order(match(present_statuses,
                                                     c("Normal", "Low", "High", "Multiplet",
                                                       "Low+Multiplet", "High+Multiplet")))]

    # Create legend labels with counts
    legend_labels <- sapply(present_statuses, function(s) {
      count <- sum(data$combined_status == s)
      pct <- round(count / nrow(data) * 100, 1)
      paste0(s, " (n=", count, ", ", pct, "%)")
    })

    legend("topright",
           legend = legend_labels,
           col = unlist(status_colors[present_statuses]),
           pch = 16,
           pt.cex = 1.2,
           cex = 0.8,
           bg = "white",
           title = "Cell Classifications")
  }

  # Add statistics box if requested
  if (show_stats) {
    total_cells <- nrow(data)

    # Build stats text based on available classifications
    stats_lines <- c(paste0("Total cells: ", format(total_cells, big.mark=",")))

    if (has_low) {
      low_count <- sum(data$low_classification == "Low", na.rm = TRUE)
      stats_lines <- c(stats_lines,
                       paste0("Low reads: ", low_count, " (",
                              round(low_count/total_cells*100, 1), "%)"))
    }

    if (has_high) {
      high_count <- sum(data$high_classification == "High", na.rm = TRUE)
      stats_lines <- c(stats_lines,
                       paste0("High reads: ", high_count, " (",
                              round(high_count/total_cells*100, 1), "%)"))
    }

    if (has_multiplet) {
      multiplet_count <- sum(data$multiplet_classification == "Multiplet", na.rm = TRUE)
      human_count <- sum(data$multiplet_classification == "Human", na.rm = TRUE)
      mouse_count <- sum(data$multiplet_classification == "Mouse", na.rm = TRUE)

      stats_lines <- c(stats_lines,
                       paste0("Multiplets: ", multiplet_count, " (",
                              round(multiplet_count/total_cells*100, 1), "%)"),
                       paste0("Human: ", human_count, " (",
                              round(human_count/total_cells*100, 1), "%)"),
                       paste0("Mouse: ", mouse_count, " (",
                              round(mouse_count/total_cells*100, 1), "%)"))
    }

    # Add overlap counts if relevant
    if (has_low && has_multiplet) {
      low_mult <- sum(data$combined_status == "Low+Multiplet")
      if (low_mult > 0) {
        stats_lines <- c(stats_lines, paste0("Low+Mult: ", low_mult))
      }
    }

    if (has_high && has_multiplet) {
      high_mult <- sum(data$combined_status == "High+Multiplet")
      if (high_mult > 0) {
        stats_lines <- c(stats_lines, paste0("High+Mult: ", high_mult))
      }
    }

    legend("bottomleft",
           legend = paste(stats_lines, collapse = "\n"),
           bty = "n",
           cex = 0.8,
           title = "Statistics")
  }

  # Add note about missing classifications if applicable
  if (n_classifications < 3 && !run_missing_classifications) {
    missing_steps <- c()
    if (!has_low) missing_steps <- c(missing_steps, "Low")
    if (!has_high) missing_steps <- c(missing_steps, "High")
    if (!has_multiplet) missing_steps <- c(missing_steps, "Multiplet")

    mtext(paste("Note: Missing", paste(missing_steps, collapse = ", "), "classification(s)"),
          side = 3, line = 0.3, cex = 0.8, col = "gray40")
  }

  # Restore par settings
  par(old_par)

  # Print summary if verbose
  if (verbose) {
    cat("\n=== Combined Classification Summary ===\n")
    cat("Classifications available:", n_classifications, "of 3\n")
    if (has_low) cat("✓ Low reads (σ =", low_sd, ")\n")
    if (has_high) cat("✓ High reads (σ =", high_sd, ")\n")
    if (has_multiplet) cat("✓ Multiplets (σ =", multiplet_sd, ")\n")
    cat("\nCell counts by status:\n")
    print(table(data$combined_status))
  }

  # Return useful information invisibly
  invisible(list(
    data = data,
    n_classifications = n_classifications,
    has_low = has_low,
    has_high = has_high,
    has_multiplet = has_multiplet,
    thresholds = thresholds,
    status_counts = table(data$combined_status)
  ))
}

# =============================================================================
# HELPER FUNCTIONS FOR THRESHOLD EXTRACTION
# =============================================================================

get_threshold_value <- function(data, type = c("low", "high"), n_sd) {
  # Extract threshold value for low or high reads classification

  type <- match.arg(type)

  # Standardize column names
  data <- standardize_column_names(data)

  # Get valid reads
  if (!"totalReads" %in% colnames(data)) {
    data$totalReads <- data$GRCh38 + data$GRCm39
  }
  valid_reads <- data$totalReads[is.finite(data$totalReads) & data$totalReads > 0]

  if (length(valid_reads) == 0) {
    warning("No valid reads for threshold calculation")
    return(NA)
  }

  # Try to fit distribution (with error handling)
  fit <- tryCatch({
    fitdist(valid_reads, "weibull")
  }, error = function(e) {
    warning("Could not fit Weibull distribution, using simple mean/sd")
    return(NULL)
  })

  if (is.null(fit)) {
    # Fallback to simple statistics
    threshold <- if (type == "low") {
      mean(valid_reads) - n_sd * sd(valid_reads)
    } else {
      mean(valid_reads) + n_sd * sd(valid_reads)
    }
  } else {
    # Use Weibull fit
    shape <- fit$estimate["shape"]
    scale <- fit$estimate["scale"]

    # Calculate parameters
    weibull_mean <- scale * gamma(1 + 1/shape)
    weibull_var <- scale^2 * (gamma(1 + 2/shape) - gamma(1 + 1/shape)^2)
    weibull_sd <- sqrt(weibull_var)

    # Calculate threshold based on type
    threshold <- if (type == "low") {
      weibull_mean - n_sd * weibull_sd
    } else {
      weibull_mean + n_sd * weibull_sd
    }
  }

  # Ensure threshold is positive for low, reasonable for high
  if (type == "low") {
    threshold <- max(threshold, 0)
  }

  return(threshold)
}

get_multiplet_thresholds <- function(data, n_sd) {
  # Extract human and mouse thresholds for multiplet classification

  # Standardize column names
  data <- standardize_column_names(data)

  # Ensure percentMouse exists
  if (!"percentMouse" %in% colnames(data)) {
    data$percentMouse <- data$GRCm39 / (data$GRCh38 + data$GRCm39)
  }

  # Get valid data
  valid_indices <- is.finite(data$percentMouse)
  valid_pm <- data$percentMouse[valid_indices]

  if (length(valid_pm) == 0) {
    warning("No valid percentMouse values for threshold calculation")
    return(list(human = NA, mouse = NA))
  }

  # Apply flipping logic
  human_pm <- valid_pm[valid_pm <= 0.5]
  mouse_pm_original <- valid_pm[valid_pm > 0.5]
  mouse_pm_flipped <- 1 - mouse_pm_original
  combined_flipped <- c(human_pm, mouse_pm_flipped)
  combined_clean <- combined_flipped[combined_flipped > 0 & combined_flipped < 1]

  if (length(combined_clean) == 0) {
    warning("No valid data for multiplet threshold calculation")
    return(list(human = 0.1, mouse = 0.9))  # Default fallback
  }

  # Try to fit log-normal distribution
  fit <- tryCatch({
    fitdist(combined_clean, "lnorm")
  }, error = function(e) {
    warning("Could not fit log-normal distribution, using simple statistics")
    return(NULL)
  })

  if (is.null(fit)) {
    # Fallback to simple percentile-based thresholds
    threshold_flipped <- quantile(combined_clean, probs = pnorm(n_sd))
  } else {
    meanlog <- fit$estimate["meanlog"]
    sdlog <- fit$estimate["sdlog"]

    # Calculate threshold in flipped space
    threshold_flipped <- exp(meanlog + n_sd * sdlog)
  }

  # Ensure threshold is within bounds
  threshold_flipped <- min(max(threshold_flipped, 0.001), 0.999)

  # Convert back to original space
  human_threshold <- threshold_flipped
  mouse_threshold <- 1 - threshold_flipped

  return(list(human = human_threshold, mouse = mouse_threshold))
}







