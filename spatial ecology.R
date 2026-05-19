main <- function() {
  
  # 1. Paths and packages
  
  # Set the project folder. The raw community and environmental data are stored here.
  project_dir <- "D:/Beetles"
  
  # Input data files
  community_file <- file.path(project_dir, "scot_beetle_community.csv")
  environment_file <- file.path(project_dir, "scot_beetle_env.csv")
  
  # All outputs from this script will be saved in one folder.
  output_dir <- file.path(project_dir, "Hmsc_JSDM_outputs")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # A small helper for saving files into the output folder.
  out_path <- function(file_name) file.path(output_dir, file_name)
  
  # Check that the two input files exist before running the analysis.
  if (!file.exists(community_file)) {
    stop("Community data file not found: ", community_file)
  }
  
  if (!file.exists(environment_file)) {
    stop("Environmental data file not found: ", environment_file)
  }
  
  # The Hmsc package is required for the Bayesian JSDM.
  if (!requireNamespace("Hmsc", quietly = TRUE)) {
    stop("The Hmsc package is required for this analysis.")
  }
  
  library(Hmsc)
  
  # Set a seed so that random parts of the workflow are reproducible.
  set.seed(123)
  
  # 2. Functions used later in the workflow
  
  # Some CSV files include an extra first column from row names.
  # This function removes that column if it looks like a simple row index.
  drop_index_column <- function(dat) {
    if (ncol(dat) == 0) return(dat)
    
    first_name <- names(dat)[1]
    first_values <- dat[[1]]
    first_numeric <- suppressWarnings(as.numeric(as.character(first_values)))
    
    looks_like_row_index <- length(first_numeric) == nrow(dat) &&
      all(!is.na(first_numeric)) &&
      all(first_numeric == seq_len(nrow(dat)))
    
    if (is.na(first_name) || first_name == "" || first_name %in% c("X", "...1")) {
      if (first_name == "" || looks_like_row_index) {
        dat <- dat[, -1, drop = FALSE]
      }
    }
    
    dat
  }
  
  # General heatmap function for matrix outputs.
  # It is used for posterior effects and species association matrices.
  make_heatmap <- function(mat, file_name, title_text, zlim = NULL,
                           value_labels = FALSE, legend_title = "Value") {
    mat <- as.matrix(mat)
    nr <- nrow(mat)
    nc <- ncol(mat)
    
    # Blue is used for negative values, white for values near zero,
    # and red for positive values.
    colours <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101)
    
    # If no colour range is given, use a symmetric range around zero.
    if (is.null(zlim)) {
      finite_values <- mat[is.finite(mat)]
      max_abs <- if (length(finite_values) == 0) 1 else max(abs(finite_values), na.rm = TRUE)
      if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
      zlim <- c(-max_abs, max_abs)
    }
    
    png(out_path(file_name), width = 1800, height = max(1200, 70 * nr + 500), res = 160)
    old_par <- par(no.readonly = TRUE)
    
    on.exit({
      par(old_par)
      dev.off()
    }, add = TRUE)
    
    par(mar = c(8, 9, 5, 6), xpd = NA)
    plot.new()
    plot.window(xlim = c(0, nc + 1.2), ylim = c(0, nr), xaxs = "i", yaxs = "i")
    
    breaks <- seq(zlim[1], zlim[2], length.out = length(colours) + 1)
    
    # Draw each matrix cell manually so that the layout is easy to control.
    for (i in seq_len(nr)) {
      for (j in seq_len(nc)) {
        y_bottom <- nr - i
        value <- mat[i, j]
        fill <- "grey85"
        
        if (is.finite(value)) {
          colour_index <- findInterval(value, breaks, all.inside = TRUE)
          fill <- colours[colour_index]
        }
        
        rect(j - 1, y_bottom, j, y_bottom + 1, col = fill, border = "white")
        
        # Add values inside cells only when the matrix is not too large.
        if (value_labels && is.finite(value)) {
          text(j - 0.5, y_bottom + 0.5, labels = round(value, 2), cex = 0.65)
        }
      }
    }
    
    axis(1, at = seq(0.5, nc - 0.5, by = 1),
         labels = colnames(mat), las = 2, tick = FALSE)
    axis(2, at = seq(nr - 0.5, 0.5, by = -1),
         labels = rownames(mat), las = 1, tick = FALSE)
    box()
    title(main = title_text, font.main = 2)
    
    # Add a simple colour bar on the right side of the heatmap.
    legend_x <- nc + 0.25
    legend_y <- seq(0.5, nr - 0.5, length.out = length(colours))
    legend_height <- if (length(legend_y) > 1) abs(diff(legend_y)[1]) else 0.2
    
    for (k in seq_along(colours)) {
      rect(
        legend_x,
        legend_y[k] - legend_height / 2,
        legend_x + 0.18,
        legend_y[k] + legend_height / 2,
        col = colours[k],
        border = NA
      )
    }
    
    text(legend_x + 0.35, 0.5, labels = format(round(zlim[1], 2), nsmall = 2), adj = 0, cex = 0.8)
    text(legend_x + 0.35, nr / 2, labels = "0", adj = 0, cex = 0.8)
    text(legend_x + 0.35, nr - 0.5, labels = format(round(zlim[2], 2), nsmall = 2), adj = 0, cex = 0.8)
    text(legend_x, nr + 0.25, labels = legend_title, adj = 0, font = 2, cex = 0.85)
  }
  
  # This heatmap shows whether posterior effects are strongly supported
  # as positive, negative, or not strongly supported.
  make_supported_heatmap <- function(supported_mat, file_name) {
    supported_mat <- as.matrix(supported_mat)
    nr <- nrow(supported_mat)
    nc <- ncol(supported_mat)
    
    colours <- c(
      "supported negative" = "#2166AC",
      "not strongly supported" = "#F7F7F7",
      "supported positive" = "#B2182B"
    )
    
    png(out_path(file_name), width = 1800, height = max(1200, 70 * nr + 500), res = 160)
    old_par <- par(no.readonly = TRUE)
    
    on.exit({
      par(old_par)
      dev.off()
    }, add = TRUE)
    
    par(mar = c(8, 9, 5, 7), xpd = NA)
    plot.new()
    plot.window(xlim = c(0, nc + 1.4), ylim = c(0, nr), xaxs = "i", yaxs = "i")
    
    for (i in seq_len(nr)) {
      for (j in seq_len(nc)) {
        y_bottom <- nr - i
        label <- supported_mat[i, j]
        fill <- colours[label]
        
        if (is.na(fill)) fill <- "grey85"
        
        rect(j - 1, y_bottom, j, y_bottom + 1, col = fill, border = "white")
      }
    }
    
    axis(1, at = seq(0.5, nc - 0.5, by = 1),
         labels = colnames(supported_mat), las = 2, tick = FALSE)
    axis(2, at = seq(nr - 0.5, 0.5, by = -1),
         labels = rownames(supported_mat), las = 1, tick = FALSE)
    box()
    title(main = "Strong posterior-supported Hmsc effects", font.main = 2)
    
    legend(
      "right",
      inset = c(-0.35, 0),
      legend = names(colours),
      fill = colours,
      border = NA,
      bty = "n",
      cex = 0.85
    )
  }
  
  # Calculate AUC from observed 0/1 values and predicted probabilities.
  # This is used as a simple species-level predictive performance measure.
  simple_auc <- function(obs, score) {
    ok <- is.finite(obs) & is.finite(score)
    obs <- obs[ok]
    score <- score[ok]
    
    if (length(unique(obs)) < 2) return(NA_real_)
    
    n_pos <- sum(obs == 1)
    n_neg <- sum(obs == 0)
    
    if (n_pos == 0 || n_neg == 0) return(NA_real_)
    
    ranks <- rank(score, ties.method = "average")
    (sum(ranks[obs == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
  }
  
  # Summarise the overall strength of residual species associations.
  # The label is only used as a compact description of the association matrix.
  association_strength_label <- function(values) {
    values <- values[is.finite(values)]
    
    if (length(values) == 0) return("not available")
    
    median_abs <- median(abs(values), na.rm = TRUE)
    max_abs <- max(abs(values), na.rm = TRUE)
    
    if (median_abs < 0.10 && max_abs < 0.40) return("generally weak")
    if (median_abs < 0.20 && max_abs < 0.60) return("weak to moderate")
    if (median_abs < 0.30) return("moderate, with some stronger pairs")
    
    "relatively strong"
  }
  
  # Extract beta samples from all MCMC chains in the fitted Hmsc model.
  # These samples are later used to calculate posterior means and intervals.
  extract_beta_samples <- function(model) {
    beta_list <- list()
    index <- 1
    
    for (chain_i in seq_along(model$postList)) {
      for (sample_i in seq_along(model$postList[[chain_i]])) {
        beta_list[[index]] <- model$postList[[chain_i]][[sample_i]]$Beta
        index <- index + 1
      }
    }
    
    simplify2array(beta_list)
  }
  
  # Fit an Hmsc model while keeping track of warnings, errors and run time.
  # This makes it easier to compare the fixed-effect and random-effect models.
  fit_hmsc_with_settings <- function(model_object, model_label, settings) {
    warnings_found <- character()
    start_time <- Sys.time()
    
    fitted <- tryCatch(
      withCallingHandlers(
        sampleMcmc(
          model_object,
          samples = settings$samples,
          transient = settings$transient,
          thin = settings$thin,
          nChains = settings$nChains,
          verbose = 1
        ),
        warning = function(w) {
          warnings_found <<- c(warnings_found, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) e
    )
    
    end_time <- Sys.time()
    
    list(
      model = fitted,
      success = inherits(fitted, "Hmsc"),
      warnings = unique(warnings_found),
      error = if (inherits(fitted, "error")) conditionMessage(fitted) else "",
      label = model_label,
      settings = settings,
      elapsed_minutes = as.numeric(difftime(end_time, start_time, units = "mins"))
    )
  }
  
  # 3. Read and prepare the beetle data
  
  comm <- read.csv(community_file, check.names = FALSE, stringsAsFactors = FALSE)
  env <- read.csv(environment_file, check.names = FALSE, stringsAsFactors = FALSE)
  
  # Remove possible row-index columns before matching the two tables.
  comm <- drop_index_column(comm)
  env <- drop_index_column(env)
  
  # If both files have a Sites column, sort by it and check that the rows match.
  # Otherwise, the script assumes that the row order in both files is already correct.
  if ("Sites" %in% names(comm) && "Sites" %in% names(env)) {
    comm <- comm[order(comm$Sites), , drop = FALSE]
    env <- env[order(env$Sites), , drop = FALSE]
    
    site_labels_match <- identical(as.character(comm$Sites), as.character(env$Sites))
    if (!site_labels_match) stop("Sites columns did not match after sorting.")
    
    site_ids <- comm$Sites
    site_matching_method <- "Sites columns were sorted and matched."
  } else {
    site_ids <- seq_len(nrow(comm))
    site_matching_method <- "No shared Sites column; row order was used."
  }
  
  # Keep only the species columns for the community response matrix.
  species_abundance <- comm
  
  if ("Sites" %in% names(species_abundance)) {
    species_abundance <- species_abundance[, names(species_abundance) != "Sites", drop = FALSE]
  }
  
  # Make sure all species columns are numeric before converting to presence-absence.
  species_abundance[] <- lapply(species_abundance, function(x) {
    suppressWarnings(as.numeric(as.character(x)))
  })
  
  # The Hmsc model is fitted to occurrence data, so abundance is converted to 0/1.
  species_pa_full <- ifelse(species_abundance > 0, 1, 0)
  species_pa_full <- as.data.frame(species_pa_full, check.names = FALSE)
  rownames(species_pa_full) <- as.character(site_ids)
  
  # Check how often each species occurs and whether it has enough variation for modelling.
  prevalence_table_full <- data.frame(
    species = names(species_pa_full),
    presences = colSums(species_pa_full == 1, na.rm = TRUE),
    absences = colSums(species_pa_full == 0, na.rm = TRUE),
    prevalence = colMeans(species_pa_full == 1, na.rm = TRUE),
    variance = vapply(species_pa_full, var, numeric(1), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  # Species with no variation cannot be fitted in a presence-absence model.
  zero_variance_species <- prevalence_table_full$species[
    !is.finite(prevalence_table_full$variance) | prevalence_table_full$variance == 0
  ]
  
  included_species <- setdiff(names(species_pa_full), zero_variance_species)
  
  if (length(included_species) == 0) {
    stop("All species had zero variance, so Hmsc cannot be fitted.")
  }
  
  species_pa_nonconstant <- species_pa_full[, included_species, drop = FALSE]
  
  excluded_species_table <- data.frame(
    species = zero_variance_species,
    exclusion_reason = rep("zero variance in presence-absence response", length(zero_variance_species)),
    stringsAsFactors = FALSE
  )
  
  # Save the response matrix and the main species-level data checks.
  write.csv(species_pa_full, out_path("hmsc_presence_absence_matrix.csv"), row.names = TRUE)
  write.csv(prevalence_table_full, out_path("hmsc_species_prevalence_table.csv"), row.names = FALSE)
  write.csv(excluded_species_table, out_path("hmsc_excluded_species_table.csv"), row.names = FALSE)
