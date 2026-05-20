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

        # Add values inside cells only when needed.
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

  # 4. Select and scale environmental predictors

  # Use the same environmental variables described in the methods section.
  predictors <- c("Management", "Moist", "pH", "Elevation", "Bryophyte")

  # Check that all required predictors are present in the environmental table.
  missing_predictors <- setdiff(predictors, names(env))

  if (length(missing_predictors) > 0) {
    stop("Missing required predictors: ", paste(missing_predictors, collapse = ", "))
  }

  env_selected <- env[, predictors, drop = FALSE]

  # Convert predictor columns to numeric before scaling and modelling.
  env_selected[] <- lapply(env_selected, function(x) {
    suppressWarnings(as.numeric(as.character(x)))
  })

  # Keep only sites with complete predictor data.
  complete_rows <- complete.cases(env_selected)
  env_complete <- env_selected[complete_rows, , drop = FALSE]
  site_ids_model <- site_ids[complete_rows]

  # Match the response matrix to the same set of complete sites.
  Y <- as.matrix(species_pa_nonconstant[complete_rows, included_species, drop = FALSE])
  storage.mode(Y) <- "numeric"

  XData <- env_complete

  # Save the original means and standard deviations so the scaling step is clear.
  scale_details <- data.frame(
    predictor = predictors,
    mean_before_scaling = vapply(env_complete, mean, numeric(1), na.rm = TRUE),
    sd_before_scaling = vapply(env_complete, sd, numeric(1), na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  # Standardise predictors so coefficients are comparable across variables.
  XData[] <- lapply(XData, function(x) as.numeric(scale(x)))

  XData_out <- data.frame(Sites = site_ids_model, XData, check.names = FALSE)
  write.csv(XData_out, out_path("hmsc_environment_scaled.csv"), row.names = FALSE)
  write.csv(scale_details, out_path("hmsc_predictor_scaling_details.csv"), row.names = FALSE)

  # Check correlation among predictors before fitting the model.
  predictor_cor <- cor(XData, use = "pairwise.complete.obs")
  write.csv(predictor_cor, out_path("hmsc_predictor_correlation_matrix.csv"), row.names = TRUE)

  high_pairs <- data.frame(
    predictor_1 = character(),
    predictor_2 = character(),
    correlation = numeric(),
    stringsAsFactors = FALSE
  )

  # Record highly correlated predictor pairs for later checking.
  if (ncol(predictor_cor) >= 2) {
    for (i in seq_len(ncol(predictor_cor) - 1)) {
      for (j in (i + 1):ncol(predictor_cor)) {
        value <- predictor_cor[i, j]
        if (is.finite(value) && abs(value) > 0.7) {
          high_pairs <- rbind(
            high_pairs,
            data.frame(
              predictor_1 = rownames(predictor_cor)[i],
              predictor_2 = colnames(predictor_cor)[j],
              correlation = value,
              stringsAsFactors = FALSE
            )
          )
        }
      }
    }
  }

  write.csv(high_pairs, out_path("hmsc_high_correlation_pairs.csv"), row.names = FALSE)

  # Save a compact data check table for reproducibility.
  data_check_table <- data.frame(
    item = c(
      "community_file",
      "environment_file",
      "original_sites",
      "sites_used_after_predictor_filtering",
      "raw_species_count",
      "species_included_in_hmsc",
      "zero_variance_species_excluded",
      "rows_dropped_missing_predictors",
      "missing_abundance_after_numeric_conversion",
      "site_matching_method"
    ),
    value = c(
      community_file,
      environment_file,
      nrow(comm),
      nrow(Y),
      ncol(species_abundance),
      ncol(Y),
      length(zero_variance_species),
      sum(!complete_rows),
      any(is.na(species_abundance)),
      site_matching_method
    ),
    stringsAsFactors = FALSE
  )

  write.csv(data_check_table, out_path("hmsc_data_check_table.csv"), row.names = FALSE)


  # 5. Fit the Hmsc models

  x_formula <- ~ Management + Moist + pH + Elevation + Bryophyte

  # MCMC settings used for both model structures.
  mcmc_settings <- list(samples = 1000, transient = 500, thin = 10, nChains = 2)

  # Model A includes environmental fixed effects only.
  fixed_model <- Hmsc(
    Y = Y,
    XData = XData,
    XFormula = x_formula,
    distr = "probit"
  )

  fixed_fit <- fit_hmsc_with_settings(fixed_model, "Model A: fixed effects only", mcmc_settings)

  if (!fixed_fit$success) {
    stop("Model A failed. Error: ", fixed_fit$error)
  }

  saveRDS(fixed_fit$model, out_path("hmsc_model_fixed.rds"))

  # Model B adds a site-level random effect so that residual species associations can be estimated.
  studyDesign <- data.frame(site = factor(seq_len(nrow(Y))))
  ranLevels <- list(site = HmscRandomLevel(units = studyDesign$site))

  random_model <- Hmsc(
    Y = Y,
    XData = XData,
    XFormula = x_formula,
    distr = "probit",
    studyDesign = studyDesign,
    ranLevels = ranLevels
  )

  random_fit <- fit_hmsc_with_settings(random_model, "Model B: site-level random effect", mcmc_settings)

  # Use the random-effect model if it fits successfully; otherwise keep the fixed-effect model.
  if (random_fit$success) {
    saveRDS(random_fit$model, out_path("hmsc_model_random.rds"))
    selected_model <- random_fit$model
    selected_fit <- random_fit
    selected_model_name <- "Model B: site-level random effect"
    selected_model_file <- "hmsc_model_random.rds"
    selected_model_reason <- "site-level random effect available"
  } else {
    selected_model <- fixed_fit$model
    selected_fit <- fixed_fit
    selected_model_name <- "Model A: fixed effects only"
    selected_model_file <- "hmsc_model_fixed.rds"
    selected_model_reason <- "site-level random effect model did not fit"
  }

  # Save the model fitting summary for checking the selected model and MCMC settings.
  model_run_summary <- data.frame(
    model = c("Model A", "Model B", "Selected model"),
    description = c(
      "fixed effects only",
      "fixed effects with site-level random effect",
      selected_model_name
    ),
    success = c(fixed_fit$success, random_fit$success, TRUE),
    elapsed_minutes = c(fixed_fit$elapsed_minutes, random_fit$elapsed_minutes, NA_real_),
    samples = c(mcmc_settings$samples, mcmc_settings$samples, mcmc_settings$samples),
    transient = c(mcmc_settings$transient, mcmc_settings$transient, mcmc_settings$transient),
    thin = c(mcmc_settings$thin, mcmc_settings$thin, mcmc_settings$thin),
    nChains = c(mcmc_settings$nChains, mcmc_settings$nChains, mcmc_settings$nChains),
    warnings = c(
      paste(fixed_fit$warnings, collapse = " | "),
      paste(random_fit$warnings, collapse = " | "),
      ""
    ),
    error = c(fixed_fit$error, random_fit$error, ""),
    selected_model_file = c("", "", selected_model_file),
    selection_reason = c("", "", selected_model_reason),
    stringsAsFactors = FALSE
  )

  write.csv(model_run_summary, out_path("hmsc_model_run_summary.csv"), row.names = FALSE)


  # 6. Summarise posterior environmental effects

  post_beta <- getPostEstimate(selected_model, parName = "Beta")
  beta_samples <- extract_beta_samples(selected_model)

  # Summarise posterior beta samples for each predictor and species.
  beta_mean <- apply(beta_samples, c(1, 2), mean, na.rm = TRUE)
  beta_sd <- apply(beta_samples, c(1, 2), sd, na.rm = TRUE)
  beta_low <- apply(beta_samples, c(1, 2), quantile, probs = 0.025, na.rm = TRUE)
  beta_high <- apply(beta_samples, c(1, 2), quantile, probs = 0.975, na.rm = TRUE)
  beta_prob_pos <- apply(beta_samples, c(1, 2), function(x) mean(x > 0, na.rm = TRUE))
  beta_prob_neg <- apply(beta_samples, c(1, 2), function(x) mean(x < 0, na.rm = TRUE))

  cov_names <- selected_model$covNames

  if (length(cov_names) != dim(beta_samples)[1]) {
    cov_names <- rownames(post_beta$mean)
  }

  if (is.null(cov_names) || length(cov_names) != dim(beta_samples)[1]) {
    cov_names <- colnames(model.matrix(x_formula, data = XData))
  }

  species_names <- selected_model$spNames

  if (is.null(species_names) || length(species_names) != dim(beta_samples)[2]) {
    species_names <- colnames(Y)
  }

  dimnames(beta_mean) <- list(cov_names, species_names)
  dimnames(beta_sd) <- list(cov_names, species_names)
  dimnames(beta_low) <- list(cov_names, species_names)
  dimnames(beta_high) <- list(cov_names, species_names)
  dimnames(beta_prob_pos) <- list(cov_names, species_names)
  dimnames(beta_prob_neg) <- list(cov_names, species_names)

  # Convert the beta summaries into a long table for easier plotting and checking.
  beta_rows <- list()
  row_index <- 1

  for (pred in predictors) {
    pred_index <- match(pred, cov_names)
    if (is.na(pred_index)) next

    for (sp in species_names) {
      beta_rows[[row_index]] <- data.frame(
        species = sp,
        predictor = pred,
        posterior_mean = beta_mean[pred, sp],
        posterior_sd = beta_sd[pred, sp],
        ci_lower_95 = beta_low[pred, sp],
        ci_upper_95 = beta_high[pred, sp],
        posterior_prob_positive = beta_prob_pos[pred, sp],
        posterior_prob_negative = beta_prob_neg[pred, sp],
        ci_excludes_zero = beta_low[pred, sp] > 0 | beta_high[pred, sp] < 0,
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1
    }
  }

  beta_long <- do.call(rbind, beta_rows)
  write.csv(beta_long, out_path("hmsc_beta_summary_long.csv"), row.names = FALSE)

  # Also save a wide table so each species has one row.
  beta_wide <- data.frame(species = species_names, stringsAsFactors = FALSE)

  for (pred in predictors) {
    sub <- beta_long[beta_long$predictor == pred, , drop = FALSE]
    beta_wide[[paste0(pred, "_mean")]] <- sub$posterior_mean[match(species_names, sub$species)]
    beta_wide[[paste0(pred, "_lower95")]] <- sub$ci_lower_95[match(species_names, sub$species)]
    beta_wide[[paste0(pred, "_upper95")]] <- sub$ci_upper_95[match(species_names, sub$species)]
    beta_wide[[paste0(pred, "_prob_positive")]] <- sub$posterior_prob_positive[match(species_names, sub$species)]
  }

  write.csv(beta_wide, out_path("hmsc_beta_summary_wide.csv"), row.names = FALSE)

  # Prepare the matrix used for the posterior mean heatmap.
  beta_heatmap <- matrix(
    NA_real_,
    nrow = length(species_names),
    ncol = length(predictors),
    dimnames = list(species_names, predictors)
  )

  for (i in seq_len(nrow(beta_long))) {
    beta_heatmap[beta_long$species[i], beta_long$predictor[i]] <- beta_long$posterior_mean[i]
  }

  make_heatmap(
    beta_heatmap,
    "hmsc_beta_heatmap.png",
    "Hmsc posterior mean environmental effects",
    value_labels = TRUE,
    legend_title = "Posterior mean"
  )


  # 7. Summarise responses by predictor

  # Count how many species have positive or negative posterior mean effects for each predictor.
  predictor_summary <- do.call(rbind, lapply(predictors, function(pred) {
    sub <- beta_long[beta_long$predictor == pred, , drop = FALSE]

    data.frame(
      predictor = pred,
      number_positive_mean = sum(sub$posterior_mean > 0, na.rm = TRUE),
      number_negative_mean = sum(sub$posterior_mean < 0, na.rm = TRUE),
      number_ci_excludes_zero = sum(sub$ci_excludes_zero, na.rm = TRUE),
      mean_posterior_coefficient = mean(sub$posterior_mean, na.rm = TRUE),
      median_posterior_coefficient = median(sub$posterior_mean, na.rm = TRUE),
      mean_absolute_posterior_coefficient = mean(abs(sub$posterior_mean), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  write.csv(predictor_summary, out_path("hmsc_predictor_response_summary.csv"), row.names = FALSE)

  png(out_path("hmsc_predictor_response_counts.png"), width = 1500, height = 1000, res = 160)
  old_par <- par(no.readonly = TRUE)
  par(mar = c(7, 5, 4, 2))

  direction_counts <- rbind(
    Positive = predictor_summary$number_positive_mean,
    Negative = predictor_summary$number_negative_mean
  )
  colnames(direction_counts) <- predictor_summary$predictor

  barplot(
    direction_counts,
    beside = TRUE,
    col = c("#B2182B", "#2166AC"),
    las = 2,
    ylab = "Number of species",
    ylim = c(0, max(direction_counts, na.rm = TRUE) + 2),
    main = "Hmsc positive and negative posterior mean effects"
  )

  legend(
    "topright",
    legend = rownames(direction_counts),
    fill = c("#B2182B", "#2166AC"),
    bty = "n"
  )

  par(old_par)
  dev.off()
# 8. Calculate fitted probabilities and five-fold cross-validation

  # First calculate fitted probabilities from the final selected model.
  # These are used for the in-sample AUC comparison.
  pred_array <- try(computePredictedValues(selected_model), silent = TRUE)

  if (!inherits(pred_array, "try-error") && length(dim(pred_array)) == 3) {
    pred_probs <- apply(pred_array, c(1, 2), mean, na.rm = TRUE)
    prediction_method <- "computePredictedValues averaged across posterior samples"
  } else {
    X_matrix <- model.matrix(x_formula, data = XData)
    beta_for_prediction <- beta_mean[match(colnames(X_matrix), rownames(beta_mean)), , drop = FALSE]
    pred_probs <- pnorm(X_matrix %*% beta_for_prediction)
    prediction_method <- "posterior mean fixed effects with probit inverse link"
  }

  pred_probs <- as.matrix(pred_probs)
  colnames(pred_probs) <- species_names
  rownames(pred_probs) <- as.character(site_ids_model)

  write.csv(
    data.frame(Sites = site_ids_model, pred_probs, check.names = FALSE),
    out_path("hmsc_predicted_probabilities.csv"),
    row.names = FALSE
  )

  observed_prev <- colMeans(Y, na.rm = TRUE)
  predicted_prev <- colMeans(pred_probs, na.rm = TRUE)

  predictive_performance <- data.frame(
    species = species_names,
    observed_prevalence = observed_prev[species_names],
    predicted_prevalence = predicted_prev[species_names],
    auc = NA_real_,
    rmse = NA_real_,
    tjur_r2 = NA_real_,
    prediction_method = prediction_method,
    stringsAsFactors = FALSE
  )

  for (sp in species_names) {
    obs <- Y[, sp]
    score <- pred_probs[, sp]

    predictive_performance$auc[predictive_performance$species == sp] <- simple_auc(obs, score)
    predictive_performance$rmse[predictive_performance$species == sp] <-
      sqrt(mean((obs - score)^2, na.rm = TRUE))
    predictive_performance$tjur_r2[predictive_performance$species == sp] <-
      mean(score[obs == 1], na.rm = TRUE) - mean(score[obs == 0], na.rm = TRUE)
  }

  write.csv(predictive_performance, out_path("hmsc_species_predictive_performance.csv"), row.names = FALSE)

  # Five-fold cross-validation.
  # The folds are fixed by the seed so the validation split can be reproduced.
  set.seed(123)
  n_folds <- 5
  fold_id <- sample(rep(seq_len(n_folds), length.out = nrow(Y)))

  cv_partition <- data.frame(
    Sites = site_ids_model,
    fold = fold_id,
    stringsAsFactors = FALSE
  )
  write.csv(cv_partition, out_path("hmsc_cv_partition.csv"), row.names = FALSE)

  cv_pred_probs <- matrix(
    NA_real_,
    nrow = nrow(Y),
    ncol = ncol(Y),
    dimnames = list(as.character(site_ids_model), species_names)
  )

  cv_model_rows <- list()

  for (fold in seq_len(n_folds)) {
    train_rows <- fold_id != fold
    test_rows <- fold_id == fold

    Y_train <- Y[train_rows, , drop = FALSE]
    X_train <- XData[train_rows, , drop = FALSE]
    X_test <- XData[test_rows, , drop = FALSE]

    # Use the same model structure as the final model where possible.
    cv_studyDesign <- data.frame(site = factor(seq_len(nrow(Y_train))))
    cv_ranLevels <- list(site = HmscRandomLevel(units = cv_studyDesign$site))

    cv_model <- Hmsc(
      Y = Y_train,
      XData = X_train,
      XFormula = x_formula,
      distr = "probit",
      studyDesign = cv_studyDesign,
      ranLevels = cv_ranLevels
    )

    cv_fit <- fit_hmsc_with_settings(
      cv_model,
      paste("CV fold", fold, "site-level random effect"),
      mcmc_settings
    )

    # If the random-effect model fails in a fold, use the fixed-effect structure.
    if (!cv_fit$success) {
      cv_model_fixed <- Hmsc(
        Y = Y_train,
        XData = X_train,
        XFormula = x_formula,
        distr = "probit"
      )

      cv_fit <- fit_hmsc_with_settings(
        cv_model_fixed,
        paste("CV fold", fold, "fixed effects only"),
        mcmc_settings
      )
    }

    if (!cv_fit$success) {
      stop("Cross-validation failed in fold ", fold, ". Error: ", cv_fit$error)
    }

    cv_beta_samples <- extract_beta_samples(cv_fit$model)
    cv_beta_mean <- apply(cv_beta_samples, c(1, 2), mean, na.rm = TRUE)

    cv_cov_names <- cv_fit$model$covNames
    if (is.null(cv_cov_names) || length(cv_cov_names) != dim(cv_beta_mean)[1]) {
      cv_cov_names <- colnames(model.matrix(x_formula, data = X_train))
    }

    cv_species_names <- cv_fit$model$spNames
    if (is.null(cv_species_names) || length(cv_species_names) != dim(cv_beta_mean)[2]) {
      cv_species_names <- colnames(Y_train)
    }

    dimnames(cv_beta_mean) <- list(cv_cov_names, cv_species_names)

    X_test_matrix <- model.matrix(x_formula, data = X_test)
    beta_for_cv_prediction <- cv_beta_mean[match(colnames(X_test_matrix), rownames(cv_beta_mean)), , drop = FALSE]

    if (any(is.na(match(colnames(X_test_matrix), rownames(cv_beta_mean))))) {
      stop("Coefficient names did not match the test design matrix in CV fold ", fold)
    }

    fold_pred <- pnorm(X_test_matrix %*% beta_for_cv_prediction)
    fold_pred <- as.matrix(fold_pred)
    colnames(fold_pred) <- cv_species_names

    cv_pred_probs[test_rows, cv_species_names] <- fold_pred

    cv_model_rows[[fold]] <- data.frame(
      fold = fold,
      model_label = cv_fit$label,
      success = cv_fit$success,
      elapsed_minutes = cv_fit$elapsed_minutes,
      samples = cv_fit$settings$samples,
      transient = cv_fit$settings$transient,
      thin = cv_fit$settings$thin,
      nChains = cv_fit$settings$nChains,
      warnings = paste(cv_fit$warnings, collapse = " | "),
      error = cv_fit$error,
      stringsAsFactors = FALSE
    )
  }

  cv_model_run <- do.call(rbind, cv_model_rows)
  write.csv(cv_model_run, out_path("hmsc_cv_model_run.csv"), row.names = FALSE)

  write.csv(
    data.frame(Sites = site_ids_model, cv_pred_probs, check.names = FALSE),
    out_path("hmsc_cv_predicted_probabilities.csv"),
    row.names = FALSE
  )

  cv_predicted_prev <- colMeans(cv_pred_probs, na.rm = TRUE)

  cv_predictive_performance <- data.frame(
    species = species_names,
    observed_prevalence = observed_prev[species_names],
    cv_predicted_prevalence = cv_predicted_prev[species_names],
    cv_auc = NA_real_,
    cv_rmse = NA_real_,
    cv_tjur_r2 = NA_real_,
    stringsAsFactors = FALSE
  )

  for (sp in species_names) {
    obs <- Y[, sp]
    score <- cv_pred_probs[, sp]

    cv_predictive_performance$cv_auc[cv_predictive_performance$species == sp] <- simple_auc(obs, score)
    cv_predictive_performance$cv_rmse[cv_predictive_performance$species == sp] <-
      sqrt(mean((obs - score)^2, na.rm = TRUE))
    cv_predictive_performance$cv_tjur_r2[cv_predictive_performance$species == sp] <-
      mean(score[obs == 1], na.rm = TRUE) - mean(score[obs == 0], na.rm = TRUE)
  }

  write.csv(
    cv_predictive_performance,
    out_path("hmsc_cv_species_predictive_performance.csv"),
    row.names = FALSE
  )

  insample_vs_cv <- merge(
    predictive_performance,
    cv_predictive_performance,
    by = "species",
    all = TRUE
  )
  write.csv(insample_vs_cv, out_path("hmsc_insample_vs_cv_performance.csv"), row.names = FALSE)

  cv_summary <- data.frame(
    n_folds = n_folds,
    mean_cv_auc = mean(cv_predictive_performance$cv_auc, na.rm = TRUE),
    median_cv_auc = median(cv_predictive_performance$cv_auc, na.rm = TRUE),
    mean_in_sample_auc = mean(predictive_performance$auc, na.rm = TRUE),
    median_in_sample_auc = median(predictive_performance$auc, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  write.csv(cv_summary, out_path("hmsc_cv_summary.csv"), row.names = FALSE)

  # Figure 4: observed prevalence against cross-validated predicted prevalence.
  png(out_path("hmsc_cv_observed_vs_predicted_prevalence.png"), width = 1300, height = 1100, res = 160)
  old_par <- par(no.readonly = TRUE)
  par(mar = c(5, 5, 4, 2))

  plot(
    cv_predictive_performance$observed_prevalence,
    cv_predictive_performance$cv_predicted_prevalence,
    pch = 19,
    col = "#2C7FB8",
    xlab = "Observed prevalence",
    ylab = "Cross-validated predicted prevalence",
    xlim = c(0, 1),
    ylim = c(0, 1)
  )
  grid(col = "grey85", lty = "dotted")
  abline(0, 1, lty = 2, col = "grey50")

  legend(
    "topleft",
    legend = c(
      "5-fold CV",
      paste0("Median AUC = ", format(round(cv_summary$median_cv_auc, 3), nsmall = 3)),
      paste0("Mean AUC = ", format(round(cv_summary$mean_cv_auc, 3), nsmall = 3))
    ),
    bty = "o",
    bg = "white",
    cex = 0.85
  )

  label_pos <- data.frame(
    species = c("sp10", "sp3", "sp6", "sp11", "sp22", "sp7", "sp8", "sp9",
                "sp19", "sp18", "sp20", "sp12", "sp17", "sp4", "sp15",
                "sp5", "sp1", "sp21", "sp14", "sp16", "sp13", "sp2"),
    label_x = c(0.035, 0.035, 0.035, 0.035, 0.035, 0.205, 0.205, 0.205,
                0.285, 0.315, 0.360, 0.365, 0.400, 0.525, 0.590,
                0.610, 0.705, 0.710, 0.755, 0.790, 0.875, 0.305),
    label_y = c(0.160, 0.132, 0.106, 0.080, 0.055, 0.198, 0.145, 0.105,
                0.272, 0.292, 0.314, 0.342, 0.392, 0.515, 0.575,
                0.570, 0.675, 0.695, 0.705, 0.748, 0.835, 0.282),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(cv_predictive_performance))) {
    sp <- cv_predictive_performance$species[i]
    px <- cv_predictive_performance$observed_prevalence[i]
    py <- cv_predictive_performance$cv_predicted_prevalence[i]

    label_row <- label_pos[label_pos$species == sp, , drop = FALSE]
    if (nrow(label_row) == 1) {
      tx <- label_row$label_x
      ty <- label_row$label_y
    } else {
      tx <- min(px + 0.025, 0.97)
      ty <- min(py + 0.015, 0.97)
    }

    if (sqrt((tx - px)^2 + (ty - py)^2) > 0.035) {
      segments(px, py, tx, ty, col = "grey45", lwd = 0.7)
    }

    text(tx, ty, labels = sp, cex = 0.75, pos = ifelse(tx < px, 2, 4))
  }

  par(old_par)
  dev.off()

  # Figure 5: in-sample AUC against cross-validated AUC.
  png(out_path("hmsc_insample_vs_cv_auc.png"), width = 1400, height = 1100, res = 160)
  old_par <- par(no.readonly = TRUE)
  par(mar = c(5, 5, 4, 2))

  plot(
    insample_vs_cv$auc,
    insample_vs_cv$cv_auc,
    pch = 19,
    col = "#4C78A8",
    xlab = "In-sample AUC",
    ylab = "Cross-validated AUC",
    xlim = c(0, 1),
    ylim = c(0, 1),
    main = "In-sample vs cross-validated AUC"
  )
  abline(0, 1, lty = 2, col = "grey55")
  abline(h = 0.5, lty = 3, col = "grey50")
  abline(v = 0.5, lty = 3, col = "grey50")

  text(
    insample_vs_cv$auc,
    insample_vs_cv$cv_auc,
    labels = insample_vs_cv$species,
    pos = 4,
    cex = 0.75
  )

  par(old_par)
  dev.off()

  # 9. Extract residual species associations

  association_attempt <- try(computeAssociations(selected_model), silent = TRUE)

  if (!inherits(association_attempt, "try-error") &&
      length(association_attempt) > 0 &&
      !is.null(association_attempt[[1]]$mean)) {
    association_matrix <- association_attempt[[1]]$mean
    association_source <- "Hmsc computeAssociations"
  } else {
    residual_matrix <- Y - pred_probs
    association_matrix <- cor(residual_matrix, use = "pairwise.complete.obs")
    association_source <- "correlation of observed minus predicted residuals"
  }

  association_matrix <- as.matrix(association_matrix)
  rownames(association_matrix) <- species_names
  colnames(association_matrix) <- species_names

  write.csv(association_matrix, out_path("hmsc_species_association_matrix.csv"), row.names = TRUE)

  make_heatmap(
    association_matrix,
    "hmsc_species_association_heatmap.png",
    "Hmsc residual species association matrix",
    zlim = c(-1, 1),
    value_labels = FALSE,
    legend_title = "Association"
  )

  assoc_values <- if (nrow(association_matrix) >= 2) {
    association_matrix[upper.tri(association_matrix)]
  } else {
    numeric(0)
  }

  association_summary <- data.frame(
    association_source = association_source,
    median_absolute_association = ifelse(length(assoc_values) == 0, NA_real_, median(abs(assoc_values), na.rm = TRUE)),
    maximum_absolute_association = ifelse(length(assoc_values) == 0, NA_real_, max(abs(assoc_values), na.rm = TRUE)),
    strength_label = association_strength_label(assoc_values),
    stringsAsFactors = FALSE
  )

  write.csv(association_summary, out_path("hmsc_species_association_summary.csv"), row.names = FALSE)

  # 10. Check MCMC diagnostics

  coda_attempt <- try(convertToCodaObject(selected_model), silent = TRUE)

  diagnostics_ok <- FALSE
  diagnostics_available <- FALSE

  diagnostics_summary <- data.frame(
    selected_model = selected_model_name,
    selected_model_file = selected_model_file,
    beta_parameter_count = NA_integer_,
    min_effective_sample_size = NA_real_,
    median_effective_sample_size = NA_real_,
    max_rhat_point_estimate = NA_real_,
    median_rhat_point_estimate = NA_real_,
    max_rhat_upper_ci = NA_real_,
    diagnostics_ok = FALSE,
    stringsAsFactors = FALSE
  )

  if (!inherits(coda_attempt, "try-error") &&
      "Beta" %in% names(coda_attempt) &&
      requireNamespace("coda", quietly = TRUE)) {
    beta_coda <- coda_attempt$Beta
    ess <- coda::effectiveSize(beta_coda)
    gelman_attempt <- try(coda::gelman.diag(beta_coda, multivariate = FALSE), silent = TRUE)

    if (!inherits(gelman_attempt, "try-error")) {
      rhat <- gelman_attempt$psrf[, "Point est."]
      rhat_upper <- gelman_attempt$psrf[, "Upper C.I."]

      diagnostics_available <- TRUE

      diagnostics_ok <- max(rhat, na.rm = TRUE) <= 1.10 &&
        median(rhat, na.rm = TRUE) <= 1.05 &&
        min(ess, na.rm = TRUE) >= 100

      diagnostics_table <- data.frame(
        parameter = names(ess),
        effective_sample_size = as.numeric(ess),
        rhat_point_estimate = as.numeric(rhat[names(ess)]),
        rhat_upper_ci = as.numeric(rhat_upper[names(ess)]),
        stringsAsFactors = FALSE
      )

      write.csv(diagnostics_table, out_path("hmsc_beta_diagnostics_table.csv"), row.names = FALSE)

      diagnostics_summary <- data.frame(
        selected_model = selected_model_name,
        selected_model_file = selected_model_file,
        beta_parameter_count = length(ess),
        min_effective_sample_size = min(ess, na.rm = TRUE),
        median_effective_sample_size = median(ess, na.rm = TRUE),
        max_rhat_point_estimate = max(rhat, na.rm = TRUE),
        median_rhat_point_estimate = median(rhat, na.rm = TRUE),
        max_rhat_upper_ci = max(rhat_upper, na.rm = TRUE),
        diagnostics_ok = diagnostics_ok,
        stringsAsFactors = FALSE
      )
    }
  }

  if (!diagnostics_available) {
    diagnostics_table <- data.frame(
      parameter = character(),
      effective_sample_size = numeric(),
      rhat_point_estimate = numeric(),
      rhat_upper_ci = numeric(),
      stringsAsFactors = FALSE
    )
    write.csv(diagnostics_table, out_path("hmsc_beta_diagnostics_table.csv"), row.names = FALSE)
  }

  write.csv(diagnostics_summary, out_path("hmsc_diagnostics_summary.csv"), row.names = FALSE)

  # 11. Save a compact list of created output files

  output_files <- sort(list.files(output_dir, recursive = TRUE, full.names = TRUE))
  output_files <- unique(c(output_files, out_path("created_files.csv")))

  created_files_table <- data.frame(
    file_path = output_files,
    file_name = basename(output_files),
    stringsAsFactors = FALSE
  )

  write.csv(created_files_table, out_path("created_files.csv"), row.names = FALSE)
}

tryCatch(
  {
    main()
    message("Hmsc JSDM analysis completed successfully.")
  },
  error = function(e) {
    project_dir <- "D:/Beetles"
    output_dir <- file.path(project_dir, "Hmsc_JSDM_outputs")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    log_file <- file.path(output_dir, "hmsc_error_log.txt")

    writeLines(
      c(
        "Hmsc JSDM analysis failed.",
        paste("Time:", as.character(Sys.time())),
        paste("Error:", conditionMessage(e)),
        "",
        "Traceback:",
        paste(capture.output(traceback()), collapse = "\n")
      ),
      con = log_file,
      useBytes = TRUE
    )

    message("Hmsc JSDM analysis failed. See hmsc_error_log.txt in the output folder.")
    stop("Analysis stopped because an error occurred. Check hmsc_error_log.txt for details.")
  }
)
