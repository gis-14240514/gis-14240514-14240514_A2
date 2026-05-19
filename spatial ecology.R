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
  "fig_01_hmsc_beta_heatmap.png",
  "Hmsc posterior mean environmental effects",
  value_labels = TRUE,
  legend_title = "Posterior mean"
)

# Classify effects with strong posterior support.
supported_mat <- matrix(
  "not strongly supported",
  nrow = length(species_names),
  ncol = length(predictors),
  dimnames = list(species_names, predictors)
)

for (i in seq_len(nrow(beta_long))) {
  if (isTRUE(beta_long$ci_excludes_zero[i]) &&
      beta_long$posterior_mean[i] > 0 &&
      beta_long$posterior_prob_positive[i] >= 0.95) {
    supported_mat[beta_long$species[i], beta_long$predictor[i]] <- "supported positive"
  }
  
  if (isTRUE(beta_long$ci_excludes_zero[i]) &&
      beta_long$posterior_mean[i] < 0 &&
      beta_long$posterior_prob_negative[i] >= 0.95) {
    supported_mat[beta_long$species[i], beta_long$predictor[i]] <- "supported negative"
  }
}

make_supported_heatmap(supported_mat, "fig_02_hmsc_supported_effects_heatmap.png")


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

png(out_path("fig_03_hmsc_predictor_response_counts.png"), width = 1500, height = 1000, res = 160)
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

png(out_path("fig_04_hmsc_coefficient_distribution_by_predictor.png"), width = 1500, height = 1000, res = 160)
old_par <- par(no.readonly = TRUE)
par(mar = c(7, 5, 4, 2))

coef_split <- split(beta_long$posterior_mean, beta_long$predictor)

boxplot(
  coef_split[predictors],
  las = 2,
  col = "#A6CEE3",
  ylab = "Posterior mean coefficient",
  main = "Hmsc coefficient distributions by predictor"
)

abline(h = 0, lty = 2, col = "grey40")

par(old_par)
dev.off()
}
