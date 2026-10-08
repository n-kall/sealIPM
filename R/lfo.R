extract_future_log_lik <- function(forecast) {
    # Stan sums the conditional likelihood over the entire forecast horizon.
    # Averaging over draws integrates over parameter and future-state uncertainty.
    forecast$draws(
        variables = "log_lik",
        format = "draws_matrix"
    )
}

log_mean_exp_cols <- function(log_lik_matrix) {
    if (ncol(log_lik_matrix) == 0L) {
        return(numeric())
    }

    matrixStats::colLogSumExps(log_lik_matrix) -
        log(nrow(log_lik_matrix))
}

elpd_lfo <- function(forecast) {
    log_lik <- extract_future_log_lik(forecast) |>
        log_mean_exp_cols()

    return(log_lik)
}

# Observation likelihoods and the matching state-year indices. Keep complete
# composition records together: the likelihood is multinomial, not cellwise.
lfo_observation_fields <- function() {
    list(
        aerial_count = c("n_aerial_years", "aerial_year", "obs_aerial_count"),
        hunting_bags_finland = c("n_hunting_bag_years_finland", "hunting_bag_year_finland", "obs_hunting_bag_finland"),
        hunting_bags_sweden = c("n_hunting_bag_years_sweden", "hunting_bag_year_sweden", "obs_hunting_bag_sweden"),
        hunting_comp_finland = c("n_hunting_comp_years_finland", "hunting_comp_year_finland", "obs_hunting_comp_finland"),
        hunting_comp_sweden = c("n_hunting_comp_years_sweden", "hunting_comp_year_sweden", "obs_hunting_comp_sweden"),
        bycatch_comp = c("n_bycatch_years", "bycatch_comp_year", "obs_bycatch_comp"),
        pregnancy_count = c("n_pregnancy_years", "pregnancy_count_year", "obs_pregnancy_count", "pregnancy_sample_size"),
        reproductive_signs = c("n_reproductive_years", "reproductive_signs_year", "obs_reproductive_signs_finland")
    )
}

lfo_training_data <- function(stan_data, cutoff) {
    for (fields in lfo_observation_fields()) {
        keep <- stan_data[[fields[2L]]] <= cutoff
        stan_data[[fields[1L]]] <- as.integer(sum(keep))
        for (field in fields[-1L]) {
            value <- stan_data[[field]]
            stan_data[[field]] <- if (is.matrix(value)) {
                value[keep, , drop = FALSE]
            } else {
                value[keep]
            }
        }
    }
    stan_data
}


# Preserve the hidden initialization year used by build_grey_stan_data().
lfo_prefix_data <- function(stan_data, n_state_years) {
    out <- lfo_training_data(stan_data, n_state_years)
    out$n_state_years <- as.integer(n_state_years)
    for (field in c("hunting_quota_sweden", "hunting_quota_finland")) {
        out[[field]] <- out[[field]][seq_len(n_state_years)]
    }
    for (field in c("herring_index_baltic_proper_gulf_finland",
                    "herring_index_gulf_bothnia")) {
        out[[field]] <- out[[field]][seq_len(n_state_years + 1L)]
    }
    out
}

get_species_lfo_model <- function(species, ...) {
    stan_file <- system.file(
        "bin",
        "stan",
        paste0(species, "_seal_IPM_lfo.stan"),
        package = "sealIPM",
        mustWork = TRUE
    )

    instantiate::stan_package_model(
        name = paste0(species, "_seal_IPM_lfo"),
        package = "sealIPM",
        stan_file = stan_file
    )
}

lfo_score_proposal <- function(fit, stan_data, gq_args) {
    species <- fit$species
    lfo_model <- get_species_lfo_model(species)

    future_years <- seq.int(max(fit$years) + 1L,
        min(fit$years) + stan_data$n_state_years - 2L)
    future_data <- build_future_stan_data(
        species = species,
        actual_data = fit$data,
        future_years = future_years,
        fitted_stan_data = fit$stan_data
    )

    lfo_variables <- names(
        lfo_model$variables()$parameters
    )
    past_draws <- fit$fit$draws(
        variables = lfo_variables,
        format = "draws_matrix"
    )

    scored <- do.call(lfo_model$generate_quantities, c(
        list(
            fitted_params = past_draws,
            data = utils::modifyList(fit$stan_data, future_data)
        ),
        gq_args
    ))
    if (any(scored$return_codes() != 0L)) {
        stop("LFO likelihood generation failed.", call. = FALSE)
    }
    future_lik <- scored$draws(
        variables = "log_lik_future_by_year",
        format = "draws_matrix"
    )
    expected_names <- sprintf("log_lik_future_by_year[%d]",
                               seq_along(future_years))
    if (!all(expected_names %in% colnames(future_lik)) ||
        nrow(future_lik) != nrow(past_draws)) {
        stop("LFO likelihoods do not match the proposal draws and future years.", call. = FALSE)
    }
    future_lik <- future_lik[, expected_names, drop = FALSE]
    if (anyNA(future_lik) || any(future_lik == Inf)) {
        stop("LFO likelihood generation produced invalid values.", call. = FALSE)
    }
    annual <- matrix(0, nrow(future_lik), stan_data$n_state_years)
    annual[, seq.int(fit$stan_data$n_state_years + 1L,
                      stan_data$n_state_years)] <- future_lik
    list(fit = fit, scored = scored,
         chain_id = fit$fit$draws(format = "draws_df")$.chain,
         log_lik = annual)
}

lfo_refit <- function(fit, stan_data, training_end, refit_args) {
    n_state <- match(training_end, seq.int(min(fit$years) - 1L,
                                          min(fit$years) + stan_data$n_state_years - 2L))
    training_data <- lfo_prefix_data(stan_data, n_state)
    model <- get_species_model(fit$species)
    if (is.null(refit_args$init)) {
        refit_args$init <- function() grey_init_fun(n_state,
            2L * training_data$n_age_classes, fit$priors)
    }
    sampled <- do.call(model$sample, c(list(data = training_data), refit_args))
    if (any(sampled$return_codes() != 0L)) {
        stop("LFO refit failed at year ", training_end, call. = FALSE)
    }
    fit$fit <- sampled
    fit$years <- seq.int(min(fit$years), training_end)
    fit$stan_data <- training_data
    fit
}

lfo_psis <- function(log_ratios, chain_id) {
    S <- length(log_ratios)
    if (S < 10L) stop("LFO requires at least 10 proposal draws.", call. = FALSE)
    if (anyNA(log_ratios) || any(log_ratios == Inf) ||
        all(log_ratios == -Inf)) {
        stop("LFO ratios are undefined or all proposal weights are zero.", call. = FALSE)
    }
    constant <- all(log_ratios == log_ratios[1L])
    if (constant || any(log_ratios == -Inf)) {
        # Zero weights give infinite reciprocal ratios, preventing an
        # autocorrelation estimate. Retain weights but flag unavailable
        # efficiency; ESS then uses loo's iid default.
        r_eff <- 1
    } else {
        inverse <- exp(-log_ratios - max(-log_ratios))
        r_eff <- loo::relative_eff(matrix(inverse, ncol = 1L), chain_id = chain_id)
        if (any(!is.finite(r_eff)) || any(r_eff <= 0)) {
            stop("Cannot estimate relative efficiency of LFO ratios.", call. = FALSE)
        }
    }
    ps <- suppressWarnings(loo::psis(log_ratios, r_eff = r_eff))
    if (constant) {
        ps$diagnostics$pareto_k[] <- 0
        ps$diagnostics$n_eff[] <- S
    }
    attr(ps, "efficiency_estimated") <- constant || all(is.finite(log_ratios))
    ps
}

#' Forward PSIS leave-future-out log weights
#'
#' Starts from the supplied initial training fit, predicts the next year,
#' then incorporates one additional observed year at a time. Raw log ratios
#' are cumulative sums of annual observation log likelihoods since the last
#' exact fit (positive sign). Refits occur only when Pareto-k exceeds the
#' threshold. Future state innovations are simulated from the fitted continuous
#' model's priors, not from the separate discrete forecast simulator.
#'
#' @param fit An initial grey-seal MCMC posterior fit, e.g. the first five years.
#' @param future_years Consecutive evaluation years beginning immediately after
#'   the final fitted year. With horizon = 1, each is predicted using observations
#'   through the preceding year.
#' @param actual_data Complete raw observations and covariates, including the
#'   evaluation years. Defaults to fit$data.
#' @param horizon Number of future years scored jointly at each training cutoff.
#' @param k_threshold Refit when Pareto-k exceeds this value (default 0.7).
#' @param refit Whether to refit when PSIS fails its threshold. If FALSE, retain
#'   the approximate weights and scores and flag them as unreliable.
#' @param refit_args Named CmdStan sample() controls for necessary refits.
#' @param gq_args Named CmdStan generate_quantities() controls.
#' @return A list with normalized log_weights and psis_objects, each indexed by
#'   the first prediction year; pointwise diagnostics and predictive scores;
#'   proposals containing the exact fit, future continuous scoring paths, and
#'   annual log likelihoods; log_predictive_lik vectors for each prediction block;
#'   and refit_years. `pointwise$proposal` identifies which proposal's scored draws
#'   correspond to each weight vector. Pass `psis_objects[[j]]` to bayesplot with
#'   a one-column prediction matrix from that proposal's scored draws.
#'   Different proposals can have different draw counts and must not be mixed.
#'
#' @details Validation is conditional on supplied covariates. The prefix of the
#'   full data must match the original training data, including covariate scaling
#'   and the hidden initialization year. The initial prediction and predictions
#'   immediately after refits use equal weights. The year being predicted is
#'   never included in its importance ratios. High-k scores are retained.
#' @export
lfo_log_weights <- function(fit, future_years, actual_data = fit$data,
                            horizon = 1L, k_threshold = 0.7, refit = TRUE,
                            refit_args = list(), gq_args = list()) {
    if (!inherits(fit, "seal_ipm") || !identical(fit$species, "grey") ||
        !identical(fit$method, "sample") || !isTRUE(fit$stan_data$prior_only == 0)) {
        stop("`fit` must be an initial grey-seal MCMC posterior fit.", call. = FALSE)
    }
    years <- validate_state_years(fit$years, "fit$years")
    future_years <- validate_state_years(future_years, "future_years")
    if (future_years[1L] != max(years) + 1L) {
        stop("`future_years` must start immediately after the final fitted year.", call. = FALSE)
    }
    if (is.null(actual_data)) stop("Supply complete `actual_data`.", call. = FALSE)
    if (!is.numeric(horizon) || length(horizon) != 1L ||
        !is.finite(horizon) || horizon < 1 || horizon != floor(horizon) ||
        horizon > length(future_years)) {
        stop("`horizon` must be a positive integer within future_years.", call. = FALSE)
    }
    if (!is.numeric(k_threshold) || length(k_threshold) != 1L ||
        !is.finite(k_threshold) || k_threshold <= 0 || k_threshold > 0.7) {
        stop("`k_threshold` must be in (0, 0.7].", call. = FALSE)
    }
    if (!is.logical(refit) || length(refit) != 1L || is.na(refit)) {
        stop("`refit` must be TRUE or FALSE.", call. = FALSE)
    }
    for (args in list(refit_args, gq_args)) {
        if (!is.list(args) || (length(args) &&
            (is.null(names(args)) || anyNA(names(args)) ||
             any(!nzchar(names(args))) || anyDuplicated(names(args)))) ||
            any(c("data", "fitted_params", "fixed_param") %in% names(args))) {
            stop("Supply named sampling/GQ controls, excluding data, fitted_params and fixed_param.", call. = FALSE)
        }
    }
    full_years <- c(years, future_years)
    stan_data <- build_stan_data(fit$species, actual_data, full_years,
                                 fit$priors, prior_only = FALSE)
    prefix <- lfo_prefix_data(stan_data, fit$stan_data$n_state_years)
    for (field in names(fit$stan_data)) {
        if (!isTRUE(all.equal(prefix[[field]], fit$stan_data[[field]],
                              check.attributes = FALSE))) {
            stop("Full data change the fitted training prefix: ", field,
                 ". Use the same observations, priors and covariate scaling.", call. = FALSE)
        }
    }
    # Keep complete raw data available to the same future-data builder used
    # by forecast_ipm(), including after an adaptive refit.
    fit$data <- actual_data
    prediction_years <- future_years[seq_len(length(future_years) - horizon + 1L)]
    training_ends <- prediction_years - 1L
    proposals <- list(lfo_score_proposal(fit, stan_data, gq_args))
    current <- proposals[[1L]]
    cutoff <- fit$stan_data$n_state_years
    weights <- psis_objects <- log_predictive_lik <- vector("list", length(prediction_years))
    names(weights) <- names(psis_objects) <- names(log_predictive_lik) <- as.character(prediction_years)
    pointwise <- data.frame(training_end = training_ends,
        prediction_year = prediction_years, prediction_end = prediction_years + horizon - 1L,
        proposal = NA_integer_, fit_end = NA_integer_, pareto_k = NA_real_,
        attempted_pareto_k = NA_real_, psis_ess = NA_real_,
        reliable = FALSE, refitted = FALSE, elpd = NA_real_)
    state_years <- seq.int(min(years) - 1L, max(future_years))
    for (j in seq_along(prediction_years)) {
        index <- match(training_ends[j], state_years)
        lr <- if (index == cutoff) rep(0, nrow(current$log_lik)) else
            rowSums(current$log_lik[, seq.int(cutoff + 1L, index), drop = FALSE])
        if (all(lr == -Inf)) {
            ps <- NULL
            k <- Inf
        } else {
            ps <- lfo_psis(lr, current$chain_id)
            k <- as.numeric(ps$diagnostics$pareto_k)
        }
        pointwise$attempted_pareto_k[j] <- k
        reliable <- !is.na(k) && k <= k_threshold
        if (!reliable && refit) {
            new_fit <- lfo_refit(fit, stan_data, training_ends[j], refit_args)
            current <- lfo_score_proposal(new_fit, stan_data, gq_args)
            proposals[[length(proposals) + 1L]] <- current
            cutoff <- index
            ps <- lfo_psis(rep(0, nrow(current$log_lik)), current$chain_id)
            k <- 0
            reliable <- TRUE
            pointwise$refitted[j] <- TRUE
        }
        if (is.null(ps)) {
            stop("All importance weights are zero at training year ",
                 training_ends[j], "; enable refits.", call. = FALSE)
        }
        lw <- as.numeric(stats::weights(ps, log = TRUE, normalize = TRUE))
        lp <- rowSums(current$log_lik[, seq.int(index + 1L, index + horizon), drop = FALSE])
        weights[[j]] <- lw
        psis_objects[[j]] <- ps
        log_predictive_lik[[j]] <- lp
        pointwise$proposal[j] <- length(proposals)
        pointwise$fit_end[j] <- max(current$fit$years)
        pointwise$pareto_k[j] <- k
        pointwise$psis_ess[j] <- as.numeric(ps$diagnostics$n_eff)
        pointwise$reliable[j] <- reliable
        pointwise$elpd[j] <- matrixStats::logSumExp(lw + lp)
    }
    list(log_weights = weights, psis_objects = psis_objects,
         pointwise = pointwise, proposals = proposals,
         log_predictive_lik = log_predictive_lik,
         elpd_lfo = sum(pointwise$elpd),
         refit_years = pointwise$training_end[pointwise$refitted])
}

#' @rdname lfo_log_weights
#' @export
psis_lfo <- function(fit, future_years, actual_data = fit$data,
                     horizon = 1L, k_threshold = 0.7, refit = TRUE,
                     refit_args = list(), gq_args = list()) {
    lfo_log_weights(fit, future_years, actual_data, horizon,
                    k_threshold, refit, refit_args, gq_args)
}
