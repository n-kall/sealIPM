# Independent continuous/discrete forecasts and comparison utilities

library(cmdstanr)
library(ggplot2)
library(patchwork)

# -----------------------------------------------------------------------------
# Compile model and expose Stan functions
# -----------------------------------------------------------------------------

stan_file <- "src/stan/grey_seal_IPM_forecast.stan"

if (!exists("m", inherits = FALSE)) {
    m <- cmdstan_model(
        stan_file,
        force_recompile = TRUE
    )
    m$expose_functions()
}

required_functions <- c(
    "create_aging_matrix",
    "update_population_from_survivors",
    "run_state_process_from_first_population",
    "run_state_process_from_survivors_rng"
)

missing_functions <- required_functions[
    !vapply(
        required_functions,
        function(x) is.function(m$functions[[x]]),
        logical(1)
    )
]

if (length(missing_functions) > 0L) {
    stop(
        "The following exposed Stan functions are unavailable: ",
        paste(missing_functions, collapse = ", "),
        ". Remove `m`, recompile the updated Stan model, and source this file again.",
        call. = FALSE
    )
}

.default_survivors_final <- c(
    180,
    240,
    310,
    400,
    500,
    3000,
    175,
    235,
    300,
    390,
    480,
    2800
)

.default_S_diag <- c(
    0.82,
    0.90,
    0.94,
    0.96,
    0.97,
    0.97,
    0.80,
    0.88,
    0.92,
    0.94,
    0.95,
    0.95
)

.default_fate_probability_tolerance <- 1e-6

.state_output_names <- c(
    "birth_rate",
    "pregnancy_rate",
    "population_total",
    "population_comp",
    "survivors",
    "deaths_or_bycatch",
    "hunted_sweden",
    "hunted_finland",
    "bycatch_expected",
    "hunting_bag_total_sweden",
    "hunting_bag_total_finland",
    "hunted_total",
    "reproductive_probs"
)

name_state_process_output <- function(x) {
    if (length(x) != length(.state_output_names)) {
        stop(
            "Unexpected state-process tuple length: ",
            length(x),
            ". Expected ",
            length(.state_output_names),
            ".",
            call. = FALSE
        )
    }
    names(x) <- .state_output_names
    x
}

validate_forecast_inputs <- function(
    n_state_years,
    n_age_classes,
    survivors_final,
    baseline_birth_rate,
    S_diag,
    hunting_selectivity_sweden,
    hunting_selectivity_finland,
    hunting_quota_sweden,
    hunting_quota_finland,
    epsilon_h_sw,
    epsilon_h_fi,
    pi_s,
    pi_c,
    t_mate_to_preg,
    t_birth_to_end_hunt,
    rel_tol,
    abs_tol,
    max_num_steps,
    fate_probability_tolerance
) {
    n_demo_groups <- 2L * n_age_classes

    stopifnot(
        n_state_years >= 1L,
        n_age_classes >= 2L,
        length(survivors_final) == n_demo_groups,
        length(S_diag) == n_demo_groups,
        length(hunting_selectivity_sweden) == n_demo_groups,
        length(hunting_selectivity_finland) == n_demo_groups,
        length(baseline_birth_rate) == n_state_years + 1L,
        length(hunting_quota_sweden) == n_state_years,
        length(hunting_quota_finland) == n_state_years,
        length(epsilon_h_sw) == n_state_years,
        length(epsilon_h_fi) == n_state_years,
        length(pi_s) == n_state_years,
        length(pi_c) == n_state_years,
        all(is.finite(survivors_final)),
        all(survivors_final >= 0),
        all(is.finite(S_diag)),
        all(S_diag > 0 & S_diag <= 1),
        all(baseline_birth_rate >= 0 & baseline_birth_rate <= 1),
        all(hunting_quota_sweden >= 0),
        all(hunting_quota_finland >= 0),
        all(pi_s >= 0 & pi_s <= 1),
        all(pi_c >= 0 & pi_c <= 1),
        t_mate_to_preg >= 0 && t_mate_to_preg <= 1,
        t_birth_to_end_hunt >= 0 && t_birth_to_end_hunt <= 1,
        is.finite(rel_tol) && rel_tol > 0,
        is.finite(abs_tol) && abs_tol > 0,
        max_num_steps >= 1L,
        is.finite(fate_probability_tolerance),
        fate_probability_tolerance > 0
    )

    invisible(n_demo_groups)
}

# -----------------------------------------------------------------------------
# Run one continuous forecast independently
# -----------------------------------------------------------------------------

run_continuous_forecast <- function(
    seed = 123L,
    n_state_years = 10L,
    n_age_classes = 6L,
    survivors_final = .default_survivors_final,
    birth_rate_first = 0.82,
    baseline_birth_rate = rep(0.86, n_state_years + 1L),
    density_dependence_intercept = 0.05,
    density_dependence_slope = 1e-5,
    S_diag = .default_S_diag,
    hunting_selectivity_sweden = rep(0, 2L * n_age_classes),
    hunting_selectivity_finland = rep(0, 2L * n_age_classes),
    hunting_quota_sweden = rep(1000L, n_state_years),
    hunting_quota_finland = rep(1000L, n_state_years),
    hunting_effort_sd_sweden = 0.2,
    hunting_effort_sd_finland = 0.2,
    epsilon_h_sw = NULL,
    epsilon_h_fi = NULL,
    t_mate_to_preg = 0.5,
    t_birth_to_end_hunt = 0.75,
    pi_s = rep(0.8, n_state_years),
    pi_c = rep(0.7, n_state_years),
    prob_of_ca = 0.2,
    rel_tol = 1e-8,
    abs_tol = 1e-8,
    max_num_steps = 100000L,
    fate_probability_tolerance = .default_fate_probability_tolerance
) {
    n_state_years <- as.integer(n_state_years)
    n_age_classes <- as.integer(n_age_classes)
    n_demo_groups <- 2L * n_age_classes
    set.seed(as.integer(seed))

    if (is.null(epsilon_h_sw)) {
        epsilon_h_sw <- abs(rnorm(n_state_years))
    }
    if (is.null(epsilon_h_fi)) {
        epsilon_h_fi <- abs(rnorm(n_state_years))
    }

    validate_forecast_inputs(
        n_state_years,
        n_age_classes,
        survivors_final,
        baseline_birth_rate,
        S_diag,
        hunting_selectivity_sweden,
        hunting_selectivity_finland,
        hunting_quota_sweden,
        hunting_quota_finland,
        epsilon_h_sw,
        epsilon_h_fi,
        pi_s,
        pi_c,
        t_mate_to_preg,
        t_birth_to_end_hunt,
        rel_tol,
        abs_tol,
        max_num_steps,
        fate_probability_tolerance
    )

    mu_m <- -log(S_diag)
    aging_matrix <- m$functions$create_aging_matrix(
        n_demo_groups,
        n_age_classes
    )

    # Construct the first population from the same survivor boundary used by
    # the discrete process. This preserves a fair year-1 comparison.
    population_comp_first <- m$functions$update_population_from_survivors(
        survivors_final,
        aging_matrix,
        birth_rate_first,
        rnorm(1L),
        rnorm(1L),
        n_age_classes
    )
    population_total_first <- sum(population_comp_first)

    # Year 1 is supplied explicitly, so the first birth and sex innovations are
    # placeholders. Elements 2:n are used for later transitions.
    epsilon_birth <- c(0, rnorm(n_state_years - 1L))
    epsilon_sex <- c(0, rnorm(n_state_years - 1L))
    transition_noise_raw <- matrix(
        rnorm(3L * n_demo_groups * n_state_years),
        nrow = 3L * n_demo_groups,
        ncol = n_state_years
    )

    result <- m$functions$run_state_process_from_first_population(
        n_state_years,
        n_age_classes,
        population_comp_first,
        birth_rate_first,
        population_total_first,
        baseline_birth_rate,
        density_dependence_intercept,
        density_dependence_slope,
        aging_matrix,
        mu_m,
        hunting_selectivity_sweden,
        hunting_selectivity_finland,
        as.integer(hunting_quota_sweden),
        as.integer(hunting_quota_finland),
        hunting_effort_sd_sweden,
        hunting_effort_sd_finland,
        epsilon_h_sw,
        epsilon_h_fi,
        t_mate_to_preg,
        t_birth_to_end_hunt,
        epsilon_birth,
        epsilon_sex,
        transition_noise_raw,
        pi_s,
        pi_c,
        prob_of_ca,
        rel_tol,
        abs_tol,
        as.integer(max_num_steps),
        fate_probability_tolerance
    )

    name_state_process_output(result)
}

# -----------------------------------------------------------------------------
# Run one discrete forecast independently
# -----------------------------------------------------------------------------

run_discrete_forecast <- function(
    seed = NULL,
    n_state_years = 10L,
    n_age_classes = 6L,
    survivors_final = .default_survivors_final,
    birth_rate_first = 0.82,
    baseline_birth_rate = rep(0.86, n_state_years + 1L),
    density_dependence_intercept = 0.05,
    density_dependence_slope = 1e-5,
    S_diag = .default_S_diag,
    hunting_selectivity_sweden = rep(0, 2L * n_age_classes),
    hunting_selectivity_finland = rep(0, 2L * n_age_classes),
    hunting_quota_sweden = rep(1000L, n_state_years),
    hunting_quota_finland = rep(1000L, n_state_years),
    hunting_effort_sd_sweden = 0.2,
    hunting_effort_sd_finland = 0.2,
    epsilon_h_sw = NULL,
    epsilon_h_fi = NULL,
    t_mate_to_preg = 0.5,
    t_birth_to_end_hunt = 0.75,
    pi_s = rep(0.8, n_state_years),
    pi_c = rep(0.7, n_state_years),
    prob_of_ca = 0.2,
    rel_tol = 1e-8,
    abs_tol = 1e-8,
    max_num_steps = 100000L,
    fate_probability_tolerance = .default_fate_probability_tolerance
) {
    n_state_years <- as.integer(n_state_years)
    n_age_classes <- as.integer(n_age_classes)

    # set.seed() controls only R-generated effort effects. CmdStanR's exposed
    # function interface does not necessarily expose a public RNG-reset helper.
    if (!is.null(seed)) {
        set.seed(as.integer(seed))
    }
    if (is.null(epsilon_h_sw)) {
        epsilon_h_sw <- abs(rnorm(n_state_years))
    }
    if (is.null(epsilon_h_fi)) {
        epsilon_h_fi <- abs(rnorm(n_state_years))
    }

    validate_forecast_inputs(
        n_state_years,
        n_age_classes,
        survivors_final,
        baseline_birth_rate,
        S_diag,
        hunting_selectivity_sweden,
        hunting_selectivity_finland,
        hunting_quota_sweden,
        hunting_quota_finland,
        epsilon_h_sw,
        epsilon_h_fi,
        pi_s,
        pi_c,
        t_mate_to_preg,
        t_birth_to_end_hunt,
        rel_tol,
        abs_tol,
        max_num_steps,
        fate_probability_tolerance
    )

    mu_m <- -log(S_diag)

    # The joint-ODE discrete signature no longer accepts S_diag,
    # ode_init_state, or ode_times.
    result <- m$functions$run_state_process_from_survivors_rng(
        n_state_years,
        n_age_classes,
        survivors_final,
        birth_rate_first,
        baseline_birth_rate,
        density_dependence_intercept,
        density_dependence_slope,
        mu_m,
        hunting_selectivity_sweden,
        hunting_selectivity_finland,
        as.integer(hunting_quota_sweden),
        as.integer(hunting_quota_finland),
        hunting_effort_sd_sweden,
        hunting_effort_sd_finland,
        epsilon_h_sw,
        epsilon_h_fi,
        t_mate_to_preg,
        t_birth_to_end_hunt,
        pi_s,
        pi_c,
        prob_of_ca,
        rel_tol,
        abs_tol,
        as.integer(max_num_steps),
        fate_probability_tolerance
    )

    result <- name_state_process_output(result)

    fate_sum <-
        result$survivors +
        result$deaths_or_bycatch +
        result$hunted_sweden +
        result$hunted_finland

    stopifnot(
        all(is.finite(result$population_total)),
        all(is.finite(result$population_comp)),
        all(is.finite(result$bycatch_expected)),
        all(result$population_total >= 0),
        all(result$population_comp >= 0),
        all(result$survivors >= 0),
        all(result$deaths_or_bycatch >= 0),
        all(result$hunted_sweden >= 0),
        all(result$hunted_finland >= 0),
        all(result$bycatch_expected >= 0),
        all(result$population_comp == round(result$population_comp)),
        all(result$survivors == round(result$survivors)),
        all(result$deaths_or_bycatch == round(result$deaths_or_bycatch)),
        all(result$hunted_sweden == round(result$hunted_sweden)),
        all(result$hunted_finland == round(result$hunted_finland)),
        isTRUE(all.equal(fate_sum, result$population_comp, tolerance = 1e-8))
    )

    result
}

# -----------------------------------------------------------------------------
# Plot one forecast independently
# -----------------------------------------------------------------------------

plot_forecast <- function(
    forecast,
    process = c("continuous", "discrete"),
    survivor_boundary = NULL
) {
    process <- match.arg(process)
    required_outputs <- c(
        "population_total",
        "birth_rate",
        "pregnancy_rate",
        "hunted_total",
        "hunting_bag_total_sweden",
        "hunting_bag_total_finland"
    )
    missing_outputs <- setdiff(required_outputs, names(forecast))
    if (length(missing_outputs) > 0L) {
        stop(
            "The forecast is missing: ",
            paste(missing_outputs, collapse = ", "),
            call. = FALSE
        )
    }

    n_state_years <- length(forecast$population_total)
    process_label <- if (process == "continuous") {
        "Continuous approximation"
    } else {
        "Discrete process"
    }
    process_colour <- if (process == "continuous") "steelblue" else "firebrick"

    population_data <- data.frame(
        year = seq_len(n_state_years),
        population = forecast$population_total
    )
    rate_data <- rbind(
        data.frame(
            year = seq_len(n_state_years),
            quantity = "Birth rate",
            value = forecast$birth_rate
        ),
        data.frame(
            year = seq_len(n_state_years),
            quantity = "Pregnancy rate",
            value = forecast$pregnancy_rate
        )
    )
    hunting_data <- rbind(
        data.frame(
            year = seq_len(n_state_years),
            quantity = "Total hunted",
            value = forecast$hunted_total
        ),
        data.frame(
            year = seq_len(n_state_years),
            quantity = "Swedish hunting bag",
            value = forecast$hunting_bag_total_sweden
        ),
        data.frame(
            year = seq_len(n_state_years),
            quantity = "Finnish hunting bag",
            value = forecast$hunting_bag_total_finland
        )
    )

    population_plot <- ggplot(population_data, aes(year, population)) +
        geom_line(colour = process_colour, linewidth = 0.9) +
        geom_point(colour = process_colour, size = 1.8) +
        scale_x_continuous(breaks = seq_len(n_state_years)) +
        labs(
            x = "Forecast year",
            y = "Total population",
            title = paste(process_label, "population trajectory")
        ) +
        theme_minimal()

    if (!is.null(survivor_boundary)) {
        population_plot <- population_plot +
            geom_point(
                data = data.frame(year = 0, population = survivor_boundary),
                aes(year, population),
                inherit.aes = FALSE,
                shape = 4,
                size = 3,
                colour = "black"
            ) +
            scale_x_continuous(breaks = 0:n_state_years) +
            labs(
                subtitle = "Year-0 cross is the fitted survivor boundary; the connected trajectory begins at year 1"
            )
    }

    rates_plot <- ggplot(rate_data, aes(year, value)) +
        geom_line(colour = process_colour, linewidth = 0.9) +
        geom_point(colour = process_colour, size = 1.5) +
        facet_wrap(vars(quantity), scales = "free_y", ncol = 1) +
        scale_x_continuous(breaks = seq_len(n_state_years)) +
        labs(
            x = "Forecast year",
            y = NULL,
            title = paste(process_label, "reproductive rates")
        ) +
        theme_minimal()

    hunting_plot <- ggplot(hunting_data, aes(year, value, colour = quantity)) +
        geom_line(linewidth = 0.9) +
        geom_point(size = 1.5) +
        scale_x_continuous(breaks = seq_len(n_state_years)) +
        labs(
            x = "Forecast year",
            y = "Count",
            colour = NULL,
            title = paste(process_label, "hunting outcomes")
        ) +
        theme_minimal() +
        theme(legend.position = "top")

    list(
        population = population_plot,
        rates = rates_plot,
        hunting = hunting_plot,
        combined = population_plot +
            rates_plot +
            hunting_plot +
            plot_layout(widths = c(1, 1, 1))
    )
}

run_and_plot_continuous_forecast <- function(...) {
    call_args <- list(...)
    forecast <- do.call(run_continuous_forecast, call_args)
    survivors <- if ("survivors_final" %in% names(call_args)) {
        call_args$survivors_final
    } else {
        .default_survivors_final
    }
    list(
        forecast = forecast,
        plots = plot_forecast(forecast, "continuous", sum(survivors))
    )
}

run_and_plot_discrete_forecast <- function(...) {
    call_args <- list(...)
    forecast <- do.call(run_discrete_forecast, call_args)
    survivors <- if ("survivors_final" %in% names(call_args)) {
        call_args$survivors_final
    } else {
        .default_survivors_final
    }
    list(
        forecast = forecast,
        plots = plot_forecast(forecast, "discrete", sum(survivors))
    )
}

# -----------------------------------------------------------------------------
# Compare repeated independent forecasts
# -----------------------------------------------------------------------------

compare_forecasts <- function(
    n_state_years = 10L,
    n_age_classes = 6L,
    n_replicates = 200L,
    master_seed = 123L,
    survivors_final = .default_survivors_final,
    birth_rate_first = 0.82,
    baseline_birth_rate = rep(0.86, n_state_years + 1L),
    density_dependence_intercept = 0.05,
    density_dependence_slope = 1e-5,
    S_diag = .default_S_diag,
    hunting_selectivity_sweden = rep(0, 2L * n_age_classes),
    hunting_selectivity_finland = rep(0, 2L * n_age_classes),
    hunting_quota_sweden = rep(1000L, n_state_years),
    hunting_quota_finland = rep(1000L, n_state_years),
    hunting_effort_sd_sweden = 0.2,
    hunting_effort_sd_finland = 0.2,
    rel_tol = 1e-8,
    abs_tol = 1e-8,
    max_num_steps = 100000L,
    fate_probability_tolerance = .default_fate_probability_tolerance
) {
    set.seed(master_seed)
    epsilon_h_sw <- abs(rnorm(n_state_years))
    epsilon_h_fi <- abs(rnorm(n_state_years))

    common_args <- list(
        n_state_years = n_state_years,
        n_age_classes = n_age_classes,
        survivors_final = survivors_final,
        birth_rate_first = birth_rate_first,
        baseline_birth_rate = baseline_birth_rate,
        density_dependence_intercept = density_dependence_intercept,
        density_dependence_slope = density_dependence_slope,
        S_diag = S_diag,
        hunting_selectivity_sweden = hunting_selectivity_sweden,
        hunting_selectivity_finland = hunting_selectivity_finland,
        hunting_quota_sweden = hunting_quota_sweden,
        hunting_quota_finland = hunting_quota_finland,
        hunting_effort_sd_sweden = hunting_effort_sd_sweden,
        hunting_effort_sd_finland = hunting_effort_sd_finland,
        epsilon_h_sw = epsilon_h_sw,
        epsilon_h_fi = epsilon_h_fi,
        rel_tol = rel_tol,
        abs_tol = abs_tol,
        max_num_steps = max_num_steps,
        fate_probability_tolerance = fate_probability_tolerance
    )

    seeds <- master_seed + seq_len(n_replicates) - 1L
    continuous_results <- lapply(seeds, function(seed) {
        do.call(run_continuous_forecast, c(list(seed = seed), common_args))
    })
    discrete_results <- lapply(seeds, function(seed) {
        do.call(run_discrete_forecast, c(list(seed = seed), common_args))
    })

    extract_matrix <- function(results, quantity) {
        vapply(results, function(x) x[[quantity]], numeric(n_state_years))
    }
    summarize_draws <- function(draws, process, quantity) {
        data.frame(
            year = seq_len(nrow(draws)),
            process = process,
            quantity = quantity,
            mean = rowMeans(draws),
            lower = apply(draws, 1L, quantile, probs = 0.05),
            median = apply(draws, 1L, median),
            upper = apply(draws, 1L, quantile, probs = 0.95)
        )
    }

    continuous_population_draws <- extract_matrix(
        continuous_results,
        "population_total"
    )
    discrete_population_draws <- extract_matrix(
        discrete_results,
        "population_total"
    )
    continuous_birth_rate_draws <- extract_matrix(
        continuous_results,
        "birth_rate"
    )
    discrete_birth_rate_draws <- extract_matrix(discrete_results, "birth_rate")
    continuous_hunted_draws <- extract_matrix(
        continuous_results,
        "hunted_total"
    )
    discrete_hunted_draws <- extract_matrix(discrete_results, "hunted_total")

    population_summary <- rbind(
        summarize_draws(
            continuous_population_draws,
            "Continuous approximation",
            "Total population"
        ),
        summarize_draws(
            discrete_population_draws,
            "Discrete process",
            "Total population"
        )
    )
    quantity_summary <- rbind(
        summarize_draws(
            continuous_birth_rate_draws,
            "Continuous approximation",
            "Birth rate"
        ),
        summarize_draws(
            discrete_birth_rate_draws,
            "Discrete process",
            "Birth rate"
        ),
        summarize_draws(
            continuous_hunted_draws,
            "Continuous approximation",
            "Total hunted"
        ),
        summarize_draws(
            discrete_hunted_draws,
            "Discrete process",
            "Total hunted"
        )
    )

    continuous_result <- continuous_results[[1L]]
    discrete_result <- discrete_results[[1L]]
    single_data <- rbind(
        data.frame(
            year = seq_len(n_state_years),
            process = "Continuous approximation",
            population = continuous_result$population_total
        ),
        data.frame(
            year = seq_len(n_state_years),
            process = "Discrete process",
            population = discrete_result$population_total
        )
    )
    mean_difference <- data.frame(
        year = seq_len(n_state_years),
        continuous_mean = rowMeans(continuous_population_draws),
        discrete_mean = rowMeans(discrete_population_draws)
    )
    mean_difference$difference <- mean_difference$discrete_mean -
        mean_difference$continuous_mean
    mean_difference$relative_difference <- mean_difference$difference /
        mean_difference$continuous_mean

    colours <- c(
        "Continuous approximation" = "steelblue",
        "Discrete process" = "firebrick"
    )
    p1 <- ggplot(single_data, aes(year, population, colour = process)) +
        geom_line(linewidth = 0.9) +
        geom_point(size = 1.8) +
        geom_point(
            data = data.frame(year = 0, population = sum(survivors_final)),
            aes(year, population),
            inherit.aes = FALSE,
            shape = 4,
            size = 3
        ) +
        scale_x_continuous(breaks = 0:n_state_years) +
        scale_colour_manual(values = colours) +
        labs(
            x = "Forecast year",
            y = "Total population",
            colour = NULL,
            title = "Example trajectories",
            subtitle = "Year-0 cross is the common fitted survivor boundary"
        ) +
        theme_minimal() +
        theme(legend.position = "top")

    p2 <- ggplot(
        quantity_summary,
        aes(year, median, colour = process, fill = process)
    ) +
        geom_ribbon(
            aes(ymin = lower, ymax = upper),
            alpha = 0.16,
            colour = NA
        ) +
        geom_line(linewidth = 0.9) +
        facet_wrap(vars(quantity), scales = "free_y", ncol = 1) +
        scale_colour_manual(values = colours) +
        scale_fill_manual(values = colours) +
        labs(
            x = "Forecast year",
            y = NULL,
            colour = NULL,
            fill = NULL,
            title = "Birth-rate and hunting distributions"
        ) +
        theme_minimal() +
        theme(legend.position = "top")

    p3 <- ggplot(
        population_summary,
        aes(year, median, colour = process, fill = process)
    ) +
        geom_ribbon(
            aes(ymin = lower, ymax = upper),
            alpha = 0.16,
            colour = NA
        ) +
        geom_line(linewidth = 0.9) +
        geom_line(aes(y = mean), linewidth = 0.65, linetype = "dashed") +
        scale_colour_manual(values = colours) +
        scale_fill_manual(values = colours) +
        labs(
            x = "Forecast year",
            y = "Total population",
            colour = NULL,
            fill = NULL,
            title = "Population distributions",
            subtitle = paste0(
                "Median and 90% intervals from ",
                n_replicates,
                " simulations; dashed lines are means"
            )
        ) +
        theme_minimal() +
        theme(legend.position = "top")

    p4 <- ggplot(mean_difference, aes(year, relative_difference)) +
        geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
        geom_line(colour = "purple4", linewidth = 0.9) +
        geom_point(colour = "purple4", size = 1.8) +
        scale_y_continuous(labels = scales::label_percent()) +
        labs(
            x = "Forecast year",
            y = "Relative mean difference",
            title = "Discrete minus continuous population mean",
            subtitle = "Difference divided by the continuous-process mean"
        ) +
        theme_minimal()

    comparison_plot <- p1 +
        p2 +
        p3 +
        p4 +
        plot_layout(guides = "collect") &
        theme(legend.position = "top")

    structure(
        list(
            continuous_result = continuous_result,
            discrete_result = discrete_result,
            continuous_results = continuous_results,
            discrete_results = discrete_results,
            continuous_population_draws = continuous_population_draws,
            discrete_population_draws = discrete_population_draws,
            population_summary = population_summary,
            quantity_summary = quantity_summary,
            mean_difference = mean_difference,
            survivor_boundary = sum(survivors_final),
            plots = list(
                single_trajectory = p1,
                quantity_distribution = p2,
                population_distribution = p3,
                mean_difference = p4,
                comparison = comparison_plot
            )
        ),
        class = "seal_ipm_forecast_comparison"
    )
}

# Examples:
# continuous <- run_and_plot_continuous_forecast(seed = 123L)
# print(continuous$plots$combined)
#
# discrete <- run_and_plot_discrete_forecast(seed = 123L)
# print(discrete$plots$combined)
#
# comparison <- compare_forecasts(n_replicates = 200L)
# print(comparison$plots$comparison)
