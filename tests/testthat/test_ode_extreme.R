library(cmdstanr)

stan_file <- "src/stan/grey_seal_IPM_forecast.stan"

model <- cmdstan_model(
    stan_file,
    force_recompile = TRUE
)
model$expose_functions()

fate_probability <- model$functions$compute_fate_probabilities(
    1.5 / 12,
    8 / 12,
    904.299,
    0,
    0.125893,
    1e-8,
    1e-10,
    100000L,
    1e-6
)

survival_to_hunting_start <- exp(-0.125893 * (1.5 / 12))
reference_death_before_hunting <- 1 - survival_to_hunting_start
reference_death_during_hunting <- survival_to_hunting_start * stats::integrate(
    function(tau) {
        0.125893 * exp(
            -904.299 * ((8 / 12) * tau - tau^2 / 2) -
                0.125893 * tau
        )
    },
    lower = 0,
    upper = 8 / 12,
    rel.tol = 1e-11,
    abs.tol = 1e-13
)$value

alive_at_end_hunting <- exp(
    -904.299 * (8 / 12)^2 / 2 - 0.125893 * (1.5 / 12 + 8 / 12)
)
reference_death_after_hunting <- alive_at_end_hunting *
    (1 - exp(-0.125893 * (1 - 1.5 / 12 - 8 / 12)))

stopifnot(
    all(is.finite(fate_probability)),
    all(fate_probability >= 0),
    isTRUE(all.equal(sum(fate_probability), 1, tolerance = 1e-12)),
    # Death is a remainder of ODE probabilities, so compare absolute error
    abs(fate_probability[4] -
        (reference_death_before_hunting + reference_death_during_hunting +
            reference_death_after_hunting)) < 1e-8
)

# With very small mortality, death is a near-zero remainder of the original
# hunting ODE solves. Small solver boundary errors are handled by tolerance.
near_zero_death_probability <-
    model$functions$compute_fate_probabilities(
        1.5 / 12,
        8 / 12,
        1e6,
        1e6,
        1e-12,
        1e-8,
        1e-10,
        100000L,
        1e-6
    )

stopifnot(
    all(is.finite(near_zero_death_probability)),
    all(near_zero_death_probability >= 0),
    isTRUE(all.equal(
        sum(near_zero_death_probability),
        1,
        tolerance = 1e-12
    ))
)

# Exercise a much more concentrated hunting hazard than the reported case.
# The restored RK45 solves must conserve probability within tolerance.
extreme_fate_probability <-
    model$functions$compute_fate_probabilities(
        1.5 / 12,
        8 / 12,
        1e6,
        1e6,
        0.5,
        1e-8,
        1e-10,
        100000L,
        1e-6
    )

stopifnot(
    all(is.finite(extreme_fate_probability)),
    all(extreme_fate_probability >= 0),
    isTRUE(all.equal(
        sum(extreme_fate_probability),
        1,
        tolerance = 1e-12
    ))
)
