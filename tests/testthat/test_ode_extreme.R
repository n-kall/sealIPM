library(cmdstanr)

stan_file <- "src/stan/grey_seal_IPM_forecast.stan"

model <- cmdstan_model(
    stan_file,
    force_recompile = TRUE
)
model$expose_functions()

fate_probability <- model$functions$compute_joint_fate_probabilities(
    8 / 12,
    904.299,
    0,
    0.125893,
    1e-8,
    1e-10,
    100000L,
    1e-6
)

reference_death_during_hunting <- stats::integrate(
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
    -904.299 * (8 / 12)^2 / 2 - 0.125893 * (8 / 12)
)
reference_death_after_hunting <- alive_at_end_hunting *
    (1 - exp(-0.125893 * (1 - 8 / 12)))

stopifnot(
    all(is.finite(fate_probability)),
    all(fate_probability >= 0),
    isTRUE(all.equal(sum(fate_probability), 1, tolerance = 1e-12)),
    isTRUE(all.equal(
        fate_probability[2],
        reference_death_during_hunting + reference_death_after_hunting,
        tolerance = 1e-8
    ))
)

# Very small non-hunting mortality produces a probability close to zero. The
# unscaled integrate_1d implementation failed here because a relative-error
# target on a near-zero integral implied an impractically small absolute error.
near_zero_death_probability <-
    model$functions$compute_joint_fate_probabilities(
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

# Exercise a much stiffer case than the reported failure. The old scalar BDF
# solve could overshoot below zero here; positive quadrature must remain on the
# probability simplex.
extreme_fate_probability <-
    model$functions$compute_joint_fate_probabilities(
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
