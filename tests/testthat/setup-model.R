skip_if_not_installed("cmdstanr")

stan_file <- system.file(
    "bin",
    "stan",
    "grey_seal_IPM.stan",
    package = "sealIPM",
    mustWork = TRUE
)

test_model <- cmdstanr::cmdstan_model(
    stan_file,
    quiet = TRUE,
    compile = TRUE,
    force_recompile = TRUE
)

test_model$expose_functions()
