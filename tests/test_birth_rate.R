library(cmdstanr)
library(tidyverse)
library(purrr)

m <- cmdstan_model("src/stan/grey_seal_IPM.stan", force_recompile = TRUE)

m$expose_functions()

pup_survivals <- seq(0.1, 0.85, by = 0.05)

birth_rate <- map(pup_survivals, function(.x) {
    mu_m <- m$functions$mortality_rates(
        phi_pup = .x,
        phi_adult = 0.9,
        c = 0.5,
        n_age_classes = 6,
        male_pup_offset = 0,
        male_adult_offset = 0
    )

    m$functions$euler_lotka_birth_rate(
        mu_m = mu_m,
        n_age_classes = 6,
        phi_a = 0.9
    )
})

tibble(pup_survival = pup_survivals, bk = unlist(birth_rate)) |>
    ggplot(aes(x = pup_survival, y = bk)) +
    geom_line() +
    geom_hline(yintercept = 1, linetype = "dashed") +
    ggtitle(
        "Birth rate at carrying capacity (b_K) as a function of pup_survival (phi_pup)",
        "other parameters set to prior means"
    )
