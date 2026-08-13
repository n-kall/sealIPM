library(sealIPM)

dat <- data(grey_seal_data)

grey_seal_ipm_fit <- fit_ipm(
  data = grey_seal_data,
  species = "grey",
  years = 2003:2024,
  iter_warmup = 1000,
  iter_sampling = 1000,
  chains = 4
)

usethis::use_data(grey_seal_ipm_fit, overwrite = TRUE)
