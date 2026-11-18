initialize_forecast_refit <- function(init, stan_data) {
  positive <- which(stan_data$obs_aerial_count > 0)
  if (!length(positive)) {
    return(init)
  }

  first <- positive[which.min(stan_data$aerial_year[positive])]
  detection <- init$aerial_count_mu
  initial_population <- stan_data$obs_aerial_count[first] / detection
  population_peak <- max(stan_data$obs_aerial_count) / detection

  # Later windows need starting populations on the scale of their observations
  init$population_init_size <- initial_population
  init$carrying_capacity <- pmax(init$carrying_capacity, population_peak)
  init
}
