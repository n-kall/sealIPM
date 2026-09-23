functions {
  #include functions/ode.stanfunctions
  #include functions/logitnormal.stanfunctions
  #include functions/observation_models/aerial_counts.stanfunctions
  #include functions/observation_models/hunting_bags.stanfunctions
  #include functions/observation_models/hunting_composition.stanfunctions
  #include functions/observation_models/bycatch_composition.stanfunctions
  #include functions/observation_models/pregnancy.stanfunctions
  #include functions/observation_models/reproductive_signs.stanfunctions
  #include functions/statespace.stanfunctions
  #include functions/statespace_discrete.stanfunctions
}

data {
  // ------------------------------
  // FORECAST HORIZON AND DIMENSIONS
  // ------------------------------
  int<lower=1> n_future_years; // number of years in the forecast horizon
  int<lower=1> n_age_classes;  // number of age classes per sex

  // Environmental covariates. The additional element is needed because
  // pregnancy in forecast year y depends on the baseline birth rate for y + 1.
  vector[n_future_years + 1]
    future_herring_index_baltic_proper_gulf_finland;
  vector[n_future_years + 1]
    future_herring_index_gulf_bothnia;

  // Country-specific annual hunting quotas.
  array[n_future_years] int<lower=0>
    future_hunting_quota_sweden;
  array[n_future_years] int<lower=0>
    future_hunting_quota_finland;

  // Hunting and reproductive timing, expressed as fractions of a year.
  real<lower=0, upper=1> t_mate_to_preg;
  real<lower=0, upper=1> t_birth_to_start_hunt;
  real<lower=0, upper=1> t_birth_to_end_hunt;
  real<lower=0, upper=1> t_hunt;

  // ODE solver controls used by the joint competing-risks fate ODE.
  real<lower=1e-15> rel_tol;
  real<lower=1e-15> abs_tol;
  int<lower=1> max_num_steps;
  // Maximum accepted numerical deviation from the probability simplex
  // returned by the joint fate ODE. Larger deviations cause rejection.
  real<lower=1e-15> fate_probability_tolerance;

  // Observation-model coefficient of variation for reported hunting bags.
  real<lower=0> harvest_bag_cv;

  // Sample sizes used for posterior predictive observation generation.
  array[n_future_years] int<lower=0>
    future_hunting_sample_size_sweden;
  array[n_future_years] int<lower=0>
    future_hunting_sample_size_finland;
  array[n_future_years] int<lower=0>
    future_bycatch_sample_size;
  array[n_future_years] int<lower=0>
    future_reproductive_signs_sample_size;
  array[n_future_years] int<lower=0>
    future_pregnancy_sample_size;

  // ------------------------------
  // LEAVE-FUTURE-OUT OBSERVATIONS
  // ------------------------------
  // Each count may be zero when that observation type is unavailable.
  int<lower=0> n_future_aerial;
  int<lower=0> n_future_hunting_bag_sweden;
  int<lower=0> n_future_hunting_bag_finland;
  int<lower=0> n_future_hunting_comp_sweden;
  int<lower=0> n_future_hunting_comp_finland;
  int<lower=0> n_future_bycatch;
  int<lower=0> n_future_pregnancy;
  int<lower=0> n_future_reproductive_signs;

  // Aerial survey observations.
  array[n_future_aerial]
    int<lower=1, upper=n_future_years> aerial_year;
  array[n_future_aerial]
    int<lower=0> future_obs_aerial_count;

  // Swedish hunting-bag observations.
  array[n_future_hunting_bag_sweden]
    int<lower=1, upper=n_future_years> hunting_bag_year_sweden;
  array[n_future_hunting_bag_sweden]
    real<lower=0> future_obs_hunting_bag_sweden;

  // Finnish hunting-bag observations.
  array[n_future_hunting_bag_finland]
    int<lower=1, upper=n_future_years> hunting_bag_year_finland;
  array[n_future_hunting_bag_finland]
    real<lower=0> future_obs_hunting_bag_finland;

  // Swedish hunting-composition observations, stored year-major.
  array[n_future_hunting_comp_sweden]
    int<lower=1, upper=n_future_years> hunting_comp_year_sweden;
  array[n_future_hunting_comp_sweden, 2 * n_age_classes]
    int<lower=0> future_obs_hunting_comp_sweden;

  // Finnish hunting-composition observations, stored year-major.
  array[n_future_hunting_comp_finland]
    int<lower=1, upper=n_future_years> hunting_comp_year_finland;
  array[n_future_hunting_comp_finland, 2 * n_age_classes]
    int<lower=0> future_obs_hunting_comp_finland;

  // Bycatch-composition observations, stored year-major.
  array[n_future_bycatch]
    int<lower=1, upper=n_future_years> bycatch_year;
  array[n_future_bycatch, 2 * n_age_classes]
    int<lower=0> future_obs_bycatch_comp;

  // Pregnancy observations.
  array[n_future_pregnancy]
    int<lower=1, upper=n_future_years> pregnancy_count_year;
  array[n_future_pregnancy]
    int<lower=0> future_obs_pregnancy_count;
  array[n_future_pregnancy]
    int<lower=0> future_obs_pregnancy_sample_size;

  // Reproductive-sign observations.
  array[n_future_reproductive_signs]
    int<lower=1, upper=n_future_years> reproductive_signs_year;
  array[n_future_reproductive_signs, 4]
    int<lower=0> future_obs_reproductive_signs_finland;
}

transformed data {
  // Consecutive state-process indices used by observation RNG functions.
  array[n_future_years] int future_year;
  int n_demo_groups = 2 * n_age_classes;

  for (year in 1:n_future_years) {
    future_year[year] = year;
  }
}

parameters {

  // Aerial-survey observation parameters.
  real<lower=0, upper=1> aerial_count_mu;
  real<lower=0> aerial_count_overdispersion;

  // Hunting-effort process parameters.
  real<lower=0> hunting_effort_sd_sweden;
  real<lower=0> hunting_effort_sd_finland;

  // Herring-dependent baseline birth-rate parameters.
  real<lower=0, upper=1> max_baseline_birth_rate;
  real<lower=0, upper=1> min_baseline_birth_rate_prop;
  real herring_intercept_scaled;
  real herring_slope;
  real<lower=0, upper=1> herring_weight;
  real<lower=0, upper=1> density_dependence_scaled;

  // Reproductive-sign reporting parameters.
  real<lower=0, upper=1> report_ca_mean;
  real<lower=0, upper=1> report_placental_mean;
  real<lower=0, upper=1> prob_of_ca;
  real<lower=0> report_placental_sd;
  real<lower=0> report_ca_sd;

  // Derived fitted demographic parameters carried into the forecast.
  real density_dependence_intercept;
  vector<lower=0>[n_demo_groups] mu_m;
  vector[n_demo_groups] S_diag;
  vector[n_demo_groups] hunting_selectivity_finland;
  vector[n_demo_groups] hunting_selectivity_sweden;
  vector[n_demo_groups] bycatch_bias;

  real<lower=0, upper=1> birth_rate_at_carrying_capacity;
  real<lower=0, upper=1> min_baseline_birth_rate_actual;
  real density_dependence_slope;

  // Final fitted states. The discrete forecast starts from survivors_final.
  // population_comp_final and population_total_final are retained for output
  // compatibility and for calculating the first forecast birth rate.
  vector<lower=0>[n_demo_groups] population_comp_final;
  real<lower=0> population_total_final;
  real<lower=0, upper=1> pregnancy_rate_final;
  vector<lower=0>[n_demo_groups] survivors_final;
}

generated quantities {
  // ------------------------------
  // BIRTH-RATE INPUT
  // ------------------------------
  vector<lower=0, upper=1>[n_future_years + 1]
    baseline_birth_rate_future = compute_baseline_birth_rate(
      min_baseline_birth_rate_actual,
      max_baseline_birth_rate,
      herring_intercept_scaled,
      herring_slope,
      herring_weight,
      future_herring_index_baltic_proper_gulf_finland,
      future_herring_index_gulf_bothnia
    );

  // ------------------------------
  // FUTURE PROCESS INNOVATIONS
  // ------------------------------
  vector<lower=0>[n_future_years] epsilon_h_sw_future;
  vector<lower=0>[n_future_years] epsilon_h_fi_future;
  vector<lower=0>[n_future_years] epsilon_placental_future;
  vector<lower=0>[n_future_years] epsilon_ca_future;

  // ------------------------------
  // FUTURE STATE-PROCESS OUTPUTS
  // ------------------------------
  vector<lower=0, upper=1>[n_future_years] birth_rate_future;
  vector<lower=0, upper=1>[n_future_years] pregnancy_rate_future;
  vector<lower=0>[n_future_years] population_total_future;
  matrix<lower=0>[n_demo_groups, n_future_years] population_comp_future;
  matrix<lower=0>[n_demo_groups, n_future_years] survivors_future;
  matrix<lower=0>[n_demo_groups, n_future_years] deaths_or_bycatch_future;
  matrix<lower=0>[n_demo_groups, n_future_years] hunted_sweden_future;
  matrix<lower=0>[n_demo_groups, n_future_years] hunted_finland_future;
  matrix<lower=0>[n_demo_groups, n_future_years] bycatch_expected_future;
  vector<lower=0>[n_future_years] hunting_bag_total_sweden_future;
  vector<lower=0>[n_future_years] hunting_bag_total_finland_future;
  vector<lower=0>[n_future_years] hunted_total_future;
  matrix<lower=0, upper=1>[4, n_future_years] reproductive_probs_future;

  for (year in 1:n_future_years) {
    epsilon_h_sw_future[year] = abs(std_normal_rng());
    epsilon_h_fi_future[year] = abs(std_normal_rng());
    epsilon_placental_future[year] = abs(std_normal_rng());
    epsilon_ca_future[year] = abs(std_normal_rng());
  }

  // Annual reporting probabilities for reproductive signs.
  vector<lower=0, upper=1>[n_future_years] pi_s_future =
    report_placental_mean
    * exp(-epsilon_placental_future * report_placental_sd);

  vector<lower=0, upper=1>[n_future_years] pi_c_future =
    report_ca_mean
    * exp(-epsilon_ca_future * report_ca_sd);

  // The first future birth rate is updated from the final fitted population.
  real<lower=0, upper=1> birth_rate_future_first = update_birth_rate(
    baseline_birth_rate_future[1],
    density_dependence_intercept,
    density_dependence_slope,
    population_total_final
  );

  // -------------------------------------------------------------------------
  // DISCRETE DEMOGRAPHIC FORECAST
  // -------------------------------------------------------------------------
  (
    birth_rate_future,
    pregnancy_rate_future,
    population_total_future,
    population_comp_future,
    survivors_future,
    deaths_or_bycatch_future,
    hunted_sweden_future,
    hunted_finland_future,
    bycatch_expected_future,
    hunting_bag_total_sweden_future,
    hunting_bag_total_finland_future,
    hunted_total_future,
    reproductive_probs_future
  ) = run_state_process_from_survivors_rng(
    n_future_years,
    n_age_classes,
    survivors_final,
    birth_rate_future_first,
    baseline_birth_rate_future,
    density_dependence_intercept,
    density_dependence_slope,
    mu_m,
    hunting_selectivity_sweden,
    hunting_selectivity_finland,
    future_hunting_quota_sweden,
    future_hunting_quota_finland,
    hunting_effort_sd_sweden,
    hunting_effort_sd_finland,
    epsilon_h_sw_future,
    epsilon_h_fi_future,
    t_mate_to_preg,
    t_birth_to_end_hunt,
    pi_s_future,
    pi_c_future,
    prob_of_ca,
    rel_tol,
    abs_tol,
    max_num_steps,
    fate_probability_tolerance
  );

  // ------------------------------
  // POSTERIOR PREDICTIVE OBSERVATIONS
  // ------------------------------
  array[n_future_years] int future_aerial_count = aerial_count_rng(
    future_year,
    population_total_future,
    aerial_count_mu,
    aerial_count_overdispersion
  );

  array[n_future_years] real future_harvest_bags_sweden =
    harvest_bags_rng(
      future_year,
      hunting_bag_total_sweden_future,
      harvest_bag_cv
    );

  array[n_future_years] real future_harvest_bags_finland =
    harvest_bags_rng(
      future_year,
      hunting_bag_total_finland_future,
      harvest_bag_cv
    );

  array[n_future_years, n_demo_groups] int future_hunting_comp_finland =
    hunting_comp_rng(
      future_year,
      hunted_finland_future,
      hunting_bag_total_finland_future,
      future_hunting_sample_size_finland
    );

  array[n_future_years, n_demo_groups] int future_hunting_comp_sweden =
    hunting_comp_rng(
      future_year,
      hunted_sweden_future,
      hunting_bag_total_sweden_future,
      future_hunting_sample_size_sweden
    );

  array[n_future_years, n_demo_groups] int future_bycatch_comp =
    bycatch_comp_rng(
      future_year,
      bycatch_expected_future,
      bycatch_bias,
      future_bycatch_sample_size
    );

  array[n_future_years] int future_pregnancy = pregnancy_rng(
    future_year,
    future_pregnancy_sample_size,
    pregnancy_rate_future
  );

  array[n_future_years, 4] int future_reproductive_signs =
    reproductive_signs_rng(
      future_year,
      reproductive_probs_future,
      future_reproductive_signs_sample_size
    );

  // ------------------------------
  // LEAVE-FUTURE-OUT LOG LIKELIHOODS
  // ------------------------------
  // Vectors are initialized to zero so observation types with no records
  // make an exact zero contribution to the joint log likelihood
  vector[n_future_aerial] log_lik_future_aerial =
    rep_vector(0.0, n_future_aerial);
  vector[n_future_hunting_bag_sweden]
    log_lik_future_harvest_bags_sweden =
      rep_vector(0.0, n_future_hunting_bag_sweden);
  vector[n_future_hunting_bag_finland]
    log_lik_future_harvest_bags_finland =
      rep_vector(0.0, n_future_hunting_bag_finland);
  vector[n_future_hunting_comp_sweden]
    log_lik_future_hunting_comp_sweden =
      rep_vector(0.0, n_future_hunting_comp_sweden);
  vector[n_future_hunting_comp_finland]
    log_lik_future_hunting_comp_finland =
      rep_vector(0.0, n_future_hunting_comp_finland);
  vector[n_future_bycatch] log_lik_future_bycatch =
    rep_vector(0.0, n_future_bycatch);
  vector[n_future_pregnancy] log_lik_future_pregnancy =
    rep_vector(0.0, n_future_pregnancy);
  vector[n_future_reproductive_signs]
    log_lik_future_reproductive_signs =
      rep_vector(0.0, n_future_reproductive_signs);

  if (n_future_aerial > 0) {
    log_lik_future_aerial = aerial_count_pointwise_log_lik(
      future_obs_aerial_count,
      aerial_year,
      population_total_future,
      aerial_count_mu,
      aerial_count_overdispersion
    );
  }

  if (n_future_hunting_bag_sweden > 0) {
    log_lik_future_harvest_bags_sweden =
      harvest_bags_pointwise_log_lik(
        future_obs_hunting_bag_sweden,
        hunting_bag_year_sweden,
        hunting_bag_total_sweden_future,
        harvest_bag_cv
      );
  }

  if (n_future_hunting_bag_finland > 0) {
    log_lik_future_harvest_bags_finland =
      harvest_bags_pointwise_log_lik(
        future_obs_hunting_bag_finland,
        hunting_bag_year_finland,
        hunting_bag_total_finland_future,
        harvest_bag_cv
      );
  }

  if (n_future_hunting_comp_sweden > 0) {
    log_lik_future_hunting_comp_sweden =
      hunting_comp_pointwise_log_lik(
        future_obs_hunting_comp_sweden,
        hunting_comp_year_sweden,
        hunted_sweden_future,
        hunting_bag_total_sweden_future
      );
  }

  if (n_future_hunting_comp_finland > 0) {
    log_lik_future_hunting_comp_finland =
      hunting_comp_pointwise_log_lik(
        future_obs_hunting_comp_finland,
        hunting_comp_year_finland,
        hunted_finland_future,
        hunting_bag_total_finland_future
      );
  }

  if (n_future_bycatch > 0) {
    log_lik_future_bycatch = bycatch_comp_pointwise_log_lik(
      future_obs_bycatch_comp,
      bycatch_year,
      bycatch_expected_future,
      bycatch_bias
    );
  }

  if (n_future_pregnancy > 0) {
    log_lik_future_pregnancy = pregnancy_pointwise_log_lik(
      future_obs_pregnancy_count,
      pregnancy_count_year,
      future_obs_pregnancy_sample_size,
      pregnancy_rate_future
    );
  }

  if (n_future_reproductive_signs > 0) {
    log_lik_future_reproductive_signs =
      reproductive_signs_pointwise_log_lik(
        future_obs_reproductive_signs_finland,
        reproductive_signs_year,
        reproductive_probs_future
      );
  }

  // Combined conditional log likelihood for the simulated future state path.
  real log_lik =
    sum(log_lik_future_aerial)
    + sum(log_lik_future_harvest_bags_sweden)
    + sum(log_lik_future_harvest_bags_finland)
    + sum(log_lik_future_hunting_comp_sweden)
    + sum(log_lik_future_hunting_comp_finland)
    + sum(log_lik_future_bycatch)
    + sum(log_lik_future_pregnancy)
    + sum(log_lik_future_reproductive_signs);
}
