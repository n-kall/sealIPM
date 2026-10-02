functions {
  #include functions/observation_models/aerial_counts.stanfunctions
  #include functions/observation_models/hunting_bags.stanfunctions
  #include functions/observation_models/hunting_composition.stanfunctions
  #include functions/observation_models/bycatch_composition.stanfunctions
  #include functions/observation_models/pregnancy.stanfunctions
  #include functions/observation_models/reproductive_signs.stanfunctions
  #include functions/pregnancy_births.stanfunctions
  #include functions/deaths.stanfunctions
  #include functions/aging.stanfunctions
  #include functions/hunting.stanfunctions
  #include functions/fates.stanfunctions
  #include functions/run_process.stanfunctions
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
  real<lower=0, upper=1> pregnancy_exposure_scaled;
  real<lower=0, upper=1> t_birth_to_start_hunt;
  real<lower=0, upper=1> hunting_duration; // shared hunting duration in years

  // ODE solver controls used by the joint competing-risks fate ODE.
  real<lower=1e-15> rel_tol;
  real<lower=1e-15> abs_tol;
  int<lower=1> max_num_steps;

  // Maximum accepted numerical deviation from the probability simplex
  // returned by the joint fate ODE. Larger deviations cause rejection.
  real<lower=1e-15> fate_probability_tolerance;

  // Observation-model coefficient of variation for reported hunting bags.
  real<lower=0> hunting_bag_cv;

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
  real<lower=0, upper=1> birth_rate_baseline_max;
  real<lower=0, upper=1> birth_rate_baseline_min_max_ratio;
  real herring_birth_rate_midpoint; // -alpha / beta in weighted herring-index units
  real herring_slope;
  real<lower=0, upper=1> herring_weight;
  real<lower=0, upper=1> density_dependence_scaled;

  // Reproductive-sign reporting parameters.
  real<lower=0, upper=1> report_ca_mean;
  real<lower=0, upper=1> report_placental_mean;
  real<lower=0, upper=1> ca_probability_without_birth;
  real<lower=0> report_placental_sd;
  real<lower=0> report_ca_sd;

  // Derived fitted demographic parameters carried into the forecast.
  real density_dependence_intercept;
  vector<lower=0>[n_demo_groups] non_hunting_mortality_rate;
  vector[n_demo_groups] survival_probability;
  vector[n_demo_groups] hunting_selectivity_finland;
  vector[n_demo_groups] hunting_selectivity_sweden;
  vector[n_demo_groups] bycatch_selectivity;

  real<lower=0, upper=1> birth_rate_at_carrying_capacity;
  real<lower=0, upper=1> birth_rate_baseline_min;
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
    birth_rate_baseline_future = compute_baseline_birth_rate(
      birth_rate_baseline_min,
      birth_rate_baseline_max,
      herring_birth_rate_midpoint,
      herring_slope,
      herring_weight,
      future_herring_index_baltic_proper_gulf_finland,
      future_herring_index_gulf_bothnia
    );

  // ------------------------------
  // FUTURE PROCESS INNOVATIONS
  // ------------------------------
  vector<lower=0>[n_future_years] hunting_effort_noise_sweden_future;
  vector<lower=0>[n_future_years] hunting_effort_noise_finland_future;
  vector<lower=0>[n_future_years] placental_scar_detection_noise_future;
  vector<lower=0>[n_future_years] ca_detection_noise_future;

  // ------------------------------
  // FUTURE STATE-PROCESS OUTPUTS
  // ------------------------------
  vector<lower=0, upper=1>[n_future_years] birth_rate_future;
  vector<lower=0, upper=1>[n_future_years] pregnancy_rate_future;
  vector<lower=0>[n_future_years] population_total_future;
  // Aerial surveys use ages 1+ of both sexes, matching the fitted model.
  vector<lower=0>[n_future_years] non_pup_population_total_future;
  matrix<lower=0>[n_demo_groups, n_future_years] population_comp_future;
  matrix<lower=0>[n_demo_groups, n_future_years] survivors_future;
  matrix<lower=0>[n_demo_groups, n_future_years] deaths_or_bycatch_future;
  matrix<lower=0>[n_demo_groups, n_future_years] hunted_sweden_future;
  matrix<lower=0>[n_demo_groups, n_future_years] hunted_finland_future;
  // Expected natural deaths plus bycatch, not bycatch counts alone
  matrix<lower=0>[n_demo_groups, n_future_years] non_hunting_deaths_expected_future;
  vector<lower=0>[n_future_years] hunting_bag_total_sweden_future;
  vector<lower=0>[n_future_years] hunting_bag_total_finland_future;
  vector<lower=0>[n_future_years] hunted_total_future;
  matrix<lower=0, upper=1>[4, n_future_years] reproductive_probs_future;

  for (year in 1:n_future_years) {
    hunting_effort_noise_sweden_future[year] = abs(std_normal_rng());
    hunting_effort_noise_finland_future[year] = abs(std_normal_rng());
    placental_scar_detection_noise_future[year] = abs(std_normal_rng());
    ca_detection_noise_future[year] = abs(std_normal_rng());
  }

  // Annual reporting probabilities for reproductive signs.
  vector<lower=0, upper=1>[n_future_years] placental_scar_detection_probability_future =
    report_placental_mean
    * exp(-placental_scar_detection_noise_future * report_placental_sd);

  vector<lower=0, upper=1>[n_future_years] ca_detection_probability_future =
    report_ca_mean
    * exp(-ca_detection_noise_future * report_ca_sd);

  // The first future birth rate is updated from the final fitted population.
  real<lower=0, upper=1> birth_rate_future_first = update_birth_rate(
    birth_rate_baseline_future[1],
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
    non_hunting_deaths_expected_future,
    hunting_bag_total_sweden_future,
    hunting_bag_total_finland_future,
    hunted_total_future,
    reproductive_probs_future
  ) = run_state_process_from_survivors_rng(
    n_future_years,
    n_age_classes,
    survivors_final,
    birth_rate_future_first,
    birth_rate_baseline_future,
    density_dependence_intercept,
    density_dependence_slope,
    non_hunting_mortality_rate,
    hunting_selectivity_sweden,
    hunting_selectivity_finland,
    future_hunting_quota_sweden,
    future_hunting_quota_finland,
    hunting_effort_sd_sweden,
    hunting_effort_sd_finland,
    hunting_effort_noise_sweden_future,
    hunting_effort_noise_finland_future,
    pregnancy_exposure_scaled,
    t_birth_to_start_hunt,
    hunting_duration,
    placental_scar_detection_probability_future,
    ca_detection_probability_future,
    ca_probability_without_birth,
    rel_tol,
    abs_tol,
    max_num_steps,
    fate_probability_tolerance
  );

  // ------------------------------
  // POSTERIOR PREDICTIVE OBSERVATIONS
  // ------------------------------
  for (year in 1:n_future_years) {
    non_pup_population_total_future[year] =
      sum(population_comp_future[2:n_age_classes, year])
      + sum(population_comp_future[(n_age_classes + 2):(2 * n_age_classes), year]);
  }

  array[n_future_years] int future_aerial_count = aerial_count_rng(
    future_year,
    non_pup_population_total_future,
    aerial_count_mu,
    aerial_count_overdispersion
  );

  array[n_future_years] real future_hunting_bags_sweden =
    hunting_bags_rng(
      future_year,
      hunting_bag_total_sweden_future,
      hunting_bag_cv
    );

  array[n_future_years] real future_hunting_bags_finland =
    hunting_bags_rng(
      future_year,
      hunting_bag_total_finland_future,
      hunting_bag_cv
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
      non_hunting_deaths_expected_future,
      bycatch_selectivity,
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
    log_lik_future_hunting_bags_sweden =
      rep_vector(0.0, n_future_hunting_bag_sweden);
  vector[n_future_hunting_bag_finland]
    log_lik_future_hunting_bags_finland =
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
      non_pup_population_total_future,
      aerial_count_mu,
      aerial_count_overdispersion
    );
  }

  if (n_future_hunting_bag_sweden > 0) {
    log_lik_future_hunting_bags_sweden =
      hunting_bags_pointwise_log_lik(
        future_obs_hunting_bag_sweden,
        hunting_bag_year_sweden,
        hunting_bag_total_sweden_future,
        hunting_bag_cv
      );
  }

  if (n_future_hunting_bag_finland > 0) {
    log_lik_future_hunting_bags_finland =
      hunting_bags_pointwise_log_lik(
        future_obs_hunting_bag_finland,
        hunting_bag_year_finland,
        hunting_bag_total_finland_future,
        hunting_bag_cv
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
      non_hunting_deaths_expected_future,
      bycatch_selectivity
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
    + sum(log_lik_future_hunting_bags_sweden)
    + sum(log_lik_future_hunting_bags_finland)
    + sum(log_lik_future_hunting_comp_sweden)
    + sum(log_lik_future_hunting_comp_finland)
    + sum(log_lik_future_bycatch)
    + sum(log_lik_future_pregnancy)
    + sum(log_lik_future_reproductive_signs);
}
