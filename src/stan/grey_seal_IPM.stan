functions {
  #include stanfunctions/ode.stanfunctions
  #include stanfunctions/logitnormal.stanfunctions
  #include stanfunctions/likelihoods.stanfunctions
  #include stanfunctions/statespace.stanfunctions
}

data {

  int<lower=1> n_state_years; // number of modeled years for the state process
  int<lower=1> n_age_classes; // number of age classes per sex

  vector[2 * n_age_classes] population_init; // initial population structure

  // Environmental covariates
  vector[n_state_years + 1] herring_index_baltic_proper_gulf_finland; // herring weight-at-age for Baltic Proper + Gulf of Finland for previous year
  vector[n_state_years + 1] herring_index_gulf_bothnia; // herring weight-at-age for Gulf of Bothnia for previous year

  // Hunting quotas
  array[n_state_years] int<lower=0> hunting_quota_sweden; // Sweden's hunting quota
  array[n_state_years] int<lower=0> hunting_quota_finland; // Finland's hunting quota


  // Hunting/reproductive timing
  real<lower=0> t_mate_to_preg; // time from mating to pregnancy
  real<lower=0> t_birth_to_start_hunt; // time from birth to start of hunting period
  real<lower=0> t_birth_to_end_hunt; // time from birth to end of hunting period
  real<lower=0> t_hunt; // length of hunting period


  // ----------------------------
  // OBSERVATIONS
  // ----------------------------
  
  // Aerial survey observations
  int<lower=0, upper=n_state_years> n_aerial_years; // total number of surveyed years
  array[n_aerial_years] int<lower=1, upper=n_state_years> aerial_year; // which year was the count (1 is first year)
  array[n_aerial_years] int<lower=0> obs_aerial_count; // the counts for each observed year

  // Hunting bag observations
  int<lower=0, upper=n_state_years> n_hunting_bag_years_sweden; // total number of observed years
  array[n_hunting_bag_years_sweden] int<lower=1> hunting_bag_year_sweden; // which year is it from (1 is first year)
  array[n_hunting_bag_years_sweden] real<lower=0> obs_hunting_bag_sweden; // reported number of hunted seals

  int<lower=0, upper=n_state_years> n_hunting_bag_years_finland; // total number of observed years
  array[n_hunting_bag_years_finland] int<lower=1> hunting_bag_year_finland; // which year is it from (1 is first year)
  array[n_hunting_bag_years_finland] real<lower=0> obs_hunting_bag_finland; // reported number of hunted seals

  // Hunting composition observations, year-major
  int<lower=0, upper=n_state_years> n_hunting_comp_years_sweden; // total number of observed years
  array[n_hunting_comp_years_sweden] int<lower=1> hunting_comp_year_sweden; // which year is it from (1 is first year)
  array[n_hunting_comp_years_sweden, 2 * n_age_classes] int<lower=0> obs_hunting_comp_sweden; // reported hunting compositions

  int<lower=0, upper=n_state_years> n_hunting_comp_years_finland; // total number of observed years
  array[n_hunting_comp_years_finland] int<lower=1> hunting_comp_year_finland; // which year is it from (1 is first year)
  array[n_hunting_comp_years_finland, 2 * n_age_classes] int<lower=0> obs_hunting_comp_finland; // reported hunting compositions

  // Bycatch comp observations
  int<lower=0, upper=n_state_years> n_bycatch_years;
  array[n_bycatch_years] int<lower=1> bycatch_comp_year; // which year is it from (1 is first year)
  array[n_bycatch_years, 2 * n_age_classes] int<lower=0> obs_bycatch_comp; // reported bycatch compositions


  // Pregnancy observations
  int<lower=0, upper=n_state_years> n_pregnancy_years; // total number of surveyed years
  array[n_pregnancy_years] int<lower=1> pregnancy_count_year; // which year was it (1 is first year)
  array[n_pregnancy_years] int<lower=0> obs_pregnancy_count; // observed count of pregnancy
  array[n_pregnancy_years] int<lower=0> pregnancy_sample_size; // total checked for pregnancy

  // Reproductive signs observations
  int<lower=0, upper=n_state_years> n_reproductive_years; // total number of surveyed years
  array[n_reproductive_years] int<lower=1> reproductive_signs_year; // which year was it (1 is first year)
  array[n_reproductive_years, 4] int<lower=0> obs_reproductive_signs_finland; // observed reproductive signs


  // ----------------------------
  // LEAVE-FUTURE-OUT
  // ----------------------------

  // int<lower=0> n_future_aerial;
  // int<lower=0> n_future_hunting_bag_sweden;
  // int<lower=0> n_future_hunting_bag_finland;
  // int<lower=0> n_future_hunting_comp_sweden;
  // int<lower=0> n_future_hunting_comp_finland;
  // int<lower=0> n_future_bycatch;
  // int<lower=0> n_future_pregnancy;
  // int<lower=0> n_future_reproductive_signs;

  // // final-year observations

  // array[n_future_aerial]
  // int<lower=0> future_obs_aerial_count;

  // array[n_future_hunting_bag_sweden]
  // real<lower=0> future_obs_hunting_bag_sweden;

  // array[n_future_hunting_bag_finland]
  // real<lower=0> future_obs_hunting_bag_finland;

  // array[n_future_hunting_comp_sweden, 2 * n_age_classes]
  // int<lower=0> future_obs_hunting_comp_sweden;

  // array[n_future_hunting_comp_finland, 2 * n_age_classes]
  // int<lower=0> future_obs_hunting_comp_finland;

  // array[n_future_bycatch, 2 * n_age_classes]
  // int<lower=0> future_obs_bycatch_comp;

  // array[n_future_pregnancy]
  // int<lower=0> future_obs_pregnancy_count;

  // array[n_future_pregnancy]
  // int<lower=0> future_pregnancy_sample_size;

  // array[n_future_reproductive_signs, 4]
  // int<lower=0> future_obs_reproductive_signs_finland;

  // ----------------------------
  // PRIOR HYPERPARAMETERS
  // ----------------------------

  // Prior covariance/Cholesky structure for selectivity biases
  matrix[2 * n_age_classes, 2 * n_age_classes] bias_cholesky_factor;

  // Initial population size: lognormal(log_mean, log_sd)
  real prior_initial_population_log_mean;
  real<lower=0> prior_initial_population_log_sd;

  // Carrying capacity: lognormal(log_mean, log_sd)
  real prior_carrying_capacity_log_mean;
  real<lower=0> prior_carrying_capacity_log_sd;

  // Aerial detection probability: beta(alpha, beta)
  real<lower=0> prior_aerial_detection_alpha;
  real<lower=0> prior_aerial_detection_beta;

  // Aerial overdispersion: lognormal(log_mean, log_sd)
  real prior_aerial_overdispersion_log_mean;
  real<lower=0> prior_aerial_overdispersion_log_sd;

  // Log-selectivity random effects: normal(mean, sd)
  real<lower=0> prior_hunting_selectivity_sd;
  real<lower=0> prior_bycatch_selectivity_sd;

  // Hunting effort SD: cauchy(location, scale), constrained positive
  real<lower=0> prior_hunting_effort_sd_location;
  real<lower=0> prior_hunting_effort_sd_scale;

  // Herring birth-rate regression
  real<lower=0> prior_herring_intercept_scaled_sd;
  real<lower=0> prior_herring_slope_sd;

  // Male mortality log offsets: cauchy(location, scale)
  real prior_male_pup_mortality_offset_location;
  real<lower=0> prior_male_pup_mortality_offset_scale;

  real prior_male_adult_mortality_offset_location;
  real<lower=0> prior_male_adult_mortality_offset_scale;

  // Reproductive-sign process SDs: normal(location, scale), constrained positive
  real prior_reproductive_process_sd_location;
  real<lower=0> prior_reproductive_process_sd_scale;

  // Harvest bag observation coefficient of variation
  real<lower=0> harvest_bag_cv;

  // ----------------------------
  // ALGORITHM SETTINGS
  // ----------------------------
  int<lower=1> population_burn_in; // number of iterations for initializing the population

  // ODE solver settings
  real<lower=0> rel_tol;
  real<lower=0> abs_tol;
  int<lower=1> max_num_steps;

}

transformed data {

  int n_demo_groups = 2 * n_age_classes;

  // aging matrix
  matrix[n_demo_groups, n_demo_groups] aging_matrix = create_aging_matrix(n_demo_groups, n_age_classes);

  // ode
  vector[1] ode_init_state;
  array[1] real ode_times;
  ode_init_state[1] = 0.0;
  ode_times[1] = t_birth_to_end_hunt;

}

parameters {
  // observations
  real<lower=0, upper=1> aerial_count_mu; // mu (maybe need to account for different probability of haul-out for pups and adults)
  real<lower=0> aerial_count_overdispersion; // r

  // Biases
  vector[n_demo_groups] hunting_selectivity_sweden_sc; // g_sw
  vector[n_demo_groups] hunting_selectivity_finland_sc; // g_fi
  vector[n_demo_groups] bycatch_selectivity_sc; // g_bc

  real<lower=0> hunting_effort_sd_sweden; // sigma_sw
  real<lower=0> hunting_effort_sd_finland; // sigma_fi
  vector<lower=0>[n_state_years] epsilon_h_sw; //
  vector<lower=0>[n_state_years] epsilon_h_fi;

  real<lower=0, upper=1> survival_shape; // c


  //Natural mortality
  real<lower=0, upper=1> phi_a_sc;
  real<lower=0, upper=1> phi_sc;
  real male_pup_survival_offset; // nu_0
  real male_adult_survival_offset; // nu_5+

  //Birth rate
  real<lower=0> carrying_capacity; // K_max
  //real log_density_dependence_slope; 

  real<lower=0, upper=1> max_baseline_birth_rate; // b0_max
  real<lower=0, upper=1> min_baseline_birth_rate; // b0_min / b0_max


  real herring_intercept_scaled; // alpha_sc
  real herring_slope;            // beta

  real<lower=0, upper=1> herring_weight; // w


  real<lower=0, upper=1> density_dependence_scaled; // theta0


  // state process
  real population_init_size; // n_0
  vector[n_state_years] epsilon_birth; //
  vector[n_state_years] epsilon_sex;
  matrix[3 * n_demo_groups, n_state_years] transition_noise_raw;

  // reproductive signs

  real<lower=0, upper=1> report_ca_mean;
  real<lower=0, upper=1> report_placental_mean;
  real<lower=0, upper=1> prob_of_ca;
  real<lower=0> report_placental_sd;
  real<lower=0> report_ca_sd;
  vector<lower=0>[n_state_years] epsilon_ca;
  vector<lower=0>[n_state_years] epsilon_placental;

}

transformed parameters {

  // ----------------------------
  // TIME-INVARIANT PARAMETERS
  // ----------------------------

  real<lower=0, upper=1> phi_a = 0.8 + 0.2 * phi_a_sc;
  real<lower=0, upper=1> phi_pup = phi_a * phi_sc;
  
  // density dependence
  real density_dependence_intercept = compute_density_dependence_intercept(
    max_baseline_birth_rate,
    density_dependence_scaled
  );

  // natural mortality
  vector<lower=0>[n_demo_groups] mu_m = mortality_rates(
    phi_pup,
    phi_a,
    survival_shape,
    n_age_classes,
    male_pup_survival_offset,
    male_adult_survival_offset
  );
  vector[n_demo_groups] S_diag = exp(-mu_m);

  // hunting biases

  vector[n_demo_groups] hunting_selectivity_finland = bias_cholesky_factor * hunting_selectivity_finland_sc;
  vector[n_demo_groups] hunting_selectivity_sweden = bias_cholesky_factor * hunting_selectivity_sweden_sc;

  // bycatch bias
  vector[n_demo_groups] bycatch_bias = bias_cholesky_factor * bycatch_selectivity_sc;

  // birth rate at carrying capacity
  real<lower=0, upper=0.8> birth_rate_at_carrying_capacity;
  real<lower=0, upper=1> min_baseline_birth_rate_actual;

  min_baseline_birth_rate_actual =
  max_baseline_birth_rate * min_baseline_birth_rate;

  real<lower=0, upper=1> reference_baseline_birth_rate =
  min_baseline_birth_rate_actual
  + (max_baseline_birth_rate - min_baseline_birth_rate_actual)
  * inv_logit(herring_intercept_scaled * herring_slope);


  birth_rate_at_carrying_capacity =
  2.0 * (1.0 - phi_a) /
  exp(sum(-mu_m[1:(n_age_classes - 1)]));


  //real density_dependence_slope = exp(log_density_dependence_slope);
  
   real density_dependence_slope =
   log(
     1.0
     - log(birth_rate_at_carrying_capacity / reference_baseline_birth_rate)
       / density_dependence_intercept
   ) / carrying_capacity;

  vector<lower=0, upper=1>[n_state_years] pi_s =
  report_placental_mean * exp(-epsilon_placental * report_placental_sd);
  vector<lower=0, upper=1>[n_state_years] pi_c =
  report_ca_mean * exp(-epsilon_ca * report_ca_sd);

  // ----------------------------
  // YEAR-TO-YEAR STATE TRANSITION
  // ----------------------------

  // precompute baseline birth rates for each year
  // baseline birth rates just depend on herring
  vector<lower=0, upper=1>[n_state_years + 1] baseline_birth_rate = compute_baseline_birth_rate(
    min_baseline_birth_rate_actual,
    max_baseline_birth_rate,
    herring_intercept_scaled,
    herring_slope,
    herring_weight,
    herring_index_baltic_proper_gulf_finland,
    herring_index_gulf_bothnia
  );

  // ----------------------------
  // STATE PROCESS OUTPUTS
  // ----------------------------

  vector<lower=0, upper=1>[n_state_years] birth_rate;
  vector<lower=0, upper=1>[n_state_years] pregnancy_rate;

  vector[n_state_years] population_total;

  matrix[n_demo_groups, n_state_years] population_comp;

  matrix[n_demo_groups, n_state_years] survivors;
  matrix[n_demo_groups, n_state_years] deaths_or_bycatch;
  matrix[n_demo_groups, n_state_years] hunted_sweden;
  matrix[n_demo_groups, n_state_years] hunted_finland;
  matrix[n_demo_groups, n_state_years] bycatch_expected;

  vector[n_state_years] hunting_bag_total_sweden;
  vector[n_state_years] hunting_bag_total_finland;

  vector<lower=0>[n_state_years] hunted_total;

  matrix[4, n_state_years] reproductive_probs;


  // ----------------------------
  // STATE PROCESS TUPLES
  // ----------------------------

  tuple(vector[n_demo_groups], real, real) init_state;

  tuple(
    vector[n_state_years],  // 1 birth_rate
    vector[n_state_years],  // 2 pregnancy_rate
    vector[n_state_years],  // 3 population_total
    matrix[n_demo_groups, n_state_years],  // 4 population_comp
    matrix[n_demo_groups, n_state_years],  // 5 survivors
    matrix[n_demo_groups, n_state_years],  // 6 deaths_or_bycatch
    matrix[n_demo_groups, n_state_years],  // 7 hunted_sweden
    matrix[n_demo_groups, n_state_years],  // 8 hunted_finland
    matrix[n_demo_groups, n_state_years],  // 9 bycatch_expected
    vector[n_state_years],  // 10 hunting_bag_total_sweden
    vector[n_state_years],  // 11 hunting_bag_total_finland
    vector[n_state_years],  // 12 hunted_total
    matrix[4, n_state_years]   // 13 reproductive_probs
  ) state_process;

  

real initial_birth_rate;

initial_birth_rate =
  update_birth_rate(
    baseline_birth_rate[1],
    density_dependence_intercept,
    density_dependence_slope,
    sum(population_init)
  );

if (
  is_nan(initial_birth_rate) ||
  is_inf(initial_birth_rate) ||
  initial_birth_rate <= 0 ||
  initial_birth_rate >= 1
) {
  reject(
    "Bad initial_birth_rate = ", initial_birth_rate,
    ", baseline_birth_rate[1] = ", baseline_birth_rate[1],
    ", density_dependence_intercept = ", density_dependence_intercept,
    ", density_dependence_slope = ", density_dependence_slope,
    ", sum(population_init) = ", sum(population_init)
  );
}

  init_state =
initialize_population_with_burnin(
  population_init,
  population_init_size,
  population_burn_in,
  initial_birth_rate,
  aging_matrix,
  S_diag,
  n_age_classes
);
  
  state_process =
  run_state_process_from_first_population(
    n_state_years,
    n_age_classes,
    init_state.1,
    init_state.2,
    init_state.3,
    baseline_birth_rate,
    density_dependence_intercept,
    density_dependence_slope,
    aging_matrix,
    S_diag,
    mu_m,
    hunting_selectivity_sweden,
    hunting_selectivity_finland,
    hunting_quota_sweden,
    hunting_quota_finland,
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
    ode_init_state,
    ode_times,
    rel_tol,
    abs_tol,
    max_num_steps
  );

  birth_rate =
  state_process.1;

  pregnancy_rate =
  state_process.2;

  population_total =
  state_process.3;

  population_comp =
  state_process.4;

  survivors =
  state_process.5;

  deaths_or_bycatch =
  state_process.6;

  hunted_sweden =
  state_process.7;

  hunted_finland =
  state_process.8;

  bycatch_expected =
  state_process.9;

  hunting_bag_total_sweden =
  state_process.10;

  hunting_bag_total_finland =
  state_process.11;

  hunted_total =
  state_process.12;

  reproductive_probs =
  state_process.13;
  
}

model {

  // ----------------------------
  // PRIORS
  // ----------------------------

  //Initial population size
  population_init_size ~ lognormal(prior_initial_population_log_mean, prior_initial_population_log_sd);

  // Natural mortality
  // phi_a_sc ~ uniform(0,1); // implied prior by constraints
  // phi_sc ~ uniform(0,1); // implied prior by constraints

  male_pup_survival_offset ~ cauchy(prior_male_pup_mortality_offset_location, prior_male_pup_mortality_offset_scale);
  male_adult_survival_offset ~ cauchy(prior_male_adult_mortality_offset_location, prior_male_adult_mortality_offset_scale);

  // Carrying capactiy
   carrying_capacity ~ lognormal(prior_carrying_capacity_log_mean, prior_carrying_capacity_log_sd);
  //log_density_dependence_slope ~ normal(log(1e-5), 1);
 
  
  // Hunting and bycatch bias

  hunting_selectivity_sweden_sc ~ normal(0, prior_hunting_selectivity_sd);
  hunting_selectivity_finland_sc ~ normal(0, prior_hunting_selectivity_sd);
  bycatch_selectivity_sc ~ normal(0, prior_bycatch_selectivity_sd);

  // Hunting effort sd
  hunting_effort_sd_sweden ~ cauchy(prior_hunting_effort_sd_location, prior_hunting_effort_sd_scale);
  hunting_effort_sd_finland ~ cauchy(prior_hunting_effort_sd_location, prior_hunting_effort_sd_scale);

  // Birth rate
  // b0max ~ uniform (0, 1); // implied prior by bounds
  // b0min_sc ~ uniform (0, 1); // implied prior by bounds

  herring_intercept_scaled ~ normal(0, prior_herring_intercept_scaled_sd);
  herring_slope ~ normal(0, prior_herring_slope_sd);
  // herring_weight ~ uniform(0, 1); // implied prior by bounds


  // Observation of aerial survey
  aerial_count_mu ~ beta(prior_aerial_detection_alpha, prior_aerial_detection_beta);
  aerial_count_overdispersion ~ lognormal (prior_aerial_overdispersion_log_mean, prior_aerial_overdispersion_log_sd);

  // Observation of reproductive signs
  // prob_ca_nonpreg ~ uniform(0, 1); // implied prior by constraints
  // pi_s_mean ~ uniform(0, 1); // implied prior by constraints
  // pi_c_mean ~ uniform(0, 1); // implied prior by constraints

  report_placental_sd ~ normal (0, 0.1);
  report_ca_sd ~ normal (0, 0.1);


  // standard normals for stochasticity
  epsilon_h_sw ~ std_normal();
  epsilon_h_fi ~ std_normal();

  epsilon_ca ~ std_normal();
  epsilon_placental ~ std_normal();

  // Stochasicity for birth process
  epsilon_birth ~ std_normal();
  epsilon_sex ~ std_normal();

  // Stochasticity for state transitions
  for (i in 1:n_state_years){
    transition_noise_raw[,i] ~ std_normal();
  }

  // ----------------------------
  // LIKELIHOODS
  // ----------------------------

  // Aerial surveys
  target += aerial_count_lpmf(
    obs_aerial_count |
    aerial_year,
    population_total,
    aerial_count_mu,
    aerial_count_overdispersion
  );

  // Harvest totals
  target += harvest_bags_lpdf(
    obs_hunting_bag_finland |
    hunting_bag_year_finland,
    hunting_bag_total_finland,
    harvest_bag_cv
  );
  target += harvest_bags_lpdf(
    obs_hunting_bag_sweden |
    hunting_bag_year_sweden,
    hunting_bag_total_sweden,
    harvest_bag_cv
  );

  // Hunting comp
  target += hunting_comp_lpmf(
    obs_hunting_comp_sweden |
    hunting_comp_year_sweden,
    hunted_sweden,
    hunting_bag_total_sweden
  );

  target += hunting_comp_lpmf(
    obs_hunting_comp_finland |
    hunting_comp_year_finland,
    hunted_finland,
    hunting_bag_total_finland
  );

  // Bycatch comp
  target += bycatch_comp_lpmf(
    obs_bycatch_comp |
    bycatch_comp_year,
    bycatch_expected,
    bycatch_bias
  );

  // Pregnancy
  target += pregnancy_lpmf(
    obs_pregnancy_count |
    pregnancy_count_year,
    pregnancy_sample_size,
    pregnancy_rate
  );

  // Reproductive signs
  target += reproductive_signs_lpmf(
    obs_reproductive_signs_finland |
    reproductive_signs_year,
    reproductive_probs
  );
}

generated quantities {

  // ----------------------------
  // IN-SAMPLE PREDICTIONS
  // ----------------------------

  // --------------------------------------------------
  // LEAVE-FUTURE-OUT LOG LIKELIHOODS
  // --------------------------------------------------

  // real log_lik_future_aerial = 0;
  // real log_lik_future_harvest_bags_sweden = 0;
  // real log_lik_future_harvest_bags_finland = 0;
  // real log_lik_future_hunting_comp_sweden = 0;
  // real log_lik_future_hunting_comp_finland = 0;
  // real log_lik_future_bycatch = 0;
  // real log_lik_future_pregnancy = 0;
  // real log_lik_future_reproductive_signs = 0;

  // if (n_future_aerial > 0) {
  //   log_lik_future_aerial =
  //     sum(
  //       aerial_count_pointwise_log_lik(
  //       future_obs_aerial_count,
  //       rep_array(n_state_years, 1),
  //       population_total,
  //       aerial_count_mu,
  //       aerial_count_overdispersion
  //     ));
  // }

  // if (n_future_hunting_bag_sweden > 0) {
  //   log_lik_future_harvest_bags_sweden =
  //     sum(
  //       harvest_bags_pointwise_log_lik(
  //       future_obs_hunting_bag_sweden,
  //       rep_array(n_state_years, 1),
  //       hunting_bag_total_sweden,
  //       harvest_bag_cv
  //     ));
  // }

  // if (n_future_hunting_bag_finland > 0) {
  //   log_lik_future_harvest_bags_finland =
  //     sum(
  //       harvest_bags_pointwise_log_lik(
  //       future_obs_hunting_bag_finland,
  //       rep_array(n_state_years, 1),
  //       hunting_bag_total_finland,
  //       harvest_bag_cv
  //       )
  //     );
  // }

  // if (n_future_hunting_comp_sweden > 0) {
  //   log_lik_future_hunting_comp_sweden =
  //     sum(
  //       hunting_comp_pointwise_log_lik(
  //       future_obs_hunting_comp_sweden,
  //       rep_array(n_state_years, 1),
  //       hunted_sweden,
  //       hunting_bag_total_sweden
  //       )
  //     );
  // }

  // if (n_future_hunting_comp_finland > 0) {
  //   log_lik_future_hunting_comp_finland =
  //     sum(
  //       hunting_comp_pointwise_log_lik(
  //       future_obs_hunting_comp_finland,
  //       rep_array(n_state_years, 1),
  //       hunted_finland,
  //       hunting_bag_total_finland
  //       )
  //     );
  // }

  // if (n_future_bycatch > 0) {
  //   log_lik_future_bycatch =
  //     sum(
  //       bycatch_comp_pointwise_log_lik(
  //       future_obs_bycatch_comp,
  //       rep_array(n_state_years, 1),
  //       bycatch_expected,
  //       bycatch_bias
  //       )
  //     );
  // }

  // if (n_future_pregnancy > 0) {
  //   log_lik_future_pregnancy =
  //     sum(
  //       pregnancy_pointwise_log_lik(
  //       future_obs_pregnancy_count,
  //       rep_array(n_state_years, 1),
  //       future_pregnancy_sample_size,
  //       pregnancy_rate
  //       )
  //     );
  // }
    
  // if (n_future_reproductive_signs > 0) {
  //   log_lik_future_reproductive_signs =
  //     sum(
  //       reproductive_signs_pointwise_log_lik(
  //       future_obs_reproductive_signs_finland,
  //       rep_array(n_state_years, 1),
  //       reproductive_probs
  //       )
  //     );
  // }

    
  // real log_lik_future_total =
  //     log_lik_future_aerial
  //   + log_lik_future_harvest_bags_sweden
  //   + log_lik_future_harvest_bags_finland
  //   + log_lik_future_hunting_comp_sweden
  //   + log_lik_future_hunting_comp_finland
  //   + log_lik_future_bycatch
  //   + log_lik_future_pregnancy
  //   + log_lik_future_reproductive_signs;
  
  // ----------------------------
  // POINTWISE LOG LIKELIHOODS
  // ----------------------------

  // Aerial surveys
  vector[n_aerial_years] log_lik_aerial;
  log_lik_aerial =
    aerial_count_pointwise_log_lik(
      obs_aerial_count,
      aerial_year,
      population_total,
      aerial_count_mu,
      aerial_count_overdispersion
    );

  // Hunting bag totals
  vector[n_hunting_bag_years_finland] log_lik_harvest_bags_finland;
  log_lik_harvest_bags_finland =
    harvest_bags_pointwise_log_lik(
      obs_hunting_bag_finland,
      hunting_bag_year_finland,
      hunting_bag_total_finland,
      harvest_bag_cv
    );

  vector[n_hunting_bag_years_sweden] log_lik_harvest_bags_sweden;
  log_lik_harvest_bags_sweden =
    harvest_bags_pointwise_log_lik(
      obs_hunting_bag_sweden,
      hunting_bag_year_sweden,
      hunting_bag_total_sweden,
      harvest_bag_cv
    );

  // Hunting composition
  vector[n_hunting_comp_years_finland] log_lik_hunting_comp_finland;
  log_lik_hunting_comp_finland =
    hunting_comp_pointwise_log_lik(
      obs_hunting_comp_finland,
      hunting_comp_year_finland,
      hunted_finland,
      hunting_bag_total_finland
    );

  vector[n_hunting_comp_years_sweden] log_lik_hunting_comp_sweden;
  log_lik_hunting_comp_sweden =
    hunting_comp_pointwise_log_lik(
      obs_hunting_comp_sweden,
      hunting_comp_year_sweden,
      hunted_sweden,
      hunting_bag_total_sweden
    );

  // Bycatch
  vector[n_bycatch_years] log_lik_bycatch;
  log_lik_bycatch =
    bycatch_comp_pointwise_log_lik(
      obs_bycatch_comp,
      bycatch_comp_year,
      bycatch_expected,
      bycatch_bias
    );

  // Pregnancy
  vector[n_pregnancy_years] log_lik_pregnancy;
  log_lik_pregnancy =
    pregnancy_pointwise_log_lik(
      obs_pregnancy_count,
      pregnancy_count_year,
      pregnancy_sample_size,
      pregnancy_rate
    );

  // Reproductive signs
  vector[n_reproductive_years] log_lik_reproductive_signs;
  log_lik_reproductive_signs =
    reproductive_signs_pointwise_log_lik(
      obs_reproductive_signs_finland,
      reproductive_signs_year,
      reproductive_probs
    );


  // Output for forecasting
  vector[n_demo_groups] population_comp_final = population_comp[, n_state_years];
  real population_total_final = population_total[n_state_years];
  real birth_rate_final = birth_rate[n_state_years];
  real pregnancy_rate_final = pregnancy_rate[n_state_years];
  
}
