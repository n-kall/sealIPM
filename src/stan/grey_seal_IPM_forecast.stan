functions {

  #include stanfunctions/ode.stanfunctions
  #include stanfunctions/logitnormal.stanfunctions
  #include stanfunctions/likelihoods.stanfunctions
  #include stanfunctions/statespace.stanfunctions

}

data {

  int<lower=1> n_future_years;
  int<lower=1> n_age_classes;

  vector[n_future_years + 1] future_herring_index_baltic_proper_gulf_finland;
  vector[n_future_years + 1] future_herring_index_gulf_bothnia;
  array[n_future_years] int future_hunting_quota_sweden;
  array[n_future_years] int future_hunting_quota_finland;

  // Hunting/reproductive timing
  real<lower=0> t_mate_to_preg; // time from mating to pregnancy
  real<lower=0> t_birth_to_start_hunt; // time from birth to start of hunting period
  real<lower=0> t_birth_to_end_hunt; // time from birth to end of hunting period
  real<lower=0> t_hunt; // length of hunting period

  // ODE solver settings
  real<lower=0> rel_tol;
  real<lower=0> abs_tol;
  int<lower=1> max_num_steps;

  real<lower=0> harvest_bag_cv;

  array[n_future_years] int future_hunting_sample_size;
  array[n_future_years] int future_bycatch_sample_size;
  array[n_future_years] int future_reproductive_signs_sample_size;
  array[n_future_years] int future_pregnancy_sample_size;

}

transformed data {

  array[n_future_years] int future_year;
  for (y in 1:n_future_years) {
    future_year[y] = y;
  }

  int n_demo_groups = n_age_classes * 2;

  // aging matrix
  matrix[n_demo_groups, n_demo_groups] aging_matrix = create_aging_matrix(n_demo_groups, n_age_classes);

  vector[1] ode_init_state;
  array[1] real ode_times;
  ode_init_state[1] = 0.0;
  ode_times[1] = t_birth_to_end_hunt;

}

parameters {

  real aerial_count_mu;
  real aerial_count_overdispersion;

  real<lower=0> hunting_effort_sd_sweden;
  real<lower=0> hunting_effort_sd_finland;

  real<lower=0, upper=1> max_baseline_birth_rate;
  real<lower=0, upper=1> min_baseline_birth_rate;

  real herring_intercept_scaled;
  real herring_slope;

  real<lower=0, upper=1> herring_weight;
  real<lower=0, upper=1> density_dependence_scaled;

  real<lower=0, upper=1> report_ca_mean;
  real<lower=0, upper=1> report_placental_mean;
  real<lower=0, upper=1> prob_of_ca;
  real<lower=0> report_placental_sd;
  real<lower=0> report_ca_sd;

  real density_dependence_intercept;

  vector<lower=0>[n_demo_groups] mu_m;
  vector[n_demo_groups] S_diag;

  vector[n_demo_groups] hunting_selectivity_finland;
  vector[n_demo_groups] hunting_selectivity_sweden;
  vector[n_demo_groups] bycatch_bias;

  real<lower=0, upper=0.8> birth_rate_at_carrying_capacity;
  real<lower=0, upper=1> min_baseline_birth_rate_actual;

  real density_dependence_slope;

  vector[n_demo_groups] population_comp_final;
  real population_total_final;
  real<lower=0, upper=1> birth_rate_final;
  real<lower=0, upper=1> pregnancy_rate_final;

}

generated quantities {

  vector<lower=0, upper=1>[n_future_years + 1] baseline_birth_rate_future = compute_baseline_birth_rate(
    min_baseline_birth_rate_actual,
    max_baseline_birth_rate,
    herring_intercept_scaled,
    herring_slope,
    herring_weight,
    future_herring_index_baltic_proper_gulf_finland,
    future_herring_index_gulf_bothnia
  );
  

  tuple(
    vector[n_future_years],  // 1 birth_rate
    vector[n_future_years],  // 2 pregnancy_rate
    vector[n_future_years],  // 3 population_total
    matrix[n_demo_groups, n_future_years],  // 4 population_comp
    matrix[n_demo_groups, n_future_years],  // 5 survivors
    matrix[n_demo_groups, n_future_years],  // 6 deaths_or_bycatch
    matrix[n_demo_groups, n_future_years],  // 7 hunted_sweden
    matrix[n_demo_groups, n_future_years],  // 8 hunted_finland
    matrix[n_demo_groups, n_future_years],  // 9 bycatch_expected
    vector[n_future_years],  // 10 hunting_bag_total_sweden
    vector[n_future_years],  // 11 hunting_bag_total_finland
    vector[n_future_years],  // 12 hunted_total
    matrix[4, n_future_years]   // 13 reproductive_probs
  ) state_process_future;

  vector<lower=0>[n_future_years] epsilon_h_sw_future;
  vector<lower=0>[n_future_years] epsilon_h_fi_future;
  vector[n_future_years] epsilon_birth_future;
  vector[n_future_years] epsilon_sex_future;
  vector<lower=0>[n_future_years] epsilon_placental_future;
  vector<lower=0>[n_future_years] epsilon_ca_future;
  matrix[3 * n_demo_groups, n_future_years] transition_noise_raw_future;

  vector<lower=0, upper=1>[n_future_years] birth_rate_future;
  vector<lower=0, upper=1>[n_future_years] pregnancy_rate_future;

  vector[n_future_years] population_total_future;

  matrix[n_demo_groups, n_future_years] population_comp_future;

  matrix[n_demo_groups, n_future_years] survivors_future;
  matrix[n_demo_groups, n_future_years] deaths_or_bycatch_future;
  matrix[n_demo_groups, n_future_years] hunted_sweden_future;
  matrix[n_demo_groups, n_future_years] hunted_finland_future;
  matrix[n_demo_groups, n_future_years] bycatch_expected_future;

  vector[n_future_years] hunting_bag_total_sweden_future;
  vector[n_future_years] hunting_bag_total_finland_future;


  vector<lower=0>[n_future_years] hunted_total_future;

  matrix[4, n_future_years] reproductive_probs_future;

    

  for (i in 1:n_future_years) {
    epsilon_h_sw_future[i] = abs(normal_rng(0, 1));
    epsilon_h_fi_future[i] = abs(normal_rng(0, 1));

    epsilon_birth_future[i] = normal_rng(0, 1);
    epsilon_sex_future[i] = normal_rng(0, 1);

    epsilon_placental_future[i] = abs(normal_rng(0, 1));
    epsilon_ca_future[i] = abs(normal_rng(0, 1));

    for (k in 1:(3 * n_demo_groups)) {
      transition_noise_raw_future[k, i] = normal_rng(0, 1);
    }
  }

  vector<lower=0, upper=1>[n_future_years] pi_s_future =
  report_placental_mean * exp(-epsilon_placental_future * report_placental_sd);
  vector<lower=0, upper=1>[n_future_years] pi_c_future =
  report_ca_mean * exp(-epsilon_ca_future * report_ca_sd);

  
  state_process_future =
  run_state_process_from_first_population(
    n_future_years,
    n_age_classes,
    population_comp_final,
    birth_rate_final,
    population_total_final,
    baseline_birth_rate_future,
    density_dependence_intercept,
    density_dependence_slope,
    aging_matrix,
    S_diag,
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
    epsilon_birth_future,
    epsilon_sex_future,
    transition_noise_raw_future,
    pi_s_future,
    pi_c_future,
    prob_of_ca,
    ode_init_state,
    ode_times,
    rel_tol,
    abs_tol,
    max_num_steps
  );

  birth_rate_future =
  state_process_future.1;

  pregnancy_rate_future =
  state_process_future.2;

  population_total_future =
  state_process_future.3;

  population_comp_future =
  state_process_future.4;

  survivors_future =
  state_process_future.5;

  deaths_or_bycatch_future =
  state_process_future.6;

  hunted_sweden_future =
  state_process_future.7;

  hunted_finland_future =
  state_process_future.8;

  bycatch_expected_future =
  state_process_future.9;

  hunting_bag_total_sweden_future =
  state_process_future.10;

  hunting_bag_total_finland_future =
  state_process_future.11;

  hunted_total_future =
  state_process_future.12;

  reproductive_probs_future =
  state_process_future.13;

  
  array[n_future_years] int future_aerial_count = aerial_count_rng(
    future_year,
    population_total_future,
    aerial_count_mu,
    aerial_count_overdispersion
  );


  array[n_future_years] real future_harvest_bags_sweden =
  harvest_bags_rng(future_year, hunting_bag_total_sweden_future, harvest_bag_cv);

  array[n_future_years] real future_harvest_bags_finland =
  harvest_bags_rng(future_year, hunting_bag_total_finland_future, harvest_bag_cv);


  array[n_future_years, n_demo_groups] int future_hunting_comp_finland =
  hunting_comp_rng(
    future_year,
    hunted_finland_future,
    hunting_bag_total_finland_future,
    future_hunting_sample_size
  );

  array[n_future_years, n_demo_groups] int future_hunting_comp_sweden =
  hunting_comp_rng(
    future_year,
    hunted_sweden_future,
    hunting_bag_total_sweden_future,
    future_hunting_sample_size
  );

  array[n_future_years, n_demo_groups] int future_bycatch_comp =
  bycatch_comp_rng(
    future_year,
    bycatch_expected_future,
    bycatch_bias,
    future_bycatch_sample_size
  );

  array[n_future_years] int future_pregnancy =
  pregnancy_rng(
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
  
}
