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

  // --------------------------------
  // PROCESS AND ALGORITHM PARAMETERS
  // --------------------------------


  int<lower=1> n_state_years; // number of modeled years for the state process
  int<lower=1> n_age_classes; // number of age classes per sex

  vector[2 * n_age_classes] population_init; // initial population structure

  int<lower=1> population_burn_in; // number of iterations for initializing the population

  // Hunting/reproductive timing
  real<lower=0> pregnancy_exposure_scaled; // Scaled pregnancy exposure coefficient, not raw elapsed time
  real<lower=0> t_birth_to_start_hunt; // time from birth to start of hunting period
  real<lower=0> hunting_duration; // shared hunting duration in years

  // ODE solver settings
  real<lower=0> rel_tol;
  real<lower=0> abs_tol;
  int<lower=1> max_num_steps;
  real<lower=0> fate_probability_tolerance;
  int<lower=0, upper=1> prior_only;

  //---------------
  // COVARIATE DATA
  //---------------

  // Environmental covariates
  vector[n_state_years + 1] herring_index_baltic_proper_gulf_finland; // herring weight-at-age for Baltic Proper + Gulf of Finland for previous year
  vector[n_state_years + 1] herring_index_gulf_bothnia; // herring weight-at-age for Gulf of Bothnia for previous year

  // Hunting quotas
  array[n_state_years] int<lower=0> hunting_quota_sweden; // Sweden's hunting quota
  array[n_state_years] int<lower=0> hunting_quota_finland; // Finland's hunting quota

  // ----------------------------
  // OBSERVATIONS
  // ----------------------------

  // Aerial survey observations
  int<lower=0, upper=n_state_years> n_aerial_years; // total number of surveyed years
  array[n_aerial_years] int<lower=1, upper=n_state_years> aerial_year; // which year was the count (1 is first year)
  array[n_aerial_years] int<lower=0> obs_aerial_count; // the counts for each observed year

  // Hunting bag observations
  int<lower=0, upper=n_state_years> n_hunting_bag_years_sweden; // total number of observed years
  array[n_hunting_bag_years_sweden] int<lower=1, upper=n_state_years> hunting_bag_year_sweden; // which year is it from (1 is first year)
  array[n_hunting_bag_years_sweden] real<lower=0> obs_hunting_bag_sweden; // reported number of hunted seals

  int<lower=0, upper=n_state_years> n_hunting_bag_years_finland; // total number of observed years
  array[n_hunting_bag_years_finland] int<lower=1, upper=n_state_years> hunting_bag_year_finland; // which year is it from (1 is first year)
  array[n_hunting_bag_years_finland] real<lower=0> obs_hunting_bag_finland; // reported number of hunted seals

  // Hunting composition observations, year-major
  int<lower=0, upper=n_state_years> n_hunting_comp_years_sweden; // total number of observed years
  array[n_hunting_comp_years_sweden] int<lower=1, upper=n_state_years> hunting_comp_year_sweden; // which year is it from (1 is first year)
  array[n_hunting_comp_years_sweden, 2 * n_age_classes] int<lower=0> obs_hunting_comp_sweden; // reported hunting compositions

  int<lower=0, upper=n_state_years> n_hunting_comp_years_finland; // total number of observed years
  array[n_hunting_comp_years_finland] int<lower=1, upper=n_state_years> hunting_comp_year_finland; // which year is it from (1 is first year)
  array[n_hunting_comp_years_finland, 2 * n_age_classes] int<lower=0> obs_hunting_comp_finland; // reported hunting compositions

  // Bycatch comp observations
  int<lower=0, upper=n_state_years> n_bycatch_years;
  array[n_bycatch_years] int<lower=1, upper=n_state_years> bycatch_comp_year; // which year is it from (1 is first year)
  array[n_bycatch_years, 2 * n_age_classes] int<lower=0> obs_bycatch_comp; // reported bycatch compositions

  // Pregnancy observations
  int<lower=0, upper=n_state_years> n_pregnancy_years; // total number of surveyed years
  array[n_pregnancy_years] int<lower=1, upper=n_state_years> pregnancy_count_year; // which year was it (1 is first year)
  array[n_pregnancy_years] int<lower=0> obs_pregnancy_count; // observed count of pregnancy
  array[n_pregnancy_years] int<lower=0> pregnancy_sample_size; // total checked for pregnancy

  // Reproductive signs observations
  int<lower=0, upper=n_state_years> n_reproductive_years; // total number of surveyed years
  array[n_reproductive_years] int<lower=1, upper=n_state_years> reproductive_signs_year; // which year was it (1 is first year)
  array[n_reproductive_years, 4] int<lower=0> obs_reproductive_signs_finland; // observed reproductive signs


  // ----------------------------
  // PRIOR HYPERPARAMETERS
  // ----------------------------

  // Prior covariance/Cholesky structure for selectivity biases
  matrix[2 * n_age_classes, 2 * n_age_classes] selectivity_cholesky;

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

  // Hunting effort SD: t4(location, scale), constrained positive
  real<lower=0> prior_hunting_effort_sd_location;
  real<lower=0> prior_hunting_effort_sd_scale;

  // Herring birth-rate regression
  real<lower=0> prior_herring_birth_rate_midpoint_sd;
  real<lower=0> prior_herring_slope_sd;

  // Male log mortality offsets: normal(location, scale)
  real prior_male_pup_mortality_log_offset_location;
  real<lower=0> prior_male_pup_mortality_log_offset_scale;

  real prior_male_adult_mortality_log_offset_location;
  real<lower=0> prior_male_adult_mortality_log_offset_scale;

  // Hunting bag observation coefficient of variation
  real<lower=0> hunting_bag_cv;
}

transformed data {

  int n_demo_groups = 2 * n_age_classes;

  // aging matrix
  matrix[n_demo_groups, n_demo_groups] aging_matrix = create_aging_matrix(n_demo_groups, n_age_classes);


    // Composition predictions use the same sample size as the corresponding
  // observed year.
  array[n_hunting_comp_years_sweden] int
    hunting_sample_size_sweden_pred;
  array[n_hunting_comp_years_finland] int
    hunting_sample_size_finland_pred;
  array[n_bycatch_years] int
    bycatch_sample_size_pred;
  array[n_reproductive_years] int
    reproductive_signs_sample_size_pred;

  for (i in 1:n_hunting_comp_years_sweden) {
    hunting_sample_size_sweden_pred[i] =
      sum(obs_hunting_comp_sweden[i]);
  }

  for (i in 1:n_hunting_comp_years_finland) {
    hunting_sample_size_finland_pred[i] =
      sum(obs_hunting_comp_finland[i]);
  }

  for (i in 1:n_bycatch_years) {
    bycatch_sample_size_pred[i] =
      sum(obs_bycatch_comp[i]);
  }

  for (i in 1:n_reproductive_years) {
    reproductive_signs_sample_size_pred[i] =
      sum(obs_reproductive_signs_finland[i]);
  }

}

parameters {
  // aerial survey
  real<lower=0, upper=1> aerial_count_mu; // mu (maybe need to account for different probability of haul-out for pups and adults)
  real<lower=0> aerial_count_overdispersion; // r

  // selectivity
  vector[n_demo_groups] hunting_selectivity_sweden_sc; // g_sw
  vector[n_demo_groups] hunting_selectivity_finland_sc; // g_fi
  vector[n_demo_groups] bycatch_selectivity_sc; // g_bc

  // hunting effort
  real<lower=0> hunting_effort_sd_sweden; // sigma_sw
  real<lower=0> hunting_effort_sd_finland; // sigma_fi
  vector<lower=0>[n_state_years] hunting_effort_noise_sweden; //
  vector<lower=0>[n_state_years] hunting_effort_noise_finland;

  // survival curve
  real<lower=0, upper=1> mortality_age_shape; // c

  // natural mortality
  real<lower=0.8, upper=1> female_adult_survival_probability;
  real<lower=0, upper=1> pup_to_adult_survival_ratio;
  real male_pup_mortality_log_offset; // nu_0, positive values increase mortality
  real male_adult_mortality_log_offset; // nu_5+, positive values increase mortality

  // birth rate
  real<lower=0> carrying_capacity; // K at reference herring conditions
  real<lower=0, upper=1> birth_rate_baseline_max; // b0_max
  real<lower=0, upper=1> birth_rate_baseline_min_max_ratio; // Minimum baseline birth rate divided by maximum
  real<lower=0, upper=1> birth_rate_at_capacity_to_baseline_ratio; // Birth rate at carrying capacity divided by reference baseline birth rate
  real herring_birth_rate_midpoint; // -alpha / beta in weighted herring-index units
  real herring_slope;            // beta
  real<lower=0, upper=1> herring_weight; // w
  real<lower=0, upper=1> density_dependence_scaled; // Scaled parameter used to derive the density dependence intercept

  // state process
  real<lower=0> population_init_size; // n_0
  vector[n_state_years] birth_count_noise; //
  vector[n_state_years] pup_sex_allocation_noise;
  matrix[3 * n_demo_groups, n_state_years] transition_noise_raw;

  // reproductive signs

  real<lower=0, upper=1> report_ca_mean;
  real<lower=0, upper=1> report_placental_mean;
  real<lower=0, upper=1> ca_probability_without_birth;
  real<lower=0> report_placental_sd;
  real<lower=0> report_ca_sd;
  vector<lower=0>[n_state_years] ca_detection_noise;
  vector<lower=0>[n_state_years] placental_scar_detection_noise;

}

transformed parameters {

  // ----------------------------
  // TIME-INVARIANT PARAMETERS
  // ----------------------------

  real female_pup_survival_probability = female_adult_survival_probability * pup_to_adult_survival_ratio;

  // density dependence
  real density_dependence_intercept = compute_density_dependence_intercept(
    birth_rate_baseline_max,
    density_dependence_scaled
  );

  // natural mortality
  vector[n_demo_groups] non_hunting_mortality_rate = mortality_rates(
    female_pup_survival_probability,
    female_adult_survival_probability,
    mortality_age_shape,
    n_age_classes,
    male_pup_mortality_log_offset,
    male_adult_mortality_log_offset
  );

  vector[n_demo_groups] survival_probability = exp(-non_hunting_mortality_rate);

  // hunting selectivity
  vector[n_demo_groups] hunting_selectivity_finland = selectivity_cholesky * hunting_selectivity_finland_sc;
  vector[n_demo_groups] hunting_selectivity_sweden = selectivity_cholesky * hunting_selectivity_sweden_sc;

  // bycatch selectivity
  vector[n_demo_groups] bycatch_selectivity = selectivity_cholesky * bycatch_selectivity_sc;

  real birth_rate_baseline_min;

  birth_rate_baseline_min =
  birth_rate_baseline_max * birth_rate_baseline_min_max_ratio;

  // Baseline birth rate at reference herring conditions before density effects
  real birth_rate_baseline_reference =
  birth_rate_baseline_min
  + (birth_rate_baseline_max - birth_rate_baseline_min)
  * inv_logit(-herring_birth_rate_midpoint * herring_slope);

  // Birth rate required for demographic replacement under the survival schedule
  real birth_rate_replacement;
  // Density-adjusted birth rate at carrying capacity under reference herring conditions
  real birth_rate_at_carrying_capacity = birth_rate_baseline_reference * birth_rate_at_capacity_to_baseline_ratio;

  birth_rate_replacement = compute_birth_rate_replacement(
    non_hunting_mortality_rate,
    n_age_classes,
    female_adult_survival_probability,
    0
  );

  real density_dependence_slope = compute_density_dependence_slope(
    birth_rate_at_capacity_to_baseline_ratio,
    density_dependence_intercept,
    carrying_capacity
  );

  // check when density dependence slope is a bad value
  if (
    is_nan(density_dependence_slope) ||
  is_inf(density_dependence_slope)
  ) {
    reject(
    "Invalid density_dependence_slope: ",
    "slope = ", density_dependence_slope,
    ", bK = ", birth_rate_at_carrying_capacity,
    ", reference birth rate = ",
      birth_rate_baseline_reference,
    ", bK / reference birth rate = ",
      birth_rate_at_carrying_capacity /
      birth_rate_baseline_reference,
    ", log ratio = ",
      log(
        birth_rate_at_carrying_capacity /
        birth_rate_baseline_reference
      ),
    ", density intercept = ",
      density_dependence_intercept,
    ", outer log argument = ",
      1.0 -
      log(
        birth_rate_at_carrying_capacity /
        birth_rate_baseline_reference
      ) / density_dependence_intercept,
    ", carrying capacity = ", carrying_capacity,
    ", max baseline birth rate = ",
      birth_rate_baseline_max,
    ", min baseline birth rate actual = ",
      birth_rate_baseline_min,
    ", density dependence scaled = ",
      density_dependence_scaled,
    ", herring intercept scaled = ",
      herring_birth_rate_midpoint,
    ", herring slope = ", herring_slope
    );
  }

  // reporting probabilities
  vector<lower=0, upper=1>[n_state_years] placental_scar_detection_probability =
  report_placental_mean * exp(-placental_scar_detection_noise * report_placental_sd);
  vector<lower=0, upper=1>[n_state_years] ca_detection_probability =
  report_ca_mean * exp(-ca_detection_noise * report_ca_sd);

  // ----------------------------
  // YEAR-TO-YEAR STATE TRANSITION
  // ----------------------------

  // precompute baseline birth rates for each year
  // baseline birth rates just depend on herring
  vector<lower=0, upper=1>[n_state_years + 1] birth_rate_baseline = compute_baseline_birth_rate(
    birth_rate_baseline_min,
    birth_rate_baseline_max,
    herring_birth_rate_midpoint,
    herring_slope,
    herring_weight,
    herring_index_baltic_proper_gulf_finland,
    herring_index_gulf_bothnia
  );

  // ----------------------------
  // STATE PROCESS OUTPUTS
  // ----------------------------

  vector<lower=0, upper=1>[n_state_years] birth_rate; // yearly birth rate
  vector<lower=0, upper=1>[n_state_years] pregnancy_rate; // yearly pregnancy rate

  vector[n_state_years] population_total; // yearly total population estimate
  vector[n_state_years] non_pup_population_total; // yearly population estimate excluding pups (ages 1+)

  matrix[n_demo_groups, n_state_years] population_comp; // yearly demographic composition

  matrix[n_demo_groups, n_state_years] survivors; // yearly survivors per demo group
  matrix[n_demo_groups, n_state_years] deaths_or_bycatch; // yearly deaths or bycatch
  matrix[n_demo_groups, n_state_years] hunted_sweden; // yearly hunted in sweden
  matrix[n_demo_groups, n_state_years] hunted_finland; // yearly hunted in finland
  matrix[n_demo_groups, n_state_years] non_hunting_deaths_expected; // expected natural deaths plus bycatch

  vector[n_state_years] hunting_bag_total_sweden; // yearly total hunting bag sweden
  vector[n_state_years] hunting_bag_total_finland; // yearly total hunting bag finland

  vector<lower=0>[n_state_years] hunted_total; // yearly total hunted

  matrix[4, n_state_years] reproductive_probs; // yearly probabilities of reproductive signs


  // starting birth rate
  real birth_rate_initial =
  update_birth_rate(
    birth_rate_baseline[1],
    density_dependence_intercept,
    density_dependence_slope,
    sum(population_init)
  );


  // ----------------------------
  // STATE PROCESS TUPLES
  // ----------------------------

  // starting demographics
  tuple(vector[n_demo_groups], real, real) init_state =
initialize_population_with_burnin(
  population_init,
  population_init_size,
  population_burn_in,
  birth_rate_initial,
  aging_matrix,
  survival_probability,
  n_age_classes
);


  (birth_rate,
  pregnancy_rate,
  population_total,
  non_pup_population_total,
  population_comp,
  survivors,
  deaths_or_bycatch,
  hunted_sweden,
  hunted_finland,
  non_hunting_deaths_expected,
  hunting_bag_total_sweden,
  hunting_bag_total_finland,
  hunted_total,
  reproductive_probs) = run_state_process_from_first_population(
    n_state_years,
    n_age_classes,
    init_state.1,
    init_state.2,
    init_state.3,
    birth_rate_baseline,
    density_dependence_intercept,
    density_dependence_slope,
    aging_matrix,
    non_hunting_mortality_rate,
    hunting_selectivity_sweden,
    hunting_selectivity_finland,
    hunting_quota_sweden,
    hunting_quota_finland,
    hunting_effort_sd_sweden,
    hunting_effort_sd_finland,
    hunting_effort_noise_sweden,
    hunting_effort_noise_finland,
    pregnancy_exposure_scaled,
    t_birth_to_start_hunt,
    hunting_duration,
    birth_count_noise,
    pup_sex_allocation_noise,
    transition_noise_raw,
    placental_scar_detection_probability,
    ca_detection_probability,
    ca_probability_without_birth,
    rel_tol,
    abs_tol,
    max_num_steps,
    fate_probability_tolerance
  );


  // ----------------------------
  // PRIORS
  // ----------------------------

  real lprior = 0;

  //Initial population size
  real lprior_population_init_size = lognormal_lpdf(population_init_size | prior_initial_population_log_mean, prior_initial_population_log_sd);
  lprior += lprior_population_init_size;

  // Natural mortality
  // female_adult_survival_probability ~ uniform(0,1); // implied prior by constraints
  real lprior_pup_to_adult_survival_ratio = beta_lpdf(pup_to_adult_survival_ratio | 1, 1);
  lprior += lprior_pup_to_adult_survival_ratio;

  real lprior_male_pup_mortality_log_offset = normal_lpdf(male_pup_mortality_log_offset | prior_male_pup_mortality_log_offset_location, prior_male_pup_mortality_log_offset_scale);
  lprior += lprior_male_pup_mortality_log_offset;

  real lprior_male_adult_mortality_log_offset = normal_lpdf(male_adult_mortality_log_offset | prior_male_adult_mortality_log_offset_location, prior_male_adult_mortality_log_offset_scale);
  lprior += lprior_male_adult_mortality_log_offset;

  // Carrying capacity
  real lprior_carrying_capacity = lognormal_lpdf(carrying_capacity | prior_carrying_capacity_log_mean, prior_carrying_capacity_log_sd);
  lprior += lprior_carrying_capacity;

  // Hunting and bycatch bias

  real lprior_hunting_selectivity_sweden_sc = normal_lpdf(hunting_selectivity_sweden_sc | 0, prior_hunting_selectivity_sd);
  lprior += lprior_hunting_selectivity_sweden_sc;

  real lprior_hunting_selectivity_finland_sc = normal_lpdf(hunting_selectivity_finland_sc | 0, prior_hunting_selectivity_sd);
  lprior += lprior_hunting_selectivity_finland_sc;

  real lprior_bycatch_selectivity_sc = normal_lpdf(bycatch_selectivity_sc | 0, prior_bycatch_selectivity_sd);
  lprior += lprior_bycatch_selectivity_sc;

  // Hunting effort sd
  real lprior_hunting_effort_sd_sweden = student_t_lpdf(hunting_effort_sd_sweden | 4, prior_hunting_effort_sd_location, prior_hunting_effort_sd_scale);
  lprior += lprior_hunting_effort_sd_sweden;

  real lprior_hunting_effort_sd_finland =  student_t_lpdf(hunting_effort_sd_finland | 4, prior_hunting_effort_sd_location, prior_hunting_effort_sd_scale);
  lprior += lprior_hunting_effort_sd_finland;

  // Birth rate
  // b0max ~ uniform (0, 1); // implied prior by bounds
  // b0min_sc ~ uniform (0, 1); // implied prior by bounds

  real lprior_birth_rate_diff = normal_lpdf(
    compute_birth_rate_replacement(non_hunting_mortality_rate, n_age_classes, female_adult_survival_probability, 1)
    - log(birth_rate_baseline_reference)
    - log(birth_rate_at_capacity_to_baseline_ratio) | 0, 0.05
  );
  lprior += lprior_birth_rate_diff;

  real lprior_herring_birth_rate_midpoint = normal_lpdf(herring_birth_rate_midpoint | 0, prior_herring_birth_rate_midpoint_sd);
  lprior += lprior_herring_birth_rate_midpoint;

  real lprior_herring_slope = normal_lpdf(herring_slope | 0, prior_herring_slope_sd);
  lprior += lprior_herring_slope;

  // herring_weight ~ uniform(0, 1); // implied prior by bounds

  // Observation of aerial survey
  real lprior_aerial_count_mu = beta_lpdf(aerial_count_mu | prior_aerial_detection_alpha, prior_aerial_detection_beta);
  lprior += lprior_aerial_count_mu;

  real lprior_aerial_count_overdispersion = lognormal_lpdf(aerial_count_overdispersion | prior_aerial_overdispersion_log_mean, prior_aerial_overdispersion_log_sd);
  lprior += lprior_aerial_count_overdispersion;

  // Observation of reproductive signs with fixed half-normal(0, 0.1) SD priors
  // ca_probability_without_birth ~ uniform(0, 1); // implied prior by constraints
  // report_placental_mean ~ uniform(0, 1); // implied prior by constraints
  // report_ca_mean ~ uniform(0, 1); // implied prior by constraints

  real lprior_report_placental_sd = normal_lpdf(report_placental_sd | 0, 0.1);
  lprior += lprior_report_placental_sd;

  real lprior_report_ca_sd = normal_lpdf(report_ca_sd | 0, 0.1);
  lprior += lprior_report_ca_sd;
}

model {

  // ----------------------------
  // PRIORS
  // ----------------------------

  target += lprior;

  // standard normals for stochasticity
  hunting_effort_noise_sweden ~ std_normal();
  hunting_effort_noise_finland ~ std_normal();

  ca_detection_noise ~ std_normal();
  placental_scar_detection_noise ~ std_normal();

  // Stochasicity for birth process
  birth_count_noise ~ std_normal();
  pup_sex_allocation_noise ~ std_normal();

  // Stochasticity for state transitions
  for (i in 1:n_state_years){
    transition_noise_raw[,i] ~ std_normal();
  }

  // ----------------------------
  // LIKELIHOODS
  // ----------------------------

  if (prior_only != 1) {

    // Aerial surveys
    target += aerial_count_lpmf(
      obs_aerial_count |
      aerial_year,
      non_pup_population_total,
      aerial_count_mu,
      aerial_count_overdispersion
    );

    // Hunting totals
    target += hunting_bags_lpdf(
      obs_hunting_bag_finland |
      hunting_bag_year_finland,
      hunting_bag_total_finland,
      hunting_bag_cv
    );
    target += hunting_bags_lpdf(
      obs_hunting_bag_sweden |
      hunting_bag_year_sweden,
      hunting_bag_total_sweden,
      hunting_bag_cv
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
      non_hunting_deaths_expected,
      bycatch_selectivity
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
}

generated quantities {

  // ----------------------------
  // IN-SAMPLE PREDICTIONS
  // ----------------------------
  // Generate one prediction for every observed aerial-survey year.
  array[n_aerial_years] int aerial_count_pred =
    aerial_count_rng(
      aerial_year,
      non_pup_population_total,
      aerial_count_mu,
      aerial_count_overdispersion
    );

  // Generate one prediction for every observed national hunting bag.
  array[n_hunting_bag_years_sweden] real hunting_bags_sweden_pred =
    hunting_bags_rng(
      hunting_bag_year_sweden,
      hunting_bag_total_sweden,
      hunting_bag_cv
    );

  array[n_hunting_bag_years_finland] real hunting_bags_finland_pred =
    hunting_bags_rng(
      hunting_bag_year_finland,
      hunting_bag_total_finland,
      hunting_bag_cv
    );

  array[n_hunting_comp_years_sweden, n_demo_groups] int
    hunting_comp_sweden_pred =
      hunting_comp_rng(
        hunting_comp_year_sweden,
        hunted_sweden,
        hunting_bag_total_sweden,
        hunting_sample_size_sweden_pred
      );

  array[n_hunting_comp_years_finland, n_demo_groups] int
    hunting_comp_finland_pred =
      hunting_comp_rng(
        hunting_comp_year_finland,
        hunted_finland,
        hunting_bag_total_finland,
        hunting_sample_size_finland_pred
      );

  array[n_bycatch_years, n_demo_groups] int bycatch_comp_pred =
    bycatch_comp_rng(
      bycatch_comp_year,
      non_hunting_deaths_expected,
      bycatch_selectivity,
      bycatch_sample_size_pred
    );

  array[n_pregnancy_years] int pregnancy_count_pred =
    pregnancy_rng(
      pregnancy_count_year,
      pregnancy_sample_size,
      pregnancy_rate
    );

  array[n_reproductive_years, 4] int
    reproductive_signs_finland_pred =
      reproductive_signs_rng(
        reproductive_signs_year,
        reproductive_probs,
        reproductive_signs_sample_size_pred
      );


  // ----------------------------
  // POINTWISE LOG LIKELIHOODS
  // ----------------------------

  // Aerial surveys
  vector[n_aerial_years] log_lik_aerial_count = 
    aerial_count_pointwise_log_lik(
      obs_aerial_count,
      aerial_year,
      non_pup_population_total,
      aerial_count_mu,
      aerial_count_overdispersion
    );

  // Hunting bag totals
  vector[n_hunting_bag_years_finland] log_lik_hunting_bags_finland =
    hunting_bags_pointwise_log_lik(
      obs_hunting_bag_finland,
      hunting_bag_year_finland,
      hunting_bag_total_finland,
      hunting_bag_cv
    );

  vector[n_hunting_bag_years_sweden] log_lik_hunting_bags_sweden =
    hunting_bags_pointwise_log_lik(
      obs_hunting_bag_sweden,
      hunting_bag_year_sweden,
      hunting_bag_total_sweden,
      hunting_bag_cv
    );

  // Hunting composition
  vector[n_hunting_comp_years_finland] log_lik_hunting_comp_finland = 
    hunting_comp_pointwise_log_lik(
      obs_hunting_comp_finland,
      hunting_comp_year_finland,
      hunted_finland,
      hunting_bag_total_finland
    );

  vector[n_hunting_comp_years_sweden] log_lik_hunting_comp_sweden = 
    hunting_comp_pointwise_log_lik(
      obs_hunting_comp_sweden,
      hunting_comp_year_sweden,
      hunted_sweden,
      hunting_bag_total_sweden
    );

  // Bycatch
  vector[n_bycatch_years] log_lik_bycatch_comp = 
    bycatch_comp_pointwise_log_lik(
      obs_bycatch_comp,
      bycatch_comp_year,
      non_hunting_deaths_expected,
      bycatch_selectivity
    );

  // Pregnancy
  vector[n_pregnancy_years] log_lik_pregnancy_count = 
    pregnancy_pointwise_log_lik(
      obs_pregnancy_count,
      pregnancy_count_year,
      pregnancy_sample_size,
      pregnancy_rate
    );

  // Reproductive signs
  vector[n_reproductive_years] log_lik_reproductive_signs = 
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
  vector[n_demo_groups] survivors_final = survivors[, n_state_years];

}
