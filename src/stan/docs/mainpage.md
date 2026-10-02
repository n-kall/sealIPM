# sealIPM Stan reference

This reference describes the Stan functions used to fit and forecast the
Baltic grey seal integrated population model. Use the **Function index** to
find a function by name, or the search box to find a parameter or topic.

## Start here

| Component | Reference | Purpose |
| :--- | :--- | :--- |
| Process runners and initialization | @ref run_process.stanfunctions | Continuous fitting, discrete forecasting, and initial population burn-in |
| Pregnancy and births | @ref pregnancy_births.stanfunctions | Baseline and density-dependent reproductive rates |
| Deaths | @ref deaths.stanfunctions | Demographic non-hunting mortality rates |
| Hunting | @ref hunting.stanfunctions | Hunting pressure, ODE derivative, and annual fate probabilities |
| Aging | @ref aging.stanfunctions | Age-class transitions and newborn population updates |
| Fates | @ref fates.stanfunctions | Logistic-normal allocation of expected fate counts |

## Observation models

| Observations | Reference |
| :--- | :--- |
| Aerial counts (ages 1+, excluding pups) | @ref aerial_counts.stanfunctions |
| Hunting bags | @ref hunting_bags.stanfunctions |
| Hunting age and sex composition | @ref hunting_composition.stanfunctions |
| Bycatch age and sex composition | @ref bycatch_composition.stanfunctions |
| Pregnancy | @ref pregnancy.stanfunctions |
| Reproductive signs | @ref reproductive_signs.stanfunctions |

## Variable meanings

| Variable | Interpretation |
| :--- | :--- |
| `t_birth_to_start_hunt` | Time from birth to hunting onset in years |
| `hunting_duration` | Shared hunting season length in years, not time from birth to the end of hunting |
| `male_pup_mortality_log_offset`, `male_adult_mortality_log_offset` | Male minus female log mortality, so a positive offset increases mortality and reduces survival |
| `herring_birth_rate_midpoint` | Weighted herring index where baseline birth rate is halfway between its minimum and maximum, using `herring_slope * (weighted_herring_index - herring_birth_rate_midpoint)` |
| `non_hunting_deaths_expected` | Expected natural deaths plus bycatch, used as the demographic basis of the bycatch composition model rather than as bycatch death counts |
| `carrying_capacity` | Carrying capacity at reference herring conditions |

The hunting endpoint relative to birth is `t_birth_to_start_hunt + hunting_duration`
The two countries retain the paper's shared hunting-season approximation

The herring midpoint is -alpha / beta in the paper's logistic predictor
It replaces the previous additive offset with its negative
The zero-centered normal prior is unchanged by this sign reversal

## Annual process timeline

For years after the initial state, update the birth rate using the previous
population total, age the previous survivors, and add newborn pups
Calculate hunting pressures and annual fate probabilities from this population
The hunting probabilities include survival between birth and hunting onset
Annual survival already includes a full year of non-hunting mortality
The survivors form the boundary for the next year's ageing and births

Pregnancy in the current year uses the next year's baseline birth rate and
the current population density
Aerial counts use both sexes aged one and older and exclude pups of the year
Fitting uses continuous demographic approximations while forecasting uses
integer binomial births and multinomial fate draws

## Prior configuration

The default configuration also contains fixed process inputs and solver controls
Only one hunting duration is supplied
The reporting SD priors remain fixed half-normal distributions with location zero
and scale 0.1, so there are no unused configurable reporting SD hyperparameters
All saved `lprior` components remain available for priorsense selection
