#' @importFrom rlang .data
NULL

build_stan_data <- function(species, data, years, prior_spec = NULL) {
  switch(
    species,
    grey = build_grey_stan_data(data, years = years, prior_spec)
  )
}

DEMOGRAPHIC_GROUPS <- c(
  paste("Female", 0:4),
  "Female 5+",
  paste("Male", 0:4),
  "Male 5+"
)

reindex_years <- function(df, year_col, years_keep, offset = 0) {

  year_map <- stats::setNames(
    seq_along(years_keep) + offset,
    years_keep
  )

  df |>
    dplyr::mutate(
      Year = as.integer(.data[[year_col]])
    ) |>
    dplyr::filter(
      .data$Year %in% years_keep
    ) |>
    dplyr::arrange(.data$Year) |>
    dplyr::mutate(
      Year = unname(
        year_map[as.character(.data$Year)]
      )
    )
}

build_grey_herring <- function(herring, years) {

  out <- list()

  compute_index <- function(dat, region, ages_5plus) {

    dat_region <- dat |>
      dplyr::filter(
        .data$Region == region,
        .data$Year >= 2001,
        .data$Age %in% ages_5plus
      )

    age_weights <- dat_region |>
      dplyr::group_by(.data$Age) |>
      dplyr::summarise(
        mean_catch = mean(.data$Catch, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        weight = .data$mean_catch / sum(.data$mean_catch)
      )

    full_index <- dat_region |>
      dplyr::left_join(
        age_weights |>
          dplyr::select(.data$Age, .data$weight),
        by = "Age"
      ) |>
      dplyr::group_by(.data$Year) |>
      dplyr::summarise(
        age5plus = sum(.data$Weight * .data$weight),
        .groups = "drop"
      ) |>
      dplyr::arrange(.data$Year)

    scaling_mean <- mean(full_index$age5plus)
    scaling_sd <- stats::sd(full_index$age5plus)

    full_index <- full_index |>
      dplyr::mutate(
        age5plus_scaled =
          (.data$age5plus - scaling_mean) /
          scaling_sd
      )

    full_index |>
      dplyr::filter(
        .data$Year %in% c(
          seq(
            min(years) - 2,
            max(years)
          )
        )
      ) |>
      dplyr::pull(.data$age5plus_scaled)
  }

  out$herring_index_baltic_proper_gulf_finland <-
    compute_index(
      herring,
      "Gulf of Finland",
      c("5", "6", "7", "8+")
    )

  out$herring_index_gulf_bothnia <-
    compute_index(
      herring,
      "Bothnian Bay",
      c(
        as.character(5:14),
        "15+"
      )
    )

  out
}

build_grey_hunting_quotas <- function(hunting_quotas, years) {

  out <- list()

  out$hunting_quota_finland <-
    hunting_quotas |>
    dplyr::filter(
      .data$Country == "Finland"
    ) |>
    reindex_years(
      year_col = "Year",
      years_keep = years,
      offset = 0
    )  |>
    dplyr::pull(.data$Quota)

  out$hunting_quota_sweden <-
    hunting_quotas |>
    dplyr::filter(
      .data$Country == "Sweden"
    ) |>
    reindex_years(
      year_col = "Year",
      years_keep = years,
      offset = 0
    ) |>
    dplyr::pull(.data$Quota)

  out
}

build_grey_aerial_counts <- function(aerial_counts, years) {

  out <- list()

  aerial_counts <-
    aerial_counts |>
    reindex_years(
      year_col = "Year",
      years_keep = years,
      offset = 1
    )

  out$n_aerial_years <- nrow(aerial_counts)
  out$aerial_year <- as.integer(aerial_counts$Year)
  out$obs_aerial_count <- aerial_counts$Count

  out
}

build_grey_hunting_bags <- function(hunting_bags, years) {

  out <- list()

  hunting_bag_finland <-
    hunting_bags |>
    dplyr::filter(
      .data$Country == "Finland"
    ) |>
    reindex_years(
      year_col = "Year",
      years_keep = years,
      offset = 0
    )

  out$n_hunting_bag_years_finland <-
    nrow(hunting_bag_finland)

  out$hunting_bag_year_finland <-
    as.integer(hunting_bag_finland$Year)

  out$obs_hunting_bag_finland <-
    hunting_bag_finland$Count

  hunting_bag_sweden <-
    hunting_bags |>
    dplyr::filter(
      .data$Country == "Sweden"
    ) |>
    reindex_years(
      year_col = "Year",
      years_keep = years,
      offset = 0
    )

  out$n_hunting_bag_years_sweden <-
    nrow(hunting_bag_sweden)

  out$hunting_bag_year_sweden <-
    as.integer(hunting_bag_sweden$Year)

  out$obs_hunting_bag_sweden <-
    hunting_bag_sweden$Count

  out
}

build_grey_hunting_samples <- function(samples, years) {

  out <- list()

  hunting_comp_finland <-
    samples |>
    dplyr::filter(
      .data$Source == "Hunt",
      .data$Country == "Finland",
      !is.na(.data$Sex),
      !is.na(.data$Age)
    ) |>
    dplyr::mutate(
      Group = factor(
        paste(.data$Sex, .data$Age),
        levels = DEMOGRAPHIC_GROUPS
      )
    ) |>
    reindex_years(
      year_col = "Year",
      years_keep = years,
      offset = 0
    ) |>
    dplyr::select(
      .data$Year,
      .data$Group
    ) |>
    table()

  out$n_hunting_comp_years_finland <-
    nrow(hunting_comp_finland)

  out$hunting_comp_year_finland <-
    as.integer(
      rownames(hunting_comp_finland)
    )

  out$obs_hunting_comp_finland <-
    unclass(hunting_comp_finland)

  hunting_comp_sweden <-
    samples |>
    dplyr::filter(
      .data$Source == "Hunt",
      .data$Country == "Sweden",
      !is.na(.data$Sex),
      !is.na(.data$Age)
    ) |>
    dplyr::mutate(
      Group = factor(
        paste(.data$Sex, .data$Age),
        levels = DEMOGRAPHIC_GROUPS
      )
    ) |>
    reindex_years(
      year_col = "Year",
      years_keep = years,
      offset = 0
    ) |>
    dplyr::select(
      .data$Year,
      .data$Group
    ) |>
    table()

  out$n_hunting_comp_years_sweden <-
    nrow(hunting_comp_sweden)

  out$hunting_comp_year_sweden <-
    as.integer(
      rownames(hunting_comp_sweden)
    )

  out$obs_hunting_comp_sweden <-
    unclass(hunting_comp_sweden)

  out
}

build_grey_bycatch <- function(samples, years) {

  out <- list()

  bycatch <-
    samples |>
    dplyr::filter(
      .data$Source == "Bycatch",
      !is.na(.data$Sex),
      !is.na(.data$Age)
    ) |>
    dplyr::mutate(
      Group = factor(
        paste(.data$Sex, .data$Age),
        levels = DEMOGRAPHIC_GROUPS
      )
    ) |>
    reindex_years(
      year_col = "Year",
      years_keep = years,
      offset = 0
    ) |>
    dplyr::select(
      .data$Year,
      .data$Group
    ) |>
    table()

  out$n_bycatch_years <-
    nrow(bycatch)

  out$bycatch_comp_year <-
    as.integer(
      rownames(bycatch)
    )

  out$obs_bycatch_comp <-
    unclass(bycatch)

  out
}

build_grey_reproductive_signs <- function(reproductive_signs, years) {

  out <- list()

  reproductive_signs <-
    reproductive_signs |>
    reindex_years(
      year_col = "Year",
      years_keep = years,
      offset = 0
    ) |>
    dplyr::mutate(
      Sign = factor(
        .data$Sign,
        levels = c(
          "none",
          "scar",
          "CA",
          "both"
        )
      )
    ) |>
    dplyr::select(
      .data$Year,
      .data$Sign
    ) |>
    table()

  out$n_reproductive_years <-
    nrow(reproductive_signs)

  out$reproductive_signs_year <-
    as.integer(
      rownames(reproductive_signs)
    )

  out$obs_reproductive_signs_finland <-
    unclass(reproductive_signs)

  out
}

build_grey_pregnancy_status <- function(pregnancy, years) {

  out <- list()

  pregnancy_status <-
    pregnancy |>
    reindex_years(
      year_col = "Year",
      years_keep = years,
      offset = 0
    ) |>
    dplyr::mutate(
      Status = factor(
        .data$Status,
        levels = c("not pregnant", "pregnant")
      )
    ) |>
    dplyr::select(
      .data$Year,
      .data$Status
    ) |>
    table()

  out$n_pregnancy_years <-
    nrow(pregnancy_status)

  out$pregnancy_count_year <-
    as.integer(
      rownames(pregnancy_status)
    )

  out$obs_pregnancy_count <-
    pregnancy_status[, "pregnant"]

  out$pregnancy_sample_size <-
    rowSums(pregnancy_status)

  out
}

build_grey_stan_data <- function(data,
                                 years,
                                 prior_spec) {

  process_years <- c(min(years) - 1, years)
  
  c(
    prior_spec,
    list(
      n_state_years = length(process_years),
      n_age_classes = 6
    ),

    build_grey_herring(
      data$herring,
      years
    ),

    build_grey_hunting_quotas(
      data$hunting_quotas,
      process_years
    ),

    build_grey_aerial_counts(
      data$aerial_counts,
      years
    ),

    build_grey_hunting_bags(
      data$hunting_bags,
      process_years
    ),

    build_grey_hunting_samples(
      data$samples,
      process_years
    ),

    build_grey_bycatch(
      data$samples,
      process_years
    ),

    build_grey_reproductive_signs(
      data$reproductive_signs,
      process_years
    ),

    build_grey_pregnancy_status(
      data$pregnancy_status,
      process_years
    )
  )
}
