library(dplyr)
library(forcats)
library(tidyr)
options(readr.show_col_types = FALSE)

# Aerial counts
aerial_counts_df <- read_csv("Aerial_count_total.csv") |>
  pivot_longer(-Area, names_to = "Year", values_to = "Count") |>
  filter(Area == "Total") |>
  filter(!is.na(Count)) |>
  select(Year, Count) |>
  mutate(Year = as.numeric(Year))

# Hunting bags
hunting_bag_fi <- read_delim("Hunting_bag_FI.csv") |>
    filter(!is.na(FINLAND)) |>
    pivot_longer(-FINLAND, names_to = "Year") |>
    filter(!str_detect(FINLAND, "Data sources:")) |>
  separate(FINLAND, into = c("Type", "Region"), sep = " -") |>
    mutate(Year = replace_when(
    Year,
    Year == "2016/17" ~ "2017",
    Year == "2017/18" ~ "2018",
    Year == "2018/19" ~ "2019",
    Year == "2019/20" ~ "2020",
    Year == "2020/21" ~ "2021",
    Year == "2021/22" ~ "2022",
    Year == "2022/23" ~ "2023",
    Year == "2023/24" ~ "2024",
    Year == "2024/25" ~ "2025")
    )

hunting_bag_total_fi <- hunting_bag_fi |>
  filter(Type == "Total bag") |>
  group_by(Year) |>
  summarise(Count = sum(value)) |>
  mutate(Country = "Finland") |>
  arrange(Year)

hunting_quota_finland <-  hunting_bag_fi |>
  filter(Type == "Quota") |>
  group_by(Year) |>
  summarise(Quota = sum(value)) |>
  mutate(Country = "Finland") |>
  arrange(Year)


hunting_bag_sw <- read_delim("Hunting_bag_SW.csv", skip = 1) |>
    rename(Region = `...1`) |>
    pivot_longer(-Region, names_to = "Year") |>
    mutate(type = if_else(Region == "Quota", "Quota", "Count"))


hunting_quota_sweden <- hunting_bag_sw |>
  filter(type == "Quota") |>
  arrange(Year) |>
  rename(Quota = value) |>
  select(Year, Quota) |>
  mutate(Country = "Sweden")

hunting_bag_total_sw <- hunting_bag_sw |>
  filter(Region == "Sum", type == "Count") |>
  rename(Count = value) |>
  select(Year, Count) |>
  mutate(Country = "Sweden")


hunting_bags_df <- bind_rows(hunting_bag_total_fi, hunting_bag_total_sw)
hunting_quotas_df <- bind_rows(hunting_quota_finland, hunting_quota_sweden)


# hunting and bycatch samples (finland)
samples_fi <- read_csv("Samples_FI.csv",
                       guess_max = 4000) |>
    filter(!is.na(Year)) |>
    mutate(
        Source = case_when(
            Source %in% c("by-catch", "By-catch", "by-caught", "by-Caught", "By-caught", "Stranded") ~ "Bycatch",
            Source %in% c("hunted", "Hunted") ~ "Hunt",
            TRUE ~ NA_character_
            ),
        Sex = case_when(
            Sex %in% c(0, "2?") ~ NA_character_,
            Sex == "1" ~ "Male",
            Sex == "2" ~ "Female",
            TRUE ~ NA_character_
        ),
        `Age class` = case_when(
            `Age class` == "adult" ~ "Adult",
            `Age class` == "juvenile" ~ "Juvenile",
            TRUE ~ as.character(`Age class`)
        ),
        Pregnant = replace_when(
            Pregnent,
            !(Pregnent %in% c(0, 1)) ~ NA
        ),
        Month = replace_when(
            Month,
            Month == "4.-5." ~ "5",
            Month == "kevät" ~ "5",
            Month == "6-7" ~ "7",
            Month == "16" ~ NA
        ),
        Age = replace_when(
            `Age (years)`,
            `Age (years)` %in% seq(0, 4) ~ `Age (years)`,
            as.numeric(`Age (years)`) >= 5 ~ "5+"
        ),
        Age.Sex = if_else(!is.na(Age) & !is.na(Sex), paste(Sex, Age), NA),
        CA = if_else(Age.Sex != "Female 5+", NA, CA),
        Scar = if_else(Age.Sex != "Female 5+", NA, Scar)
    ) |>
  mutate(
    Source = case_when(
      `Mortality reason` %in% c(4, 6) ~ "Bycatch",
      `Mortality reason` %in% c(1,2,3,5) ~ "Hunt",
      TRUE ~ Source
    )
  )


samples_comp_fi <- samples_fi |>
  filter(!is.na(Sex),
         !is.na(Age),
         !is.na(Source)
         ) |>
  select(Year, Age, Sex, Source) |>
  mutate(Country = "Finland")



# reproductive signs
reproductive_signs_df <- samples_fi |>
    filter(
        Month %in% c(4, 5, 6),
        !is.na(Scar),
        !is.na(CA)
    ) |>
  mutate(Sign = case_when(
    Scar == 1 & CA == 1 ~ "both",
    Scar == 1 ~ "scar",
    CA == 1 ~ "CA",
    .default = "none"
  )
  ) |>
  select(Year, Sign)
    

# hunting and bycatch samples (sweden)
col_types <- readr::cols(.default = readr::col_character())

samples_sw <- read_csv("Samples_SW.csv", skip = 1, col_types = col_types)

samples_sw_24 <- read_csv("Samples_SW_2024.csv", skip = 1, col_types = col_types)

samples_sw <- samples_sw |>
  bind_rows(samples_sw_24 |> rename("Ålder (formel)" = "Ålder")) |>
    rename(
        Year = "Fyndår",
        Sex = "Kön",
        Found = "Fyndomständighet",
        Age = "Ålder (formel)",
        Month = "Fyndmånad",
        Weight = "Kroppsvikt kg",
        "Uterus status" = "Uterusstatus",
        "Pregnancy status" = "Dräktighetsstatus HELCOM",
        Length = "Kroppslängd cm",
        "Blubber thickness" = "Späck (mm)"
    ) |>
  select(
    Year,
    Sex,
    Found,
    Age,
    Month,
    Weight,
    `Uterus status`,
    `Pregnancy status`,
    Length,
    `Blubber thickness`
  ) |>
    mutate(
        Sex = replace_when(
            Sex,
            Sex == "hane" ~ "Male",
            Sex == "Hane" ~ "Male",
            Sex == "Hona" ~ "Female",
            Sex == "hona" ~ "Female",
            Sex == "Okänt" ~ NA,
            Sex == "okänt" ~ NA,
            Sex == "?" ~ NA
        ),
        Source = case_when(
            Found %in% c("Jakt", "jakt", "Jakt?", "Hunt") ~ "Hunt",
            Found %in% c("", NA_character_) ~ NA_character_,
            TRUE ~ "Bycatch"
        ),
        `Age class` = case_when(
          Age %in% c("0", "1", "2", "3", "4") ~ Age,
          as.numeric(Age) >= 5 ~ "5+",
          Age %in% c(">25", "16+") ~ "5+"
        ),
        Weight = replace_when(
            Weight,
            Weight == "150200" ~ "175",
            Weight == "200225" ~ "212",
            Weight == "200250" ~ "225",
            as.numeric(Weight) > 500 ~ NA
        ),
        `Pregnancy status` = replace_when(
            `Pregnancy status`,
            Sex == "Male" ~ "Male"
        ),
        `Uterus status` = replace_when(
            `Uterus status`,
            `Uterus status` %in% c("Juvenil", "Juvenile", "juvenil ", "juvenil","juvenil, hormonpåverkan") ~ "Juvenile"),
        Age.Sex = paste(Sex, `Age class`)
    )


samples_comp_sw <- samples_sw |>
  filter(!is.na(Sex),
         !is.na(`Age class`),
         !is.na(Source)
         ) |>
  select(Year, Sex, `Age class`, Source) |>
  rename(Age = "Age class") |>
  select(Year, Age, Sex, Source) |>
  mutate(Country = "Sweden",
         Year = as.numeric(Year))


# composition
samples_df <- bind_rows(samples_comp_fi, samples_comp_sw)

# pregnancy signs

pregnancy_status_df <- samples_sw |>
  filter(
    `Pregnancy status` %in% c("Dräktig", "Ej dräktig"),
    `Uterus status` != "Juvenile",
    (`Age class` == "5+" | `Age class` == "Unknown" |  is.na(`Age class`)),
    as.numeric(Month) >= 8
  ) |>
  mutate(Status = case_when(
    `Pregnancy status` == "Dräktig" ~ "pregnant",
    `Pregnancy status` == "Ej dräktig" ~ "not pregnant"
  )
  ) |>
  select(Year, Status)

pregnancy_status_df <- samples_sw |>
  dplyr::mutate(
    Year = as.integer(Year),
    Month = as.numeric(Month),
    Source = dplyr::coalesce(Source, "Unknown"),
    `Uterus status` = dplyr::coalesce(`Uterus status`, ""),
    `Age class` = dplyr::coalesce(`Age class`, "Unknown"),
    `Pregnancy status` = dplyr::coalesce(`Pregnancy status`, "Unknown")
  ) |>
  dplyr::filter(
    `Pregnancy status` %in% c("Dräktig", "Ej dräktig"),
    Source == "Hunt",
    Month >= 8,
    `Uterus status` != "Juvenile",
    `Age class` %in% c("5+", "Unknown")
  ) |>
  dplyr::mutate(
    Status = dplyr::case_when(
      `Pregnancy status` == "Dräktig" ~ "pregnant",
      `Pregnancy status` == "Ej dräktig" ~ "not pregnant"
    )
  ) |>
  dplyr::select(Year, Status)

# Herring

herring_weight_BP_GoF <- read_delim(
  "Herring_mean_weight_25-29_32.csv",
  delim = " "
) |>
  pivot_longer(
    -Year,
    names_to = "Age",
    values_to = "Weight"
  ) |>
  mutate(
    Region = "Gulf of Finland",
    Weight = Weight * 1000
  )

herring_weight_BP_BB <- read_csv(
  "Herring_mean_weight_GoB.csv"
) |>
  pivot_longer(
    -Year,
    names_to = "Age",
    values_to = "Weight"
  ) |>
  mutate(
    Region = "Bothnian Bay"
  )

herring_catch_BP_GoF <- read_delim(
  "Herring_catch_in_numbers_25-29_32.csv"
) |>
  select(-`SOPCOF %`) |>
  pivot_longer(
    -Year,
    names_to = "Age",
    values_to = "Catch"
  ) |>
  mutate(
    Region = "Gulf of Finland",
    Age = str_remove(Age, "Age ")
  )

herring_catch_BP_GoB <- read_csv(
  "Herring_catch_in_numbers_GoB.csv"
) |>
  select(-`...1`) |>
  pivot_longer(
    -Year,
    names_to = "Age",
    values_to = "Catch"
  ) |>
  mutate(
    Region = "Bothnian Bay"
  )

herring_df <- inner_join(
  herring_weight_BP_GoF,
  herring_catch_BP_GoF,
  by = c("Year", "Age", "Region")
) |>
  bind_rows(
    inner_join(
      herring_weight_BP_BB,
      herring_catch_BP_GoB,
      by = c("Year", "Age", "Region")
    )
  )

# Reproduce original Age5plus calculation exactly

herring_index_baltic_proper_gulf_finland <- {

  catch_weights <- herring_df |>
    filter(
      Region == "Gulf of Finland",
      Age %in% c("5", "6", "7", "8+")
    ) |>
    group_by(Age) |>
    summarise(
      mean_catch = mean(Catch, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      w = mean_catch / sum(mean_catch)
    )

  herring_df |>
    filter(
      Region == "Gulf of Finland",
      Age %in% c("5", "6", "7", "8+")
    ) |>
    left_join(
      catch_weights |> select(Age, w),
      by = "Age"
    ) |>
    group_by(Year) |>
    summarise(
      age5plus = sum(Weight * w),
      .groups = "drop"
    ) |>
    arrange(Year) |>
    mutate(
      mean_weight_scaled = as.numeric(scale(age5plus))
    ) |>
    pull(mean_weight_scaled)
}

herring_index_gulf_bothnia <- {

  catch_weights <- herring_df |>
    filter(
      Region == "Bothnian Bay",
      Age %in% c(
        as.character(5:14),
        "15+"
      )
    ) |>
    group_by(Age) |>
    summarise(
      mean_catch = mean(Catch, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      w = mean_catch / sum(mean_catch)
    )

  herring_df |>
    filter(
      Region == "Bothnian Bay",
      Age %in% c(
        as.character(5:14),
        "15+"
      )
    ) |>
    left_join(
      catch_weights |> select(Age, w),
      by = "Age"
    ) |>
    group_by(Year) |>
    summarise(
      age5plus = sum(Weight * w),
      .groups = "drop"
    ) |>
    arrange(Year) |>
    mutate(
      mean_weight_scaled = as.numeric(scale(age5plus))
    ) |>
    pull(mean_weight_scaled)
}


grey_seal_data <- list(
  hunting_quotas = hunting_quotas_df,
  herring = herring_df,
  aerial_counts = aerial_counts_df,
  hunting_bags = hunting_bags_df,
  samples = samples_df,
  reproductive_signs = reproductive_signs_df,
  pregnancy_status = pregnancy_status_df
)

usethis::use_data(grey_seal_data, overwrite = TRUE)
