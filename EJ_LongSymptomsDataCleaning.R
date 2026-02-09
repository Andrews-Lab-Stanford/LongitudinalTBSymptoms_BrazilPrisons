################################################################################
#              DATA CLEANING FOR LONGITUDINAL SYMPTOMS ANALYSES
# Last Updated: August 12
# By: Esther Jung
################################################################################

#### Load Libraries =====
library(dplyr)
library(stringr)
library(tidyr)
library(networkD3)
library(htmlwidgets)
library(ggplot2)
library(haven)

#### Load Data from Most Current REDCap Pull =====
widesymptoms <- definition1cohort %>%
  mutate(
    # Fix specific death dates
    death_date = case_when(
      record_id == 224 ~ as.Date("2023-06-16"),
      record_id == 2038 ~ as.Date("2024-06-02"),
      record_id == 3627 ~ as.Date("2024-03-27"),
      TRUE ~ death_date
    ),
    # Fix specific follow-up visit dates
    fu_date1 = case_when(
      record_id == 1460 ~ as.Date("2024-03-19"),
      record_id == 2637 ~ as.Date("2023-10-10"),
      record_id == 3607 ~ as.Date("2024-03-19"),
      TRUE ~ fu_date1
    )
  ) %>%
  select(c(
    record_id,
    round,
    enroll_date,
    visit_date,
    dob,
    prison,
    age,
    prison_time,
    prev_arrest,
    timeleft,
    preexisting_tb,
    race,
    education,
    prev_tb,
    smoking,
    csmoke,
    any_drugs,
    hiv1,
    symptomatic_cellmates,
    tb_contact,
    LunitTB,
    LunitTB_1,
    LunitTB_2,
    LunitTB_3,
    LunitTB_4,
    LunitTB_5,
    LunitTB_6,
    visit_date,
    symptom_date,
    any_symptoms, # symptoms at baseline
    sum_symptoms,
    cough,
    sputum,
    blood,
    fever,
    appetite,
    weightloss,
    sweats,
    chestpain,
    breathing,
    expectoracao_sangue1,
    visitdate_1,
    symptom_date_1,
    symptoms_1,
    sum_symptoms_1,
    cough_1,
    sputum_1,
    blood_1,
    fever_1,
    appetite_1,
    weightloss_1,
    sweats_1,
    chestpain_1,
    breathing_1,
    bloodysputum_1,
    tbfu1_1,
    tbfu1_2,
    visitdate_2,
    symptom_date_2,
    symptoms_2,
    sum_symptoms_2,
    cough_2,
    sputum_2,
    blood_2,
    fever_2,
    appetite_2,
    weightloss_2,
    sweats_2,
    chestpain_2,
    breathing_2,
    bloodysputum_2,
    tbfu2_1,
    tbfu2_2,
    visitdate_3,
    symptom_date_3,
    symptoms_3,
    sum_symptoms_3,
    cough_3,
    sputum_3,
    blood_3,
    fever_3,
    appetite_3,
    weightloss_3,
    sweats_3,
    chestpain_3,
    breathing_3,
    bloodysputum_3,
    tbfu3_1,
    tbfu3_2,
    visitdate_4,
    symptom_date_4,
    symptoms_4,
    sum_symptoms_4,
    cough_4,
    sputum_4,
    blood_4,
    fever_4,
    appetite_4,
    weightloss_4,
    sweats_4,
    chestpain_4,
    breathing_4,
    bloodysputum_4,
    tbfu4_1,
    tbfu4_2,
    visitdate_5,
    symptom_date_5,
    symptoms_5,
    sum_symptoms_5,
    cough_5,
    sputum_5,
    blood_5,
    fever_5,
    appetite_5,
    weightloss_5,
    sweats_5,
    chestpain_5,
    breathing_5,
    bloodysputum_5,
    tbfu5_1,
    tbfu5_2,
    visitdate_6,
    symptom_date_6,
    symptoms_6,
    sum_symptoms_6,
    cough_6,
    sputum_6,
    blood_6,
    fever_6,
    appetite_6,
    weightloss_6,
    sweats_6,
    chestpain_6,
    breathing_6,
    bloodysputum_6,
    tbfu6_1,
    tbfu6_2,
    TBincidence1,
    TBincidence2,
    fu_date1,
    fu_date2,
    death_date,
    fu_days,
    fu_visit
  )) %>%
  mutate(symptoms_0 = any_symptoms) %>%
  mutate(bloodysputum = expectoracao_sangue1) # Add FU 7 when complete


# Clean and classify baseline status
widesymptoms_clean <- widesymptoms %>%
  filter(symptoms_0 %in% c("Any", "None")) %>%
  select(record_id, fu_date2, symptoms_0, everything()) %>%
  mutate(
    id = row_number(),
    baseline_status = case_when(
      symptoms_0 == "Any" ~ "Symptomatic",
      symptoms_0 == "None" ~ "Asymptomatic"
    )
  )

# Pivot to long and classify visit symptom status
symptom_long <- widesymptoms_clean %>%
  select(id, TBincidence1, TBincidence2, baseline_status, starts_with("symptoms_")) %>%
  mutate(across(starts_with("symptoms_"), ~ as.character(.))) %>%
  pivot_longer(
    cols = starts_with("symptoms_"),
    names_to = "visit",
    values_to = "symptom_status"
  ) %>%
  mutate(
    visit_num = str_extract(visit, "\\d+") %>% as.integer(),
    visit_status = case_when(
      symptom_status == "Any" ~ "+",
      symptom_status == "None" ~ "-",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(visit_status))

symptom_binary <- widesymptoms_clean %>%
  mutate(across(starts_with("symptoms_"), ~ case_when(
    . == "Any" ~ 1,
    . == "None" ~ 0,
    TRUE ~ NA_real_
  )))

count_symptoms <- symptom_binary %>%
  rowwise() %>%
  mutate(
    symptom_vec = list(c_across(starts_with("symptoms_"))),

    # Remove trailing NAs
    last_non_na = max(which(!is.na(symptom_vec)), na.rm = TRUE),
    trimmed_vec = list(symptom_vec[1:last_non_na]),
    obs_vec = list(trimmed_vec[!is.na(trimmed_vec)]),
    switch_count = if (length(obs_vec) > 1) sum(diff(obs_vec) != 0) else 0,
    n_obs = length(obs_vec[[1]]),
    n_symp = sum(obs_vec[[1]] == 1),

    # Intermittent NA: NA inside the trimmed (non-trailing) portion
    is_intermittent_na = any(is.na(trimmed_vec[2:(length(trimmed_vec) - 1)]))
  ) %>%
  ungroup()

sum(count_symptoms$is_intermittent_na)
intermittents <- count_symptoms %>% filter(is_intermittent_na == TRUE)

count_clean <- count_symptoms %>%
  filter(is_intermittent_na == FALSE)

hist(count_clean$switch_count)
table(count_clean$switch_count)

#### Checking ==================================================================

# STEP 1: Pivot long and label symptom status
symptom_long <- widesymptoms_clean %>%
  select(id, fu_date2, starts_with("symptom_date_"), starts_with("symptoms_")) %>%
  pivot_longer(
    cols = matches("^symptom_date_\\d+|^symptoms_\\d+"),
    names_to = c(".value", "visit_num"),
    names_pattern = "(symptom_date|symptoms)_(\\d+)"
  ) %>%
  mutate(
    visit_num = as.integer(visit_num),
    visit_status = case_when(
      symptoms == "Any" ~ "Symptomatic",
      symptoms == "None" ~ "Asymptomatic",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(visit_status))

# STEP 2: Identify switch per individual
switch_info <- symptom_long %>%
  arrange(id, visit_num) %>%
  group_by(id) %>%
  mutate(baseline_status = first(visit_status)) %>%
  filter(visit_status != baseline_status, !is.na(symptom_date)) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    id,
    switch_direction = paste(baseline_status, "→", visit_status),
    switch_date = symptom_date
  )

# STEP 3: Merge with full dataset, classify switch timing
switch_summary <- widesymptoms_clean %>%
  select(id, fu_date2) %>%
  left_join(switch_info, by = "id") %>%
  mutate(
    switch_direction = replace_na(switch_direction, "No Switch"),
    rel_days_to_tb = as.numeric(switch_date - fu_date2),
    switch_timing_class = case_when(
      switch_direction == "No Switch" ~ "No Switch",
      is.na(fu_date2) ~ "No TB diagnosis",
      rel_days_to_tb < -365 ~ "≥1 year before TB",
      rel_days_to_tb >= -365 & rel_days_to_tb < -30 ~ "≥30 days before TB",
      between(rel_days_to_tb, -30, 30) ~ "±30 days around TB",
      TRUE ~ "Other" # Catch-all for unexpected cases
    )
  ) %>%
  left_join(select(widesymptoms_clean, id, TBincidence1, TBincidence2), by = "id")

# STEP 4: Tabulate results
switch_table <- switch_summary %>%
  count(TBincidence1, switch_direction, switch_timing_class) %>%
  pivot_wider(names_from = switch_timing_class, values_from = n, values_fill = 0)

# Check if any switches occurred after TB incidence

switch_after_tb <- switch_summary %>%
  filter(
    TBincidence1 == 1, # Have TB incidence
    switch_direction != "No Switch", # Had a switch
    !is.na(rel_days_to_tb), # Valid relative days
    rel_days_to_tb > 0
  ) # Switch after TB diagnosis

# Check how many such cases exist
n_after_tb <- nrow(switch_after_tb)
n_after_tb


#### Data setup ================================================================

# Analyze without individuals who have intermittent missingness
clean_ids <- count_clean$record_id

# Step 1: Pivot symptoms into long format
symptom_long2 <- widesymptoms_clean %>%
  select(record_id, starts_with("symptoms_"), starts_with("TBincidence")) %>%
  pivot_longer(
    cols = matches("^symptoms_\\d+"),
    names_to = "visit",
    values_to = "symptom_status"
  ) %>%
  mutate(
    visit_num = as.integer(str_extract(visit, "\\d+")),
    symptom_binary = case_when(
      symptom_status == "Any" ~ 1,
      symptom_status == "None" ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  arrange(record_id, visit_num)

# Step 2: Group and create lagged symptoms
symptom_long2 <- symptom_long2 %>%
  group_by(record_id) %>%
  mutate(
    symptom_lag = lag(symptom_binary),
    transition = case_when(
      symptom_lag == 0 & symptom_binary == 0 ~ "A-A",
      symptom_lag == 0 & symptom_binary == 1 ~ "A-S",
      symptom_lag == 1 & symptom_binary == 0 ~ "S-A",
      symptom_lag == 1 & symptom_binary == 1 ~ "S-S",
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup()

symptom_long2 <- symptom_long2 %>%
  filter(record_id %in% clean_ids)

tb_long <- widesymptoms_clean %>%
  select(record_id, matches("^tbfu\\d+_1$")) %>%
  pivot_longer(
    cols = -record_id,
    names_to = "tb_visit",
    values_to = "tb_dx"
  ) %>%
  mutate(
    visit_num = as.integer(str_extract(tb_visit, "\\d+"))
  )

symptom_long2 <- symptom_long2 %>%
  left_join(tb_long, by = c("record_id", "visit_num"))

model_data <- symptom_long2 %>%
  mutate(
    transition = factor(transition, levels = c(
      "A-A", "A-S", "S-A", "S-S"
    )),
    record_id = as.integer(record_id)
  )

# Add Longitudinal Lunit Scores (Time-Varying Covariate)
lunit_long <- widesymptoms_clean %>%
  select(record_id, LunitTB, starts_with("LunitTB_")) %>%
  mutate(LunitTB_0 = LunitTB) %>%
  select(-LunitTB) %>%
  pivot_longer(
    cols = starts_with("LunitTB_"),
    names_to = "visit",
    values_to = "LunitTB_score"
  ) %>%
  mutate(
    visit_num = as.integer(str_extract(visit, "\\d+")),
    record_id = as.integer(record_id)
  ) %>%
  select(record_id, visit_num, LunitTB_score)

model_data <- model_data %>%
  left_join(lunit_long, by = c("record_id", "visit_num")) %>%
  rename(LunitTB_end = LunitTB_score) %>%
  mutate(visit_num_start = visit_num - 1) %>%
  left_join(lunit_long, by = c("record_id", "visit_num_start" = "visit_num")) %>%
  rename(LunitTB_start = LunitTB_score)

covariates <- widesymptoms_clean %>%
  select(record_id, age, prison, prev_tb, smoking, any_drugs, prison_time, race, education, symptomatic_cellmates, tb_contact, LunitTB) %>%
  mutate(record_id = as.integer(record_id))

model_data <- model_data %>%
  merge(covariates, by = "record_id")

model_data <- model_data %>%
  mutate(record_id = as.integer(as.character(record_id))) %>%
  filter(!is.na(transition))

# Step A: Pivot all individual symptoms + sum_symptoms into long
symptom_details_long <- widesymptoms_clean %>%
  # include baseline (visit 0) symptoms too
  select(
    record_id,
    cough, sputum, blood, fever, appetite,
    weightloss, sweats, chestpain, breathing, bloodysputum, sum_symptoms,
    starts_with("cough_"), starts_with("sputum_"),
    starts_with("blood_"), starts_with("fever_"), starts_with("appetite_"),
    starts_with("weightloss_"), starts_with("sweats_"), starts_with("chestpain_"),
    starts_with("breathing_"), starts_with("bloodysputum_"),
    starts_with("sum_symptoms_")
  ) %>%
  # rename baseline to *_0 so they align with visit_num logic
  rename_with(
    ~ paste0(.x, "_0"),
    c(
      cough, sputum, blood, fever, appetite,
      weightloss, sweats, chestpain, breathing, bloodysputum, sum_symptoms
    )
  ) %>%
  # standardize labels
  mutate(across(-record_id, ~ as_factor(.) %>% as.character())) %>%
  pivot_longer(
    cols = -record_id,
    names_to = c("symptom", "visit_num"),
    names_pattern = "(.*)_(\\d+)$",
    values_to = "symptom_value"
  ) %>%
  mutate(
    visit_num = as.integer(visit_num),
    record_id = as.integer(record_id)
  ) %>%
  pivot_wider(
    names_from = symptom,
    values_from = symptom_value
  )
