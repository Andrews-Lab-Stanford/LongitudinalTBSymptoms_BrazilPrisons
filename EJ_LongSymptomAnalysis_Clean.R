################################################################################
#                         LONGITUDINAL SYMPTOMS ANALYSES
# Last Updated: August 12
# By: Esther Jung
################################################################################

#### Load Libraries and Data =====
library(lme4)
library(geepack)
library(rstanarm)
library(gtsummary)
library(sandwich)
library(lmtest)
library(dplyr)
library(broom)

# load("data/symptom_details_long1.RData")
# load("data/model_data1.RData")

table(model_data$TBincidence1)
table(model_data$tb_dx)

#### Table 1 ===================================================================

table1 <- model_data %>%
  mutate(
    TBincidence1 = factor(TBincidence1, levels = c("No", "Yes"), labels = c("No TB", "TB"))
  )

table1_baseline <- table1 %>%
  group_by(record_id) %>%
  slice(1) %>%
  ungroup()

table1_baseline <- table1_baseline %>%
  mutate(
    smoking_current = ifelse(smoking == "Current Smoker", "Yes", "No") %>% factor(),
    LunitTB = ifelse(LunitTB < 50, "Normal", "Abnormal") %>% factor()
  )

table1 <- table1_baseline %>%
  select(tb_dx, race, transition, age, prison, prev_tb, smoking_current, any_drugs, prison_time, education, symptomatic_cellmates, tb_contact, LunitTB) %>%
  tbl_summary(
    by = transition,
    type = list(
      prison ~ "categorical",
      prev_tb ~ "categorical",
      any_drugs ~ "categorical",
      smoking_current ~ "categorical",
      age ~ "continuous",
      prison_time ~ "continuous",
      tb_dx ~ "categorical",
      race ~ "categorical",
      education ~ "categorical",
      symptomatic_cellmates ~ "categorical",
      tb_contact ~ "categorical",
      LunitTB ~ "categorical"
    )
  ) %>%
  add_p(
    test = where(is.factor) ~ "chisq.test"
  ) %>%
  add_overall()

kruskal.test(age ~ transition, table1_baseline)
kruskal.test(prison_time ~ transition, table1_baseline)

chisq.test(table1_baseline$race, table1_baseline$transition)
chisq.test(table1_baseline$LunitTB, table1_baseline$transition)
chisq.test(table1_baseline$prison, table1_baseline$transition)
chisq.test(table1_baseline$education, table1_baseline$transition)
chisq.test(table1_baseline$prev_tb, table1_baseline$transition)
chisq.test(table1_baseline$tb_contact, table1_baseline$transition)
chisq.test(table1_baseline$smoking, table1_baseline$transition)
chisq.test(table1_baseline$any_drugs, table1_baseline$transition)
chisq.test(table1_baseline$symptomatic_cellmates, table1_baseline$transition)

chisq.test(table1_baseline$tb_dx, table1_baseline$transition)
table(table1_baseline$tb_dx, table1_baseline$transition)

library(gtsummary)
library(flextable)
library(officer)

table1_flex <- as_flex_table(table1)
doc <- read_docx() %>%
  body_add_flextable(table1_flex) %>%
  body_add_par(" ")
print(doc, target = "table1.docx")

#### Predictors of the change ==================================================
library(dplyr)
library(sandwich)
library(lmtest)
library(broom)

## 1. Predictors of becoming symptomatic (A-S vs A-A)
asymp_data <- model_data %>%
  filter(symptom_lag == 0) %>% # started asymptomatic
  mutate(
    change_binary = ifelse(transition == "A-S", 1, 0),
    LunitTB = ifelse(LunitTB < 50, "Normal", "Abnormal"),
    LunitTB = factor(LunitTB, levels = c("Normal", "Abnormal"))
  )

model_asymp <- glm(
  change_binary ~ age + prison + prev_tb + smoking + any_drugs +
    prison_time + symptomatic_cellmates + tb_contact + LunitTB,
  data = asymp_data,
  family = poisson(link = "log")
)

# Cluster-robust SEs
vcov_asymp <- vcovCL(model_asymp, cluster = ~record_id)
robust_asymp <- coeftest(model_asymp, vcov = vcov_asymp)

# Tidy results into a table
table_asymp <- broom::tidy(robust_asymp) %>%
  mutate(
    RR = exp(estimate),
    lower_CI = exp(estimate - 1.96 * std.error),
    upper_CI = exp(estimate + 1.96 * std.error),
    `95% CI` = paste0(round(lower_CI, 2), "-", round(upper_CI, 2)),
    RR = round(RR, 2),
    p.value = round(p.value, 3)
  ) %>%
  select(term, RR, `95% CI`, p.value)

## 2. Predictors of becoming asymptomatic (S-A vs S-S)
symp_data <- model_data %>%
  filter(symptom_lag == 1) %>% # started symptomatic
  mutate(
    change_binary = ifelse(transition == "S-A", 1, 0),
    LunitTB = ifelse(LunitTB < 50, "Normal", "Abnormal"),
    LunitTB = factor(LunitTB, levels = c("Normal", "Abnormal"))
  )

model_symp <- glm(
  change_binary ~ age + prison + prev_tb + smoking + any_drugs +
    prison_time + symptomatic_cellmates + tb_contact + LunitTB,
  data = symp_data,
  family = poisson(link = "log")
)

# Cluster-robust SEs
vcov_symp <- vcovCL(model_symp, cluster = ~record_id)
robust_symp <- coeftest(model_symp, vcov = vcov_symp)

# Tidy results into a table
table_symp <- broom::tidy(robust_symp) %>%
  mutate(
    RR = exp(estimate),
    lower_CI = exp(estimate - 1.96 * std.error),
    upper_CI = exp(estimate + 1.96 * std.error),
    `95% CI` = paste0(round(lower_CI, 2), "-", round(upper_CI, 2)),
    RR = round(RR, 2),
    p.value = round(p.value, 3)
  ) %>%
  select(term, RR, `95% CI`, p.value)

# Output
table_asymp
table_symp

#### Symptoms ==================================================================

# Recode symptoms to binary inside symptom_details_long
symptom_details_long <- symptom_details_long %>%
  mutate(across(
    c(
      cough, sputum, blood, fever, appetite,
      weightloss, sweats, chestpain, breathing, bloodysputum
    ),
    ~ case_when(
      . %in% c("Yes", "Any", 1) ~ 1,
      . %in% c("No", "None", 0) ~ 0,
      TRUE ~ NA_real_
    )
  )) %>%
  mutate(sum_symptoms = as.numeric(sum_symptoms)) %>%
  group_by(record_id) %>%
  arrange(visit_num, .by_group = TRUE) %>%
  # Lag variables (baseline = visit 0 now included!)
  mutate(across(
    c(
      cough, sputum, blood, fever, appetite,
      weightloss, sweats, chestpain, breathing, bloodysputum
    ),
    list(
      lag = ~ lag(.),
      change = ~ case_when(
        lag(.) == 0 & . == 1 ~ "Onset",
        lag(.) == 1 & . == 0 ~ "Resolution",
        lag(.) == 1 & . == 1 ~ "Persist",
        lag(.) == 0 & . == 0 ~ "Absent",
        TRUE ~ NA_character_
      )
    ),
    .names = "{.col}_{.fn}"
  )) %>%
  mutate(
    sum_symptoms_lag = lag(sum_symptoms),
    sum_symptoms_change = sum_symptoms - sum_symptoms_lag
  ) %>%
  ungroup()

transition_lookup <- model_data %>%
  select(record_id, visit_num, transition)

# Symptom change summary stratified by transition
symptom_change_summary <- symptom_details_long %>%
  filter(visit_num > 0) %>%
  select(record_id, visit_num, ends_with("_change")) %>%
  select(-sum_symptoms_change) %>%
  pivot_longer(
    cols = ends_with("_change"),
    names_to = "symptom",
    values_to = "change"
  ) %>%
  left_join(transition_lookup, by = c("record_id", "visit_num")) %>%
  group_by(transition, symptom, change) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(transition, symptom, change) %>%
  filter(!is.na(change))

denoms <- tibble::tibble(
  transition = c(
    "A-S",
    "S-S",
    "S-A"
  ),
  denom = c(1092, 755, 1025)
)

symptom_change_pct_clean <- symptom_change_summary %>%
  filter(change %in% c("Onset", "Resolution")) %>%
  group_by(transition, symptom) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  left_join(denoms, by = "transition") %>%
  mutate(pct = 100 * n / denom) %>%
  ungroup() %>%
  # clean labels
  mutate(symptom = gsub("_change", "", symptom)) %>%
  filter(symptom != "blood") %>%
  filter(!is.na(transition)) %>%
  mutate(
    symptom = recode(symptom,
      weightloss = "Weight loss",
      sweats = "Night sweats",
      fever = "Fever",
      cough = "Cough",
      sputum = "Productive cough",
      chestpain = "Chest pain",
      breathing = "Breathing trouble",
      bloodysputum = "Blood in sputum",
      appetite = "Appetite"
    )
  )
symptom_change_pct_clean

symptom_order <- c(
  "Blood in sputum",
  "Breathing trouble",
  "Cough",
  "Productive cough",
  "Chest pain",
  "Fever",
  "Night sweats",
  "Weight loss",
  "Appetite"
)

transition_order <- c(
  "S-A",
  "S-S",
  "A-S"
)

symptom_change_pct_clean <- symptom_change_pct_clean %>%
  mutate(symptom = factor(symptom, levels = symptom_order)) %>%
  mutate(transition = factor(transition, levels = transition_order))

symptom_heatmap <- ggplot(symptom_change_pct_clean, aes(x = symptom, y = transition, fill = pct)) +
  geom_tile(color = "white", ) +
  geom_text(
    aes(label = ifelse(pct >= 1, paste0(round(pct), "%"), "")),
    color = "black", size = 3, fontface = "bold"
  ) +
  scale_fill_gradient2(
    low = "white", mid = "lightsalmon", high = "firebrick3",
    midpoint = 30, limits = c(0, max(symptom_change_pct_clean$pct, na.rm = TRUE))
  ) +
  scale_y_discrete(labels = c(
    "S-A" = "Symptomatic to Asymptomatic",
    "S-S" = "Symptomatic to Symptomatic",
    "A-S" = "Asymptomatic to Symptomatic"
  )) +
  labs(
    x = "Symptom",
    y = "Symptom transition",
    fill = "% change"
  ) +
  coord_fixed(ratio = 0.8) +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(color = "black"), # all text black
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    legend.text = element_text(color = "black", size = 14),
    legend.title = element_text(color = "black", size = 14),
    panel.grid = element_blank()
  )
symptom_heatmap

ggsave(
  filename = "~/Dropbox/TB-R01-Brazil/figures/longsymptoms_symptomheatmap.png",
  plot = symptom_heatmap,
  width = 10, # inches
  height = 4, # inches
  dpi = 800
)


# Symptom count change stratified by transition
sum_change_summary <- symptom_details_long %>%
  filter(visit_num > 0) %>%
  left_join(transition_lookup, by = c("record_id", "visit_num")) %>%
  group_by(transition) %>%
  summarise(
    avg_change = mean(sum_symptoms_change, na.rm = TRUE),
    median_change = median(sum_symptoms_change, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

table(model_data$transition)

#### Main: Modified Poisson model ====================================================

table(model_data$transition)
table(model_data$transition, model_data$tb_dx)

poisson_model <- glm(
  tb_dx ~ transition + age + prison + prev_tb + smoking + any_drugs + prison_time + symptomatic_cellmates + tb_contact + LunitTB,
  data = model_data,
  family = poisson(link = "log")
)

# Compute cluster-robust standard errors (clustered by individual ID)
clustered_se <- vcovCL(poisson_model, cluster = ~record_id)

# Tidy the model with robust SEs
tidy_results <- coeftest(poisson_model, vcov = clustered_se) %>%
  broom::tidy() %>%
  mutate(
    RR = exp(estimate),
    lower_CI = exp(estimate - 1.96 * std.error),
    upper_CI = exp(estimate + 1.96 * std.error)
  )

print(tidy_results)

# Visualize
plot_data <- tidy_results %>%
  filter(str_detect(term, "transition")) %>%
  mutate(
    transition = case_when(
      term == "transitionA-S" ~ "Asymptomatic → Symptomatic",
      term == "transitionS-S" ~ "Symptomatic → Symptomatic",
      term == "transitionS-A" ~ "Symptomatic → Asymptomatic",
      TRUE ~ term
    ),
    transition = factor(transition, levels = c("Asymptomatic → Symptomatic", "Symptomatic → Symptomatic", "Symptomatic → Asymptomatic"))
  )

dummy_row <- data.frame(
  term       = NA,
  estimate   = NA,
  std.error  = NA,
  statistic  = NA,
  p.value    = NA,
  RR         = NA,
  lower_CI   = NA,
  upper_CI   = NA,
  transition = "Reference"
)

plot_data_dummy <- rbind(dummy_row, plot_data)

plot_data_dummy$transition <- factor(
  plot_data_dummy$transition,
  levels = c(
    "Symptomatic → Asymptomatic",
    "Symptomatic → Symptomatic",
    "Asymptomatic → Symptomatic",
    "Reference"
  )
)
dark_colors <- c(
  "Asymptomatic → Symptomatic" = "#1C7E50",
  "Symptomatic → Symptomatic" = "#7570B3",
  "Symptomatic → Asymptomatic" = "#D95F02"
)

plot_tb_main <- ggplot(plot_data_dummy, aes(x = RR, y = transition, color = transition)) +
  geom_point(shape = 15, size = 3, na.rm = TRUE) +
  geom_errorbarh(aes(xmin = lower_CI, xmax = upper_CI), height = 0.2, na.rm = TRUE) +
  geom_vline(xintercept = 1, color = "gray40", alpha = 0.3, linetype = "dashed") +
  scale_x_continuous(
    limits = c(0.5, 3.5),
    breaks = seq(0.5, 3.5, by = .5) # set tick positions
  ) +
  labs(x = "Adjusted Risk Ratio (95% CI)", y = NULL) +
  scale_color_manual(values = dark_colors) + # apply darker colors
  theme_minimal(base_size = 12, base_family = "Arial") +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.length = unit(2, "mm"),
    legend.position = "none"
  )

plot_tb_main

ggsave(
  filename = "~/Dropbox/TB-R01-Brazil/figures/longsymptoms_forestplot.png",
  plot = plot_tb_main,
  width = 5, # inches
  height = 2, # inches
  dpi = 1000
)

#### Kaplan Meier Curves =======================================================

baseline_data <- model_data %>%
  group_by(record_id) %>%
  slice(1) %>% # take first row of each record_id
  transmute(
    record_id,
    TBincidence1,
    TBincidence2,
    visit = "symptoms_0",
    symptom_status = ifelse(symptom_lag == 1, "Any", "None"),
    visit_num = 0,
    symptom_binary = symptom_lag, # baseline = lagged symptom status
    symptom_lag = NA,
    transition = NA
  ) %>%
  mutate(record_id = as.character(record_id))

model_data_long <- bind_rows(
  baseline_data %>% mutate(across(everything(), as.character)),
  model_data %>% mutate(across(everything(), as.character))
) %>%
  arrange(record_id, visit_num) %>%
  mutate(
    visit_num = as.numeric(visit_num),
    symptom_binary = as.numeric(symptom_binary)
  )

baseline_status <- model_data_long %>%
  filter(visit_num == 0) %>%
  select(record_id, baseline_symptom = symptom_binary)

long_with_baseline <- model_data_long %>%
  left_join(baseline_status, by = "record_id") %>%
  arrange(record_id, visit_num) %>%
  mutate(switched = ifelse(symptom_binary != baseline_symptom, 1, 0))

first_switch <- long_with_baseline %>%
  group_by(record_id) %>%
  summarise(
    baseline_symptom = first(baseline_symptom),
    first_change_visit = visit_num[switched == 1][1],
    n_visits = length(visit_num),
  )

visit_props <- first_switch %>%
  group_by(baseline_symptom) %>%
  do({
    bs <- .
    max_visit <- max(long_with_baseline$visit_num, na.rm = TRUE)
    tibble(visit_num = 0:max_visit) %>%
      rowwise() %>%
      mutate(
        total_pop = sum(bs$n_visits > visit_num),
        n_at_risk = sum((is.na(bs$first_change_visit) | bs$first_change_visit > visit_num) & bs$n_visits > visit_num),
        n_events = ifelse(visit_num == 0, 0, sum(bs$first_change_visit == visit_num, na.rm = TRUE)),
        prop_event = ifelse(visit_num == 0, 0, n_events / (n_events + n_at_risk)),
        sp = ifelse(visit_num == 0, 1, (n_at_risk) / (total_pop))
      )
  }) %>%
  ungroup()

km_curves <- visit_props %>%
  group_by(baseline_symptom) %>%
  arrange(visit_num) %>%
  mutate(
    surv_prob = cumprod(1 - prop_event)
  ) %>%
  ungroup()

# Plot for asymptomatic
p_asym <- ggplot(
  km_curves %>% filter(baseline_symptom == 0),
  aes(x = visit_num, y = surv_prob)
) +
  geom_step(color = "blue3", size = 1.5) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_x_continuous(breaks = unique(km_curves$visit_num)) +
  labs(
    title = "Started Asymptomatic",
    x = "Visit Number",
    y = "Probability of Remaining Asymptomatic"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black")
  )

# How many people in the cohort started asymptomatic and developed new symptoms at any point?
ever_symptomatic <- long_with_baseline %>%
  group_by(record_id) %>%
  summarise(
    baseline_symptom = first(baseline_symptom),
    ever_symptom = any(symptom_binary == 1, na.rm = TRUE) # ever had symptoms
  ) %>%
  filter(baseline_symptom == 0) # only those starting asymptomatic

# Count and proportion
table_asym <- ever_symptomatic %>%
  summarise(
    n_asym = n(),
    n_developed = sum(ever_symptom),
    prop_developed = mean(ever_symptom)
  )

table_asym

# Plot for symptomatic
p_sym <- ggplot(
  km_curves %>% filter(baseline_symptom == 1),
  aes(x = visit_num, y = surv_prob)
) +
  geom_step(color = "red3", size = 1.5) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_x_continuous(breaks = unique(km_curves$visit_num)) +
  labs(
    title = "Started Symptomatic",
    x = "Visit Number",
    y = "Probability of Remaining Symptomatic"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black")
  )

library(cowplot)
combined <- plot_grid(p_asym, p_sym, ncol = 2)
combined
ggsave(
  filename = "~/Dropbox/TB-R01-Brazil/figures/longsymptoms_kmcurves.png",
  plot = combined,
  width = 8, height = 4, dpi = 1000
)

km_table <- km_curves %>%
  arrange(baseline_symptom, visit_num) %>%
  mutate(
    surv_prob = round(surv_prob, 3),
    prop_event = round(prop_event, 3)
  )

km_table

#### Next transition ===========================================================

post_as <- model_data %>%
  group_by(record_id) %>%
  arrange(visit_num, .by_group = TRUE) %>%
  mutate(next_transition = lead(transition)) %>%
  filter(transition == "A-S" & !is.na(next_transition)) %>%
  ungroup() %>%
  count(next_transition) %>%
  mutate(
    prop = n / sum(n),
    prior_transition = "A-S"
  )

post_sa <- model_data %>%
  group_by(record_id) %>%
  arrange(visit_num, .by_group = TRUE) %>%
  mutate(next_transition = lead(transition)) %>%
  filter(transition == "S-A" & !is.na(next_transition)) %>%
  ungroup() %>%
  count(next_transition) %>%
  mutate(
    prop = n / sum(n),
    prior_transition = "S-A"
  )

post_combined <- bind_rows(post_as, post_sa)

# Plot together
ggplot(post_combined, aes(x = next_transition, y = prop, fill = next_transition)) +
  geom_col() +
  facet_wrap(~prior_transition) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Next Transitions After A-S and S-A",
    x = "Next Transition",
    y = "Proportion"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

#### Other Poisson models ====================================================

model_data_smoke <- model_data %>%
  mutate(smoking = case_when(
    smoking == "Never Smoker" ~ "Never",
    TRUE ~ "Ever"
  )) %>%
  mutate(smoking = factor(smoking, levels = c("Never", "Ever")))

#### Interaction with smoking
poisson_model_smoking <- glm(
  tb_dx ~ transition * smoking + age + prison + prev_tb + any_drugs + prison_time,
  data = model_data_smoke,
  family = poisson(link = "log")
)

clustered_se_smoking <- vcovCL(poisson_model_smoking, cluster = ~record_id)

tidy_smoking <- coeftest(poisson_model_smoking, vcov = clustered_se_smoking) %>%
  broom::tidy() %>%
  mutate(
    RR = exp(estimate),
    lower_CI = exp(estimate - 1.96 * std.error),
    upper_CI = exp(estimate + 1.96 * std.error)
  )

print(tidy_smoking)

#### Interaction with drugs
poisson_model_drugs <- glm(
  tb_dx ~ transition * any_drugs + age + prison + prev_tb + smoking + prison_time,
  data = model_data,
  family = poisson(link = "log")
)

# Cluster-robust SEs
clustered_se <- vcovCL(poisson_model_drugs, cluster = ~record_id)

# Tidy output while keeping term names
tidy_drugs <- tidy(poisson_model_drugs, conf.int = TRUE) %>%
  mutate(
    std.error = sqrt(diag(clustered_se)),
    estimate = coef(poisson_model_drugs),
    statistic = estimate / std.error,
    p.value = 2 * (1 - pnorm(abs(statistic))),
    RR = exp(estimate),
    lower_CI = exp(estimate - 1.96 * std.error),
    upper_CI = exp(estimate + 1.96 * std.error)
  )

print(tidy_drugs)

#### Interaction: transition × prev_tb
poisson_model_tb <- glm(
  tb_dx ~ transition * prev_tb + age + prison + smoking + any_drugs + prison_time,
  data = model_data,
  family = poisson(link = "log")
)

# Cluster-robust SEs
clustered_se <- vcovCL(poisson_model_tb, cluster = ~record_id)

# Tidy output while keeping term names
tidy_prevtb <- tidy(poisson_model_tb, conf.int = TRUE) %>%
  mutate(
    std.error = sqrt(diag(clustered_se)),
    estimate = coef(poisson_model_drugs),
    statistic = estimate / std.error,
    p.value = 2 * (1 - pnorm(abs(statistic))),
    RR = exp(estimate),
    lower_CI = exp(estimate - 1.96 * std.error),
    upper_CI = exp(estimate + 1.96 * std.error)
  )

print(tidy_prevtb)

#### LME4 GLMM Logistic ========================================================
# # Crude
# summary(glmer(tb_dx ~ transition + (1|record_id), data = model_data, family = binomial()))
#
# # Adjusted
# model_mixed <- glmer(
#   tb_dx ~ transition + age + prison + prev_tb + smoking + any_drugs + prison_time + (1 | record_id),
#   data = model_data_complete,
#   family = binomial()
# )
# summary(model_mixed)
# exp(cbind(OR = fixef(model_mixed), confint(model_mixed, method = "Wald")))
#
# model_data_complete <- model_data %>%
#   select(tb_dx, transition, age, prison, prev_tb, smoking, any_drugs, prison_time, record_id) %>%
#   na.omit() %>%
#   mutate(
#     transition = factor(transition),
#     prison = factor(prison),
#     prev_tb = factor(prev_tb),
#     smoking = factor(smoking),
#     any_drugs = factor(any_drugs)
#   )

#### Bayesian GLMM Logistic ====================================================

# model_stan <- stan_glmer(
#   tb_dx ~ transition + age + prison + prev_tb + smoking + any_drugs + prison_time + (1 | record_id),
#   data = model_data_complete,
#   family = binomial(link = "logit"),
#   prior = normal(0, 2.5),
#   prior_intercept = normal(0, 5),
#   chains = 4,
#   iter = 2000,
#   cores = 4,
#   seed = 123
# # )
#
# library(broom.mixed)
# library(dplyr)
#
# # Extract tidy summary of fixed effects
# tidy_results <- tidy(model_stan, effects = "fixed")
#
# # Compute OR, Wald 95% CI, and pseudo p-value
# tidy_results <- tidy_results %>%
#   mutate(
#     p_value = 2 * (1 - pnorm(abs(estimate / std.error))),  # Two-tailed Wald p-value
#     OR = exp(estimate),
#     lower_CrI = exp(estimate - 1.96 * std.error),  # Wald CI lower
#     upper_CrI = exp(estimate + 1.96 * std.error)   # Wald CI upper
#   )
#
# print(tidy_results)
#### Added Feb 2026: Model adjusted for LunitTB at start of interval (time-varying) ====

# Ensure libraries are loaded
library(sandwich)
library(lmtest)
library(broom)
library(dplyr)

print("Running Poisson model with adjustment for LunitTB_start (time-varying X-ray)...")

# Create categorical version of time-varying X-ray (Cutoff 50)
model_data <- model_data %>%
  mutate(
    LunitTB_start_cat = case_when(
      LunitTB_start >= 50 ~ "Abnormal",
      LunitTB_start < 50 ~ "Normal",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("Normal", "Abnormal"))
  )

# Poisson model replacing baseline 'LunitTB' with 'LunitTB_start_cat'
poisson_model_tv_xray <- glm(
  tb_dx ~ transition + age + prison + prev_tb + smoking + any_drugs + prison_time + symptomatic_cellmates + tb_contact + LunitTB_start_cat,
  data = model_data,
  family = poisson(link = "log")
)

# Compute cluster-robust standard errors
clustered_se_tv <- vcovCL(poisson_model_tv_xray, cluster = ~record_id)

# Tidy results
tidy_results_tv <- coeftest(poisson_model_tv_xray, vcov = clustered_se_tv) %>%
  broom::tidy() %>%
  mutate(
    RR = exp(estimate),
    lower_CI = exp(estimate - 1.96 * std.error),
    upper_CI = exp(estimate + 1.96 * std.error)
  )

print("Results for Time-Varying X-ray Adjustment (Categorized Normal/Abnormal):")
print(tidy_results_tv)

# Visualize Time-Varying X-ray Adjusted Model
library(ggplot2)
library(stringr)

plot_data_tv <- tidy_results_tv %>%
  filter(str_detect(term, "transition")) %>%
  mutate(
    transition = case_when(
      term == "transitionA-S" ~ "Asymptomatic → Symptomatic",
      term == "transitionS-S" ~ "Symptomatic → Symptomatic",
      term == "transitionS-A" ~ "Symptomatic → Asymptomatic",
      TRUE ~ term
    ),
    transition = factor(transition, levels = c("Asymptomatic → Symptomatic", "Symptomatic → Symptomatic", "Symptomatic → Asymptomatic"))
  )

dummy_row <- data.frame(
  term       = NA,
  estimate   = NA,
  std.error  = NA,
  statistic  = NA,
  p.value    = NA,
  RR         = NA,
  lower_CI   = NA,
  upper_CI   = NA,
  transition = "Reference"
)

plot_data_dummy_tv <- rbind(dummy_row, plot_data_tv)

plot_data_dummy_tv$transition <- factor(
  plot_data_dummy_tv$transition,
  levels = c(
    "Symptomatic → Asymptomatic",
    "Symptomatic → Symptomatic",
    "Asymptomatic → Symptomatic",
    "Reference"
  )
)

dark_colors <- c(
  "Asymptomatic → Symptomatic" = "#1C7E50",
  "Symptomatic → Symptomatic" = "#7570B3",
  "Symptomatic → Asymptomatic" = "#D95F02"
)

plot_tb_main_tv <- ggplot(plot_data_dummy_tv, aes(x = RR, y = transition, color = transition)) +
  geom_point(shape = 15, size = 3, na.rm = TRUE) +
  geom_errorbarh(aes(xmin = lower_CI, xmax = upper_CI), height = 0.2, na.rm = TRUE) +
  geom_vline(xintercept = 1, color = "gray40", alpha = 0.3, linetype = "dashed") +
  scale_x_continuous(
    limits = c(0.5, 3.5),
    breaks = seq(0.5, 3.5, by = .5) # set tick positions
  ) +
  labs(x = "Adjusted Risk Ratio (95% CI)", y = NULL) +
  scale_color_manual(values = dark_colors) + # apply darker colors
  theme_minimal(base_size = 12, base_family = "Arial") +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.length = unit(2, "mm"),
    legend.position = "none"
  )

print(plot_tb_main_tv)

ggsave(
  filename = "~/Dropbox/TB-R01-Brazil/figures/longsymptoms_forestplot_TV_adj.png",
  plot = plot_tb_main_tv,
  width = 5, # inches
  height = 2, # inches
  dpi = 1000
)
