################################################################################
#         LONGITUDINAL SYMPTOMS ANALYSES - STRATIFIED BY X-RAY STATUS
# Reviewer Request: Stratify by chest radiographic status at first visit
#                   of each 4-month interval
# Last Updated: [Current Date]
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
library(ggplot2)
library(stringr)
library(tidyr)

#### Prepare Data with Radiographic Status at Interval Start =====

# The radiographic status at the "first visit of the given 4-month interval"
# corresponds to the visit where symptom_lag is defined (the start of the interval)
# For each transition from visit N-1 to visit N, we need LunitTB at visit N-1

# Check if LunitTB_start is already in model_data (added by cleaning script)
if ("LunitTB_start" %in% names(model_data)) {
  cat("Using pre-calculated LunitTB_start from model_data\n")
  model_data_stratified <- model_data %>%
    mutate(
      LunitTB_for_strat = ifelse(!is.na(LunitTB_start), LunitTB_start, LunitTB)
    )
} else {
  # Fallback: Attempt to extract visit-specific LunitTB data if available
  # Check if we can load the raw data to get longitudinal LunitTB
  tryCatch(
    {
      # Try to load raw data - adjust path/object name as needed
      # This assumes the raw data object might be available
      if (exists("definition1cohort") || exists("widesymptoms")) {
        # Extract longitudinal LunitTB data
        if (exists("definition1cohort")) {
          raw_data <- definition1cohort
        } else {
          raw_data <- widesymptoms
        }

        # Check if visit-specific LunitTB columns exist
        lunit_cols <- names(raw_data)[grepl("^LunitTB", names(raw_data))]

        if (length(lunit_cols) > 1) {
          # Create longitudinal LunitTB dataset
          lunit_long <- raw_data %>%
            select(record_id, all_of(lunit_cols)) %>%
            pivot_longer(
              cols = starts_with("LunitTB"),
              names_to = "visit",
              values_to = "LunitTB_score"
            ) %>%
            mutate(
              visit_num = case_when(
                visit == "LunitTB" ~ 0L,
                TRUE ~ as.integer(str_extract(visit, "\\d+"))
              )
            ) %>%
            filter(!is.na(visit_num)) %>%
            select(record_id, visit_num, LunitTB_score) %>%
            mutate(record_id = as.integer(record_id))

          # Merge with model_data to get LunitTB at the start of each interval
          # For transition at visit N, we need LunitTB at visit N-1 (where symptom_lag is from)
          model_data_with_lunit <- model_data %>%
            mutate(
              # The visit where symptom_lag is defined (start of interval)
              interval_start_visit = visit_num - 1
            ) %>%
            left_join(
              lunit_long %>%
                rename(LunitTB_interval_start = LunitTB_score),
              by = c("record_id", "interval_start_visit" = "visit_num")
            ) %>%
            # Use interval-start LunitTB if available, otherwise fall back to baseline
            mutate(
              LunitTB_for_stratification = ifelse(
                !is.na(LunitTB_interval_start),
                LunitTB_interval_start,
                LunitTB
              )
            )

          cat("Using visit-specific LunitTB data for stratification\n")
          model_data_stratified <- model_data_with_lunit
        } else {
          cat("Only baseline LunitTB available, using baseline for stratification\n")
          model_data_stratified <- model_data
        }
      } else {
        cat("Raw data not found, using baseline LunitTB from model_data\n")
        model_data_stratified <- model_data
      }
    },
    error = function(e) {
      cat("Could not extract visit-specific LunitTB, using baseline:\n", e$message, "\n")
      model_data_stratified <- model_data
    }
  )
}

# Create radiographic status variable for stratification
# Normal < 50, Abnormal >= 50
# Check if visit-specific LunitTB was successfully created
if ("LunitTB_for_stratification" %in% names(model_data_stratified)) {
  # Use visit-specific LunitTB at interval start (fallback path)
  model_data_stratified <- model_data_stratified %>%
    mutate(LunitTB_for_strat = LunitTB_for_stratification)
}
# If neither path set LunitTB_for_strat (e.g. fallback failed), use baseline
if (!"LunitTB_for_strat" %in% names(model_data_stratified)) {
  model_data_stratified <- model_data_stratified %>%
    mutate(LunitTB_for_strat = LunitTB)
}

# Convert to categorical
model_data_stratified <- model_data_stratified %>%
  mutate(
    # Convert to categorical (handle both numeric and character/factor inputs)
    LunitTB_cat = case_when(
      is.numeric(LunitTB_for_strat) & LunitTB_for_strat < 50 ~ "Normal",
      is.numeric(LunitTB_for_strat) & LunitTB_for_strat >= 50 ~ "Abnormal",
      is.character(LunitTB_for_strat) & LunitTB_for_strat == "Normal" ~ "Normal",
      is.character(LunitTB_for_strat) & LunitTB_for_strat == "Abnormal" ~ "Abnormal",
      is.factor(LunitTB_for_strat) & as.character(LunitTB_for_strat) == "Normal" ~ "Normal",
      is.factor(LunitTB_for_strat) & as.character(LunitTB_for_strat) == "Abnormal" ~ "Abnormal",
      TRUE ~ NA_character_
    ),
    LunitTB_cat = factor(LunitTB_cat, levels = c("Normal", "Abnormal")),
    # Ensure LunitTB variable is properly formatted for models (baseline/control)
    LunitTB = case_when(
      is.numeric(LunitTB) & LunitTB < 50 ~ "Normal",
      is.numeric(LunitTB) & LunitTB >= 50 ~ "Abnormal",
      TRUE ~ as.character(LunitTB)
    ),
    LunitTB = factor(LunitTB, levels = c("Normal", "Abnormal"))
  ) %>%
  filter(!is.na(LunitTB_cat)) # Remove missing radiographic status

# Check sample sizes by stratum
table(model_data_stratified$transition)
table(model_data_stratified$LunitTB_cat)
table(model_data_stratified$LunitTB_cat, model_data_stratified$transition)
table(model_data_stratified$LunitTB_cat, model_data_stratified$transition, model_data_stratified$tb_dx)

#### Stratified Main Analysis: Modified Poisson Model ========================

# Normal X-ray stratum
model_data_normal <- model_data_stratified %>%
  filter(LunitTB_cat == "Normal")

poisson_model_normal <- glm(
  tb_dx ~ transition + age + prison + prev_tb + smoking + any_drugs +
    prison_time + symptomatic_cellmates + tb_contact,
  data = model_data_normal,
  family = poisson(link = "log")
)

# Compute cluster-robust standard errors
clustered_se_normal <- vcovCL(poisson_model_normal, cluster = ~record_id)

# Tidy results with robust SEs
tidy_results_normal <- coeftest(poisson_model_normal, vcov = clustered_se_normal) %>%
  broom::tidy() %>%
  mutate(
    RR = exp(estimate),
    lower_CI = exp(estimate - 1.96 * std.error),
    upper_CI = exp(estimate + 1.96 * std.error),
    stratum = "Normal X-ray"
  )

print("=== Normal X-ray Stratum ===")
print(tidy_results_normal)

# Abnormal X-ray stratum
model_data_abnormal <- model_data_stratified %>%
  filter(LunitTB_cat == "Abnormal")

poisson_model_abnormal <- glm(
  tb_dx ~ transition + age + prison + prev_tb + smoking + any_drugs +
    prison_time + symptomatic_cellmates + tb_contact,
  data = model_data_abnormal,
  family = poisson(link = "log")
)

# Compute cluster-robust standard errors
clustered_se_abnormal <- vcovCL(poisson_model_abnormal, cluster = ~record_id)

# Tidy results with robust SEs
tidy_results_abnormal <- coeftest(poisson_model_abnormal, vcov = clustered_se_abnormal) %>%
  broom::tidy() %>%
  mutate(
    RR = exp(estimate),
    lower_CI = exp(estimate - 1.96 * std.error),
    upper_CI = exp(estimate + 1.96 * std.error),
    stratum = "Abnormal X-ray"
  )

print("=== Abnormal X-ray Stratum ===")
print(tidy_results_abnormal)

#### Combined Results Table ===================================================

# Combine results for comparison
combined_results <- bind_rows(
  tidy_results_normal,
  tidy_results_abnormal
) %>%
  filter(str_detect(term, "transition")) %>%
  mutate(
    transition = case_when(
      term == "transitionA-S" ~ "Asymptomatic → Symptomatic",
      term == "transitionS-S" ~ "Symptomatic → Symptomatic",
      term == "transitionS-A" ~ "Symptomatic → Asymptomatic",
      TRUE ~ term
    ),
    transition = factor(transition, levels = c(
      "Asymptomatic → Symptomatic",
      "Symptomatic → Symptomatic",
      "Symptomatic → Asymptomatic"
    ))
  )

# Create summary table
summary_table <- combined_results %>%
  select(stratum, transition, RR, lower_CI, upper_CI, p.value) %>%
  mutate(
    `95% CI` = paste0(round(lower_CI, 2), "-", round(upper_CI, 2)),
    RR = round(RR, 2),
    p.value = round(p.value, 3)
  ) %>%
  arrange(stratum, transition)

print("=== Combined Stratified Results ===")
print(summary_table)

#### Stratified Figure 1: Forest Plot ========================================

# Prepare plot data
plot_data_normal <- tidy_results_normal %>%
  filter(str_detect(term, "transition")) %>%
  mutate(
    transition = case_when(
      term == "transitionA-S" ~ "Asymptomatic → Symptomatic",
      term == "transitionS-S" ~ "Symptomatic → Symptomatic",
      term == "transitionS-A" ~ "Symptomatic → Asymptomatic",
      TRUE ~ term
    ),
    transition = factor(transition, levels = c(
      "Asymptomatic → Symptomatic",
      "Symptomatic → Symptomatic",
      "Symptomatic → Asymptomatic"
    ))
  )

plot_data_abnormal <- tidy_results_abnormal %>%
  filter(str_detect(term, "transition")) %>%
  mutate(
    transition = case_when(
      term == "transitionA-S" ~ "Asymptomatic → Symptomatic",
      term == "transitionS-S" ~ "Symptomatic → Symptomatic",
      term == "transitionS-A" ~ "Symptomatic → Asymptomatic",
      TRUE ~ term
    ),
    transition = factor(transition, levels = c(
      "Asymptomatic → Symptomatic",
      "Symptomatic → Symptomatic",
      "Symptomatic → Asymptomatic"
    ))
  )

# Add reference rows
dummy_row_normal <- data.frame(
  term       = NA,
  estimate   = NA,
  std.error  = NA,
  statistic  = NA,
  p.value    = NA,
  RR         = NA,
  lower_CI   = NA,
  upper_CI   = NA,
  stratum    = "Normal X-ray",
  transition = "Reference"
)

dummy_row_abnormal <- data.frame(
  term       = NA,
  estimate   = NA,
  std.error  = NA,
  statistic  = NA,
  p.value    = NA,
  RR         = NA,
  lower_CI   = NA,
  upper_CI   = NA,
  stratum    = "Abnormal X-ray",
  transition = "Reference"
)

plot_data_normal_dummy <- rbind(dummy_row_normal, plot_data_normal)
plot_data_abnormal_dummy <- rbind(dummy_row_abnormal, plot_data_abnormal)

# Combine for plotting
plot_data_combined <- bind_rows(
  plot_data_normal_dummy,
  plot_data_abnormal_dummy
) %>%
  mutate(
    transition = factor(
      transition,
      levels = c(
        "Symptomatic → Asymptomatic",
        "Symptomatic → Symptomatic",
        "Asymptomatic → Symptomatic",
        "Reference"
      )
    ),
    stratum = factor(stratum, levels = c("Normal X-ray", "Abnormal X-ray"))
  )

# Color scheme (matching original)
dark_colors <- c(
  "Asymptomatic → Symptomatic" = "#1C7E50",
  "Symptomatic → Symptomatic" = "#7570B3",
  "Symptomatic → Asymptomatic" = "#D95F02"
)

# Create stratified forest plot
plot_tb_stratified <- ggplot(
  plot_data_combined,
  aes(x = RR, y = transition, color = transition)
) +
  geom_point(shape = 15, size = 3, na.rm = TRUE) +
  geom_errorbarh(aes(xmin = lower_CI, xmax = upper_CI), height = 0.2, na.rm = TRUE) +
  geom_vline(xintercept = 1, color = "gray40", alpha = 0.3, linetype = "dashed") +
  facet_wrap(~stratum, ncol = 1) +
  scale_x_continuous(
    limits = c(0, 5.5),
    breaks = seq(0,5.5, by = .5)
  ) +
  labs(
    x = "Adjusted Risk Ratio (95% CI)",
    y = NULL,
  ) +
  scale_color_manual(values = dark_colors) +
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
    legend.position = "none",
    strip.text = element_text(size = 12),
  )

plot_tb_stratified

# Save plot
ggsave(
  filename = "~/Dropbox/TB-R01-Brazil/figures/longsymptoms_forestplot_stratified_xray.png",
  plot = plot_tb_stratified,
  width = 6, # inches
  height = 5, # inches
  dpi = 1000
)

#### Alternative: Side-by-side comparison plot ================================

# Create a version with side-by-side panels
plot_tb_stratified_side <- ggplot(
  plot_data_combined,
  aes(x = RR, y = transition, color = transition)
) +
  geom_point(shape = 15, size = 3, na.rm = TRUE) +
  geom_errorbarh(aes(xmin = lower_CI, xmax = upper_CI), height = 0.2, na.rm = TRUE) +
  geom_vline(xintercept = 1, color = "gray40", alpha = 0.3, linetype = "dashed") +
  facet_wrap(~stratum, ncol = 2) +
  scale_x_continuous(
    limits = c(0.5, 4),
    breaks = seq(0.5, 3.5, by = .5)
  ) +
  labs(
    x = "Adjusted Risk Ratio (95% CI)",
    y = NULL,
    title = "TB Risk by Symptom Transition, Stratified by Chest Radiographic Status"
  ) +
  scale_color_manual(values = dark_colors) +
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
    legend.position = "none",
    strip.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

plot_tb_stratified_side

# Save side-by-side plot
ggsave(
  filename = "~/Dropbox/TB-R01-Brazil/figures/longsymptoms_forestplot_stratified_xray_sidebyside.png",
  plot = plot_tb_stratified_side,
  width = 10, # inches
  height = 4, # inches
  dpi = 1000
)

#### Descriptive Statistics by Stratum ========================================

# Summary statistics for each stratum
desc_normal <- model_data_normal %>%
  group_by(transition) %>%
  summarise(
    n = n(),
    n_tb = sum(tb_dx == 1, na.rm = TRUE),
    pct_tb = 100 * mean(tb_dx == 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(stratum = "Normal X-ray")

desc_abnormal <- model_data_abnormal %>%
  group_by(transition) %>%
  summarise(
    n = n(),
    n_tb = sum(tb_dx == 1, na.rm = TRUE),
    pct_tb = 100 * mean(tb_dx == 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(stratum = "Abnormal X-ray")

desc_combined <- bind_rows(desc_normal, desc_abnormal)

print("=== Descriptive Statistics by Stratum ===")
print(desc_combined)

#### Test for Interaction (Optional) ==========================================

# Test if the effect of transition differs by radiographic status
# This is a formal test of whether stratification is necessary

poisson_model_interaction <- glm(
  tb_dx ~ transition * LunitTB_cat + age + prison + prev_tb + smoking +
    any_drugs + prison_time + symptomatic_cellmates + tb_contact,
  data = model_data_stratified,
  family = poisson(link = "log")
)

clustered_se_interaction <- vcovCL(poisson_model_interaction, cluster = ~record_id)

tidy_interaction <- coeftest(poisson_model_interaction, vcov = clustered_se_interaction) %>%
  broom::tidy() %>%
  mutate(
    RR = exp(estimate),
    lower_CI = exp(estimate - 1.96 * std.error),
    upper_CI = exp(estimate + 1.96 * std.error)
  )

# Extract interaction terms
interaction_terms <- tidy_interaction %>%
  filter(str_detect(term, "transition.*LunitTB|LunitTB.*transition"))

print("=== Interaction Terms (Transition × X-ray Status) ===")
print(interaction_terms)

# Likelihood ratio test for interaction
poisson_model_no_interaction <- glm(
  tb_dx ~ transition + LunitTB_cat + age + prison + prev_tb + smoking +
    any_drugs + prison_time + symptomatic_cellmates + tb_contact,
  data = model_data_stratified,
  family = poisson(link = "log")
)

lr_test <- lrtest(poisson_model_no_interaction, poisson_model_interaction)
print("=== Likelihood Ratio Test for Interaction ===")
print(lr_test)

#### Export Results ============================================================

# Save results to CSV for easy reference
write.csv(
  summary_table,
  file = "~/Dropbox/TB-R01-Brazil/2-analysis/Longitudinal Symptoms/stratified_results_xray.csv",
  row.names = FALSE
)

write.csv(
  desc_combined,
  file = "~/Dropbox/TB-R01-Brazil/2-analysis/Longitudinal Symptoms/stratified_descriptives_xray.csv",
  row.names = FALSE
)

print("=== Analysis Complete ===")
print("Results saved to:")
print("  - stratified_results_xray.csv")
print("  - stratified_descriptives_xray.csv")
print("  - longsymptoms_forestplot_stratified_xray.png")
print("  - longsymptoms_forestplot_stratified_xray_sidebyside.png")
