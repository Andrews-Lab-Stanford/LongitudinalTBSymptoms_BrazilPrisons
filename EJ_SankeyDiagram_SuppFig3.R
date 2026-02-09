################################################################################
#         SANKEY DIAGRAM: SYMPTOM STATE TRANSITIONS
# Reviewer Request: Visualize symptom transitions for participants with ≥3 visits
# Last Updated: February 2026
# By: Esther Jung
################################################################################

#### Load Libraries ============================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggalluvial) # For Sankey/alluvial diagrams
library(stringr)

#### Descriptive: Participants by Number of Visits ============================

# Count number of symptom assessments per participant
# NOTE: Each row in model_data represents a TRANSITION (interval between 2 visits)
# So n_transitions rows = n_transitions + 1 actual visits (including baseline)
visit_counts <- model_data %>%
    group_by(record_id) %>%
    summarise(
        n_transitions = sum(!is.na(transition)),
        n_visits = n_transitions + 1, # Actual number of visits = transitions + 1
        has_tb = any(tb_dx == "TB", na.rm = TRUE)
    ) %>%
    ungroup()

# Summary statistics
visit_summary <- visit_counts %>%
    summarise(
        total_participants = n(),
        mean_visits = mean(n_visits),
        median_visits = median(n_visits),
        min_visits = min(n_visits),
        max_visits = max(n_visits),
        n_with_3plus = sum(n_visits >= 3),
        pct_with_3plus = 100 * mean(n_visits >= 3),
        n_with_4plus = sum(n_visits >= 4),
        pct_with_4plus = 100 * mean(n_visits >= 4)
    )

print("=== Summary of Participant Visit Counts ===")
print(visit_summary)

# Distribution of visit counts
visit_distribution <- visit_counts %>%
    count(n_visits) %>%
    mutate(
        pct = 100 * n / sum(n),
        cumulative_pct = cumsum(pct)
    )

print("=== Distribution of Number of Visits per Participant ===")
print(visit_distribution)

# Histogram
ggplot(visit_counts, aes(x = n_visits)) +
    geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
    labs(
        title = "Distribution of Number of Symptom Assessments per Participant",
        x = "Number of Symptom Assessments",
        y = "Number of Participants"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold")
    )

ggsave(
    filename = "~/Dropbox/TB-R01-Brazil/figures/longsymptoms_visit_distribution.png",
    width = 6,
    height = 4,
    dpi = 800
)

#### Prepare Data for Sankey Diagram ==========================================

# Filter to participants with at least 3 symptom assessments (2+ transitions)
# Since each row is a transition, we need participants with at least 2 rows
sankey_data <- model_data %>%
    group_by(record_id) %>%
    filter(n() >= 2) %>% # At least 2 transitions = 3 visits total
    arrange(visit_num) %>%
    ungroup()

# Count participants included
n_participants_sankey <- length(unique(sankey_data$record_id))
print(paste("Participants with ≥3 visits (≥2 transitions) included in Sankey:", n_participants_sankey))

# Create flow data for alluvial diagram
# We need to reconstruct the visit sequence including baseline
# Baseline status = symptom_lag at visit_num = 1
# Visit 1 status = symptom_binary at visit_num = 1
# Visit 2 status = symptom_binary at visit_num = 2, etc.

alluvial_prep <- sankey_data %>%
    select(record_id, visit_num, symptom_lag, symptom_binary) %>%
    mutate(
        # Create state labels
        state_current = case_when(
            symptom_binary == 1 ~ "Symptomatic",
            symptom_binary == 0 ~ "Asymptomatic",
            TRUE ~ NA_character_
        ),
        state_previous = case_when(
            symptom_lag == 1 ~ "Symptomatic",
            symptom_lag == 0 ~ "Asymptomatic",
            TRUE ~ NA_character_
        )
    )

# For each participant, get baseline and first 2 follow-up visits
alluvial_data <- alluvial_prep %>%
    filter(visit_num <= 2) %>%
    group_by(record_id) %>%
    arrange(visit_num) %>%
    summarise(
        # Baseline is symptom_lag from first transition
        Baseline = first(state_previous),
        # Visit 1 is symptom_binary from first transition
        Visit_1 = first(state_current),
        # Visit 2 is symptom_binary from second transition (if exists)
        Visit_2 = if (n() >= 2) state_current[2] else NA_character_,
        .groups = "drop"
    ) %>%
    filter(!is.na(Baseline) & !is.na(Visit_1) & !is.na(Visit_2)) %>% # Require baseline + 2 follow-ups
    count(Baseline, Visit_1, Visit_2, name = "Freq") %>%
    arrange(desc(Freq))

# Calculate final N for participants with complete Baseline + 2 follow-ups
n_final <- sum(alluvial_data$Freq)
print(paste("Total participants in Sankey (Baseline + 2 follow-ups):", n_final))

# Convert to lodes (long) format for more control over labels
alluvial_lodes <- alluvial_data %>%
    mutate(baseline_status = Baseline) %>% # Preserve baseline status as an ID for coloring
    to_lodes_form(key = "Visit", axes = 1:3) %>%
    group_by(Visit, stratum) %>%
    mutate(
        # Calculate sub-n for each stratum at each visit
        stratum_n = sum(Freq),
        # Create final label with whole-number n
        label_with_n = paste0(stratum, "\n(n = ", stratum_n, ")")
    ) %>%
    ungroup()

# Create Sankey diagram using ggalluvial
sankey_plot <- ggplot(
    alluvial_lodes,
    aes(x = Visit, stratum = stratum, alluvium = alluvium, y = Freq)
) +
    geom_alluvium(aes(fill = baseline_status), width = 1 / 12, alpha = 0.7, curve_type = "sigmoid") +
    geom_stratum(width = 1 / 12, fill = "gray90", color = "gray50") +
    # Use stat = "stratum" to correctly position and aggregate labels per box
    geom_label(
        stat = "stratum",
        aes(label = label_with_n),
        size = 3
    ) +
    scale_x_discrete(
        limits = c("Baseline", "Visit_1", "Visit_2"),
        labels = c("Baseline", "Visit 1", "Visit 2"),
        expand = c(0.05, 0.05)
    ) +
    scale_fill_manual(
        values = c(
            "Asymptomatic" = "#4575B4",
            "Symptomatic" = "#D73027"
        ),
        name = "Baseline Status"
    ) +
    labs(
        title = NULL,
        subtitle = NULL,
        x = NULL,
        y = "Number of Participants"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        legend.position = "bottom",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 11, face = "bold")
    )

print(sankey_plot)

ggsave(
    filename = "~/Dropbox/TB-R01-Brazil/figures/longsymptoms_sankey_diagram.png",
    plot = sankey_plot,
    width = 10,
    height = 6,
    dpi = 800
)

#### Create Summary Table Instead of Sankey ===================================

# The Sankey diagram is too messy with this many transitions
# Instead, create a clean summary table showing visit patterns

# Summary by baseline status and number of visits
visit_pattern_summary <- model_data %>%
    group_by(record_id) %>%
    summarise(
        n_transitions = n(),
        n_visits = n_transitions + 1,
        baseline_status = first(ifelse(symptom_lag == 1, "Symptomatic", "Asymptomatic")),
        has_tb = any(tb_dx == "TB", na.rm = TRUE)
    ) %>%
    ungroup() %>%
    group_by(baseline_status, n_visits) %>%
    summarise(
        n_participants = n(),
        n_tb = sum(has_tb),
        pct_tb = 100 * mean(has_tb),
        .groups = "drop"
    ) %>%
    arrange(baseline_status, n_visits)

print("=== Visit Pattern Summary by Baseline Status ===")
print(visit_pattern_summary)

# Create a cleaner visualization: bar chart of visit counts by baseline status
visit_bar_plot <- model_data %>%
    group_by(record_id) %>%
    summarise(
        n_visits = n() + 1,
        baseline_status = first(ifelse(symptom_lag == 1, "Symptomatic", "Asymptomatic"))
    ) %>%
    ggplot(aes(x = n_visits, fill = baseline_status)) +
    geom_bar(position = "dodge", alpha = 0.8) +
    scale_fill_manual(
        values = c("Asymptomatic" = "#4575B4", "Symptomatic" = "#D73027"),
        name = "Baseline Status"
    ) +
    labs(
        title = "Distribution of Follow-up Visits by Baseline Symptom Status",
        x = "Number of Symptom Assessments",
        y = "Number of Participants"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom"
    )

print(visit_bar_plot)

ggsave(
    filename = "~/Dropbox/TB-R01-Brazil/figures/longsymptoms_visit_distribution_by_baseline.png",
    plot = visit_bar_plot,
    width = 8,
    height = 5,
    dpi = 800
)

#### Export Results ============================================================

# Save summary statistics
write.csv(
    visit_summary,
    file = "~/Dropbox/TB-R01-Brazil/2-analysis/Longitudinal Symptoms/visit_summary_statistics.csv",
    row.names = FALSE
)

write.csv(
    visit_distribution,
    file = "~/Dropbox/TB-R01-Brazil/2-analysis/Longitudinal Symptoms/visit_distribution.csv",
    row.names = FALSE
)

write.csv(
    visit_pattern_summary,
    file = "~/Dropbox/TB-R01-Brazil/2-analysis/Longitudinal Symptoms/visit_pattern_by_baseline.csv",
    row.names = FALSE
)

print("=== Analysis Complete ===")
print("Generated files:")
print("  - longsymptoms_visit_distribution.png")
print("  - longsymptoms_visit_distribution_by_baseline.png")
print("  - visit_summary_statistics.csv")
print("  - visit_distribution.csv")
print("  - visit_pattern_by_baseline.csv")
