# ============================================================
# 03 Period comparisons
# Synthetic veterinary pharmacovigilance work sample
# ============================================================

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Paths ----------------------------------------------------

prepared_file <- here(
  "output",
  "01_data_preparation",
  "01_prepared_data.rds"
)

monthly_file <- here(
  "output",
  "02_monthly_trends",
  "02_monthly_trends.rds"
)

output_dir <- here(
  "output",
  "03_period_comparison"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(
  file.exists(prepared_file),
  file.exists(monthly_file)
)

# 2. Load data ------------------------------------------------

prepared <- readRDS(prepared_file)
monthly <- readRDS(monthly_file)

cases <- prepared$cases
monthly_strata <- monthly$monthly_strata

stopifnot(
  !anyDuplicated(cases$case_id),
  !anyNA(cases$initial_received_date)
)

# 3. Define comparison periods -------------------------------

assign_period <- function(x) {
  case_when(
    x >= as.Date("2024-07-01") &
      x < as.Date("2025-01-01") ~ "2024 H2",
    
    x >= as.Date("2025-01-01") &
      x < as.Date("2025-07-01") ~ "2025 H1",
    
    x >= as.Date("2025-07-01") &
      x < as.Date("2026-01-01") ~ "2025 H2",
    
    TRUE ~ NA_character_
  )
}

period_levels <- c(
  "2024 H2",
  "2025 H1",
  "2025 H2"
)

cases_period <- cases %>%
  mutate(
    period = assign_period(initial_received_date)
  ) %>%
  filter(!is.na(period))

monthly_period <- monthly_strata %>%
  mutate(
    period = assign_period(month)
  ) %>%
  filter(!is.na(period))

# 4. Summarize product-level counts and denominators -----------
# Denominators are summed from the monthly distribution grid.
# They are never summed from individual case rows.

period_product <- monthly_period %>%
  group_by(product_id, period) %>%
  summarise(
    n_months = n_distinct(month),
    n_reports = sum(n_reports),
    doses_distributed = sum(doses_distributed),
    .groups = "drop"
  )

stopifnot(
  all(period_product$n_months == 6),
  all(period_product$doses_distributed > 0),
  nrow(period_product) ==
    n_distinct(monthly_strata$product_id) * 3,
  sum(period_product$n_reports) == nrow(cases_period)
)

# 5. Calculate exact Poisson intervals for reporting ratios ----

poisson_lower <- function(n, confidence = 0.95) {
  alpha <- 1 - confidence
  
  result <- numeric(length(n))
  positive <- n > 0
  
  result[positive] <- 0.5 * qchisq(
    alpha / 2,
    df = 2 * n[positive]
  )
  
  result
}

poisson_upper <- function(n, confidence = 0.95) {
  alpha <- 1 - confidence
  
  0.5 * qchisq(
    1 - alpha / 2,
    df = 2 * (n + 1)
  )
}

period_product <- period_product %>%
  mutate(
    reports_per_100k_doses =
      100000 * n_reports / doses_distributed,
    
    reporting_ratio_lower =
      100000 * poisson_lower(n_reports) /
      doses_distributed,
    
    reporting_ratio_upper =
      100000 * poisson_upper(n_reports) /
      doses_distributed
  )

# 6. Construct product-event-period counts --------------------
# Every primary event receives the same product-period
# denominator. Do not sum denominators across event categories.

event_levels <- c(
  "Vomiting",
  "Diarrhea",
  "Lethargy",
  "Pruritus",
  "Decreased appetite",
  "Ataxia",
  "Tremor",
  "Hypersensitivity"
)

stopifnot(
  !anyNA(cases_period$primary_event),
  all(cases_period$primary_event %in% event_levels)
)

event_counts <- cases_period %>%
  count(
    product_id,
    period,
    primary_event,
    name = "n_reports"
  )

event_grid <- expand_grid(
  product_id = sort(unique(period_product$product_id)),
  period = period_levels,
  primary_event = event_levels
)

period_event <- event_grid %>%
  left_join(
    event_counts,
    by = c(
      "product_id",
      "period",
      "primary_event"
    )
  ) %>%
  mutate(
    n_reports = replace_na(n_reports, 0L)
  ) %>%
  left_join(
    period_product %>%
      select(
        product_id,
        period,
        doses_distributed
      ),
    by = c("product_id", "period")
  ) %>%
  mutate(
    reports_per_100k_doses =
      100000 * n_reports / doses_distributed
  )

stopifnot(
  !anyNA(period_event$doses_distributed),
  all(period_event$doses_distributed > 0),
  sum(period_event$n_reports) == nrow(cases_period)
)

# Each case has one primary event in this synthetic dataset.

event_reconciliation <- period_event %>%
  group_by(product_id, period) %>%
  summarise(
    event_total = sum(n_reports),
    .groups = "drop"
  ) %>%
  left_join(
    period_product %>%
      select(product_id, period, n_reports),
    by = c("product_id", "period")
  )

stopifnot(
  all(
    event_reconciliation$event_total ==
      event_reconciliation$n_reports
  )
)

# 7. Define the reporting ratio comparison --------------------
# Ratio = current-period reports per distributed dose divided
# by reference-period reports per distributed dose.
#
# Exact conditional Poisson intervals are used.
# These intervals assume independent Poisson counts and fixed
# denominators. They do not account for reporting bias,
# overdispersion, confounding, or denominator uncertainty.

compare_reporting_ratios <- function(
    current_n,
    reference_n,
    current_doses,
    reference_doses
) {
  
  if (current_n == 0 && reference_n == 0) {
    return(
      tibble(
        ratio = NA_real_,
        ratio_lower = NA_real_,
        ratio_upper = NA_real_,
        interval_status = "Not estimable: both counts are zero"
      )
    )
  }
  
  fit <- poisson.test(
    x = c(current_n, reference_n),
    T = c(current_doses, reference_doses),
    conf.level = 0.95
  )
  
  estimate <- (current_n / current_doses) /
    (reference_n / reference_doses)
  
  tibble(
    ratio = estimate,
    ratio_lower = unname(fit$conf.int[1]),
    ratio_upper = unname(fit$conf.int[2]),
    interval_status = if (
      current_n == 0 || reference_n == 0
    ) {
      "Boundary estimate: one count is zero"
    } else {
      "Finite estimate"
    }
  )
}

make_comparison <- function(
    data,
    reference_period,
    comparison_label
) {
  
  current <- data %>%
    filter(period == "2025 H2") %>%
    select(
      product_id,
      primary_event,
      current_n = n_reports,
      current_doses = doses_distributed
    )
  
  reference <- data %>%
    filter(period == reference_period) %>%
    select(
      product_id,
      primary_event,
      reference_n = n_reports,
      reference_doses = doses_distributed
    )
  
  paired <- current %>%
    left_join(
      reference,
      by = c("product_id", "primary_event")
    )
  
  stopifnot(
    !anyNA(paired),
    all(paired$current_doses > 0),
    all(paired$reference_doses > 0)
  )
  
  estimates <- lapply(
    seq_len(nrow(paired)),
    function(i) {
      compare_reporting_ratios(
        current_n = paired$current_n[i],
        reference_n = paired$reference_n[i],
        current_doses = paired$current_doses[i],
        reference_doses = paired$reference_doses[i]
      )
    }
  ) %>%
    bind_rows()
  
  bind_cols(paired, estimates) %>%
    mutate(
      current_period = "2025 H2",
      reference_period = reference_period,
      comparison = comparison_label,
      
      current_per_100k =
        100000 * current_n / current_doses,
      
      reference_per_100k =
        100000 * reference_n / reference_doses
    )
}

# 8. Compare overall reporting ratios -------------------------

overall_input <- period_product %>%
  mutate(primary_event = "All primary events")

overall_comparisons <- bind_rows(
  make_comparison(
    overall_input,
    reference_period = "2025 H1",
    comparison_label = "2025 H2 vs 2025 H1"
  ),
  make_comparison(
    overall_input,
    reference_period = "2024 H2",
    comparison_label = "2025 H2 vs 2024 H2"
  )
)

print(
  overall_comparisons %>%
    select(
      product_id,
      comparison,
      current_n,
      reference_n,
      ratio,
      ratio_lower,
      ratio_upper
    )
)

# 9. Compare event-specific reporting ratios ------------------

event_comparisons <- bind_rows(
  make_comparison(
    period_event,
    reference_period = "2025 H1",
    comparison_label = "2025 H2 vs 2025 H1"
  ),
  make_comparison(
    period_event,
    reference_period = "2024 H2",
    comparison_label = "2025 H2 vs 2024 H2"
  )
)

# 10. Figure 1: product-level reporting ratios by period -------

period_colors <- c(
  "2024 H2" = "#8B98A5",
  "2025 H1" = "#355C7D",
  "2025 H2" = "#C06C45"
)

period_plot_data <- period_product %>%
  mutate(
    period = factor(
      period,
      levels = period_levels
    )
  )

p_period <- ggplot(
  period_plot_data,
  aes(
    x = period,
    y = reports_per_100k_doses,
    color = period
  )
) +
  geom_errorbar(
    aes(
      ymin = reporting_ratio_lower,
      ymax = reporting_ratio_upper
    ),
    width = 0.15,
    linewidth = 0.7
  ) +
  geom_point(size = 3) +
  facet_wrap(
    ~ product_id,
    nrow = 1
  ) +
  scale_color_manual(values = period_colors) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    title = "Reporting ratios across comparison periods",
    subtitle = "Synthetic data with pointwise 95% Poisson intervals",
    x = NULL,
    y = "Reports per 100,000 distributed doses",
    caption = paste(
      "Each period covers six months. Cases are assigned by initial receipt date.",
      "Ratios are pooled across species and countries without adjustment.",
      "Distributed doses are not a measure of unique animals exposed.",
      sep = "\n"
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(
      hjust = 0,
      size = 9
    ),
    legend.position = "none"
  )

# 11. Figure 2: event-specific period comparisons --------------
# A log scale cannot display zero or infinite estimates.
# Such comparisons are retained in the output table and
# explicitly counted in the figure caption.

event_comparisons <- event_comparisons %>%
  mutate(
    plot_eligible =
      is.finite(ratio) &
      is.finite(ratio_lower) &
      is.finite(ratio_upper) &
      ratio > 0 &
      ratio_lower > 0 &
      ratio_upper > 0
  )

n_not_plotted <- sum(!event_comparisons$plot_eligible)

comparison_levels <- c(
  "2025 H2 vs 2025 H1",
  "2025 H2 vs 2024 H2"
)

forest_data <- event_comparisons %>%
  filter(plot_eligible) %>%
  mutate(
    primary_event = factor(
      primary_event,
      levels = rev(event_levels)
    ),
    comparison = factor(
      comparison,
      levels = comparison_levels
    )
  )

if (nrow(forest_data) == 0) {
  stop("No event comparisons can be displayed on a log scale.")
}

p_event <- ggplot(
  forest_data,
  aes(
    x = primary_event,
    y = ratio,
    color = product_id
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "gray45"
  ) +
  geom_errorbar(
    aes(
      ymin = ratio_lower,
      ymax = ratio_upper
    ),
    width = 0.2,
    linewidth = 0.6
  ) +
  geom_point(size = 2.3) +
  facet_grid(
    comparison ~ product_id,
    drop = FALSE
  ) +
  coord_flip() +
  scale_x_discrete(drop = FALSE) +
  scale_y_log10() +
  scale_color_manual(
    values = c(
      SYN_A = "#355C7D",
      SYN_B = "#C06C45",
      SYN_C = "#54856A"
    )
  ) +
  labs(
    title = "Changes in event-specific reporting ratios",
    subtitle = paste(
      "Current period: 2025 H2;",
      "pointwise 95% exact conditional Poisson intervals"
    ),
    x = NULL,
    y = "Reporting ratio ratio (log scale)",
    caption = paste(
      "Values above 1 indicate more reports per distributed dose in 2025 H2.",
      "Comparisons are exploratory, unadjusted, and do not establish causality.",
      paste0(
        n_not_plotted,
        " boundary or non-estimable comparisons are not plotted; ",
        "all comparisons are retained in the CSV."
      ),
      sep = "\n"
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(
      hjust = 0,
      size = 9
    ),
    legend.position = "none"
  )

print(p_period)
print(p_event)

# 12. Save figures --------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "01_product_reporting_ratios_by_period.png"
  ),
  plot = p_period,
  width = 10,
  height = 5.5,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(
    output_dir,
    "02_event_reporting_ratio_comparisons.png"
  ),
  plot = p_event,
  width = 13,
  height = 9,
  dpi = 300,
  bg = "white"
)

# 13. Save tables and analysis objects -------------------------

tables <- list(
  "period_product" = period_product,
  "period_event" = period_event,
  "overall_period_comparisons" = overall_comparisons,
  "event_period_comparisons" = event_comparisons
)

for (table_name in names(tables)) {
  write.csv(
    tables[[table_name]],
    file.path(
      output_dir,
      paste0(table_name, ".csv")
    ),
    row.names = FALSE
  )
}

saveRDS(
  list(
    period_product = period_product,
    period_event = period_event,
    overall_comparisons = overall_comparisons,
    event_comparisons = event_comparisons
  ),
  file.path(output_dir, "03_period_comparison.rds")
)

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "session_info.txt")
)