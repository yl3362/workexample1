# ============================================================
# 04 Product-event disproportionality analysis
# Synthetic veterinary pharmacovigilance work sample
# ============================================================

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Paths ----------------------------------------------------

input_file <- here(
  "output",
  "01_data_preparation",
  "01_prepared_data.rds"
)

output_dir <- here(
  "output",
  "04_disproportionality"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(file.exists(input_file))

# 2. Load prepared cases --------------------------------------

prepared <- readRDS(input_file)
cases <- prepared$cases

product_levels <- c(
  "SYN_A",
  "SYN_B",
  "SYN_C"
)

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

required_columns <- c(
  "case_id",
  "product_id",
  "primary_event",
  "initial_received_date",
  "species",
  "country",
  "serious",
  "outcome",
  "batch_id"
)

stopifnot(
  all(required_columns %in% names(cases)),
  !anyNA(cases$case_id),
  !anyDuplicated(cases$case_id),
  !anyNA(cases$product_id),
  !anyNA(cases$primary_event),
  !anyNA(cases$initial_received_date),
  all(cases$product_id %in% product_levels),
  all(cases$primary_event %in% event_levels)
)

cases <- cases %>%
  mutate(
    initial_received_date = as.Date(initial_received_date)
  )

stopifnot(
  all(
    cases$initial_received_date >= as.Date("2023-01-01") &
      cases$initial_received_date <= as.Date("2025-12-31")
  )
)

# 3. Define analysis windows ----------------------------------
# The recent window and the full window overlap.
# Their results are not independent replications.

window_definitions <- tibble(
  analysis_window = c(
    "2025 H2",
    "2023-2025"
  ),
  start_date = as.Date(c(
    "2025-07-01",
    "2023-01-01"
  )),
  end_date = as.Date(c(
    "2026-01-01",
    "2026-01-01"
  ))
)

# 4. Construct product-event 2 x 2 tables -----------------------
#
#                         Target event    Other events
# Target product               a               b
# Other two products           c               d
#
# ROR = (a / b) / (c / d) = a * d / (b * c)
#
# Each case contributes exactly one primary event in this
# synthetic dataset. No distribution denominator is used.

make_contingency_counts <- function(
    data,
    window_name,
    start_date,
    end_date
) {
  
  selected <- data %>%
    filter(
      initial_received_date >= start_date,
      initial_received_date < end_date
    )
  
  total_reports <- nrow(selected)
  
  if (total_reports == 0) {
    stop("No reports found for analysis window: ", window_name)
  }
  
  product_counts <- selected %>%
    count(
      product_id,
      name = "product_total"
    )
  
  event_counts <- selected %>%
    count(
      primary_event,
      name = "event_total"
    )
  
  pair_counts <- selected %>%
    count(
      product_id,
      primary_event,
      name = "a"
    )
  
  result <- expand_grid(
    product_id = product_levels,
    primary_event = event_levels
  ) %>%
    left_join(
      pair_counts,
      by = c("product_id", "primary_event")
    ) %>%
    left_join(
      product_counts,
      by = "product_id"
    ) %>%
    left_join(
      event_counts,
      by = "primary_event"
    ) %>%
    mutate(
      across(
        c(a, product_total, event_total),
        ~ replace_na(.x, 0L)
      ),
      b = product_total - a,
      c = event_total - a,
      d = total_reports - a - b - c,
      analysis_window = window_name,
      total_reports = total_reports,
      comparator_total = c + d
    )
  
  stopifnot(
    all(result$a >= 0),
    all(result$b >= 0),
    all(result$c >= 0),
    all(result$d >= 0),
    all(
      result$a + result$b + result$c + result$d ==
        total_reports
    )
  )
  
  result
}

contingency_counts <- lapply(
  seq_len(nrow(window_definitions)),
  function(i) {
    
    make_contingency_counts(
      data = cases,
      window_name = window_definitions$analysis_window[i],
      start_date = window_definitions$start_date[i],
      end_date = window_definitions$end_date[i]
    )
  }
) %>%
  bind_rows()

# 5. Estimate ROR and uncertainty ------------------------------
# For positive cells, use the uncorrected cross-product ROR
# with a log-Wald confidence interval.
#
# If any cell is zero, add 0.5 to all four cells and label
# the result as continuity corrected.
#
# Tables with an empty row or column margin are not estimable.
# Fisher's exact test provides a separate exploratory p-value.

estimate_ror <- function(a, b, c, d) {
  
  margins <- c(
    a + b,
    c + d,
    a + c,
    b + d
  )
  
  if (any(margins == 0)) {
    return(
      tibble(
        ror = NA_real_,
        ror_lower = NA_real_,
        ror_upper = NA_real_,
        p_value = NA_real_,
        zero_cell_correction = FALSE,
        estimation_method = "Not estimable: empty margin"
      )
    )
  }
  
  cells <- as.numeric(c(a, b, c, d))
  corrected <- any(cells == 0)
  
  adjusted <- if (corrected) {
    cells + 0.5
  } else {
    cells
  }
  
  log_ror <- log(adjusted[1]) +
    log(adjusted[4]) -
    log(adjusted[2]) -
    log(adjusted[3])
  
  standard_error <- sqrt(
    sum(1 / adjusted)
  )
  
  z <- qnorm(0.975)
  
  contingency_table <- matrix(
    cells,
    nrow = 2,
    byrow = TRUE
  )
  
  fisher_result <- fisher.test(
    contingency_table,
    alternative = "two.sided"
  )
  
  tibble(
    ror = exp(log_ror),
    ror_lower = exp(log_ror - z * standard_error),
    ror_upper = exp(log_ror + z * standard_error),
    p_value = fisher_result$p.value,
    zero_cell_correction = corrected,
    estimation_method = if (corrected) {
      "0.5-corrected ROR with log-Wald interval"
    } else {
      "Uncorrected ROR with log-Wald interval"
    }
  )
}

estimates <- lapply(
  seq_len(nrow(contingency_counts)),
  function(i) {
    
    estimate_ror(
      a = contingency_counts$a[i],
      b = contingency_counts$b[i],
      c = contingency_counts$c[i],
      d = contingency_counts$d[i]
    )
  }
) %>%
  bind_rows()

ror_results <- bind_cols(
  contingency_counts,
  estimates
) %>%
  mutate(
    target_event_percent = if_else(
      a + b > 0,
      100 * a / (a + b),
      NA_real_
    ),
    comparator_event_percent = if_else(
      c + d > 0,
      100 * c / (c + d),
      NA_real_
    )
  )

# Adjust across the 24 product-event tests within each window.
# Windows overlap, and tests share reports and comparators.
# Adjusted p-values remain exploratory rather than a
# validated safety-signal decision rule.

ror_results <- ror_results %>%
  group_by(analysis_window) %>%
  mutate(
    p_adjusted_bh = p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%
  ungroup()

# 6. Create an exploratory review flag ------------------------
# This is a transparent portfolio exercise rule.
# It is not an official regulatory or company threshold.
#
# Require:
# - at least 5 target-product reports of the target event;
# - no zero-cell correction;
# - lower 95% ROR bound above 1.
#
# The flag is for organizing case review, not declaring
# causality, statistical confirmation, or clinical risk.

minimum_pair_reports <- 5L

ror_results <- ror_results %>%
  mutate(
    minimum_pair_reports = minimum_pair_reports,
    
    exploratory_review_flag =
      a >= minimum_pair_reports &
      !zero_cell_correction &
      !is.na(ror_lower) &
      ror_lower > 1,
    
    sparse_cell_flag = pmin(a, b, c, d) < 5,
    
    review_reason = case_when(
      is.na(ror) ~
        "Not estimable",
      
      exploratory_review_flag ~
        "Meets the prespecified exploratory review rule",
      
      zero_cell_correction ~
        "Zero cell: continuity-corrected estimate",
      
      a < minimum_pair_reports ~
        "Fewer than 5 target product-event reports",
      
      TRUE ~
        "Does not meet the exploratory review rule"
    )
  ) %>%
  arrange(
    analysis_window,
    product_id,
    desc(ror)
  )

print(
  ror_results %>%
    filter(analysis_window == "2025 H2") %>%
    select(
      product_id,
      primary_event,
      a,
      b,
      c,
      d,
      ror,
      ror_lower,
      ror_upper,
      exploratory_review_flag
    )
)

# 7. Figure: product-event ROR estimates -----------------------
# All estimable pairs are shown, not only flagged pairs.
# A dagger marks estimates requiring a zero-cell correction.

window_levels <- c(
  "2025 H2",
  "2023-2025"
)

ror_results <- ror_results %>%
  mutate(
    plot_eligible =
      is.finite(ror) &
      is.finite(ror_lower) &
      is.finite(ror_upper) &
      ror > 0 &
      ror_lower > 0 &
      ror_upper > 0
  )

n_not_plotted <- sum(!ror_results$plot_eligible)

plot_data <- ror_results %>%
  filter(plot_eligible) %>%
  mutate(
    primary_event = factor(
      primary_event,
      levels = rev(event_levels)
    ),
    analysis_window = factor(
      analysis_window,
      levels = window_levels
    ),
    correction_marker = if_else(
      zero_cell_correction,
      "\u2020",
      ""
    )
  )

if (nrow(plot_data) == 0) {
  stop("No ROR estimates can be displayed on a log scale.")
}

p_ror <- ggplot(
  plot_data,
  aes(
    x = primary_event,
    y = ror,
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
      ymin = ror_lower,
      ymax = ror_upper
    ),
    width = 0.2,
    linewidth = 0.6
  ) +
  geom_point(size = 2.3) +
  geom_text(
    aes(label = correction_marker),
    nudge_x = 0.25,
    size = 3,
    show.legend = FALSE
  ) +
  facet_grid(
    analysis_window ~ product_id,
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
    title = "Product-event reporting disproportionality",
    subtitle = paste(
      "Synthetic data;",
      "reporting odds ratios with pointwise 95% intervals"
    ),
    x = NULL,
    y = "Reporting odds ratio (log scale)",
    caption = paste(
      "Each product is compared with the other two products in the same window.",
      "Results are unadjusted; overlapping windows are not independent replications.",
      "\u2020 indicates a 0.5 continuity correction. Intervals are not multiplicity adjusted.",
      paste0(
        n_not_plotted,
        " non-estimable or non-finite results are omitted from the plot ",
        "and retained in the CSV."
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

print(p_ror)

ggsave(
  filename = file.path(
    output_dir,
    "01_product_event_ror.png"
  ),
  plot = p_ror,
  width = 13,
  height = 9,
  dpi = 300,
  bg = "white"
)

# 8. Prepare a review table for the recent window --------------

recent_review_table <- ror_results %>%
  filter(
    analysis_window == "2025 H2",
    exploratory_review_flag
  ) %>%
  arrange(
    desc(ror_lower),
    desc(a)
  )

# Retain all matching cases, without selecting by outcome.
# Empty output is valid if no pairs meet the review rule.

review_pairs <- recent_review_table %>%
  select(
    product_id,
    primary_event
  ) %>%
  distinct()

recent_review_cases <- cases %>%
  filter(
    initial_received_date >= as.Date("2025-07-01"),
    initial_received_date < as.Date("2026-01-01")
  ) %>%
  semi_join(
    review_pairs,
    by = c("product_id", "primary_event")
  ) %>%
  arrange(
    product_id,
    primary_event,
    initial_received_date,
    case_id
  )

stopifnot(
  !anyDuplicated(recent_review_cases$case_id)
)

# 9. Summarize the exploratory review workload -----------------

review_summary <- ror_results %>%
  group_by(analysis_window) %>%
  summarise(
    n_product_event_pairs = n(),
    n_estimable_pairs = sum(!is.na(ror)),
    n_corrected_pairs = sum(zero_cell_correction),
    n_sparse_pairs = sum(sparse_cell_flag),
    n_pairs_for_exploratory_review =
      sum(exploratory_review_flag),
    .groups = "drop"
  )

print(review_summary)

# 10. Save tables and analysis objects -------------------------

tables <- list(
  "product_event_ror_results" = ror_results,
  "recent_exploratory_review_table" = recent_review_table,
  "recent_exploratory_review_cases" = recent_review_cases,
  "review_summary" = review_summary
)

for (table_name in names(tables)) {
  write.csv(
    tables[[table_name]],
    file.path(
      output_dir,
      paste0(table_name, ".csv")
    ),
    row.names = FALSE,
    na = ""
  )
}

saveRDS(
  list(
    ror_results = ror_results,
    recent_review_table = recent_review_table,
    recent_review_cases = recent_review_cases,
    review_summary = review_summary,
    window_definitions = window_definitions,
    minimum_pair_reports = minimum_pair_reports
  ),
  file.path(output_dir, "04_disproportionality.rds")
)

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "session_info.txt")
)