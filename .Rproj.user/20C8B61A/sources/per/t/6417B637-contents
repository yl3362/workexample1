# ============================================================
# 02 Monthly reporting trends and distribution volume
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
  "02_monthly_trends"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(file.exists(input_file))

# 2. Load prepared data ---------------------------------------

prepared <- readRDS(input_file)

cases <- prepared$cases
distribution <- prepared$distribution

stopifnot(
  !anyDuplicated(cases$case_id),
  !anyNA(cases$initial_received_date)
)

cases <- cases %>%
  mutate(
    reporting_month = as.Date(
      format(initial_received_date, "%Y-%m-01")
    )
  )

distribution <- distribution %>%
  mutate(month = as.Date(month))

join_keys <- c(
  "month",
  "product_id",
  "species",
  "country"
)

# 3. Validate the distribution table --------------------------

stopifnot(
  all(join_keys %in% names(distribution)),
  "doses_distributed" %in% names(distribution),
  !anyNA(distribution[join_keys]),
  is.numeric(distribution$doses_distributed),
  all(is.finite(distribution$doses_distributed)),
  all(distribution$doses_distributed > 0),
  all(
    distribution$doses_distributed ==
      floor(distribution$doses_distributed)
  ),
  all(format(distribution$month, "%d") == "01")
)

duplicate_keys <- distribution %>%
  count(
    across(all_of(join_keys)),
    name = "n_rows"
  ) %>%
  filter(n_rows > 1)

if (nrow(duplicate_keys) > 0) {
  stop("Duplicate monthly distribution keys were found.")
}

expected_grid <- expand_grid(
  month = seq(
    as.Date("2023-01-01"),
    as.Date("2025-12-01"),
    by = "month"
  ),
  product_id = c("SYN_A", "SYN_B", "SYN_C"),
  species = c("Dog", "Cat"),
  country = c("US", "GB", "DE")
)

missing_distribution <- expected_grid %>%
  anti_join(distribution, by = join_keys)

unexpected_distribution <- distribution %>%
  anti_join(expected_grid, by = join_keys)

if (
  nrow(missing_distribution) > 0 ||
  nrow(unexpected_distribution) > 0
) {
  stop("The distribution table does not match the expected grid.")
}

# 4. Aggregate cases before joining denominators --------------

case_counts <- cases %>%
  count(
    month = reporting_month,
    product_id,
    species,
    country,
    name = "n_reports"
  )

unmatched_cases <- case_counts %>%
  anti_join(
    distribution,
    by = join_keys
  )

if (nrow(unmatched_cases) > 0) {
  stop("Some case strata have no matching distribution records.")
}

# The complete distribution grid retains months with zero reports.
# A zero is assigned only when the matching case count is absent.

monthly_strata <- distribution %>%
  left_join(
    case_counts,
    by = join_keys
  ) %>%
  mutate(
    n_reports = replace_na(n_reports, 0L),
    reports_per_100k_doses =
      100000 * n_reports / doses_distributed
  ) %>%
  arrange(
    product_id,
    species,
    country,
    month
  )

stopifnot(
  nrow(monthly_strata) == nrow(distribution),
  sum(monthly_strata$n_reports) == nrow(cases),
  sum(monthly_strata$doses_distributed) ==
    sum(distribution$doses_distributed)
)

# 5. Calculate product-level monthly summaries ----------------
# Sum counts and denominators before calculating each ratio.
# Do not average the stratum-specific ratios.

monthly_product <- monthly_strata %>%
  group_by(month, product_id) %>%
  summarise(
    n_reports = sum(n_reports),
    doses_distributed = sum(doses_distributed),
    .groups = "drop"
  ) %>%
  mutate(
    reports_per_100k_doses =
      100000 * n_reports / doses_distributed
  ) %>%
  arrange(product_id, month)

# 6. Add exact Poisson intervals ------------------------------
# These intervals describe count uncertainty conditional on
# the supplied denominator and a Poisson count model.
# They do not account for underreporting, distribution-to-use
# differences, confounding, or uncertainty in the denominator.

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

monthly_product <- monthly_product %>%
  mutate(
    reporting_ratio_lower =
      100000 * poisson_lower(n_reports) /
      doses_distributed,
    reporting_ratio_upper =
      100000 * poisson_upper(n_reports) /
      doses_distributed
  )

stopifnot(
  all(
    monthly_product$reporting_ratio_lower <=
      monthly_product$reports_per_100k_doses
  ),
  all(
    monthly_product$reporting_ratio_upper >=
      monthly_product$reports_per_100k_doses
  )
)

# 7. Overall descriptive summary ------------------------------

product_summary <- monthly_product %>%
  group_by(product_id) %>%
  summarise(
    n_months = n(),
    total_reports = sum(n_reports),
    total_doses_distributed = sum(doses_distributed),
    .groups = "drop"
  ) %>%
  mutate(
    reports_per_100k_doses =
      100000 * total_reports / total_doses_distributed
  )

print(product_summary)

# 8. Plot settings --------------------------------------------

product_colors <- c(
  SYN_A = "#355C7D",
  SYN_B = "#C06C45",
  SYN_C = "#54856A"
)

trend_theme <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.text = element_text(
        face = "bold",
        hjust = 0
      ),
      plot.title = element_text(face = "bold"),
      plot.caption = element_text(
        hjust = 0,
        size = 9
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      legend.position = "none"
    )
}

month_axis <- function() {
  scale_x_date(
    date_breaks = "6 months",
    date_labels = "%b\n%Y"
  )
}

count_labels <- function(x) {
  format(
    x,
    big.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}

# 9. Figure 1: monthly report counts ---------------------------

p_counts <- ggplot(
  monthly_product,
  aes(
    x = month,
    y = n_reports,
    color = product_id
  )
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.5) +
  facet_wrap(
    ~ product_id,
    ncol = 1,
    scales = "fixed"
  ) +
  scale_color_manual(values = product_colors) +
  month_axis() +
  scale_y_continuous(
    limits = c(0, NA),
    labels = count_labels,
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    title = "Monthly adverse event report counts",
    subtitle = "Synthetic data, January 2023 to December 2025",
    x = NULL,
    y = "Unique reports",
    caption = paste(
      "Each case is counted once using its initial receipt month.",
      "The latest available case version is retained."
    )
  ) +
  trend_theme()

# 10. Figure 2: monthly distribution volume --------------------

p_distribution <- ggplot(
  monthly_product,
  aes(
    x = month,
    y = doses_distributed,
    color = product_id
  )
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.5) +
  facet_wrap(
    ~ product_id,
    ncol = 1,
    scales = "fixed"
  ) +
  scale_color_manual(values = product_colors) +
  month_axis() +
  scale_y_continuous(
    limits = c(0, NA),
    labels = count_labels,
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    title = "Monthly product distribution",
    subtitle = "Synthetic data, January 2023 to December 2025",
    x = NULL,
    y = "Doses distributed",
    caption = paste(
      "Distribution is summed across species and countries.",
      "Distributed doses do not represent unique animals treated."
    )
  ) +
  trend_theme()

# 11. Figure 3: monthly reporting ratios -----------------------

p_ratio <- ggplot(
  monthly_product,
  aes(
    x = month,
    y = reports_per_100k_doses,
    color = product_id,
    fill = product_id
  )
) +
  geom_ribbon(
    aes(
      ymin = reporting_ratio_lower,
      ymax = reporting_ratio_upper
    ),
    alpha = 0.15,
    color = NA
  ) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.5) +
  facet_wrap(
    ~ product_id,
    ncol = 1,
    scales = "fixed"
  ) +
  scale_color_manual(values = product_colors) +
  scale_fill_manual(values = product_colors) +
  month_axis() +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    title = "Monthly reports relative to distribution volume",
    subtitle = "Synthetic data with pointwise 95% Poisson intervals",
    x = NULL,
    y = "Reports per 100,000 distributed doses",
    caption = paste(
      "Descriptive reporting ratios, not incidence or causal risk.",
      "Intervals assume Poisson counts and fixed denominators.",
      "Results are pooled across species and countries without adjustment.",
      sep = "\n"
    )
  ) +
  trend_theme()

print(p_counts)
print(p_distribution)
print(p_ratio)

# 12. Save figures --------------------------------------------

plots <- list(
  "01_monthly_report_counts" = p_counts,
  "02_monthly_distribution" = p_distribution,
  "03_monthly_reporting_ratios" = p_ratio
)

for (plot_name in names(plots)) {
  ggsave(
    filename = file.path(
      output_dir,
      paste0(plot_name, ".png")
    ),
    plot = plots[[plot_name]],
    width = 10,
    height = 8,
    dpi = 300,
    bg = "white"
  )
}

# 13. Save tables and analysis objects -------------------------

tables <- list(
  "monthly_strata" = monthly_strata,
  "monthly_product" = monthly_product,
  "product_summary" = product_summary
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
    monthly_strata = monthly_strata,
    monthly_product = monthly_product,
    product_summary = product_summary
  ),
  file.path(output_dir, "02_monthly_trends.rds")
)

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "session_info.txt")
)