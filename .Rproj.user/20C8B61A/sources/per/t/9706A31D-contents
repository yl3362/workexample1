# ============================================================
# 01 Data preparation and quality checks
# Synthetic veterinary pharmacovigilance work sample
# ============================================================



library(here)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Paths ----------------------------------------------------

input_file <- here(
  "data",
  "Synthetic_Veterinary_PV_Data.xlsx"
)

output_dir <- here(
  "output",
  "01_data_preparation"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(file.exists(input_file))

# 2. Import ---------------------------------------------------

raw <- read_excel(
  input_file,
  sheet = "cases_raw",
  na = c("", "NA")
)

distribution <- read_excel(
  input_file,
  sheet = "monthly_distribution",
  na = c("", "NA")
)

required_columns <- c(
  "record_id", "case_id", "case_version",
  "initial_received_date", "version_received_date",
  "administration_date", "onset_date",
  "product_id", "species", "country",
  "age_years", "weight_kg", "batch_id",
  "sex", "neutered", "primary_event",
  "serious", "seriousness_reason", "outcome",
  "reporter_type", "concomitant_medication",
  "off_label_use"
)

missing_columns <- setdiff(
  required_columns,
  names(raw)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

date_columns <- c(
  "initial_received_date",
  "version_received_date",
  "administration_date",
  "onset_date"
)

raw <- raw %>%
  mutate(
    across(all_of(date_columns), as.Date)
  )

distribution <- distribution %>%
  mutate(month = as.Date(month))

# 3. Check identifiers before selecting versions --------------

stopifnot(
  !anyNA(raw$record_id),
  !anyNA(raw$case_id),
  !anyNA(raw$case_version),
  !anyDuplicated(raw$record_id)
)

duplicate_versions <- raw %>%
  count(case_id, case_version) %>%
  filter(n > 1)

if (nrow(duplicate_versions) > 0) {
  stop("Repeated case_id + case_version combinations found.")
}

# Initial receipt date should be consistent within each case.
inconsistent_receipt <- raw %>%
  group_by(case_id) %>%
  summarise(
    n_dates = n_distinct(initial_received_date),
    .groups = "drop"
  ) %>%
  filter(n_dates != 1)

stopifnot(nrow(inconsistent_receipt) == 0)

# 4. Retain the latest version of each case --------------------
# Follow-up versions are not independent new cases.
# The original 'raw' object remains unchanged.

cases <- raw %>%
  group_by(case_id) %>%
  slice_max(
    order_by = case_version,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  arrange(initial_received_date, case_id) %>%
  mutate(
    reporting_month = as.Date(
      format(initial_received_date, "%Y-%m-01")
    ),
    reporting_year = as.integer(
      format(initial_received_date, "%Y")
    ),
    reporting_half = if_else(
      as.integer(format(initial_received_date, "%m")) <= 6,
      "H1",
      "H2"
    ),
    reporting_period = paste0(
      reporting_year,
      "_",
      reporting_half
    )
  )

stopifnot(!anyDuplicated(cases$case_id))

# 5. Data quality checks --------------------------------------

quality_checks <- tibble(
  check = c(
    "Missing initial receipt date",
    "Missing version receipt date",
    "Initial receipt outside study period",
    "Version received before initial receipt",
    "Version received after snapshot",
    "Administration after event onset",
    "Event onset after initial receipt",
    "Administration after initial receipt",
    "Nonpositive age",
    "Nonpositive weight",
    "Inconsistent seriousness flag",
    "Inconsistent death outcome"
  ),
  n_flagged = c(
    sum(is.na(cases$initial_received_date)),
    sum(is.na(cases$version_received_date)),
    
    sum(
      cases$initial_received_date < as.Date("2023-01-01") |
        cases$initial_received_date > as.Date("2025-12-31"),
      na.rm = TRUE
    ),
    
    sum(
      cases$version_received_date < cases$initial_received_date,
      na.rm = TRUE
    ),
    
    sum(
      cases$version_received_date > as.Date("2025-12-31"),
      na.rm = TRUE
    ),
    
    sum(
      cases$administration_date > cases$onset_date,
      na.rm = TRUE
    ),
    
    sum(
      cases$onset_date > cases$initial_received_date,
      na.rm = TRUE
    ),
    
    sum(
      cases$administration_date > cases$initial_received_date,
      na.rm = TRUE
    ),
    
    sum(cases$age_years <= 0, na.rm = TRUE),
    sum(cases$weight_kg <= 0, na.rm = TRUE),
    
    sum(
      (cases$serious == "No") !=
        (cases$seriousness_reason == "Not applicable"),
      na.rm = TRUE
    ),
    
    sum(
      (cases$outcome == "Died") !=
        (cases$seriousness_reason == "Death"),
      na.rm = TRUE
    )
  )
)

write.csv(
  quality_checks,
  file.path(output_dir, "01_quality_checks.csv"),
  row.names = FALSE
)

print(quality_checks)

if (any(quality_checks$n_flagged > 0)) {
  stop(
    "Data quality checks flagged records. ",
    "Review 01_quality_checks.csv before continuing."
  )
}

# 6. Describe version selection -------------------------------

selection_summary <- tibble(
  stage = c(
    "Raw export rows",
    "Unique cases",
    "Superseded versions excluded",
    "Cases with follow-up history",
    "Final analysis cases"
  ),
  n = c(
    nrow(raw),
    n_distinct(raw$case_id),
    nrow(raw) - nrow(cases),
    sum(table(raw$case_id) > 1),
    nrow(cases)
  )
)

write.csv(
  selection_summary,
  file.path(output_dir, "02_case_selection_summary.csv"),
  row.names = FALSE
)

print(selection_summary)

# 7. Missingness in the final analysis dataset -----------------
# Count both blank cells (NA) and explicit "Unknown" responses.
# Keep their original values unchanged in 'cases'.

is_missing_or_unknown <- function(x) {
  is.na(x) |
    trimws(tolower(as.character(x))) %in% c("", "unknown")
}

missingness <- cases %>%
  select(all_of(required_columns)) %>%
  summarise(
    across(
      everything(),
      ~ sum(is_missing_or_unknown(.x))
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "n_missing"
  ) %>%
  mutate(
    n_cases = nrow(cases),
    percent_missing = 100 * n_missing / n_cases
  ) %>%
  arrange(desc(percent_missing))

write.csv(
  missingness,
  file.path(output_dir, "03_missingness_summary.csv"),
  row.names = FALSE
)

# 8. Figure 1: missing or unknown information ------------------

plot_labels <- c(
  age_years = "Age",
  weight_kg = "Weight",
  administration_date = "Administration date",
  onset_date = "Event onset date",
  batch_id = "Batch",
  sex = "Sex",
  neutered = "Neuter status",
  outcome = "Outcome",
  concomitant_medication = "Concomitant medication",
  off_label_use = "Off-label use"
)

plot_data <- missingness %>%
  filter(variable %in% names(plot_labels)) %>%
  mutate(
    variable_label = unname(plot_labels[variable]),
    variable_label = reorder(
      variable_label,
      percent_missing
    )
  )

p_missing <- ggplot(
  plot_data,
  aes(
    x = percent_missing,
    y = variable_label
  )
) +
  geom_col(
    fill = "#355C7D",
    width = 0.7
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.1f%% (n = %s)",
        percent_missing,
        format(n_missing, big.mark = ",", trim = TRUE)
      )
    ),
    hjust = -0.1,
    size = 3.6
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.4))
  ) +
  labs(
    title = "Completeness of case information",
    subtitle = paste0(
      "Synthetic data; latest version of each case; N = ",
      format(nrow(cases), big.mark = ",", trim = TRUE)
    ),
    x = "Missing or unknown (%)",
    y = NULL,
    caption = paste(
      "Missing includes blank cells and explicit Unknown responses.",
      "Percentages use all unique cases as the denominator."
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0)
  )

print(p_missing)

ggsave(
  filename = file.path(
    output_dir,
    "01_case_information_missingness.png"
  ),
  plot = p_missing,
  width = 9,
  height = 6,
  dpi = 300,
  bg = "white"
)

# 9. Save data for the next analysis step ----------------------

saveRDS(
  list(
    raw = raw,
    cases = cases,
    distribution = distribution,
    quality_checks = quality_checks,
    selection_summary = selection_summary,
    missingness = missingness
  ),
  file.path(output_dir, "01_prepared_data.rds")
)

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "session_info.txt")
)