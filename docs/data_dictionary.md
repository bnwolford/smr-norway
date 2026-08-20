# Data dictionary

## Input: cohort file (one row per participant)

Column names are configurable via the `COL_*` settings in each analysis script;
the names below are the defaults used in `smr_analysis_dates.R` and the shipped
example.

| Column | Type | Required | Description |
|--------|------|----------|-------------|
| `id` | any | yes | Unique participant identifier. |
| `sex` | string/int | yes | Sex, coded however you like; mapped to Male/Female via `SEX_MAP`. |
| `start_obs` | date | yes | Date observation begins (cohort entry). |
| `end_obs` | date | yes | Date observation ends (administrative censoring date). |
| `death_date` | date | no | Date of death; blank/NA for survivors. A death after `end_obs` is treated as censoring. |
| `age_entry` | numeric | yes* | Age at `start_obs`. *Required only for the age-input and three-date scripts (no DOB); omit if you supply a true date of birth.* |
| `birth_date` | date | yes* | Date of birth. *Required only for `smr_analysis.R`; set `COL_BIRTH` and leave `age_entry` unused.* |

Dates may be in any format listed in the `DATE_FORMATS` setting (ISO
`YYYY-MM-DD` recommended).

## Reference rates: `data/norway_reference_mortality_rates.csv`

Source: Statistics Norway StatBank table 07902 (see `data-raw/fetch_ssb_rates.R`).

| Column | Type | Description |
|--------|------|-------------|
| `sex` | string | "Male" or "Female". |
| `age` | integer | Single year of age, 0–106. |
| `year` | integer | Calendar year. |
| `qx` | numeric | Probability of death between age *x* and *x*+1 (0–1). |
| `mx` | numeric | Central mortality rate = `qx / (1 − qx/2)`. **This is multiplied by person-years to get expected deaths.** |
| `qx_per_1000` | numeric | Raw SSB value (probability of death per 1,000), retained for traceability. |

## Output tables (written by every analysis script)

`smr_overall.csv`, `smr_by_sex.csv`, `smr_by_age.csv` share these columns:

| Column | Description |
|--------|-------------|
| `sex` / `age_band` | Stratifying variables (absent in the overall table). |
| `observed` | Observed deaths in the cohort (integer). |
| `expected` | Expected deaths = Σ(person-years × `mx`). |
| `person_years` | Person-time contributed by the stratum. |
| `smr` | Standardized mortality ratio = observed / expected. |
| `smr_lo`, `smr_hi` | Exact Poisson (Garwood) 95% confidence limits for the SMR. |
| `p_value` | Two-sided Poisson p-value testing H₀: SMR = 1. |

`smr_person_time_cells.csv` — the full Lexis table, one row per
age × sex × calendar-year cell: `sex`, `age`, `year`, `py` (person-years),
`ev` (deaths in the cell), `mx`, `expected`, and `age_band`.

`smr_excluded_rows.csv` — any input rows dropped for missing/invalid values or
`exit < entry` (written only if such rows exist).
