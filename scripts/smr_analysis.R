################################################################################
## Standardized Mortality Ratio (SMR) vs. the Norwegian general population
## ---------------------------------------------------------------------------
## Compares OBSERVED deaths in your cohort to EXPECTED deaths, where expected
## deaths come from applying age-, sex- and calendar-year-specific mortality
## rates (Statistics Norway life tables) to the person-time your cohort
## actually contributed.
##
## Method: exact Lexis expansion (person-time split across age-band x calendar-
## year cells as each subject ages), then  SMR = O / E  with Poisson-based CIs.
##
## Dependencies: data.table only (base R otherwise).
################################################################################

library(data.table)

## ============================================================================
## 1. USER SETTINGS  --  edit everything in this block to fit your data
## ============================================================================

## ---- 1a. File paths --------------------------------------------------------
COHORT_FILE <- "cohort.csv"                          # your individual-level data
RATES_FILE  <- "norway_reference_mortality_rates.csv" # reference rates (from SSB)
OUT_DIR     <- "."                                    # where results are written

## ---- 1b. Column names in YOUR cohort file ----------------------------------
##   Rename the right-hand strings to match your file's headers.
COL_ID       <- "id"            # unique subject identifier
COL_SEX      <- "sex"           # sex variable
COL_BIRTH    <- "dob"           # date of birth
COL_ENTRY    <- "entry_date"    # start of follow-up (cohort entry)
COL_EXIT     <- "exit_date"     # end of follow-up (death / censor / study end)
COL_STATUS   <- "status"        # vital status at exit (see SEX/STATUS coding below)

## ---- 1c. How sex is coded in YOUR file -------------------------------------
##   Map YOUR codes -> "Male"/"Female" (which is how the rates file is coded).
##   Examples: list("Male"=c("M","1","male"), "Female"=c("F","2","female"))
SEX_MAP <- list(
  Male   = c("Male",   "M", "1", "male",   "m"),
  Female = c("Female", "F", "2", "female", "f")
)

## ---- 1d. How the death event is coded in COL_STATUS ------------------------
##   Which value(s) in COL_STATUS mean "died". Everything else = censored/alive.
DEATH_CODES <- c("1", "dead", "death", "died", "yes", "TRUE")

## ---- 1e. Date format -------------------------------------------------------
##   Format string(s) used by as.Date(). If dates are already Date objects or
##   ISO "YYYY-MM-DD", the default is fine. Add more formats if mixed.
DATE_FORMATS <- c("%Y-%m-%d", "%d.%m.%Y", "%d/%m/%Y", "%m/%d/%Y")

## ---- 1f. Age banding for reporting -----------------------------------------
##   The reference rates are single-year-of-age, so person-time is always split
##   exactly. AGE_BREAKS only controls how results are GROUPED for the strata
##   table. Standard 5-year bands with an open-ended top group:
AGE_BREAKS <- c(seq(0, 85, by = 5), Inf)   # 0-4,5-9,...,80-84,85+

## ---- 1g. Optional: restrict the analysis -----------------------------------
##   Clamp follow-up to the calendar years covered by your rates file, or to a
##   study window. Set to NULL to use the full span present in the data.
STUDY_START <- NULL   # e.g. as.Date("2005-01-01")
STUDY_END   <- NULL   # e.g. as.Date("2020-12-31")

## ============================================================================
## 2. HELPERS  (no editing needed below unless you want to change the method)
## ============================================================================

## Robust multi-format date parser
parse_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  x <- as.character(x)
  out <- as.Date(rep(NA_real_, length(x)), origin = "1970-01-01")
  for (fmt in DATE_FORMATS) {
    miss <- is.na(out)
    if (!any(miss)) break
    out[miss] <- as.Date(x[miss], format = fmt)
  }
  out
}

## Map a raw sex vector to "Male"/"Female"
recode_sex <- function(x) {
  x <- as.character(x)
  out <- rep(NA_character_, length(x))
  for (target in names(SEX_MAP)) out[x %in% SEX_MAP[[target]]] <- target
  out
}

## Exact Poisson 95% CI for an observed count O (Garwood method)
pois_ci <- function(O, conf = 0.95) {
  a <- (1 - conf) / 2
  lo <- ifelse(O == 0, 0, qgamma(a,     O,       rate = 1))
  hi <-                    qgamma(1 - a, O + 1,   rate = 1)
  c(lower = lo, upper = hi)
}

## ============================================================================
## 3. LOAD DATA
## ============================================================================

rates <- fread(RATES_FILE)
## expected columns: sex, age, year, mx   (mx = central mortality rate / person-year)
stopifnot(all(c("sex", "age", "year", "mx") %in% names(rates)))
setkey(rates, sex, age, year)

raw <- fread(COHORT_FILE)

coh <- data.table(
  id     = raw[[COL_ID]],
  sex    = recode_sex(raw[[COL_SEX]]),
  dob    = parse_date(raw[[COL_BIRTH]]),
  entry  = parse_date(raw[[COL_ENTRY]]),
  exit   = parse_date(raw[[COL_EXIT]]),
  died   = as.integer(tolower(as.character(raw[[COL_STATUS]])) %in% tolower(DEATH_CODES))
)

## Optional study-window clamp (right-censor at STUDY_END, delay entry to STUDY_START)
if (!is.null(STUDY_START)) coh[entry < STUDY_START, entry := STUDY_START]
if (!is.null(STUDY_END)) {
  coh[exit > STUDY_END, `:=`(exit = STUDY_END, died = 0L)]  # death after window = censored
}

## ---- Data-quality checks ---------------------------------------------------
prob <- coh[is.na(sex) | is.na(dob) | is.na(entry) | is.na(exit) | exit < entry]
if (nrow(prob) > 0) {
  warning(sprintf("Dropping %d rows with bad/missing dates, unmapped sex, or exit<entry.",
                  nrow(prob)))
  fwrite(prob, file.path(OUT_DIR, "smr_excluded_rows.csv"))
}
coh <- coh[!is.na(sex) & !is.na(dob) & !is.na(entry) & !is.na(exit) & exit >= entry]
stopifnot(nrow(coh) > 0)

## ---- Same-day follow-up (exit == entry) ------------------------------------
## Two participants may have observation entry and exit on the SAME day (e.g. a
## death dated the day observation starts). Such rows produce zero-length
## follow-up, so the Lexis split drops them AND their death silently vanishes.
## SAMEDAY controls the handling (documented in README + methods):
##   "keep"    -> credit 1 day of follow-up so the death is counted (default,
##                recommended for an SMR: never silently drop a known death).
##   "exclude" -> drop these rows entirely (use only if same-day dates are a
##                data artifact / landmark-survival requirement, and report it).
SAMEDAY <- "keep"
n_sameday <- coh[exit == entry, .N]
if (n_sameday > 0L) {
  if (SAMEDAY == "keep") {
    coh[exit == entry, exit := entry + 1]
    cat(sprintf("Same-day follow-up: %d row(s) given 1 day so the death is counted.\n", n_sameday))
  } else if (SAMEDAY == "exclude") {
    coh <- coh[exit > entry]
    cat(sprintf("Same-day follow-up: %d row(s) EXCLUDED (SAMEDAY='exclude').\n", n_sameday))
  }
}

## ============================================================================
## 4. LEXIS EXPANSION  --  split each subject's person-time into
##    (age-in-completed-years x calendar-year) cells
## ============================================================================
## For every subject we walk the timeline from entry to exit, cutting at each
## birthday and at each Jan-1 boundary. Each resulting segment contributes
## person-years to exactly one (age, year) cell; the death (if any) is credited
## to the final segment.

DAY <- 1 / 365.25

split_one <- function(id, sex, dob, entry, exit, died) {
  # candidate cut points: entry, exit, every birthday, every Jan 1st
  yrs  <- seq(as.integer(format(entry, "%Y")), as.integer(format(exit, "%Y")))
  jan1 <- as.Date(paste0(yrs, "-01-01"), format = "%Y-%m-%d")
  # birthdays within the interval
  bdays <- as.Date(paste0(yrs, "-", format(dob, "%m-%d")), format = "%Y-%m-%d")
  bdays <- bdays[!is.na(bdays)]   # Feb-29 birthdays are NA in non-leap years -> drop
  cuts  <- sort(unique(c(entry, exit, jan1, bdays)))
  cuts  <- cuts[cuts >= entry & cuts <= exit]
  n <- length(cuts) - 1L
  if (n < 1L) return(NULL)
  seg_start <- cuts[1:n]
  seg_end   <- cuts[2:(n + 1L)]
  py   <- as.numeric(seg_end - seg_start) * DAY
  age  <- as.integer(floor(as.numeric(seg_start - dob) * DAY))
  year <- as.integer(format(seg_start, "%Y"))
  ev   <- integer(n); if (died == 1L) ev[n] <- 1L   # death credited to last segment
  data.table(id, sex, age, year, py, ev)
}

message("Expanding person-time (Lexis split)...")
lex <- rbindlist(lapply(seq_len(nrow(coh)), function(i)
  split_one(coh$id[i], coh$sex[i], coh$dob[i], coh$entry[i], coh$exit[i], coh$died[i])
))
lex <- lex[py > 0]

## ============================================================================
## 5. MERGE REFERENCE RATES  ->  EXPECTED DEATHS
## ============================================================================
## Clamp ages above the max in the rates table to that max (open-ended top age).
maxage <- rates[, max(age)]
lex[age > maxage, age := maxage]
## Clamp calendar years to the range covered by the rates table (nearest year).
yrng <- rates[, range(year)]
lex[year < yrng[1], year := yrng[1]]
lex[year > yrng[2], year := yrng[2]]

lex <- merge(lex, rates[, .(sex, age, year, mx)],
             by = c("sex", "age", "year"), all.x = TRUE)
if (anyNA(lex$mx)) warning(sprintf("%d person-time cells had no matching rate (set to 0 expected).",
                                   sum(is.na(lex$mx))))
lex[is.na(mx), mx := 0]
lex[, expected := py * mx]                 # expected deaths in this cell

## ============================================================================
## 6. SMR  --  overall, by sex, and by age band
## ============================================================================

smr_from <- function(dt, by = NULL) {
  agg <- dt[, .(observed = sum(ev), expected = sum(expected),
                person_years = sum(py)), by = by]
  agg[, smr := observed / expected]
  ci <- t(mapply(pois_ci, agg$observed))
  agg[, `:=`(smr_lo = ci[, 1] / expected,
             smr_hi = ci[, 2] / expected)]
  # two-sided mid-p-ish exact Poisson p-value vs SMR = 1
  agg[, p_value := mapply(function(O, E) {
        if (E <= 0) return(NA_real_)
        p <- if (O >= E) ppois(O - 1, E, lower.tail = FALSE) else ppois(O, E)
        min(1, 2 * p)
      }, observed, expected)]
  agg[]
}

overall  <- smr_from(lex)
by_sex   <- smr_from(lex, "sex")

lex[, age_band := cut(age, breaks = AGE_BREAKS, right = FALSE,
                      include.lowest = TRUE)]
by_age   <- smr_from(lex, c("sex", "age_band"))

## ============================================================================
## 7. OUTPUT
## ============================================================================

fmt <- function(d) {
  d <- copy(d)
  numcols <- c("expected","smr","smr_lo","smr_hi","person_years","p_value")
  for (c in intersect(numcols, names(d))) d[[c]] <- round(d[[c]], 3)
  d
}

cat("\n================ OVERALL SMR ================\n")
print(fmt(overall))
cat("\n================ SMR BY SEX =================\n")
print(fmt(by_sex))
cat("\n============ SMR BY SEX x AGE BAND ==========\n")
print(fmt(by_age))

fwrite(fmt(overall), file.path(OUT_DIR, "smr_overall.csv"))
fwrite(fmt(by_sex),  file.path(OUT_DIR, "smr_by_sex.csv"))
fwrite(fmt(by_age),  file.path(OUT_DIR, "smr_by_age.csv"))
fwrite(lex,          file.path(OUT_DIR, "smr_person_time_cells.csv"))

cat("\nWrote: smr_overall.csv, smr_by_sex.csv, smr_by_age.csv, smr_person_time_cells.csv\n")
