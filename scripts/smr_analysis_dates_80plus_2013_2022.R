################################################################################
## SMR vs. the Norwegian general population
##   ---- RESTRICTED VARIANT: attained age >= 80 AND calendar years 2013-2022
##
## This is smr_analysis_dates.R with two restrictions applied:
##
##   (1) CALENDAR PERIOD 2013-2022  -- via STUDY_START / STUDY_END.
##       Each person's entry/exit is clamped to [2013-01-01, 2022-12-31].
##       A death occurring AFTER 2022-12-31 is converted to a censoring
##       event (died := 0) rather than a counted death. This is the correct
##       handling of left-truncation / right-censoring at the window edges.
##
##   (2) ATTAINED AGE >= 80         -- via a filter on the Lexis table.
##       Applied AFTER person-time expansion, so only the age-80+ slices of
##       each person's follow-up (and only deaths occurring at age 80+) are
##       kept. A person who enters at 78 and dies at 83 contributes their
##       80-83 person-time and their death, but NOT their 78-79 time.
##       (Filtering the INPUT rows instead would wrongly keep under-80 time.)
##
## "age >= 80" = attained age 80 and older (the [80,81) band and up), because
## `age` is the integer floor of exact age at the start of each segment.
## For strictly-over-80 (81+), set AGE_MIN <- 81.
##
## Cohort input still has THREE dates (start obs, end obs, death) plus age at
## entry; DOB is reconstructed from entry date and age at entry.
################################################################################

library(data.table)

## ============================================================================
## 1. USER SETTINGS
## ============================================================================
COHORT_FILE <- "cohort_dates.csv"
RATES_FILE  <- "norway_reference_mortality_rates.csv"
OUT_DIR     <- "."

## ---- Column names in YOUR file ---------------------------------------------
COL_ID        <- "id"
COL_SEX       <- "sex"
COL_ENTRY     <- "start_obs"    # date observation starts           (REQUIRED)
COL_ENDOBS    <- "end_obs"      # date observation ends / censored   (REQUIRED)
COL_DEATH     <- "death_date"   # date of death; blank/NA if alive   (REQUIRED)
COL_AGE_ENTRY <- "age_entry"    # age at start of observation        (REQUIRED, no DOB)
## If you HAVE date of birth instead of age, set COL_AGE_ENTRY <- NA and give:
COL_BIRTH     <- NA             # date of birth column name, or NA

## ---- Age handling (only used when reconstructing DOB from age) -------------
AGE_OFFSET <- 0.5   # completed years -> mid-year birthday; set 0 for exact age

## ---- RESTRICTIONS ----------------------------------------------------------
AGE_MIN      <- 80                        # keep only attained age >= AGE_MIN
STUDY_START  <- as.Date("2013-01-01")     # calendar window start (inclusive)
STUDY_END    <- as.Date("2022-12-31")     # calendar window end   (inclusive)

## ---- Sex / dates -----------------------------------------------------------
SEX_MAP <- list(
  Male   = c("Male",   "M", "1", "male",   "m"),
  Female = c("Female", "F", "2", "female", "f")
)
DATE_FORMATS <- c("%Y-%m-%d", "%d.%m.%Y", "%d/%m/%Y", "%m/%d/%Y")
## Age bands: default tops out at 85. Because we restrict to 80+, the only
## bands that appear are [80,85) and [85,Inf). For finer resolution among the
## very old, use e.g. AGE_BREAKS <- c(80, 85, 90, 95, Inf)
AGE_BREAKS   <- c(seq(0, 85, by = 5), Inf)

## ============================================================================
## 2. HELPERS
## ============================================================================
parse_date <- function(x) {
  if (inherits(x, "Date")) return(as.Date(x))
  x <- as.character(x); x[x %in% c("", "NA", "NaN", ".")] <- NA
  out <- as.Date(rep(NA_real_, length(x)), origin = "1970-01-01")
  for (fmt in DATE_FORMATS) { miss <- is.na(out) & !is.na(x); if (!any(miss)) break
    out[miss] <- as.Date(x[miss], format = fmt) }
  out
}
recode_sex <- function(x) {
  x <- as.character(x); out <- rep(NA_character_, length(x))
  for (target in names(SEX_MAP)) out[x %in% SEX_MAP[[target]]] <- target
  out
}
pois_ci <- function(O, conf = 0.95) {
  a <- (1 - conf) / 2
  lo <- ifelse(O == 0, 0, qgamma(a, O, rate = 1)); hi <- qgamma(1 - a, O + 1, rate = 1)
  c(lower = lo, upper = hi)
}
YEAR <- 365.25

## ============================================================================
## 3. LOAD + DERIVE exit/event + RECONSTRUCT dob
## ============================================================================
rates <- fread(RATES_FILE); setkey(rates, sex, age, year)
raw   <- fread(COHORT_FILE)

entry   <- parse_date(raw[[COL_ENTRY]])
end_obs <- parse_date(raw[[COL_ENDOBS]])
death   <- parse_date(raw[[COL_DEATH]])

## ---- exit date + event flag ------------------------------------------------
## died only if a death date exists AND occurs on/before the observation close.
died <- as.integer(!is.na(death) & death <= end_obs)
## exit is the death date for those who died within the window, else obs close.
exit <- as.Date(ifelse(died == 1L, death, end_obs), origin = "1970-01-01")

## ---- date of birth ---------------------------------------------------------
if (!is.na(COL_BIRTH) && COL_BIRTH %in% names(raw)) {
  dob <- parse_date(raw[[COL_BIRTH]])
} else {
  age_entry <- as.numeric(raw[[COL_AGE_ENTRY]])
  dob <- entry - (age_entry + AGE_OFFSET) * YEAR
}

coh <- data.table(
  id = raw[[COL_ID]], sex = recode_sex(raw[[COL_SEX]]),
  dob = dob, entry = entry, exit = exit, died = died
)

## ---- Apply CALENDAR-PERIOD restriction (2013-2022) -------------------------
## Clamp entry forward and exit back to the window. Deaths that fall after the
## window close are converted to censoring (died := 0). Anyone whose entire
## follow-up lies outside the window drops out at the exit < entry filter below.
if (!is.null(STUDY_START)) coh[entry < STUDY_START, entry := STUDY_START]
if (!is.null(STUDY_END))   coh[exit  > STUDY_END, `:=`(exit = STUDY_END, died = 0L)]

prob <- coh[is.na(sex) | is.na(dob) | is.na(entry) | is.na(exit) | exit < entry]
if (nrow(prob) > 0) {
  warning(sprintf("Dropping %d rows with bad/missing values or exit<entry.", nrow(prob)))
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
cat(sprintf("Cohort after calendar restriction: %d subjects, %d deaths, %.0f person-years span.\n",
            nrow(coh), coh[, sum(died)], coh[, sum(as.numeric(exit - entry)/YEAR)]))

## ============================================================================
## 4-7. IDENTICAL LEXIS / SMR ENGINE
## ============================================================================
DAY <- 1 / YEAR
split_one <- function(id, sex, dob, entry, exit, died) {
  yrs   <- seq(as.integer(format(entry, "%Y")), as.integer(format(exit, "%Y")))
  jan1  <- as.Date(paste0(yrs, "-01-01"), format = "%Y-%m-%d")
  bdays <- as.Date(paste0(yrs, "-", format(dob, "%m-%d")), format = "%Y-%m-%d")
  bdays <- bdays[!is.na(bdays)]                       # Feb-29 in non-leap years -> drop
  cuts  <- sort(unique(c(entry, exit, jan1, bdays)))
  cuts  <- cuts[cuts >= entry & cuts <= exit]
  n <- length(cuts) - 1L; if (n < 1L) return(NULL)
  ss <- cuts[1:n]; se <- cuts[2:(n + 1L)]
  data.table(id, sex, age = as.integer(floor(as.numeric(ss - dob) * DAY)),
             year = as.integer(format(ss, "%Y")), py = as.numeric(se - ss) * DAY,
             ev = `if`(died == 1L, c(rep(0L, n - 1L), 1L), integer(n)))
}
message("Expanding person-time (Lexis split)...")
lex <- rbindlist(lapply(seq_len(nrow(coh)), function(i)
  split_one(coh$id[i], coh$sex[i], coh$dob[i], coh$entry[i], coh$exit[i], coh$died[i])))
lex <- lex[py > 0]

maxage <- rates[, max(age)]; lex[age > maxage, age := maxage]; lex[age < 0, age := 0]
yrng <- rates[, range(year)]
lex[year < yrng[1], year := yrng[1]]; lex[year > yrng[2], year := yrng[2]]
lex <- merge(lex, rates[, .(sex, age, year, mx)], by = c("sex","age","year"), all.x = TRUE)
lex[is.na(mx), mx := 0]; lex[, expected := py * mx]

## ---- Apply ATTAINED-AGE restriction (>= AGE_MIN) ---------------------------
## Applied AFTER expansion so only the age-80+ person-time slices (and deaths
## occurring at age 80+) are retained.
lex <- lex[age >= AGE_MIN]
stopifnot(nrow(lex) > 0)
cat(sprintf("After age >= %d restriction: %.1f person-years, %d deaths at age %d+.\n",
            AGE_MIN, lex[, sum(py)], lex[, sum(ev)], AGE_MIN))

smr_from <- function(dt, by = NULL) {
  agg <- dt[, .(observed = sum(ev), expected = sum(expected),
                person_years = sum(py)), by = by]
  agg[, smr := observed / expected]
  ci <- t(mapply(pois_ci, agg$observed))
  agg[, `:=`(smr_lo = ci[,1]/expected, smr_hi = ci[,2]/expected)]
  agg[, p_value := mapply(function(O,E){ if(E<=0) return(NA_real_)
        p <- if(O>=E) ppois(O-1,E,lower.tail=FALSE) else ppois(O,E); min(1,2*p)},
        observed, expected)]
  agg[]
}
overall <- smr_from(lex); by_sex <- smr_from(lex, "sex")
lex[, age_band := cut(age, breaks = AGE_BREAKS, right = FALSE, include.lowest = TRUE)]
by_age <- smr_from(lex, c("sex","age_band"))

fmt <- function(d){ d<-copy(d)
  for(c in intersect(c("expected","smr","smr_lo","smr_hi","person_years","p_value"),names(d)))
    d[[c]] <- round(d[[c]],3); d }
cat("\n================ OVERALL SMR (age 80+, 2013-2022) ================\n"); print(fmt(overall))
cat("\n================ SMR BY SEX =================\n"); print(fmt(by_sex))
cat("\n============ SMR BY SEX x AGE BAND ==========\n"); print(fmt(by_age))
fwrite(fmt(overall), file.path(OUT_DIR,"smr_overall.csv"))
fwrite(fmt(by_sex),  file.path(OUT_DIR,"smr_by_sex.csv"))
fwrite(fmt(by_age),  file.path(OUT_DIR,"smr_by_age.csv"))
fwrite(lex,          file.path(OUT_DIR,"smr_person_time_cells.csv"))
cat("\nDone.\n")
