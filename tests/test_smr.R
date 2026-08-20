################################################################################
## test_smr.R  --  self-checks for the SMR toolkit
##
## Run from the repository root:  source("tests/test_smr.R")
##
## Checks:
##   1. Person-time conservation: sum of Lexis-segment person-years equals the
##      naive (undivided) follow-up time, to within floating-point tolerance.
##   2. Death conservation: every death flagged in the cleaned cohort is counted
##      exactly once in the Lexis table (with SAMEDAY = "keep").
##   3. Same-day handling: a forced entry==exit death is retained.
##   4. Output sanity: SMR == observed / expected; CI brackets the point estimate.
################################################################################

library(data.table)
ok <- function(cond, msg) {
  if (isTRUE(cond)) cat(sprintf("  PASS  %s\n", msg))
  else { cat(sprintf("  FAIL  %s\n", msg)); assign(".TEST_FAILED", TRUE, envir = .GlobalEnv) }
}
.TEST_FAILED <- FALSE

## ---- Build a small deterministic cohort with a forced same-day death --------
set.seed(42)
n <- 200
entry <- as.Date("2010-01-01") + sample(0:2000, n, replace = TRUE)
dur   <- sample(200:3000, n, replace = TRUE)
endob <- entry + dur
died  <- runif(n) < 0.4
deathd <- entry + sample(50:2800, n, replace = TRUE)
death_chr <- ifelse(died, format(deathd, "%Y-%m-%d"), "")
death_chr[1] <- format(entry[1], "%Y-%m-%d")   # force a same-day death
endob[1]     <- entry[1]                        # ensure it counts (within window)
coh_in <- data.table(
  id = 1:n, sex = ifelse(runif(n) < 0.5, "M", "F"),
  start_obs = format(entry, "%Y-%m-%d"),
  end_obs   = format(endob, "%Y-%m-%d"),
  death_date = death_chr,
  age_entry = sample(40:90, n, replace = TRUE)
)
fwrite(coh_in, "cohort_dates.csv")
file.copy("data/norway_reference_mortality_rates.csv",
          "norway_reference_mortality_rates.csv", overwrite = TRUE)

## ---- Run the analysis ------------------------------------------------------
source("scripts/smr_analysis_dates.R")

## ---- 1. Person-time conservation -------------------------------------------
naive_py <- coh[, sum(as.numeric(exit - entry) / 365.25)]
lex_py   <- lex[, sum(py)]
ok(abs(naive_py - lex_py) < 1e-6,
   sprintf("person-time conserved (naive %.4f == Lexis %.4f)", naive_py, lex_py))

## ---- 2 & 3. Death conservation and same-day retention ----------------------
ok(coh[, sum(died)] == lex[, sum(ev)],
   sprintf("deaths conserved (coh %d == lex %d)", coh[, sum(died)], lex[, sum(ev)]))
ok(1L %in% lex[ev == 1L, id], "forced same-day death (id 1) is retained")

## ---- 4. Output sanity ------------------------------------------------------
ok(abs(overall$smr - overall$observed / overall$expected) < 1e-9,
   "overall SMR == observed / expected")
ok(overall$smr_lo <= overall$smr && overall$smr <= overall$smr_hi,
   "overall 95% CI brackets the point estimate")

## ---- Cleanup ---------------------------------------------------------------
file.remove("cohort_dates.csv", "norway_reference_mortality_rates.csv")

if (isTRUE(.TEST_FAILED)) stop("Some tests FAILED (see above).") else
  cat("\nAll tests passed.\n")
