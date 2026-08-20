## Paste in AFTER sourcing the analysis script (uses coh, lex, and re-reads raw).
## Reconciles raw death dates -> counted deaths, stage by stage.
raw2   <- fread(COHORT_FILE)
d_raw  <- parse_date(raw2[[COL_DEATH]])
e_obs  <- parse_date(raw2[[COL_ENDOBS]])
ent    <- parse_date(raw2[[COL_ENTRY]])

n_datepresent <- sum(!is.na(d_raw))
n_afterclose  <- sum(!is.na(d_raw) & !is.na(e_obs) & d_raw > e_obs)
n_flagged     <- sum(coh$died)                     # deaths surviving row-clean + calendar clamp
n_lex         <- lex[, sum(ev)]                    # deaths surviving Lexis + py>0 (+ age filter)

cat(sprintf("death dates present in file      : %d\n", n_datepresent))
cat(sprintf("  - fall AFTER observation close : %d  (censored by design)\n", n_afterclose))
cat(sprintf("died==1 after clean + calendar    : %d\n", n_flagged))
cat(sprintf("observed in lex (final)           : %d\n", n_lex))
cat(sprintf("  >> lost between coh and lex      : %d\n", n_flagged - n_lex))

## Identify the same-day (entry==exit) decedents that split_one drops:
sameday <- coh[died == 1L & exit == entry]
cat(sprintf("\nsame-day entry==exit decedents (dropped by split_one): %d\n", nrow(sameday)))
if (nrow(sameday)) print(sameday[, .(id, sex, entry, exit)])

## If using the age-restricted script, deaths removed by the age floor:
if (exists("AGE_MIN")) {
  ## rebuild unfiltered lex would be needed for exact count; approximate via coh ages
  cat(sprintf("\n(age-restricted script: AGE_MIN = %d — deaths below this age are excluded)\n", AGE_MIN))
}
