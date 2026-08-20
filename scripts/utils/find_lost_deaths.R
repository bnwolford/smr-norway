## Paste in AFTER sourcing the analysis script. Finds the deaths in `coh`
## that never made it into `lex` (the 565 -> 563 loss at n_lex).

died_ids <- coh[died == 1L, id]
lex_death_ids <- lex[ev == 1L, unique(id)]
lost <- setdiff(died_ids, lex_death_ids)
cat(sprintf("deaths flagged in coh : %d\n", length(died_ids)))
cat(sprintf("deaths landing in lex : %d\n", length(lex_death_ids)))
cat(sprintf("LOST deaths           : %d\n\n", length(lost)))

if (length(lost)) {
  info <- coh[id %in% lost, .(id, sex, dob, entry, exit, died,
                              follow_days = as.numeric(exit - entry))]
  print(info)
  ## classify each lost death. Only two real causes at this step:
  ##  (a) zero-length follow-up: exit == entry -> split_one() returns NULL
  ##  (b) age floor (restricted script): death at attained age < AGE_MIN
  age_at_exit <- info[, floor(as.numeric(exit - dob) / 365.25)]
  info[, age_exit := age_at_exit]
  info[, reason := fifelse(exit == entry, "zero follow-up (entry==exit) -> NULL",
                    fifelse(exists("AGE_MIN") && age_at_exit < AGE_MIN,
                            "death at age < AGE_MIN -> removed by age filter",
                            "unexpected -- inspect this row"))]
  cat("\nclassification:\n"); print(info[, .(id, follow_days, age_exit, reason)])
}
