################################################################################
## run_example.R  --  end-to-end demonstration on synthetic data
##
## Runs the three-date SMR analysis on the shipped synthetic cohort, prints the
## SMR tables, and saves a forest plot. Run from the repository root:
##
##     source("run_example.R")
##
## Produces (in the repo root / docs/figures):
##   smr_overall.csv, smr_by_sex.csv, smr_by_age.csv, smr_person_time_cells.csv
##   docs/figures/example_output_gg.png
################################################################################

## The analysis scripts read file paths from their USER SETTINGS block. For the
## demo we point those settings at the shipped example files by overriding them
## via a tiny in-place edit of the environment: we set the globals the script
## expects BEFORE sourcing is not possible (the script hard-codes them), so
## instead we temporarily copy the example inputs to the working directory names
## the script uses, then restore.

stopifnot(file.exists("scripts/smr_analysis_dates.R"))

## Stage the inputs under the names the script expects.
file.copy("examples/example_cohort_dates.csv", "cohort_dates.csv", overwrite = TRUE)
file.copy("data/norway_reference_mortality_rates.csv",
          "norway_reference_mortality_rates.csv", overwrite = TRUE)

## Run the analysis (writes smr_*.csv to the working directory).
source("scripts/smr_analysis_dates.R")

## Plot from the by-age output.
if (requireNamespace("ggplot2", quietly = TRUE)) {
  source("scripts/smr_plot_gg.R")
  p <- plot_smr_gg("smr_by_age.csv")
  dir.create("docs/figures", showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave("docs/figures/example_output_gg.png", p,
                  width = 8, height = 6, dpi = 150)
  cat("\nForest plot saved to docs/figures/example_output_gg.png\n")
} else {
  message("ggplot2 not installed; skipping plot. install.packages('ggplot2')")
}

## Clean up the staged input copies (outputs are left in place for inspection).
file.remove("cohort_dates.csv", "norway_reference_mortality_rates.csv")

cat("\nDemo complete. See smr_overall.csv / smr_by_sex.csv / smr_by_age.csv.\n")
