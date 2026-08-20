# Methods

## Standardized mortality ratio

Mortality in the cohort was compared with that expected in the general
Norwegian population using the standardized mortality ratio (SMR), defined as
the ratio of the number of observed deaths to the number of expected deaths,
SMR = O/E. Expected deaths were obtained by the indirect standardization (person-years)
method, applying age-, sex- and calendar-year-specific mortality rates from the
general population to the person-time accrued by cohort members.

## Reference mortality rates

Sex-, single-year-of-age- and calendar-year-specific mortality rates for the
Norwegian general population were obtained from the official life tables
published by Statistics Norway (Statistisk sentralbyrå; StatBank table 07902).
Rates were retrieved for single years of age from 0 to 106 years and for each
calendar year spanning the follow-up period. The published annual probability
of death at age *x*, q(x), was converted to the central mortality rate m(x)
used for the expected-count calculation under the standard assumption of a
uniform distribution of deaths within each year of age,
m(x) = q(x) / (1 − q(x)/2).

## Follow-up and person-time

Each participant contributed person-time from the date of entry into the cohort
until the date of exit, defined as the earliest of death, loss to follow-up, or
the end of the observation period; the latter two were treated as censoring.
Deaths recorded after the close of observation were censored at the close date.
Person-time was partitioned into age-, sex- and calendar-year-specific cells by
exact Lexis expansion: each participant's follow-up interval was subdivided at
every attained birthday and at every 1 January boundary, so that each resulting
segment fell entirely within a single combination of one-year age band and one
calendar year. Each segment contributed its exact duration (in years) of
person-time to the corresponding cell, and any death was assigned to the final
segment of the participant's follow-up. This procedure accounts correctly for
delayed entry (left truncation), for participants ageing through multiple age
bands during follow-up, and for secular change in population mortality across
calendar years.

## Expected deaths and the SMR

Within each age × sex × calendar-year cell, the expected number of deaths was
computed as the product of the accrued person-years and the corresponding
Norwegian reference rate. Expected deaths were summed, and the SMR was
calculated as the ratio of total observed to total expected deaths, overall and
within strata defined by sex and by five-year age band. An SMR greater than 1.0
indicates higher mortality in the cohort than in the age- and sex-matched
general population, whereas an SMR below 1.0 indicates lower mortality.

## Statistical inference

The observed number of deaths was assumed to follow a Poisson distribution.
Exact 95% confidence intervals for each SMR were derived from the exact Poisson
confidence limits for the observed count, divided by the expected count
(Garwood method). Two-sided *p*-values for the null hypothesis SMR = 1 were
obtained from the Poisson distribution with mean equal to the expected number of
deaths.

Analyses were performed in R version 4.5 (R Foundation for Statistical
Computing, Vienna, Austria).
