################################################################################
## fetch_ssb_rates.R
##
## Reproducibly (re)build data/norway_reference_mortality_rates.csv from the
## authoritative source: Statistics Norway (SSB) StatBank table 07902,
## "Life tables, by sex, age x, contents and year".
##
##   https://www.ssb.no/en/statbank/table/07902
##   API:  https://data.ssb.no/api/v0/en/table/07902   (JSON-stat2)
##
## Requires network access to data.ssb.no and packages: httr, jsonlite, data.table.
## The committed CSV covers years 2000-2024, ages 0-106, Male & Female. The table
## itself goes back to 1966 -- change YEARS below to extend it.
##
## Run from the repository root:  Rscript data-raw/fetch_ssb_rates.R
################################################################################

library(httr)
library(jsonlite)
library(data.table)

## ---- What to pull ----------------------------------------------------------
YEARS   <- as.character(2000:2024)   # extend back to "1966" if your cohort predates 2000
OUTFILE <- "data/norway_reference_mortality_rates.csv"
TABLE   <- "07902"
URL     <- sprintf("https://data.ssb.no/api/v0/en/table/%s", TABLE)

## ---- Build the JSON-stat2 query --------------------------------------------
## Kjonn: 1=Males, 2=Females.  ContentsCode: Dodssannsynlighet = probability of
## death per 1000 (qx); Dode = deaths per 100k (dx), kept for traceability.
query <- list(
  query = list(
    list(code = "Kjonn",        selection = list(filter = "item", values = list("1","2"))),
    list(code = "AlderX",       selection = list(filter = "all",  values = list("*"))),
    list(code = "ContentsCode", selection = list(filter = "item",
                                                 values = list("Dodssannsynlighet","Dode"))),
    list(code = "Tid",          selection = list(filter = "item", values = as.list(YEARS)))
  ),
  response = list(format = "json-stat2")
)

message("Querying SSB table ", TABLE, " for years ", YEARS[1], "-", tail(YEARS,1), " ...")
resp <- POST(URL, body = toJSON(query, auto_unbox = TRUE), content_type_json(),
             timeout(120))
stop_for_status(resp)
js <- fromJSON(content(resp, "text", encoding = "UTF-8"), simplifyVector = FALSE)

## ---- Parse JSON-stat2 into a tidy table ------------------------------------
## Values are a flat array in row-major order over the dimensions listed in js$id,
## each of size js$size. Reconstruct the index combinations in the same order.
dims  <- unlist(js$id)
sizes <- unlist(js$size)

## For each dimension, the category codes in position order.
codes_in_order <- function(dim_name) {
  idx <- js$dimension[[dim_name]]$category$index      # code -> position (0-based)
  nm  <- names(idx)
  nm[order(unlist(idx))]
}
ord <- lapply(dims, codes_in_order); names(ord) <- dims

## Cartesian product in row-major order (last dimension varies fastest).
grid <- do.call(CJ, c(rev(ord), sorted = FALSE))       # data.table cross join
setcolorder(grid, rev(names(grid)))                    # restore dims order
setnames(grid, dims)
grid[, value := unlist(js$value)]

## ---- Reshape to sex x age x year with qx and mx ----------------------------
grid[, sex  := c("1"="Male","2"="Female")[Kjonn]]
grid[, age  := as.integer(AlderX)]
grid[, year := as.integer(Tid)]
wide <- dcast(grid, sex + age + year ~ ContentsCode, value.var = "value")
setnames(wide, c("Dodssannsynlighet","Dode"), c("qx_per_1000","dx_per_100k"))

wide[, qx := qx_per_1000 / 1000]                       # probability of death in year
wide[, mx := qx / (1 - qx/2)]                          # central rate (uniform deaths)
out <- wide[, .(sex, age, year, qx, mx, qx_per_1000)][order(sex, year, age)]

dir.create(dirname(OUTFILE), showWarnings = FALSE, recursive = TRUE)
fwrite(out, OUTFILE)
message(sprintf("Wrote %d rows to %s (ages %d-%d, years %d-%d, sexes %s).",
                nrow(out), OUTFILE, min(out$age), max(out$age),
                min(out$year), max(out$year), paste(unique(out$sex), collapse=", ")))
