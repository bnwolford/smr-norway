################################################################################
## plot_smr()  --  forest plot + observed/expected panel from SMR output
## ---------------------------------------------------------------------------
## Plug this at the END of any of the SMR scripts. It reads the by-age result
## the script already produced (the `by_age` object in memory, or the saved
## smr_by_age.csv) and writes a two-panel PNG:
##    A) forest plot of SMR by age band & sex, log axis, reference line at 1
##    B) observed vs expected deaths by age band
##
## Base R graphics only -- no extra packages.
##
## USAGE (after the script has run and `by_age` exists):
##     source("smr_plot.R")
##     plot_smr(by_age, outfile = "smr_forest.png")
## or from the saved file:
##     plot_smr("smr_by_age.csv", outfile = "smr_forest.png")
################################################################################

plot_smr <- function(by_age,
                     outfile   = "smr_forest.png",
                     sex_col   = "sex",
                     band_col  = "age_band",
                     cols      = c(Female = "#C44E52", Male = "#4C72B0"),
                     obs_col   = "#55A868",   # observed-deaths bar colour
                     exp_col   = "grey70",    # expected-deaths bar colour
                     floor_smr = 0.04,        # left edge for zero-death arrows
                     width     = 2760, height = 1560, res = 300) {

  ## ---- accept a data.frame/data.table OR a path to the CSV -----------------
  if (is.character(by_age) && length(by_age) == 1L && file.exists(by_age))
    by_age <- utils::read.csv(by_age, stringsAsFactors = FALSE)
  by_age <- as.data.frame(by_age, stringsAsFactors = FALSE)

  need <- c(sex_col, band_col, "observed", "expected", "smr", "smr_lo", "smr_hi")
  miss <- setdiff(need, names(by_age))
  if (length(miss)) stop("by_age is missing column(s): ", paste(miss, collapse = ", "))

  ## ---- order age bands by their lower bound --------------------------------
  low  <- as.numeric(sub("^[^0-9-]*(-?[0-9.]+).*", "\\1", as.character(by_age[[band_col]])))
  ord_bands <- unique(by_age[[band_col]][order(low)])
  nb   <- length(ord_bands)
  ypos <- rev(seq_len(nb))                       # top row = youngest band
  names(ypos) <- as.character(ord_bands)

  ## pretty band labels:  "[70,75)" -> "70-75",  "[85,Inf]" -> "85+"
  pretty_band <- function(b) {
    b <- as.character(b)
    lo <- sub("^[^0-9-]*(-?[0-9.]+).*", "\\1", b)
    hi <- sub(".*,\\s*([0-9.A-Za-z]+)[^0-9.A-Za-z]*$", "\\1", b)
    ifelse(grepl("Inf", hi), paste0(lo, "+"), paste0(lo, "-", hi))
  }
  band_lab <- pretty_band(ord_bands)

  sexes <- intersect(names(cols), unique(by_age[[sex_col]]))
  if (!length(sexes)) { sexes <- unique(by_age[[sex_col]])
    cols <- setNames(c("#C44E52", "#4C72B0", "#55A868")[seq_along(sexes)], sexes) }
  offs <- if (length(sexes) == 2) c(0.16, -0.16) else rep(0, length(sexes))
  names(offs) <- sexes

  ## ---- device --------------------------------------------------------------
  grDevices::png(outfile, width = width, height = height, res = res)
  on.exit(grDevices::dev.off(), add = TRUE)
  op <- graphics::par(mfrow = c(1, 2), mar = c(4, 6, 3, 1), mgp = c(2.2, 0.6, 0),
                      tcl = -0.3, cex.axis = 0.8, cex.lab = 0.95, xaxs = "i")
  on.exit(graphics::par(op), add = TRUE)

  ## ===== Panel A: forest plot (log x) =======================================
  xmax <- max(4, ceiling(max(by_age$smr_hi[is.finite(by_age$smr_hi)], na.rm = TRUE)))
  xlim <- c(floor_smr * 0.8, max(xmax, 8))
  graphics::plot(NA, xlim = xlim, ylim = c(0.4, nb + 0.6), log = "x",
                 xlab = "Standardized mortality ratio (log scale)",
                 ylab = "", yaxt = "n", xaxt = "n")
  ticks <- c(0.1, 0.25, 0.5, 1, 2, 4, 8, 16)
  ticks <- ticks[ticks >= xlim[1] & ticks <= xlim[2]]
  graphics::axis(1, at = ticks, labels = ticks)
  graphics::axis(2, at = ypos, labels = band_lab, las = 1)
  graphics::mtext("Age band (years)", side = 2, line = 4.2, cex = 0.95)
  graphics::abline(v = 1, lty = 2, col = "grey35")

  for (sx in sexes) {
    d <- by_age[by_age[[sex_col]] == sx, ]
    d <- d[match(as.character(ord_bands), as.character(d[[band_col]])), ]
    for (i in seq_len(nb)) {
      if (is.na(d$smr[i])) next
      y <- ypos[i] + offs[sx]; cc <- cols[[sx]]
      if (isTRUE(d$observed[i] == 0)) {
        graphics::arrows(d$smr_hi[i], y, floor_smr, y, length = 0.05,
                         col = cc, lwd = 1.2)
        graphics::points(d$smr_hi[i], y, pch = 21, bg = "white", col = cc, cex = 0.8)
      } else {
        graphics::segments(d$smr_lo[i], y, d$smr_hi[i], y, col = cc, lwd = 1.4)
        graphics::points(d$smr[i], y, pch = 19, col = cc,
                         cex = 0.5 + 0.35 * sqrt(d$observed[i]))
      }
    }
  }
  graphics::text(1, 0.55, "SMR = 1", pos = 4, offset = 0.2, cex = 0.7, col = "grey35")
  graphics::legend("bottomright", legend = sexes, pch = 19,
                   col = unlist(cols[sexes]), bty = "n", cex = 0.85)
  graphics::mtext("SMR vs. reference population, by age and sex",
                  side = 3, line = 0.6, adj = 0, cex = 1.0, font = 1)
  graphics::mtext("a", side = 3, line = 1.2, adj = -0.28, cex = 1.1, font = 2)

  ## ===== Panel B: observed vs expected deaths ===============================
  agg <- stats::aggregate(cbind(observed, expected) ~ get(band_col), data = by_age, FUN = sum)
  names(agg)[1] <- band_col
  agg <- agg[match(as.character(ord_bands), as.character(agg[[band_col]])), ]
  graphics::par(mar = c(4, 1, 3, 1))
  xb <- c(0, max(agg$expected, agg$observed, na.rm = TRUE) * 1.05)
  graphics::plot(NA, xlim = xb, ylim = c(0.4, nb + 0.6),
                 xlab = "Deaths", ylab = "", yaxt = "n")
  h <- 0.34
  for (i in seq_len(nb)) {
    graphics::rect(0, ypos[i] + 0.03, agg$expected[i], ypos[i] + 0.03 + h,
                   col = exp_col, border = NA)
    graphics::rect(0, ypos[i] - 0.03 - h, agg$observed[i], ypos[i] - 0.03,
                   col = obs_col, border = NA)
  }
  graphics::legend("bottomright", legend = c("Expected", "Observed"),
                   fill = c(exp_col, obs_col), border = NA, bty = "n", cex = 0.85)
  graphics::mtext("Observed vs. expected deaths", side = 3, line = 0.6,
                  adj = 0, cex = 1.0, font = 1)
  graphics::mtext("b", side = 3, line = 1.2, adj = -0.05, cex = 1.1, font = 2)

  invisible(outfile)
}
