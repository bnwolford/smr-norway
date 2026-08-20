################################################################################
## plot_smr_gg()  --  ggplot2 forest plot from SMR output
## ---------------------------------------------------------------------------
## ggplot2 version of the SMR forest plot. Points are drawn semi-transparent
## (alpha = 0.5 by default) so where markers and confidence intervals overlap
## between sexes you can still see the interval underneath.
##
## Plug at the END of any SMR script (after `by_age` exists), or run standalone
## against the saved CSV.
##
## USAGE:
##     source("smr_plot_gg.R")
##     p <- plot_smr_gg(by_age)                      # from in-memory object
##     p <- plot_smr_gg("smr_by_age.csv")            # or from the saved file
##     ggsave("smr_forest_gg.png", p, width = 8, height = 6, dpi = 300)
##  (the function also writes the file directly if you pass `outfile=`)
################################################################################

plot_smr_gg <- function(by_age,
                        outfile   = NULL,        # if given, ggsave() to this path
                        sex_col   = "sex",
                        band_col  = "age_band",
                        cols      = c(Female = "#C44E52", Male = "#4C72B0"),
                        point_alpha = 0.5,       # <-- transparency so CIs show through
                        line_alpha  = 0.6,
                        floor_smr = 0.04,        # left edge for zero-death arrows
                        size_range = c(1.5, 7),  # marker size range (~ observed deaths)
                        width = 8, height = 6, dpi = 300) {

  suppressPackageStartupMessages({ library(ggplot2) })

  ## ---- accept a data.frame/data.table OR a path to the CSV -----------------
  if (is.character(by_age) && length(by_age) == 1L && file.exists(by_age))
    by_age <- utils::read.csv(by_age, stringsAsFactors = FALSE)
  d <- as.data.frame(by_age, stringsAsFactors = FALSE)

  need <- c(sex_col, band_col, "observed", "expected", "smr", "smr_lo", "smr_hi")
  miss <- setdiff(need, names(d))
  if (length(miss)) stop("by_age is missing column(s): ", paste(miss, collapse = ", "))

  ## ---- order age bands by lower bound, build pretty labels -----------------
  low <- as.numeric(sub("^[^0-9-]*(-?[0-9.]+).*", "\\1", as.character(d[[band_col]])))
  band_order <- unique(d[[band_col]][order(low)])
  pretty_band <- function(b) {
    b <- as.character(b)
    lo <- sub("^[^0-9-]*(-?[0-9.]+).*", "\\1", b)
    hi <- sub(".*,\\s*([0-9.A-Za-z]+)[^0-9.A-Za-z]*$", "\\1", b)
    ifelse(grepl("Inf", hi), paste0(lo, "+"), paste0(lo, "-", hi))
  }
  lab_map <- setNames(pretty_band(band_order), as.character(band_order))
  # factor with youngest at TOP  (reverse, because ggplot y grows upward)
  d$band_f <- factor(lab_map[as.character(d[[band_col]])],
                     levels = rev(unname(lab_map[as.character(band_order)])))
  d$sex_f  <- d[[sex_col]]

  ## ---- split zero-death strata (drawn as arrows) from the rest -------------
  d$zero <- !is.na(d$observed) & d$observed == 0
  pts  <- d[!d$zero & !is.na(d$smr), ]
  zero <- d[d$zero, ]
  if (nrow(zero)) zero$smr_lo_arrow <- floor_smr    # arrow tail

  sexes <- intersect(names(cols), unique(d$sex_f))
  if (!length(sexes)) { sexes <- unique(d$sex_f)
    cols <- setNames(c("#C44E52", "#4C72B0", "#55A868")[seq_along(sexes)], sexes) }
  dodge <- ggplot2::position_dodge(width = 0.55)

  brks <- c(0.1, 0.25, 0.5, 1, 2, 4, 8, 16)
  xmax <- max(4, max(d$smr_hi[is.finite(d$smr_hi)], na.rm = TRUE))
  brks <- brks[brks >= floor_smr * 0.8 & brks <= xmax * 1.2]

  p <- ggplot2::ggplot(pts, ggplot2::aes(x = smr, y = band_f, colour = sex_f)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", colour = "grey35") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = smr_lo, xmax = smr_hi),
                           orientation = "y", width = 0, linewidth = 0.7,
                           alpha = line_alpha, position = dodge) +
    ggplot2::geom_point(ggplot2::aes(size = observed),
                        alpha = point_alpha, position = dodge) +
    ggplot2::scale_size_area(max_size = size_range[2], guide = "none")

  ## zero-death arrows (open points at the upper CI, arrow toward the floor)
  if (nrow(zero)) {
    p <- p +
      ggplot2::geom_segment(data = zero,
                            ggplot2::aes(x = smr_hi, xend = smr_lo_arrow,
                                         y = band_f, yend = band_f, colour = sex_f),
                            arrow = grid::arrow(length = grid::unit(0.06, "inches")),
                            linewidth = 0.6, alpha = line_alpha,
                            position = dodge, inherit.aes = FALSE) +
      ggplot2::geom_point(data = zero,
                          ggplot2::aes(x = smr_hi, y = band_f, colour = sex_f),
                          shape = 21, fill = "white", alpha = 1,
                          position = dodge, inherit.aes = FALSE)
  }

  p <- p +
    ggplot2::scale_x_log10(breaks = brks, labels = brks,
                           limits = c(floor_smr * 0.8, max(xmax, 8))) +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::labs(x = "Standardized mortality ratio (log scale)",
                  y = "Age band (years)",
                  title = "SMR vs. reference population, by age and sex") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = c(0.9, 0.12),
                   legend.background = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 12, hjust = 0),
                   axis.ticks.length = ggplot2::unit(-2, "pt"))

  if (!is.null(outfile))
    ggplot2::ggsave(outfile, p, width = width, height = height, dpi = dpi)
  p
}
