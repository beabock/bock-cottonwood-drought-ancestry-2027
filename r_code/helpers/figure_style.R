# ==============================================================================
# Shared figure style for every figure in the manuscript and appendices (sourced by r_code/28_figures.R).
#   - theme_bw() base, Arial, 8 pt text (9 pt axis titles, 10 pt bold panel letters), AJB sizes.
#   - One color rule: color always means watering treatment (blue = well-watered, gold = drought), with shape and
#     line type as redundant cues (circle and solid line = well-watered; triangle and dashed line = drought), so the
#     figures read the same in color, under color-vision deficiency and in grayscale. Anything else is gray or black.
#   - The blue/gold pair stays distinct under simulated deuteranopia, protanopia and tritanopia (CIELAB distance 79-106).
# ==============================================================================
CM <- 1 / 2.54
W_1COL <- 8.9 * CM; W_2COL <- 18.4 * CM; H_MAX <- 23 * CM
FONT_PT <- 8; FONT <- "Arial"
X_ANC <- "Fremont ancestry (%)"
PT_SIZE <- 1.3; PT_ALPHA <- 0.5; PT_STROKE <- 0.4; PT_EDGE <- "#444444"; LN_W <- 0.8; BAND_ALPHA <- 0.2
TRT_LEVELS <- c("Well-watered", "Drought")
TRT_COL <- c("Well-watered" = "#3B9AB2", "Drought" = "#F2AD00")
TRT_SHAPE <- c("Well-watered" = 21, "Drought" = 24)
TRT_LTY <- c("Well-watered" = "solid", "Drought" = "22")
NEUTRAL_FILL <- "grey80"; NEUTRAL_LINE <- "black"; NEUTRAL_BAND <- "grey40"

theme_fig <- function() {
  theme_bw(base_size = FONT_PT, base_family = FONT) +
    theme(text = element_text(size = FONT_PT, colour = "black"),
          axis.text = element_text(size = FONT_PT, colour = "black"),
          axis.title = element_text(size = FONT_PT + 1, colour = "black"),
          legend.text = element_text(size = FONT_PT), legend.title = element_text(size = FONT_PT),
          legend.position = "bottom", legend.key.size = unit(0.4, "cm"), legend.margin = margin(0, 0, 0, 0),
          strip.text = element_text(size = FONT_PT), strip.background = element_rect(fill = "grey92", colour = "grey50"),
          panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey88", linewidth = 0.3),
          plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(4, 6, 4, 4),
          plot.tag = element_text(size = 10, face = "bold", family = FONT),
          plot.subtitle = element_text(size = FONT_PT, margin = margin(0, 0, 2, 0)))
}
trt_scales <- function() list(
  scale_fill_manual(values = TRT_COL, name = "Treatment"), scale_colour_manual(values = TRT_COL, name = "Treatment"),
  scale_shape_manual(values = TRT_SHAPE, name = "Treatment"), scale_linetype_manual(values = TRT_LTY, name = "Treatment"))

# bold panel letters (and an optional plain label) inside faceted panels
panel_letters <- function(d, facet_var, label_var = NULL, x_label = NULL) {
  d <- d %>% mutate(lab = LETTERS[row_number()])
  out <- list(geom_text(data = d, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE, hjust = -0.5, vjust = 1.4,
                        size = 10 / ggplot2::.pt, fontface = "bold", family = FONT))
  if (!is.null(label_var)) out <- c(out, list(geom_text(data = d, aes(x = x_label, y = Inf, label = .data[[label_var]]), inherit.aes = FALSE,
                                                        vjust = 1.6, size = 9 / ggplot2::.pt, family = FONT)))
  out }

# writes <name>.pdf (vector, cairo) and <name>.png (600 dpi) at final size; draws the gtable inside each device
save_fig <- function(p, name, w, h, dir = "figures") {
  stopifnot(h <= H_MAX + 1e-9)
  draw <- function() { grid::grid.newpage(); grid::grid.draw(if (inherits(p, "patchwork")) patchwork::patchworkGrob(p) else ggplot2::ggplotGrob(p)) }
  grDevices::cairo_pdf(file.path(dir, paste0(name, ".pdf")), width = w, height = h, bg = "white"); draw(); invisible(dev.off())
  ragg::agg_png(file.path(dir, paste0(name, ".png")), width = w, height = h, units = "in", res = 600, background = "white"); draw(); invisible(dev.off())
  cat(sprintf("wrote %s (.pdf, .png): %.1f x %.1f cm\n", file.path(dir, name), w * 2.54, h * 2.54)) }
