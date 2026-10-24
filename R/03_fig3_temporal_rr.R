## Figure 3: monthly relative risk of outbreak size, Baseline (main) and national-trend model (sensitivity); colour-blind-safe (Okabe-Ito)
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
infile  <- file.path(TAB_DIR, "Table_temporalRR_monthly.csv")
outdir  <- FIG_DIR
MERGE_GAP <- 4L
d <- fread(infile, encoding = "UTF-8")
d[, `:=`(t = .I, month = match(Month, month.abb))]
d[, date := as.Date(sprintf("%d-%02d-15", Year, month))]
long <- rbind(
  d[, .(t, date, Year, month, RR = `RR (primary)`,  ci = `95% CrI (primary)`,  sig = `Signal (primary)`,  model = "Baseline model (main analysis)")],
  d[, .(t, date, Year, month, RR = `RR (national)`, ci = `95% CrI (national)`, sig = `Signal (national)`, model = "National-trend model without iid month effect (sensitivity)")])
long[, c("lo", "hi") := tstrsplit(ci, "-", fixed = TRUE, type.convert = TRUE)]
long[, model := factor(model, levels = unique(model))]
long[, status := fifelse(grepl("elevated", sig), "Credibly higher", fifelse(grepl("lower|suppress", sig), "Credibly lower", NA_character_))]
## periods: credibly higher months <= MERGE_GAP apart merged
per <- long[status == "Credibly higher"][order(model, t)][, grp := cumsum(c(1, diff(t)) > MERGE_GAP), by = model][
  , .(start = min(date) - 15, end = max(date) + 15, pk = date[which.max(RR)], pkRR = max(RR),
      pklo = lo[which.max(RR)], pkhi = hi[which.max(RR)]), by = .(model, grp)]
per[, lab := sprintf("Peak %s %d\nRR %.2f (%.2f to %.2f)", month.abb[as.integer(format(pk, "%m"))], as.integer(format(pk, "%Y")), pkRR, pklo, pkhi)]
print(per)
## Ethiopian seasons as background bands
seas <- data.table(date = seq(as.Date("2015-01-01"), as.Date("2024-12-01"), by = "month"))
seas[, m := as.integer(format(date, "%m"))]
seas[, season := fifelse(m %in% c(12, 1, 2), "Bega (dry)", fifelse(m %in% 3:5, "Belg (short rains)", fifelse(m %in% 6:8, "Kiremt (main rains)", "Tseday (post-rains)")))]
seas[, `:=`(xmin = date, xmax = seq(as.Date("2015-02-01"), as.Date("2025-01-01"), by = "month"))]
seas_cols <- c("Bega (dry)" = "#F0E442", "Belg (short rains)" = "#E69F00", "Kiremt (main rains)" = "#56B4E9", "Tseday (post-rains)" = "#009E73")   # Okabe-Ito
ymax <- 8
p <- ggplot() +
  geom_rect(data = seas, aes(xmin = xmin, xmax = xmax, ymin = 0.05, ymax = Inf, fill = season), alpha = 0.13) +
  geom_rect(data = per, aes(xmin = start, xmax = end, ymin = 0.05, ymax = Inf), fill = NA, colour = "grey20", linewidth = 0.25, linetype = "dotted") +
  geom_ribbon(data = long, aes(x = date, ymin = lo, ymax = pmin(hi, ymax)), fill = "grey45", alpha = 0.25) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30", linewidth = 0.4) +
  geom_line(data = long, aes(x = date, y = RR), colour = "black", linewidth = 0.8) +
  geom_point(data = long[!is.na(status)], aes(x = date, y = RR, colour = status, shape = status), size = 1.8) +
  geom_text(data = per, aes(x = pk, y = 6.2, label = lab), size = 2.6, lineheight = 0.9, family = "Liberation Serif", colour = "grey15") +
  scale_fill_manual(values = seas_cols, name = "Ethiopian season") +
  scale_colour_manual(values = c("Credibly higher" = "#0072B2", "Credibly lower" = "#D55E00"), name = "95% CrI excludes 1", na.translate = FALSE) +
  scale_shape_manual(values = c("Credibly higher" = 16, "Credibly lower" = 17), name = "95% CrI excludes 1", na.translate = FALSE) +
  scale_y_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 8), labels = c("0.1", "0.25", "0.5", "1", "2", "4", "8")) +
  coord_cartesian(ylim = c(0.08, ymax)) +
  scale_x_date(breaks = seq(as.Date("2015-01-01"), as.Date("2024-01-01"), by = "year"), date_labels = "%Y", limits = c(as.Date("2015-01-01"), as.Date("2024-12-31")), expand = c(0, 0)) +
  facet_wrap(~ model, ncol = 1) +
  labs(x = NULL, y = "Rate ratio of outbreak size (log scale)",
       caption = "RR = rate ratio of the temporal trend relative to the 2015-2024 average month, after seasonality, conditional on an outbreak. Ribbon: 95% CrI. Dotted frames: periods of credibly larger outbreak size (credibly higher months up to 4 months apart merged).") +
  theme_bw(base_size = 10, base_family = "Liberation Serif") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "grey94"), strip.text = element_text(size = 10, hjust = 0),
        legend.position = "bottom", legend.box = "vertical", legend.margin = margin(0, 0, 0, 0),
        plot.caption = element_text(size = 7.5, colour = "grey30", hjust = 0)) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(alpha = 0.35)), colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))
ggsave(file.path(outdir, "Fig3_temporal_RR.png"),  p, width = 10, height = 6.4, dpi = 300)
ggsave(file.path(outdir, "Fig3_temporal_RR.tiff"), p, width = 10, height = 6.4, dpi = 300, compression = "lzw")
