## Supplementary figure: the four HMM indicators by latent profile (each dot = one zone)
##   Rscript Fig_S_profile_indicators.R send_stage3/Table_S_zone_profile_community_class.csv send_stage3
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
args <- commandArgs(TRUE); infile <- args[1]; outdir <- if (length(args) >= 2) args[2] else "."
z <- fread(infile)
lab <- c("Higher burden,\nlarge peaks", "Low burden,\nsmall peaks", "Little or no\nrecorded activity")
z[, prof_n := paste0(lab[profile], "\n(n = ", .N, ")"), by = profile]
z[, prof_n := factor(prof_n, unique(prof_n[order(profile)]))]
z[, peak_log := log10(pmax(peak_reached, 0.01))]
L <- melt(z, id.vars = c("zone_name", "prof_n"), measure.vars = c("burden", "entry", "persistence_active", "peak_log"))
L[, variable := factor(variable, c("burden", "entry", "persistence_active", "peak_log"),
                       c("Burden\n(share of months in high state)", "Entry\n(monthly probability of becoming active)",
                         "Persistence of active months\n(probability of staying active)", "Peak reached\n(log10 cases per month)"))]
p <- ggplot(L, aes(prof_n, value, fill = prof_n)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6) +
  geom_jitter(width = 0.15, height = 0, size = 1.3, alpha = 0.8, shape = 21, colour = "grey20") +
  facet_wrap(~variable, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("#BD0026", "#FD8D3C", "#FFEDA0"), guide = "none") +
  labs(x = NULL, y = NULL, title = "HMM indicators by latent profile",
       subtitle = "Profiles were formed from burden and peak reached; entry and persistence are shown for description only",
       caption = "Each dot is one zone (posterior mean). Boxes show the median and interquartile range.") +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 9.5), strip.background = element_rect(fill = "grey95"),
        plot.title = element_text(face = "bold"), axis.text.x = element_text(size = 8.5))
ggsave(file.path(outdir, "Fig_S_profile_indicators.png"), p, width = 13, height = 5, dpi = 600, bg = "white")
z[, .(n = .N, burden = median(burden), entry = median(entry), persistence = median(persistence_active), peak = median(peak_reached)), by = .(profile)][order(profile)] |> print()
cat("Saved\n")
