## Figure 4: outbreak occurrence and expected outbreak size by zone, Baseline hurdle model (style of earlier Fig 4)
suppressPackageStartupMessages({library(data.table); library(sf); library(ggplot2); library(patchwork)})
sf_use_s2(FALSE)
infile   <- file.path(TAB_DIR, "TableS3_zone_summary_OCCPOP.csv")
shp_path <- "F:/Cholera/shape_22/Adm_2/eth_admbnda_adm2_csa_bofedb_2021.shp"
outdir   <- getOption("cholera.proj", ".")
FONT <- "sans"                                                                           # Windows: "Arial" also works

z   <- fread(infile)
shp <- st_read(shp_path, quiet = TRUE) |> st_make_valid() |> st_transform(20137)
shp$key <- toupper(gsub("[^A-Za-z0-9]", "", shp$ADM2_PCODE)); shp <- st_simplify(shp, dTolerance = 1000)
m <- merge(shp[, c("key", "geometry")], z, by.x = "key", by.y = "admin2_clean")
m$no_outbreak <- m$N_outbreak_months == 0
m$mu_plot <- ifelse(m$no_outbreak, NA, m$mu)          # grey in panel B
n_none <- sum(m$no_outbreak)

## north arrow and scale bar (drawn by hand, no ggspatial needed)
bb <- st_bbox(shp); W <- unname(bb["xmax"] - bb["xmin"]); H <- unname(bb["ymax"] - bb["ymin"])
## scale bar in the empty bottom-right corner, clear of the country outline
x0 <- unname(bb["xmax"]) - 4e5 - 0.03 * W; y0 <- unname(bb["ymin"]) + 0.02 * H; h <- 0.016 * H
sb <- data.frame(xmin = x0 + (0:3) * 1e5, xmax = x0 + (1:4) * 1e5, ymin = y0, ymax = y0 + h,
                 f = c("black", "white", "black", "white"))
## north arrow (orienteering style): ring with a half black / half white arrow and N above
ax <- unname(bb["xmin"]) + 0.07 * W; ay <- unname(bb["ymax"]) - 0.14 * H; s <- 0.06 * H
ring <- data.frame(x = ax + 0.5 * s * cos(seq(0, 2 * pi, length.out = 100)),
                   y = ay + 0.5 * s * sin(seq(0, 2 * pi, length.out = 100)))
arL <- data.frame(x = c(ax, ax - 0.32 * s, ax), y = c(ay + 0.62 * s, ay - 0.42 * s, ay - 0.18 * s))
arR <- data.frame(x = c(ax, ax + 0.32 * s, ax), y = c(ay + 0.62 * s, ay - 0.42 * s, ay - 0.18 * s))
deco <- list(
  geom_rect(data = sb, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = sb$f,
            colour = "black", linewidth = 0.25, inherit.aes = FALSE),
  annotate("text", x = x0 + c(0, 2e5, 4e5), y = y0 + h + 0.03 * H, label = c("0", "200", "400 km"), size = 3, family = FONT),
  geom_path(data = ring, aes(x, y), colour = "black", linewidth = 0.5, inherit.aes = FALSE),
  geom_polygon(data = arR, aes(x, y), fill = "white", colour = "black", linewidth = 0.4, inherit.aes = FALSE),
  geom_polygon(data = arL, aes(x, y), fill = "black", colour = "black", linewidth = 0.4, inherit.aes = FALSE),
  annotate("text", x = ax, y = ay + 0.9 * s, label = "N", size = 4.2, family = FONT),
  coord_sf(crs = 20137, datum = 4326, expand = TRUE),
  scale_x_continuous(breaks = seq(34, 48, 2)), scale_y_continuous(breaks = seq(4, 14, 2)),
  theme_bw(base_size = 11, base_family = FONT),
  theme(panel.border = element_rect(colour = "black", linewidth = 0.6),
        panel.grid = element_line(colour = "grey88", linewidth = 0.3), axis.title = element_blank(),
        axis.text = element_text(size = 8, colour = "grey30"),
        plot.title = element_text(size = 13, face = "bold"),
        legend.title = element_text(size = 11), legend.text = element_text(size = 9),
        legend.key.height = unit(1.4, "cm"), legend.key.width = unit(0.45, "cm")))

pA <- ggplot(m) + geom_sf(aes(fill = pi * 100), colour = "grey40", linewidth = 0.15) +
  scale_fill_distiller(palette = "Greens", direction = 1, limits = c(0, 40), breaks = seq(10, 40, 10),
                       labels = function(x) paste0(x, "%"), name = "Outbreak\nprobability (\u03c0)") +
  deco + labs(title = "A. Probability of an outbreak per month")
pB <- ggplot(m) + geom_sf(aes(fill = mu_plot), colour = "grey40", linewidth = 0.15) +
  scale_fill_distiller(palette = "Reds", direction = 1, limits = c(0, 550), breaks = seq(100, 500, 100),
                       na.value = "grey85", name = "Expected cases\nper outbreak\nmonth (\u03bc)") +
  deco + labs(title = "B. Expected outbreak size")

p <- (pA | pB) + plot_annotation(
  title = "Zone-level cholera outbreak probability and size in Ethiopia",
  subtitle = "Bayesian hurdle model (INLA, BYM2 spatial effect), 2015\u20132024",
  caption = paste0("Source: EPHI national surveillance data. Posterior means from the Baseline hurdle model: (A) probability of at least one case in a month; ",
                   "(B) expected cases in a month with an outbreak.\nGrey zones in panel B (n = ", n_none, ") reported no outbreak month during the study period. ",
                   "Administrative boundaries: OCHA Ethiopia (2021). Projection: Adindan / UTM zone 37N."),
  theme = theme(plot.title = element_text(size = 17, face = "bold", hjust = 0.5, family = FONT),
                plot.subtitle = element_text(size = 12, face = "italic", hjust = 0.5, colour = "grey30", family = FONT),
                plot.caption = element_text(size = 8.5, colour = "grey30", hjust = 0, family = FONT)))
ggsave(file.path(outdir, "Fig4_occurrence_size_maps_v2.png"),  p, width = 13, height = 6, dpi = 600, bg = "white")
ggsave(file.path(outdir, "Fig4_occurrence_size_maps_v2.tiff"), p, width = 13, height = 6, dpi = 600, compression = "lzw", bg = "white")
