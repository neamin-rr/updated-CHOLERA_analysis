## Figure 4: zone-level outbreak probability and expected outbreak size, Baseline hurdle model (ColorBrewer, colour-blind safe)
suppressPackageStartupMessages({library(data.table); library(sf); library(ggplot2); library(patchwork); library(RColorBrewer)})
sf_use_s2(FALSE)
infile   <- file.path(PROJ)
shp_path <- "F:/Cholera/shape_22/Adm_2/eth_admbnda_adm2_csa_bofedb_2021.shp"
outdir   <- file.path(PROJ)
FONT <- "serif"   # or "Times New Roman"
z <- fread(infile)
shp <- st_read(shp_path, quiet = TRUE) |> st_make_valid() |> st_transform(20137)
shp$key <- toupper(gsub("[^A-Za-z0-9]", "", shp$ADM2_PCODE)); shp <- st_simplify(shp, dTolerance = 1200)
m <- merge(shp[, c("key", "ADM1_EN", "geometry")], z, by.x = "key", by.y = "admin2_clean")
m$no_outbreak <- m$N_outbreak_months == 0
adm1 <- aggregate(shp["geometry"], by = list(R = shp$ADM1_EN), FUN = function(x) x[1])
cent <- suppressWarnings(st_point_on_surface(m[m$no_outbreak, ]))
bb <- st_bbox(shp); W <- bb["xmax"] - bb["xmin"]; H <- bb["ymax"] - bb["ymin"]
x0 <- unname(bb["xmin"] + 0.03 * W); y0 <- unname(bb["ymin"] + 0.02 * H); h <- 0.018 * H
sb  <- data.frame(xmin = x0 + c(0, 1e5, 2e5, 3e5), xmax = x0 + c(1e5, 2e5, 3e5, 4e5), ymin = y0, ymax = y0 + h, f = c("black", "white", "black", "white"))
sbt <- data.frame(x = x0 + c(0, 2e5, 4e5), y = y0 + 2.4 * h, lab = c("0", "200", "400 km"))
ax <- unname(bb["xmin"] + 0.08 * W); ay <- unname(bb["ymax"] - 0.16 * H); s <- 0.06 * H
na <- data.frame(x = c(ax - 0.4 * s, ax, ax + 0.4 * s, ax), y = c(ay, ay + 1.6 * s, ay, ay + 0.45 * s))
deco <- list(
  geom_sf(data = adm1, fill = NA, colour = "grey25", linewidth = 0.25),
  geom_sf(data = cent, shape = 1, size = 1.4, stroke = 0.5, colour = "black"),
  geom_rect(data = sb, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = sb$f, colour = "black", linewidth = 0.2, inherit.aes = FALSE),
  geom_text(data = sbt, aes(x, y, label = lab), size = 2.4, family = FONT, inherit.aes = FALSE),
  geom_polygon(data = na, aes(x, y), fill = "black", inherit.aes = FALSE),
  annotate("text", x = ax, y = ay + 2.1 * s, label = "N", fontface = "bold", size = 3.4, family = FONT),
  coord_sf(crs = 20137, datum = 4326, expand = FALSE),
  theme_bw(base_size = 10, base_family = FONT),
  theme(panel.grid = element_line(colour = "grey90", linewidth = 0.25), axis.title = element_blank(),
        axis.text = element_text(size = 7, colour = "grey30"), legend.position = "right",
        plot.title = element_text(size = 11, face = "bold"), legend.title = element_text(size = 9)))
pA <- ggplot(m) + geom_sf(aes(fill = pi * 100), colour = "white", linewidth = 0.1) +
  scale_fill_distiller(palette = "YlGnBu", direction = 1, limits = c(0, 40), breaks = seq(0, 40, 10),
                       name = "Outbreak\nprobability\nper month (%)") + deco + labs(title = "A. Probability of an outbreak")
pB <- ggplot(m) + geom_sf(aes(fill = mu), colour = "white", linewidth = 0.1) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, trans = "log10", limits = c(5, 600),
                       breaks = c(5, 10, 25, 50, 100, 250, 500), oob = scales::squish,
                       name = "Expected cases\nper outbreak\nmonth (log scale)") + deco + labs(title = "B. Expected outbreak size")
p <- (pA | pB) + plot_annotation(caption = paste0(
  "Posterior means from the Baseline hurdle model (population-adjusted), 2015-2024. Open circles: zones with no reported outbreak month; values are model-based ",
  "(borrowed from neighbouring zones and other months).\nDark lines: regional boundaries. Boundaries: OCHA admin-2 (2021). Projection: Adindan / UTM zone 37N."),
  theme = theme(plot.caption = element_text(size = 7.5, colour = "grey30", hjust = 0, family = FONT)))
ggsave(file.path(outdir, "Fig4_occurrence_size_maps.png"),  p, width = 11, height = 4.7, dpi = 600, bg = "white")
ggsave(file.path(outdir, "Fig4_occurrence_size_maps.tiff"), p, width = 11, height = 4.7, dpi = 600, compression = "lzw", bg = "white")
