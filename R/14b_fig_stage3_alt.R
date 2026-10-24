## Stage 3 figure: (A) latent profiles, (B) spatial risk communities with central zones, (C) classes, uncertain zones hatched
## Same style as Fig 4 (north arrow, scale bar, graticule, colour-blind safe colours, 600 dpi)
suppressPackageStartupMessages({library(data.table); library(sf); library(ggplot2); library(patchwork)})
sf_use_s2(FALSE)
args <- commandArgs(TRUE)
infile   <- if (length(args) >= 1) args[1] else "Table_S_zone_profile_community_class.csv"
shp_path <- if (length(args) >= 2) args[2] else "eth_admbnda_adm2_csa_bofedb_2021.shp"
outdir   <- if (length(args) >= 3) args[3] else "."
FONT <- "sans"; N_LABEL <- 5                                   # number of central zones labelled in panel B

z <- fread(infile)
if (!"profile" %in% names(z)) stop("column 'profile' missing")
shp <- st_read(shp_path, quiet = TRUE) |> st_make_valid() |> st_transform(20137) |> st_simplify(dTolerance = 1000)
m <- merge(shp[, c("ADM2_PCODE", "geometry")], z, by.x = "ADM2_PCODE", by.y = "zone")
stopifnot(nrow(m) == 92)
m$uncertain <- as.logical(m$uncertain)

## labels and classes (classes: communities grouped by their mean indicators, as in the Stage 3 script)
m$prof_f <- factor(m$profile, 1:3, c("Higher estimated activity", "Intermediate", "Lowest estimated activity"))
m$comm_f <- factor(m$community)
low_comm <- c(1, 3)                                            # communities whose class is "Lower-activity" in the Stage 3 output
m$class_f <- factor(ifelse(m$community %in% low_comm, "Lower estimated activity", "Higher estimated peaks"),
                    c("Higher estimated peaks", "Lower estimated activity"))

## hatching for uncertain zones (diagonal lines clipped to the zones)
bb <- st_bbox(shp); W <- unname(bb["xmax"] - bb["xmin"]); H <- unname(bb["ymax"] - bb["ymin"])
step <- 18000; offs <- seq(-H, W, by = step)
hl <- lapply(offs, function(o) st_linestring(rbind(c(bb["xmin"] + o, bb["ymin"]), c(bb["xmin"] + o + H, bb["ymax"]))))
hatch <- st_intersection(st_sfc(hl, crs = 20137), st_union(m[m$uncertain, ]))

## central zones (highest eigenvector centrality)
lab <- m[order(-m$centrality), ][1:N_LABEL, ]
lab_xy <- cbind(st_drop_geometry(lab), st_coordinates(st_point_on_surface(lab)))

## north arrow and scale bar
x0 <- unname(bb["xmax"]) - 4e5 - 0.03 * W; y0 <- unname(bb["ymin"]) + 0.02 * H; h <- 0.016 * H
sb <- data.frame(xmin = x0 + (0:3) * 1e5, xmax = x0 + (1:4) * 1e5, ymin = y0, ymax = y0 + h)
ax <- unname(bb["xmin"]) + 0.07 * W; ay <- unname(bb["ymax"]) - 0.14 * H; s <- 0.06 * H
ring <- data.frame(x = ax + 0.5 * s * cos(seq(0, 2 * pi, length.out = 100)), y = ay + 0.5 * s * sin(seq(0, 2 * pi, length.out = 100)))
arL <- data.frame(x = c(ax, ax - 0.32 * s, ax), y = c(ay + 0.62 * s, ay - 0.42 * s, ay - 0.18 * s))
arR <- data.frame(x = c(ax, ax + 0.32 * s, ax), y = c(ay + 0.62 * s, ay - 0.42 * s, ay - 0.18 * s))
deco <- list(
  geom_rect(data = sb, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = c("black", "white", "black", "white"),
            colour = "black", linewidth = 0.25, inherit.aes = FALSE),
  annotate("text", x = x0 + c(0, 2e5, 4e5), y = y0 + h + 0.03 * H, label = c("0", "200", "400 km"), size = 2.8, family = FONT),
  geom_path(data = ring, aes(x, y), colour = "black", linewidth = 0.5, inherit.aes = FALSE),
  geom_polygon(data = arR, aes(x, y), fill = "white", colour = "black", linewidth = 0.4, inherit.aes = FALSE),
  geom_polygon(data = arL, aes(x, y), fill = "black", colour = "black", linewidth = 0.4, inherit.aes = FALSE),
  annotate("text", x = ax, y = ay + 0.9 * s, label = "N", size = 4, family = FONT),
  coord_sf(crs = 20137, datum = 4326, expand = TRUE),
  scale_x_continuous(breaks = seq(34, 48, 4)), scale_y_continuous(breaks = seq(4, 14, 2)),
  theme_bw(base_size = 11, base_family = FONT),
  theme(panel.border = element_rect(colour = "black", linewidth = 0.6), panel.grid = element_line(colour = "grey88", linewidth = 0.3),
        axis.title = element_blank(), axis.text = element_text(size = 7.5, colour = "grey30"),
        plot.title = element_text(size = 12, face = "bold"), legend.position = "bottom", legend.direction = "vertical",
        legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 9)))

pA <- ggplot(m) + geom_sf(aes(fill = prof_f), colour = "grey40", linewidth = 0.15) +
  scale_fill_manual(values = c("#CB181D", "#FC9272", "#FEE5D9"), name = "Latent profile", drop = FALSE,
                    labels = paste0(levels(m$prof_f), " (n = ", as.vector(table(m$prof_f)), ")")) +
  deco + labs(title = "A. Latent profiles")
pB <- ggplot(m) + geom_sf(aes(fill = comm_f), colour = "grey40", linewidth = 0.15) +
  geom_point(data = lab_xy, aes(X, Y), shape = 21, fill = "white", colour = "black", size = 1.8, stroke = 0.6) +
  geom_label(data = lab_xy, aes(X, Y, label = zone_name), size = 2.6, family = FONT, nudge_y = 45000,
             label.padding = unit(0.1, "lines"), label.size = 0.15, fill = alpha("white", 0.85)) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#CC79A7"), name = "Spatial risk community",
                    labels = paste0("Community ", levels(m$comm_f), " (n = ", as.vector(table(m$comm_f)), ")")) +
  guides(fill = guide_legend(ncol = 2)) + deco + labs(title = "B. Spatial risk communities")
pC <- ggplot(m) + geom_sf(aes(fill = class_f), colour = "grey40", linewidth = 0.15) +
  geom_sf(data = hatch, colour = "grey15", linewidth = 0.25) +
  scale_fill_manual(values = c("#D55E00", "#9ECAE1"), name = "Class",
                    labels = paste0(levels(m$class_f), " (n = ", as.vector(table(m$class_f)), ")")) +
  deco + labs(title = "C. Classes (hatched: uncertain)")

n_unc <- sum(m$uncertain)
p <- (pA | pB | pC) + plot_annotation(
  title = "Latent profiles, spatial risk communities and classes of cholera activity in Ethiopia",
  subtitle = "Population-adjusted hidden Markov model indicators, 2015–2024",
  caption = paste0("Profiles: latent profile analysis of burden and peak reached. Communities: Louvain partition of the zone adjacency network, ",
                   "edges weighted by similarity in burden and peak reached;\nlabelled points mark the ", N_LABEL,
                   " zones with the highest eigenvector centrality. Hatched zones (n = ", n_unc, ") kept their profile, community and class in fewer than 70% ",
                   "of 500 posterior draws. Adjacency does not represent transmission.\nAdministrative boundaries: OCHA Ethiopia (2021). Projection: Adindan / UTM zone 37N."),
  theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5, family = FONT),
                plot.subtitle = element_text(size = 11.5, face = "italic", hjust = 0.5, colour = "grey30", family = FONT),
                plot.caption = element_text(size = 8.5, colour = "grey30", hjust = 0, family = FONT)))
ggsave(file.path(outdir, "Fig_Stage3_profile_community_class.png"),  p, width = 16, height = 7.5, dpi = 600, bg = "white")
ggsave(file.path(outdir, "Fig_Stage3_profile_community_class.tiff"), p, width = 16, height = 7.5, dpi = 600, compression = "lzw", bg = "white")
cat("Saved\n")
