## Figure 5 (profiles, communities, classes) and Figure 6 (central zones), updated results, style of the earlier figures
##   Rscript Fig5_Fig6_stage3_updated.R <zone table csv> <shapefile> <output folder>
suppressPackageStartupMessages({library(data.table); library(sf); library(ggplot2); library(patchwork)})
sf_use_s2(FALSE)
args <- commandArgs(TRUE)
infile   <- if (length(args) >= 1) args[1] else "Table_S_zone_profile_community_class.csv"
shp_path <- if (length(args) >= 2) args[2] else "eth_admbnda_adm2_csa_bofedb_2021.shp"
outdir   <- if (length(args) >= 3) args[3] else "."
FONT <- "sans"; N_LABEL <- 5; CENT_MIN <- 0.10

z <- fread(infile)
shp <- st_read(shp_path, quiet = TRUE) |> st_make_valid() |> st_transform(20137) |> st_simplify(dTolerance = 1000)
m <- merge(shp[, c("ADM2_PCODE", "geometry")], z, by.x = "ADM2_PCODE", by.y = "zone"); stopifnot(nrow(m) == 92)
m$uncertain <- as.logical(m$uncertain)

## profiles, communities, classes
m$prof_f <- factor(m$profile, 1:3, c("Higher burden, large peaks", "Low burden, small peaks", "Little or no recorded activity"))
low_comm <- c(1, 3)                                              # communities in the lower-activity class (Stage 3 output)
m$class_f <- factor(ifelse(m$community %in% low_comm, "Lower estimated activity", "Higher estimated peaks"),
                    c("Higher estimated peaks", "Lower estimated activity"))
## community colours shaded by class: blues = lower-activity class, oranges/reds = higher-peak class
comm_col <- c("1" = "#0072B2", "2" = "#E69F00", "3" = "#009E73", "4" = "#F0E442", "5" = "#CC79A7", "6" = "#56B4E9")   # Okabe-Ito, colour-blind safe
m$comm_f <- factor(m$community, 1:6)
## class x certainty (replaces hatching)
m$cc_f <- factor(paste0(m$class_f, ifelse(m$uncertain, " - uncertain", " - stable")),
  c("Higher estimated peaks - stable", "Higher estimated peaks - uncertain", "Lower estimated activity - stable", "Lower estimated activity - uncertain"))

## hatching for uncertain zones
bb <- st_bbox(shp); W <- unname(bb["xmax"] - bb["xmin"]); H <- unname(bb["ymax"] - bb["ymin"])
offs <- seq(-H, W, by = 18000)
hl <- lapply(offs, function(o) st_linestring(rbind(c(bb["xmin"] + o, bb["ymin"]), c(bb["xmin"] + o + H, bb["ymax"]))))
hatch <- st_intersection(st_sfc(hl, crs = 20137), st_union(m[m$uncertain, ]))

## north arrow (top right) and scale bar (bottom left)
x0 <- unname(bb["xmin"]) + 0.01 * W; y0 <- unname(bb["ymin"]) - 0.07 * H; h <- 0.016 * H   # below the country outline
sb <- data.frame(xmin = x0 + (0:2) * 1e5, xmax = x0 + (1:3) * 1e5, ymin = y0, ymax = y0 + h)
YLIM <- c(unname(bb["ymin"]) - 0.09 * H, unname(bb["ymax"]) + 0.01 * H)
ax <- unname(bb["xmax"]) - 0.06 * W; ay <- unname(bb["ymax"]) - 0.08 * H; s <- 0.06 * H
ring <- data.frame(x = ax + 0.5 * s * cos(seq(0, 2 * pi, length.out = 100)), y = ay + 0.5 * s * sin(seq(0, 2 * pi, length.out = 100)))
arL <- data.frame(x = c(ax, ax - 0.32 * s, ax), y = c(ay + 0.62 * s, ay - 0.42 * s, ay - 0.18 * s))
arR <- data.frame(x = c(ax, ax + 0.32 * s, ax), y = c(ay + 0.62 * s, ay - 0.42 * s, ay - 0.18 * s))
deco <- function(base = 11, legend_pos = "bottom") list(
  geom_rect(data = sb, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = c("black", "white", "black"),
            colour = "black", linewidth = 0.25, inherit.aes = FALSE),
  annotate("text", x = x0 + 3e5 + 0.012 * W, y = y0 + h / 2, label = "300 km", hjust = 0, size = base / 3.6, family = FONT),
  geom_path(data = ring, aes(x, y), colour = "black", linewidth = 0.5, inherit.aes = FALSE),
  geom_polygon(data = arR, aes(x, y), fill = "white", colour = "black", linewidth = 0.4, inherit.aes = FALSE),
  geom_polygon(data = arL, aes(x, y), fill = "black", colour = "black", linewidth = 0.4, inherit.aes = FALSE),
  annotate("text", x = ax, y = ay + 0.9 * s, label = "N", size = base / 2.8, family = FONT),
  coord_sf(crs = 20137, datum = 4326, ylim = YLIM, expand = FALSE),
  scale_x_continuous(breaks = seq(34, 48, 2)), scale_y_continuous(breaks = seq(4, 14, 2)),
  theme_bw(base_size = base, base_family = FONT),
  theme(panel.border = element_rect(colour = "black", linewidth = 0.6), panel.grid = element_line(colour = "grey88", linewidth = 0.3),
        axis.title = element_blank(), axis.text = element_text(size = base * 0.68, colour = "grey30"),
        plot.title = element_text(size = base * 1.15, face = "bold", hjust = 0.5), legend.position = legend_pos,
        legend.title = element_text(size = base * 1.05, face = "bold"), legend.text = element_text(size = base * 0.98),
        legend.key.size = unit(0.5, "cm"), plot.margin = margin(4, 10, 4, 10)))
bord <- list(colour = "white", linewidth = 0.25)

## community labels drawn on panel B (largest contiguous piece of each community)
comm_lab <- do.call(rbind, lapply(levels(m$comm_f), function(k) {
  g <- st_union(st_geometry(m[m$community == as.integer(k), ]))
  parts <- suppressWarnings(st_cast(st_make_valid(g), "POLYGON"))
  big <- parts[which.max(st_area(parts))]
  p <- st_coordinates(st_point_on_surface(big))
  data.frame(comm = k, X = p[1, 1], Y = p[1, 2]) }))

## ---------------- Figure 5 ----------------
nlab <- function(f) paste0(levels(f), " (", as.vector(table(f)), ")")
pA <- ggplot(m) + geom_sf(aes(fill = prof_f), colour = "white", linewidth = 0.25) +
  scale_fill_manual(values = c("#99000D", "#E6550D", "#FEE391"), labels = nlab(m$prof_f), name = "Latent profile") +
  guides(fill = guide_legend(ncol = 1)) + deco() + labs(title = "A. Latent profiles (3)")
pB <- ggplot(m) + geom_sf(aes(fill = comm_f), colour = "white", linewidth = 0.25) +
  lapply(c(-1, 1), function(dx) lapply(c(-1, 1), function(dy)
    geom_text(data = comm_lab, aes(X + dx * 8000, Y + dy * 8000, label = comm), colour = "white", size = 6, fontface = "bold", inherit.aes = FALSE))) +
  geom_text(data = comm_lab, aes(X, Y, label = comm), colour = "grey10", size = 6, fontface = "bold", inherit.aes = FALSE) +
  scale_fill_manual(values = comm_col[levels(m$comm_f)], labels = paste0(levels(m$comm_f), " (", as.vector(table(m$comm_f)), ")"),
                    name = "Spatial risk community") +
  guides(fill = guide_legend(ncol = 3)) + deco() + labs(title = "B. Spatial risk communities (6)")
pC <- ggplot(m) + geom_sf(aes(fill = cc_f), colour = "white", linewidth = 0.25) +
  scale_fill_manual(values = c("#D55E00", "#F5C29C", "#0072B2", "#A6CDEB"), drop = FALSE, name = "Class",
                    labels = paste0(levels(m$cc_f), " (", as.vector(table(m$cc_f)), ")")) +
  guides(fill = guide_legend(ncol = 1)) + deco() + labs(title = "C. Classes and stability")
n_unc <- sum(m$uncertain)
p5 <- (pA | pB | pC) + plot_annotation(
  title = "Latent profiles, spatial risk communities and classes of cholera activity in Ethiopia (2015-2024)",
  subtitle = "Population-adjusted hidden Markov model indicators (burden and peak reached)",
  caption = paste0("(A) Latent profile analysis of burden and peak reached. (B) Louvain communities of the zone adjacency network, edges weighted by similarity in burden and peak reached. ",
                   "\n(C) Colour denotes class (orange = higher estimated peaks, blue = lower estimated activity) and shading denotes stability: dark = stable zones (n = ", 92 - n_unc,
                   ") that kept their profile, community and class in at least 70% of 500 posterior draws; pale = uncertain zones (n = ", n_unc, ").\nAll stable zones lie in ", paste0("Community 3 (", sum(!m$uncertain & m$community == 3), "), Community 5 (", sum(!m$uncertain & m$community == 5), ") and Community 6 (", sum(!m$uncertain & m$community == 6), ")"), ".",
                   "\nCommunity numbers are shown on the map in (B); numbers in brackets in each legend give the number of zones. Classes in (C) are not the profiles of (A): class assignment uses the community means of burden and peak reached.\nAdjacency does not represent transmission. Boundaries: OCHA Ethiopia (2021); projection Adindan / UTM 37N."),
  theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5, family = FONT),
                plot.subtitle = element_text(size = 11.5, face = "italic", hjust = 0.5, colour = "grey30", family = FONT),
                plot.caption = element_text(size = 8.5, colour = "grey30", hjust = 0, family = FONT)))
ggsave(file.path(outdir, "Fig5_profiles_communities_classes.png"),  p5, width = 16, height = 7.5, dpi = 600, bg = "white")
ggsave(file.path(outdir, "Fig5_profiles_communities_classes.tiff"), p5, width = 16, height = 7.5, dpi = 600, compression = "lzw", bg = "white")

## ---------------- Figure 6 ----------------
pts <- m[m$centrality >= CENT_MIN, ]
pts_xy <- cbind(st_drop_geometry(pts), st_coordinates(st_point_on_surface(st_geometry(pts))))
lab_xy <- pts_xy[order(-pts_xy$centrality), ][1:N_LABEL, ]
## label offsets (metres) so names do not cover each other; edit if needed
lab_xy$dx <- 0; lab_xy$dy <- 0
off <- list("Bale" = c(150000, -10000), "West Arsi" = c(-160000, 20000), "Guji" = c(-140000, -40000),
            "East Bale" = c(170000, 10000), "Arsi" = c(60000, 75000))
for (k in names(off)) { i <- lab_xy$zone_name == k; lab_xy$dx[i] <- off[[k]][1]; lab_xy$dy[i] <- off[[k]][2] }
lab_xy$lx <- lab_xy$X + lab_xy$dx; lab_xy$ly <- lab_xy$Y + lab_xy$dy
p6 <- ggplot(m) + geom_sf(aes(fill = comm_f), colour = "white", linewidth = 0.25) +
  scale_fill_manual(values = comm_col[levels(m$comm_f)], name = "Spatial risk\ncommunity",
                    labels = paste0("Community ", levels(m$comm_f), " (", as.vector(table(m$comm_f)), ")")) +
  geom_segment(data = lab_xy, aes(x = X, y = Y, xend = lx, yend = ly), colour = "grey15", linewidth = 0.35) +
  geom_point(data = pts_xy, aes(X, Y, size = centrality, colour = centrality), alpha = 0.95) +
  geom_point(data = pts_xy, aes(X, Y, size = centrality), shape = 21, colour = "white", fill = NA, stroke = 0.5, show.legend = FALSE) +
  scale_size_continuous(range = c(2, 9), limits = c(CENT_MIN, 1), breaks = c(0.2, 0.4, 0.6, 0.8, 1.0), name = "Eigenvector\ncentrality") +
  scale_colour_gradient(low = "#D9D9D9", high = "#000000", limits = c(CENT_MIN, 1), breaks = c(0.2, 0.4, 0.6, 0.8, 1.0), name = "Eigenvector\ncentrality") +
  geom_label(data = lab_xy, aes(lx, ly, label = paste0(zone_name, " (", sprintf("%.2f", centrality), ")")), size = 3.3, family = FONT,
             label.padding = unit(0.15, "lines"), label.size = 0.2, fill = alpha("white", 0.9)) +
  guides(size = guide_legend(order = 1), colour = guide_legend(order = 1), fill = guide_legend(order = 2)) +
  deco(base = 12, legend_pos = "right") +
  labs(title = "Central zones in the zone adjacency network",
       subtitle = paste0("Six spatial risk communities; bubbles (centrality >= ", sprintf("%.2f", CENT_MIN),
                         ") sized and shaded by eigenvector centrality; the ", N_LABEL, " most central zones are labelled"),
       caption = "Centrality describes position in the adjacency network weighted by similarity in burden and peak reached; it does not show that a zone drives transmission.\nBoundaries: OCHA Ethiopia (2021); projection Adindan / UTM 37N.") +
  theme(plot.subtitle = element_text(size = 10, colour = "grey30", hjust = 0.5), plot.title = element_text(size = 15),
        plot.caption = element_text(size = 8.5, colour = "grey30", hjust = 0))
ggsave(file.path(outdir, "Fig6_central_zones.png"),  p6, width = 11, height = 8.5, dpi = 600, bg = "white")
ggsave(file.path(outdir, "Fig6_central_zones.tiff"), p6, width = 11, height = 8.5, dpi = 600, compression = "lzw", bg = "white")
cat("Saved Fig 5 and Fig 6\n")
