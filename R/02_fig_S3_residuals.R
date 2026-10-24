## Fig. S3: zone-level residuals of the Baseline hurdle model and residual Moran's I
suppressPackageStartupMessages({library(data.table); library(sf); library(spdep); library(ggplot2); library(patchwork)})
S   <- readRDS(file.path(PROJ))   # laptop: file.path(PROJ)
shp_path <- "F:/Cholera/shape_22/Adm_2/eth_admbnda_adm2_csa_bofedb_2021.shp"
outdir   <- getOption("cholera.proj", ".")
dat <- as.data.table(S$data); fo <- S$fits$occ; fi <- S$fits$int; NN <- nrow(dat)
A <- as.matrix(S$graph_Q); diag(A) <- 0; A <- (A != 0) * 1
nbq <- mat2listw(A, style = "B", zero.policy = TRUE)$neighbours
lw  <- nb2listw(nbq, style = "W", zero.policy = TRUE)
# occurrence: zone-level standardized residual, sum(y - p) / sqrt(sum p(1-p))
p  <- fo$summary.fitted.values$mean[1:NN]
ro <- tapply(dat$Y_bin - p, dat$ID_area, sum) / sqrt(tapply(p * (1 - p), dat$ID_area, sum))
# intensity: Pearson residual under the zero-truncated NB, zone sum / sqrt(n), zones with >= 3 outbreak months
lp <- fi$summary.linear.predictor[1:NN, ]; mu <- exp(lp$mean + 0.5 * lp$sd^2)
sz <- fi$summary.hyperpar[grep("^size", rownames(fi$summary.hyperpar))[1], "mean"]
q  <- -expm1(-sz * log1p(mu / sz)); m <- mu / q; v <- (mu * (1 + mu / sz) + mu^2) / q - m^2
r  <- ifelse(is.na(dat$O_int), NA_real_, (dat$O_int - m) / sqrt(v))
ri <- tapply(r, dat$ID_area, function(x) { x <- x[!is.na(x)]; if (length(x) >= 3) sum(x) / sqrt(length(x)) else NA_real_ })
keep <- !is.na(ri)
mo <- moran.test(as.numeric(ro), lw, zero.policy = TRUE)
mi <- moran.test(as.numeric(ri[keep]), subset.listw(lw, keep, zero.policy = TRUE), zero.policy = TRUE)
cat(sprintf("Occurrence Moran's I = %.3f, p = %.3f | Intensity Moran's I = %.3f, p = %.3f (%d zones)\n",
            mo$estimate[1], mo$p.value, mi$estimate[1], mi$p.value, sum(keep)))
# maps
shp <- st_read(shp_path, quiet = TRUE); shp <- st_transform(st_make_valid(shp), 20137)
shp$key <- toupper(gsub("[^A-Za-z0-9]", "", shp$ADM2_PCODE)); shp <- st_simplify(shp, dTolerance = 1500)
zo <- unique(dat[, .(ID_area, admin2_clean)])[order(ID_area)]
zo[, `:=`(r_occ = as.numeric(ro), r_int = as.numeric(ri))]
g <- merge(shp[, c("key", "geometry")], zo, by.x = "key", by.y = "admin2_clean")
lim <- max(abs(c(zo$r_occ, zo$r_int)), na.rm = TRUE)
mapf <- function(var, ttl, I) ggplot(g) + geom_sf(aes(fill = .data[[var]]), colour = "white", linewidth = 0.08) +
  scale_fill_gradient2(low = "#2a78d6", mid = "#f1f0ec", high = "#d95926", midpoint = 0, limits = c(-lim, lim),
                       na.value = "grey85", name = "Zone residual") +
  labs(title = ttl, subtitle = I) + theme_void(base_size = 10, base_family = "Liberation Serif") +
  theme(panel.border = element_rect(fill = NA, colour = "grey40", linewidth = 0.3), plot.title = element_text(size = 10.5),
        plot.subtitle = element_text(size = 9, colour = "grey25"))
# Moran scatterplots (residual vs spatial lag)
scat <- function(x, w, ttl) { d <- data.frame(x = x, lag = lag.listw(w, x, zero.policy = TRUE))
  ggplot(d, aes(x, lag)) + geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) + geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_point(colour = "#2a78d6", size = 1.6, alpha = 0.8) + geom_smooth(method = "lm", se = FALSE, colour = "#123a6b", linewidth = 0.6, formula = y ~ x) +
    labs(title = ttl, x = "Zone residual", y = "Mean residual of neighbours") + theme_bw(base_size = 10, base_family = "Liberation Serif") +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(size = 10.5)) }
lab_o <- sprintf("Moran's I = %.2f, p = %.2f", mo$estimate[1], mo$p.value)
lab_i <- sprintf("Moran's I = %.2f, p = %.2f (%d zones with at least 3 outbreak months)", mi$estimate[1], mi$p.value, sum(keep))
p <- (mapf("r_occ", "A. Occurrence residuals", lab_o) | mapf("r_int", "B. Outbreak-size residuals", lab_i)) /
     (scat(as.numeric(ro), lw, "C. Moran scatterplot, occurrence") | scat(as.numeric(ri[keep]), subset.listw(lw, keep, zero.policy = TRUE), "D. Moran scatterplot, outbreak size")) +
     plot_layout(guides = "collect", heights = c(1.15, 1)) & theme(legend.position = "right")
ggsave(file.path(outdir, "Fig_S3_residual_spatial.png"),  p, width = 10, height = 8, dpi = 300, bg = "white")
ggsave(file.path(outdir, "Fig_S3_residual_spatial.tiff"), p, width = 10, height = 8, dpi = 300, compression = "lzw", bg = "white")
