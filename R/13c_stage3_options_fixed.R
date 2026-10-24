## Stability with the profiles held fixed: profiles are fitted once on the posterior means; each draw (or refit) is then
## assigned to those profiles (predict). This never fails on degenerate draws and gives each zone's membership probability.
suppressPackageStartupMessages({ library(data.table); library(mclust) })
set.seed(2026)
proj <- path.expand("~/Documents/spatial_anlysis/update_HMM"); OUT <- file.path(proj, "test3_options")
D <- readRDS(file.path(proj, "final_pop/improve_hmm/Draws_updated_indicators_tau2_2_hi170_pop_rmin0.5_hm_FINAL.rds")); I <- D$ind; Z <- length(D$zone_order); nd <- nrow(I$burden)
R <- readRDS(file.path(proj, "joint_20draws_long_pop/Draws_joint20_long_indicators.rds"))$ind
tab_of <- function(SRC, d) data.table(burden = SRC$burden[d, ], entry = SRC$entry[d, ], persistence_active = SRC$persistence_active[d, ], peak_reached = SRC$peak_reached[d, ])
mean_tab <- data.table(burden = colMeans(I$burden), entry = colMeans(I$entry), persistence_active = colMeans(I$persistence_active), peak_reached = colMeans(I$peak_reached))
fitprep <- function(M, v) { pr <- lapply(v, function(j) { x <- M[[j]]; sk <- mean((x - mean(x))^3) / sd(x)^3; lg <- is.finite(sk) && sk > 1; off <- 0.5 * min(x[x > 0])
  y <- if (lg) log(x + off) else x; list(lg = lg, off = off, m = mean(y), s = sd(y)) }); names(pr) <- v; pr }
apply_prep <- function(M, pr) sapply(names(pr), function(j) { p <- pr[[j]]; x <- M[[j]]; y <- if (p$lg) log(pmax(x, 0) + p$off) else x; (y - p$m) / p$s })
ds <- unique(round(seq(1, nd, length.out = 500)))
sets <- list("burden+peak" = c("burden", "peak_reached"), "burden+entry+peak" = c("burden", "entry", "peak_reached"), "all four" = c("burden", "entry", "persistence_active", "peak_reached"))
res <- rbindlist(lapply(names(sets), function(sn) rbindlist(lapply(2:5, function(k) { v <- sets[[sn]]; pr <- fitprep(mean_tab, v)
  m0 <- tryCatch(Mclust(apply_prep(mean_tab, pr), G = k, modelNames = "VVI", verbose = FALSE), error = function(e) NULL)
  if (is.null(m0)) return(data.table(indicators = sn, K = k))
  main <- m0$classification; p <- pmax(m0$z, 1e-12); en <- 1 + sum(p * log(p)) / (nrow(p) * log(ncol(p)))
  CL <- sapply(ds, function(d) predict(m0, newdata = apply_prep(tab_of(I, d), pr))$classification)            # Z x draws
  CR <- sapply(seq_len(nrow(R$burden)), function(j) predict(m0, newdata = apply_prep(tab_of(R, j), pr))$classification)
  a <- apply(CL, 2, function(x) adjustedRandIndex(main, x)); ar <- apply(CR, 2, function(x) adjustedRandIndex(main, x))
  keep <- rowMeans(CL == main)                                                                                   # each zone: share of draws in its main profile
  data.table(indicators = sn, K = k, BIC = round(m0$bic, 1), entropy = round(en, 3), smallest = min(table(main)),
             ARI_draws = round(mean(a), 3), ARI_draws_lo = round(quantile(a, .025), 3), ARI_refits = round(mean(ar), 3),
             zones_70 = sum(keep >= 0.7), zones_90 = sum(keep >= 0.9), median_keep = round(median(keep), 2)) }))))
cat("Profiles fitted once on the posterior means; 500 draws and 20 refits assigned to them\n(zones_70 / zones_90 = zones whose draws fall in their main profile at least 70% / 90% of the time)\n\n")
print(res, row.names = FALSE); fwrite(res, file.path(OUT, "Table_options_LPA_fixed_profiles.csv"))
