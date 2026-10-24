## Compare Stage 3 options by stability: indicator sets x K (LPA), network weights (Louvain), and a consensus grouping.
## Run: cd ~/Documents/spatial_anlysis/update_HMM ; nice -n 19 Rscript Stage3_stability_options.R > test3_options/log.txt 2>&1
suppressPackageStartupMessages({ library(data.table); library(mclust); library(igraph) })
set.seed(2026)
proj <- getOption("cholera.proj", path.expand("~/Documents/spatial_anlysis/update_HMM"))
DR <- getOption("cholera.draws", file.path(proj, "final_pop/improve_hmm/Draws_updated_indicators_tau2_2_hi170_pop_rmin0.5_hm_FINAL.rds"))
RF <- getOption("cholera.refits", file.path(proj, "joint_20draws_long_pop/Draws_joint20_long_indicators.rds"))
OUT <- file.path(proj, "test3_options"); dir.create(OUT, showWarnings = FALSE)
ND <- getOption("cholera.ndraws", 200)                                   # draws per option (200 is enough to compare options)
D <- readRDS(DR); I <- D$ind; zo <- D$zone_order; Z <- length(zo); nd <- nrow(I$burden)
R <- if (file.exists(RF)) readRDS(RF)$ind else NULL
S1 <- readRDS(file.path(proj, "INLA_Rebuilt_From_Zero_OCCPOP.rds"))
Qm <- as.matrix(S1$graph_Q); diag(Qm) <- 0; A <- ((Qm != 0) + t(Qm != 0) > 0) * 1; g0 <- graph_from_adjacency_matrix(A, mode = "undirected"); el <- as_edgelist(g0, names = FALSE)
tab_of <- function(SRC, d) data.table(burden = SRC$burden[d, ], entry = SRC$entry[d, ], persistence_active = SRC$persistence_active[d, ], peak_reached = SRC$peak_reached[d, ])
mean_tab <- data.table(burden = colMeans(I$burden), entry = colMeans(I$entry), persistence_active = colMeans(I$persistence_active), peak_reached = colMeans(I$peak_reached))
prep <- function(M, v) { X <- as.matrix(M[, ..v]); for (j in v) { x <- X[, j]; sk <- mean((x - mean(x))^3) / sd(x)^3; if (is.finite(sk) && sk > 1) X[, j] <- log(x + 0.5 * min(x[x > 0])) }; scale(X) }
lpa <- function(M, v, k) { m <- tryCatch(Mclust(prep(M, v), G = k, modelNames = "VVI", verbose = FALSE), error = function(e) NULL); if (is.null(m)) rep(NA_integer_, Z) else m$classification }
mm <- function(x) (x - min(x)) / (max(x) - min(x) + 1e-12)
louv <- function(M, w) { a <- mm(M[[w[1]]]); b <- mm(M[[w[2]]]); g <- g0; E(g)$weight <- 0.5 * ((a[el[, 1]] + a[el[, 2]]) / 2 + (b[el[, 1]] + b[el[, 2]]) / 2) + 1e-3; as.integer(membership(cluster_louvain(g))) }
ari <- function(main, L) sapply(L, function(x) if (anyNA(x)) NA_real_ else adjustedRandIndex(main, x))
ds <- unique(round(seq(1, nd, length.out = min(ND, nd))))
cat(sprintf("%d draws, %d zones; %d draws used per option; refits: %s\n", nd, Z, length(ds), if (is.null(R)) "not found" else nrow(R$burden)))
cat("\nBetween-zone / within-zone SD ratio (higher = indicator separates zones better than its own uncertainty):\n")
print(round(sapply(I, function(M) sd(colMeans(M)) / mean(apply(M, 2, sd))), 2))
## ---- A. LPA: indicator sets x K
sets <- list("burden+peak" = c("burden", "peak_reached"), "burden+entry+peak" = c("burden", "entry", "peak_reached"), "all four" = c("burden", "entry", "persistence_active", "peak_reached"))
resA <- rbindlist(lapply(names(sets), function(sn) rbindlist(lapply(2:5, function(k) { v <- sets[[sn]]
  m0 <- tryCatch(Mclust(prep(mean_tab, v), G = k, modelNames = "VVI", verbose = FALSE), error = function(e) NULL)
  if (is.null(m0)) return(data.table(indicators = sn, K = k, BIC = NA, entropy = NA, smallest = NA, ARI_draws = NA, ARI_draws_lo = NA, ARI_refits = NA))
  main <- m0$classification; p <- pmax(m0$z, 1e-12); en <- 1 + sum(p * log(p)) / (nrow(p) * log(ncol(p)))
  a <- ari(main, lapply(ds, function(d) lpa(tab_of(I, d), v, k)))
  ar <- if (is.null(R)) NA else mean(ari(main, lapply(seq_len(nrow(R$burden)), function(j) lpa(tab_of(R, j), v, k))), na.rm = TRUE)
  data.table(indicators = sn, K = k, BIC = round(m0$bic, 1), entropy = round(en, 3), smallest = min(table(main)), ARI_draws = round(mean(a, na.rm = TRUE), 3),
             ARI_draws_lo = round(quantile(a, .025, na.rm = TRUE), 3), ARI_refits = round(ar, 3)) }))))
cat("\nA. LPA options (ARI = agreement with the main solution; 1 = identical)\n"); print(resA, row.names = FALSE); fwrite(resA, file.path(OUT, "Table_options_LPA.csv"))
## ---- B. Louvain: edge weights
ws <- list("burden+entry" = c("burden", "entry"), "burden+peak" = c("burden", "peak_reached"))
resB <- rbindlist(lapply(names(ws), function(wn) { w <- ws[[wn]]; main <- louv(mean_tab, w)
  a <- ari(main, lapply(ds, function(d) louv(tab_of(I, d), w))); ar <- if (is.null(R)) NA else mean(ari(main, lapply(seq_len(nrow(R$burden)), function(j) louv(tab_of(R, j), w))))
  data.table(weights = wn, communities = max(main), ARI_draws = round(mean(a), 3), ARI_draws_lo = round(quantile(a, .025), 3), ARI_refits = round(ar, 3)) }))
cat("\nB. Louvain options\n"); print(resB, row.names = FALSE); fwrite(resB, file.path(OUT, "Table_options_Louvain.csv"))
## ---- C. consensus grouping for the most stable LPA option with K >= 3 (or best overall)
ok <- resA[!is.na(ARI_draws) & smallest >= 5]; best <- if (nrow(ok[K >= 3])) ok[K >= 3][which.max(ARI_draws)] else ok[which.max(ARI_draws)]
v <- sets[[best$indicators]]; k <- best$K; cat(sprintf("\nC. Consensus grouping for: %s, K = %d\n", best$indicators, k))
CO <- matrix(0, Z, Z); nn <- 0; for (d in ds) { cl <- lpa(tab_of(I, d), v, k); if (anyNA(cl)) next; CO <- CO + outer(cl, cl, "=="); nn <- nn + 1 }; CO <- CO / nn
hc <- hclust(as.dist(1 - CO), method = "average"); cons <- cutree(hc, k)
share <- sapply(1:Z, function(z) { m <- setdiff(which(cons == cons[z]), z); if (!length(m)) NA else mean(CO[z, m]) })
ct <- data.table(zone = zo, consensus_group = cons, mean_coassignment = round(share, 3), mean_burden = round(mean_tab$burden, 4), mean_peak = round(mean_tab$peak_reached, 1))
cs <- ct[, .(zones = .N, mean_coassignment = round(mean(mean_coassignment, na.rm = TRUE), 3), stable_zones_70 = sum(mean_coassignment >= 0.7, na.rm = TRUE),
             burden = round(mean(mean_burden), 4), peak = round(mean(mean_peak), 1)), by = consensus_group][order(-burden)]
print(cs, row.names = FALSE); cat(sprintf("Zones co-assigned with their group in >= 70%% of draws: %d of %d\n", sum(share >= 0.7, na.rm = TRUE), Z))
fwrite(ct, file.path(OUT, "Table_consensus_groups.csv")); fwrite(cs, file.path(OUT, "Table_consensus_summary.csv"))
cat("\nDone. Tables in", OUT, "\n")
