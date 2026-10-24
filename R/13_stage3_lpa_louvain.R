# ==============================================================================
# Stage2_LPA_Louvain_final.R  -  profiles (LPA) and communities (Louvain) from the FINAL HMM
# Lean version (v2): final draws file, 500 draws spread evenly over the chains, split-half check, stability over the 20 Stage 1 refits.
#
# Input  (made by HMM_check_ALL.R, rerun after the long fit; it reloads the cached fit in minutes):
#   Draws_updated_indicators_<TAG>.rds  = list(ind = list(burden, entry, persistence_active, peak_reached), zone_order)
#                                         each element: draws x 92 matrix
#   plus the Stage 1 file (for the adjacency graph) and the shapefile (for names and maps).
# Indicators: burden, entry (= gam12 + gam13), persistence_active, peak_reached.
# Not included on purpose: Local Moran, SaTScan, joint Stage 1 propagation (optional later), recovery test, K-states, hierarchical and seasonal runs.
# Wording: profiles / communities / classes of zones; "higher / lower estimated activity"; no transmission, hub, hotspot, corridor language.
#
# Run from R (VS Code terminal):  Rscript Stage2_LPA_Louvain_final.R
# Paths: edit the block below, or set options() before sourcing.
# ==============================================================================
need <- c("data.table", "mclust", "igraph", "ggplot2", "sf"); miss <- need[!sapply(need, requireNamespace, quietly = TRUE)]
if (length(miss)) stop("Missing R packages: ", paste(miss, collapse = ", "), ". Install them first (install.packages(...)).")
suppressPackageStartupMessages({ library(data.table); library(mclust); library(igraph); library(ggplot2); library(sf) })
set.seed(2026)

# ----------------------------- USER SETTINGS ----------------------------------
proj      <- getOption("cholera.proj", ".")
out       <- getOption("cholera.out",  file.path(proj, "improve_hmm"))
DRAWS_RDS <- getOption("cholera.draws", file.path(out, "Draws_updated_indicators_tau2_2_hi170_rmin0.5_hm_FINAL.rds"))   # main final HMM (made by Joint_20draws_long.R)
REFITS_RDS <- getOption("cholera.refits", { f <- file.path(proj, "joint_20draws_long", "Draws_joint20_long_indicators.rds"); if (file.exists(f)) f else file.path(proj, "joint_20draws", "Draws_joint_indicators.rds") })   # 20 Stage 1 refits
STAGE1    <- getOption("cholera.stage1", file.path(proj, "INLA_Rebuilt_From_Zero_OCCPOP.rds"))
SHP       <- getOption("cholera.shp", { f <- Sys.glob(file.path(proj, "*", "eth_admbnda_adm2*.shp")); if (length(f)) f[1] else "" })   # zone shapefile; "" skips maps
N_DRAWS_STAB <- getOption("cholera.ndraws", 500)   # posterior draws used for stability
N_NULL <- 1000; N_BOOT <- 1000; N_LABEL <- 10
K_RANGE <- 2:8; K_CHOSEN <- NA                     # NA = highest BIC; set 4 to force four profiles
UNCERTAIN_CUT <- 0.70                              # a zone keeping its profile / community / class in fewer than this share of draws is flagged uncertain
SAVE_FIGS <- TRUE
# -------------------------------------------------------------------------------
OUT2 <- file.path(out, "stage2_final"); dir.create(OUT2, showWarnings = FALSE, recursive = TRUE)
hd <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
CHK <- list(); chk <- function(step, check, pass, detail = "") { CHK[[length(CHK) + 1]] <<- data.table(step, check, result = ifelse(is.na(pass), "SKIP", ifelse(pass, "PASS", "FAIL")), detail); cat(sprintf("  [%s] %s: %s %s\n", ifelse(is.na(pass), "SKIP", ifelse(pass, "PASS", "FAIL")), step, check, detail)) }

# ---- load
stopifnot(file.exists(DRAWS_RDS), file.exists(STAGE1))
D <- readRDS(DRAWS_RDS); zone_order <- D$zone_order; Z <- length(zone_order); IND4 <- D$ind
vars <- getOption("cholera.vars", c("burden", "peak_reached")); stopifnot(all(vars %in% names(IND4)))   # indicators used for LPA (set options(cholera.vars = ...) to drop a redundant one)
cat("Indicators used for the LPA:", paste(vars, collapse = ", "), "\n")
nd <- nrow(IND4$burden); cat(sprintf("%d posterior draws, %d zones\n", nd, Z))
S1 <- readRDS(STAGE1); stopifnot(identical(S1$zone_order, zone_order))
shp0 <- if (nzchar(SHP)) tryCatch(st_read(SHP, quiet = TRUE), error = function(e) NULL) else NULL
zname <- if (!is.null(shp0)) setNames(as.character(if ("ADM2_EN" %in% names(shp0)) shp0$ADM2_EN else shp0$ADM2_PCODE), toupper(gsub("[^A-Z0-9]", "", shp0$ADM2_PCODE))) else NULL
nm_of <- function(z) { r <- if (!is.null(zname)) unname(zname[z]) else z; ifelse(is.na(r), z, r) }
draw_tab <- function(d, SRC = IND4) data.table(zone = zone_order, burden = SRC$burden[d, ], entry = SRC$entry[d, ], persistence_active = SRC$persistence_active[d, ], peak_reached = SRC$peak_reached[d, ])
ind_mean <- data.table(zone = zone_order, burden = colMeans(IND4$burden), entry = colMeans(IND4$entry), persistence_active = colMeans(IND4$persistence_active), peak_reached = colMeans(IND4$peak_reached))
chk("Input", "92 zones, finite indicators", Z == 92 && all(is.finite(as.matrix(ind_mean[, ..vars]))), sprintf("%d zones", Z))
chk("Input", "burden between 0 and 1", all(ind_mean$burden >= 0 & ind_mean$burden <= 1))

# ---- 5. LPA
hd("5. LATENT PROFILE ANALYSIS (Gaussian mixture, diagonal covariance with varying volume, K = 2-8)")
PR <- NULL   # transformation fixed on the posterior means (first call) and applied unchanged to every draw and refit
prep <- function(M) { if (is.null(PR)) { PR <<- lapply(vars, function(v) { x <- M[[v]]; sk <- mean((x - mean(x))^3) / sd(x)^3; lg <- is.finite(sk) && sk > 1
      off <- 0.5 * min(x[x > 0]); y <- if (lg) log(x + off) else x; list(lg = lg, off = off, m = mean(y), s = sd(y)) }); names(PR) <<- vars }
  sapply(vars, function(v) { p <- PR[[v]]; x <- M[[v]]; y <- if (p$lg) log(pmax(x, 0) + p$off) else x; (y - p$m) / p$s }) }
Xz <- prep(ind_mean); cat("Variables with skew > 1 are logged (offset = half the smallest positive value), then z-scored.\nCorrelation of the four indicators:\n"); print(round(cor(Xz), 2))
vif <- sapply(seq_along(vars), function(k) 1 / (1 - summary(lm(Xz[, k] ~ Xz[, -k]))$r.squared)); names(vif) <- vars; cat("VIF:", round(vif, 2), "\n")
ent <- function(m) { if (is.null(m) || is.null(m$z) || NCOL(m$z) < 2) return(NA_real_); p <- pmax(m$z, 1e-12); 1 + sum(p * log(p)) / (nrow(p) * log(ncol(p))) }
lpa <- lapply(K_RANGE, function(k) tryCatch(Mclust(Xz, G = k, modelNames = "VVI", verbose = FALSE), error = function(e) NULL))   # NULL when a K cannot be fitted
tabK <- data.table(K = K_RANGE, BIC = sapply(lpa, function(m) if (is.null(m)) NA_real_ else m$bic), entropy = sapply(lpa, ent),
                   smallest_n = sapply(lpa, function(m) if (is.null(m)) NA_integer_ else min(table(m$classification))), fitted_K = sapply(lpa, function(m) if (is.null(m)) NA_integer_ else m$G))
if (all(is.na(tabK$BIC))) stop("No LPA solution could be fitted with these indicators (check the collinearity above)")
print(tabK); fwrite(tabK, file.path(OUT2, "Table_S_LPA_selection.csv"))
Kbest <- if (is.na(K_CHOSEN)) tabK$K[which.max(tabK$BIC)] else K_CHOSEN; if (is.null(lpa[[match(Kbest, K_RANGE)]])) stop("chosen K could not be fitted"); cat("Chosen K:", Kbest, "(highest BIC in the mclust sign convention)\n")
lp <- lpa[[match(Kbest, K_RANGE)]]
chk("LPA", "entropy of chosen solution >= 0.8", tabK$entropy[match(Kbest, K_RANGE)] >= 0.8, sprintf("entropy %.3f", tabK$entropy[match(Kbest, K_RANGE)]))
chk("LPA", "smallest profile has >= 5 zones", tabK$smallest_n[match(Kbest, K_RANGE)] >= 5, sprintf("%d zones", tabK$smallest_n[match(Kbest, K_RANGE)]))
prof <- copy(ind_mean)[, `:=`(profile = lp$classification, max_prob = apply(lp$z, 1, max), zone_name = nm_of(zone))]
pm <- prof[, .(n = .N, burden = mean(burden), entry = mean(entry), persistence_active = mean(persistence_active), peak_reached = mean(peak_reached)), by = profile][order(-burden)]
pm[, profile_label := { lab <- rep("Intermediate", .N); lab[1] <- "Highest estimated activity"; lab[which.max(peak_reached)] <- "Largest peaks"; lab[which.min(burden)] <- "Lowest estimated activity"; lab }]
prof[pm, profile_label := i.profile_label, on = "profile"]
print(pm[, .(profile, profile_label, n, burden = round(burden, 3), entry = round(entry, 4), persistence_active = round(persistence_active, 3), peak_reached = round(peak_reached, 0))])
cat("Labels are automatic (rank on burden and peak); read the means above and edit the wording if the pattern differs.\n")
fwrite(prof, file.path(OUT2, "Table_S_LPA_zone_profiles.csv")); fwrite(pm, file.path(OUT2, "Table_S_LPA_profile_means.csv"))
for (pf in pm$profile) cat(sprintf("\n--- profile %d: %s (%d zones) ---\n%s\n", pf, pm[profile == pf, profile_label], pm[profile == pf, n], paste(prof[profile == pf, zone_name], collapse = ", ")))

# ---- 6. Network and Louvain
hd("6. ADJACENCY NETWORK, LOUVAIN COMMUNITIES, CENTRALITY")
Qm <- as.matrix(S1$graph_Q); diag(Qm) <- 0; Adj <- (Qm != 0) * 1; Adj <- (Adj + t(Adj) > 0) * 1
g0 <- graph_from_adjacency_matrix(Adj, mode = "undirected"); el <- as_edgelist(g0, names = FALSE)
cat(sprintf("Network: %d zones, %d borders, %d connected components\n", vcount(g0), ecount(g0), components(g0)$no))
mm <- function(x) (x - min(x)) / (max(x) - min(x) + 1e-12)
ew <- function(M) { b <- mm(M$burden); v <- mm(M$peak_reached); 0.5 * ((b[el[, 1]] + b[el[, 2]]) / 2 + (v[el[, 1]] + v[el[, 2]]) / 2) + 1e-3 }   # border weight = mean of burden and peak reached of the two zones
mk_g <- function(w) { g <- g0; E(g)$weight <- w; g }
g <- mk_g(ew(ind_mean)); lv <- cluster_louvain(g, weights = E(g)$weight); memb <- membership(lv); Qmod <- modularity(lv)
cat(sprintf("Louvain: %d communities, modularity Q = %.3f\n", max(memb), Qmod)); chk("Louvain", "modularity Q > 0.3", Qmod > 0.3, sprintf("Q = %.3f", Qmod))
Qnull <- replicate(N_NULL, { gn <- rewire(g0, keeping_degseq(loops = FALSE, niter = 10 * ecount(g0))); E(gn)$weight <- sample(E(g)$weight); modularity(cluster_louvain(gn, weights = E(gn)$weight)) })
cat(sprintf("Null networks (degree-preserving, weights shuffled): mean Q %.3f (SD %.3f) | Z = %.1f | p = %.4f\n", mean(Qnull), sd(Qnull), (Qmod - mean(Qnull)) / sd(Qnull), (1 + sum(Qnull >= Qmod)) / (N_NULL + 1)))
ari_edge <- replicate(N_BOOT, adjustedRandIndex(memb, membership(cluster_louvain(mk_g(E(g)$weight * runif(ecount(g), .8, 1.2))))))
cat(sprintf("Edge-weight perturbation (+-20%%, %d runs): mean ARI %.3f (95%% interval %.3f-%.3f)\n", N_BOOT, mean(ari_edge), quantile(ari_edge, .025), quantile(ari_edge, .975)))
cen <- eigen_centrality(g, weights = E(g)$weight, scale = TRUE)$vector
ctab <- data.table(zone = zone_order, zone_name = nm_of(zone_order), community = as.integer(memb), centrality = round(cen, 3))[order(-centrality)]
cat("\nMost central zones (position in the border network weighted by estimated activity; this does not show that a zone drives cases elsewhere):\n"); print(ctab[1:N_LABEL])

# 6 -> 3 classes (explicit rule on community means of z-scored log burden and log peak)
classify <- function(mb, M) {
  zb <- as.numeric(scale(log(M$burden + .5 * min(M$burden[M$burden > 0])))); zp <- as.numeric(scale(log(M$peak_reached + .5 * min(M$peak_reached[M$peak_reached > 0]))))
  cmn <- data.table(community = as.integer(mb), zb, zp)[, .(burden_z = mean(zb), peak_z = mean(zp)), by = community]
  cmn[, class := fifelse(peak_z > 0, "Higher-peak communities", fifelse(burden_z > 0, "Higher-burden communities", "Lower-activity communities"))]
  cmn$class[match(as.integer(mb), cmn$community)] }
cls_main <- classify(memb, ind_mean)
cat("Merging rule: community mean peak z > 0 -> Higher-peak; else burden z > 0 -> Higher-burden; else Lower-activity.\n"); print(table(cls_main))
comm <- data.table(zone = zone_order, zone_name = nm_of(zone_order), community = as.integer(memb), class = cls_main, centrality = cen)

# ---- 6a. stability across posterior draws (label-free ARI + zone-by-zone share)
hd("6a. STABILITY OVER POSTERIOR DRAWS")
ds <- unique(round(seq(1, nd, length.out = min(N_DRAWS_STAB, nd)))); n <- length(ds)   # evenly spaced, so all chains are represented
run_set <- function(rows, SRC) { k <- length(rows); LP <- CM <- matrix(NA_integer_, k, Z); CL <- matrix(NA_character_, k, Z)
  for (b in seq_len(k)) { M <- draw_tab(rows[b], SRC)
    LP[b, ] <- tryCatch(predict(lp, newdata = prep(M))$classification, error = function(e) NA_integer_)   # profiles held fixed; each draw/refit assigned to them
    mb <- as.integer(membership(cluster_louvain(mk_g(ew(M))))); CM[b, ] <- mb; CL[b, ] <- classify(mb, M) }
  list(LP = LP, CM = CM, CL = CL) }
RS <- run_set(ds, IND4); LP <- RS$LP; CM <- RS$CM; CL <- RS$CL
ari <- function(Mx, main) apply(Mx, 1, function(v) if (anyNA(v)) NA_real_ else adjustedRandIndex(main, v))
co <- function(Mx, main) sapply(1:Z, function(z) { mem <- setdiff(which(main == main[z]), z); if (!length(mem)) NA_real_ else mean(apply(Mx, 1, function(v) if (anyNA(v)) NA_real_ else mean(v[mem] == v[z])), na.rm = TRUE) })
a_lpa <- ari(LP, lp$classification); a_com <- ari(CM, as.integer(memb)); a_cls <- apply(CL, 1, function(v) adjustedRandIndex(cls_main, v))
sm <- function(x) sprintf("mean ARI %.2f (95%% interval %.2f-%.2f)", mean(x, na.rm = TRUE), quantile(x, .025, na.rm = TRUE), quantile(x, .975, na.rm = TRUE))
cat(sprintf("Profiles    over %d draws: %s\nCommunities over %d draws: %s\nThree classes over %d draws: %s\n", n, sm(a_lpa), n, sm(a_com), n, sm(a_cls)))
zs <- data.table(zone = zone_order, share_profile = co(LP, lp$classification), share_community = co(CM, as.integer(memb)), share_class = colMeans(sweep(CL, 2, cls_main, "=="), na.rm = TRUE))
zs[, uncertain := pmin(share_profile, share_community, share_class, na.rm = TRUE) < UNCERTAIN_CUT]
cat(sprintf("Zones keeping profile, community and class in >= %.0f%% of draws: %d of %d; flagged uncertain: %d\n", 100 * UNCERTAIN_CUT, sum(!zs$uncertain), Z, sum(zs$uncertain)))
cat("Uncertain zones:", paste(nm_of(zs[uncertain == TRUE, zone]), collapse = ", "), "\n")
h1 <- seq_len(floor(n / 2)); h2 <- setdiff(seq_len(n), h1)
half <- data.table(measure = c("profiles ARI", "communities ARI", "classes ARI"),
  first_half = c(mean(a_lpa[h1], na.rm = TRUE), mean(a_com[h1]), mean(a_cls[h1])), second_half = c(mean(a_lpa[h2], na.rm = TRUE), mean(a_com[h2]), mean(a_cls[h2])))
cat("\nSplit-half check (are 500 draws enough? the two halves should agree):\n"); print(half[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
cat(sprintf("Zone share keeping its class: correlation between halves %.2f\n", cor(colMeans(sweep(CL[h1, , drop = FALSE], 2, cls_main, "==")), colMeans(sweep(CL[h2, , drop = FALSE], 2, cls_main, "==")))))
fwrite(half, file.path(OUT2, "Table_S_stability_split_half.csv"))
chk("Stability", "split halves agree (difference in mean ARI < 0.05)", max(abs(half$first_half - half$second_half)) < 0.05, sprintf("max difference %.3f", max(abs(half$first_half - half$second_half))))

## ---- 6b. stability over the 20 refits on joint Stage 1 draws (first-stage uncertainty)
hd("6b. STABILITY OVER THE 20 STAGE 1 REFITS")
if (file.exists(REFITS_RDS)) {
  R20 <- readRDS(REFITS_RDS); stopifnot(identical(as.character(R20$zone_order), as.character(zone_order))); I20 <- R20$ind; nr <- nrow(I20$burden)
  RR <- run_set(seq_len(nr), I20)
  r_lpa <- ari(RR$LP, lp$classification); r_com <- ari(RR$CM, as.integer(memb)); r_cls <- apply(RR$CL, 1, function(v) adjustedRandIndex(cls_main, v))
  cat(sprintf("Source: %s\nProfiles    over %d refits: %s\nCommunities over %d refits: %s\nThree classes over %d refits: %s\n", REFITS_RDS, nr, sm(r_lpa), nr, sm(r_com), nr, sm(r_cls)))
  zs[, `:=`(refit_share_profile = co(RR$LP, lp$classification), refit_share_community = co(RR$CM, as.integer(memb)), refit_share_class = colMeans(sweep(RR$CL, 2, cls_main, "=="), na.rm = TRUE))]
  REF <- data.table(item = c("ARI_lpa_refits_mean", "ARI_lpa_refits_lo", "ARI_lpa_refits_hi", "ARI_comm_refits_mean", "ARI_comm_refits_lo", "ARI_comm_refits_hi", "ARI_class_refits_mean", "n_refits"),
    value = c(mean(r_lpa, na.rm = TRUE), quantile(r_lpa, .025, na.rm = TRUE), quantile(r_lpa, .975, na.rm = TRUE), mean(r_com), quantile(r_com, .025), quantile(r_com, .975), mean(r_cls), nr))
  fwrite(REF, file.path(OUT2, "Table_S_stability_refits.csv"))
  REFTXT <- sprintf("Across the %d refits on joint first-stage draws, profiles were reproduced with %s, communities with %s and the three classes with mean ARI %.2f.", nr, sm(r_lpa), sm(r_com), mean(r_cls))
} else { cat("Refit file not found:", REFITS_RDS, "- run Joint_20draws_long.R first; skipped.\n"); REFTXT <- character(0) }
chk("Stability", "no missing ARI for communities and classes", !anyNA(a_com) && !anyNA(a_cls), sprintf("%d draws", n))
comm <- merge(merge(comm, zs, by = "zone"), prof[, .(zone, profile, profile_label, max_prob, burden, entry, persistence_active, peak_reached)], by = "zone")
fwrite(comm, file.path(OUT2, "Table_S_zone_profile_community_class.csv"))
fwrite(data.table(item = c("Q", "null_mean", "null_sd", "Z", "p_null", "ARI_edge_mean", "ARI_edge_lo", "ARI_edge_hi", "ARI_comm_draws_mean", "ARI_comm_lo", "ARI_comm_hi", "ARI_lpa_mean", "ARI_lpa_lo", "ARI_lpa_hi", "ARI_class_mean", "n_uncertain"),
  value = c(Qmod, mean(Qnull), sd(Qnull), (Qmod - mean(Qnull)) / sd(Qnull), (1 + sum(Qnull >= Qmod)) / (N_NULL + 1), mean(ari_edge), quantile(ari_edge, .025), quantile(ari_edge, .975),
            mean(a_com, na.rm = TRUE), quantile(a_com, .025, na.rm = TRUE), quantile(a_com, .975, na.rm = TRUE), mean(a_lpa, na.rm = TRUE), quantile(a_lpa, .025, na.rm = TRUE), quantile(a_lpa, .975, na.rm = TRUE), mean(a_cls), sum(zs$uncertain))),
  file.path(OUT2, "Table_S_stability_summary.csv"))

# ---- 8. maps and figure
hd("8. MAPS")
if (!is.null(shp0)) {
  key <- toupper(gsub("[^A-Z0-9]", "", shp0$ADM2_PCODE)); shp <- shp0[match(zone_order, key), ]
  if (anyNA(match(zone_order, key))) cat("Some zone codes were not found in the shapefile; those zones are blank.\n")
  shp$profile <- factor(comm$profile_label[match(zone_order, comm$zone)]); shp$class <- factor(comm$class[match(zone_order, comm$zone)])
  shp$community <- factor(comm$community[match(zone_order, comm$zone)]); shp$uncertain <- comm$uncertain[match(zone_order, comm$zone)]
  mp <- function(var, ttl) ggplot(shp) + geom_sf(aes(fill = .data[[var]]), colour = "grey30", linewidth = 0.15) + geom_sf(data = shp[which(shp$uncertain), ], fill = NA, colour = "black", linewidth = 0.6) + labs(title = ttl, fill = NULL, caption = "Black outline = zone not stable across posterior draws") + theme_void() + theme(legend.position = "bottom")
  for (v in list(c("profile", "Latent profiles"), c("community", "Louvain communities"), c("class", "Three community classes"))) {
    p <- mp(v[1], v[2]); print(p)
    if (SAVE_FIGS) { ggsave(file.path(OUT2, paste0("Fig_map_", v[1], ".png")), p, width = 8, height = 8, dpi = 600); ggsave(file.path(OUT2, paste0("Fig_map_", v[1], ".pdf")), p, width = 8, height = 8) } }
} else cat("No shapefile given (set options(cholera.shp = '.../file.shp')); maps skipped.\n")

# ---- checklist and sentences
hd("9. CHECKLIST")
CH <- rbindlist(CHK); print(CH, row.names = FALSE); fwrite(CH, file.path(OUT2, "Table_S_stage2_checklist.csv"))
txt <- c(sprintf("LPA: %d profiles chosen by BIC (entropy %.2f); across %d posterior draws the profiles were reproduced with %s.", Kbest, tabK$entropy[match(Kbest, K_RANGE)], n, sm(a_lpa)),
         sprintf("Louvain: %d communities, modularity %.3f versus %.3f (SD %.3f) in degree-preserving null networks (p = %.4f); across %d posterior draws %s; %d of %d zones kept their profile, community and class in at least %.0f%% of draws.",
                 max(memb), Qmod, mean(Qnull), sd(Qnull), (1 + sum(Qnull >= Qmod)) / (N_NULL + 1), n, sm(a_com), sum(!zs$uncertain), Z, 100 * UNCERTAIN_CUT))
txt <- c(txt, sprintf("Split-half check: mean ARI for classes %.2f and %.2f in the two halves of the draws.", half$first_half[3], half$second_half[3]), REFTXT)
writeLines(txt, file.path(OUT2, "Stage2_sentences.txt")); cat("\n", paste(txt, collapse = "\n"), "\n")
cat("\nDone. Outputs in:", OUT2, "\n")
