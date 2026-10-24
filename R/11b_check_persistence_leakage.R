# ==============================================================================
# Check_forward_persistence_leakage.R
# Purpose: (1) check that the forward one-step-ahead test (train 2015-2021, test 2022-2024) uses only
#          information available at each prediction time; (2) rerun it with the lagged covariates
#          standardized on the TRAINING years only; (3) compare with the saved result.
# Self-contained: rebuilds the main Baseline formulas and priors. Needs the data (dat) and the queen-contiguity
# matrix Q: taken from this session if present, otherwise from INLA_Rebuilt_From_Zero_OCCPOP.rds.
# Not run by its author (no R in his sandbox): expect a small first-run fix. Time: 4 INLA fits, ~5-15 min.
# ==============================================================================
suppressPackageStartupMessages({ library(INLA); library(data.table); library(pROC) })

# ---- 0. SET-UP: rebuild what the main script defines (only if it is not already in this session)
output_path <- if (exists("output_path")) output_path else file.path(PROJ)
S1FILE <- file.path(output_path, "INLA_Rebuilt_From_Zero_OCCPOP.rds")
cols <- c("ID_area", "ID_time", "year", "O", "E", "Y_bin", "O_int", "logE", "sin_month", "cos_month", "sin2_month", "cos2_month",
          "log_pop_s", "log_pop_zone_s", "ID_area_bym", "ID_time_rw", "ID_time_iid")
if (!exists("dat") || !all(cols %in% names(dat)) || !exists("Q")) {
  cat("Loading data and graph from", S1FILE, "\n"); S0 <- readRDS(S1FILE)
  if (!exists("dat") || !all(cols %in% names(dat))) dat <- as.data.table(S0$data)
  if (!exists("Q")) { qn <- intersect(c("Q", "graph", "W", "adj"), names(S0)); if (length(qn)) Q <- S0[[qn[1]]] }
}
dat <- as.data.table(dat); setorder(dat, ID_time, ID_area)
if (!"Y_bin" %in% names(dat)) dat[, Y_bin := as.integer(O > 0)]
if (!"O_int" %in% names(dat)) dat[, O_int := fifelse(O > 0, O, NA_integer_)]
if (!"logE"  %in% names(dat)) dat[, logE := log(E)]
miss <- setdiff(cols, names(dat))
if (length(miss)) stop("dat lacks columns: ", paste(miss, collapse = ", "), ". Run Steps 1-5 of the main script (data set-up) first.")
if (!exists("Q")) stop("The adjacency matrix Q is not in this session or in the Stage 1 file. Run the main script up to the graph step (queen contiguity) first.")
n_time <- max(dat$ID_time); n_area <- max(dat$ID_area); stopifnot(nrow(dat) == n_time * n_area)

# main-model priors (Selected_priors.csv: the defaults were retained)
prior_bym2     <- list(theta1 = list(prior = "pc.prec", param = c(1, 0.01)), theta2 = list(prior = "pc", param = c(0.5, 0.5)))
prior_bym2_int <- prior_bym2
prior_rw2      <- list(theta = list(prior = "pc.prec", param = c(20, 0.01)))     # occurrence RW2, U = 20
prior_rw2_int  <- list(theta = list(prior = "pc.prec", param = c(10, 0.01)))     # intensity RW2, U = 10
prior_iid      <- list(theta = list(prior = "pc.prec", param = c(0.5, 0.01)))
prior_iid_int  <- prior_iid
INT_FAM <- "zeroinflatednbinomial0"                                               # zero-truncated NB (main)
INT_CF  <- list(hyper = list(theta2 = list(initial = -15, fixed = TRUE)))         # zero probability fixed at ~0, default size prior
is_int_fam <- function(f) f %in% c("nbinomial", "zeroinflatednbinomial0")
size_of <- function(f) f$summary.hyperpar[grep("^size", rownames(f$summary.hyperpar))[1], "mean"]
mu_nb   <- function(f, idx = NULL) { lp <- f$summary.linear.predictor; m <- exp(lp$mean + 0.5 * lp$sd^2); if (is.null(idx)) m else m[idx] }

base <- paste("sin_month + cos_month + sin2_month + cos2_month",
              "+ f(ID_area_bym, model='bym2', graph=Q, scale.model=TRUE, constr=TRUE, hyper=prior_bym2)",
              "+ f(ID_time_rw, model='rw2', scale.model=TRUE, constr=TRUE, hyper=prior_rw2)",
              "+ f(ID_time_iid, model='iid', hyper=prior_iid)")
base_int <- gsub("prior_iid", "prior_iid_int", gsub("prior_bym2", "prior_bym2_int", gsub("prior_rw2", "prior_rw2_int", base, fixed = TRUE), fixed = TRUE), fixed = TRUE)
forms_occ <- list(Baseline = as.formula(paste("Y_bin ~", base, "+ log_pop_s")))
forms_int <- list(Baseline = as.formula(paste("O_int ~ offset(logE) +", base_int, "+ log_pop_zone_s")))
cat("Set-up done:", nrow(dat), "zone-months,", n_area, "zones,", n_time, "months.\n")

# same fast fit as the saved forward test (Step 3C / 8a)
ho_fit <- function(form, family, d) inla(form, family = family, data = d,
                                         control.family = if (is_int_fam(family)) INT_CF else list(),
                                         control.predictor = list(compute = TRUE, link = 1),
                                         control.compute = list(waic = FALSE, cpo = FALSE, dic = FALSE),
                                         control.inla = list(strategy = "gaussian", int.strategy = "eb"))

# ---- 1. build the lagged covariates exactly as Step 8, but keep the raw (unscaled) values
dP <- copy(dat); setorder(dP, ID_area, ID_time)
dP[, lag_own := as.numeric(shift(Y_bin, 1L, fill = 0L)), by = ID_area]
dP[, prev_cases := shift(O, 1L, fill = 0L), by = ID_area]
setorder(dP, ID_time, ID_area)                                   # back to the order of dat (time-major)
stopifnot(identical(dP$ID_area, dat$ID_area), identical(dP$ID_time, dat$ID_time))
Ym <- matrix(dP$Y_bin, nrow = n_time, ncol = n_area, byrow = TRUE)   # rows = months, cols = zones
Yb <- rbind(0, Ym[-n_time, , drop = FALSE])                      # previous month
nbshare <- (Yb %*% t(as.matrix(Q))) / matrix(rowSums(as.matrix(Q)), n_time, n_area, byrow = TRUE)
dP[, lag_nb_raw := as.vector(t(nbshare))]
dP[, ls_raw := log1p(prev_cases)]
dP[, prev_occ := as.numeric(prev_cases > 0)]

# Here dP is time-major (month, then zone), as dat; so byrow = TRUE gives rows = months, cols = zones.
# (Step 8 builds the same matrix while dP is still zone-major, with byrow = FALSE: same result.) CHECK 1c confirms it.

te <- which(dP$year >= 2022); tr <- setdiff(seq_len(nrow(dP)), te)
cat(sprintf("Training zone-months: %d (2015-2021) | test zone-months: %d (2022-2024)\n", length(tr), length(te)))

# ---- 2. LEAKAGE CHECKS (all should print TRUE)
chk <- list()
o <- dP[order(ID_area, ID_time)]
o[, O_prev_true := shift(O, 1L, fill = 0L), by = ID_area]
chk[["1a prev_cases = observed count of the PREVIOUS month (all rows)"]] <- isTRUE(all.equal(o$prev_cases, o$O_prev_true))
chk[["1b lag_own = outbreak in the PREVIOUS month (all rows)"]]          <- isTRUE(all.equal(o$lag_own, as.numeric(o$O_prev_true > 0)))
# neighbour share for one test row, recomputed by hand from the previous month
i  <- te[which(dP$ID_time[te] == min(dP$ID_time[te]))[1]]       # January 2022, first zone
z  <- dP$ID_area[i]; tt <- dP$ID_time[i]
nb <- which(as.matrix(Q)[z, ] != 0)
hand <- mean(dat$Y_bin[dat$ID_time == tt - 1 & dat$ID_area %in% nb])
chk[["1c lag_nb for Jan 2022 = share of neighbours with an outbreak in Dec 2021"]] <- isTRUE(all.equal(dP$lag_nb_raw[i], hand))
chk[["1d no test-period outcome is used as a lag for a TRAINING month"]] <- all(dP$ID_time[tr] < min(dP$ID_time[te]))
chk[["1e lags of the first test month come from the last training month"]] <- (tt - 1) == max(dP$ID_time[tr])
print(data.table(check = names(chk), ok = unlist(chk)))
cat("Reading: lags use the observed previous month, which is known at prediction time = no outcome leakage.\n",
    "Parameters are estimated once on 2015-2021 (one-step-ahead with fixed parameters).\n")

# ---- 3. two standardizations: full period (as saved) and TRAINING period only (no leakage)
z_full  <- function(x) (x - mean(x)) / sd(x)
z_train <- function(x) (x - mean(x[tr])) / sd(x[tr])
mk <- function(zf) { d <- copy(dP); d[, `:=`(lag_nb = zf(lag_nb_raw), prev_size = zf(ls_raw))]; d }
d_full <- mk(z_full); d_train <- mk(z_train)
cat(sprintf("\nScaling constants  lag_nb: full mean %.4f sd %.4f | train mean %.4f sd %.4f\n",
            mean(dP$lag_nb_raw), sd(dP$lag_nb_raw), mean(dP$lag_nb_raw[tr]), sd(dP$lag_nb_raw[tr])))
cat(sprintf("                   log(prev+1): full mean %.4f sd %.4f | train mean %.4f sd %.4f\n",
            mean(dP$ls_raw), sd(dP$ls_raw), mean(dP$ls_raw[tr]), sd(dP$ls_raw[tr])))

f_occ_P <- update(forms_occ[["Baseline"]], . ~ . + lag_own + lag_nb)
f_int_P <- update(forms_int[["Baseline"]], . ~ . + prev_occ + prev_size)
yb <- dat$Y_bin; NR <- nrow(dat)

fw_occ <- function(form, d0) {
  d <- copy(d0); d[te, Y_bin := NA_integer_]                      # outcome withheld in the test period
  f <- ho_fit(form, "binomial", d)
  p <- pmin(pmax(f$summary.fitted.values$mean[seq_len(NR)], 1e-9), 1 - 1e-9); prev <- mean(yb[tr])
  data.table(AUC = as.numeric(roc(yb[te], p[te], quiet = TRUE)$auc),
             logloss = -mean(yb[te] * log(p[te]) + (1 - yb[te]) * log(1 - p[te])),
             Brier = mean((yb[te] - p[te])^2), Brier_climatology = mean((yb[te] - prev)^2))
}
fw_int <- function(form, d0) {
  d <- copy(d0); d[te, O_int := NA_integer_]
  g <- ho_fit(form, INT_FAM, d)
  mu <- mu_nb(g, seq_len(NR)); sz <- size_of(g); tp <- te[!is.na(dat$O_int[te])]
  yi <- dat$O_int[tp]; mm <- mu[tp]; p0 <- (sz / (sz + mm))^sz
  data.table(logscore_truncNB = -mean(dnbinom(yi, size = sz, mu = mm, log = TRUE) - log1p(-p0)),
             MAE_log = median(abs(log(yi) - log(mm / (1 - p0)))), n_test_positive = length(tp))
}

cat("\nFitting 4 forward models (Gaussian approximation)...\n")
RES <- rbind(
  cbind(part = "occurrence", scaling = "full period (as saved)", fw_occ(f_occ_P, d_full)),
  cbind(part = "occurrence", scaling = "training period only",   fw_occ(f_occ_P, d_train)),
  cbind(part = "intensity",  scaling = "full period (as saved)", fw_int(f_int_P, d_full)),
  cbind(part = "intensity",  scaling = "training period only",   fw_int(f_int_P, d_train)), fill = TRUE)
print(RES[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)])

# ---- 4. compare with the saved result
sv <- file.path(output_path, "Table_S_persistence_forward_validation.csv")
if (file.exists(sv)) { cat("\nSaved forward result (Baseline + persistence row is the one to compare):\n"); print(fread(sv)) }
fwrite(RES, file.path(output_path, "Table_S_persistence_forward_leakage_check.csv"))
cat("\nREAD: if the 'training period only' AUC / Brier / log score differ from the full-period ones only in the\n",
    "third decimal, the scaling leak is negligible; report the training-only values or state the full-period scaling.\n",
    "Saved: Table_S_persistence_forward_leakage_check.csv\n")
