# ==============================================================================
# HMM_states_2_3_4.R  -  how many hidden states? Compares the population HMM with 2, 3 and 4 states
# Same model as the main analysis (NB likelihood, offset E/mean(E), b_pop * standardised log population, BYM2-informed zone deviations,
# positive-gap ordering, data-informed Dirichlet transition priors with concentration 10), only the number of states changes.
# For each K: chains on all months (fit, WAIC) + chains with 2024 hidden (hold-out), one chain per worker. With CORES=40: 8 + 5 per K = 39 chains.
# Each chain: 20,000 burn-in + 20,000 kept (thin 20) = 1,000 draws. Chains are saved when done (resumable: rerun the same command).
# Run:   cd ~/Documents/spatial_anlysis/update_HMM
#        Rscript -e 'options(cholera.check = TRUE); source("HMM_states_2_3_4.R")'      # quick check (seconds)
#        Rscript -e 'options(cholera.smoke = TRUE); source("HMM_states_2_3_4.R")'      # end-to-end test (a few minutes), results meaningless
#        mkdir -p test_states_K234 ; CORES=40 setsid nohup nice -n 19 Rscript HMM_states_2_3_4.R > test_states_K234/log.txt 2>&1 < /dev/null &
#        tail -f test_states_K234/log.txt
# Output (update_HMM/test_states_K234): Table_S_number_of_states.csv, States_sentences.txt, chain files
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(parallel); library(coda) })
proj <- getOption("cholera.proj", ".")
OUT  <- file.path(proj, "test_states_K234")
KS <- c(2, 3, 4); BURN <- 20000; NIT <- 20000; THIN <- 20; T_HOLD <- 12
SMOKE <- isTRUE(getOption("cholera.smoke", FALSE))
if (SMOKE) { BURN <- 100; NIT <- 100; THIN <- 2; OUT <- file.path(proj, "test_states_K234_smoke") }
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
NW <- as.integer(getOption("cholera.cores", Sys.getenv("CORES", 12)))
PER_K <- max(4, NW %/% length(KS))                               # chains per K (40 cores -> 13)
NCH_FULL <- as.integer(Sys.getenv("NCH_FULL", max(2, ceiling(0.6 * PER_K)))); NCH_HOLD <- as.integer(Sys.getenv("NCH_HOLD", max(2, PER_K - NCH_FULL)))
say <- function(...) cat(format(Sys.time(), "%d %b %H:%M:%S"), ..., "\n")
set.seed(2026)

## ---- data (same files as the main analysis)
S0 <- readRDS(file.path(proj, "INLA_Rebuilt_From_Zero_OCCPOP.rds")); dat <- as.data.table(S0$data); setorder(dat, ID_area, ID_time)
zone_order <- S0$zone_order; Z <- length(zone_order); TT <- max(dat$ID_time); stopifnot(nrow(dat) == Z * TT)
Y <- matrix(dat$O, Z, TT, byrow = TRUE); Em <- matrix(dat$E, Z, TT, byrow = TRUE); OFF <- Em / mean(Em)
zpop <- as.numeric(scale(log(as.numeric(tapply(dat$monthly_population, dat$ID_area, mean)))))
## Stage 1 zone fields: posterior mean and SD (minimum SD 0.05) from the joint Stage 1 draws
JF <- readRDS(file.path(proj, "joint_fields_20.rds")); stopifnot(identical(as.character(JF$zone_order), as.character(zone_order)))
d_occ <- rowMeans(JF$DO); d_int <- rowMeans(JF$DI)
prec_occ <- 1 / pmax(apply(JF$DO, 1, sd), 0.05)^2; prec_int <- 1 / pmax(apply(JF$DI, 1, sd), 0.05)^2
d_occ <- d_occ - log(mean(exp(d_occ))); d_int <- d_int - log(mean(exp(d_int)))              # mean of exp(delta) = 1

## ---- state anchors, empirical starting states and transition priors for each K
pos <- Y[Y > 0]; q95 <- as.numeric(quantile(pos, 0.95)); q75 <- as.numeric(quantile(pos, 0.75))
setup <- function(K) {
  st <- matrix(1L, Z, TT)
  if (K == 2) { st[Y > 0] <- 2L; anch <- c(0.01, mean(pos)) }
  if (K == 3) { st[Y > 0 & Y < q95] <- 2L; st[Y >= q95] <- 3L; anch <- c(0.01, median(pos[pos < q95]), 170) }
  if (K == 4) { st[Y > 0 & Y < q75] <- 2L; st[Y >= q75 & Y < q95] <- 3L; st[Y >= q95] <- 4L
                a3 <- median(pos[pos >= q75 & pos < q95]); anch <- c(0.01, median(pos[pos < q75]), a3, max(170, 2 * a3)) }   # top anchor kept above the third
  trm <- matrix(0, K, K); for (z in 1:Z) for (t in 2:TT) trm[st[z, t - 1], st[z, t]] <- trm[st[z, t - 1], st[z, t]] + 1
  self <- diag(trm) / pmax(rowSums(trm), 1)
  for (k in 3:K) if (K >= 3 && sum(trm[k, ]) < 100) self[k] <- self[k - 1]                     # few transitions: borrow from the state below (as in the main model)
  CONC <- 10; alpha <- matrix(0, K, K)
  for (i in 1:K) { alpha[i, ] <- (1 - self[i]) * CONC / (K - 1); alpha[i, i] <- max(self[i], 0.05) * CONC }
  list(K = K, st = st, anch = anch, alpha = alpha, self = self, ntr = rowSums(trm))
}
SET <- lapply(setNames(KS, paste0("K", KS)), setup)

## ---- JAGS model text for K states
model_text <- function(K, hold) {
  gaps <- if (K >= 3) paste0(sapply(3:K, function(k) sprintf("  eta[%d] ~ dnorm(log(anch[%d] - anch[%d]), 1)\n  tau[%d] ~ dgamma(2, 2)\n", k, k, k - 1, k)), collapse = "") else ""
  zg <- if (K >= 3) paste0(sapply(3:K, function(k) sprintf("    lg[z, %d] ~ dnorm(eta[%d], tau[%d])\n    mu[z, %d] <- mu[z, %d] + exp(lg[z, %d])\n", k, k, k, k, k - 1, k)), collapse = "") else ""
  paste0("model {\n  pi0[1:K] ~ ddirch(pi_a[1:K])\n  b_pop ~ dnorm(0, 1)\n  sd_g12 ~ dunif(1, 2 * (anch[2] - anch[1]))\n  tau_g12 <- pow(sd_g12, -2)\n", gaps,
    "  for (z in 1:Z) {\n    mu1[z] ~ dnorm(anch[1], 1e4)T(1e-4, )\n    g12[z] ~ dnorm(anch[2] - anch[1], tau_g12)T(0, )\n    mu[z, 1] <- mu1[z]\n    mu[z, 2] <- mu1[z] + g12[z]\n", zg,
    "    delta[z, 1] ~ dnorm(d_occ[z], prec_occ[z])\n    for (k in 2:K) { delta[z, k] ~ dnorm(d_int[z], prec_int[z]) }\n",
    "    for (k in 1:K) { lam_s[z, k] <- mu[z, k] * exp(delta[z, k]) }\n    r[z] ~ dgamma(2, 0.5)T(0.5, )\n",
    "    for (i in 1:K) { gam[z, i, 1:K] ~ ddirch(alpha[i, 1:K]) }\n    S[z, 1] ~ dcat(pi0[1:K])\n    for (t in 2:TT) { S[z, t] ~ dcat(gam[z, S[z, t - 1], 1:K]) }\n",
    "    for (t in 1:TT) {\n      lam[z, t] <- max(lam_s[z, S[z, t]] * OFF[z, t] * exp(b_pop * zpop[z]), 1e-6)\n      p[z, t] <- r[z] / (r[z] + lam[z, t])\n      Y[z, t] ~ dnegbin(p[z, t], r[z])\n    }\n",
    if (hold) "    for (j in 1:NH) { Yh[z, j] <- Y[z, TT - NH + j] }\n" else "", "  }\n}\n")
}

lse <- function(M) { m <- apply(M, 1, max); m + log(rowSums(exp(M - m))) }                    # row-wise log-sum-exp

## ---- one chain; heavy summaries are computed here so that only small results are saved
one_chain <- function(task) {
  K <- task$K; hold <- task$hold; ch <- task$ch; s <- SET[[paste0("K", K)]]
  f <- file.path(OUT, sprintf("K%d_%s_chain%d.rds", K, if (hold) "holdout" else "full", ch))
  if (file.exists(f)) { r <- tryCatch(readRDS(f), error = function(e) NULL); if (!is.null(r)) return(r) }
  suppressPackageStartupMessages({ library(rjags); library(coda) })
  mf <- file.path(OUT, sprintf("model_K%d_%s_chain%d.jags", K, if (hold) "holdout" else "full", ch)); writeLines(model_text(K, hold), mf)
  Yd <- Y; if (hold) Yd[, (TT - T_HOLD + 1):TT] <- NA
  jd <- list(Y = Yd, OFF = OFF, Z = Z, TT = TT, K = K, alpha = s$alpha, anch = s$anch, pi_a = c(8, rep(1, K - 1)), zpop = zpop,
             d_occ = d_occ, d_int = d_int, prec_occ = prec_occ, prec_int = prec_int)
  if (hold) jd$NH <- T_HOLD
  m <- jags.model(mf, data = jd, inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 9100 + 100 * K + 10 * hold + ch, S = s$st),
                  n.chains = 1, quiet = TRUE, n.adapt = 1000)
  update(m, BURN)
  mon <- c("S", "gam", "mu", "delta", "b_pop", "r", "pi0"); if (hold) mon <- c(mon, "Yh")
  x <- as.matrix(coda.samples(m, mon, n.iter = NIT, thin = THIN)[[1]]); nd <- nrow(x)
  col <- function(fmt, ...) x[, sprintf(fmt, ...), drop = FALSE]
  A  <- array(col("S[%d,%d]", rep(1:Z, TT), rep(1:TT, each = Z)), c(nd, Z, TT))
  bp <- x[, "b_pop"]; rr <- col("r[%d]", 1:Z)
  LS <- array(NA_real_, c(nd, Z, K)); for (k in 1:K) LS[, , k] <- col("mu[%d,%d]", 1:Z, k) * exp(col("delta[%d,%d]", 1:Z, k)) * exp(outer(bp, zpop))
  G  <- array(NA_real_, c(nd, Z, K, K)); for (i in 1:K) for (j in 1:K) G[, , i, j] <- col("gam[%d,%d,%d]", 1:Z, i, j)
  P0 <- col("pi0[%d]", 1:K)
  res <- list(K = K, hold = hold, ch = ch, nd = nd,
              par = cbind(b_pop = bp, sapply(1:K, function(k) rowMeans(G[, , k, k])), sapply(1:K, function(k) rowMeans(log(LS[, , k])))),
              mu_nat = sapply(1:K, function(k) apply(col("mu[%d,%d]", 1:Z, k), 1, median)),
              self = sapply(1:K, function(k) rowMeans(G[, , k, k])),
              occ = sapply(1:K, function(k) apply(A == k, 1, mean)),
              top_months = apply(A == K, c(1, 2), sum))                                               # draws x zones: months in the top state
  if (!hold) {
    ## forward algorithm (states summed out): zone log-likelihood and one-step-ahead predictive log density for every zone-month
    llz <- matrix(0, nd, Z); lp1 <- array(0, c(nd, Z, TT))
    for (z in 1:Z) {
      lf <- function(t) sapply(1:K, function(k) dnbinom(Y[z, t], size = rr[, z], mu = pmax(LS[, z, k] * OFF[z, t], 1e-6), log = TRUE))
      a <- log(P0) + lf(1); prev <- lse(a); lp1[, z, 1] <- prev
      for (t in 2:TT) {
        e <- lf(t); an <- sapply(1:K, function(j) lse(a + log(G[, z, , j]))) + e
        cur <- lse(an); lp1[, z, t] <- cur - prev; prev <- cur; a <- an }
      llz[, z] <- prev }
    res$llz <- llz; res$lp1 <- matrix(lp1, nd)                                                    # nd x (Z*TT)
    ## fit: posterior predictive check, randomized quantile residuals, in-sample error
    LAM <- array(NA_real_, c(nd, Z, TT)); offA <- array(rep(OFF, each = nd), c(nd, Z, TT))
    for (k in 1:K) { idx <- A == k; LAM[idx] <- (array(LS[, , k], c(nd, Z, TT)) * offA)[idx] }
    LAM <- pmax(LAM, 1e-6)
    sub <- sort(sample(nd, min(200, nd))); RSZ <- array(rr, c(nd, Z, TT))
    rep_stats <- t(sapply(sub, function(d) { yr <- rnbinom(Z * TT, size = RSZ[d, , ], mu = LAM[d, , ])
      c(zeros = mean(yr == 0), total = sum(yr), max = max(yr), meanpos = mean(yr[yr > 0]), high = mean(yr >= q95)) }))
    res$ppc <- rep_stats
    res$rqr_sd <- mean(sapply(sub[1:min(100, length(sub))], function(d) { mu_ <- LAM[d, , ]; sz <- RSZ[d, , ]
      lo <- pnbinom(Y - 1, size = sz, mu = mu_); hi <- pnbinom(Y, size = sz, mu = mu_); u <- lo + runif(length(Y)) * (hi - lo)
      sd(qnorm(pmin(pmax(u, 1e-10), 1 - 1e-10))) }))
    res$lam_mean <- apply(LAM, c(2, 3), mean)
  } else {
    res$yh <- col("Yh[%d,%d]", rep(1:Z, T_HOLD), rep(1:T_HOLD, each = Z))                         # nd x (Z*T_HOLD), predictive draws for 2024
  }
  saveRDS(res, paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f); res
}

if (isTRUE(getOption("cholera.check", FALSE))) {
  for (K in KS) { s <- SET[[paste0("K", K)]]
    say(sprintf("K = %d: anchors %s | prior self-transitions %s | transitions out of each state %s", K, paste(round(s$anch, 1), collapse = " / "),
                paste(round(s$self, 3), collapse = " / "), paste(s$ntr, collapse = " / "))) }
  say("rjags:", requireNamespace("rjags", quietly = TRUE), "| loo:", requireNamespace("loo", quietly = TRUE), "| cores:", detectCores(), "| workers:", NW)
  say("Output folder:", OUT); stop("check finished (this is not an error): now start the run.", call. = FALSE) }

tasks <- unlist(lapply(KS, function(K) unlist(lapply(c(FALSE, TRUE), function(h) lapply(seq_len(if (h) NCH_HOLD else NCH_FULL), function(ch) list(K = K, hold = h, ch = ch))), recursive = FALSE)), recursive = FALSE)
say(sprintf("Running %d chains (K = %s; per K %d full + %d hold-out; %d burn-in + %d kept, thin %d) on %d workers", length(tasks),
            paste(KS, collapse = ", "), NCH_FULL, NCH_HOLD, BURN, NIT, THIN, min(NW, length(tasks))))
cl <- makeCluster(min(NW, length(tasks)), port = 23000L + as.integer(Sys.getpid()) %% 15000L)
clusterExport(cl, c("one_chain", "OUT", "SET", "model_text", "lse", "Y", "OFF", "Z", "TT", "zpop", "d_occ", "d_int", "prec_occ", "prec_int",
                    "BURN", "NIT", "THIN", "T_HOLD", "q95"), envir = environment())
RES <- parLapplyLB(cl, tasks, function(t) tryCatch(one_chain(t), error = function(e) list(error = conditionMessage(e), task = t))); stopCluster(cl)
bad <- Filter(function(r) !is.null(r$error), RES)
if (length(bad)) { for (b in bad) say(sprintf("FAILED K=%d hold=%s chain %d: %s", b$task$K, b$task$hold, b$task$ch, b$error)); stop("some chains failed") }
say("All chains finished. Summarising ...")

## ---- summary for each K
crps <- function(dr, y) { n <- length(dr); s <- sort(dr); mean(abs(dr - y)) - sum((2 * seq_len(n) - n - 1) * s) / n^2 }
get <- function(K, hold) Filter(function(r) r$K == K && r$hold == hold, RES)
yobs <- c(mean(Y == 0), sum(Y), max(Y), mean(Y[Y > 0]), mean(Y >= q95))
rows <- list(); PW <- list(); LZ <- list()
for (K in KS) {
  fl <- get(K, FALSE); ho <- get(K, TRUE)
  rh <- tryCatch(max(gelman.diag(mcmc.list(lapply(fl, function(r) mcmc(r$par))), autoburnin = FALSE, multivariate = FALSE)$psrf[, 1], na.rm = TRUE), error = function(e) NA)
  RB <- function(L, v) do.call(rbind, lapply(L, `[[`, v)); llz <- RB(fl, "llz"); lp1 <- RB(fl, "lp1"); LZ[[paste0("K", K)]] <- llz; PW[[paste0("K", K)]] <- lp1
  waic_pts <- function(L) { m <- apply(L, 2, max); lppd <- m + log(colMeans(exp(sweep(L, 2, m)))); pw <- apply(L, 2, var); -2 * (lppd - pw) }
  wz <- waic_pts(llz); w1 <- waic_pts(lp1)
  ppc <- RB(fl, "ppc"); pv <- sapply(1:5, function(i) mean(ppc[, i] >= yobs[i]))
  inside <- sum(sapply(1:5, function(i) yobs[i] >= quantile(ppc[, i], 0.025) && yobs[i] <= quantile(ppc[, i], 0.975)))
  lm <- Reduce(`+`, lapply(fl, `[[`, "lam_mean")) / length(fl)
  yh <- RB(ho, "yh"); yt <- as.vector(Y[, (TT - T_HOLD + 1):TT])
  cr <- mean(sapply(seq_along(yt), function(i) crps(yh[, i], yt[i])))
  cov <- mean(sapply(seq_along(yt), function(i) { q <- quantile(yh[, i], c(0.025, 0.975)); yt[i] >= q[1] && yt[i] <= q[2] }))
  occ <- colMeans(RB(fl, "occ")); self <- colMeans(RB(fl, "self")); mun <- colMeans(RB(fl, "mu_nat"))
  topm <- colMeans(RB(fl, "top_months"))
  rows[[length(rows) + 1]] <- data.table(K = K, max_rhat = round(rh, 3),
    WAIC_zone = round(sum(wz), 1), WAIC_zone_se = round(sqrt(Z * var(wz)), 1),
    WAIC_onestep = round(sum(w1), 1), WAIC_onestep_se = round(sqrt(length(w1) * var(w1)), 1),
    ppc_inside_of_5 = inside, ppc_p = paste(sprintf("%.2f", pv), collapse = " / "), rqr_sd = round(mean(sapply(fl, `[[`, "rqr_sd")), 2),
    MAE = round(mean(abs(lm - Y)), 2), RMSE = round(sqrt(mean((lm - Y)^2)), 1),
    holdout_CRPS = round(cr, 2), holdout_cov95 = round(100 * cov, 1),
    state_means = paste(signif(mun, 3), collapse = " / "), self_transitions = paste(sprintf("%.3f", self), collapse = " / "),
    expected_stay = paste(sprintf("%.1f", 1 / (1 - self)), collapse = " / "), share_zone_months = paste(sprintf("%.3f", occ), collapse = " / "),
    smallest_state_share = round(min(occ), 4), zones_lt3_top_months = sum(topm < 3))
}
TAB <- rbindlist(rows)
TAB[, dWAIC_onestep := round(WAIC_onestep - min(WAIC_onestep), 1)]; TAB[, dWAIC_zone := round(WAIC_zone - min(WAIC_zone), 1)]
## SE of differences against K = 3 (paired)
pd <- function(La, Lb) { f <- function(L) { m <- apply(L, 2, max); lppd <- m + log(colMeans(exp(sweep(L, 2, m)))); -2 * (lppd - apply(L, 2, var)) }
  d <- f(La) - f(Lb); c(diff = sum(d), se = sqrt(length(d) * var(d))) }
TAB[, diff_vs_K3_onestep := sapply(K, function(k) { v <- pd(PW[[paste0("K", k)]], PW[["K3"]]); sprintf("%.1f (SE %.1f)", v[1], v[2]) })]
if (requireNamespace("loo", quietly = TRUE)) {
  TAB[, LOOIC_zone := sapply(K, function(k) round(loo::loo(LZ[[paste0("K", k)]])$estimates["looic", 1], 1))]
  TAB[, pareto_k_gt_0.7 := sapply(K, function(k) sum(loo::pareto_k_values(loo::loo(LZ[[paste0("K", k)]])) > 0.7))]
}
fwrite(TAB, file.path(OUT, "Table_S_number_of_states.csv"))
print(TAB)

## ---- decision rule and sentences
b <- TAB[which.min(WAIC_onestep)]
txt <- c(sprintf("Models with 2, 3 and 4 hidden states were compared (%d chains of %s iterations each, burn-in %s, thinning %d; %d further chains with 2024 hidden).", NCH_FULL, format(BURN + NIT, big.mark = ","), format(BURN, big.mark = ","), THIN, NCH_HOLD),
  sprintf("K = %d: one-step-ahead WAIC %s (difference from K = 3 %s), hold-out CRPS %.2f, 95%% coverage %.1f%%, %d of 5 predictive statistics inside their intervals, residual SD %.2f, smallest state %.1f%% of zone-months, %d zones with fewer than 3 months in the top state, maximum R-hat %.3f.",
          TAB$K, format(TAB$WAIC_onestep, big.mark = ","), TAB$diff_vs_K3_onestep, TAB$holdout_CRPS, TAB$holdout_cov95, TAB$ppc_inside_of_5, TAB$rqr_sd,
          100 * TAB$smallest_state_share, TAB$zones_lt3_top_months, TAB$max_rhat),
  "Reading: a lower WAIC or CRPS is better. A difference smaller than about twice its SE is not a clear difference; then the model with fewer states is preferred.",
  "A state that holds very few zone-months, or a top state that most zones never reach, is weakly identified and is a reason not to add that state.")
writeLines(txt, file.path(OUT, "States_sentences.txt")); cat("\n", paste(txt, collapse = "\n"), "\n")
say("Done. Send: cat", file.path(OUT, "States_sentences.txt"), "and", file.path(OUT, "Table_S_number_of_states.csv"))
