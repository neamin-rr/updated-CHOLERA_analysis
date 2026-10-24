DEF <- list(cholera.tau23 = c(2, 2), cholera.anchhi = NA, cholera.hiermean = TRUE, cholera.cores = 4)   # defaults for this run; any option already set before sourcing (for example by Run_all_variants.R) is kept
for (nm in names(DEF)) if (is.null(getOption(nm))) do.call(options, setNames(list(DEF[[nm]]), nm))
# ==============================================================================
# HMM_check_small.R : one self-contained script. Small run (about 5-10 minutes) that answers: is the HMM good?
# Paste the whole thing into R (or source it). Needs the Stage 1 .rds, JAGS, rjags, coda, data.table.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(rjags); library(coda); library(Matrix) })
set.seed(2026)
proj <- getOption("cholera.proj", "."); RDS <- "INLA_Rebuilt_From_Zero_OCCPOP.rds"
out  <- file.path(proj, "improve_hmm"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
ANCHHI <- getOption("cholera.anchhi", NA)   # centre of the prior for the HIGH-state mean. NA = 95th percentile of positive counts (Methods, about 389). Set a number (e.g. 170) to centre it on the observed mean of high months; the 95th percentile is still used as the state threshold
TAU23 <- getOption("cholera.tau23", c(10, 3))   # prior on the precision of the log gap moderate -> high state mean. c(10, 3) = tight (sd 0.55); c(2, 2) = wide. Set with options(cholera.tau23 = c(2, 2))
NCH <- getOption("cholera.nch", 4); NITER <- ceiling(12000 / NCH); BURN <- 2000; THIN <- 6; ADAPT <- 1000   # more chains = more cores used; the total number of kept draws stays the same (12000 iterations shared over the chains), each chain keeps its own burn-in   # JAGS adaptation iterations before burn-in (the smoke test uses only 20, so it finishes in a couple of minutes)
if (isTRUE(getOption("cholera.smoke", FALSE))) { NITER <- 150; BURN <- 100; THIN <- 1; ADAPT <- 20 }   # SMOKE TEST (about 1 minute): only checks the code runs end to end. Results are meaningless; do not read the verdict.          # 4 chains as in the main analysis; 4 x 500 = 2,000 draws (same size as before)          # small on purpose; raise NITER to 10000 for a better look
FINAL <- isTRUE(getOption("cholera.final", FALSE)); FTEST <- isTRUE(getOption("cholera.finaltest", FALSE)); SEG <- 5000L   # FINAL: the long run. 4 chains, 300,000 iterations each = 200,000 burn-in + 100,000 kept, thinning 50 -> 2,000 draws per chain, 8,000 in total. Each chain is saved every SEG iterations and can be resumed after a stop.
if (FINAL) { NCH <- 4; NITER <- 100000; BURN <- 200000; THIN <- 50; ADAPT <- 1000
  if (FTEST) { NITER <- 2000; BURN <- 2000; SEG <- 500L; ADAPT <- 100 } }   # FTEST: a 15-minute end-to-end test of the long-run machinery (results meaningless)

# ---- data (same as the main script) --------------------------------------------------------------
S0 <- readRDS(file.path(proj, RDS)); dat <- as.data.table(S0$data); setorder(dat, ID_area, ID_time)
zone_order <- S0$zone_order; Z <- length(zone_order); T <- max(dat$ID_time); stopifnot(nrow(dat) == Z * T)
Y   <- matrix(dat$O, Z, T, byrow = TRUE)
Em  <- matrix(dat$E, Z, T, byrow = TRUE); OFF <- Em / mean(Em)                  # JUAN 2: population enters through the expected count
popz <- tapply(dat$monthly_population, dat$ID_area, mean)
bym <- function(f) { r <- f$summary.random$ID_area_bym; list(m = r$mean[1:Z], s = r$sd[1:Z]) }
fo <- bym(S0$fits$occ); fi <- bym(S0$fits$int)
d_occ <- fo$m - mean(fo$m); d_int <- fi$m - mean(fi$m); sd_occ <- pmax(fo$s, 0.05); sd_int <- pmax(fi$s, 0.05)
cat(sprintf("Zone multipliers before centring: mean exp(d_int) %.3f | mean exp(d_occ) %.3f\n", mean(exp(d_int)), mean(exp(d_occ))))
d_occ <- d_occ - log(mean(exp(d_occ))); d_int <- d_int - log(mean(exp(d_int)))   # centre so the average zone multiplier is exactly 1 (otherwise every state mean is inflated by mean(exp(delta)))
thr <- as.numeric(quantile(Y[Y > 0], 0.95)); anch3 <- thr; anch2 <- median(Y[Y > 0 & Y < thr]); anch1 <- 0.01
st_emp <- matrix(1L, Z, T); st_emp[Y > 0 & Y < thr] <- 2L; st_emp[Y >= thr] <- 3L
trm <- matrix(0, 3, 3); for (z in 1:Z) for (t in 2:T) { a <- st_emp[z, t - 1]; b <- st_emp[z, t]; trm[a, b] <- trm[a, b] + 1 }
self_emp <- diag(trm) / rowSums(trm); CONC <- 10
a_self <- pmax(self_emp, 0.05) * CONC; a_self[3] <- a_self[2]
a_oth  <- (1 - c(self_emp[1], self_emp[2], self_emp[2])) * CONC / 2
alpha <- matrix(0, 3, 3); for (i in 1:3) { alpha[i, ] <- a_oth[i]; alpha[i, i] <- a_self[i] }
cat(sprintf("Zones %d | months %d | cases %s | zeros %.1f%% | min E %.3g | anchors %.0f / %.0f\n", Z, T, format(sum(Y), big.mark = ","), 100 * mean(Y == 0), min(dat$E), anch2, anch3))

# ---- OPTIONAL EXTENSIONS (all OFF by default; each is switched on with an option, so every concern can be tested one by one) ----------
EXT_POP  <- isTRUE(getOption("cholera.pop", FALSE))      # zone population as a covariate on the expected count (the offset may not remove size completely)
EXT_SEAS <- isTRUE(getOption("cholera.season", FALSE))   # a common calendar-phase effect (12 values, centred) on the expected count
EXT_NB   <- isTRUE(getOption("cholera.nb", FALSE))       # association with the neighbouring zones' counts in the PREVIOUS month (an association, not a mechanism)
KAPPA    <- getOption("cholera.kappa", NA)               # NA = each zone has its own transition rows (as in the Methods). A number (e.g. 30) = partial pooling: zone rows are drawn around a national row with that strength
HIERMEAN <- isTRUE(getOption("cholera.hiermean", FALSE))   # TRUE = the national high-state level is LEARNED from all zones (hyper-mean with a wide prior) instead of being fixed at the anchor; zones with few high months then borrow from the others
RMIN     <- getOption("cholera.rmin", 0.05)              # lower bound on the NB size r. 0.05 allows extreme over-dispersion; raising it (0.5) tightens the variance
ph_t <- ((1:T) - 1) %% 12 + 1                            # calendar phase, assuming the time index starts in January
zpop <- as.numeric(scale(log(as.numeric(popz))))         # standardised log zone population
Qg <- tryCatch(S0$graph_Q, error = function(e) NULL); NBm <- matrix(0, Z, T)
if (!is.null(Qg) && nrow(Qg) == Z) {
  Nb <- (abs(as.matrix(Qg)) > 0) * 1; diag(Nb) <- 0
  vals <- log1p((Nb %*% Y)[, -T] / pmax(rowSums(Nb), 1))                       # log(1 + mean count of the neighbours in the previous month)
  NBm <- cbind(0, (vals - mean(vals)) / sd(vals))                              # standardised; month 1 has no previous month (set to 0)
} else if (EXT_NB) { cat("Neighbour term switched off: graph_Q not found in the Stage 1 file.\n"); EXT_NB <- FALSE }
EXT <- paste0("", if (EXT_POP) "_pop", if (EXT_SEAS) "_seas", if (EXT_NB) "_nb", if (!is.na(KAPPA)) paste0("_k", KAPPA), if (RMIN != 0.05) paste0("_rmin", RMIN), if (HIERMEAN) "_hm")
TAGSTR <- sprintf("tau%s_%s_hi%s%s", TAU23[1], TAU23[2], ifelse(is.na(ANCHHI), "389", ANCHHI), EXT)
TAG_FIT <- TAGSTR; if (FINAL) TAGSTR <- paste0(TAGSTR, if (FTEST) "_FINALTEST" else "_FINAL")   # saved fits use TAG_FIT (shared with the short checks); output tables use TAGSTR
cat(sprintf("Model version: %s\n", if (EXT == "") "as in the Methods (no extensions)" else paste("extensions", EXT)))

# ---- the HMM (negative binomial emission with the population offset) -------------------------------
# The model text is built from the options: terms that are switched off are NOT put into the model at all.
# (Switched-off terms that stay in the graph still cost time on every iteration, because they touch all 11,040 zone-months.)
mult_terms <- paste0("", if (EXT_POP) " + b_pop * zpop[z]", if (EXT_SEAS) " + seas[ph[t]]", if (EXT_NB) " + b_nb * NB[z, t]")
mult_expr  <- if (nzchar(mult_terms)) paste0(" * exp(0", mult_terms, ")") else ""
trans_txt  <- if (!is.na(KAPPA)) "for (j in 1:3) { alphaZ[z, i, j] <- max(KAP * gam_nat[i, j], 0.01) }\n      gam[z, i, 1:3] ~ ddirch(alphaZ[z, i, 1:3])" else "gam[z, i, 1:3] ~ ddirch(alpha[i, 1:3])"
model_txt <- paste0("
model {
  pi0[1:3] ~ ddirch(c(8, 1, 1))
", if (!is.na(KAPPA)) "  for (i in 1:3) { gam_nat[i, 1:3] ~ ddirch(alpha[i, 1:3]) }\n", if (EXT_POP) "  b_pop ~ dnorm(0, 1)\n", if (EXT_NB) "  b_nb ~ dnorm(0, 1)\n",
if (EXT_SEAS) "  for (k in 1:12) { sraw[k] ~ dnorm(0, 4) }\n  for (k in 1:12) { seas[k] <- sraw[k] - mean(sraw[1:12]) }\n", "
  eta0 <- log(anch3 - anch2)
  eta23 ~ dnorm(eta0, prec_eta)
  tau23 ~ dgamma(tau_a, tau_b)
  sd_g12 ~ dunif(1, 2 * (anch2 - anch1))
  tau_g12 <- pow(sd_g12, -2)
  for (z in 1:Z) {
    mu1[z] ~ dnorm(anch1, 1e4)T(1e-4, )
    g12[z] ~ dnorm(anch2 - anch1, tau_g12)T(0, )
    lg23[z] ~ dnorm(eta23, tau23)
    mu[z, 1] <- mu1[z]
    mu[z, 2] <- mu1[z] + g12[z]
    mu[z, 3] <- mu[z, 2] + exp(lg23[z])
    delta[z, 1] ~ dnorm(d_occ[z], prec_occ[z])
    delta[z, 2] ~ dnorm(d_int[z], prec_int[z])
    delta[z, 3] ~ dnorm(d_int[z], prec_int[z])
    for (s in 1:3) { lam_s[z, s] <- mu[z, s] * exp(delta[z, s]) }
    r[z] ~ dgamma(2, 0.5)T(rmin, )
    for (i in 1:3) {
      ", trans_txt, "
    }
    S[z, 1] ~ dcat(pi0[1:3])
    for (t in 2:T) { S[z, t] ~ dcat(gam[z, S[z, t - 1], 1:3]) }
    for (t in 1:T) {
      lam[z, t] <- max(lam_s[z, S[z, t]] * OFF[z, t]", mult_expr, ", 1e-6)
      p[z, t] <- r[z] / (r[z] + lam[z, t])
      Y[z, t] ~ dnegbin(p[z, t], r[z])
    }
  }
}")
mf <- file.path(out, sprintf("hmm_check_model_%s.jags", TAGSTR)); writeLines(model_txt, mf)   # one model file per variant, so parallel runs never overwrite each other
mk_data <- function(OFFm) c(list(Y = Y, OFF = OFFm, Z = Z, T = T, alpha = alpha, d_occ = d_occ, d_int = d_int,
                               prec_occ = 1 / sd_occ^2, prec_int = 1 / sd_int^2, anch1 = anch1, anch2 = anch2, anch3 = ifelse(is.na(ANCHHI), anch3, ANCHHI), tau_a = TAU23[1], tau_b = TAU23[2],
                               prec_eta = ifelse(HIERMEAN, 1, 1e6), rmin = RMIN),
                               if (!is.na(KAPPA)) list(KAP = KAPPA), if (EXT_POP) list(zpop = zpop), if (EXT_SEAS) list(ph = ph_t), if (EXT_NB) list(NB = NBm))
mk_init <- function(ch) list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 1000 + ch, S = st_emp)
library(parallel)
NCORES <- max(1, min(getOption("cholera.cores", 4), NCH))   # one core per chain: 4 chains can use at most 4 cores, so 64 cores gives no speed-up here (it only starts idle workers)
if (exists("cl_hmm")) try(stopCluster(cl_hmm), silent = TRUE)
cl_hmm <- makeCluster(NCORES, port = 20000L + as.integer(Sys.getpid()) %% 20000L)   # own port per run: runs started together would otherwise all pick the same random port (the script sets a fixed seed)
run_chain <- function(a) {                                   # one chain on one core
  if (!is.null(a$ckpt) && file.exists(a$ckpt)) { old <- tryCatch(readRDS(a$ckpt), error = function(e) NULL); if (!is.null(old)) return(old) }   # this chain was finished before a stop: reuse it
  suppressPackageStartupMessages({ library(rjags); library(coda) })
  if (is.null(a$seg)) {                                      # short runs: one go
    m <- jags.model(a$mf, data = a$data, inits = a$init, n.chains = 1, quiet = TRUE, n.adapt = a$adapt); update(m, a$burn)
    r <- coda.samples(m, a$mon, n.iter = a$n, thin = a$thin)[[1]]
    if (!is.null(a$ckpt)) { saveRDS(r, paste0(a$ckpt, ".tmp")); file.rename(paste0(a$ckpt, ".tmp"), a$ckpt) }   # each chain is saved the moment it finishes
    return(r) }
  # ---- the long run: done in blocks of a$seg iterations; after each block the state of the chain and the kept draws are saved, so a stop loses at most one block ----
  say <- function(...) cat(format(Sys.time(), "%d %b %H:%M:%S"), ..., "\n", file = paste0(a$ckpt, "_progress.txt"), append = TRUE)
  sf <- paste0(a$ckpt, "_state.rds"); bf <- function(i) sprintf("%s_block%03d.rds", a$ckpt, i)
  st <- if (file.exists(sf)) tryCatch(readRDS(sf), error = function(e) NULL) else NULL
  done <- if (is.null(st)) 0L else st$done; nblk <- if (is.null(st)) 0L else st$nblk
  init <- if (is.null(st) || is.null(st$state)) a$init else st$state
  m <- jags.model(a$mf, data = a$data, inits = init, n.chains = 1, quiet = TRUE, n.adapt = if (done == 0) a$adapt else min(a$adapt, 500))
  if (done > 0) say("resumed from iteration", done)
  total <- a$burn + a$n
  say("iteration", done, "of", total)   # shows at once on the dashboard; the next line comes after the first block of iterations
  while (done < total) {
    k <- min(a$seg, if (done < a$burn) a$burn - done else total - done)
    if (done < a$burn) update(m, k) else {
      blk <- coda.samples(m, a$mon, n.iter = k, thin = a$thin)[[1]]; nblk <- nblk + 1L
      saveRDS(blk, paste0(bf(nblk), ".tmp")); file.rename(paste0(bf(nblk), ".tmp"), bf(nblk)) }
    done <- done + k
    stt <- tryCatch(m$state(internal = TRUE)[[1]], error = function(e) NULL)   # the full state of the chain, including its random-number generator
    saveRDS(list(done = done, nblk = nblk, state = stt), paste0(sf, ".tmp")); file.rename(paste0(sf, ".tmp"), sf)
    say("iteration", done, "of", total)
  }
  r <- coda::mcmc(do.call(rbind, lapply(seq_len(nblk), function(i) as.matrix(readRDS(bf(i))))), start = a$burn + a$thin, thin = a$thin)
  saveRDS(r, paste0(a$ckpt, ".tmp")); file.rename(paste0(a$ckpt, ".tmp"), a$ckpt)
  unlink(c(sf, vapply(seq_len(nblk), bf, "")))
  r
}
fit <- function(OFFm, tag, key) {                           # each finished fit is saved, so a crash or power cut does not lose it: starting the script again reloads it
  cf <- file.path(out, sprintf("cache_%s_%s_n%d_b%d.rds", TAG_FIT, key, NITER, BURN))
  if (file.exists(cf)) { cached <- tryCatch(readRDS(cf), error = function(e) NULL)
    if (!is.null(cached)) { cat("Loaded the saved fit:", basename(cf), "\n"); return(cached) }
    cat("The saved fit", basename(cf), "is damaged (the machine probably stopped while it was being written); fitting again\n"); unlink(cf) }
  cat("Fitting:", tag, "-", NCH, "chains on", NCORES, "cores\n")
  tasks <- lapply(1:NCH, function(ch) list(mf = mf, data = mk_data(OFFm), init = mk_init(ch), burn = BURN, adapt = ADAPT, n = NITER, thin = THIN, seg = if (FINAL) SEG else NULL, ckpt = file.path(out, sprintf("chain_%s_%s_n%d_b%d_ch%d.rds", TAG_FIT, key, NITER, BURN, ch)), mon = c("gam", "mu", "delta", "r", "S", "eta23", if (EXT_POP) "b_pop", if (EXT_NB) "b_nb", if (EXT_SEAS) "seas")))
  res_fit <- as.mcmc.list(parLapplyLB(cl_hmm, tasks, run_chain)); saveRDS(res_fit, paste0(cf, ".tmp")); file.rename(paste0(cf, ".tmp"), cf); unlink(list.files(out, pattern = sprintf("^chain_%s_%s_n%d_b%d_ch", TAG_FIT, key, NITER, BURN), full.names = TRUE)); res_fit   # written under a temporary name and renamed, so a stop in the middle never leaves a half-written file
}
mc  <- fit(OFF, "with the population offset (main)", "main")
mc0 <- if (FINAL) { f0 <- file.path(out, sprintf("cache_%s_nooff_n1091_b2000.rds", TAG_FIT))
  if (!file.exists(f0)) stop("Final run: the saved short-run fit without offset (", basename(f0), ") is missing. Run the 11-chain check for this version first.")
  cat("Offset comparison uses the saved short-run fit without offset:", basename(f0), "\n"); readRDS(f0)
} else fit(matrix(1, Z, T), "without offset: absolute counts, as in the old code", "nooff")

# ---- helper: pull draws ---------------------------------------------------------------------------
X  <- as.matrix(mc); nd <- nrow(X)
gam_self <- sapply(1:3, function(i) rowMeans(X[, sprintf("gam[%d,%d,%d]", 1:Z, i, i)]))          # draws x 3, national mean over zones
stn <- c("Baseline", "Moderate", "High-intensity")
trm3 <- rbindlist(lapply(1:3, function(i) rbindlist(lapply(1:3, function(j) { v <- rowMeans(X[, sprintf("gam[%d,%d,%d]", 1:Z, i, j)])    # mean over zones, one value per draw
  data.table(from = stn[i], to = stn[j], mean = round(mean(v), 3), lo = round(quantile(v, .025), 3), hi = round(quantile(v, .975), 3)) }))))
cat("\nTable 4: transition probabilities (mean over zones; 95% interval over draws)\n"); print(dcast(trm3[, .(from, to, txt = sprintf("%.3f [%.3f, %.3f]", mean, lo, hi))], from ~ to, value.var = "txt")[match(stn, from)], row.names = FALSE)
fwrite(trm3, file.path(out, sprintf("Table_4_transition_matrix_%s.csv", TAGSTR)))
Smat <- function(x) array(x[, sprintf("S[%d,%d]", rep(1:Z, T), rep(1:T, each = Z))], c(nrow(x), Z, T))
EFFM <- function(k) { e <- matrix(0, Z, T)                                                         # multiplier from the optional terms, one draw (1 everywhere when all are off)
  if (EXT_POP) e <- e + X[k, "b_pop"] * zpop
  if (EXT_SEAS) e <- e + matrix(X[k, sprintf("seas[%d]", ph_t)], Z, T, byrow = TRUE)
  if (EXT_NB) e <- e + X[k, "b_nb"] * NBm
  exp(e) }
if (EXT != "") { cat("\nExtension terms (posterior mean and 95% interval):\n")
  if (EXT_POP) cat(sprintf("  population (per SD of log population): %.3f (%.3f to %.3f)\n", mean(X[, "b_pop"]), quantile(X[, "b_pop"], .025), quantile(X[, "b_pop"], .975)))
  if (EXT_NB)  cat(sprintf("  neighbours' previous-month counts (per SD): %.3f (%.3f to %.3f)\n", mean(X[, "b_nb"]), quantile(X[, "b_nb"], .025), quantile(X[, "b_nb"], .975)))
  if (EXT_SEAS) { sm <- colMeans(X[, sprintf("seas[%d]", 1:12)]); cat(sprintf("  calendar phase effects (log scale) range %.3f to %.3f\n", min(sm), max(sm))) }
  if (HIERMEAN) { lv <- anch2 + exp(X[, "eta23"]); cat(sprintf("  national high-state level learned from the data: %.0f (%.0f to %.0f); the anchor was %.0f\n", mean(lv), quantile(lv, .025), quantile(lv, .975), ifelse(is.na(ANCHHI), anch3, ANCHHI))) }
  if (!is.na(KAPPA)) cat(sprintf("  transition rows pooled toward the national row with strength %s\n", KAPPA)) }
burden <- function(x) { A <- Smat(x); apply(A >= 2, 2, mean) }                                    # share of months in moderate/high, per zone

# ---- CHECKS -------------------------------------------------------------------------------------
res <- list(); add <- function(q, ok, detail) { res[[length(res) + 1]] <<- data.table(question = q, ok = ok, detail = detail); cat(sprintf("[%s] %s  %s\n", if (is.na(ok)) "SKIP" else if (ok) "YES" else "NO ", q, detail)) }
cat("\n============ IS THE HMM GOOD? ============\n")
# 1 convergence (small run: threshold 1.1 here; the FULL run uses 1.01)
rh <- gelman.diag(mc[, sprintf("gam[%d,%d,%d]", 1:Z, 1, 1)], multivariate = FALSE)$psrf[, 1]
add("1. Chains agree (R-hat < 1.1 for baseline self-transition in at least 90% of zones)", mean(rh < 1.1, na.rm = TRUE) >= 0.9, sprintf("%.0f%% of zones, max R-hat %.2f%s", 100 * mean(rh < 1.1, na.rm = TRUE), max(rh, na.rm = TRUE), if (FINAL) "" else " (small run)"))
# 2 mechanics
ordok <- mean(sapply(1:Z, function(z) mean(X[, sprintf("mu[%d,1]", z)] < X[, sprintf("mu[%d,2]", z)] & X[, sprintf("mu[%d,2]", z)] < X[, sprintf("mu[%d,3]", z)])))
add("2. States ordered (baseline < moderate < high) in every draw", ordok == 1, sprintf("%.3f of draws", ordok))
rows <- X[, sprintf("gam[%d,%d,%d]", 1, 1, 1:3)]; add("3. Transition rows sum to 1", all(abs(rowSums(rows) - 1) < 1e-6), "zone 1, baseline row")
# 3 expected stay (JUAN 5)
cat("\nExpected stay in each state = 1 / (1 - p_ii), national mean over zones:\n")
for (i in 1:3) { d <- 1 / (1 - gam_self[, i]); cat(sprintf("  state %d: p_ii %.3f | %.1f months (95%% interval %.1f to %.1f)\n", i, mean(gam_self[, i]), mean(d), quantile(d, .025), quantile(d, .975))) }
add("4. High state is neither absorbing nor fleeting (stay 1.5 to 8 months)", mean(1 / (1 - gam_self[, 3])) > 1.5 & mean(1 / (1 - gam_self[, 3])) < 8, sprintf("%.1f months", mean(1 / (1 - gam_self[, 3]))))
# 4 posterior predictive check
thrk <- anch3; nsim <- min(100, nd); keep <- sample(nd, nsim)
obs <- c(zero = mean(Y == 0), total = sum(Y), high = mean(Y >= thrk), max = max(Y), meanpos = mean(Y[Y > 0]))
sim <- t(sapply(keep, function(k) { A <- Smat(X[k, , drop = FALSE])[1, , ]
  lam_s <- matrix(X[k, sprintf("mu[%d,%d]", rep(1:Z, 3), rep(1:3, each = Z))] * exp(X[k, sprintf("delta[%d,%d]", rep(1:Z, 3), rep(1:3, each = Z))]), Z, 3)
  rz <- X[k, sprintf("r[%d]", 1:Z)]; lam <- matrix(lam_s[cbind(rep(1:Z, T), as.vector(A))], Z, T) * OFF * EFFM(k)
  ys <- matrix(rnbinom(Z * T, size = rep(rz, T), mu = as.vector(lam)), Z, T)
  c(zero = mean(ys == 0), total = sum(ys), high = mean(ys >= thrk), max = max(ys), meanpos = mean(ys[ys > 0])) }))
inside <- sapply(names(obs), function(n) obs[n] >= quantile(sim[, n], .025) & obs[n] <= quantile(sim[, n], .975))
cat("\nPosterior predictive check (observed vs simulated 95% interval):\n"); print(round(rbind(observed = obs, lo = apply(sim, 2, quantile, .025), hi = apply(sim, 2, quantile, .975)), 3))
add("5. Model reproduces the data (at least 4 of 5 statistics inside the 95% interval)", sum(inside) >= 4, sprintf("%d of 5 inside", sum(inside)))
ppp <- sapply(names(obs), function(n) mean(sim[, n] >= obs[n]))
cat("Predictive p-values (0.5 is ideal; below 0.05 or above 0.95 = misfit):\n"); print(round(ppp, 3))
add("5b. Predictive p-values for zeros, total, high-count share and mean positive count all between 0.05 and 0.95", all(ppp[c("zero", "total", "high", "meanpos")] > 0.05 & ppp[c("zero", "total", "high", "meanpos")] < 0.95), paste(names(ppp), round(ppp, 2), collapse = " | "))
# 5 population (JUAN 2): do states follow population or the cholera pattern?
b1 <- burden(X); b0 <- burden(as.matrix(mc0))
cat(sprintf("\nPopulation check: Spearman(zone burden, zone population) with offset %.2f | without offset %.2f\n", cor(b1, popz, method = "spearman"), cor(b0, popz, method = "spearman")))
cat(sprintf("Spearman(zone burden with offset, without offset) = %.2f\n", cor(b1, b0, method = "spearman")))
add("6. Population offset matters little for the zone ranking (Spearman with vs without >= 0.8)", cor(b1, b0, method = "spearman") >= 0.8, sprintf("%.2f", cor(b1, b0, method = "spearman")))
add("7. States do not just copy zone population (|Spearman| with population < 0.5 with offset)", abs(cor(b1, popz, method = "spearman")) < 0.5, sprintf("%.2f", cor(b1, popz, method = "spearman")))
# 6 uncertainty (JUAN 3): does each zone keep its burden group across posterior draws?
Aall <- Smat(X); cut_hi <- quantile(b1, .75); grp <- function(k) apply(Aall[k, , ] >= 2, 1, mean) >= cut_hi
G <- sapply(sample(nd, min(100, nd)), grp); stab <- rowMeans(G); stab <- pmax(stab, 1 - stab)
add("8. Zones keep their high-burden / not-high group across posterior draws (>= 70% of zones stable in at least 70% of draws)", mean(stab >= 0.7) >= 0.7, sprintf("%.0f%% of zones", 100 * mean(stab >= 0.7)))
# ---- HOLD-OUT PREDICTION TEST: how big is the prediction error, and how much smaller than simple methods? -----------------------
# The last 12 months are hidden from the model (set to missing), the HMM is refitted, and its predictions for those 12 months are
# compared with the real counts and with simple baselines that use only the earlier 108 months. This is the fair test, because the
# in-sample fit above has already seen the data it is explaining.
if (isTRUE(getOption("cholera.holdout", TRUE))) {
  HOLD <- 12; hold <- (T - HOLD + 1):T; train <- setdiff(1:T, hold)
  cat(sprintf("\n============ HOLD-OUT PREDICTION: refit with the last %d months hidden (adds about 50%% to the run time) ============\n", HOLD))
  dh <- mk_data(OFF); dh$Y[, hold] <- NA; if (EXT_NB) dh$NB[, hold[-1]] <- 0   # months after the first hidden month would use hidden counts, so the neighbour term is set to 0 there                                      # the model never sees these months
  ynodes <- sprintf("Y[%d,%d]", rep(1:Z, HOLD), rep(hold, each = Z))          # predictions for the hidden months (zone varies fastest)
  tasks <- lapply(1:NCH, function(ch) list(mf = mf, data = dh, init = mk_init(ch), burn = BURN, adapt = ADAPT, n = NITER, thin = THIN, ckpt = file.path(out, sprintf("chain_%s_holdout_n%d_b%d_ch%d.rds", TAGSTR, NITER, BURN, ch)), mon = ynodes))
  cf_h <- file.path(out, sprintf("cache_%s_holdout_n%d_b%d.rds", TAGSTR, NITER, BURN))
  Yp <- if (file.exists(cf_h)) tryCatch(readRDS(cf_h), error = function(e) NULL) else NULL
  if (is.null(Yp)) { Yp <- as.matrix(as.mcmc.list(parLapplyLB(cl_hmm, tasks, run_chain))); saveRDS(Yp, paste0(cf_h, ".tmp")); file.rename(paste0(cf_h, ".tmp"), cf_h) } else cat("Loaded the saved hold-out fit\n")        # draws x (Z*HOLD) predicted counts
  yobs <- as.vector(Y[, hold])                                                # same order as ynodes
  # baselines, all built from the 108 training months only
  ph_all <- ((1:T) - 1) %% 12 + 1
  b_zero  <- rep(0, Z * HOLD)                                                  # always predict no cases (hard to beat when 89.5% are zeros)
  b_nat   <- as.vector(sum(Y[, train]) / sum(OFF[, train]) * OFF[, hold])      # same national rate everywhere, scaled by expected count
  b_zone  <- rep(rowMeans(Y[, train]), HOLD)                                   # each zone's own average, no time pattern
  b_zph   <- as.vector(sapply(hold, function(t) rowMeans(Y[, train[ph_all[train] == ph_all[t]], drop = FALSE])))   # zone average for the same calendar phase
  b_last  <- rep(Y[, min(hold) - 1], HOLD)                                     # carry the last observed month forward
  crps_cell <- function(x, y) { n <- length(x); xs <- sort(x); mean(abs(x - y)) - sum((2 * seq_len(n) - n - 1) * xs) / n^2 }   # CRPS from draws; equals the absolute error for a single number
  crps_m <- mean(sapply(seq_along(yobs), function(j) crps_cell(Yp[, j], yobs[j])))
  pm_mean <- colMeans(Yp); pm_med <- apply(Yp, 2, median)
  bl <- list("always zero" = b_zero, "national rate x expected" = b_nat, "zone average" = b_zone, "zone average, same calendar phase" = b_zph, "last observed month" = b_last)
  tab <- data.table(method = c("HMM (hold-out)", names(bl)),
                    MAE   = c(mean(abs(yobs - pm_med)), sapply(bl, function(b) mean(abs(yobs - b)))),
                    RMSE  = c(sqrt(mean((yobs - pm_mean)^2)), sapply(bl, function(b) sqrt(mean((yobs - b)^2)))),
                    CRPS  = c(crps_m, sapply(bl, function(b) mean(abs(yobs - b)))))
  tab[, CRPS_reduction_pct := round(100 * (1 - tab$CRPS[1] / CRPS), 1)]; tab[1, CRPS_reduction_pct := NA]
  tab[, `:=`(MAE = round(MAE, 3), RMSE = round(RMSE, 2), CRPS = round(CRPS, 3))]
  cat("Prediction error for the hidden months (smaller is better). CRPS_reduction_pct = how much smaller the HMM's CRPS is than that method's (positive = HMM better):\n"); print(tab)
  lo <- apply(Yp, 2, quantile, .025); hi <- apply(Yp, 2, quantile, .975); cover <- mean(yobs >= lo & yobs <= hi)
  tot_obs <- sum(yobs); tot_pred <- apply(matrix(Yp, nrow(Yp)), 1, sum)
  cat(sprintf("Coverage of the 95%% prediction interval: %.1f%% of hidden zone-months | hidden-month total: observed %s, predicted %s (95%% interval %s to %s)\n",
              100 * cover, format(tot_obs, big.mark = ","), format(round(mean(tot_pred)), big.mark = ","), format(round(quantile(tot_pred, .025)), big.mark = ","), format(round(quantile(tot_pred, .975)), big.mark = ",")))
  fwrite(tab, file.path(out, sprintf("Table_S_HMM_holdout_prediction_%s.csv", TAGSTR)))
  best_simple <- min(tab$CRPS[tab$method %in% c("national rate x expected", "zone average", "zone average, same calendar phase", "last observed month")])
  add("9. Hold-out: the HMM predicts the hidden 12 months better than every simple baseline (CRPS)", crps_m < best_simple, sprintf("HMM %.3f vs best simple %.3f (%.1f%% smaller)", crps_m, best_simple, 100 * (1 - crps_m / best_simple)))
  add("10. Hold-out: the 95% prediction interval covers 90 to 99% of hidden zone-months", cover >= 0.90 & cover <= 0.99, sprintf("%.1f%%", 100 * cover))
  add("11. Hold-out: predicted total for the hidden months is within the 95% interval", tot_obs >= quantile(tot_pred, .025) & tot_obs <= quantile(tot_pred, .975), sprintf("observed %s", format(tot_obs, big.mark = ",")))
} else cat("\nHold-out prediction test skipped (options(cholera.holdout = FALSE)).\n")
# ---- THE FOUR INDICATORS taken from the HMM (same definitions as Stage 2) ---------------------------------------------------------
# burden      = share of the 120 months a zone spends in the HIGH state
# onset       = probability of moving from the baseline state straight to the high state in a month (gam[z,1,3]); "velocity" in Stage 2
# persistence = probability of staying in the high state from one month to the next (gam[z,3,3])
# peak        = high-state mean level, mu_high x exp(delta) (cases per zone-month at the average expected count)
# These go into the LPA and Louvain steps, so the check is: can they be extracted, do they differ between zones by more than their uncertainty,
# and for how many zones are they driven by the prior (zones that almost never reach the high state have little information on onset or persistence).
cat("\n============ THE FOUR INDICATORS (extracted per zone, with uncertainty) ============\n")
A3 <- Smat(X); bur_d <- apply(A3 == 3, c(1, 2), mean); mx_d <- apply(A3, c(1, 2), max); rm(A3)   # mx_d = highest state each zone reaches, per draw                                    # draws x zones
on_d  <- X[, sprintf("gam[%d,%d,%d]", 1:Z, 1, 3)]; per_d <- X[, sprintf("gam[%d,%d,%d]", 1:Z, 3, 3)]
pk_d  <- X[, sprintf("mu[%d,3]", 1:Z)] * exp(X[, sprintf("delta[%d,3]", 1:Z)])
IND4 <- list(burden = bur_d, onset = on_d, persistence = per_d, peak = pk_d)
ind_tab <- rbindlist(lapply(names(IND4), function(n) { M <- IND4[[n]]; pm <- colMeans(M); psd <- apply(M, 2, sd)
  data.table(indicator = n, zone_median = signif(median(pm), 3), zone_min = signif(min(pm), 3), zone_max = signif(max(pm), 3),
             between_zone_sd = signif(sd(pm), 3), mean_within_zone_sd = signif(mean(psd), 3), ratio = round(sd(pm) / mean(psd), 2)) }))
print(ind_tab)
en_d <- X[, sprintf("gam[%d,%d,%d]", 1:Z, 1, 2)] + X[, sprintf("gam[%d,%d,%d]", 1:Z, 1, 3)]                 # ENTRY: chance of leaving the baseline state for moderate or high (far more events than the direct jump to high)
en_ratio <- sd(colMeans(en_d)) / mean(apply(en_d, 2, sd))
cat(sprintf("Alternative to onset: ENTRY into any active state (baseline to moderate or high): median %.4f, between/within ratio %.2f, Spearman with burden %.2f\n", median(colMeans(en_d)), en_ratio, cor(colMeans(en_d), colMeans(bur_d), method = "spearman")))
few_high <- colMeans(bur_d) * T < 3                                                              # zones with fewer than 3 expected high months
cat(sprintf("Zones with fewer than 3 expected high months: %d of %d. For these, onset and persistence rest mostly on the prior; read them with care.\n", sum(few_high), Z))
cp <- cor(sapply(IND4, colMeans), method = "spearman"); cat("Spearman between the four indicators (zone level):\n"); print(round(cp, 2))
# ---- the old-code correlation diagnostics (Step 5.1 of the first paper), applied to the new indicators before any clustering ----------
IM <- sapply(IND4, colMeans)                                                                      # zones x 4 posterior means
cat("\nPearson correlation between the four indicators (linear):\n"); pe <- cor(IM, method = "pearson"); print(round(pe, 3))
cat("Pearson p-values for each pair:\n")
for (pr in combn(colnames(IM), 2, simplify = FALSE)) { ct <- cor.test(IM[, pr[1]], IM[, pr[2]]); cat(sprintf("  %-12s vs %-12s r = %+.3f  p = %.4f  %s\n", pr[1], pr[2], ct$estimate, ct$p.value, ifelse(ct$p.value < 0.05, "(significant)", "(n.s.)"))) }
max_cor <- max(abs(pe[upper.tri(pe)]))
cat(sprintf("Highest absolute Pearson correlation: %.3f. %s\n", max_cor, if (max_cor > 0.85) "WARNING: above 0.85, two indicators are likely redundant." else if (max_cor > 0.70) "CAUTION: 0.70 to 0.85, they share substantial information." else "GOOD: all below 0.70, the indicators capture distinct dimensions."))
cond_num <- sqrt(max(eigen(cor(scale(IM)))$values) / min(eigen(cor(scale(IM)))$values))
cat(sprintf("Condition number: %.1f (%s)\n", cond_num, if (cond_num < 10) "GOOD, below 10" else if (cond_num < 30) "CAUTION, 10 to 30" else "WARNING, above 30"))
cat(sprintf("Ranges  burden %.4f to %.4f | onset %.4f to %.4f | persistence %.3f to %.3f (old code: about 0.4 to 0.8) | peak %.0f to %.0f (old code: about 350 to 4,000; now per average expected count, so lower)\n",
            min(IM[, 1]), max(IM[, 1]), min(IM[, 2]), max(IM[, 2]), min(IM[, 3]), max(IM[, 3]), min(IM[, 4]), max(IM[, 4])))
cat("Ten highest-burden zones (check they are recognisable):\n"); print(data.table(zone = zone_order, burden = round(IM[, 1], 4), onset = round(IM[, 2], 4), persistence = round(IM[, 3], 3), peak = round(IM[, 4], 0))[order(-burden)][1:10])
zone_ind <- data.table(zone = zone_order, burden = colMeans(bur_d), onset = colMeans(on_d), persistence = colMeans(per_d), peak = colMeans(pk_d),
  burden_lo = apply(bur_d, 2, quantile, .025), burden_hi = apply(bur_d, 2, quantile, .975), peak_lo = apply(pk_d, 2, quantile, .025), peak_hi = apply(pk_d, 2, quantile, .975), few_high_months = few_high)
fwrite(zone_ind, file.path(out, sprintf("Table_S_HMM_indicators_by_zone_%s.csv", TAGSTR)))
add("12. All four indicators can be extracted for every zone (finite, probabilities between 0 and 1)", all(sapply(IND4, function(M) all(is.finite(M)))) & all(bur_d >= 0 & bur_d <= 1) & all(on_d >= 0 & on_d <= 1) & all(per_d >= 0 & per_d <= 1), sprintf("%d zones x 4 indicators", Z))
add("13. At least 3 of 4 indicators differ between zones by more than their own uncertainty (between-zone SD / mean within-zone SD > 1)", sum(ind_tab$ratio > 1) >= 3, paste(ind_tab$indicator, ind_tab$ratio, collapse = " | "))
# PEAK REACHED: the level of the highest state the zone actually visits (per draw). Defined for every zone, so zones that rarely or never reach the high state
# are not forced onto a prior-driven high-state level: they get their moderate (or baseline) level, with its own uncertainty.
lvl <- function(s) X[, sprintf("mu[%d,%d]", 1:Z, s)] * exp(X[, sprintf("delta[%d,%d]", 1:Z, s)])
pr_d <- lvl(3) * (mx_d == 3) + lvl(2) * (mx_d == 2) + lvl(1) * (mx_d == 1)
pr_ratio <- sd(colMeans(pr_d)) / mean(apply(pr_d, 2, sd)); obs_max <- apply(Y, 1, max)
cat(sprintf("PEAK REACHED (level of the highest state visited): median %.0f, range %.1f to %.0f, between/within ratio %.2f (high-state level alone: %.2f)\n", median(colMeans(pr_d)), min(colMeans(pr_d)), max(colMeans(pr_d)), pr_ratio, ind_tab$ratio[ind_tab$indicator == "peak"]))
cat(sprintf("  Spearman with the observed highest monthly count (no model): %.2f | with burden: %.2f | share of zones that visit the high state in at least half of the draws: %.0f%%\n", cor(colMeans(pr_d), obs_max, method = "spearman"), cor(colMeans(pr_d), colMeans(bur_d), method = "spearman"), 100 * mean(colMeans(mx_d == 3) >= 0.5)))
zone_ind[, `:=`(entry = colMeans(en_d), peak_reached = colMeans(pr_d), peak_reached_lo = apply(pr_d, 2, quantile, .025), peak_reached_hi = apply(pr_d, 2, quantile, .975), obs_max_count = obs_max)]
fwrite(zone_ind, file.path(out, sprintf("Table_S_HMM_indicators_by_zone_%s.csv", TAGSTR)))   # rewritten with the new columns
add("17. Peak reached differs between zones by more than its uncertainty (ratio above 1) and follows the observed highest count (Spearman at least 0.7)", pr_ratio > 1 & cor(colMeans(pr_d), obs_max, method = "spearman") >= 0.7, sprintf("ratio %.2f | Spearman with observed max %.2f", pr_ratio, cor(colMeans(pr_d), obs_max, method = "spearman")))
add("16. The alternative onset measure (entry into any active state) differs between zones by more than its uncertainty (ratio above 1)", en_ratio > 1, sprintf("ratio %.2f", en_ratio))
add("14. Few zones lean on the prior for onset and persistence (fewer than 25% of zones with under 3 expected high months)", mean(few_high) < 0.25, sprintf("%.0f%% of zones", 100 * mean(few_high)))
add("15. Indicators are distinct (old-code rule: highest |Pearson| below 0.85 and condition number below 30)", max_cor < 0.85 & cond_num < 30, sprintf("max |r| %.2f | condition number %.1f", max_cor, cond_num))
# ---- CALIBRATION: is the high-state level right, and is the spread of the model right? (the two problems found in the first small run) -----------------
kc <- sample(nd, min(30, nd)); hi_r <- c(); rq_sd <- c(); pe_act <- c(); LAMc <- matrix(0, Z, T); RZc <- rep(0, Z)
act_m <- (Y > 0) | cbind(FALSE, Y[, -T] > 0)                      # zone-months with cases this month or the month before
for (k in kc) { A <- Smat(X[k, , drop = FALSE])[1, , ]
  ls <- matrix(X[k, sprintf("mu[%d,%d]", rep(1:Z, 3), rep(1:3, each = Z))] * exp(X[k, sprintf("delta[%d,%d]", rep(1:Z, 3), rep(1:3, each = Z))]), Z, 3)
  lam_k <- matrix(ls[cbind(rep(1:Z, T), as.vector(A))], Z, T) * OFF * EFFM(k); LAMc <- LAMc + lam_k / length(kc); RZc <- RZc + X[k, sprintf("r[%d]", 1:Z)] / length(kc)
  sel <- A == 3; if (any(sel) && mean(Y[sel]) > 0) hi_r <- c(hi_r, mean(lam_k[sel]) / mean(Y[sel]))
  rzk <- matrix(X[k, sprintf("r[%d]", 1:Z)], Z, T)               # this draw's own size r, same for every month of a zone
  p_hi <- pnbinom(Y, size = rzk, mu = lam_k); p_lo <- pnbinom(Y - 1, size = rzk, mu = lam_k)   # randomized quantile residual: exact for counts, about N(0,1) if the model is right, at any number of zeros
  zq <- qnorm(pmin(pmax(p_lo + runif(Z * T) * (p_hi - p_lo), 1e-12), 1 - 1e-12)); rq_sd <- c(rq_sd, sd(zq))
  pe_act <- c(pe_act, sd(((Y - lam_k) / sqrt(lam_k + lam_k^2 / rzk))[act_m])) }
ratio_hi <- mean(hi_r); res_sd_pearson <- sd(as.vector((Y - LAMc) / sqrt(LAMc + LAMc^2 / RZc))); res_sd <- mean(rq_sd); res_sd_act <- mean(pe_act)   # res_sd is now the randomized-quantile SD; the old Pearson SD over all cells is kept for comparison (it is pulled towards 0 by the 89.5% zero cells)
cat(sprintf("\nCalibration: model mean / observed mean in months assigned to the high state = %.2f (1.00 is exact; first run with the tight prior: 1.23) | residual SD = %.2f (about 1 is right; first run: 0.44)\n", ratio_hi, res_sd))
add("18. High-state level is right: model mean within 10% of the observed mean in months assigned to the high state", is.finite(ratio_hi) & abs(ratio_hi - 1) <= 0.10, sprintf("ratio %.2f", ratio_hi))
cat(sprintf("  Spread, three ways: randomized quantile residual SD %.2f (the check; about 1 is right at any number of zeros) | Pearson SD in months with cases %.2f | Pearson SD over all cells %.2f (old definition; shrinks with the zeros)\n", res_sd, res_sd_act, res_sd_pearson))
add("19. Spread of the model is right: randomized quantile residual SD between 0.8 and 1.25", is.finite(res_sd) & res_sd >= 0.8 & res_sd <= 1.25, sprintf("%.2f | Pearson, months with cases %.2f | Pearson, all cells (old) %.2f", res_sd, res_sd_act, res_sd_pearson))
V <- rbindlist(res); fwrite(V, file.path(out, sprintf("Table_S_HMM_small_check_%s.csv", TAGSTR)))
cat(sprintf("\n=== VERDICT (4 chains, prior on high-state gap: Gamma(%s, %s)) ===\n", TAU23[1], TAU23[2]))
if (all(V$ok, na.rm = TRUE)) cat("The HMM passes every check in this small run. Next: FAST, then FULL with the Stage 2 script.\n") else cat("NOT YET GOOD. Failed:", paste(V$question[which(!V$ok)], collapse = " || "), "\n")
cat("A small run can pass or fail only for lack of iterations (R-hat). Judge convergence on the FULL run.\n")

# ---- where does any over-prediction come from? (one draw, by state)
k <- sample(nd, 1); A <- Smat(X[k, , drop = FALSE])[1, , ]
lam_s <- matrix(X[k, sprintf("mu[%d,%d]", rep(1:Z, 3), rep(1:3, each = Z))] * exp(X[k, sprintf("delta[%d,%d]", rep(1:Z, 3), rep(1:3, each = Z))]), Z, 3)
lam <- matrix(lam_s[cbind(rep(1:Z, T), as.vector(A))], Z, T) * OFF * EFFM(k); st <- factor(as.vector(A), 1:3)
cat("\nObserved vs model mean count by state (one posterior draw):\n")
print(data.frame(state = 1:3, months = as.vector(table(st)), observed_mean = tapply(as.vector(Y), st, mean), model_mean = tapply(as.vector(lam), st, mean)))
cat("sum of model means:", round(sum(lam)), "| observed total:", sum(Y), "\n")

# ---- WHAT IS LEFT OVER? (can prediction still improve?) -------------------------------------------------------------
# Idea: if the model has captured the pattern, the leftover (observed minus expected, scaled) should show NO association with
# season, last month, neighbouring zones or population. A clear association shows where a correlation could be added later.
# This is a diagnostic only. It describes association in the residuals; it says nothing about cause.
cat("\n============ LEFT-OVER ASSOCIATIONS (residual check) ============\n")
kk <- sample(nd, min(50, nd)); LAM <- matrix(0, Z, T); RZ <- rep(0, Z)
for (k in kk) {                                                  # average the model's expected count over 50 posterior draws
  A <- Smat(X[k, , drop = FALSE])[1, , ]
  ls <- matrix(X[k, sprintf("mu[%d,%d]", rep(1:Z, 3), rep(1:3, each = Z))] * exp(X[k, sprintf("delta[%d,%d]", rep(1:Z, 3), rep(1:3, each = Z))]), Z, 3)
  LAM <- LAM + matrix(ls[cbind(rep(1:Z, T), as.vector(A))], Z, T) * OFF * EFFM(k) / length(kk); RZ <- RZ + X[k, sprintf("r[%d]", 1:Z)] / length(kk)
}
RES <- (Y - LAM) / sqrt(LAM + LAM^2 / RZ)                        # Pearson residual for a negative binomial; about mean 0, sd 1 if the model is right
left <- list()
rc <- function(name, a, b, what) { ok <- is.finite(a) & is.finite(b); r <- suppressWarnings(cor(a[ok], b[ok], method = "spearman")); left[[length(left) + 1]] <<- data.table(source = name, spearman = round(r, 3), reading = what)
  cat(sprintf("  %-34s Spearman %6.3f   %s\n", name, r, what)) }
# (a) same zone, previous month: is there association over time that the 3 states do not capture?
act <- as.vector(Y[, -1] > 0 | Y[, -T] > 0)                      # only zone-months where this or the previous month had cases: zero cells make residuals copy the model's own expected count
rc("previous month (same zone, active months)", as.vector(RES[, -1])[act], as.vector(RES[, -T])[act], "if clearly positive: states miss month-to-month persistence")
# (b) season: do residuals follow the calendar phase (month index modulo 12)?
ph <- rep(((1:T) - 1) %% 12 + 1, each = Z); mph <- tapply(as.vector(RES), ph, mean)
cat(sprintf("  %-34s range of mean residual across 12 phases: %.2f to %.2f   %s\n", "season (phase of the year)", min(mph), max(mph), "if beyond about +/-0.15: add season to the transitions"))
left[[length(left) + 1]] <- data.table(source = "season (phase of the year)", spearman = NA, reading = sprintf("phase means %.2f to %.2f", min(mph), max(mph)))
# (c) neighbouring zones: do zones next to each other have similar residuals in the same month?
Qg <- tryCatch(S0$graph_Q, error = function(e) NULL)
if (!is.null(Qg) && nrow(Qg) == Z) {
  Nb <- (abs(as.matrix(Qg)) > 0) * 1; diag(Nb) <- 0
  nbm <- (Nb %*% RES) / pmax(rowSums(Nb), 1)                    # average residual of the neighbours, same month
  rc("neighbouring zones (same month)", as.vector(RES), as.vector(nbm), "if clearly positive: nearby zones share leftover pattern (association only)")
} else cat("  neighbouring zones: graph_Q not found in the Stage 1 file; skipped\n")
# (d) population and zone size
rc("zone population (zone mean residual)", rowMeans(RES), as.numeric(popz), "if clearly non-zero: the offset does not fully remove size")
# (e) is the leftover too wide? (overdispersion check)
cat(sprintf("  residual sd %.2f (about 1 is right; clearly above 1 = still under-dispersed model)\n", sd(as.vector(RES), na.rm = TRUE)))
L <- rbindlist(left); fwrite(L, file.path(out, sprintf("Table_S_HMM_leftover_associations_%s.csv", TAGSTR)))
cat("Rule of thumb: |Spearman| below about 0.10 = nothing important left. Larger values tell us which term to add in the second paper.\n")


# ---- in-sample prediction error (for information; the hold-out above is the fair test) ----------------------------------------
ins <- data.table(method = c("HMM (in-sample)", "national rate x expected", "zone average"),
  MAE  = c(mean(abs(Y - LAM)), mean(abs(Y - sum(Y) / sum(OFF) * OFF)), mean(abs(Y - rowMeans(Y)))),
  RMSE = c(sqrt(mean((Y - LAM)^2)), sqrt(mean((Y - sum(Y) / sum(OFF) * OFF)^2)), sqrt(mean((Y - rowMeans(Y))^2))))
ins[, RMSE_reduction_pct := round(100 * (1 - RMSE[1] / RMSE), 1)]; ins[1, RMSE_reduction_pct := NA]
cat("\nIn-sample prediction error (the HMM has seen these data, so this flatters it):\n"); print(ins[, `:=`(MAE = round(MAE, 3), RMSE = round(RMSE, 2))])
fwrite(ins, file.path(out, sprintf("Table_S_HMM_insample_error_%s.csv", TAGSTR)))
if (!exists("V") || !exists("L")) stop("coverage table needs V (checks) and L (left-over) - an earlier block did not finish")
# ---- CONCERN COVERAGE: every concern and gap found so far, the check that tests it, and its result in THIS run ----------------------------------
chk_status <- function(ids) { r <- sapply(ids, function(i) { w <- which(startsWith(V$question, paste0(i, ". "))); if (!length(w)) NA else V$ok[w[1]] })
  if (all(is.na(r))) "not run" else if (all(r, na.rm = TRUE)) "PASS" else paste("FAIL:", paste(ids[which(!r)], collapse = ",")) }
left_val <- function(nm) { w <- grep(nm, L$source, fixed = TRUE); if (length(w)) sprintf("%.2f", L$spearman[w[1]]) else "n/a" }
COV <- data.table(
  concern = c("Population and offsets (Juan 2): do states follow population rather than cholera?", "Zone multipliers inflated every state mean (mean 1.30)", "High state over-predicted (206 vs 168)",
    "Model variance too wide (residual SD 0.44)", "Uncertainty carried from the first stage (Juan 3)", "Expected stay in each state (Juan 5)", "Chains agree (4 chains)",
    "Prediction: better than simple methods, honest intervals", "Four indicators can be extracted, are distinct and differ between zones", "Onset weak (rare direct jumps)",
    "Peak weak (few high months)", "Persistence prior-driven", "Left-over: previous month (active months)", "Left-over: neighbouring zones", "Left-over: population", "Left-over: season",
    "Recovery test, number of states, LPA, Louvain, classes, joint first-stage draws"),
  how_tested = c("checks 6, 7", "centring + checks 5b", "checks 5, 5b, 18", "check 19", "check 8", "check 4", "check 1", "checks 9, 10, 11", "checks 12, 13, 15", "check 16 (entry into any active state)",
    "check 17 (peak reached); option cholera.hiermean", "check 14; option cholera.kappa (pooling)", "left-over block", "left-over block", "left-over block", "left-over block", "Stage 2 (TEST, FAST, FULL)"),
  result = c(chk_status(c("6", "7")), if (isTRUE(all(abs(colMeans(exp(cbind(d_occ, d_int))) - 1) < 1e-6))) "PASS (centred)" else "applied", chk_status(c("5", "5b", "18")), chk_status("19"), chk_status("8"), chk_status("4"), chk_status("1"), chk_status(c("9", "10", "11")), chk_status(c("12", "13", "15")), chk_status("16"),
    chk_status("17"), chk_status("14"), paste("Spearman", left_val("previous month")), paste("Spearman", left_val("neighbouring")), paste("Spearman", left_val("population")), if (any(grepl("season", L$source))) L$reading[grep("season", L$source)[1]] else "n/a", "not in this script"))
cat("\n============ CONCERN COVERAGE (this run) ============\n"); print(COV[, .(concern = substr(concern, 1, 70), how_tested, result)], right = FALSE)
fwrite(COV, file.path(out, sprintf("Table_S_HMM_concern_coverage_%s.csv", TAGSTR)))
cat("Left-over rule of thumb: |Spearman| below about 0.10 = nothing important left. A FAIL or a large left-over value says which option to try next (prior, cholera.hiermean, cholera.kappa, cholera.pop, cholera.season, cholera.nb).\n")

# ================= UPDATED INDICATORS + REMAINING GAPS =================
cat("\n============ UPDATED INDICATOR SET: burden | entry | persistence of active months | peak reached ============\n")
dd <- dim(Aall)
act_prev <- Aall[, , -dd[3], drop = FALSE] >= 2; act_now <- Aall[, , -1, drop = FALSE] >= 2
pers_act <- apply(act_prev & act_now, c(1, 2), sum) / pmax(apply(act_prev, c(1, 2), sum), 1)
m_act <- apply(Aall >= 2, c(1, 2), sum)
NEW4 <- list(burden = bur_d, entry = en_d, persistence_active = pers_act, peak_reached = pr_d)
nt <- rbindlist(lapply(names(NEW4), function(n) { M <- NEW4[[n]]; pm <- colMeans(M); ps <- apply(M, 2, sd)
  data.table(indicator = n, median = signif(median(pm), 3), min = signif(min(pm), 3), max = signif(max(pm), 3), ratio = round(sd(pm) / mean(ps), 2),
             spearman_population = round(cor(pm, as.numeric(popz), method = "spearman"), 2)) }))
print(nt)
IMn <- sapply(NEW4, colMeans)
pen <- cor(IMn, method = "pearson"); cat("Pearson between the updated indicators:\n"); print(round(pen, 3))
maxn <- max(abs(pen[upper.tri(pen)])); condn <- sqrt(max(eigen(cor(scale(IMn)))$values) / min(eigen(cor(scale(IMn)))$values))
cat(sprintf("Highest |Pearson| %.3f | condition number %.1f | old set (burden, onset, persistence, peak): %.3f / %.1f\n", maxn, condn, max_cor, cond_num))
h1 <- 1:floor(nd / 2); h2 <- (floor(nd / 2) + 1):nd
stab_half <- sapply(NEW4, function(M) cor(colMeans(M[h1, , drop = FALSE]), colMeans(M[h2, , drop = FALSE]), method = "spearman"))
cat("Split-half stability of the zone ranking (Spearman, chains 1-2 vs 3-4):\n"); print(round(stab_half, 3))
bm <- colMeans(bur_d)
anch_sp <- sapply(c(0.90, 0.95, 0.975), function(p) cor(bm, rowMeans(Y >= quantile(Y[Y > 0], p)), method = "spearman"))
names(anch_sp) <- c("p90", "p95", "p975"); cat("Anchor sensitivity (Spearman, model burden vs observed share of months above the threshold):\n"); print(round(anch_sp, 3))
thin_act <- colMeans(m_act) < 3
cat(sprintf("Zones with fewer than 3 expected active months: %d of %d | with fewer than 3 high months: %d of %d\n", sum(thin_act), Z, sum(few_high), Z))
if (exists("sat_zones")) { inn <- zone_order %in% sat_zones
  cat(sprintf("SaTScan comparison: zones inside clusters (n = %d) vs outside, median burden %.4f vs %.4f, median peak reached %.0f vs %.0f\n", sum(inn), median(IMn[inn, 1]), median(IMn[!inn, 1]), median(IMn[inn, 4]), median(IMn[!inn, 4]))) } else cat("SaTScan comparison skipped (create sat_zones <- c(...zone names...) to run it).\n")
zt <- data.table(zone = zone_order, burden = IMn[, 1], entry = IMn[, 2], persistence_active = IMn[, 3], peak_reached = IMn[, 4],
                 expected_active_months = colMeans(m_act), thin_active = thin_act, thin_high = few_high)
zt[, `:=`(z_burden = as.numeric(scale(log(burden + 1e-4))), z_entry = as.numeric(scale(log(entry + 1e-4))), z_persistence = as.numeric(scale(persistence_active)), z_peak = as.numeric(scale(log(peak_reached + 1))))]
for (n in names(NEW4)) { zt[[paste0(n, "_lo")]] <- apply(NEW4[[n]], 2, quantile, .025); zt[[paste0(n, "_hi")]] <- apply(NEW4[[n]], 2, quantile, .975) }
saveRDS(list(ind = NEW4, zone_order = zone_order), file.path(out, sprintf("Draws_updated_indicators_%s.rds", TAGSTR)))   # every posterior draw of the four indicators, for the stability analysis in Stage 2
fwrite(zt, file.path(out, sprintf("Table_S_HMM_updated_indicators_for_LPA_%s.csv", TAGSTR)))
fwrite(nt, file.path(out, sprintf("Table_S_HMM_updated_indicator_summary_%s.csv", TAGSTR)))
add("20. Updated set: at least 3 of 4 indicators differ between zones by more than their uncertainty (ratio above 1)", sum(nt$ratio > 1) >= 3, paste(nt$indicator, nt$ratio, collapse = " | "))
add("21. Updated set is distinct (highest |Pearson| below 0.85 and condition number below 30)", maxn < 0.85 & condn < 30, sprintf("max |r| %.2f | condition number %.1f", maxn, condn))
add("22. No updated indicator just follows zone population (|Spearman| below 0.5 for all four)", all(abs(nt$spearman_population) < 0.5), paste(nt$indicator, nt$spearman_population, collapse = " | "))
add("23. Zone burden ranking does not depend on the threshold (Spearman at least 0.8 for the 90th, 95th and 97.5th percentile)", all(anch_sp >= 0.8), paste(names(anch_sp), round(anch_sp, 2), collapse = " | "))
add("24. Zone rankings are stable between chain halves (Spearman at least 0.9 for all four)", all(stab_half >= 0.9), paste(names(stab_half), round(stab_half, 2), collapse = " | "))
add("25. Few zones lean on the prior for the updated set (under 25% of zones with fewer than 3 expected active months)", mean(thin_act) < 0.25, sprintf("%.0f%% of zones", 100 * mean(thin_act)))
V <- rbindlist(res); fwrite(V, file.path(out, sprintf("Table_S_HMM_small_check_%s.csv", TAGSTR)))
COV2 <- rbind(COV, data.table(concern = c("Persistence prior-driven: use persistence of active months", "Indicators correlated with each other (updated set)", "Indicators follow population", "Threshold (anchor) dependence",
    "Stability of zone rankings", "Zones with too few events"),
  how_tested = c("check 20 (ratio)", "check 21", "check 22", "check 23", "check 24", "check 25"),
  result = c(chk_status("20"), chk_status("21"), chk_status("22"), chk_status("23"), chk_status("24"), chk_status("25"))))
cat("\n============ CONCERN COVERAGE (with updated indicators) ============\n"); print(COV2[, .(concern = substr(concern, 1, 70), how_tested, result)], right = FALSE)
fwrite(COV2, file.path(out, sprintf("Table_S_HMM_concern_coverage_%s.csv", TAGSTR)))

# ================= END-TO-END SUMMARY =================
fix_map <- list(
  "1"  = "Chains disagree: run longer (FAST/FULL). In a small run this is expected; judge on FULL.",
  "2"  = "States not ordered: check the ordering constraint in the JAGS model (mu1 < mu2 < mu3).",
  "4"  = "High-state stay is out of range: try cholera.kappa = 30 (pool transition rows) or a different CONC.",
  "5"  = "Model does not reproduce the data: try cholera.hiermean = TRUE and cholera.rmin = 0.5.",
  "5b" = "Predictive p-values extreme: same as check 5; look at which statistic is off in the table above.",
  "6"  = "Offset changes the ranking a lot: report both rankings; population is influencing the states.",
  "7"  = "States follow population: add cholera.pop = TRUE and report the population coefficient.",
  "8"  = "Zones switch groups between draws: report intervals; do not label single zones as high or low.",
  "9"  = "HMM does not beat simple baselines: make NO forecasting claim; show the hold-out in the supplement.",
  "10" = "Prediction intervals too narrow or wide: check the NB size r (cholera.rmin) and the prior on the high state.",
  "11" = "Hidden-month total outside the interval: the model misses a national trend; try cholera.season = TRUE.",
  "12" = "Indicators not extractable: a state is never visited; check the thresholds and the data for that zone.",
  "13" = "Indicators do not separate zones: use entry and active persistence; report intervals.",
  "14" = "Many zones with few high months: use active-month indicators; flag these zones in the paper.",
  "15" = "Indicators too correlated: drop one (burden or entry) or use the updated set.",
  "16" = "Entry does not separate zones: keep onset out of the LPA; use burden, active persistence and peak reached.",
  "17" = "Peak reached does not follow the observed maximum: try cholera.hiermean = TRUE or a centred prior (cholera.anchhi = 170).",
  "18" = "High-state level off by more than 10%: try cholera.hiermean = TRUE or a centred prior (cholera.anchhi = 170 or 250).",
  "19" = "Residual SD outside 0.8 to 1.25: try cholera.rmin = 0.5 (too wide) or cholera.nb / cholera.season (too narrow).",
  "20" = "Updated indicators do not separate zones: report intervals and cluster with fewer indicators.",
  "21" = "Updated indicators too correlated: drop one before the LPA.",
  "22" = "An indicator follows population: report this and test the LPA with and without it.",
  "23" = "Burden depends on the threshold: add a sensitivity run with the 90th and 97.5th percentile anchors.",
  "24" = "Zone rankings unstable: raise the number of iterations; do not cluster on this indicator yet.",
  "25" = "Many zones with few active months: flag them and report class membership probabilities, not single labels.")
num <- sub("\\..*$", "", V$question)
V[, id := num]; V[, status := ifelse(is.na(ok), "SKIP", ifelse(ok, "PASS", "FAIL"))]
fails <- V[status == "FAIL"]
cat("\n==================== SUMMARY ====================\n")
cat(sprintf("Model version: %s | prior Gamma(%s, %s) | checks: %d pass, %d fail, %d skipped\n", if (EXT == "") "Methods default" else EXT, TAU23[1], TAU23[2], sum(V$status == "PASS"), nrow(fails), sum(V$status == "SKIP")))
cat("\nKey numbers\n")
cat(sprintf("  high-state model/observed mean %.2f (target 1.00) | residual SD %.2f (target about 1)\n", ratio_hi, res_sd))
cat(sprintf("  expected stay: baseline %.1f, moderate %.1f, high %.1f months\n", mean(1 / (1 - gam_self[, 1])), mean(1 / (1 - gam_self[, 2])), mean(1 / (1 - gam_self[, 3]))))
if (exists("tab")) cat(sprintf("  hold-out CRPS: HMM %.2f | always zero %.2f | best simple %.2f | 95%% interval coverage %.1f%%\n", tab$CRPS[1], tab$CRPS[tab$method == "always zero"], best_simple, 100 * cover))
cat(sprintf("  updated indicators, between/within ratio: %s\n", paste(nt$indicator, nt$ratio, collapse = " | ")))
cat(sprintf("  zones with few high months %d of %d | few active months %d of %d\n", sum(few_high), Z, sum(thin_act), Z))
if (nrow(fails)) {
  cat("\nGAPS and what to do next\n")
  for (i in seq_len(nrow(fails))) cat(sprintf("  [%s] %s\n       detail: %s\n       next:   %s\n", fails$id[i], substr(fails$question[i], 1, 110), fails$detail[i], if (!is.null(fix_map[[fails$id[i]]])) fix_map[[fails$id[i]]] else "see the check text"))
} else cat("\nNo gaps found in this run. Next: add the chosen options and the updated indicators to Stage 2, then TEST, FAST, FULL.\n")
cat("\nOptions to try one at a time at the top of the script (then rerun):\n  options(cholera.hiermean = TRUE)   # learn the high-state level\n  options(cholera.anchhi = 170)      # centre the prior near the observed high-month mean\n  options(cholera.rmin = 0.5)        # tighten the spread\n  options(cholera.kappa = 30)        # pool the transition rows\n  options(cholera.season = TRUE); options(cholera.nb = TRUE); options(cholera.pop = TRUE)\n")
rep_txt <- c(sprintf("HMM check summary | %s | %s", TAGSTR, format(Sys.time(), "%Y-%m-%d %H:%M")),
  sprintf("Checks: %d pass, %d fail, %d skipped", sum(V$status == "PASS"), nrow(fails), sum(V$status == "SKIP")),
  sprintf("High-state model/observed %.2f | residual SD %.2f", ratio_hi, res_sd),
  paste0("FAIL ", fails$id, ": ", fails$detail),
  capture.output(print(COV2[, .(concern = substr(concern, 1, 70), how_tested, result)], right = FALSE)),
  capture.output(print(nt)))
writeLines(rep_txt, file.path(out, sprintf("HMM_check_summary_%s.txt", TAGSTR)))
fwrite(V[, .(id, question, status, detail)], file.path(out, sprintf("Table_S_HMM_check_status_%s.csv", TAGSTR)))
cat("\nSummary saved to:", file.path(out, sprintf("HMM_check_summary_%s.txt", TAGSTR)), "\n")

# ================= JUAN'S COMMENTS: where each one is answered, and what is still outside this script =================
JC <- data.table(
  comment = c("C1 Language: association and clusters of estimated risk, not transmission", "C2 Offset and population (E, smoothing, E = 0, population in the HMM)",
              "C3 Uncertainty carried between stages (stability of profiles, classes, communities)", "C4 Data clarity (true zeros, case definition, missing data, boundaries)",
              "C5 Interpretation (expected stay, waves, beta, model scores)"),
  answered_by = c("Wording in Methods, Results and Discussion", "Checks 6, 7, 22 (offset changes little, states do not copy population, indicators do not follow population); optional population term (variant J)",
                  "First-stage zone effects enter with their own SD; checks 8, 24, 25", "Your data paragraph (not a code check)", "Check 4 and the expected-stay table (1 / (1 - p_ii))"),
  result_here = c("Text only", chk_status(c("6", "7", "22")), chk_status(c("8", "24", "25")), "Text only", chk_status("4")),
  evidence = c("-",
    sprintf("min E %.3f | Spearman burden with vs without offset %.2f | burden vs population %.2f", min(dat$E), cor(b1, b0, method = "spearman"), cor(b1, popz, method = "spearman")),
    sprintf("%.0f%% of zones keep their burden group | split-half Spearman min %.2f | %.0f%% of zones with under 3 active months", 100 * mean(stab >= 0.7), min(stab_half), 100 * mean(thin_act)),
    sprintf("zeros %.1f%% | %s cases | %d zones x %d months", 100 * mean(Y == 0), format(sum(Y), big.mark = ","), Z, T),
    sprintf("expected stay: baseline %.1f, moderate %.1f, high %.1f months", mean(1 / (1 - gam_self[, 1])), mean(1 / (1 - gam_self[, 2])), mean(1 / (1 - gam_self[, 3])))),
  still_outside_this_script = c("Check every sentence for transmission, hub, hotspot, corridor, gateway, amplifier wording",
    "Stage 1: pooled national-rate E (no smoothing, no floor); report the sensitivity with unsmoothed or population-based E; explain E = 0",
    "Stage 2: LPA and Louvain stability across posterior draws (replaces the plus or minus 20% edge perturbation)",
    "True zeros vs non-reporting; confirmed vs suspected vs epi-linked; missing-data handling; how 2015 to 2024 boundary changes were harmonised to 92 zones",
    "Waves come from the national sensitivity model (say so in the abstract); beta as scaling factor; WAIC, DIC and log-CPO disagreement; log scale for the RR figure"))
cat("\n============ JUAN'S COMMENTS (this run) ============\n"); print(JC[, .(comment = substr(comment, 1, 60), result_here, evidence)], right = FALSE)
cat("\nStill to do outside this script:\n"); for (i in seq_len(nrow(JC))) cat(sprintf("  %s\n      -> %s\n", substr(JC$comment[i], 1, 60), JC$still_outside_this_script[i]))
fwrite(JC, file.path(out, sprintf("Table_S_HMM_juan_comment_coverage_%s.csv", TAGSTR)))

# ---- CONVERGENCE of the zone-specific parameters, and the sentences for Methods and Results (every number is taken from this run) -------------------------
cat("\n============ CONVERGENCE OF THE ZONE-SPECIFIC PARAMETERS ============\n")
cvr <- function(nm, cols) { m <- mc[, cols]; x <- gelman.diag(m, multivariate = FALSE, autoburnin = FALSE)$psrf[, 1]; e <- effectiveSize(m)
  data.table(parameter = nm, n = length(cols), n_rhat_below_1.1 = sum(x < 1.1, na.rm = TRUE), pct_rhat_below_1.1 = round(100 * mean(x < 1.1, na.rm = TRUE), 1), n_rhat_below_1.05 = sum(x < 1.05, na.rm = TRUE), max_rhat = round(max(x, na.rm = TRUE), 3), min_ess = round(min(e))) }
mu3c <- sprintf("mu[%d,3]", 1:Z)
CV <- rbindlist(list(cvr("zone high-intensity mean, mu[z,3]", mu3c), cvr("zone moderate mean, mu[z,2]", sprintf("mu[%d,2]", 1:Z)), cvr("zone baseline mean, mu[z,1]", sprintf("mu[%d,1]", 1:Z)),
  cvr("zone dispersion, r[z]", sprintf("r[%d]", 1:Z)), cvr("self-transition, baseline state", sprintf("gam[%d,1,1]", 1:Z)), cvr("self-transition, moderate state", sprintf("gam[%d,2,2]", 1:Z)), cvr("self-transition, high state", sprintf("gam[%d,3,3]", 1:Z)),
  if (HIERMEAN) cvr("national high-state level, eta23", "eta23")))
print(CV, row.names = FALSE); fwrite(CV, file.path(out, sprintf("Table_S_HMM_convergence_%s.csv", TAGSTR)))
mu3 <- mc[, mu3c]
rh3 <- gelman.diag(mu3, multivariate = FALSE, autoburnin = FALSE)$psrf[, 1]
rh3l <- gelman.diag(do.call(coda::mcmc.list, lapply(mu3, function(ch) { mp <- attr(ch, "mcpar"); coda::mcmc(log(as.matrix(ch)), start = mp[1], thin = mp[3]) })), multivariate = FALSE, autoburnin = FALSE)$psrf[, 1]   # R-hat on the log scale as well
ess3 <- effectiveSize(mu3); if (!exists("Aall")) Aall <- Smat(X); eh3 <- apply(Aall == 3, 2, mean) * T
zc <- data.table(zone = zone_order, rhat_mu_high = round(rh3, 3), rhat_log_mu_high = round(rh3l, 3), ess_mu_high = round(ess3), mean_mu_high = round(colMeans(X[, mu3c]), 1), expected_high_months = round(eh3, 1), converged = rh3 < 1.1)
fwrite(zc, file.path(out, sprintf("Table_S_HMM_convergence_by_zone_%s.csv", TAGSTR)))
k3 <- sum(rh3 < 1.1, na.rm = TRUE)
if (k3 < Z) { cat(sprintf("\n%d of %d zones have R-hat of 1.1 or more for the high-intensity mean. They are listed below with how often they visit the high state:\n", Z - k3, Z)); print(zc[converged == FALSE | is.na(converged)], row.names = FALSE)
  cat(sprintf("  median expected high months: zones not converged %.1f | zones converged %.1f (a zone that rarely visits the high state gives the data little to learn its high-state mean from)\n", median(zc[converged == FALSE, expected_high_months]), median(zc[converged == TRUE, expected_high_months]))) } else cat("\nAll zone-specific high-intensity means have R-hat below 1.1 (100%).\n")
fmt <- function(x) format(round(x), big.mark = ",", trim = TRUE, scientific = FALSE)
sen <- c("METHODS", sprintf("Parameter estimation was performed using Bayesian MCMC simulation across %d independent parallel chains of %s iterations each, with a %s-iteration burn-in period and a thinning interval of %d, retaining %s draws per chain (%s in total).", NCH, fmt(BURN + NITER), fmt(BURN), THIN, fmt(NITER / THIN), fmt(NCH * NITER / THIN)),
  "", "RESULTS (convergence)",
  if (k3 == Z) sprintf("All %d zone-specific high-intensity mean parameters achieved convergence (100%%; R-hat below 1.1, maximum %.3f; smallest effective sample size %s).", Z, max(rh3, na.rm = TRUE), fmt(min(ess3)))
  else sprintf("%.1f%% (%d of %d) of the zone-specific high-intensity mean parameters achieved convergence (R-hat below 1.1; maximum %.3f; smallest effective sample size %s). The %d zones that did not converge had a median of %.1f expected high-state months, against %.1f for the converged zones.", 100 * k3 / Z, k3, Z, max(rh3, na.rm = TRUE), fmt(min(ess3)), Z - k3, median(zc[converged == FALSE, expected_high_months]), median(zc[converged == TRUE, expected_high_months])),
  sprintf("For the zone-specific moderate and baseline means, %.1f%% and %.1f%% of zones had R-hat below 1.1; for the three self-transition probabilities, %.1f%%, %.1f%% and %.1f%%.", CV$pct_rhat_below_1.1[2], CV$pct_rhat_below_1.1[3], CV$pct_rhat_below_1.1[5], CV$pct_rhat_below_1.1[6], CV$pct_rhat_below_1.1[7]),
  "", "RESULTS (fit and expected stay)",
  sprintf("The expected stay in a state is 1 / (1 - p) months for self-transition probability p: baseline %.1f, moderate %.1f and high %.1f months (posterior means over zones).", mean(1 / (1 - gam_self[, 1])), mean(1 / (1 - gam_self[, 2])), mean(1 / (1 - gam_self[, 3]))),
  sprintf("In months assigned to the high state, the model mean was %.2f times the observed mean (1.00 is exact); the randomized quantile residual SD was %.2f (about 1 when the spread is right).", ratio_hi, res_sd))
writeLines(sen, file.path(out, sprintf("Methods_Results_sentences_%s.txt", TAGSTR))); cat("\n"); cat(sen, sep = "\n"); cat("\n(Also saved as Methods_Results_sentences_", TAGSTR, ".txt)\n", sep = "")

try(stopCluster(cl_hmm), silent = TRUE)
if (nzchar(Sys.getenv("TAGFILE"))) writeLines(TAGSTR, Sys.getenv("TAGFILE"))   # lets Run_all_variants.R find this run's files
options(cholera.last_tag = TAGSTR)
