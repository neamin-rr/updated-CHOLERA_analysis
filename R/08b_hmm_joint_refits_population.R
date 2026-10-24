# ==============================================================================
# Joint_20draws_long_pop.R  -  sensitivity analysis (4) for the NEW main HMM (with log population): 20 joint Stage 1 draws, 2 chains per draw, R-hat per refit
# Same HMM as the final population run (prior centre 170, wide prior, shared national high-state gap, r >= 0.5, b_pop * standardised log population).
# Can start before the final run has finished: the refits run now; the comparison with the main run is done when the script is run again after the final run
# (finished chains are reloaded in seconds).
# Needs: rjags, coda, data.table, parallel (no INLA). Uses joint_fields_20.rds (made by make_joint_fields.R) in update_HMM.
# Per draw: 2 chains x (20,000 burn-in + 20,000 kept, thin 20) = 2,000 kept draws per refit.
# Run:   cd ~/Documents/spatial_anlysis/update_HMM
#        Rscript -e 'options(cholera.check = TRUE); source("Joint_20draws_long_pop.R")'          # quick check (seconds)
#        Rscript -e 'options(cholera.smoke = TRUE); source("Joint_20draws_long_pop.R")'          # end-to-end test (a few minutes), results meaningless
#        mkdir -p joint_20draws_long_pop ; setsid nohup nice -n 19 Rscript Joint_20draws_long_pop.R > joint_20draws_long_pop/log.txt 2>&1 &
#        tail -f joint_20draws_long_pop/log.txt
# Workers: options(cholera.cores = N) or env CORES=N; default = number of cores - 4, at most 40 (one chain per worker). Each chain is saved when done (resumable).
# Output (folder update_HMM/joint_20draws_long_pop): chain files, Table_S_joint20_long_per_draw.csv, Table_S_joint20_long_comparison.csv,
#        Table_S_joint20_long_indicators_by_zone.csv, Draws_joint20_long_indicators.rds, Joint20_long_sentences.txt
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(parallel); library(coda) })
proj  <- getOption("cholera.proj", ".")
OUTJ  <- file.path(proj, "joint_20draws_long_pop"); dir.create(OUTJ, showWarnings = FALSE, recursive = TRUE)
FIELDS <- file.path(proj, "joint_fields_20.rds")
MAINFIT <- file.path(proj, "final_pop", "improve_hmm", "cache_tau2_2_hi170_pop_rmin0.5_hm_main_n100000_b200000.rds")   # the final population run
N_JOINT <- 20; NCHD <- 2; BURN <- 20000; NIT <- 20000; THIN <- 20
if (isTRUE(getOption("cholera.smoke", FALSE))) { N_JOINT <- 2; BURN <- 50; NIT <- 50; THIN <- 1; OUTJ <- file.path(proj, "joint_20draws_smoke"); dir.create(OUTJ, showWarnings = FALSE) }   # 2-minute test of the whole script; results meaningless
NW <- as.integer(getOption("cholera.cores", Sys.getenv("CORES", min(40, max(1, detectCores() - 4)))))
set.seed(2026)
say <- function(...) cat(format(Sys.time(), "%d %b %H:%M:%S"), ..., "\n")

S0 <- readRDS(file.path(proj, "INLA_Rebuilt_From_Zero_OCCPOP.rds")); dat <- as.data.table(S0$data); setorder(dat, ID_area, ID_time)
zone_order <- S0$zone_order; Z <- length(zone_order); T <- max(dat$ID_time); stopifnot(nrow(dat) == Z * T)
Y <- matrix(dat$O, Z, T, byrow = TRUE); Em <- matrix(dat$E, Z, T, byrow = TRUE); OFF <- Em / mean(Em)
zpop <- as.numeric(scale(log(as.numeric(tapply(dat$monthly_population, dat$ID_area, mean)))))   # standardised log zone population (as in HMM_check_ALL.R)
if (!file.exists(FIELDS)) stop("joint_fields_20.rds not found in ", proj, ": run make_joint_fields.R first")
JF <- readRDS(FIELDS); stopifnot(identical(as.character(JF$zone_order), as.character(zone_order)))
DO <- JF$DO[, 1:N_JOINT]; DI <- JF$DI[, 1:N_JOINT]
if (isTRUE(getOption("cholera.check", FALSE))) {
  say("CHECK OK: fields file found (", N_JOINT, "draws )")
  say("Main final fit found:", file.exists(MAINFIT), "(not needed yet: the comparison runs when the script is run again after the final run)")
  say("rjags:", requireNamespace("rjags", quietly = TRUE), "| cores:", detectCores(), "| workers to use:", NW, "| chains in total:", N_JOINT * NCHD)
  say("Free memory (GB):", tryCatch(system("free -g | awk 'NR==2{print $7}'", intern = TRUE), error = function(e) "?"))
  say("Output folder:", OUTJ); stop("check finished (this is not an error): now start the full run.", call. = FALSE) }
thr <- as.numeric(quantile(Y[Y > 0], 0.95)); anch2 <- median(Y[Y > 0 & Y < thr]); anch1 <- 0.01; ANCHHI <- 170
st_emp <- matrix(1L, Z, T); st_emp[Y > 0 & Y < thr] <- 2L; st_emp[Y >= thr] <- 3L
trm <- matrix(0, 3, 3); for (z in 1:Z) for (t in 2:T) { a <- st_emp[z, t - 1]; b <- st_emp[z, t]; trm[a, b] <- trm[a, b] + 1 }
self_emp <- diag(trm) / rowSums(trm); CONC <- 10
a_self <- pmax(self_emp, 0.05) * CONC; a_self[3] <- a_self[2]; a_oth <- (1 - c(self_emp[1], self_emp[2], self_emp[2])) * CONC / 2
alpha <- matrix(0, 3, 3); for (i in 1:3) { alpha[i, ] <- a_oth[i]; alpha[i, i] <- a_self[i] }
say(sprintf("Priors: self-transition rates %s; transitions out of high: %d", paste(round(self_emp, 3), collapse = " / "), sum(trm[3, ])))

# ---- the same JAGS model as the final run
model_txt <- "
model {
  pi0[1:3] ~ ddirch(c(8, 1, 1))
  b_pop ~ dnorm(0, 1)
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
    for (i in 1:3) { gam[z, i, 1:3] ~ ddirch(alpha[i, 1:3]) }
    S[z, 1] ~ dcat(pi0[1:3])
    for (t in 2:T) { S[z, t] ~ dcat(gam[z, S[z, t - 1], 1:3]) }
    for (t in 1:T) {
      lam[z, t] <- max(lam_s[z, S[z, t]] * OFF[z, t] * exp(b_pop * zpop[z]), 1e-6)
      p[z, t] <- r[z] / (r[z] + lam[z, t])
      Y[z, t] ~ dnegbin(p[z, t], r[z])
    }
  }
}"

mf <- file.path(OUTJ, "hmm_model_joint.jags"); writeLines(model_txt, mf)
jdata <- list(Y = Y, OFF = OFF, Z = Z, T = T, alpha = alpha, anch1 = anch1, anch2 = anch2, anch3 = ANCHHI, tau_a = 2, tau_b = 2, prec_eta = 1, rmin = 0.5, zpop = zpop)

## one chain of one draw; returns per-iteration zone indicators (draws x zones) and the parameters needed for R-hat
one_chain <- function(task) {
  j <- task$j; ch <- task$ch
  f <- file.path(OUTJ, sprintf("draw%02d_chain%d.rds", j, ch)); if (file.exists(f)) { r <- tryCatch(readRDS(f), error = function(e) NULL); if (!is.null(r)) return(r) }
  suppressPackageStartupMessages({ library(rjags); library(coda) })
  jd <- jdata; jd$d_occ <- DO[, j]; jd$d_int <- DI[, j]; jd$prec_occ <- rep(1e4, Z); jd$prec_int <- rep(1e4, Z)
  m <- jags.model(mf, data = jd, inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 7000 + 10 * j + ch, S = st_emp), n.chains = 1, quiet = TRUE, n.adapt = 1000)
  update(m, BURN)
  x <- as.matrix(coda.samples(m, c("S", "gam", "mu", "delta", "b_pop"), n.iter = NIT, thin = THIN)[[1]]); nd <- nrow(x)
  A <- array(x[, sprintf("S[%d,%d]", rep(1:Z, T), rep(1:T, each = Z))], c(nd, Z, T))
  g <- function(i, k) x[, sprintf("gam[%d,%d,%d]", 1:Z, i, k)]
  lvl <- function(s) x[, sprintf("mu[%d,%d]", 1:Z, s)] * exp(x[, sprintf("delta[%d,%d]", 1:Z, s)]) * exp(outer(x[, "b_pop"], zpop))   # zone level incl. its population term
  mx <- apply(A, c(1, 2), max); act_prev <- A[, , -T, drop = FALSE] >= 2; act_now <- A[, , -1, drop = FALSE] >= 2
  r <- list(burden = apply(A == 3, c(1, 2), mean), entry = g(1, 2) + g(1, 3),
            persistence_active = apply(act_prev & act_now, c(1, 2), sum) / pmax(apply(act_prev, c(1, 2), sum), 1),
            peak_reached = lvl(3) * (mx == 3) + lvl(2) * (mx == 2) + lvl(1) * (mx == 1),
            gam11 = g(1, 1), gam22 = g(2, 2), gam33 = g(3, 3), mu3 = x[, sprintf("mu[%d,3]", 1:Z)], b_pop = x[, "b_pop"])
  saveRDS(r, paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f); r }

tasks <- unlist(lapply(1:N_JOINT, function(j) lapply(1:NCHD, function(ch) list(j = j, ch = ch))), recursive = FALSE)
say(sprintf("Running %d draws x %d chains = %d chains (%d burn-in + %d kept, thin %d) on %d workers; results in %s",
            N_JOINT, NCHD, length(tasks), BURN, NIT, THIN, min(NW, length(tasks)), OUTJ))
cl <- makeCluster(min(NW, length(tasks)), port = 22000L + as.integer(Sys.getpid()) %% 15000L)
clusterExport(cl, c("one_chain", "OUTJ", "jdata", "DO", "DI", "mf", "st_emp", "Z", "T", "BURN", "NIT", "THIN", "zpop"), envir = environment())
RES <- parLapplyLB(cl, tasks, one_chain); stopCluster(cl)
say("All chains finished.")

## ---- per draw: R-hat (2 chains) and pooled zone indicators
vars <- c("burden", "entry", "persistence_active", "peak_reached")
key <- sapply(tasks, function(t) t$j)
rh <- function(a, b) { v <- tryCatch(gelman.diag(mcmc.list(mcmc(a), mcmc(b)), autoburnin = FALSE, multivariate = FALSE)$psrf[, 1], error = function(e) rep(NA, ncol(a))); v }
bp <- sapply(1:N_JOINT, function(j) mean(c(RES[key == j][[1]]$b_pop, RES[key == j][[2]]$b_pop)))
say(sprintf("Population effect b_pop across refits: mean %.3f (range %.3f to %.3f)", mean(bp), min(bp), max(bp)))
per_draw <- rbindlist(lapply(1:N_JOINT, function(j) { r <- RES[key == j]
  rhs <- c(rh(r[[1]]$gam11, r[[2]]$gam11), rh(r[[1]]$gam22, r[[2]]$gam22), rh(r[[1]]$gam33, r[[2]]$gam33), rh(log(r[[1]]$mu3), log(r[[2]]$mu3)))
  data.table(draw = j, rhat_max = max(rhs, na.rm = TRUE), share_rhat_below_1.1 = mean(rhs < 1.1, na.rm = TRUE),
             gam11 = mean(c(r[[1]]$gam11, r[[2]]$gam11)), gam33 = mean(c(r[[1]]$gam33, r[[2]]$gam33))) }))
IND <- lapply(setNames(vars, vars), function(v) t(sapply(1:N_JOINT, function(j) { r <- RES[key == j]; colMeans(rbind(r[[1]][[v]], r[[2]][[v]])) })))   # N_JOINT x Z
saveRDS(list(ind = IND, zone_order = zone_order), file.path(OUTJ, "Draws_joint20_long_indicators.rds"))
fwrite(per_draw, file.path(OUTJ, "Table_S_joint20_long_per_draw.csv"))
say(sprintf("R-hat per refit: maximum over all refits %.3f; share of parameters below 1.1 = %.1f%% (lowest refit %.1f%%)",
            max(per_draw$rhat_max), 100 * mean(per_draw$share_rhat_below_1.1), 100 * min(per_draw$share_rhat_below_1.1)))

if (!file.exists(MAINFIT)) { say("Refits done. The final population run is not finished yet (", basename(MAINFIT), " not found):")
  say("run this script again when it is - the finished chains are reloaded in seconds and the comparison is made."); quit(save = "no", status = 0) }
## ---- main final run, same indicator definitions (from the saved fit; 2,000 draws are enough for zone means)
say("Computing the main-run indicators from the final fit (a few minutes) ...")
X <- as.matrix(readRDS(MAINFIT)); X <- X[sort(sample(nrow(X), min(2000, nrow(X)))), ]; nd <- nrow(X)
A <- array(X[, sprintf("S[%d,%d]", rep(1:Z, T), rep(1:T, each = Z))], c(nd, Z, T))
g <- function(i, k) X[, sprintf("gam[%d,%d,%d]", 1:Z, i, k)]
lvl <- function(s) X[, sprintf("mu[%d,%d]", 1:Z, s)] * exp(X[, sprintf("delta[%d,%d]", 1:Z, s)]) * exp(outer(X[, "b_pop"], zpop))
mx <- apply(A, c(1, 2), max); act_prev <- A[, , -T, drop = FALSE] >= 2; act_now <- A[, , -1, drop = FALSE] >= 2
MAIN <- list(burden = apply(A == 3, c(1, 2), mean), entry = g(1, 2) + g(1, 3),
             persistence_active = apply(act_prev & act_now, c(1, 2), sum) / pmax(apply(act_prev, c(1, 2), sum), 1),
             peak_reached = lvl(3) * (mx == 3) + lvl(2) * (mx == 2) + lvl(1) * (mx == 1))
saveRDS(list(ind = MAIN, zone_order = zone_order), file.path(proj, "final_pop", "improve_hmm", "Draws_updated_indicators_tau2_2_hi170_pop_rmin0.5_hm_FINAL.rds"))   # also the file Stage 3 needs
mm <- sapply(MAIN, colMeans); pooled <- sapply(IND, colMeans)
sp <- sapply(vars, function(v) cor(pooled[, v], mm[, v], method = "spearman"))
sp_each <- sapply(vars, function(v) apply(IND[[v]], 1, function(r) cor(r, mm[, v], method = "spearman")))   # N_JOINT x 4
top10 <- function(a, b) length(intersect(order(-a)[1:10], order(-b)[1:10]))
t10 <- sapply(vars, function(v) top10(pooled[, v], mm[, v]))
hb <- function(M) { b <- colMeans(M); b >= quantile(b, 0.75) }                                    # top quarter of zones by burden
stab <- mean(sapply(1:N_JOINT, function(j) mean(hb(matrix(IND$burden[j, ], 1)) == hb(MAIN$burden))))
cmp <- data.table(measure = c(paste("Spearman, pooled refits vs main:", vars), paste("Spearman, single refits vs main (median, min):", vars),
                              paste("Top-10 zones shared (of 10):", vars), "Zones in the same burden group (top quarter or not), mean over refits",
                              "Baseline self-transition: refits mean (range) | main", "High-state self-transition: refits mean (range) | main",
                              "National high-state burden: refits | main", "R-hat per refit: max | share below 1.1"),
  value = c(sprintf("%.3f", sp), sprintf("%.3f (%.3f)", apply(sp_each, 2, median), apply(sp_each, 2, min)), as.character(t10), sprintf("%.1f%%", 100 * stab),
            sprintf("%.3f (%.3f to %.3f) | %.3f", mean(per_draw$gam11), min(per_draw$gam11), max(per_draw$gam11), mean(g(1, 1))),
            sprintf("%.3f (%.3f to %.3f) | %.3f", mean(per_draw$gam33), min(per_draw$gam33), max(per_draw$gam33), mean(g(3, 3))),
            sprintf("%.4f | %.4f", mean(pooled[, "burden"]), mean(mm[, "burden"])),
            sprintf("%.3f | %.1f%%", max(per_draw$rhat_max), 100 * mean(per_draw$share_rhat_below_1.1))))
print(cmp, row.names = FALSE); fwrite(cmp, file.path(OUTJ, "Table_S_joint20_long_comparison.csv"))
fwrite(data.table(zone = zone_order, setNames(as.data.table(pooled), paste0("refits_", vars)), setNames(as.data.table(mm), paste0("main_", vars)),
                  sd_between_refits_burden = apply(IND$burden, 2, sd)), file.path(OUTJ, "Table_S_joint20_long_indicators_by_zone.csv"))
txt <- c(sprintf("Twenty joint draws of the first-stage BYM2 fields were taken from the Gaussian approximation of the Stage 1 INLA posterior, and the HMM was refitted to each (2 chains of %s iterations, %s burn-in, thinning %d; %s retained draws per refit).",
                 format(BURN + NIT, big.mark = ","), format(BURN, big.mark = ","), THIN, format(NCHD * NIT / THIN, big.mark = ",")),
         sprintf("Across refits, R-hat was below 1.1 for %.1f%% of transition and high-state parameters (maximum %.2f).", 100 * mean(per_draw$share_rhat_below_1.1), max(per_draw$rhat_max)),
         sprintf("Zone indicators from the refits agreed with the main analysis (Spearman: burden %.2f, entry %.2f, persistence %.2f, peak %.2f), and %.0f%% of zones stayed in the same burden group.",
                 sp[1], sp[2], sp[3], sp[4], 100 * stab))
writeLines(txt, file.path(OUTJ, "Joint20_long_sentences.txt")); cat("\n", txt, sep = "\n")
say("Done. Next (Stage 3 on these refits): options(cholera.draws = '", file.path(OUTJ, "Draws_joint20_long_indicators.rds"), "', cholera.ndraws = 20, cholera.out = '", OUTJ, "'); source('Stage2_LPA_Louvain_final.R')")
