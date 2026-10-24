# ==============================================================================
# Joint_stage1_propagation.R  -  carries the JOINT first-stage uncertainty into the HMM
# Idea: instead of giving the HMM only the mean and SD of each zone's first-stage spatial effect, draw N complete sets of
# BYM2 spatial fields (all 92 zones together, so the correlation between zones is kept) from the two first-stage models,
# fit one shorter HMM to each set (the field is held fixed at the draw), and look at how much the zone indicators move.
# Same HMM as the final run (version F: prior centre 170, shared national high-state gap, r >= 0.5).
# NOT RUN by its author (no R in his sandbox): expect a small first-run fix.
#
# Needs: R packages INLA, rjags, coda, data.table, parallel; JAGS 4.3.2.
#        The Stage 1 file INLA_Rebuilt_From_Zero_OCCPOP.rds with the two main fits saved with config = TRUE
#        (check first:  S <- readRDS(".../INLA_Rebuilt_From_Zero_OCCPOP.rds"); !is.null(S$fits$occ$misc$configs); !is.null(S$fits$int$misc$configs) ).
#        If either is NULL, the two main models must be refitted with  control.compute = list(config = TRUE)  before this can run.
# Run AFTER the long final run has finished (it needs the cores):
#        cd ~/Documents/spatial_anlysis/update_HMM ; nohup setsid Rscript Joint_stage1_propagation.R > joint_log.txt 2>&1 &
#        tail -f joint_log.txt
# Time: one short run takes about 40 minutes; with 10 workers and 20 draws about 1.5 to 2 hours. Pauses at 80 C are NOT built in: check sensors.
# Output (folder improve_hmm/joint_stage1): one RDS per run (resumable), Table_S_joint_propagation_*.csv, Draws_joint_indicators.rds, Joint_sentences.txt
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(parallel) })
proj   <- getOption("cholera.proj", ".")
S1FILE <- getOption("cholera.stage1", file.path(proj, "INLA_Rebuilt_From_Zero_OCCPOP.rds"))
out    <- getOption("cholera.out", file.path(proj, "improve_hmm")); OUTJ <- file.path(out, "joint_stage1"); dir.create(OUTJ, showWarnings = FALSE, recursive = TRUE)
MAIN_DRAWS <- getOption("cholera.draws", file.path(out, "Draws_updated_indicators_tau2_2_hi170_rmin0.5_hm_main_n1091_b2000_FINAL.rds"))   # main run, for comparison (edit to the real name)
N_JOINT <- getOption("cholera.njoint", 20); NW <- getOption("cholera.cores", 10)      # number of first-stage draws; parallel workers (one run per worker)
BURN <- 5000; NIT <- 5000; THIN <- 5                                                   # per short run: 1,000 kept draws
set.seed(2026)
say <- function(...) cat(format(Sys.time(), "%d %b %H:%M:%S"), ..., "\n")

# ---- data (identical to HMM_check_ALL.R)
S0 <- readRDS(S1FILE); dat <- as.data.table(S0$data); setorder(dat, ID_area, ID_time)
zone_order <- S0$zone_order; Z <- length(zone_order); T <- max(dat$ID_time); stopifnot(nrow(dat) == Z * T)
Y <- matrix(dat$O, Z, T, byrow = TRUE); Em <- matrix(dat$E, Z, T, byrow = TRUE); OFF <- Em / mean(Em)
have <- requireNamespace("INLA", quietly = TRUE) && !is.null(S0$fits$occ$misc$configs) && !is.null(S0$fits$int$misc$configs)
if (!have) stop("INLA is not installed, or the two main Stage 1 fits were not saved with config = TRUE. See the header of this script.")
thr <- as.numeric(quantile(Y[Y > 0], 0.95)); anch2 <- median(Y[Y > 0 & Y < thr]); anch1 <- 0.01; ANCHHI <- 170
st_emp <- matrix(1L, Z, T); st_emp[Y > 0 & Y < thr] <- 2L; st_emp[Y >= thr] <- 3L
trm <- matrix(0, 3, 3); for (z in 1:Z) for (t in 2:T) { a <- st_emp[z, t - 1]; b <- st_emp[z, t]; trm[a, b] <- trm[a, b] + 1 }
self_emp <- diag(trm) / rowSums(trm); CONC <- 10
a_self <- pmax(self_emp, 0.05) * CONC; a_self[3] <- a_self[2]; a_oth <- (1 - c(self_emp[1], self_emp[2], self_emp[2])) * CONC / 2
alpha <- matrix(0, 3, 3); for (i in 1:3) { alpha[i, ] <- a_oth[i]; alpha[i, i] <- a_self[i] }
say(sprintf("Anchors: moderate %.0f, upper threshold %.0f (high-state prior centre %d); self-transition rates %s; transitions out of high: %d",
            anch2, thr, ANCHHI, paste(round(self_emp, 3), collapse = " / "), sum(trm[3, ])))

# ---- the same JAGS model as the final run (hierarchical high-state gap, no extensions)
model_txt <- "
model {
  pi0[1:3] ~ ddirch(c(8, 1, 1))
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
      lam[z, t] <- max(lam_s[z, S[z, t]] * OFF[z, t], 1e-6)
      p[z, t] <- r[z] / (r[z] + lam[z, t])
      Y[z, t] ~ dnegbin(p[z, t], r[z])
    }
  }
}"
mf <- file.path(OUTJ, "hmm_model_joint.jags"); writeLines(model_txt, mf)

# ---- N_JOINT joint draws of the two BYM2 fields (all zones at once), centred the same way as the main run
smp <- function(fit, J, sd0) { x <- INLA::inla.posterior.sample(J, fit, selection = list(ID_area_bym = 1:Z), seed = sd0); sapply(x, function(s) as.numeric(s$latent[, 1])[1:Z]) }
DO <- smp(S0$fits$occ, N_JOINT, 2026L); DI <- smp(S0$fits$int, N_JOINT, 2027L)               # Z x N_JOINT; the two models are treated as independent
cen <- function(M) { M <- sweep(M, 2, colMeans(M)); sweep(M, 2, log(colMeans(exp(M)))) }     # each draw: mean zero, then average multiplier exp(delta) = 1
DO <- cen(DO); DI <- cen(DI)
stopifnot(all(dim(DO) == c(Z, N_JOINT)), all(is.finite(DO)), all(is.finite(DI)))
fo <- S0$fits$occ$summary.random$ID_area_bym$mean[1:Z]; fi <- S0$fits$int$summary.random$ID_area_bym$mean[1:Z]
say(sprintf("Mean of the sampled fields vs the first-stage posterior mean (Pearson): occurrence %.2f, intensity %.2f (should be high, above 0.5)", cor(rowMeans(DO), fo), cor(rowMeans(DI), fi)))

jdata <- list(Y = Y, OFF = OFF, Z = Z, T = T, alpha = alpha, anch1 = anch1, anch2 = anch2, anch3 = ANCHHI, tau_a = 2, tau_b = 2, prec_eta = 1, rmin = 0.5)

# ---- one short HMM per draw, on parallel workers; each result is saved the moment it is done (a stop loses only the unfinished runs)
one_run <- function(j) {
  f <- file.path(OUTJ, sprintf("run%02d.rds", j)); if (file.exists(f)) { r <- tryCatch(readRDS(f), error = function(e) NULL); if (!is.null(r)) return(r) }
  suppressPackageStartupMessages({ library(rjags); library(coda) })
  jd <- jdata; jd$d_occ <- DO[, j]; jd$d_int <- DI[, j]; jd$prec_occ <- rep(1e4, Z); jd$prec_int <- rep(1e4, Z)   # field fixed at the draw (SD 0.01): its uncertainty is carried by the draws
  m <- jags.model(mf, data = jd, inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 5000 + j, S = st_emp), n.chains = 1, quiet = TRUE, n.adapt = 1000); update(m, BURN)
  x <- as.matrix(coda.samples(m, c("S", "gam", "mu", "delta"), n.iter = NIT, thin = THIN)[[1]]); nd <- nrow(x)
  A <- array(x[, as.vector(outer(1:Z, 1:T, function(z, t) sprintf("S[%d,%d]", z, t)))], c(nd, Z, T))
  g <- function(i, k) x[, sprintf("gam[%d,%d,%d]", 1:Z, i, k)]
  lvl <- function(s) x[, sprintf("mu[%d,%d]", 1:Z, s)] * exp(x[, sprintf("delta[%d,%d]", 1:Z, s)])
  mx <- apply(A, c(1, 2), max)
  act_prev <- A[, , -T, drop = FALSE] >= 2; act_now <- A[, , -1, drop = FALSE] >= 2
  r <- list(burden = colMeans(apply(A == 3, c(1, 2), mean)), entry = colMeans(g(1, 2) + g(1, 3)),
            persistence_active = colMeans(apply(act_prev & act_now, c(1, 2), sum) / pmax(apply(act_prev, c(1, 2), sum), 1)),
            peak_reached = colMeans(lvl(3) * (mx == 3) + lvl(2) * (mx == 2) + lvl(1) * (mx == 1)),
            gam11 = mean(g(1, 1)), gam33 = mean(g(3, 3)))
  saveRDS(r, paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f); r }
say(sprintf("Running %d short HMMs (%d iterations burn-in + %d kept, thin %d) on %d workers; results in %s", N_JOINT, BURN, NIT, THIN, NW, OUTJ))
cl <- makeCluster(min(NW, N_JOINT), port = 21000L + as.integer(Sys.getpid()) %% 15000L)
clusterExport(cl, c("one_run", "OUTJ", "jdata", "DO", "DI", "mf", "st_emp", "Z", "T", "BURN", "NIT", "THIN"), envir = environment())
RES <- parLapplyLB(cl, 1:N_JOINT, one_run); stopCluster(cl)
say("All runs finished.")

# ---- summarise
vars <- c("burden", "entry", "persistence_active", "peak_reached")
IND <- lapply(setNames(vars, vars), function(v) t(sapply(RES, `[[`, v)))                        # N_JOINT x Z: zone indicators from each run
chk <- all(sapply(IND, function(M) all(is.finite(M)) && ncol(M) == Z)); say("All runs returned complete indicators:", chk)
saveRDS(list(ind = IND, zone_order = zone_order), file.path(OUTJ, "Draws_joint_indicators.rds"))   # same layout as the main draws file: usable by Stage2_LPA_Louvain_final.R
pooled <- sapply(IND, colMeans)
tab <- data.table(zone = zone_order, pooled, between_run_sd_burden = apply(IND$burden, 2, sd))
txt <- c(sprintf("Joint first-stage uncertainty: %d draws of the whole BYM2 fields from the occurrence and intensity models (treated as independent), one HMM of %d retained draws fitted to each.", N_JOINT, NIT / THIN))
if (file.exists(MAIN_DRAWS)) {
  M <- readRDS(MAIN_DRAWS)$ind; mm <- sapply(vars, function(v) colMeans(M[[v]]))
  sp <- sapply(vars, function(v) cor(pooled[, v], mm[, v], method = "spearman")); names(sp) <- vars
  say("Spearman of the pooled joint-run indicators with the main run:", paste(sprintf("%s %.2f", vars, sp), collapse = " | "))
  rel <- median(apply(IND$burden, 2, sd)) / median(pooled[, "burden"])
  say(sprintf("Between-run SD of zone burden (median over zones) %.4f against a median burden of %.4f", median(apply(IND$burden, 2, sd)), median(pooled[, "burden"])))
  tab <- cbind(tab, setNames(as.data.frame(mm), paste0("main_", vars)))
  txt <- c(txt, sprintf("The zone indicators from the joint runs agreed with those from the main run (Spearman %s).", paste(sprintf("%s %.2f", vars, sp), collapse = ", ")))
} else say("Main draws file not found: set options(cholera.draws = '...') to compare with the main run.")
fwrite(tab, file.path(OUTJ, "Table_S_joint_propagation_indicators.csv")); writeLines(txt, file.path(OUTJ, "Joint_sentences.txt"))
cat("\nNext: stability of profiles, communities and classes across these runs:\n",
    "  options(cholera.draws = '", file.path(OUTJ, "Draws_joint_indicators.rds"), "', cholera.ndraws = ", N_JOINT, ", cholera.out = '", OUTJ, "'); source('Stage2_LPA_Louvain_final.R')\n", sep = "")
say("Done.")
