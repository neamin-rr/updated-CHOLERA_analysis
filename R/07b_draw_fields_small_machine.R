## Run on the WINDOWS LAPTOP (has INLA). Draws 20 joint BYM2 fields from the two Baseline Stage 1 fits and saves them to one small file.
## Then copy  F:/Ub/Paper_one_revision/joint_fields_20.rds  to the Linux PC folder  ~/Documents/spatial_anlysis/update_HMM/
library(INLA)
S0 <- readRDS(STAGE1_RDS)
Z <- length(S0$zone_order); N_JOINT <- 20
stopifnot(!is.null(S0$fits$occ$misc$configs), !is.null(S0$fits$int$misc$configs))
smp <- function(fit, J, sd0) { x <- inla.posterior.sample(J, fit, selection = list(ID_area_bym = 1:Z), seed = sd0)
  sapply(x, function(s) as.numeric(s$latent[, 1])[1:Z]) }
DO <- smp(S0$fits$occ, N_JOINT, 2026L); DI <- smp(S0$fits$int, N_JOINT, 2027L)          # Z x 20, all zones drawn together
cen <- function(M) { M <- sweep(M, 2, colMeans(M)); sweep(M, 2, log(colMeans(exp(M)))) } # same centring as the main HMM run
DO <- cen(DO); DI <- cen(DI)
stopifnot(all(dim(DO) == c(Z, N_JOINT)), all(is.finite(DO)), all(is.finite(DI)))
fo <- S0$fits$occ$summary.random$ID_area_bym$mean[1:Z]; fi <- S0$fits$int$summary.random$ID_area_bym$mean[1:Z]
cat(sprintf("Agreement of the drawn fields with the Stage 1 means (Pearson): occurrence %.2f, intensity %.2f (should be above 0.5)\n",
            cor(rowMeans(DO), fo), cor(rowMeans(DI), fi)))
saveRDS(list(DO = DO, DI = DI, zone_order = S0$zone_order, made = Sys.time()), file.path(PROJ))
cat("Saved", JOINT_RDS, "- copy it to the machine that runs the HMM stages\n")
