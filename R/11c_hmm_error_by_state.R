## In-sample prediction error of the final HMM, overall and by state (modal state of each zone-month)
suppressPackageStartupMessages({ library(coda); library(data.table) })
set.seed(2026)
proj <- getOption("cholera.proj", ".")
S0  <- readRDS(file.path(proj, "INLA_Rebuilt_From_Zero_OCCPOP.rds")); dat <- as.data.table(S0$data); setorder(dat, ID_area, ID_time)
Z <- length(S0$zone_order); T <- max(dat$ID_time)
Y   <- matrix(dat$O, Z, T, byrow = TRUE)
Em  <- matrix(dat$E, Z, T, byrow = TRUE); OFF <- Em / mean(Em)
X   <- as.matrix(readRDS(file.path(proj, "improve_hmm", "cache_tau2_2_hi170_rmin0.5_hm_main_n100000_b200000.rds")))
scol <- sprintf("S[%d,%d]", rep(1:Z, T), rep(1:T, each = Z))          # column order: zone fastest, then month
mcol <- sprintf("mu[%d,%d]", rep(1:Z, 3), rep(1:3, each = Z)); dcol <- sprintf("delta[%d,%d]", rep(1:Z, 3), rep(1:3, each = Z))
kk <- sample(nrow(X), 1000)                                            # 1,000 posterior draws
LAM <- matrix(0, Z, T); CNT <- array(0, c(Z, T, 3))
for (k in kk) {
  A  <- matrix(X[k, scol], Z, T)
  ls <- matrix(X[k, mcol] * exp(X[k, dcol]), Z, 3)
  LAM <- LAM + matrix(ls[cbind(rep(1:Z, T), as.vector(A))], Z, T) * OFF / length(kk)
  for (s in 1:3) CNT[, , s] <- CNT[, , s] + (A == s)
}
modal <- apply(CNT, c(1, 2), which.max)
st <- c("Baseline", "Moderate", "High-intensity")
row <- function(lab, i) data.table(state = lab, n = sum(i), observed_mean = mean(Y[i]), predicted_mean = mean(LAM[i]),
  MAE = mean(abs(Y[i] - LAM[i])), RMSE = sqrt(mean((Y[i] - LAM[i])^2)), bias = mean(LAM[i] - Y[i]))
tab <- rbind(rbindlist(lapply(1:3, function(s) row(st[s], modal == s))), row("All zone-months", modal > 0))
tab[, (3:7) := lapply(.SD, round, 2), .SDcols = 3:7]
print(tab, row.names = FALSE)
cat(sprintf("\nTotal predicted %.0f vs observed %d cases\n", sum(LAM), sum(Y)))
fwrite(tab, file.path(proj, "improve_hmm", "Table_S_HMM_error_by_state_FINAL.csv"))
