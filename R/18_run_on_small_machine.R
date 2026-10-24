# ==============================================================================
# Run_on_my_PC.R : runs ONE HMM version on your Windows PC (RStudio or R console).
# Needs in the folder below:  HMM_check_ALL.R  and  INLA_Rebuilt_From_Zero_OCCPOP.rds
# Use: open this file in RStudio, change VERSION if you want, then Source it (or select all and Run).
# It uses 4 cores (one per chain) and takes a few hours. Results go to <folder>/improve_hmm.
# ==============================================================================
proj <- file.path(PROJ)          # <- folder that holds the two files above
VERSION <- "D"                                # D = learned level (wide prior, rmin 0.5); F = centred on 170; A = old tight prior
SMOKE   <- FALSE                              # TRUE = 2-minute test of the code only (results meaningless)

stopifnot(file.exists(file.path(proj, "HMM_check_ALL.R")), file.exists(file.path(proj, "INLA_Rebuilt_From_Zero_OCCPOP.rds")))
for (p in c("data.table", "rjags", "coda", "Matrix")) if (!requireNamespace(p, quietly = TRUE)) stop("Install the package first: install.packages('", p, "')   (rjags also needs JAGS from https://mcmc-jags.sourceforge.io)")
setwd(proj)
WIDE <- c(2, 2)
opt <- switch(VERSION,
  A = list(cholera.tau23 = c(10, 3), cholera.anchhi = NA,  cholera.hiermean = FALSE),
  D = list(cholera.tau23 = WIDE,     cholera.anchhi = NA,  cholera.hiermean = TRUE, cholera.rmin = 0.5),
  F = list(cholera.tau23 = WIDE,     cholera.anchhi = 170, cholera.hiermean = TRUE, cholera.rmin = 0.5),
  stop("VERSION must be A, D or F"))
do.call(options, c(opt, list(cholera.proj = proj, cholera.cores = 4, cholera.nch = 4, cholera.holdout = !SMOKE, cholera.smoke = SMOKE)))
cat("Version", VERSION, "| folder", proj, "| smoke test:", SMOKE, "| started", format(Sys.time(), "%H:%M:%S"), "\n")
source(file.path(proj, "HMM_check_ALL.R"))
cat("Finished", format(Sys.time(), "%H:%M:%S"), "- results are in", file.path(proj, "improve_hmm"), "\n")
