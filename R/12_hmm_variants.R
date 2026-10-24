# ==============================================================================
# Run_all_variants.R : runs the HMM small check for many model versions AT THE SAME TIME on a many-core machine,
# then puts every check of every version side by side in one table.
# Use:   cd <folder with HMM_check_ALL.R and the .rds>
#        nohup Rscript Run_all_variants.R > variants_log.txt 2>&1 &
# Step 0 is a 1-minute code test; if it fails nothing else is started.
# ==============================================================================
proj   <- getOption("cholera.proj", ".")
script <- file.path(proj, "HMM_check_ALL.R"); out <- file.path(proj, "improve_hmm"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
setwd(proj)
stopifnot(file.exists(script), file.exists(file.path(proj, "INLA_Rebuilt_From_Zero_OCCPOP.rds")))
library(parallel)

# ---- the versions to test (one line each; every version runs 4 chains on 4 cores) --------------------------------
WIDE <- "cholera.tau23 = c(2, 2)"
VAR <- list(
  A_tight_prior_as_before        = "cholera.tau23 = c(10, 3), cholera.anchhi = NA, cholera.hiermean = FALSE",
  B_wide_prior                   = paste(WIDE, ", cholera.anchhi = NA, cholera.hiermean = FALSE"),
  C_wide_learned_level           = paste(WIDE, ", cholera.anchhi = NA, cholera.hiermean = TRUE"),
  D_learned_level_rmin05         = paste(WIDE, ", cholera.anchhi = NA, cholera.hiermean = TRUE, cholera.rmin = 0.5"),
  E_centred_170                  = paste(WIDE, ", cholera.anchhi = 170, cholera.hiermean = TRUE"),
  F_centred_170_rmin05           = paste(WIDE, ", cholera.anchhi = 170, cholera.hiermean = TRUE, cholera.rmin = 0.5"),
  G_learned_pooled_rows_k30      = paste(WIDE, ", cholera.anchhi = NA, cholera.hiermean = TRUE, cholera.rmin = 0.5, cholera.kappa = 30"),
  H_learned_season               = paste(WIDE, ", cholera.anchhi = NA, cholera.hiermean = TRUE, cholera.rmin = 0.5, cholera.season = TRUE"),
  I_learned_neighbours           = paste(WIDE, ", cholera.anchhi = NA, cholera.hiermean = TRUE, cholera.rmin = 0.5, cholera.nb = TRUE"),
  J_learned_population           = paste(WIDE, ", cholera.anchhi = NA, cholera.hiermean = TRUE, cholera.rmin = 0.5, cholera.pop = TRUE"),
  K_everything_on                = paste(WIDE, ", cholera.anchhi = NA, cholera.hiermean = TRUE, cholera.rmin = 0.5, cholera.kappa = 30, cholera.season = TRUE, cholera.nb = TRUE, cholera.pop = TRUE"),
  L_very_wide_centred_170        = "cholera.tau23 = c(1, 1), cholera.anchhi = 170, cholera.hiermean = FALSE",
  M_learned_season_and_neighbours = paste(WIDE, ", cholera.anchhi = NA, cholera.hiermean = TRUE, cholera.rmin = 0.5, cholera.season = TRUE, cholera.nb = TRUE"),
  N_learned_rmin1                = paste(WIDE, ", cholera.anchhi = NA, cholera.hiermean = TRUE, cholera.rmin = 1"),
  O_centred_250_rmin05           = paste(WIDE, ", cholera.anchhi = 250, cholera.hiermean = TRUE, cholera.rmin = 0.5"))
ONLY <- Sys.getenv("ONLY")   # for example ONLY=A,D,F runs just those versions (the letter at the start of each name)
if (nzchar(ONLY)) VAR <- VAR[substr(names(VAR), 1, 1) %in% trimws(strsplit(ONLY, ",")[[1]])]
stopifnot(length(VAR) > 0)
CH <- as.integer(Sys.getenv("CHAINS", 4))   # chains (= cores) per version; CHAINS=20 lets ONE version use 20 cores
HOLDOUT <- TRUE        # FALSE saves about a third of the time and skips checks 9 to 11

run_one <- function(name, opts, smoke = FALSE, holdout = HOLDOUT) {
  tagfile <- file.path(out, paste0("tag_", name, ".txt")); unlink(tagfile)
  if (!smoke && exists("NCONC")) Sys.sleep(10 * ((match(name, names(VAR)) - 1) %% NCONC))   # stagger the starts by 10 s, so the versions do not all read the 344 MB data file at the same moment
  expr <- sprintf("options(%s, cholera.cores = %d, cholera.nch = %d, cholera.holdout = %s, cholera.smoke = %s); source('%s')", opts, CH, CH, holdout && !smoke, smoke, script)
  log <- file.path(out, paste0("run_log_", name, if (smoke) "_SMOKE", ".txt"))
  st <- system2("Rscript", c("-e", shQuote(expr)), stdout = log, stderr = log,
                env = c(paste0("TAGFILE=", tagfile), "R_DATATABLE_NUM_THREADS=1", "OMP_NUM_THREADS=1"))
  cat(sprintf("[%s] %s finished with status %s\n", format(Sys.time(), "%H:%M:%S"), name, st))
  st }

# ---- step 0: code test (about 1 minute) ------------------------------------------------------------------------
unlink(list.files(out, pattern = "^(run_log_[A-O]_|tag_[A-O]_)", full.names = TRUE))   # remove logs left over from earlier attempts, so the progress check only shows this run
if (Sys.getenv("SKIPTEST") == "1") { cat("Step 0 skipped (SKIPTEST=1)\n"); s0 <- 0L } else { cat("Step 0: code test\n"); s0 <- run_one("codetest", VAR[[3]], smoke = TRUE, holdout = FALSE) }
if (!identical(as.integer(s0), 0L)) { cat("CODE TEST FAILED. Last lines of its log:\n"); print(tail(readLines(file.path(out, "run_log_codetest_SMOKE.txt")), 25)); stop("Fix this first; nothing else was started.") }

# ---- step 1: all versions at once -------------------------------------------------------------------------------
ncore <- detectCores(); ram <- suppressWarnings(as.numeric(system("free -g | awk '/Mem:/{print $2}'", intern = TRUE))); if (!length(ram) || is.na(ram)) ram <- 32
CORES_TOTAL <- as.integer(Sys.getenv("CORES_TOTAL", 60))   # cores this launcher may use (60 of your 64); each version needs 4
NCONC <- as.integer(Sys.getenv("NCONC", NA)); if (is.na(NCONC)) NCONC <- max(1, min(floor(min(CORES_TOTAL, ncore) / CH), floor(ram / 8), length(VAR)))   # 4 cores per version; about 8 GB of memory per version
cat(sprintf("Step 1: %d versions, %d at a time = %d cores in use (machine: %d cores, %.0f GB memory)\n", length(VAR), NCONC, CH * NCONC, ncore, ram))
t0 <- Sys.time()
cl <- makeCluster(min(NCONC, length(VAR)))                 # PSOCK: works on Windows
on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
clusterExport(cl, c("run_one", "VAR"), envir = environment())
st <- parLapply(cl, names(VAR), function(n) try(run_one(n, VAR[[n]]), silent = TRUE))
stopCluster(cl)
cat(sprintf("All versions done in %.1f minutes\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---- step 2: side-by-side table of every check ---------------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
rows <- list(); det <- list()
for (n in names(VAR)) {
  tf <- file.path(out, paste0("tag_", n, ".txt"))
  if (!file.exists(tf)) { cat("No result for", n, "(see run_log_", n, ".txt)\n", sep = ""); next }
  tag <- readLines(tf)[1]; f <- file.path(out, sprintf("Table_S_HMM_check_status_%s.csv", tag))
  if (!file.exists(f)) { cat("No check table for", n, "\n"); next }
  d <- fread(f); d[, variant := n]; d[, tag := tag]; rows[[n]] <- d }
if (!length(rows)) stop("No variant finished. Open the run_log files in ", out)
A <- rbindlist(rows, fill = TRUE)
A[, id := factor(id, levels = unique(c(as.character(1:25), "5b")))]
M <- dcast(A, id ~ variant, value.var = "status", fun.aggregate = function(x) substr(x[1], 1, 1))   # P / F / S
M <- M[order(as.integer(sub("b", "", as.character(id))), as.character(id))]
num <- function(s, pat) suppressWarnings(as.numeric(sub(pat, "\\1", s)))
summ <- A[, .(passed = sum(status == "PASS"), failed = sum(status == "FAIL"),
              failed_excl_check1 = sum(status == "FAIL" & id != "1"),
              failing = paste(id[status == "FAIL"], collapse = ","),
              high_state_ratio = num(detail[id == "18"][1], "^ratio ([0-9.]+).*"),
              residual_sd = num(detail[id == "19"][1], "^([0-9.]+).*")), by = variant][order(failed_excl_check1, failed)]
cat("\n==================== ALL VERSIONS SIDE BY SIDE (P = pass, F = fail, S = skipped) ====================\n"); print(M, row.names = FALSE)
cat("\n==================== RANKING (check 1 is the convergence check; a small run can fail it for lack of iterations, so the ranking ignores it) ====================\n"); print(summ, row.names = FALSE)
fwrite(M, file.path(out, "Table_S_HMM_variant_checks_side_by_side.csv")); fwrite(summ, file.path(out, "Table_S_HMM_variant_ranking.csv"))
best <- summ$variant[1]; cat(sprintf("\nFewest failed checks (excluding check 1): %s (failing: %s)\n", best, summ$failing[1]))
cat("Look at it with the high-state ratio (target 1.00 +/- 0.10) and residual SD (target 0.8 to 1.25). Send me Table_S_HMM_variant_ranking.csv and the side-by-side table.\n")
jf <- file.path(out, sprintf("Table_S_HMM_juan_comment_coverage_%s.csv", A[variant == best, tag][1]))
if (file.exists(jf)) { cat("\n==================== JUAN'S COMMENTS for the top-ranked version ====================\n"); print(fread(jf)[, .(comment = substr(comment, 1, 60), result_here, evidence)], right = FALSE) }
