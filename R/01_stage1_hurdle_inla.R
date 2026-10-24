# ==============================================================================
# CHOLERA ETHIOPIA - ALL CODE IN ONE FILE  (run top to bottom in ONE R session: Source / Ctrl+Shift+Enter)
#
# PART 1  Main pipeline (Steps 0-11): raw-data checks, prior grid (81 + 243), hurdle models (zero-truncated NB intensity),
#         joint association model, waves, checks, persistence, gt tables, maps.
# PART 2  Extras that use the fitted models in memory: residual variograms, seasonal tables (absolute + by year),
#         descriptive Figure 1 / Figure 2 / Table 1, numbers for the yellow [XX] placeholders.
# PART 3  (not in this file) Stage 2 is run with HMM_check_ALL.R and Stage2_LPA_Louvain_final.R
#
# WHAT IS RELOADED AND WHAT IS REFITTED (nothing is refitted that is still valid)
#   * Every fitted model is cached in a folder whose name contains its settings:
#       occurrence : fits_OCCPOP_rw20      (RW2 bound U_RW_OCC = 20, set in STEP 1)   -> NEW folder: refitted the first time
#       intensity  : fits_INTPOP_rw10_ZT   (zero-truncated NB)                          -> reloaded if present ("loaded from file (date)")
#   * A saved fit is reloaded only if its likelihood and number of rows match the current setup; otherwise it is refitted.
#   * The intensity prior grid (243) is reused; the occurrence grid (81) is rerun once because U_RW_OCC changed (old files are renamed _OLD_<date>).
#   * Old results written before this run started are ignored by the gt tables (stale-file guard).
#   * To force a full refit of everything: REFIT_ALL <- TRUE in STEP 1.
#
# SWITCHES (survive the rm(list = ls()) inside STEP 1 because they are R options)
# ==============================================================================
options(cholera.run_part1 = TRUE,    # FALSE = skip the main pipeline (only if fit_o, fit_i, dat ... are already in memory)
        cholera.run_part2 = TRUE,    # variograms, seasonal tables, descriptive figures, yellow-item numbers
        cholera.run_part3 = FALSE)   # Stage 2 (HMM / LPA / Louvain) is NOT in this file: use HMM_check_ALL.R and Stage2_LPA_Louvain_final.R
RUN <- function(k) isTRUE(getOption(paste0("cholera.run_", k), FALSE))

# ##############################################################################
# PART 1 - MAIN PIPELINE
# ##############################################################################
if (RUN("part1")) {
  # ==============================================================================
  # CHOLERA ETHIOPIA — COMPLETE PIPELINE, FROM RAW-DATA CHECKS TO MODELLING AND CHECKS
  # One file, run top to bottom in ONE R session (Source / Ctrl+Shift+Enter).
  #   STEP 0  raw-data checks (no INLA, seconds)
  #   STEP 1  modelling: data build, PRIOR GRID (81 occurrence + 243 intensity Baseline fits with PC priors; the best set is selected and used),  E and population, hurdle models (occurrence + intensity),
  #           model selection, joint model, waves, calibration, decoupling checks D1-D5,
  #           silent months and persistence, expected-caseload waves, summary report
  #   STEP 2  spatial autocorrelation and clusters (Moran's I, LISA, Gi*)
  #   STEP 3  pre-deployment checks (fit health, priors, zero-truncation, predictive checks,
  #           residuals, out-of-time validation, rank stability)
  #   STEP 3B overfitting check: Baseline vs interaction models, in-sample vs held-out (random 20% and forward 2022-2024)
  #   STEP 3C PIT (calibration) histograms and summary for all five candidate models, both parts
  #   STEP 7  structure checks: hurdle vs zero-inflated models; population adjustment / spatial confounding (OCC_POP sensitivity)
  #   STEP 5  complementary analyses: STL decomposition, month effects (NB mixed), Gi* 90/95/99%, observed vs predicted + uncertainty map
  #   STEP 8  persistence extension (predictive model): forward one-step-ahead validation, residual autocorrelation, effect on zone results
  #   STEP 9  maps; zones without intensity data are grey
  #   STEP 10 manuscript tables with gt: each shown in the Viewer and saved as .docx (SAVE_TABLES)
  #   STEP 11 publication maps (shown in the Plots pane; saved only if SAVE_MAPS <- TRUE; also available as the stand-alone Stage_maps_publication.R)
  # Settings (OCC_POP, INT_POP, U_RW_OCC, U_RW_INT, RUN_* flags, NDRAW) are at the top of STEP 1.
  # FIGURES: drawn in the RStudio Plots pane only; nothing is saved unless SAVE_FIGS <- TRUE (STEP 1 settings).
  # Fitted models are cached in fits*/ folders, so a rerun after a crash is fast.
  # Edit the three paths (data_path, shp_path, output_path) if your folders differ.
  # ==============================================================================
  
  
  ################################################################################
  # STEP 0 — RAW-DATA CHECKS
  ################################################################################
  # ==============================================================================
  # STAGE 0 — RAW DATA CHECKS (run BEFORE the models; takes seconds, no INLA)
  # Answers: (1) is "monthly" population really one value per year?
  #          (2) zero vs non-reporting, missing cells, duplicates, outliers
  #          (3) population jumps (boundary changes), cases > population
  # Everything is printed AND saved as CSV in <output_path>/raw_checks/
  # ==============================================================================
  library(data.table)
  
  data_path   <- "F:/phd/DATA_source/correted_E_inala/Cholera_updated_data.Rdata"
  output_path <- getOption("cholera.proj", ".")
  chk_dir <- file.path(output_path, "raw_checks"); dir.create(chk_dir, showWarnings = FALSE, recursive = TRUE)
  save_chk <- function(x, name) fwrite(x, file.path(chk_dir, paste0(name, ".csv")))
  
  load(data_path)
  raw <- as.data.table(final_comprehensive_dataset_corrected)
  if ("geometry" %in% names(raw)) raw[, geometry := NULL]
  raw[, zone := toupper(gsub("[^A-Z0-9]", "", admin2Pcode))]
  cat("Columns available:\n"); print(names(raw))
  cat("\nRows:", nrow(raw), "| zones:", uniqueN(raw$zone), "| years:",
      paste(range(raw$year, na.rm = TRUE), collapse = "-"), "\n")
  
  # ------------------------------------------------------------------------------
  # 1. STRUCTURE: duplicates, missing cells, NA, impossible values
  # ------------------------------------------------------------------------------
  cat("\n=== 1. STRUCTURE ===\n")
  cat("Duplicate zone-month rows:", sum(duplicated(raw[, .(zone, year, month)])), "\n")
  cat("NA total_cases:", sum(is.na(raw$total_cases)),
      "| NA monthly_population:", sum(is.na(raw$monthly_population)),
      "| NA total_population:", if ("total_population" %in% names(raw)) sum(is.na(raw$total_population)) else "n/a", "\n")
  cat("Negative cases:", sum(raw$total_cases < 0, na.rm = TRUE),
      "| non-integer cases:", sum(raw$total_cases != round(raw$total_cases), na.rm = TRUE),
      "| population <= 0:", sum(raw$monthly_population <= 0, na.rm = TRUE), "\n")
  grid <- CJ(zone = unique(raw$zone), year = sort(unique(raw$year)), month = 1:12)
  missing_cells <- grid[!raw, on = .(zone, year, month)]
  cat("Cells missing from the zone x month grid:", nrow(missing_cells), "\n")
  save_chk(missing_cells, "missing_cells")
  
  # ------------------------------------------------------------------------------
  # 2. IS "MONTHLY" POPULATION ONE VALUE PER YEAR?
  # ------------------------------------------------------------------------------
  cat("\n=== 2. POPULATION: monthly vs yearly ===\n")
  d <- raw[!is.na(monthly_population)]
  py <- d[, .(n_distinct = uniqueN(monthly_population),
              pop_min = min(monthly_population), pop_max = max(monthly_population)), by = .(zone, year)]
  cat(sprintf("Zone-years with ONE population value for all 12 months: %d of %d (%.1f%%)\n",
              sum(py$n_distinct == 1), nrow(py), 100 * mean(py$n_distinct == 1)))
  if (mean(py$n_distinct == 1) > 0.95) {
    cat("  -> monthly_population is the SAME in every month of a year (an annual figure repeated).\n",
        "     It steps up each January. Check the size of the steps below.\n")
  } else cat("  -> population changes within the year (a true monthly series or interpolated).\n")
  if ("total_population" %in% names(d)) {
    cat(sprintf("monthly_population == total_population in %.1f%% of rows\n",
                100 * mean(d$monthly_population == d$total_population)))
  }
  
  # January step sizes: change from December (year y-1) to January (year y)
  setorder(d, zone, year, month)
  d[, pop_prev := shift(monthly_population), by = zone]
  d[, chg := monthly_population / pop_prev - 1]
  jan <- d[month == 1 & !is.na(chg)]
  dec_within <- d[month != 1 & !is.na(chg)]
  cat(sprintf("Median change Dec->Jan: %.2f%% | median change within year: %.2f%%\n",
              100 * median(jan$chg), 100 * median(dec_within$chg)))
  cat("Distribution of Dec->Jan change (%):\n"); print(round(100 * quantile(jan$chg, c(0, .01, .05, .25, .5, .75, .95, .99, 1)), 2))
  jumps <- d[!is.na(chg) & abs(chg) > 0.20, .(zone, year, month, pop_prev, monthly_population, change_pct = round(100 * chg, 1))]
  cat("Month-to-month population changes > 20% (possible boundary changes / data errors):", nrow(jumps), "\n")
  if (nrow(jumps)) print(head(jumps[order(-abs(change_pct))], 20))
  save_chk(jumps, "population_jumps_over20pct")
  save_chk(py, "population_values_per_zone_year")
  
  # national population by year (sanity: should rise smoothly ~2.5-3% per year)
  natpop <- d[month == 6, .(national_population = sum(monthly_population)), by = year][order(year)]
  natpop[, growth_pct := round(100 * (national_population / shift(national_population) - 1), 2)]
  cat("\nNational population by year (June):\n"); print(natpop)
  save_chk(natpop, "national_population_by_year")
  
  # ------------------------------------------------------------------------------
  # 3. ZEROS vs NON-REPORTING
  # ------------------------------------------------------------------------------
  cat("\n=== 3. ZEROS ===\n")
  raw[, O := total_cases]
  cat(sprintf("Zero cells: %d of %d (%.1f%%)\n", sum(raw$O == 0, na.rm = TRUE), nrow(raw),
              100 * mean(raw$O == 0, na.rm = TRUE)))
  
  # 3a. months where (almost) every zone is zero: suggests no reporting that month
  mon <- raw[, .(zones = .N, share_zero = mean(O == 0, na.rm = TRUE), cases = sum(O, na.rm = TRUE)), by = .(year, month)][order(year, month)]
  cat("Months with ALL zones zero (national silence):", sum(mon$cases == 0), "\n")
  print(mon[cases == 0])
  save_chk(mon, "national_monthly_summary")
  
  # 3b. zones that never report a case, and zones with cases in only 1-2 months
  zn <- raw[, .(months = .N, positive_months = sum(O > 0, na.rm = TRUE), total_cases = sum(O, na.rm = TRUE),
                mean_pop = mean(monthly_population, na.rm = TRUE)), by = zone][order(positive_months)]
  cat("Zones with zero cases in all months:", sum(zn$positive_months == 0), "of", nrow(zn), "\n")
  cat("Zones with 1-2 positive months:", sum(zn$positive_months %in% 1:2), "\n")
  cat("Distribution of positive months per zone:\n"); print(summary(zn$positive_months))
  save_chk(zn, "zone_summary")
  
  # 3c. zeros by year: a jump in the share of zeros in particular years = reporting changes
  yr <- raw[, .(share_zero = round(mean(O == 0, na.rm = TRUE), 3), positive_cells = sum(O > 0, na.rm = TRUE),
                cases = sum(O, na.rm = TRUE)), by = year][order(year)]
  cat("\nBy year:\n"); print(yr)
  save_chk(yr, "year_summary")
  
  # 3d. long runs of zeros inside zones that DO have cases (gaps inside an outbreak period)
  max_run0 <- function(x) { r <- rle(x == 0); if (any(r$values)) max(r$lengths[r$values]) else 0L }
  setorder(raw, zone, year, month)
  runs <- raw[, .(longest_zero_run = max_run0(O), positive_months = sum(O > 0)), by = zone]
  cat("\nLongest run of consecutive zero months per zone (summary):\n"); print(summary(runs$longest_zero_run))
  
  # ------------------------------------------------------------------------------
  # 4. OUTLIERS AND PLAUSIBILITY
  # ------------------------------------------------------------------------------
  cat("\n=== 4. OUTLIERS ===\n")
  raw[, rate100k := 1e5 * O / monthly_population]
  cat("Largest 15 zone-months by cases:\n")
  print(raw[order(-O)][1:15, .(zone, year, month, cases = O, population = monthly_population, per100k = round(rate100k, 1))])
  cat("Cases greater than population:", sum(raw$O > raw$monthly_population, na.rm = TRUE), "\n")
  cat(sprintf("Top 1%% of positive zone-months hold %.1f%% of all cases\n",
              100 * { p <- sort(raw[O > 0, O], decreasing = TRUE); sum(p[1:ceiling(0.01 * length(p))]) / sum(p) }))
  cat(sprintf("Share of positive cells with a single case: %.1f%%\n", 100 * mean(raw[O > 0, O] == 1)))
  save_chk(raw[order(-O)][1:50, .(zone, year, month, cases = O, population = monthly_population, per100k = rate100k)], "top50_zone_months")
  
  # ------------------------------------------------------------------------------
  # 5. ADMIN CODES (boundary harmonisation)
  # ------------------------------------------------------------------------------
  cat("\n=== 5. ADMIN CODES ===\n")
  zy <- raw[, .(zones = uniqueN(zone)), by = year][order(year)]
  cat("Number of distinct zones per year (should be constant):\n"); print(zy)
  first_last <- raw[, .(first_year = min(year), last_year = max(year)), by = zone]
  cat("Zones not present in every year:", sum(first_last$first_year > min(raw$year) | first_last$last_year < max(raw$year)), "\n")
  if ("admin2Name_en" %in% names(raw) || "admin2Name" %in% names(raw)) {
    nm <- intersect(c("admin2Name_en", "admin2Name"), names(raw))[1]
    per_code_names <- raw[, .(n_names = uniqueN(get(nm))), by = zone]
    cat("Codes with more than one name:", sum(per_code_names$n_names > 1), "\n")
  }
  
  # ------------------------------------------------------------------------------
  # 6. DOES THE POPULATION CHOICE MATTER FOR E?  (compare 3 denominators)
  # ------------------------------------------------------------------------------
  cat("\n=== 6. E under 3 population choices (zone-level Spearman with observed counts) ===\n")
  d[, pop_year_const := monthly_population]                        # as in the data
  d[, t := (year - min(year)) * 12 + month]
  # linear interpolation between annual June values (smooth monthly series)
  interp <- d[, {
    jj <- .SD[month == 6, .(t, p = monthly_population)]
    if (nrow(jj) >= 2) list(t = t, pop_interp = approx(jj$t, jj$p, xout = t, rule = 2)$y) else list(t = t, pop_interp = monthly_population)
  }, by = zone]
  d <- merge(d, interp, by = c("zone", "t"))
  r_nat <- sum(d$total_cases) / sum(d$monthly_population)
  zz <- d[, .(obs = sum(total_cases),
              E_data   = sum(monthly_population * r_nat),
              E_interp = sum(pop_interp * r_nat)), by = zone]
  cat(sprintf("rho(E_data, obs) = %.3f | rho(E_interp, obs) = %.3f | rho(E_data, E_interp) = %.4f\n",
              cor(zz$E_data, zz$obs, method = "spearman"), cor(zz$E_interp, zz$obs, method = "spearman"),
              cor(zz$E_data, zz$E_interp, method = "spearman")))
  cat("Largest relative difference in a zone-month E (data vs interpolated):",
      round(100 * max(abs(d$pop_interp / d$monthly_population - 1)), 1), "%\n")
  cat("\nDone. Files in", chk_dir, "\n")
  
  
  ################################################################################
  # STEP 1 — MODELLING (data build to summary report)
  ################################################################################
  # ==============================================================================
  # ETHIOPIA CHOLERA — FULL STAGE 1 PIPELINE (corrected, ONE file, final)
  #   PART A: hurdle-INLA with new E (INT_POP: intensity also scaled by log zone population); Baseline + Knorr-Held Types I-IV; selection by
  #           leave-one-out log-CPO; posterior draws; joint model; national-trend
  #           model -> INLA_Rebuilt_From_Zero.rds (+ one file per fitted model in fits/)
  #   PART B: temporal relative risk, epidemic waves, Figure 3
  #   PART C: leave-one-out AUC/calibration (occurrence), PIT (intensity),
  #           prior sensitivity of the Baseline model
  #   PART D: five checks of decoupling (occurrence vs intensity zone effects)
  #   PART E: silent-month sensitivity, persistence, last-month effect on intensity
  #   PART B2: waves in the national expected caseload (occurrence x intensity)
  #   PART F: one Markdown summary report of every table (Stage1_summary_report.md)
  # Run top to bottom. MAIN specification: OCC_POP <- TRUE (population in both parts). OCC_POP <- FALSE gives the old occurrence-without-population run.
  # Models already saved in fits/ are loaded, not refitted.
  # ==============================================================================
  
  # ##############################################################################
  # PART A — HURDLE-INLA MODELS
  # ##############################################################################
  # ==============================================================================
  # ETHIOPIA CHOLERA — HURDLE-INLA PIPELINE, REBUILT FROM ZERO (Stage 1) — CHECKED
  #
  # Lines marked "# FIX" changed after the line-by-line check.
  # Lines marked "# [C1]".."# [C5]" show which of Juan's comments they answer.
  #
  #   [C2] E = population x ONE pooled national rate; no smoothing, no floor.
  #   [C2] Missing zone-months stop the script (no population = 1).
  #   [C2] OCC_POP = TRUE refits occurrence with log(population) (sensitivity).
  #   [C3] Uncertainty from 1,000 posterior draws; 500 saved for Stage 2.
  #   [C5] Intensity predicted for EVERY zone-month; B = pi x mu (no 2nd truncation).
  #   [C5] WAIC / DIC / log-CPO disagreement reported automatically.
  #   [C5] Joint model: beta (scaling) AND correlation (from draws), new E.
  #   [C5] National-trend sensitivity model (waves) refitted with new E.
  #   [C1] No exceedance probabilities; no Hotspot/Amplifier/Gateway tiers.
  # ==============================================================================
  
  rm(list = ls())
  options(error = NULL, warn = 0)   # no automatic debugger: if an error occurs you get the message and the script stops
  options(stringsAsFactors = FALSE)
  library(INLA); library(data.table); library(sf); library(spdep)
  library(Matrix); library(pROC); library(dplyr)
  inla.setOption(num.threads = 4)
  
  data_path   <- "F:/phd/DATA_source/correted_E_inala/Cholera_updated_data.Rdata"
  shp_path    <- "F:/Cholera/shape_22/Adm_2/eth_admbnda_adm2_csa_bofedb_2021.shp"
  output_path <- getOption("cholera.proj", ".")                         # FIX: revision folder, keeps old results apart
  
  NDRAW     <- 1000    # posterior draws for uncertainty               [C3]
  NSAVE     <- 500     # draws saved for Stage 2                       [C3]
  OCC_POP   <- TRUE    # TRUE = occurrence + scaled log(population)    [C2]  MAIN SPECIFICATION: population in BOTH parts (D7: association is hidden without it; coefficient 0.49, CrI 0.24-0.74)
  RUN_START <- Sys.time()   # stale-file guard: Step 10 ignores result files written before this run started
  REFIT_ALL <- FALSE   # FALSE = reuse the fits saved by your full refit (fast); set TRUE only to delete them and refit everything again
  SAVE_FIGS <- FALSE   # FALSE = every figure is only DRAWN in the RStudio Plots pane (nothing written); TRUE = also save PNG/PDF/TIFF files
  INT_POP   <- TRUE    # TRUE = intensity + scaled log(zone population) [C2]: the data show outbreak size is NOT proportional to population (exponent ~0.1)
  # PC prior for the temporal RW2 sd, Pr(sd > U) = 0.01, on the SCALED RW2 (scale.model = TRUE).
  # U = 0.3 was far too tight for intensity: with the old (unscaled) prior the national intensity curve
  # corresponds to a scaled sd of about 3 (1.9-7.2), and WAIC was 35 lower (11,998 vs 12,033) with three waves.
  U_RW_OCC  <- 20      # occurrence: Step 3 check 2 showed the posterior sd of the scaled RW2 is about 8.5 (95% CrI 5.9-12.3), above the earlier bound U = 5 (prior-data conflict); U = 20 is weakly informative (as U = 10 for intensity). New fit folders (..._rw20); nothing cached is reused.
  U_RW_INT  <- 10      # intensity: weakly informative on the scale of the data
  # PC prior for the space-time INTERACTION sd (Types I-IV only; the Baseline has no interaction), Pr(sd > U_INT) = 0.01.
  # MAIN RUN: U_INT = 0.5 (the pre-specified value). The fitted Type IV interaction sd was about 5, so 0.5 is tight and may hold the
  # interaction models back in Table 2. SENSITIVITY: set U_INT <- 5, run again (new fit folders are created, nothing is overwritten),
  # and report both. The Baseline conclusions do not depend on this setting.
  U_INT     <- 0.5
  # INTENSITY LIKELIHOOD. TRUE (MAIN) = zero-truncated negative binomial for the positive months (INLA 'zeroinflatednbinomial0' with the zero
  # probability fixed at ~0): the correct likelihood for counts that are >= 1 by construction; WAIC was 224 lower than the ordinary NB and the
  # zone effects, population exponent and conclusions were the same (Step 7C). FALSE = ordinary NB fitted to positive months (the earlier main model;
  # Step 7C then reports the zero-truncated NB as the sensitivity analysis). With TRUE the intensity fits go to new folders (..._ZT): nothing cached is reused.
  INT_ZT    <- TRUE
  RUN_JOINT <- TRUE    # joint shared-component model                  [C5]
  RUN_NAT   <- TRUE    # national-trend sensitivity model (waves)      [C5]
  RUN_DECOUPLE <- TRUE # PART D: five checks of occurrence-vs-intensity decoupling
  RUN_SILENT   <- TRUE # PART E1: refit with the 17 silent months masked
  RUN_PERSIST  <- TRUE # PART E2: persistence extension (last month in zone / neighbours)
  RUN_LAGINT   <- TRUE # PART E3: does an outbreak last month predict LARGER counts this month?
  
  OUT_RDS <- file.path(output_path,                                  # FIX: sensitivity run gets its own file
                       if (OCC_POP) "INLA_Rebuilt_From_Zero_OCCPOP.rds" else "INLA_Rebuilt_From_Zero.rds")
  SUFFIX  <- paste0(if (OCC_POP) "_OCCPOP" else "", if (U_INT != 0.5) paste0("_int", U_INT) else "")   # the interaction-prior sensitivity run writes its own files
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  
  # ==============================================================================
  # 1. LOAD AND CLEAN                                                  [C4]
  # ==============================================================================
  load(data_path)
  raw <- as.data.table(final_comprehensive_dataset_corrected)
  if ("geometry" %in% names(raw)) raw[, geometry := NULL]
  raw[, admin2_clean := toupper(gsub("[^A-Z0-9]", "", admin2Pcode))]
  
  n_raw <- nrow(raw)
  dat <- raw[!is.na(total_cases) & !is.na(monthly_population) & monthly_population > 0 &
               !is.na(year) & !is.na(month)]
  cat("Rows in raw data:", n_raw, "| dropped (missing cases/population/date):", n_raw - nrow(dat), "\n")
  
  # FIX: population check — decides the per-100,000 scaling (section 10)
  pr <- dat[, total_population / monthly_population]
  cat(sprintf("total_population / monthly_population: min %.2f | median %.2f | max %.2f\n",
              min(pr), median(pr), max(pr)))
  cat("  ~1  -> monthly_population = people living in the zone that month (expected)\n")
  cat("  ~12 -> monthly_population = annual population / 12 -> STOP and tell me\n")
  
  # ==============================================================================
  # 2. SPATIAL STRUCTURE
  # ==============================================================================
  shape <- st_read(shp_path, quiet = TRUE) %>%
    mutate(ADM2_PCODE_clean = toupper(gsub("[^A-Z0-9]", "", ADM2_PCODE))) %>%
    st_make_valid()
  common <- intersect(unique(dat$admin2_clean), shape$ADM2_PCODE_clean)
  shape  <- shape %>% filter(ADM2_PCODE_clean %in% common) %>% arrange(ADM2_PCODE_clean)
  dat    <- dat[admin2_clean %in% common]
  n_area <- nrow(shape)
  
  nb <- poly2nb(shape, queen = TRUE)
  islands <- which(card(nb) == 0)
  if (length(islands)) {
    xy <- st_coordinates(st_centroid(shape))
    for (i in islands) { d <- sqrt((xy[,1]-xy[i,1])^2 + (xy[,2]-xy[i,2])^2); d[i] <- Inf; nb[[i]] <- order(d)[1:3] }
  }
  # FIX: island links were one-way (island -> neighbour only). forceSymmetric() keeps
  # only the UPPER triangle, so a link stored below the diagonal was silently lost.
  # make.sym.nb() adds the reverse links, so the graph is truly symmetric.
  nb <- make.sym.nb(nb)
  Q  <- as(nb2mat(nb, style = "B"), "sparseMatrix")
  stopifnot(isSymmetric(as.matrix(Q)), all(rowSums(Q) > 0))
  cat("Zones:", n_area, "| islands linked:", length(islands),
      "| connected components:", n.comp.nb(nb)$nc, " <- should be 1\n")
  
  # ==============================================================================
  # 3. SPACE-TIME INDICES (complete grid required)                     [C2]
  # ==============================================================================
  dat[, time_index := (year - min(year)) * 12 + month]
  dat[, ID_area := match(admin2_clean, shape$ADM2_PCODE_clean)]
  dat[, ID_time := time_index - min(time_index) + 1]
  n_time <- max(dat$ID_time)
  
  stopifnot(!anyDuplicated(dat[, .(ID_area, ID_time)]))
  if (nrow(dat) != n_area * n_time) {
    miss <- CJ(ID_area = 1:n_area, ID_time = 1:n_time)[!dat, on = .(ID_area, ID_time)]
    print(head(miss, 20))
    stop("Grid incomplete: ", nrow(miss), " zone-months missing. Resolve before modelling (do NOT set population = 1).")
  }
  setorder(dat, ID_time, ID_area)
  N <- nrow(dat)
  cat("Zone-months:", N, "(", n_area, "x", n_time, ")\n")
  
  # ==============================================================================
  # 4. EXPECTED COUNTS — single pooled national rate                   [C2]
  # ==============================================================================
  r_nat <- dat[, sum(total_cases) / sum(monthly_population)]
  dat[, `:=`(O = as.integer(total_cases), E = monthly_population * r_nat)]
  dat[, logE := log(E)]
  dat[, Y_bin := as.integer(O > 0)]
  dat[, O_int := fifelse(O > 0, O, NA_integer_)]
  
  cat(sprintf("Pooled national rate: %.3e per person-month | min E = %.3f | E <= 0: %d\n",
              r_nat, min(dat$E), sum(dat$E <= 0)))
  cat(sprintf("Zero zone-months: %d (%.1f%%) | months with cases: %d | total cases: %d\n",
              sum(dat$O == 0), 100 * mean(dat$O == 0), sum(dat$O > 0), sum(dat$O)))
  z_chk <- dat[, .(obs = sum(O), E = sum(E), pop = mean(monthly_population)), by = ID_area]
  cat("Zone SIR (obs/E) range:", paste(round(range(z_chk[obs > 0, obs / E]), 3), collapse = " - "),
      "  <- should vary widely (not all ~1)\n")
  # FIX: direct evidence for the response letter (old E gave rho = 0.999)
  cat(sprintf("Spearman rho(E, observed) zone level = %.3f  <- should now be well below 0.999\n",
              cor(z_chk$E, z_chk$obs, method = "spearman")))
  
  # ==============================================================================
  # 5. COVARIATES AND RANDOM-EFFECT INDICES
  # ==============================================================================
  dat[, `:=`(sin_month  = sin(2 * pi * month / 12),  cos_month  = cos(2 * pi * month / 12),
             sin2_month = sin(4 * pi * month / 12),  cos2_month = cos(4 * pi * month / 12))]
  SEAS <- c("sin_month", "cos_month", "sin2_month", "cos2_month")
  dat[, (SEAS) := lapply(.SD, function(z) as.numeric(scale(z))), .SDcols = SEAS]
  dat[, log_pop_s := as.numeric(scale(log(monthly_population)))]
  # zone-level (time-constant) log population for the intensity model: a BETWEEN-zone scaling term,
  # so it cannot pick up the small within-zone population growth and interfere with the national trend
  zl <- dat[, .(lp = log(mean(monthly_population))), by = ID_area]
  POP_SD <- sd(zl$lp)
  dat[zl, on = "ID_area", lp := i.lp]
  dat[, log_pop_zone_s := (lp - mean(zl$lp)) / POP_SD]
  dat[, lp := NULL]
  
  dat[, `:=`(ID_area_bym = ID_area, ID_time_iid = ID_time, ID_time_rw = ID_time,
             ID_st = (ID_area - 1) * n_time + ID_time)]   # FIX: delta_it index, month within zone (H_space (x) H_time)
  
  # ==============================================================================
  # 6. PRIORS AND CANDIDATE MODELS
  #    BYM2: Pr(sigma > 1) = 0.01 ; Pr(phi < 0.5) = 0.5  (minor comment: prior is on phi)
  # ==============================================================================
  prior_iid  <- list(theta = list(prior = "pc.prec", param = c(0.5, 0.01)))
  prior_rw2  <- list(theta = list(prior = "pc.prec", param = c(U_RW_OCC, 0.01)))
  prior_rw2_int <- list(theta = list(prior = "pc.prec", param = c(U_RW_INT, 0.01)))   # intensity RW2
  prior_int  <- list(theta = list(prior = "pc.prec", param = c(U_INT, 0.01)))   # space-time interaction (setting U_INT at the top of Step 1)
  # ---- intensity-likelihood helpers (ordinary NB or zero-truncated NB, switch INT_ZT at the top of Step 1) ----
  if (INT_ZT && !("zeroinflatednbinomial0" %in% names(INLA::inla.models()$likelihood)))
    stop("INT_ZT = TRUE needs the INLA family 'zeroinflatednbinomial0' (update INLA), or set INT_ZT <- FALSE.")
  INT_FAM    <- if (INT_ZT) "zeroinflatednbinomial0" else "nbinomial"
  is_int_fam <- function(f) f %in% c("nbinomial", "zeroinflatednbinomial0")
  # control.family for the intensity part: the NB size prior (NULL = INLA default) and, for the truncated NB, the zero probability fixed at ~0
  mk_int_cf  <- function(size_hyper = NULL, zt = INT_ZT) {
    if (!zt) return(if (is.null(size_hyper)) list() else list(hyper = list(theta = size_hyper)))
    list(hyper = c(if (!is.null(size_hyper)) list(theta1 = size_hyper), list(theta2 = list(initial = -15, fixed = TRUE))))
  }
  INT_CF     <- mk_int_cf()   # set again by the prior grid if a non-default size prior wins
  size_of    <- function(f) f$summary.hyperpar[grep("^size", rownames(f$summary.hyperpar))[1], "mean"]
  is_zt_fit  <- function(f) { fam <- tryCatch(f$.args$family[1], error = function(e) NULL); if (is.null(fam)) INT_ZT else identical(fam, "zeroinflatednbinomial0") }
  # NB mean parameter mu of a fitted intensity model (the posterior mean of exp(linear predictor)); idx = rows
  mu_nb      <- function(f, idx = NULL) {
    m <- if (is_zt_fit(f)) { lp <- f$summary.linear.predictor; exp(lp$mean + 0.5 * lp$sd^2) } else f$summary.fitted.values$mean
    if (is.null(idx)) m else m[idx]
  }
  zt_q       <- function(mu, sz) -expm1(-sz * log1p(mu / sz))                  # 1 - P0 of the NB (stable for small mu)
  zt_mean    <- function(mu, sz) mu / zt_q(mu, sz)                            # E(Y | Y > 0)
  zt_var     <- function(mu, sz) { q <- zt_q(mu, sz); (mu * (1 + mu / sz) + mu^2) / q - (mu / q)^2 }   # Var(Y | Y > 0)
  # Pearson residual of a positive count y under the family of the fit f
  pearson_int <- function(f, y, idx) { mu <- mu_nb(f, idx); sz <- size_of(f)
  if (is_zt_fit(f)) (y - zt_mean(mu, sz)) / sqrt(zt_var(mu, sz)) else (y - mu) / sqrt(mu * (1 + mu / sz)) }
  # Randomised in-sample PIT of the zero-truncated NB (used only if INLA gives no leave-one-out PIT for this family); NA at zero months
  pit_zt     <- function(f, seed = 12) { out <- rep(NA_real_, nrow(dat)); pos <- which(!is.na(dat$O_int)); y <- dat$O_int[pos]
  sz <- size_of(f); mu <- mu_nb(f, pos); q <- zt_q(mu, sz)
  Fk <- function(k) ifelse(k < 1, 0, (pnbinom(k, size = sz, mu = mu) - (1 - q)) / q)
  set.seed(seed); out[pos] <- runif(length(y), Fk(y - 1), Fk(y)); out }
  prior_bym2 <- list(theta1 = list(prior = "pc.prec", param = c(1, 0.01)),
                     theta2 = list(prior = "pc",      param = c(0.5, 0.5)))
  prior_bym2_int <- prior_bym2; prior_iid_int <- prior_iid    # intensity gets its own copies, so the grid can select different values per part
  int_prior_swap <- function(k) gsub("prior_iid", "prior_iid_int", gsub("prior_bym2", "prior_bym2_int", gsub("prior_rw2", "prior_rw2_int", k, fixed = TRUE), fixed = TRUE), fixed = TRUE)
  
  # ------------------------------------------------------------------------------
  # FIX: BASELINE + KNORR-HELD SPACE-TIME INTERACTIONS (notation as in the Methods)
  #
  #   Baseline (main effects only):
  #     eta_it = alpha + seasonality_t + (u_i + v_i) + gamma_t + phi_t
  #       u_i     structured spatial (ICAR)      } fitted jointly as BYM2
  #       v_i     unstructured spatial (iid)     } (Riebler et al. 2016)
  #       gamma_t structured temporal (RW2)
  #       phi_t   unstructured temporal (iid)
  #
  #   Interaction models: baseline + delta_it, delta ~ N(0, (tau_delta * H_delta)^-)
  #     Type I  : H_delta = H_v (x) H_phi   = I (x) I = I          rank nT
  #     Type II : H_delta = H_v (x) H_gamma = I (x) RW2            rank n(T-2)
  #     Type III: H_delta = H_u (x) H_phi   = ICAR (x) I           rank T(n-1)
  #     Type IV : H_delta = H_u (x) H_gamma = ICAR (x) RW2         rank (n-1)(T-2)
  #   delta_it is ordered month within zone, index (i-1)T + t, so the Kronecker
  #   product is written space (x) time throughout.
  #   Sum-to-zero constraints on delta span the null space of H_delta, so delta
  #   cannot absorb the main effects (identifiability; Goicoa et al. 2018).
  #
  # The OLD code did not do this: it had no interaction terms, and its "RW2 grouped
  # by zone" used INLA's default exchangeable group (one correlation rho_group shared
  # by all zones), not independent or spatially structured interactions.
  # ------------------------------------------------------------------------------
  n <- n_area; T_ <- n_time
  I_n <- Diagonal(n); I_T <- Diagonal(T_)
  H_u <- Diagonal(n, x = rowSums(Q)) - Q                                    # ICAR (CAR) structure
  H_u <- inla.scale.model(H_u, constr = list(A = matrix(1, 1, n), e = 0))
  D2  <- diff(Diagonal(T_), differences = 2)
  H_gamma <- inla.scale.model(crossprod(D2), constr = list(A = rbind(rep(1, T_), 1:T_), e = c(0, 0)))   # RW2
  one <- function(k) matrix(1, 1, k)
  C_month_sum <- kronecker(one(n), diag(T_))          # each month t: sum over zones = 0         (T rows)
  C_zone_sum  <- kronecker(diag(n), one(T_))          # each zone i : sum over months = 0        (n rows)
  C_zone_lin  <- kronecker(diag(n), t(1:T_))          # each zone i : no linear trend in time    (n rows)
  ec <- function(A) list(A = as.matrix(A), e = rep(0, nrow(A)))
  
  H_II  <- kronecker(I_n, H_gamma);  EC_II  <- ec(rbind(C_zone_sum, C_zone_lin));  rd_II  <- 2 * n          # rank n(T-2)
  H_III <- kronecker(H_u, I_T);      EC_III <- ec(C_month_sum);                     rd_III <- T_             # rank T(n-1)
  H_IV  <- kronecker(H_u, H_gamma);  EC_IV  <- ec(rbind(C_month_sum, C_zone_sum[-n, ], C_zone_lin[-n, ]))
  rd_IV <- 2 * n + T_ - 2                                                                                   # rank (n-1)(T-2)
  stopifnot(nrow(EC_II$A) == rd_II, nrow(EC_III$A) == rd_III, nrow(EC_IV$A) == rd_IV,
            n * T_ - rd_IV == (n - 1) * (T_ - 2))
  cat(sprintf("Interaction ranks: II %d = n(T-2) | III %d = T(n-1) | IV %d = (n-1)(T-2)\n",
              n * T_ - rd_II, n * T_ - rd_III, n * T_ - rd_IV))
  
  base <- paste(SEAS, collapse = " + ")
  baseline <- paste(base,
                    "+ f(ID_area_bym, model='bym2', graph=Q, scale.model=TRUE, constr=TRUE, hyper=prior_bym2)",   # u_i + v_i
                    "+ f(ID_time_rw, model='rw2', scale.model=TRUE, constr=TRUE, hyper=prior_rw2)",              # gamma_t
                    "+ f(ID_time_iid, model='iid', hyper=prior_iid)")                                            # phi_t
  g0 <- function(H, EC, rd) sprintf(
    "+ f(ID_st, model='generic0', Cmatrix=%s, rankdef=%s, constr=FALSE, extraconstr=%s, diagonal=1e-6, hyper=prior_int)",
    H, rd, EC)
  KH <- list(
    Baseline = baseline,                                                     # main effects, no interaction
    Type_I   = paste(baseline, "+ f(ID_st, model='iid', hyper=prior_int)"),
    Type_II  = paste(baseline, g0("H_II",  "EC_II",  "rd_II")),
    Type_III = paste(baseline, g0("H_III", "EC_III", "rd_III")),
    Type_IV  = paste(baseline, g0("H_IV",  "EC_IV",  "rd_IV")))
  
  occ_extra <- if (OCC_POP) "+ log_pop_s" else ""
  int_extra <- if (INT_POP) "+ log_pop_zone_s" else ""   # exponent of population = 1 + coefficient / POP_SD (1 = proportional)
  build_forms <- function() {
    forms_occ <<- lapply(KH, function(k) as.formula(paste("Y_bin ~", k, occ_extra)))
    forms_int <<- lapply(KH, function(k) as.formula(paste("O_int ~ offset(logE) +", int_prior_swap(k), int_extra)))
  }
  build_forms()
  
  ################################################################################
  # 6A. PRIOR GRID SEARCH FOR THE BASELINE MODEL, THEN FIT EVERYTHING WITH THE SELECTED PRIORS
  #   Adapts the grid-search idea to your model. Every Baseline refit uses one combination of:
  #     U_phi   BYM2 mixing parameter: prior median of phi  (PC prior, Pr(phi < U_phi) = 0.5)      0.3 | 0.5 | 0.7
  #     U_bym   BYM2 spatial SD:       Pr(sigma > U_bym) = 0.01                                     0.5 | 1   | 2
  #     U_rw    RW2 time trend SD:     Pr(sigma > U_rw)  = 0.01   (0.5x | 1x | 2x of the main-model bound: U_RW_OCC_EFF for occurrence, U_RW_INT_EFF for intensity)
  #     U_iid   iid month SD:          Pr(sigma > U_iid) = 0.01                                     0.25| 0.5 | 1
  #     size    NB overdispersion PC prior (intensity only): INLA default | loose (lambda 2) | tight (lambda 20)
  #   Occurrence: 3^4 = 81 fits.  Intensity: 3^4 x 3 = 243 fits (fast: positive months only).
  #   Stage 1 = SCREENING with a fast INLA setting (Gaussian approximation). Stage 2 = the best 3 + the
  #   default are refitted with the full setting used in the main analysis.
  #   HOW TO USE: this is a SENSITIVITY analysis, not a way to pick priors by the data. Priors stay as
  #   specified in advance (PC priors) unless the grid shows a clear problem. Read:
  #     - zone-effect Spearman vs the default (close to 1 = the zone results do not depend on the prior)
  #     - top-10 zones shared with the default (8-10 = stable)
  #     - range of mean log-CPO across the grid (differences < 0.001 = the data do not distinguish the priors)
  #   The run resumes: finished combinations are read back from the CSV files.
  ################################################################################
  RUN_GRID  <- TRUE            # FALSE = skip the grid and use the pre-specified PC priors above
  SELECT_RULE <- "agree"      # "agree" = change the default only if mean log-CPO AND WAIC both favour the new set; "logcpo" = log-CPO alone
  GRID_MODE <- "full"          # "full" = all combinations; "oat" = default + one change at a time (about 9-11 fits per part)
  cat("\n================ 6A: prior grid search (Baseline, PC priors) ================\n")
  # The occurrence grid files saved under an earlier specification (no population, or another RW2 bound U_RW_OCC) must not be reused:
  # move them aside once (a flag file records that this was done). The intensity grid is unchanged and is reused.
  flag_grid <- file.path(output_path, paste0("grid_occurrence_reset_for_OCCPOP_rw", U_RW_OCC, ".flag"))   # keyed to the RW2 bound: a new bound means a new occurrence grid
  if (OCC_POP && !file.exists(flag_grid)) {
    for (f in c("Table_S_prior_grid_occurrence.csv", "prior_grid_occurrence_effects.rds")) { fp <- file.path(output_path, f)
    if (file.exists(fp)) file.rename(fp, file.path(output_path, sub("(\\.[A-Za-z]+)$", paste0("_OLD_", format(Sys.time(), "%Y%m%d_%H%M"), "\\1"), f))) }
    writeLines("occurrence grid files moved aside for the population-in-both-parts specification", flag_grid)
  }
  PTAG <- ""; U_RW_OCC_EFF <- U_RW_OCC; U_RW_INT_EFF <- U_RW_INT      # effective RW2 bounds (changed below if the grid selects other values)
  if (RUN_GRID) {
    library(ggplot2)
    if (REFIT_ALL) unlink(file.path(output_path, c("Table_S_prior_grid_occurrence.csv", "Table_S_prior_grid_intensity.csv", "prior_grid_occurrence_effects.rds", "prior_grid_intensity_effects.rds")))
    levels_g <- list(U_phi = c(0.3, 0.5, 0.7), U_bym = c(0.5, 1, 2), rw_mult = c(0.5, 1, 2), U_iid = c(0.25, 0.5, 1))
    size_levels <- c("default", "loose", "tight")
    def <- list(U_phi = 0.5, U_bym = 1, rw_mult = 1, U_iid = 0.5, size = "default")
    make_design <- function(with_size) {
      lv <- levels_g; if (with_size) lv$size <- size_levels
      if (GRID_MODE == "full") return(as.data.table(do.call(expand.grid, c(lv, stringsAsFactors = FALSE))))
      rows <- list(as.data.table(def[names(lv)]))
      for (nm in names(lv)) for (v in setdiff(lv[[nm]], def[[nm]])) { r <- def[names(lv)]; r[[nm]] <- v; rows[[length(rows) + 1]] <- as.data.table(r) }
      rbindlist(rows)
    }
    grid_fit <- function(response, family, g, strategy = "gaussian") {
      pri_bym <- list(theta1 = list(prior = "pc.prec", param = c(g$U_bym, 0.01)),
                      theta2 = list(prior = "pc",      param = c(g$U_phi, 0.5)))
      base_rw <- if (response == "O_int") U_RW_INT else U_RW_OCC
      pri_rw  <- list(theta = list(prior = "pc.prec", param = c(base_rw * g$rw_mult, 0.01)))
      pri_iid <- list(theta = list(prior = "pc.prec", param = c(g$U_iid, 0.01)))
      rhs <- paste(paste(SEAS, collapse = " + "),
                   "+ f(ID_area_bym, model='bym2', graph=Q, scale.model=TRUE, constr=TRUE, hyper=pri_bym)",
                   "+ f(ID_time_rw, model='rw2', scale.model=TRUE, constr=TRUE, hyper=pri_rw)",
                   "+ f(ID_time_iid, model='iid', hyper=pri_iid)")
      lhs <- if (response == "O_int") "O_int ~ offset(logE) +" else "Y_bin ~"
      f <- as.formula(paste(lhs, rhs, if (response == "O_int") int_extra else occ_extra), env = environment())
      cf <- list()
      if (response == "O_int" && !is.null(g$size) && g$size != "default")      # PC prior on the NB overdispersion (shrinks towards Poisson); larger lambda = stronger shrinkage
        cf <- list(hyper = list(theta = list(prior = "pc.mgamma", param = if (g$size == "loose") 2 else 20)))
      inla(f, family = family, data = dat, control.family = cf,
           control.predictor = list(compute = TRUE, link = 1),
           control.compute   = list(dic = FALSE, waic = TRUE, cpo = TRUE),
           control.inla      = if (strategy == "gaussian") list(strategy = "gaussian", int.strategy = "eb") else list(strategy = "simplified.laplace"))
    }
    g_zone <- function(f) f$summary.random$ID_area_bym$mean[1:n_area]
    g_time <- function(f) f$summary.random$ID_time_rw$mean
    g_phi  <- function(f) f$summary.hyperpar[grep("^Phi", rownames(f$summary.hyperpar))[1], "mean"]
    g_cpo  <- function(f) { ok <- is.na(f$cpo$failure) | f$cpo$failure == 0; mean(log(f$cpo$cpo[ok]), na.rm = TRUE) }
    g_size <- function(f) { i <- grep("^size", rownames(f$summary.hyperpar))[1]; if (is.na(i)) NA_real_ else f$summary.hyperpar[i, "mean"] }
    
    run_grid <- function(part, response, family, with_size) {
      des <- make_design(with_size); if (!with_size) des[, size := "default"]
      des[, id := sprintf("%s_%s", part, seq_len(.N))]
      fcsv <- file.path(output_path, sprintf("Table_S_prior_grid_%s.csv", part)); feff <- file.path(output_path, sprintf("prior_grid_%s_effects.rds", part))
      done <- if (file.exists(fcsv)) fread(fcsv) else NULL
      eff  <- if (file.exists(feff)) readRDS(feff) else list()
      res  <- if (!is.null(done)) done else NULL
      t_all <- Sys.time()
      for (k in seq_len(nrow(des))) {
        g <- des[k]
        combo <- paste(g$U_phi, g$U_bym, g$rw_mult, g$U_iid, g$size)
        if (!is.null(res) && combo %in% res$combo) next
        t0 <- Sys.time()
        f <- tryCatch(grid_fit(response, family, g), error = function(e) { message("  fit failed: ", conditionMessage(e)); NULL })
        if (is.null(f)) next
        eff[[combo]] <- list(zone = g_zone(f), time = g_time(f))
        row <- data.table(combo = combo, part = part, U_phi = g$U_phi, U_bym = g$U_bym, rw_mult = g$rw_mult, U_iid = g$U_iid, size_prior = g$size,
                          WAIC = f$waic$waic, mean_logCPO = g_cpo(f), CPO_failures = sum(f$cpo$failure > 0, na.rm = TRUE),
                          phi = g_phi(f), size = g_size(f), mode_status = f$mode$mode.status)
        res <- rbind(res, row, fill = TRUE)
        fwrite(res, fcsv); saveRDS(eff, feff)
        el <- as.numeric(difftime(Sys.time(), t0, units = "mins")); left <- nrow(des) - k
        cat(sprintf("  %s %d/%d done in %.1f min (about %.0f min left)\n", part, k, nrow(des), el, el * left))
      }
      # comparison with the default combination (same screening setting)
      kd <- paste(def$U_phi, def$U_bym, def$rw_mult, def$U_iid, "default")
      ref <- eff[[kd]]
      res[, `:=`(zone_spearman_vs_default = vapply(combo, function(k) if (is.null(eff[[k]]) || is.null(ref)) NA_real_ else cor(eff[[k]]$zone, ref$zone, method = "spearman"), numeric(1)),
                 top10_shared = vapply(combo, function(k) if (is.null(eff[[k]]) || is.null(ref)) NA_integer_ else length(intersect(order(-eff[[k]]$zone)[1:10], order(-ref$zone)[1:10])), integer(1)),
                 trend_cor_vs_default = vapply(combo, function(k) if (is.null(eff[[k]]) || is.null(ref)) NA_real_ else cor(eff[[k]]$time, ref$time), numeric(1)))]
      res[, is_default := combo == kd]
      fwrite(res, fcsv)
      res
    }
    cat("\nOccurrence grid...\n"); grid_occ <- run_grid("occurrence", "Y_bin", "binomial", FALSE)
    cat("\nIntensity grid (screening with the ordinary NB; the selected priors are then used for the main intensity likelihood)...\n");  grid_int <- run_grid("intensity", "O_int", "nbinomial", TRUE)
    
    # ---------------- summary ----------------
    for (G in list(grid_occ, grid_int)) {
      part <- G$part[1]; d0 <- G[is_default == TRUE]
      cat(sprintf("\n=== %s: %d fits | default mean log-CPO %.4f | best %.4f | range %.4f ===\n", toupper(part), nrow(G),
                  d0$mean_logCPO, max(G$mean_logCPO), diff(range(G$mean_logCPO))))
      cat(sprintf("Zone-effect Spearman vs default: min %.3f | median %.3f | share >= 0.95: %.0f%%\n",
                  min(G$zone_spearman_vs_default, na.rm = TRUE), median(G$zone_spearman_vs_default, na.rm = TRUE),
                  100 * mean(G$zone_spearman_vs_default >= 0.95, na.rm = TRUE)))
      cat(sprintf("Top-10 zones shared with default: min %d | median %.0f | share >= 8: %.0f%%\n",
                  min(G$top10_shared, na.rm = TRUE), median(G$top10_shared, na.rm = TRUE), 100 * mean(G$top10_shared >= 8, na.rm = TRUE)))
      cat(sprintf("Posterior phi: min %.2f | median %.2f | max %.2f\n", min(G$phi), median(G$phi), max(G$phi)))
      fac <- c("U_phi", "U_bym", "rw_mult", "U_iid", if (part == "intensity") "size_prior")
      me <- rbindlist(lapply(fac, function(v) G[, .(factor = v, level = as.character(get(v)[1]), n = .N,
                                                    mean_logCPO = round(mean(mean_logCPO), 4), mean_phi = round(mean(phi), 3),
                                                    mean_spearman = round(mean(zone_spearman_vs_default, na.rm = TRUE), 3)), by = .(lv = get(v))][, lv := NULL]))
      cat("Main effect of each prior setting (average over the other settings):\n"); print(me)
      fwrite(me, file.path(output_path, sprintf("Table_S_prior_grid_main_effects_%s.csv", part)))
      # ---- rank all combinations: best mean log-CPO first (ties within 0.001 are NOT distinguishable) ----
      best <- max(G$mean_logCPO); tie <- G[mean_logCPO >= best - 0.001]
      topk <- G[order(-mean_logCPO)][1:min(10, .N), .(U_phi, U_bym, rw_mult, U_iid, size_prior, mean_logCPO = round(mean_logCPO, 4),
                                                      diff_vs_best = round(mean_logCPO - best, 4), WAIC = round(WAIC, 1), phi = round(phi, 3),
                                                      zone_spearman_vs_default = round(zone_spearman_vs_default, 3), top10_shared, is_default)]
      cat(sprintf("\nTop 10 of %d combinations by mean log-CPO (%s):\n", nrow(G), part)); print(topk)
      cat(sprintf("Combinations within 0.001 of the best (a tie): %d of %d | default inside the tie: %s\n", nrow(tie), nrow(G),
                  if (any(tie$is_default)) "YES -> keep the default (pre-specified, simplest)" else "NO -> report the best as a sensitivity result"))
      fwrite(G[order(-mean_logCPO)], file.path(output_path, sprintf("Table_S_prior_grid_ranked_%s.csv", part)))
      print(ggplot(G, aes(mean_logCPO, zone_spearman_vs_default, colour = is_default)) + geom_point(alpha = 0.7) +
              scale_colour_manual(values = c("FALSE" = "grey40", "TRUE" = "#B2182B"), guide = "none") +
              labs(x = "Mean log-CPO (higher = better prediction)", y = "Zone-effect Spearman vs default",
                   title = sprintf("Prior grid, %s part (%d fits; red = default)", part, nrow(G))) + theme_minimal())
    }
    
    # ---------------- stage 2: best 3 + default refitted with the full setting ----------------
    cat("\nStage 2: refit the default and the 3 best (by mean log-CPO) of each part with the main-analysis setting...\n")
    conf <- list()
    for (G in list(grid_occ, grid_int)) {
      part <- G$part[1]; resp <- if (part == "occurrence") "Y_bin" else "O_int"; fam <- if (part == "occurrence") "binomial" else "nbinomial"
      pick <- unique(rbind(G[is_default == TRUE], G[order(-mean_logCPO)][1:3]))
      ref_eff <- NULL
      out <- rbindlist(lapply(seq_len(nrow(pick)), function(i) {
        g <- pick[i]; f <- grid_fit(resp, fam, list(U_phi = g$U_phi, U_bym = g$U_bym, rw_mult = g$rw_mult, U_iid = g$U_iid, size = g$size_prior), strategy = "full")
        data.table(part = part, is_default = g$is_default, U_phi = g$U_phi, U_bym = g$U_bym, rw_mult = g$rw_mult, U_iid = g$U_iid,
                   size_prior = g$size_prior, WAIC = round(f$waic$waic, 1), mean_logCPO = round(g_cpo(f), 4), phi = round(g_phi(f), 3),
                   top10 = paste(order(-g_zone(f))[1:10], collapse = ","), zone = I(list(g_zone(f))))
      }))
      rf <- out$zone[[which(out$is_default)[1]]]
      out[, `:=`(zone_spearman_vs_default = vapply(zone, function(z) round(cor(z, rf, method = "spearman"), 3), numeric(1)))]
      out[, zone := NULL]; conf[[part]] <- out
      cat(sprintf("\n--- %s: default vs best 3 (full setting) ---\n", part)); print(out[, !"top10"])
    }
    fwrite(rbindlist(conf, fill = TRUE), file.path(output_path, "Table_S_prior_grid_stage2_confirmation.csv"))
    cat("\nREAD: if the zone-effect Spearman stays about 0.95 or higher, 8-10 of the top-10 zones are shared, and the\n",
        "mean log-CPO differs by less than about 0.001, the results do not depend on the prior. Report the range in the\n",
        "Supplement and keep the pre-specified PC priors. If one setting is clearly better, report it as a sensitivity\n",
        "result and say that the choice used the same data (it is not an independent test).\n")
    
    # ---------------- Methods + Results text, numbers filled in from this run ----------------
    rng <- function(x, d = 2) sprintf("%.*f to %.*f", d, min(x, na.rm = TRUE), d, max(x, na.rm = TRUE))
    txt_grid <- paste0(
      "**Prior sensitivity (hyperparameter grid search).** To assess whether the results depend on the prior specification, we refitted the Baseline models over a grid of penalized-complexity (PC) prior settings. For the occurrence model, the grid crossed (i) the prior median of the BYM2 mixing parameter phi (0.3, 0.5, 0.7; Pr(phi < U) = 0.5), (ii) the upper bound of the spatial standard deviation (U = 0.5, 1, 2; Pr(sigma > U) = 0.01), (iii) the upper bound of the RW2 temporal trend standard deviation (0.5x, 1x and 2x the main-analysis value), and (iv) the upper bound of the iid month standard deviation (U = 0.25, 0.5, 1), giving 3^4 = 81 combinations. For the intensity model, the same grid was crossed with three PC priors on the negative binomial overdispersion (INLA default, a weaker one with lambda = 2 and a stronger one with lambda = 20), giving 243 combinations. The intensity grid was fitted with the ordinary negative binomial; the main intensity model is the zero-truncated negative binomial, which uses the selected priors. All combinations were first fitted with a fast Gaussian approximation; the default specification and the three best-scoring combinations of each part were then refitted with the simplified Laplace approximation used in the main analysis. Models were compared by the mean leave-one-out log conditional predictive ordinate (log-CPO) and WAIC. Robustness was summarised as the Spearman correlation of the zone effects with those of the default specification, the number of the ten highest-ranked zones shared with the default, and the range of phi. Because the same data were used for fitting and comparison, the grid is reported as a sensitivity analysis. The pre-specified default was replaced only if a combination improved the mean log-CPO by more than 0.001 and also had a lower WAIC in the full-setting refit; otherwise the default PC priors were retained for the final models.\n\n",
      sprintf("**Results.** Across the %d occurrence combinations, zone effects correlated with the default at rho = %s (median %.2f), with %d to %d of the top-10 zones shared with the default; mean log-CPO ranged over %.4f (default %.4f, best %.4f) and phi over %s. For intensity (%d combinations), rho = %s (median %.2f), top-10 overlap %d to %d, mean log-CPO range %.4f (default %.4f, best %.4f), and phi %s.\n",
              nrow(grid_occ), rng(grid_occ$zone_spearman_vs_default), median(grid_occ$zone_spearman_vs_default, na.rm = TRUE),
              min(grid_occ$top10_shared, na.rm = TRUE), max(grid_occ$top10_shared, na.rm = TRUE),
              diff(range(grid_occ$mean_logCPO)), grid_occ[is_default == TRUE, mean_logCPO], max(grid_occ$mean_logCPO), rng(grid_occ$phi),
              nrow(grid_int), rng(grid_int$zone_spearman_vs_default), median(grid_int$zone_spearman_vs_default, na.rm = TRUE),
              min(grid_int$top10_shared, na.rm = TRUE), max(grid_int$top10_shared, na.rm = TRUE),
              diff(range(grid_int$mean_logCPO)), grid_int[is_default == TRUE, mean_logCPO], max(grid_int$mean_logCPO), rng(grid_int$phi)),
      "\n*(Numbers come from the screening grid; the stage-2 confirmation table gives the full-setting refits. Check the wording against the printed tables before use.)*\n")
    writeLines(txt_grid, file.path(output_path, "Methods_and_Results_prior_grid.md"))
    cat("\n---------------- Text for the manuscript (saved as Methods_and_Results_prior_grid.md) ----------------\n", txt_grid)
    
    # ---------------- SELECTION: best of {default + 3 best}, tie rule, then USE it for all models below ----------------
    # SELECT_RULE "agree" (default): leave the pre-specified default only if a combination is better by more than 0.001 in mean
    # log-CPO AND has a WAIC no higher than the default (both criteria agree). "logcpo": mean log-CPO alone.
    pick_best <- function(out) {
      d0 <- out[is_default == TRUE][1]
      out[, n_changed := (U_phi != def$U_phi) + (U_bym != def$U_bym) + (rw_mult != def$rw_mult) + (U_iid != def$U_iid) + (size_prior != def$size)]
      cand <- out[is_default == FALSE & mean_logCPO > d0$mean_logCPO + 0.001 & (SELECT_RULE == "logcpo" | WAIC <= d0$WAIC)]
      if (!nrow(cand)) return(out[is_default == TRUE][1])                       # no clear improvement: keep the default
      tie <- cand[mean_logCPO >= max(cand$mean_logCPO) - 0.001]                  # ties among the improvements: fewest changes from the default
      tie[order(n_changed, -mean_logCPO)][1]
    }
    sel_occ <- pick_best(conf$occurrence); sel_int <- pick_best(conf$intensity)
    selp <- rbind(sel_occ[, .(part, U_phi, U_bym, rw_mult, U_iid, size_prior, mean_logCPO, phi, selected_is_default = n_changed == 0)],
                  sel_int[, .(part, U_phi, U_bym, rw_mult, U_iid, size_prior, mean_logCPO, phi, selected_is_default = n_changed == 0)])
    cat("\n=== SELECTED PRIORS (used for all models below) ===\n"); print(selp)
    fwrite(selp, file.path(output_path, "Selected_priors.csv"))
    prior_bym2 <- list(theta1 = list(prior = "pc.prec", param = c(sel_occ$U_bym, 0.01)), theta2 = list(prior = "pc", param = c(sel_occ$U_phi, 0.5)))
    prior_iid  <- list(theta = list(prior = "pc.prec", param = c(sel_occ$U_iid, 0.01)))
    prior_rw2  <- list(theta = list(prior = "pc.prec", param = c(U_RW_OCC * sel_occ$rw_mult, 0.01)))
    prior_bym2_int <- list(theta1 = list(prior = "pc.prec", param = c(sel_int$U_bym, 0.01)), theta2 = list(prior = "pc", param = c(sel_int$U_phi, 0.5)))
    prior_iid_int  <- list(theta = list(prior = "pc.prec", param = c(sel_int$U_iid, 0.01)))
    prior_rw2_int  <- list(theta = list(prior = "pc.prec", param = c(U_RW_INT * sel_int$rw_mult, 0.01)))
    INT_CF <- mk_int_cf(if (sel_int$size_prior == "default") NULL else list(prior = "pc.mgamma", param = if (sel_int$size_prior == "loose") 2 else 20))
    build_forms()
    U_RW_OCC_EFF <- U_RW_OCC * sel_occ$rw_mult; U_RW_INT_EFF <- U_RW_INT * sel_int$rw_mult
    PTAG <- if (sel_occ$n_changed == 0 && sel_int$n_changed == 0) "" else
      paste0("_selP_o", paste(sel_occ$U_phi, sel_occ$U_bym, sel_occ$rw_mult, sel_occ$U_iid, sep = "-"),
             "_i", paste(sel_int$U_phi, sel_int$U_bym, sel_int$rw_mult, sel_int$U_iid, sel_int$size_prior, sep = "-"))
    cat("Fit folders get the tag:", if (nzchar(PTAG)) PTAG else "(none: selected = default priors)", "\n")
  } else cat("RUN_GRID is FALSE: prior grid skipped.\n")
  
  
  
  # ==============================================================================
  # 7. FIT ALL CANDIDATES
  # ==============================================================================
  # FIX (speed): inla.cpo() refits the model once per CPO failure. With 11,040
  # interaction terms each refit is slow (one model ran overnight), so it is OFF.
  # CPO failures are counted and reported instead, and mean log-CPO is also given
  # with failed observations excluded (see section 8).
  RECOMPUTE_CPO <- FALSE
  
  safe_inla <- function(f, family) {
    for (s in c("simplified.laplace", "adaptive")) {
      fit <- tryCatch(inla(f, family = family, data = dat, control.family = if (is_int_fam(family)) INT_CF else list(),
                           control.predictor = list(compute = TRUE, link = 1),
                           control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                           control.inla      = list(strategy = s)),
                      error = function(e) { message("  ", s, " failed: ", e$message); NULL })
      if (!is.null(fit)) {
        if (RECOMPUTE_CPO && sum(fit$cpo$failure > 0, na.rm = TRUE) > 0)
          fit <- tryCatch(inla.cpo(fit), error = function(e) { message("  inla.cpo failed"); fit })
        return(fit)
      }
    }
    stop("INLA failed: ", deparse(f))
  }
  
  # FIX (resume): each model is saved to its own file as soon as it finishes.
  # If R stops or crashes, re-running skips every model already saved.
  # occurrence and intensity fits live in separate folders, so changing the intensity specification
  # does not force the (slow) occurrence models to be refitted
  UTAG <- if (U_INT != 0.5) paste0("_int", U_INT) else ""      # a different interaction prior gets its own fit folders (no cached fit is reused)
  fit_dir_occ <- file.path(output_path, paste0("fits", SUFFIX, if (U_RW_OCC != 0.3) paste0("_rw", U_RW_OCC) else "", PTAG))     # SUFFIX already carries the U_INT tag
  ZTAG <- if (INT_ZT) "_ZT" else ""                            # zero-truncated intensity fits live in their own folder
  fit_dir_int <- file.path(output_path, paste0("fits", if (INT_POP) "_INTPOP" else "", "_rw", U_RW_INT, PTAG, UTAG, ZTAG))
  fit_dir <- fit_dir_occ                                   # kept for older references
  if (REFIT_ALL) { unlink(c(fit_dir_occ, fit_dir_int), recursive = TRUE); cat("REFIT_ALL = TRUE: cached fits deleted, every model is refitted.\n") }
  dir.create(fit_dir_occ, showWarnings = FALSE); dir.create(fit_dir_int, showWarnings = FALSE)
  fit_or_load <- function(forms, family, part, dir) {
    out <- list()
    for (nm in names(forms)) {
      fp <- file.path(dir, paste0(part, "_", nm, ".rds"))
      if (file.exists(fp)) {
        x_old <- suppressWarnings(readRDS(fp))                      # saved fits can mention a package namespace that is not loaded: harmless
        fam_ok <- tryCatch({ fm <- as.character(x_old$.args$family); length(fm) == 0 || identical(fm[1], family[1]) }, error = function(e) TRUE)
        n_ok   <- tryCatch(nrow(x_old$summary.fitted.values) == nrow(dat), error = function(e) FALSE)
        if (isTRUE(fam_ok) && isTRUE(n_ok)) { cat(sprintf("  %s %-8s loaded from file (%s)\n", part, nm, format(file.info(fp)$mtime, "%Y-%m-%d %H:%M"))); out[[nm]] <- x_old; next }
        cat(sprintf("  %s %-8s saved fit does not match the current data/likelihood -> refitted\n", part, nm))
      }
      t0 <- Sys.time()
      out[[nm]] <- safe_inla(forms[[nm]], family)
      saveRDS(out[[nm]], fp)
      cat(sprintf("  %s %-8s done in %.1f min | CPO failures: %d\n", part, nm,
                  as.numeric(difftime(Sys.time(), t0, units = "mins")),
                  sum(out[[nm]]$cpo$failure > 0, na.rm = TRUE)))
    }
    out
  }
  cat("\nFitting occurrence models...\n"); fits_occ <- fit_or_load(forms_occ, "binomial",  "occ", fit_dir_occ)
  cat("Fitting intensity models (", if (INT_ZT) "zero-truncated NB" else "ordinary NB", ")...\n", sep = "");   fits_int <- fit_or_load(forms_int, INT_FAM, "int", fit_dir_int)
  
  # ==============================================================================
  # 8. MODEL SELECTION (best leave-one-out log-CPO; WAIC/DIC reported) [C5]
  # ==============================================================================
  # Logarithmic score: LS = -sum(log CPO_i); smaller = better prediction.
  # mean log-CPO = -LS / n; larger (closer to 0) = better. Same information.
  metrics <- function(fits, auc = FALSE) {
    cpo_ok <- function(x) is.na(x$cpo$failure) | x$cpo$failure == 0
    # FIX: observations whose CPO is reliable in EVERY candidate model. Models are
    # compared on this common set so that each is scored on the same observations
    # (a model with more failures is otherwise scored on an easier subset).
    ok_all <- Reduce(`&`, lapply(fits, cpo_ok))
    m <- data.table(Model = names(fits),
                    WAIC = sapply(fits, function(x) x$waic$waic),
                    DIC  = sapply(fits, function(x) x$dic$dic),
                    mean_logCPO        = sapply(fits, function(x) mean(log(x$cpo$cpo), na.rm = TRUE)),
                    CPO_failures       = sapply(fits, function(x) sum(x$cpo$failure > 0, na.rm = TRUE)),
                    mean_logCPO_noFail = sapply(fits, function(x) mean(log(x$cpo$cpo[cpo_ok(x)]), na.rm = TRUE)),
                    mean_logCPO_common = sapply(fits, function(x) mean(log(x$cpo$cpo[ok_all]), na.rm = TRUE)),
                    LS_common          = sapply(fits, function(x) -sum(log(x$cpo$cpo[ok_all]), na.rm = TRUE)),
                    n_common           = sapply(fits, function(x) sum(ok_all & !is.na(x$cpo$cpo))))
    # FIX: reliability of WAIC (Vehtari et al. 2017) and leave-one-out deviance
    m[, `:=`(p_WAIC = sapply(fits, function(x) x$waic$p.eff),
             p_DIC  = sapply(fits, function(x) x$dic$p.eff),
             n_obs_pWAIC_gt_0.4 = sapply(fits, function(x) sum(x$waic$local.p.eff > 0.4, na.rm = TRUE)),
             LOO_deviance = sapply(fits, function(x) -2 * sum(log(x$cpo$cpo), na.rm = TRUE)))]
    if (auc) m[, AUC_insample := sapply(fits, function(x)
      as.numeric(roc(dat$Y_bin, x$summary.fitted.values$mean[1:N], quiet = TRUE)$auc))]
    # FIX: SELECTION RULE = out-of-sample prediction: mean log-CPO on the common
    # set of reliable observations. WAIC/DIC favoured the interaction models but
    # were unreliable for them (hundreds of observations with p_WAIC > 0.4) and
    # their leave-one-out prediction was worse. Models within 0.001 of the best
    # score are treated as tied and the simplest (earliest in the list) is chosen.
    m[, order_in_list := seq_len(.N)]
    best <- max(m$mean_logCPO_common)
    m[, selected := FALSE]
    m[mean_logCPO_common >= best - 0.001][order(order_in_list)][1, Model] -> sel
    m[Model == sel, selected := TRUE]
    m[order(-mean_logCPO_common)][, order_in_list := NULL][]
  }
  met_occ <- metrics(fits_occ, auc = TRUE)
  met_int <- metrics(fits_int)
  cat("\n--- Occurrence (Table 2) ---\n"); print(met_occ, digits = 5)
  cat("\n--- Intensity (Table 2) ---\n");  print(met_int, digits = 6)
  fwrite(met_occ, file.path(output_path, paste0("Table2_occurrence_metrics", SUFFIX, ".csv")))   # FIX: save tables
  fwrite(met_int, file.path(output_path, paste0("Table2_intensity_metrics",  SUFFIX, ".csv")))
  
  sel_o <- met_occ[selected == TRUE, Model]; sel_i <- met_int[selected == TRUE, Model]
  fit_o <- fits_occ[[sel_o]]
  fit_i <- fits_int[[sel_i]]
  cat("\nSelected (best leave-one-out log-CPO): occurrence =", sel_o, "| intensity =", sel_i, "\n")
  for (nm in c("occ", "int")) {
    m <- if (nm == "occ") met_occ else met_int
    w <- m[which.min(WAIC), Model]
    if (w != m[selected == TRUE, Model])
      cat("NOTE (", nm, "): best WAIC is", w, "-> report the WAIC/DIC vs log-CPO disagreement in the text.\n")
  }
  
  # ==============================================================================
  # 9. REPORTABLE PARAMETERS (Table 3)                                 [C1][C5]
  # ==============================================================================
  tab <- function(f) round(f[, c("mean", "0.025quant", "0.975quant")], 3)
  cat("\n--- Occurrence fixed effects ---\n");  print(tab(fit_o$summary.fixed))
  cat("\n--- Occurrence hyperparameters ---\n"); print(tab(fit_o$summary.hyperpar))
  cat("\n--- Intensity fixed effects ---\n");   print(tab(fit_i$summary.fixed))
  cat("\n--- Intensity hyperparameters ---\n");  print(tab(fit_i$summary.hyperpar))
  cat("Phi = structured share of the spatial random-effect VARIANCE (not 'spillover').\n")
  cat("Precision for ID_st = space-time interaction (Knorr-Held type of the selected model).\n")
  
  # FIX: bimodal-season check (minor comment), printed explicitly
  for (nm in c("occurrence", "intensity")) {
    f <- if (nm == "occurrence") fit_o$summary.fixed else fit_i$summary.fixed
    h2 <- f[c("sin2_month", "cos2_month"), c("0.025quant", "0.975quant")]
    cat(sprintf("Second harmonic (%s): %s\n", nm,
                if (any(h2[,1] > 0 | h2[,2] < 0)) "at least one CrI excludes 0 -> model supports a second peak"
                else "both CrIs include 0 -> do NOT claim a bimodal season from this model"))
  }
  
  if (INT_POP && "log_pop_zone_s" %in% rownames(fit_i$summary.fixed)) {
    bp <- fit_i$summary.fixed["log_pop_zone_s", ]
    cat(sprintf("Population exponent of outbreak size (1 = proportional): %.2f (95%% CrI %.2f to %.2f)\n",
                1 + bp$mean / POP_SD, 1 + bp[["0.025quant"]] / POP_SD, 1 + bp[["0.975quant"]] / POP_SD))
  }
  size_row <- grep("^size", rownames(fit_i$summary.hyperpar))
  if (length(size_row)) {
    k <- fit_i$summary.hyperpar[size_row, ]
    cat(sprintf("Overdispersion alpha = 1/size = %.2f (95%% CrI %.2f-%.2f)\n",
                1 / k$mean, 1 / k[["0.975quant"]], 1 / k[["0.025quant"]]))
  }
  tab3 <- rbind(                                                       # FIX: save Table 3
    data.table(component = "occurrence", type = "fixed", term = rownames(fit_o$summary.fixed),    tab(fit_o$summary.fixed)),
    data.table(component = "occurrence", type = "hyper", term = rownames(fit_o$summary.hyperpar), tab(fit_o$summary.hyperpar)),
    data.table(component = "intensity",  type = "fixed", term = rownames(fit_i$summary.fixed),    tab(fit_i$summary.fixed)),
    data.table(component = "intensity",  type = "hyper", term = rownames(fit_i$summary.hyperpar), tab(fit_i$summary.hyperpar)))
  fwrite(tab3, file.path(output_path, paste0("Table3_parameters", SUFFIX, ".csv")))
  
  # ---- SEASONAL CURVE BY CALENDAR MONTH: which months are highest? (posterior draws, with uncertainty) ----
  # The Fourier coefficients cannot be read month by month. This rebuilds the seasonal curve for months 1-12 from the
  # posterior draws of the four coefficients and shows it relative to the average month. Only the SEASONAL shape is shown
  # (the national wave trend and the zone effects are not included). Occurrence: odds ratio of an outbreak vs the average
  # month. Intensity: ratio of the expected outbreak size vs the average month.
  season_curve <- function(fit, part, ratio_name) {
    sel <- setNames(as.list(rep(1, length(SEAS))), SEAS)
    sm  <- inla.posterior.sample(1000, fit, selection = sel)
    rn  <- sub(":1$", "", rownames(sm[[1]]$latent))
    Bd  <- t(sapply(sm, function(z) z$latent[, 1]))[, match(SEAS, rn), drop = FALSE]
    bm  <- unique(dat[, c("month", SEAS), with = FALSE])[order(month)]
    E   <- Bd %*% t(as.matrix(bm[, ..SEAS]))                      # draws x 12, model scale
    Ec  <- E - rowMeans(E)                                        # relative to the average month
    pk  <- tabulate(apply(E, 1, which.max), nbins = 12) / nrow(E)
    data.table(part = part, month = bm$month, month_name = month.abb[bm$month],
               ratio = exp(apply(Ec, 2, median)), lo = exp(apply(Ec, 2, quantile, 0.025)), hi = exp(apply(Ec, 2, quantile, 0.975)),
               P_highest_month = round(pk, 3), measure = ratio_name)
  }
  sc <- rbind(season_curve(fit_o, "Occurrence (outbreak yes/no)", "odds ratio vs average month"),
              season_curve(fit_i, "Intensity (outbreak size)",     "size ratio vs average month"))
  sc[, `:=`(ratio = round(ratio, 2), lo = round(lo, 2), hi = round(hi, 2))]
  cat("\n--- Seasonal curve by month (relative to the average month; 95% CrI) ---\n"); print(sc[, .(part, month_name, ratio, lo, hi, P_highest_month)])
  for (pp in unique(sc$part)) { d <- sc[part == pp][order(-ratio)]
  cat(sprintf("%s: highest %s (%.2f), %s (%.2f), %s (%.2f); lowest %s (%.2f)\n", pp, d$month_name[1], d$ratio[1], d$month_name[2], d$ratio[2],
              d$month_name[3], d$ratio[3], d$month_name[12], d$ratio[12])) }
  # ---- RISK WINDOW AND PEAK (seasonal): consecutive months above average, wrapping across December-January ----
  runs_text <- function(d, flag) {                       # d ordered by month; flag = logical, length 12
    if (!any(flag)) return("none")
    if (all(flag)) return("all months")
    st <- which(!flag)[1]; idx <- c(st:12, if (st > 1) 1:(st - 1)); f <- flag[idx]; m <- d$month_name[idx]
    r <- rle(f); en <- cumsum(r$lengths); sta <- en - r$lengths + 1
    paste(vapply(which(r$values), function(k) if (sta[k] == en[k]) m[sta[k]] else paste0(m[sta[k]], "-", m[en[k]]), character(1)), collapse = ", ")
  }
  win <- rbindlist(lapply(unique(sc$part), function(pp) { d <- sc[part == pp][order(month)]; pk <- d[which.max(ratio)]
  data.table(part = pp, peak_month = pk$month_name, peak_ratio = pk$ratio, peak_ratio_lo = pk$lo, peak_ratio_hi = pk$hi, P_peak_month = pk$P_highest_month,
             window_above_average = runs_text(d, d$ratio > 1), window_credibly_above_average = runs_text(d, d$lo > 1),
             months_in_credible_window = sum(d$lo > 1)) }))
  cat("\n--- Seasonal RISK WINDOW and PEAK ---\n"); print(win)
  fwrite(win, file.path(output_path, paste0("Table_S_seasonal_risk_window", SUFFIX, ".csv")))
  cat("window_above_average = months whose central estimate is above the average month; window_credibly_above_average =\n",
      "months whose whole 95% CrI is above 1 (the defensible 'risk window'). Word it as the 'period of highest outbreak\n",
      "probability / expected size', not as a transmission window. It is the SEASONAL shape only; epidemic waves differ by year.\n")
  fwrite(sc, file.path(output_path, paste0("Table_S_seasonal_curve_by_month", SUFFIX, ".csv")))
  print(ggplot(sc, aes(factor(month_name, levels = month.abb), ratio, group = part)) + geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
          geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey80", alpha = 0.6) + geom_line() + geom_point() + facet_wrap(~part) + scale_y_log10() +
          labs(x = NULL, y = "Ratio vs average month (log scale)", title = "Seasonal pattern by calendar month (model-based, 95% CrI)") + theme_minimal())
  cat("READ: a month is credibly above average when its whole interval is above 1. P_highest_month is the share of posterior\n",
      "draws in which that month is the seasonal peak. Wide intervals or a flat curve = weak seasonality; the wave trend,\n",
      "not the season, drives most of the variation.\n")
  
  # ==============================================================================
  # 10. EXPECTED CASELOAD FROM POSTERIOR DRAWS                         [C3][C5]
  #     pi = P(outbreak), mu = E(cases | outbreak), B = pi * mu (no extra truncation)
  #     The two hurdle parts share no parameters, so their posteriors are
  #     independent and combining independent draws gives the correct joint draw.
  # ==============================================================================
  cat("\nDrawing", NDRAW, "posterior samples from each model...\n")
  draw_eta <- function(fit) {
    s <- inla.posterior.sample(NDRAW, fit, selection = list(Predictor = 1:N))
    sapply(s, function(x) x$latent[, 1])
  }
  set.seed(2024); eta_o <- draw_eta(fit_o)
  set.seed(2025); smp_i <- inla.posterior.sample(NDRAW, fit_i, selection = list(Predictor = 1:N))
  eta_i <- sapply(smp_i, function(x) x$latent[, 1])
  k_sz  <- grep("^size", names(smp_i[[1]]$hyperpar))[1]                  # size (1/overdispersion) of the intensity part, one value per draw
  sz_draw <- if (is.na(k_sz)) rep(size_of(fit_i), NDRAW) else vapply(smp_i, function(x) unname(x$hyperpar[k_sz]), numeric(1))
  rm(smp_i); invisible(gc())
  PI <- plogis(eta_o); MU_NB <- exp(eta_i)                              # MU_NB = mean parameter of the NB (before truncation)
  cat("Check: draws vs INLA means of the NB mean parameter (should be ~1):",
      round(cor(rowMeans(MU_NB), mu_nb(fit_i, 1:N)), 3), "\n")
  MU <- MU_NB                                                           # MU = expected cases GIVEN an outbreak, E(Y | Y > 0)
  if (INT_ZT) { for (b in seq_len(NDRAW)) MU[, b] <- zt_mean(MU_NB[, b], sz_draw[b]) }   # zero-truncated NB: mean = mu / (1 - P0), per draw
  MU_NB <- MU_NB[, seq_len(min(500, NDRAW)), drop = FALSE]              # kept for the predictive checks (Step 3, check 4)
  B <- PI * MU
  rm(eta_o, eta_i); invisible(gc())                                    # FIX: free memory
  
  qs <- function(M) t(apply(M, 1, quantile, probs = c(0.025, 0.975)))
  q_pi <- qs(PI); q_mu <- qs(MU); q_B <- qs(B)
  res <- dat[, .(admin2_clean, zone, ID_area, ID_time, year, month, O, E, monthly_population, total_population)]
  res[, `:=`(pi = rowMeans(PI), pi_lo = q_pi[,1], pi_hi = q_pi[,2],
             mu = rowMeans(MU), mu_lo = q_mu[,1], mu_hi = q_mu[,2],
             B  = rowMeans(B),  B_lo  = q_B[,1],  B_hi  = q_B[,2],
             predicted_only = is.na(dat$O_int))]
  fwrite(res, file.path(output_path, paste0("zone_month_results", SUFFIX, ".csv")))
  
  # Zone summaries: aggregate WITHIN each draw, then take quantiles (proper CrIs)
  zsum <- function(M, scale = NULL) {
    X <- if (is.null(scale)) M else M * scale
    Z <- rowsum(X, dat$ID_area) / as.vector(table(dat$ID_area))
    data.table(ID_area = as.integer(rownames(Z)), mean = rowMeans(Z),
               lo = apply(Z, 1, quantile, 0.025), hi = apply(Z, 1, quantile, 0.975))
  }
  # FIX: ONE denominator for every rate = monthly_population (the same one used in E).
  # The old code mixed total_population (burden) and monthly_population (intensity).
  per100k_month <- 1 / dat$monthly_population * 1e5        # cases per 100,000 per month
  per100k_year  <- 12 * per100k_month                      # annualised per 100,000
  zone <- data.table(ID_area = 1:n_area, admin2_clean = shape$ADM2_PCODE_clean,
                     zone = dat[match(1:n_area, ID_area), zone])
  zone <- zone[zsum(PI)[, .(ID_area, pi = mean, pi_lo = lo, pi_hi = hi)], on = "ID_area"]
  zone <- zone[zsum(MU)[, .(ID_area, mu = mean, mu_lo = lo, mu_hi = hi)], on = "ID_area"]
  zone <- zone[zsum(MU, per100k_month)[, .(ID_area, mu_100k = mean, mu_100k_lo = lo, mu_100k_hi = hi)], on = "ID_area"]
  zone <- zone[zsum(B,  per100k_year)[,  .(ID_area, burden_100k = mean, burden_lo = lo, burden_hi = hi)], on = "ID_area"]
  zone <- zone[dat[, .(N_outbreak_months = sum(O > 0), total_cases = sum(O)), by = ID_area], on = "ID_area"]
  fwrite(zone, file.path(output_path, paste0("TableS3_zone_summary", SUFFIX, ".csv")))
  
  cat("\n--- Top 10 zones: outbreak probability ---\n");                 print(zone[order(-pi)][1:10, .(zone, pi = round(pi, 3), pi_lo = round(pi_lo, 3), pi_hi = round(pi_hi, 3))])
  cat("\n--- Top 10 zones: cases per 100k per outbreak month ---\n");    print(zone[order(-mu_100k)][1:10, .(zone, mu_100k = round(mu_100k, 1), mu_100k_lo = round(mu_100k_lo, 1), mu_100k_hi = round(mu_100k_hi, 1))])
  cat("\n--- Top 10 zones: annual expected caseload per 100k ---\n");    print(zone[order(-burden_100k)][1:10, .(zone, burden_100k = round(burden_100k, 1), burden_lo = round(burden_lo, 1), burden_hi = round(burden_hi, 1))])
  
  # ==============================================================================
  # 11. JOINT SHARED-COMPONENT MODEL (decoupling test) — new E         [C5]
  #     Occurrence: logit(pi) = seasonal + w_i
  #     Intensity : log(mu)   = log(E) + seasonal + beta * w_i + v_i
  #     w, v = BYM2 fields. beta = scaling factor; correlation(w, beta*w + v)
  #     is derived separately from posterior draws. Temporal terms are omitted
  #     here (simplified relative to the selected models) -> say so in Methods.
  # ==============================================================================
  joint <- NULL
  if (RUN_JOINT) {
    cat("\nFitting joint shared-component model...\n")
    di <- dat[O > 0]; n_o <- N; n_i <- nrow(di); NN <- n_o + n_i
    zo <- function(v) c(as.numeric(v), rep(0, n_i)); zi <- function(v) c(rep(0, n_o), as.numeric(v))
    jd <- list(                                                        # FIX: list, so the 2-column Y is safe
      int_o = c(rep(1, n_o), rep(0, n_i)), int_i = c(rep(0, n_o), rep(1, n_i)),
      s1o = zo(dat$sin_month), c1o = zo(dat$cos_month), s2o = zo(dat$sin2_month), c2o = zo(dat$cos2_month),
      s1i = zi(di$sin_month),  c1i = zi(di$cos_month),  s2i = zi(di$sin2_month),  c2i = zi(di$cos2_month),
      off = c(rep(0, n_o), di$logE),
      lp_i = c(rep(0, n_o), di$log_pop_zone_s),                       # zone population, intensity part
      lp_o = zo(dat$log_pop_s),                                       # zone population, occurrence part (used when OCC_POP)
      s_shared = c(dat$ID_area, rep(NA, n_i)),
      s_copy   = c(rep(NA, n_o), di$ID_area),
      s_own    = c(rep(NA, n_o), di$ID_area))
    Yj <- matrix(NA_real_, NN, 2); Yj[1:n_o, 1] <- dat$Y_bin; Yj[(n_o + 1):NN, 2] <- di$O; jd$Y <- Yj
    hb <- list(beta = list(prior = "normal", param = c(0, 0.1)))     # mean 0, precision 0.1 (variance 10)
    fj <- Y ~ -1 + int_o + s1o + c1o + s2o + c2o + int_i + s1i + c1i + s2i + c2i + offset(off) +
      f(s_shared, model = "bym2", graph = Q, scale.model = TRUE, constr = TRUE, hyper = prior_bym2) +
      f(s_copy, copy = "s_shared", fixed = FALSE, hyper = hb) +
      f(s_own, model = "bym2", graph = Q, scale.model = TRUE, constr = TRUE, hyper = prior_bym2)
    if (INT_POP) fj <- update(fj, . ~ . + lp_i)                       # same population term as the intensity model
    if (OCC_POP) fj <- update(fj, . ~ . + lp_o)                       # and as the occurrence model
    mj <- inla(fj, family = c("binomial", INT_FAM), data = jd, control.family = list(list(), INT_CF),
               Ntrials = c(rep(1, n_o), rep(NA, n_i)),
               control.predictor = list(compute = TRUE, link = c(rep(1, n_o), rep(2, n_i))),
               control.compute = list(waic = TRUE, config = TRUE),
               control.inla = list(strategy = "simplified.laplace"))
    br <- grep("^Beta", rownames(mj$summary.hyperpar))[1]
    beta <- mj$summary.hyperpar[br, ]
    set.seed(42); sm <- inla.posterior.sample(NDRAW, mj)
    rn <- rownames(sm[[1]]$latent)
    iw <- grep("^s_shared:", rn)[1:n_area]; iv <- grep("^s_own:", rn)[1:n_area]   # first n = total BYM2 effect
    bi <- grep("^Beta", names(sm[[1]]$hyperpar))[1]
    rr <- vapply(sm, function(s) { w <- s$latent[iw]; cor(w, s$hyperpar[bi] * w + s$latent[iv]) }, numeric(1))
    qr <- quantile(rr, c(.025, .5, .975))
    cat(sprintf("Joint model: beta = %.3f (95%% CrI %.3f to %.3f) | spatial correlation = %.3f (95%% CrI %.3f to %.3f) | WAIC %.1f\n",
                beta$mean, beta[["0.025quant"]], beta[["0.975quant"]], qr[2], qr[1], qr[3], mj$waic$waic))
    joint <- list(fit = mj, beta = beta, cor_draws = rr, cor_q = qr)
    fwrite(data.table(beta = beta$mean, beta_lo = beta[["0.025quant"]], beta_hi = beta[["0.975quant"]],   # FIX: save
                      cor = qr[2], cor_lo = qr[1], cor_hi = qr[3], WAIC = mj$waic$waic),
           file.path(output_path, paste0("joint_model_beta_correlation", SUFFIX, ".csv")))
    rm(sm); invisible(gc())
  }
  
  # ==============================================================================
  # 12. NATIONAL-TREND SENSITIVITY MODEL (waves, Figure 3)             [C5]  # FIX: was missing
  #     The waves in the paper came from this model, not the area-specific one,
  #     so it must be refitted with the new E as well.
  #     Spec (confirmed against your original fit_nat): seasonal + BYM2 + ONE
  #     national RW2 (not grouped by zone).
  # ==============================================================================
  nat <- NULL
  if (RUN_NAT) {
    cat("\nFitting national-trend sensitivity model...\n")
    KH_nat <- paste(base, "+ f(ID_area_bym, model='bym2', graph=Q, scale.model=TRUE, constr=TRUE, hyper=prior_bym2)",
                    "+ f(ID_time_rw, model='rw2', scale.model=TRUE, constr=TRUE, hyper=prior_rw2)")
    nat_o <- safe_inla(as.formula(paste("Y_bin ~", KH_nat, occ_extra)), "binomial")
    nat_i <- safe_inla(as.formula(paste("O_int ~ offset(logE) +", int_prior_swap(KH_nat), int_extra)), INT_FAM)
    tl <- unique(dat[, .(ID_time, year, month)])[order(ID_time)]
    rw <- function(fit, pre, expo) {
      r <- fit$summary.random$ID_time_rw[, c("mean", "0.025quant", "0.975quant")]
      if (expo) r <- exp(r)                                  # intensity: relative rate vs average
      setnames(as.data.table(r), paste0(pre, c("", "_lo", "_hi")))
    }
    wave <- cbind(tl, rw(nat_o, "occ_logodds", FALSE), rw(nat_i, "int_RR", TRUE))
    fwrite(wave, file.path(output_path, paste0("Fig3_national_trend", SUFFIX, ".csv")))
    cat("National model WAIC: occurrence", round(nat_o$waic$waic, 1), "| intensity", round(nat_i$waic$waic, 1), "\n")
    nat <- list(occ = nat_o, int = nat_i, trend = wave)
  }
  
  # ==============================================================================
  # 13. SAVE (Stage 2 reads this file)                                 [C3]
  # ==============================================================================
  set.seed(7); keep <- sample(NDRAW, NSAVE)                            # FIX: seed -> reproducible draws
  saveRDS(list(
    data             = dat,
    zone_month       = res,
    zone_summary     = zone,
    fits             = list(occ = fit_o, int = fit_i),
    metrics          = list(occ = met_occ, int = met_int),
    draws_for_stage2 = list(pi = PI[, keep], mu = MU[, keep], burden = B[, keep]),   # N x NSAVE, rows = dat order
    joint            = joint,
    national         = nat,
    graph_Q          = Q,
    zone_order       = shape$ADM2_PCODE_clean,
    settings         = list(r_nat = r_nat, min_E = min(dat$E), OCC_POP = OCC_POP, INT_POP = INT_POP, INT_ZT = INT_ZT, POP_SD = POP_SD,
                            row_order = "setorder(ID_time, ID_area)",
                            E_definition = "monthly_population x pooled national rate (2015-2024); no smoothing, no floor"),
    saved_at         = Sys.time()), OUT_RDS)
  # read_S(): read the saved state; if the file is unreadable (full disk, file locked by another session),
  # rebuild the same list from the objects in memory, so the later parts never depend on a damaged file.
  read_S <- function(path) {
    S_try <- tryCatch(readRDS(path), error = function(e) { message("Cannot read ", basename(path), ": ", conditionMessage(e), " -> using the objects in memory"); NULL })
    if (!is.null(S_try)) return(S_try)
    stopifnot(all(sapply(c("dat", "zone", "fit_o", "fit_i", "met_occ", "met_int", "Q", "shape", "PI", "MU", "B", "keep"), exists)))
    list(data = dat, zone_month = res, zone_summary = zone, fits = list(occ = fit_o, int = fit_i), metrics = list(occ = met_occ, int = met_int),
         draws_for_stage2 = list(pi = PI[, keep], mu = MU[, keep], burden = B[, keep]), joint = joint, national = nat,
         graph_Q = Q, zone_order = shape$ADM2_PCODE_clean)
  }
  cat("Saved state file:", OUT_RDS, "| readable:", !is.null(tryCatch(readRDS(OUT_RDS), error = function(e) NULL)), "\n")
  # all 10 candidate fits are saved one per file in the fits*/ folder (kept out of OUT_RDS to keep it small)
  cat("\nSaved:", OUT_RDS, "\n")
  
  
  # ##############################################################################
  # PART B — TEMPORAL RELATIVE RISK AND EPIDEMIC WAVES (main run only)
  # ##############################################################################
  if (FALSE) { cat("(unused)\n") } else {   # population in both parts is now the MAIN run, so Part B always runs
    
    # ==============================================================================
    # STAGE 1b — TEMPORAL RELATIVE RISK AND EPIDEMIC WAVES (new E)
    #
    # Reads the RDS saved at the end of PART A (runs automatically after it).
    #
    #   Primary   : selected intensity model (Baseline or Type I-IV)
    #   National  : sensitivity model, ONE shared RW2 (seasonal + BYM2 + RW2, no group)
    #
    # Changes from the old temporal_RR_FINAL.R (marked # FIX):
    #   - Both models use the full grid with NA response, so the RW2 is estimated
    #     for all 120 months. No 103-month time-ID mapping is needed.
    #   - CPO "caution years" REMOVED. A CPO failure flag means INLA's leave-one-out
    #     approximation was unreliable for that observation; it is not a sign of
    #     overfitting. Failures are counted in Stage 1 (table column CPO_failures).
    #   - "significant" -> "credibly elevated / credibly lower" (95% CrI excludes 1).
    #   - Log scale for RR (Juan: vertical scale compressed the estimates).
    #   - Wave labels "Wave 1", not "Wave: Wave 1" (minor comment, Figure 3).
    #   - Wave merge rule (gap <= 4 months) reported with a sensitivity check (2, 6).
    #   - Posterior sampling restricted to the time term (much less memory).
    #   - Season label: Belg = short rains (not harvest).
    # ==============================================================================
    
    library(INLA); library(data.table); library(ggplot2); library(scales)
    
    if (!exists("output_path")) output_path <- getOption("cholera.proj", ".")
    IN_RDS      <- if (exists("OUT_RDS")) OUT_RDS else file.path(output_path, "INLA_Rebuilt_From_Zero_OCCPOP.rds")   # the CURRENT run (never the old file without _OCCPOP)
    N_DRAWS     <- 2000
    MERGE_GAP   <- 4L          # elevated runs separated by <= 4 months = one wave
    GAP_SENS    <- c(2L, 6L)   # sensitivity of the wave count to this rule
    
    S       <- read_S(IN_RDS)
    dat     <- S$data
    fit_i   <- S$fits$int
    stopifnot(!is.null(S$national))
    fit_nat <- S$national$int
    n_time  <- max(dat$ID_time)
    tl      <- unique(dat[, .(t = ID_time, year, month)])[order(t)]
    cat("Primary intensity model:", S$metrics$int$Model[1], "\n")
    
    # FIX: CPO check is now a plain report, not a reliability filter
    cat(sprintf("CPO failures: primary %d | national %d\n",
                sum(fit_i$cpo$failure > 0, na.rm = TRUE), sum(fit_nat$cpo$failure > 0, na.rm = TRUE)))
    
    season_of <- function(m) fifelse(m %in% c(12, 1, 2), "Bega (dry)",
                                     fifelse(m %in% 3:5,         "Belg (short rains)",
                                             fifelse(m %in% 6:8,         "Kiremt (main rains)",
                                                     "Tseday (post-rains)")))
    SEASON_LEVELS <- c("Bega (dry)", "Belg (short rains)", "Kiremt (main rains)", "Tseday (post-rains)")
    season_cols   <- setNames(c("#FFF3C4", "#FCD9C8", "#D3E5F5", "#DCEFDC"), SEASON_LEVELS)
    
    # ==============================================================================
    # 1. MONTHLY RR FROM THE TEMPORAL TERM
    #    RR = exp(temporal effect) = intensity in that month relative to the
    #    2015-2024 average month, after seasonality, conditional on an outbreak.
    #    Primary (selected space-time model): per draw, main temporal effects
    #    (RW2 + iid) + the 92 zones' interaction averaged within each month,
    #    then exponentiated. For Types III/IV the interaction averages to 0 by
    #    construction, so the curve is the main temporal effect.
    # ==============================================================================
    n_area_b <- length(unique(dat$ID_area))
    get_rr <- function(fit, label) {
      sr     <- names(fit$summary.random)
      main_t <- intersect(c("ID_time_rw", "ID_time_iid"), sr)
      int_t  <- intersect("ID_st", sr)
      if (!length(main_t)) stop(label, ": model has no temporal term")
      if (length(main_t) == 1 && !length(int_t)) {               # one smooth curve: use marginals
        rw <- as.data.table(fit$summary.random[[main_t]])
        cat(sprintf("[%s] single curve (%s), %d months\n", label, main_t, nrow(rw)))
        return(data.table(t = rw$ID, RR = exp(rw$mean),
                          lo = exp(rw$`0.025quant`), hi = exp(rw$`0.975quant`)))
      }
      terms <- c(main_t, int_t)
      cat(sprintf("[%s] terms %s, %d draws...\n", label, paste(terms, collapse = " + "), N_DRAWS))
      sel <- setNames(lapply(terms, function(x) seq_len(nrow(fit$summary.random[[x]]))), terms)
      set.seed(42)
      sm  <- inla.posterior.sample(N_DRAWS, fit, selection = sel)
      rn  <- rownames(sm[[1]]$latent)
      idx <- lapply(setNames(terms, terms), function(x) grep(paste0("^", x, ":"), rn))
      lt  <- vapply(sm, function(s) {
        v <- Reduce(`+`, lapply(main_t, function(x) s$latent[idx[[x]], 1]))
        if (length(int_t)) v <- v + rowMeans(matrix(s$latent[idx$ID_st, 1], n_time, n_area_b))   # month within zone
        v
      }, numeric(n_time))
      RRd <- exp(lt)
      data.table(t = 1:n_time, RR = rowMeans(RRd),
                 lo = apply(RRd, 1, quantile, 0.025), hi = apply(RRd, 1, quantile, 0.975))
    }
    
    build <- function(fit, label) {
      d <- merge(tl, get_rr(fit, label), by = "t")
      d[, `:=`(model = label, season = factor(season_of(month), levels = SEASON_LEVELS),
               direction = fifelse(lo > 1, "credibly elevated",
                                   fifelse(hi < 1, "credibly lower", "includes 1")))]
      d[order(t)]
    }
    LAB_P <- "Selected model (primary)"
    LAB_N <- "National trend (sensitivity)"
    rr_p <- build(fit_i,   LAB_P)
    rr_n <- build(fit_nat, LAB_N)
    cat(sprintf("Months with 95%% CrI excluding 1: primary %d | national %d\n",
                rr_p[direction != "includes 1", .N], rr_n[direction != "includes 1", .N]))
    cat(sprintf("Median CrI width (hi/lo): primary %.1f | national %.1f  <- Juan: primary CrIs very wide\n",
                median(rr_p$hi / rr_p$lo), median(rr_n$hi / rr_n$lo)))
    
    # ==============================================================================
    # 2. TABLES: monthly (both models) and annual peaks
    # ==============================================================================
    fmt <- function(d, tag) d[, .(t, RR = sprintf("%.2f", RR), CrI = sprintf("%.2f-%.2f", lo, hi), Signal = direction)] |>
      setnames(c("RR", "CrI", "Signal"), paste0(c("RR", "95% CrI", "Signal"), " (", tag, ")"))
    monthly <- merge(tl[, .(t, Year = year, Month = month.abb[month])], fmt(rr_p, "primary"), by = "t") |>
      merge(fmt(rr_n, "national"), by = "t")
    fwrite(monthly[order(t)][, !"t"], file.path(output_path, "Table_temporalRR_monthly.csv"), bom = TRUE)
    
    peaks <- function(d) d[, .SD[which.max(RR)], by = year][
      , .(Year = year, `Peak month` = month.abb[month], `Peak RR` = sprintf("%.2f", RR),
          `95% CrI` = sprintf("%.2f-%.2f", lo, hi), `Excludes 1` = fifelse(direction == "includes 1", "No", "Yes"))]
    annual <- merge(peaks(rr_p), peaks(rr_n), by = "Year", suffixes = c(" (primary)", " (national)"))
    print(annual, row.names = FALSE)
    fwrite(annual, file.path(output_path, "Table_annualRR_combined.csv"), bom = TRUE)
    
    # ==============================================================================
    # 3. WAVES — national model only, detected from the data (not hard-coded)
    # ==============================================================================
    find_waves <- function(d, gap) {
      wv <- d[direction == "credibly elevated"][order(t)]
      if (!nrow(wv)) return(data.table())
      wv[, grp := cumsum(c(1, diff(t)) > gap)]
      w <- wv[, .(t_start = min(t), t_end = max(t),
                  From = sprintf("%s %d", month.abb[month[1]], year[1]),
                  To   = sprintf("%s %d", month.abb[month[.N]], year[.N]),
                  `Elevated months` = .N,
                  `Peak RR` = round(max(RR), 2),
                  `Peak 95% CrI` = sprintf("%.2f-%.2f", lo[which.max(RR)], hi[which.max(RR)]),
                  Peak = sprintf("%s %d", month.abb[month[which.max(RR)]], year[which.max(RR)])), by = grp]
      w[, Wave := .I][, grp := NULL][]
    }
    waves <- find_waves(rr_n, MERGE_GAP)
    cat("\n--- Epidemic waves, national model (gap <=", MERGE_GAP, "months) ---\n")
    print(waves[, !c("t_start", "t_end")], row.names = FALSE)
    fwrite(waves[, !c("t_start", "t_end")], file.path(output_path, "Table_waves.csv"), bom = TRUE)
    for (g in GAP_SENS) cat(sprintf("Sensitivity: gap <= %d months -> %d wave(s)\n", g, nrow(find_waves(rr_n, g))))
    
    # ==============================================================================
    # 4. FIGURE 3 — stacked panels, log scale, seasons, waves on national panel
    # ==============================================================================
    rr_plot <- rbind(rr_p, rr_n)[, model := factor(model, levels = c(LAB_P, LAB_N))]
    YMIN <- 0.05; YMAX <- 50                                     # outside -> squished to the edge
    bands <- rr_p[, runid := rleid(season)][, .(xmin = min(t) - 0.5, xmax = max(t) + 0.5, season = season[1]), by = runid]
    qlab  <- tl[month %in% c(1, 4, 7, 10)][, label := paste0(year, " Q", ceiling(month / 3))]
    WAVE_Y <- min(YMAX / 1.6, max(rr_n$hi) * 1.6)
    wave_dt <- if (nrow(waves)) waves[, .(t_start, t_end, t_mid = (t_start + t_end) / 2,
                                          lab = sprintf("Wave %d\npeak RR %.2f", Wave, `Peak RR`),
                                          model = factor(LAB_N, levels = levels(rr_plot$model)))] else NULL
    
    p <- ggplot(rr_plot, aes(x = t)) +
      geom_rect(data = bands, inherit.aes = FALSE,
                aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = season), alpha = 0.35) +   # FIX: log scale cannot show 0
      scale_fill_manual(values = season_cols, name = "Ethiopian season", drop = FALSE) +
      geom_vline(xintercept = tl[month == 1, t] - 0.5, colour = "grey70", linewidth = 0.4) +
      geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30") +
      geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.25) +
      geom_line(aes(y = RR), colour = "#08306B", linewidth = 0.8) +
      geom_point(data = rr_plot[direction != "includes 1"], aes(y = RR, colour = direction), size = 1.6) +
      scale_colour_manual(values = c("credibly elevated" = "#0072B2", "credibly lower" = "#D55E00"),
                          name = "95% CrI excludes 1") +
      facet_wrap(~ model, ncol = 1) +
      scale_x_continuous(breaks = qlab$t, labels = qlab$label, expand = c(0.01, 0)) +
      scale_y_log10(limits = c(YMIN, YMAX), oob = squish,
                    breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 8, 16, 32)) +           # FIX: log scale
      labs(x = NULL, y = "Relative risk of outbreak intensity (log scale)",
           caption = paste0("RR: intensity in each month relative to the 2015-2024 average month, after seasonality, ",
                            "conditional on an outbreak.\nPrimary: temporal main effects (RW2 + iid) of the selected model, plus the zone-averaged ",
                            "interaction if one is selected. Waves are identified from the national model only.")) +
      theme_bw(base_size = 11) +
      theme(strip.text = element_text(face = "bold"),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.5),
            legend.position = "bottom", plot.caption = element_text(hjust = 0, size = 7, colour = "grey40"))
    if (!is.null(wave_dt))
      p <- p +
      geom_segment(data = wave_dt, inherit.aes = FALSE,
                   aes(x = t_start, xend = t_end, y = WAVE_Y, yend = WAVE_Y), colour = "#7A0000", linewidth = 0.8,
                   arrow = arrow(ends = "both", length = unit(0.12, "cm"), type = "closed")) +
      geom_text(data = wave_dt, inherit.aes = FALSE, aes(x = t_mid, y = WAVE_Y * 1.5, label = lab),
                size = 2.8, colour = "#7A0000", lineheight = 0.85)                 # FIX: "Wave 1"
    print(p)
    
    f3 <- file.path(output_path, "Figure3_temporalRR")
    if (SAVE_FIGS) {
      ggsave(paste0(f3, ".png"),  p, width = 13, height = 9.5, dpi = 600, bg = "white")
      ggsave(paste0(f3, ".pdf"),  p, width = 13, height = 9.5, bg = "white", device = PDF_DEVICE)
      ggsave(paste0(f3, ".tiff"), p, width = 13, height = 9.5, dpi = 600, bg = "white", compression = "lzw")
      cat("Saved Figure 3 (PNG, PDF, TIFF) to", output_path, "\n")
    } else cat("Figure 3 drawn in the Plots pane (not saved; set SAVE_FIGS <- TRUE to save)\n")
    
  }  # end PART B
  
  
  # ##############################################################################
  # PART C — DIAGNOSTICS AND PRIOR SENSITIVITY FOR THE SELECTED MODEL (main run only)
  # ##############################################################################
  if (TRUE) {
    
    library(INLA); library(data.table); library(pROC)
    
    S      <- read_S(OUT_RDS)
    dat    <- S$data
    Q      <- S$graph_Q
    N      <- nrow(dat)
    n_area <- length(S$zone_order)
    
    sel_o <- S$metrics$occ[selected == TRUE, Model]
    sel_i <- S$metrics$int[selected == TRUE, Model]
    cat("Selected models: occurrence =", sel_o, "| intensity =", sel_i, "\n")
    fit_o <- readRDS(file.path(fit_dir_occ, paste0("occ_", sel_o, ".rds")))
    fit_i <- readRDS(file.path(fit_dir_int, paste0("int_", sel_i, ".rds")))
    
    # ==============================================================================
    # 1. OCCURRENCE: LEAVE-ONE-OUT PREDICTION OF OUTBREAK / NO OUTBREAK
    #    CPO_i = P(y_i | all other data). For a binary outcome the leave-one-out
    #    probability of an outbreak is CPO_i if y_i = 1 and 1 - CPO_i if y_i = 0.
    #    ("Leave one OBSERVATION out", not leave one zone or one month out.)
    # ==============================================================================
    cpo_o  <- fit_o$cpo$cpo
    fail_o <- fit_o$cpo$failure
    p_loo  <- ifelse(dat$Y_bin == 1, cpo_o, 1 - cpo_o)
    ok_o   <- (is.na(fail_o) | fail_o == 0) & is.finite(p_loo)
    cat(sprintf("\nOccurrence CPO failures: %d of %d\n", sum(fail_o > 0, na.rm = TRUE), N))
    
    auc_loo <- as.numeric(roc(dat$Y_bin[ok_o], p_loo[ok_o], quiet = TRUE)$auc)
    auc_in  <- as.numeric(roc(dat$Y_bin, fit_o$summary.fitted.values$mean[1:N], quiet = TRUE)$auc)
    prev    <- mean(dat$Y_bin[ok_o])
    brier   <- mean((dat$Y_bin[ok_o] - p_loo[ok_o])^2)
    brier0  <- mean((dat$Y_bin[ok_o] - prev)^2)             # always predict the overall frequency
    cat(sprintf("AUC in-sample %.3f | AUC leave-one-out %.3f\n", auc_in, auc_loo))
    cat(sprintf("Brier (LOO) %.4f vs %.4f for a model that only knows the overall frequency (skill %.1f%%)\n",
                brier, brier0, 100 * (1 - brier / brier0)))
    
    brk   <- unique(quantile(p_loo[ok_o], probs = seq(0, 1, 0.1)))
    calib <- data.table(bin = cut(p_loo[ok_o], brk, include.lowest = TRUE),
                        p = p_loo[ok_o], y = dat$Y_bin[ok_o])[
                          , .(n = .N, predicted = mean(p), observed = mean(y)), by = bin][order(predicted)]
    cat("\nCalibration of leave-one-out probabilities (predicted should match observed):\n")
    print(calib[, .(bin, n, predicted = round(predicted, 3), observed = round(observed, 3))])
    fwrite(calib, file.path(output_path, "Table_S_LOO_calibration_occurrence.csv"))
    fwrite(data.table(AUC_insample = auc_in, AUC_LOO = auc_loo, Brier = brier, Brier_null = brier0),
           file.path(output_path, "Table_S_LOO_AUC_occurrence.csv"))
    
    # ==============================================================================
    # 2. INTENSITY: PIT (probability integral transform) — calibration of the NB model
    #    Well-calibrated predictions give a roughly flat histogram. A pile-up near 1
    #    means observed counts are larger than predicted (the very large outbreaks).
    #    Observations with unreliable CPO are excluded.
    # ==============================================================================
    pit_i  <- fit_i$cpo$pit
    fail_i <- fit_i$cpo$failure
    if (is.null(pit_i) || all(is.na(pit_i))) {                    # INLA gives no leave-one-out PIT for this family
      cat("(INLA gives no leave-one-out PIT for this likelihood: using the randomised in-sample PIT of the zero-truncated NB.)\n")
      pit_i <- pit_zt(fit_i); fail_i <- rep(0, length(pit_i))
    }
    ok_i   <- !is.na(pit_i) & (is.na(fail_i) | fail_i == 0)
    cat(sprintf("\nIntensity PIT (n = %d reliable): share < 0.025 = %.3f | share > 0.975 = %.3f (0.025 expected each)\n",
                sum(ok_i), mean(pit_i[ok_i] < 0.025), mean(pit_i[ok_i] > 0.975)))
    draw_pit <- function() { hist(pit_i[ok_i], breaks = 20, freq = FALSE, col = "grey80", border = "white",
                                  main = "PIT, intensity model (Baseline)", xlab = "PIT value", ylab = "Density"); abline(h = 1, lty = 2) }
    draw_pit()                                                  # shown in the Plots pane
    if (SAVE_FIGS) { png(file.path(output_path, "Fig_S_PIT_intensity.png"), width = 1600, height = 1000, res = 200); draw_pit(); dev.off() }
    
    # ==============================================================================
    # 3. PRIOR SENSITIVITY — Baseline refitted with tighter / looser PC priors
    #    PC prior: Pr(sigma > U) = 0.01. "default" is the main analysis (loaded, not refitted).
    # ==============================================================================
    SEAS <- c("sin_month", "cos_month", "sin2_month", "cos2_month")
    fit_base <- function(response, family, U_bym, U_rw, U_iid) {
      pri_bym <- list(theta1 = list(prior = "pc.prec", param = c(U_bym, 0.01)),
                      theta2 = list(prior = "pc",      param = c(0.5, 0.5)))
      pri_rw  <- list(theta = list(prior = "pc.prec", param = c(U_rw, 0.01)))
      pri_iid <- list(theta = list(prior = "pc.prec", param = c(U_iid, 0.01)))
      rhs <- paste(paste(SEAS, collapse = " + "),
                   "+ f(ID_area_bym, model='bym2', graph=Q, scale.model=TRUE, constr=TRUE, hyper=pri_bym)",
                   "+ f(ID_time_rw, model='rw2', scale.model=TRUE, constr=TRUE, hyper=pri_rw)",
                   "+ f(ID_time_iid, model='iid', hyper=pri_iid)")
      lhs <- if (response == "O_int") "O_int ~ offset(logE) +" else "Y_bin ~"
      rhs <- paste(rhs, if (response == "O_int") int_extra else occ_extra)
      f <- as.formula(paste(lhs, rhs), env = environment())
      inla(f, family = family, data = dat, control.family = if (is_int_fam(family)) INT_CF else list(),
           control.predictor = list(compute = TRUE, link = 1),
           control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE),
           control.inla      = list(strategy = "simplified.laplace"))
    }
    settings <- data.table(setting = c("tight", "default", "loose"),
                           U_bym = c(0.5, 1, 2), U_rw = c(0.15, 0.3, 0.6), U_iid = c(0.25, 0.5, 1))
    zone_eff <- function(f) f$summary.random$ID_area_bym$mean[1:n_area]     # u_i + v_i
    phi_of   <- function(f) f$summary.hyperpar[grep("^Phi", rownames(f$summary.hyperpar))[1], "mean"]
    mlcpo    <- function(f) { ok <- is.na(f$cpo$failure) | f$cpo$failure == 0
    mean(log(f$cpo$cpo[ok]), na.rm = TRUE) }
    
    run_part <- function(part, response, family, fit_def) {
      base_eff <- zone_eff(fit_def)
      rbindlist(lapply(seq_len(nrow(settings)), function(k) {
        s <- settings[k]
        cat(sprintf("  %s | %s priors...\n", part, s$setting))
        U_rw_k <- (if (response == "O_int") U_RW_INT_EFF else U_RW_OCC_EFF) * s$U_rw / 0.3    # tight/default/loose = 0.5x/1x/2x of the RW2 bound used in the main model of that part
        f <- if (s$setting == "default") fit_def else fit_base(response, family, s$U_bym, U_rw_k, s$U_iid)
        e <- zone_eff(f)
        data.table(part = part, setting = s$setting, U_bym = s$U_bym, U_rw = U_rw_k, U_iid = s$U_iid,
                   phi = round(phi_of(f), 3), WAIC = round(f$waic$waic, 1), mean_logCPO = round(mlcpo(f), 4),
                   zone_effect_spearman_vs_default = round(cor(e, base_eff, method = "spearman"), 3),
                   top10_zones_shared_with_default = length(intersect(order(-e)[1:10], order(-base_eff)[1:10])))
      }))
    }
    if (sel_o == "Baseline" && sel_i == "Baseline") {
      cat("\nPrior sensitivity (refits the Baseline; a few minutes per fit)...\n")
      sens <- rbind(run_part("occurrence", "Y_bin", "binomial",   fit_o),
                    run_part("intensity",  "O_int", INT_FAM,      fit_i))
      cat("\n--- Prior sensitivity ---\n"); print(sens)
      fwrite(sens, file.path(output_path, "Table_S_prior_sensitivity.csv"))
      cat("\nIf Spearman is close to 1 and 8-10 of the top-10 zones are shared, the zone results\n",
          "do not depend on the prior. Report the range in the Supplement.\n")
    } else {
      cat("\nPrior sensitivity skipped: it is written for the Baseline model (selected:", sel_o, "/", sel_i, ").\n")
    }
    
  }  # end PART C
  
  
  # ##############################################################################
  # PART D — IS WHERE OUTBREAKS OCCUR DECOUPLED FROM HOW LARGE THEY ARE?
  #   Five complementary checks. No single one proves decoupling; read them together.
  #   D1  posterior correlation of the two zone effects (all draws, all zones)
  #   D2  same, only zones with enough positive months (shrinkage guard)
  #   D3  raw-data correlation (frequency vs relative size), bootstrap CI
  #   D4  ranking agreement: do the same zones top both lists? (vs chance)
  #   D5  joint-model beta under tight / default / loose priors
  # ##############################################################################
  if (RUN_DECOUPLE) {
    cat("\n================ PART D: decoupling ================\n")
    K_MIN  <- 5                                   # min. positive months for a zone's intensity to be informative
    npos   <- dat[, sum(O > 0), by = ID_area][order(ID_area)]$V1
    enough <- npos >= K_MIN
    cat(sprintf("Zones with >= %d positive months: %d of %d (the other %d are shrunk to their neighbours)\n",
                K_MIN, sum(enough), n_area, sum(!enough)))
    
    # ---- D1/D2: posterior draws of the zone effects (u_i + v_i) of the two SEPARATE fits.
    # The hurdle parts share no parameters, so their posteriors are independent and
    # pairing draw b of one with draw b of the other is a valid joint draw.
    draw_zone <- function(fit, seed) {
      set.seed(seed)
      s <- inla.posterior.sample(NDRAW, fit, selection = list(ID_area_bym = 1:n_area))
      sapply(s, function(x) x$latent[, 1])                    # n_area x NDRAW
    }
    Zo <- draw_zone(fit_o, 101); Zi <- draw_zone(fit_i, 102)
    rho_draw <- function(idx) vapply(seq_len(NDRAW), function(b) cor(Zo[idx, b], Zi[idx, b], method = "spearman"), numeric(1))
    rho_all <- rho_draw(rep(TRUE, n_area)); rho_k <- rho_draw(enough)
    summ <- function(r) data.table(median = round(median(r), 3), lo = round(unname(quantile(r, .025)), 3),
                                   hi = round(unname(quantile(r, .975)), 3),
                                   P_positive = round(mean(r > 0), 3), P_abs_below_0.3 = round(mean(abs(r) < 0.3), 3))
    D12 <- rbind(cbind(check = "D1 posterior Spearman, all zones",                 summ(rho_all)),
                 cbind(check = sprintf("D2 posterior Spearman, zones with >=%d positive months", K_MIN), summ(rho_k)))
    
    # ---- D3: raw data (no model). frequency = share of months with a case;
    # relative size = cases / E in the positive months only.
    zz <- dat[, .(freq = mean(O > 0), npos = sum(O > 0), rel_size = sum(O[O > 0]) / sum(E[O > 0])), by = ID_area][npos >= K_MIN]
    rho_emp <- cor(zz$freq, zz$rel_size, method = "spearman")
    set.seed(303)
    bt <- replicate(2000, { i <- sample(nrow(zz), replace = TRUE); cor(zz$freq[i], zz$rel_size[i], method = "spearman") })
    D3 <- data.table(check = sprintf("D3 raw data Spearman, %d zones (bootstrap CI)", nrow(zz)),
                     median = round(rho_emp, 3), lo = round(unname(quantile(bt, .025)), 3), hi = round(unname(quantile(bt, .975)), 3),
                     P_positive = round(mean(bt > 0), 3), P_abs_below_0.3 = round(mean(abs(bt) < 0.3), 3))
    
    # ---- D4: ranking agreement (posterior-mean zone effects). Expected overlap by chance
    # is q^2/m; p-value is hypergeometric.
    ze_o <- rowMeans(Zo); ze_i <- rowMeans(Zi)
    overlap <- function(idx, q, label) {
      ids <- which(idx); m <- length(ids)
      a <- ids[order(-ze_o[ids])[1:q]]; b <- ids[order(-ze_i[ids])[1:q]]; X <- length(intersect(a, b))
      data.table(zones = label, top_q = q, shared = X, expected_by_chance = round(q^2 / m, 1),
                 p_value = signif(phyper(X - 1, q, m - q, q, lower.tail = FALSE), 3))
    }
    D4 <- rbind(overlap(rep(TRUE, n_area), 10, "all"), overlap(rep(TRUE, n_area), 23, "all"),
                overlap(enough, 10, sprintf(">=%d positive months", K_MIN)))
    
    cat("\n--- D1-D3: correlation between 'how often' and 'how large' (zone level) ---\n"); print(rbind(D12, D3))
    cat("\n--- D4: top-q zones shared between the occurrence and intensity rankings ---\n"); print(D4)
    
    # ---- D6: the lists you SEE in the top-10 tables. Zone outbreak probability (pi) vs outbreak size (mu), with uncertainty ----
    # Three measures of outbreak size are compared with pi: (i) expected cases per outbreak month, absolute; (ii) the same per 100,000
    # people; (iii) annual expected caseload per 100,000 (pi x mu). NOTE: size per 100,000 is mechanically larger in SMALL zones when outbreak
    # size hardly depends on population (exponent about 0.1), so (ii) mostly ranks small-population zones.
    D6 <- NULL
    if (exists("PI") && exists("MU")) {
      zmean <- function(M, w = NULL) { X <- if (is.null(w)) M else M * w; rowsum(X, dat$ID_area) / as.vector(table(dat$ID_area)) }
      Pz <- zmean(PI); Mz <- zmean(MU); Mz100 <- zmean(MU, per100k_month)
      sp_d <- function(A, B) vapply(seq_len(ncol(A)), function(b) cor(A[, b], B[, b], method = "spearman"), numeric(1))
      ov_d <- function(A, B, q = 10) vapply(seq_len(ncol(A)), function(b) length(intersect(order(-A[, b])[1:q], order(-B[, b])[1:q])), numeric(1))
      one <- function(lab, sp_x, ov_x) data.table(comparison = lab, spearman_median = round(median(sp_x), 2), spearman_lo = round(quantile(sp_x, .025), 2),
                                                  spearman_hi = round(quantile(sp_x, .975), 2), top10_overlap_median = median(ov_x),
                                                  top10_overlap_lo = quantile(ov_x, .025), top10_overlap_hi = quantile(ov_x, .975))
      D6 <- rbind(one("pi vs expected cases per outbreak month (absolute)", sp_d(Pz, Mz), ov_d(Pz, Mz)),
                  one("pi vs expected cases per outbreak month per 100,000", sp_d(Pz, Mz100), ov_d(Pz, Mz100)),
                  one("pi vs size ranking by absolute cases: top 23 (25%)", sp_d(Pz, Mz), ov_d(Pz, Mz, 23)))
      cat("\n--- D6: zone outbreak probability vs outbreak size (posterior draws; top-10 overlap expected by chance = 1.09, P(0 overlap) = 0.30) ---\n")
      print(D6)
      cat(sprintf("Population vs size rank: Spearman(zone population, size per outbreak month) = %.2f ; (zone population, size per 100,000) = %.2f\n",
                  cor(dat[, .(p = mean(monthly_population)), by = ID_area][order(ID_area)]$p, rowMeans(Mz), method = "spearman"),
                  cor(dat[, .(p = mean(monthly_population)), by = ID_area][order(ID_area)]$p, rowMeans(Mz100), method = "spearman")))
      fwrite(D6, file.path(output_path, "Table_S_decoupling_descriptive_lists.csv"))
      cat("READ: zero overlap between two top-10 lists is NOT evidence of decoupling by itself: with 92 zones, two unrelated top-10 lists\n",
          "have no zone in common about 30% of the time. Judge by the correlation with its interval and by the overlap relative to 1.09.\n",
          "A per-100,000 size list is dominated by small zones with very wide intervals (see the mu_100k intervals), so it is not a\n",
          "like-for-like comparison with the occurrence list.\n")
    }
    
    # ---- D5: joint model, beta under three priors
    D5 <- NULL
    if (RUN_JOINT && exists("jd") && exists("fj")) {
      jfit_prior <- function(prec, ndraw = 300) {
        hb <- list(beta = list(prior = "normal", param = c(0, prec)))
        f2 <- fj; environment(f2) <- environment()            # so f() sees THIS hb
        m2 <- inla(f2, family = c("binomial", INT_FAM), data = jd, control.family = list(list(), INT_CF),
                   Ntrials = c(rep(1, n_o), rep(NA, n_i)),
                   control.predictor = list(compute = TRUE, link = c(rep(1, n_o), rep(2, n_i))),
                   control.compute = list(waic = TRUE, config = TRUE),
                   control.inla = list(strategy = "simplified.laplace"))
        b <- m2$summary.hyperpar[grep("^Beta", rownames(m2$summary.hyperpar))[1], ]
        set.seed(42); s2 <- inla.posterior.sample(ndraw, m2)
        rn <- rownames(s2[[1]]$latent)
        iw <- grep("^s_shared:", rn)[1:n_area]; iv <- grep("^s_own:", rn)[1:n_area]
        bi <- grep("^Beta", names(s2[[1]]$hyperpar))[1]
        r  <- vapply(s2, function(s) { w <- s$latent[iw]; cor(w, s$hyperpar[bi] * w + s$latent[iv]) }, numeric(1))
        q  <- quantile(r, c(.025, .5, .975))
        data.table(beta_prior_precision = prec, beta = round(b$mean, 3), beta_lo = round(b[["0.025quant"]], 3),
                   beta_hi = round(b[["0.975quant"]], 3), cor = round(q[2], 3), cor_lo = round(q[1], 3), cor_hi = round(q[3], 3))
      }
      cat("\nD5: refitting the joint model with beta prior precision 0.01, 0.1, 1 ...\n")
      D5 <- rbindlist(lapply(c(0.01, 0.1, 1), jfit_prior))
      cat("\n--- D5: joint-model beta and spatial correlation under three priors ---\n"); print(D5)
      fwrite(D5, file.path(output_path, "Table_S_decoupling_beta_prior.csv"))
    }
    
    # ---- D7: does the joint-model association depend on how POPULATION is handled? (2 x 2: occurrence with/without, intensity with/without) ----
    # Zone population drives occurrence (more people, more chances of a case) but barely drives outbreak size (exponent about 0.1). If the occurrence field
    # is left unadjusted, it carries population, and the shared-component beta can pick up population. The clean reading is the (population in BOTH parts)
    # row. Compare beta and the spatial correlation across the four rows; the exponent is the population exponent of outbreak size in that fit.
    D7 <- NULL
    if (RUN_JOINT && exists("jd") && exists("n_i")) {
      joint_pop <- function(occ_pop, int_pop, ndraw = 500) {
        jd2 <- jd; jd2$lp_o <- c(as.numeric(dat$log_pop_s), rep(0, n_i))
        f2 <- Y ~ -1 + int_o + s1o + c1o + s2o + c2o + int_i + s1i + c1i + s2i + c2i + offset(off) +
          f(s_shared, model = "bym2", graph = Q, scale.model = TRUE, constr = TRUE, hyper = prior_bym2) +
          f(s_copy, copy = "s_shared", fixed = FALSE, hyper = hb) +
          f(s_own, model = "bym2", graph = Q, scale.model = TRUE, constr = TRUE, hyper = prior_bym2)
        if (occ_pop) f2 <- update(f2, . ~ . + lp_o)
        if (int_pop) f2 <- update(f2, . ~ . + lp_i)
        environment(f2) <- environment()
        m2 <- inla(f2, family = c("binomial", INT_FAM), data = jd2, control.family = list(list(), INT_CF), Ntrials = c(rep(1, n_o), rep(NA, n_i)),
                   control.predictor = list(compute = TRUE, link = c(rep(1, n_o), rep(2, n_i))),
                   control.compute = list(waic = TRUE, config = TRUE), control.inla = list(strategy = "simplified.laplace"))
        b <- m2$summary.hyperpar[grep("^Beta", rownames(m2$summary.hyperpar))[1], ]
        set.seed(42); s2 <- inla.posterior.sample(ndraw, m2)
        rn <- rownames(s2[[1]]$latent); iw <- grep("^s_shared:", rn)[1:n_area]; iv <- grep("^s_own:", rn)[1:n_area]
        bi <- grep("^Beta", names(s2[[1]]$hyperpar))[1]
        r <- vapply(s2, function(z) { w <- z$latent[iw]; cor(w, z$hyperpar[bi] * w + z$latent[iv]) }, numeric(1)); q <- quantile(r, c(.025, .5, .975))
        ex <- if (int_pop && "lp_i" %in% rownames(m2$summary.fixed)) round(1 + m2$summary.fixed["lp_i", "mean"] / POP_SD, 2) else NA_real_
        data.table(occurrence_population = occ_pop, intensity_population = int_pop, beta = round(b$mean, 3), beta_lo = round(b[["0.025quant"]], 3),
                   beta_hi = round(b[["0.975quant"]], 3), cor = round(q[2], 3), cor_lo = round(q[1], 3), cor_hi = round(q[3], 3),
                   population_exponent_of_size = ex, WAIC = round(m2$waic$waic, 1))
      }
      cat("\nD7: joint model with and without a population term in each part (4 fits)...\n")
      D7 <- rbindlist(lapply(list(c(FALSE, FALSE), c(FALSE, TRUE), c(TRUE, FALSE), c(TRUE, TRUE)), function(v) joint_pop(v[1], v[2])))
      cat("\n--- D7: joint-model beta and spatial correlation by population adjustment ---\n"); print(D7)
      fwrite(D7, file.path(output_path, "Table_S_decoupling_population_adjustment.csv"))
      cat("READ: if beta and the correlation stay clearly above 0 in the (TRUE, TRUE) row, zones with more frequent outbreaks also tend to have\n",
          "larger outbreaks AFTER accounting for population: report an association, not decoupling. If the association appears only when population\n",
          "is handled in one part, say that the result depends on the population adjustment and report all four rows. Wording: association, not cause.\n")
    }
    fwrite(rbind(D12, D3), file.path(output_path, "Table_S_decoupling_correlations.csv"))
    fwrite(D4, file.path(output_path, "Table_S_decoupling_rank_overlap.csv"))
    
    cat("\nHOW TO READ PART D (decide from the numbers, do not pre-write the conclusion)\n",
        " * D1/D2/D3 intervals all include 0, D5 beta near 0, D4 overlap near chance -> 'no clear association'.\n",
        " * D1/D2/D5 intervals above 0 and D4 overlap above chance (check the top-23 row, not only top-10) -> 'modest positive association'.\n",
        " * If the model-based (D1, D2, D5) and raw-data (D3) results disagree, say so: the raw correlation is noisy for zones with few positive months.\n",
        " * Check D7 before writing: an association that appears only under one population adjustment depends on that adjustment.\n",
        " * Do not write 'decoupled' or 'independent', and do not write 'cause'.\n")
  }  # end PART D
  
  
  # ##############################################################################
  # PART E — SENSITIVITY TO THE 17 SILENT MONTHS AND A PERSISTENCE EXTENSION
  #   (occurrence Baseline; main analysis keeps the silent months)
  # ##############################################################################
  if ((RUN_SILENT || RUN_PERSIST || RUN_LAGINT) && sel_o == "Baseline" && sel_i == "Baseline") {
    cat("\n================ PART E ================\n")
    run_occ <- function(f, d) inla(f, family = "binomial", data = d,
                                   control.predictor = list(compute = TRUE, link = 1),
                                   control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                                   control.inla      = list(strategy = "simplified.laplace"))
    phi_e  <- function(x) x$summary.hyperpar[grep("^Phi", rownames(x$summary.hyperpar))[1], "mean"]
    zone_e <- function(x) x$summary.random$ID_area_bym$mean[1:n_area]
    nat_e  <- function(x) x$summary.random$ID_time_rw$mean
    mlcpo_e <- function(x) { ok <- is.na(x$cpo$failure) | x$cpo$failure == 0; mean(log(x$cpo$cpo[ok]), na.rm = TRUE) }
    
    if (RUN_SILENT) {
      silent <- dat[, .(cases = sum(O)), by = ID_time][cases == 0, ID_time]
      cat("Silent months (no case in any zone):", length(silent), "\n")
      dA <- copy(dat); dA[ID_time %in% silent, Y_bin := NA_integer_]
      cat("Fitting occurrence Baseline with silent months masked...\n")
      fitA <- run_occ(forms_occ[["Baseline"]], dA)
      sensA <- data.table(model = c("main (silent months kept)", "silent months masked"),
                          phi = round(c(phi_e(fit_o), phi_e(fitA)), 3))
      cat("\n--- E1. Silent months ---\n"); print(sensA)
      cat(sprintf("Zone effects: Spearman = %.3f | top-10 shared = %d of 10 | national trend (other months) r = %.3f\n",
                  cor(zone_e(fit_o), zone_e(fitA), method = "spearman"),
                  length(intersect(order(-zone_e(fit_o))[1:10], order(-zone_e(fitA))[1:10])),
                  cor(nat_e(fit_o)[setdiff(seq_len(n_time), silent)], nat_e(fitA)[setdiff(seq_len(n_time), silent)])))
      fwrite(sensA, file.path(output_path, "Table_S_silent_months_sensitivity.csv"))
    }
    
    if (RUN_PERSIST) {
      dB <- copy(dat); setorder(dB, ID_area, ID_time)
      dB[, lag_own := shift(Y_bin, 1L, fill = 0L), by = ID_area]
      Ymat <- matrix(dB$Y_bin, nrow = n_time, ncol = n_area)
      Ybef <- rbind(0, Ymat[-n_time, , drop = FALSE])
      nbshare <- as.matrix(Ybef %*% t(as.matrix(Q))) / rep(rowSums(as.matrix(Q)), each = n_time)
      dB[, lag_nb := as.vector(nbshare)]
      dB[, `:=`(lag_own = as.numeric(lag_own), lag_nb = as.numeric(scale(lag_nb)))]
      setorder(dB, ID_time, ID_area)
      stopifnot(identical(dB$ID_area, dat$ID_area), identical(dB$ID_time, dat$ID_time))
      cat("\nFitting occurrence Baseline + persistence terms...\n")
      fitB <- run_occ(update(forms_occ[["Baseline"]], . ~ . + lag_own + lag_nb), dB)
      cmpB <- data.table(model = c("Baseline", "Baseline + persistence"),
                         WAIC = round(c(fit_o$waic$waic, fitB$waic$waic), 1),
                         mean_logCPO = round(c(mlcpo_e(fit_o), mlcpo_e(fitB)), 5),
                         CPO_failures = c(sum(fit_o$cpo$failure > 0, na.rm = TRUE), sum(fitB$cpo$failure > 0, na.rm = TRUE)),
                         phi = round(c(phi_e(fit_o), phi_e(fitB)), 3))
      cat("\n--- E2. Persistence ---\n"); print(cmpB)
      print(round(fitB$summary.fixed[c("lag_own", "lag_nb"), c("mean", "0.025quant", "0.975quant")], 3))
      cat(sprintf("Odds ratio, outbreak in the zone last month: %.1f (95%% CrI %.1f-%.1f)\n",
                  exp(fitB$summary.fixed["lag_own", "mean"]), exp(fitB$summary.fixed["lag_own", "0.025quant"]),
                  exp(fitB$summary.fixed["lag_own", "0.975quant"])))
      fwrite(cmpB, file.path(output_path, "Table_S_persistence_model.csv"))
      cat("Better mean log-CPO by > 0.001 -> persistence improves PREDICTION; keep the Baseline as the\n",
          "descriptive model and report this as an extension (association / persistence, not transmission).\n")
    }
    
    # ---- E3. Does having an outbreak (or a big one) LAST month predict a LARGER count THIS month?
    # Intensity Baseline + (i) outbreak last month in the same zone (continuing vs newly started)
    # and (ii) size last month (log1p cases, standardised). Coefficients are on the log scale:
    # exp(coef) = ratio of expected cases. Association/prediction only, not cause.
    if (RUN_LAGINT) {
      dC <- copy(dat); setorder(dC, ID_area, ID_time)
      dC[, prev_cases := shift(O, 1L, fill = 0L), by = ID_area]
      dC[, `:=`(prev_occ = as.numeric(prev_cases > 0), prev_size = log1p(prev_cases))]
      dC[, prev_size := as.numeric(scale(prev_size))]
      setorder(dC, ID_time, ID_area)
      stopifnot(identical(dC$ID_area, dat$ID_area), identical(dC$ID_time, dat$ID_time))
      cat("\nFitting intensity Baseline + last-month terms (positive months only)...\n")
      fitC <- inla(update(forms_int[["Baseline"]], . ~ . + prev_occ + prev_size), family = INT_FAM, data = dC, control.family = INT_CF,
                   control.predictor = list(compute = TRUE, link = 1),
                   control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                   control.inla      = list(strategy = "simplified.laplace"))
      okc <- (is.na(fit_i$cpo$failure) | fit_i$cpo$failure == 0) & (is.na(fitC$cpo$failure) | fitC$cpo$failure == 0)
      cmpC <- data.table(model = c("Intensity Baseline", "Baseline + last-month terms"),
                         WAIC = round(c(fit_i$waic$waic, fitC$waic$waic), 1),
                         mean_logCPO_common = round(c(mean(log(fit_i$cpo$cpo[okc]), na.rm = TRUE),
                                                      mean(log(fitC$cpo$cpo[okc]), na.rm = TRUE)), 5),
                         n_positive_months = sum(!is.na(dat$O_int)),
                         n_continuing = sum(dC$prev_occ == 1 & !is.na(dC$O_int)))
      cat("\n--- E3. Last-month terms in the intensity model ---\n"); print(cmpC)
      kk <- fitC$summary.fixed[c("prev_occ", "prev_size"), c("mean", "0.025quant", "0.975quant")]
      print(round(kk, 3))
      cat(sprintf("Ratio of expected cases, continuing (outbreak last month) vs newly started: %.2f (95%% CrI %.2f-%.2f)\n",
                  exp(kk["prev_occ", "mean"]), exp(kk["prev_occ", "0.025quant"]), exp(kk["prev_occ", "0.975quant"])))
      fwrite(cbind(cmpC,
                   prev_occ_ratio = exp(kk["prev_occ", "mean"]), prev_occ_lo = exp(kk["prev_occ", "0.025quant"]), prev_occ_hi = exp(kk["prev_occ", "0.975quant"]),
                   prev_size_ratio_per_SD = exp(kk["prev_size", "mean"]), prev_size_lo = exp(kk["prev_size", "0.025quant"]), prev_size_hi = exp(kk["prev_size", "0.975quant"])),
             file.path(output_path, "Table_S_lagged_intensity.csv"))
      cat(sprintf("Size last month: ratio of expected cases per SD = %.2f (95%% CrI %.2f-%.2f)\n",
                  exp(kk["prev_size", "mean"]), exp(kk["prev_size", "0.025quant"]), exp(kk["prev_size", "0.975quant"])))
      cat("READ: 'continuing vs newly started' (prev_occ) and 'size last month' (prev_size) are separate questions.\n",
          "A ratio CrI around 1 for prev_occ = no extra size from merely having had an outbreak; a prev_size CrI above 1 with a log-CPO gain =\n",
          "larger counts follow larger counts in the same zone (within-zone persistence of size; association, not transmission).\n")
    }
  }  # end PART E
  
  
  # ##############################################################################
  # PART B2 — WAVES IN THE NATIONAL EXPECTED CASELOAD (occurrence x intensity)
  #   Part B finds waves from the national INTENSITY term only (size of an outbreak once it
  #   occurs). Years with many zones having outbreaks (e.g. 2016-17, 2020) can be missed by it.
  #   Here a wave = months in which the model's expected national caseload (sum over the
  #   92 zones of P(outbreak) x expected cases) is credibly above the 2015-2024 average month,
  #   computed draw by draw so the 95% CrI is proper. Same merge rule as Part B.
  # ##############################################################################
  if (exists("find_waves") && exists("MERGE_GAP")) {
    cat("\n================ PART B2: national expected-caseload waves ================\n")
    Bm <- if (exists("B")) B else readRDS(OUT_RDS)$draws_for_stage2$burden
    stopifnot(nrow(Bm) == nrow(dat))
    tot <- rowsum(Bm, dat$ID_time); tot <- tot[order(as.integer(rownames(tot))), , drop = FALSE]   # months x draws
    RRb <- sweep(tot, 2, colMeans(tot), "/")                          # relative to the average month, per draw
    tl2 <- unique(dat[, .(t = ID_time, year, month)])[order(t)]
    nb  <- data.table(tl2, RR = rowMeans(RRb), lo = apply(RRb, 1, quantile, .025), hi = apply(RRb, 1, quantile, .975),
                      observed_cases = dat[, .(c = sum(O)), by = ID_time][order(ID_time)]$c)
    nb[, direction := fifelse(lo > 1, "credibly elevated", fifelse(hi < 1, "credibly lower", "includes 1"))]
    fwrite(nb[, !"direction"], file.path(output_path, "Table_national_caseload_RR_monthly.csv"))
    wB <- find_waves(nb, MERGE_GAP)
    if (nrow(wB)) {
      wB[, `Observed cases in wave` := sapply(seq_len(.N), function(i) sum(nb[t >= t_start[i] & t <= t_end[i], observed_cases]))]
      wB[, `Share of all observed cases (%)` := round(100 * `Observed cases in wave` / sum(nb$observed_cases), 1)]
    }
    cat("\n--- Waves in national expected caseload (gap <=", MERGE_GAP, "months) ---\n")
    print(wB[, !c("t_start", "t_end")], row.names = FALSE)
    fwrite(wB[, !c("t_start", "t_end")], file.path(output_path, "Table_waves_national_caseload.csv"), bom = TRUE)
    for (g in GAP_SENS) cat(sprintf("Sensitivity: gap <= %d months -> %d wave(s)\n", g, nrow(find_waves(nb, g))))
    cat("Part B (intensity-only) asks when outbreaks were larger once they occurred; B2 asks when national\n",
        "caseload was high (many zones affected and/or larger outbreaks). Report B2 as the epidemic waves.\n")
  }  # end PART B2
  
  
  # ##############################################################################
  # PART F — SUMMARY REPORT (one file to share: Stage1_summary_report.md)
  #   Collects every table saved by Parts A-E into one Markdown report with the key
  #   numbers and a plain-language reading. Missing tables are skipped, so it works
  #   after a partial run and after the OCC_POP sensitivity run (file name has SUFFIX).
  # ##############################################################################
  cat("\n================ PART F: summary report ================\n")
  rd <- function(fname) { fp <- file.path(output_path, fname); if (file.exists(fp)) fread(fp) else NULL }
  md_table <- function(x) {
    if (is.null(x) || !nrow(x)) return("_not available (part not run)_\n")
    x <- as.data.frame(x); x[] <- lapply(x, function(v) if (is.numeric(v)) format(round(v, 4), trim = TRUE, scientific = FALSE) else as.character(v))
    paste0("| ", paste(names(x), collapse = " | "), " |\n",
           "|", paste(rep("---", ncol(x)), collapse = "|"), "|\n",
           paste(apply(x, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")), collapse = "\n"), "\n")
  }
  sec <- function(title, body) paste0("\n## ", title, "\n\n", body, "\n")
  
  R <- c(sprintf("# Stage 1 summary report\n\nGenerated: %s | run: %s | zones: %d | months: %d | zone-months: %d\n",
                 format(Sys.time(), "%Y-%m-%d %H:%M"), if (OCC_POP) "MAIN (population in both parts)" else "MAIN (population in intensity only)",
                 n_area, n_time, N))
  
  # 1 data and exposure
  R <- c(R, sec("1. Data and expected counts (E)", paste0(
    sprintf("- Cases: %d | zero zone-months: %.1f%% | months with cases: %d\n", sum(dat$O), 100 * mean(dat$O == 0), sum(dat$O > 0)),
    sprintf("- E = annual population projection (repeated in the 12 months of a year) x pooled national rate %.3e per person-month (about %.1f per 100,000 per year); no smoothing, no floor; minimum E = %.3f\n",
            r_nat, r_nat * 12 * 1e5, min(dat$E)),
    "- Zeros mean *no reported case*; 17 months with no case in any zone are kept in the main analysis (sensitivity in section 8).\n")))
  
  # 2 model comparison
  sel_txt <- sprintf("Selected by best leave-one-out log-CPO: occurrence = **%s**, intensity = **%s**. Where WAIC/DIC prefer an interaction model, p_WAIC > 0.4 flags WAIC as unreliable; the leave-one-out score is used.\n", sel_o, sel_i)
  R <- c(R, sec("2. Model comparison (Table 2)", paste0(sel_txt, "\n**Occurrence**\n\n", md_table(rd(paste0("Table2_occurrence_metrics", SUFFIX, ".csv"))),
                                                        "\n**Intensity**\n\n", md_table(rd(paste0("Table2_intensity_metrics", SUFFIX, ".csv"))))))
  
  # 3 parameters
  t3 <- rd(paste0("Table3_parameters", SUFFIX, ".csv"))
  R <- c(R, sec("3. Parameters of the selected models (Table 3)", paste0(md_table(t3),
                                                                         "\nphi = structured share of the spatial variance (not transmission or spillover). Intensity time-trend precision is imprecisely estimated.\n")))
  
  # 4 zone results
  zs <- rd(paste0("TableS3_zone_summary", SUFFIX, ".csv"))
  if (!is.null(zs)) {
    top <- function(col, k = 10) { z <- zs[order(-get(col))][1:k]; z[, .(zone, estimate = round(get(col), 3))] }
    R <- c(R, sec("4. Zone results (top 10)", paste0(
      "**Outbreak probability (per zone-month)**\n\n", md_table(top("pi")),
      "\n**Cases per 100,000 per outbreak month**\n\n", md_table(top("mu_100k")),
      "\n**Annual expected caseload per 100,000**\n\n", md_table(top("burden_100k")),
      sprintf("\nZones with no reported case in any month: %d; with 1-2 positive months: %d (estimates for these rely on neighbouring zones).\n",
              sum(zs$N_outbreak_months == 0), sum(zs$N_outbreak_months %in% 1:2)))))
  }
  
  # 5 decoupling
  dj <- rd(paste0("joint_model_beta_correlation", SUFFIX, ".csv"))
  R <- c(R, sec("5. Occurrence vs intensity (decoupling checks)", paste0(
    "**Joint shared-component model**\n\n", md_table(dj),
    "\n**D1-D3 zone-level correlations**\n\n", md_table(rd("Table_S_decoupling_correlations.csv")),
    "\n**D4 ranking overlap**\n\n", md_table(rd("Table_S_decoupling_rank_overlap.csv")),
    "\n**D5 beta under three priors**\n\n", md_table(rd("Table_S_decoupling_beta_prior.csv")),
    "\n**E3 outbreak last month and size this month**\n\n", md_table(rd("Table_S_lagged_intensity.csv")),
    "\nWording: read each row from its own interval. Intervals that include 0 mean *no clear association*; intervals above 0 mean a *positive association*. Check the population-adjustment table (D7) before writing, and never write decoupled or independent. Associations only, no causal claim.\n")))
  
  # 6 waves and time
  R <- c(R, sec("6. National temporal pattern and waves", paste0(
    "**Waves in national expected caseload (occurrence x intensity)**\n\n", md_table(rd("Table_waves_national_caseload.csv")),
    "\n**Waves in national outbreak size only (intensity term)**\n\n", md_table(rd("Table_waves.csv")),
    "\n**Annual relative risk**\n\n", md_table(rd("Table_annualRR_combined.csv")),
    if (!is.null(nat)) sprintf("\nNational-trend model WAIC: occurrence %.1f | intensity %.1f (compare with the selected models in section 2).\n",
                               nat$occ$waic$waic, nat$int$waic$waic) else "")))
  
  # 7 diagnostics
  R <- c(R, sec("7. Predictive diagnostics and prior sensitivity", paste0(
    "**Occurrence leave-one-out AUC and Brier**\n\n", md_table(rd("Table_S_LOO_AUC_occurrence.csv")),
    "\n**Calibration of leave-one-out probabilities**\n\n", md_table(rd("Table_S_LOO_calibration_occurrence.csv")),
    "\n**Prior sensitivity (tight / default / loose PC priors)**\n\n", md_table(rd("Table_S_prior_sensitivity.csv")),
    "\nIntensity PIT histogram: shown in the Plots pane (saved as Fig_S_PIT_intensity.png only if SAVE_FIGS is TRUE).\n")))
  
  # 8 sensitivity extensions
  R <- c(R, sec("8. Silent months and persistence", paste0(
    "**17 silent months masked**\n\n", md_table(rd("Table_S_silent_months_sensitivity.csv")),
    "\n**Persistence extension (occurrence)**\n\n", md_table(rd("Table_S_persistence_model.csv")))))
  
  R <- c(R, sec("9. Limits to state in the paper", paste0(
    "- Surveillance data: zero = no reported case, not proven absence.\n",
    "- All results are associations; no transmission, hub or spillover claims.\n",
    "- Intensity rests on the positive zone-months only (about 10% of zone-months), so its trend and phi are weakly identified.\n",
    "- No space-time interaction was selected by leave-one-out log-CPO. The Baseline does not model month-to-month persistence within zones (residual lag-1 correlation is large; see the persistence extension), so it describes average structure and is not a forecasting model.\n",
    "- Population is an annual projection repeated across months.\n")))
  
  R <- c(R, sec("10. Files", paste0("All tables and figures are in `", output_path, "`; one file per fitted model in `", basename(fit_dir_occ), "/` (occurrence) and `", basename(fit_dir_int), "/` (intensity).\n")))
  
  rep_file <- file.path(output_path, paste0("Stage1_summary_report", SUFFIX, ".md"))
  writeLines(paste(R, collapse = "\n"), rep_file, useBytes = TRUE)
  cat("Report written:", rep_file, "\n")
  
  
  ################################################################################
  # STEP 2 — SPATIAL AUTOCORRELATION AND LOCAL CLUSTERS
  ################################################################################
  # ==============================================================================
  # STAGE 1e — SPATIAL AUTOCORRELATION AND LOCAL CLUSTERS (replaces ArcMap)
  #   Global Moran's I, Local Moran's I (LISA) and Getis-Ord Gi* for the 92 zones,
  #   using the SAME neighbour graph as the INLA models (Q from the saved RDS).
  #
  # Run AFTER Stage1_ALL_corrected.R (needs INLA_Rebuilt_From_Zero.rds,
  # TableS3_zone_summary.csv and, for the model effects, the fits/ folder).
  #
  # Variables analysed (each is a different question):
  #   crude_inc   observed cases per 100,000 per year (raw; unstable in small zones)
  #   burden_100k model-based expected caseload per 100,000 per year (smoothed)
  #   pi          model-based outbreak probability per zone-month
  #   eff_occ     occurrence zone effect u_i+v_i (log-odds scale, Baseline fit)
  #   eff_int     intensity zone effect u_i+v_i (log scale, Baseline fit)
  # WORDING: "clusters of higher/lower estimated risk" - not transmission hotspots.
  # Local tests use 92 zones and many tests -> false-discovery-rate (FDR) adjusted.
  # ==============================================================================
  library(data.table); library(spdep); library(sf); library(ggplot2)
  
  output_path <- getOption("cholera.proj", ".")
  shp_path    <- "F:/Cholera/shape_22/Adm_2/eth_admbnda_adm2_csa_bofedb_2021.shp"
  NSIM <- 9999; set.seed(2026)
  
  # Inputs: dat (zone-month data), Q (neighbour graph), zone_order, zone summary.
  # 1) from the saved RDS; 2) if that file cannot be read, from the objects already in this
  #    R session (run this script in the SAME session right after Stage1_ALL_corrected.R).
  rds <- if (exists("OUT_RDS")) OUT_RDS else file.path(output_path, "INLA_Rebuilt_From_Zero_OCCPOP.rds")      # the CURRENT run
  cat("RDS size (MB):", if (file.exists(rds)) round(file.size(rds) / 1e6, 1) else "file not found", "\n")
  S <- tryCatch(readRDS(rds), error = function(e) { message("Cannot read RDS: ", conditionMessage(e)); NULL })
  if (!is.null(S)) {
    dat <- S$data; Q <- S$graph_Q; zone_order <- S$zone_order
    zs  <- if (exists("zone")) copy(zone) else fread(file.path(output_path, paste0("TableS3_zone_summary", if (exists("SUFFIX")) SUFFIX else "_OCCPOP", ".csv")))
  } else if (all(sapply(c("dat", "Q", "shape", "zone"), exists))) {
    cat("Using objects already in memory (dat, Q, shape, zone).\n")
    zone_order <- shape$ADM2_PCODE_clean
    zs <- copy(zone)
  } else stop("RDS unreadable and the main script's objects are not in memory. Re-run Stage1_ALL_corrected.R to the end.")
  n_area <- length(zone_order)
  setorder(zs, ID_area)
  stopifnot(nrow(zs) == n_area, nrow(Q) == n_area)
  
  # ---- neighbours = the model graph (queen contiguity + island links)
  nb    <- mat2listw(as.matrix(Q), style = "B")$neighbours
  lw_W  <- nb2listw(nb, style = "W", zero.policy = TRUE)             # Moran's I, LISA (row-standardised)
  lw_B  <- nb2listw(include.self(nb), style = "B", zero.policy = TRUE) # Gi* (binary, includes the zone itself)
  
  # ---- variables
  n_time <- max(dat$ID_time)
  crude <- dat[, .(cases = sum(O), pm = sum(monthly_population)), by = ID_area][order(ID_area)]
  X <- data.table(ID_area = 1:n_area,
                  crude_inc   = crude$cases / crude$pm * 12 * 1e5,
                  burden_100k = zs$burden_100k,
                  pi          = zs$pi)
  # model effects: the fits of THIS run (in memory), else the fit folders of the current specification (never the old "fits" folder)
  get_eff <- function(obj, path) {
    f <- if (!is.null(obj)) obj else if (file.exists(path)) readRDS(path) else NULL
    if (is.null(f)) NA_real_ else f$summary.random$ID_area_bym$mean[1:n_area] }
  fo_path <- if (exists("fit_dir_occ")) file.path(fit_dir_occ, paste0("occ_", if (exists("sel_o")) sel_o else "Baseline", ".rds")) else file.path(output_path, "fits_OCCPOP_rw20", "occ_Baseline.rds")
  fi_path <- if (exists("fit_dir_int")) file.path(fit_dir_int, paste0("int_", if (exists("sel_i")) sel_i else "Baseline", ".rds")) else file.path(output_path, "fits_INTPOP_rw10_ZT", "int_Baseline.rds")
  cat("Occurrence fit used:", if (exists("fit_o")) "in memory (this run)" else fo_path, "\nIntensity fit used:", if (exists("fit_i")) "in memory (this run)" else fi_path, "\n")
  X[, eff_occ := get_eff(if (exists("fit_o")) fit_o else NULL, fo_path)]
  X[, eff_int := get_eff(if (exists("fit_i")) fit_i else NULL, fi_path)]
  vars <- setdiff(names(X), "ID_area")
  
  # ==============================================================================
  # 1. GLOBAL MORAN'S I (is there spatial autocorrelation at all?)
  # ==============================================================================
  glob <- rbindlist(lapply(vars, function(v) {
    x  <- X[[v]]
    mt <- moran.test(x, lw_W, randomisation = TRUE, zero.policy = TRUE)
    mc <- moran.mc(x, lw_W, nsim = NSIM, zero.policy = TRUE)
    data.table(variable = v, Morans_I = round(unname(mt$estimate["Moran I statistic"]), 3),
               expected_I = round(unname(mt$estimate["Expectation"]), 3),
               z = round(unname(mt$statistic), 2), p_analytic = signif(mt$p.value, 3),
               p_permutation = signif((sum(mc$res >= mc$statistic) + 1) / (NSIM + 1), 3))
  }))
  cat("\n--- Global Moran's I (zone level, 92 zones) ---\n"); print(glob)
  fwrite(glob, file.path(output_path, "Table_S_global_morans_I.csv"))
  
  # Moran's I by year (is the pattern stable?), crude incidence per 100,000
  by_year <- rbindlist(lapply(sort(unique(dat$year)), function(y) {
    d <- dat[year == y, .(cases = sum(O), pm = sum(monthly_population)), by = ID_area][order(ID_area)]
    x <- d$cases / d$pm * 12 * 1e5
    if (sd(x) == 0) return(data.table(year = y, Morans_I = NA_real_, p_permutation = NA_real_))
    mc <- moran.mc(x, lw_W, nsim = 999, zero.policy = TRUE)
    data.table(year = y, Morans_I = round(unname(mc$statistic), 3),
               p_permutation = signif((sum(mc$res >= mc$statistic) + 1) / 1000, 3))
  }))
  cat("\n--- Moran's I by year (crude incidence) ---\n"); print(by_year)
  fwrite(by_year, file.path(output_path, "Table_S_morans_I_by_year.csv"))
  
  # ==============================================================================
  # 2. LOCAL MORAN'S I (LISA) AND GETIS-ORD Gi*
  # ==============================================================================
  # Gi* p-values by conditional permutation (the analytic z is unreliable for skewed variables with many zeros):
  # the zone's own value is held fixed and its neighbours' values are drawn at random from the other zones.
  gi_perm <- function(x, nb_self, nsim = 4999) {
    tot <- sum(x); n <- length(x)
    vapply(seq_len(n), function(i) {
      ids <- nb_self[[i]]; k <- length(setdiff(ids, i)); g_obs <- sum(x[ids]) / tot; xs <- x[-i]
      g_sim <- vapply(seq_len(nsim), function(s) (x[i] + sum(sample(xs, k))) / tot, numeric(1))
      (1 + sum(abs(g_sim - mean(g_sim)) >= abs(g_obs - mean(g_sim)))) / (nsim + 1)
    }, numeric(1))
  }
  local_one <- function(v) {
    x  <- X[[v]]; z <- as.numeric(scale(x)); lagz <- lag.listw(lw_W, z, zero.policy = TRUE)
    lm <- localmoran_perm(x, lw_W, nsim = NSIM, zero.policy = TRUE, iseed = 2026)
    pc <- grep("^Pr\\(z.*Sim", colnames(lm), value = TRUE)[1]
    if (is.na(pc)) pc <- grep("^Pr", colnames(lm), value = TRUE)[1]
    p  <- lm[, pc]; q <- p.adjust(p, "fdr")
    quad <- fifelse(z > 0 & lagz > 0, "High-High", fifelse(z < 0 & lagz < 0, "Low-Low",
                                                           fifelse(z > 0 & lagz < 0, "High-Low", "Low-High")))
    gz <- as.numeric(localG(x, lw_B, zero.policy = TRUE))           # Gi* z-scores (self included)
    gp <- gi_perm(x, lw_B$neighbours); gq <- p.adjust(gp, "fdr")
    data.table(ID_area = 1:n_area, zone = zone_order, variable = v, value = x,
               Ii = lm[, 1], p_local = p, q_local_FDR = q,
               lisa = fifelse(q < 0.05, quad, "Not significant"),
               Gi_star_z = gz, p_Gi = gp, q_Gi_FDR = gq,
               gi = fifelse(gq < 0.05 & gz > 0, "Hot (higher-value cluster)",
                            fifelse(gq < 0.05 & gz < 0, "Cold (lower-value cluster)", "Not significant")))
  }
  loc <- rbindlist(lapply(vars, local_one))
  fwrite(loc, file.path(output_path, "Table_S_local_clusters_LISA_Gi.csv"))
  
  cat("\n--- LISA clusters (FDR < 0.05) ---\n");  print(dcast(loc[, .N, by = .(variable, lisa)], variable ~ lisa, value.var = "N", fill = 0))
  cat("\n--- Gi* clusters (FDR < 0.05) ---\n");   print(dcast(loc[, .N, by = .(variable, gi)],   variable ~ gi,   value.var = "N", fill = 0))
  cat("\nZones in a High-High LISA cluster AND a Gi* hot cluster (burden_100k):\n")
  print(loc[variable == "burden_100k" & lisa == "High-High" & gi == "Hot (higher-value cluster)", zone])
  
  # overlap of occurrence and intensity clusters (supports the decoupling analysis)
  if (all(c("eff_occ", "eff_int") %in% vars)) {
    hh <- function(v) loc[variable == v & gi == "Hot (higher-value cluster)", ID_area]
    a <- hh("eff_occ"); b <- hh("eff_int")
    cat(sprintf("\nGi* hot zones: occurrence effect %d | intensity effect %d | in both %d\n",
                length(a), length(b), length(intersect(a, b))))
  }
  
  # ==============================================================================
  # 3. MAPS
  # ==============================================================================
  shape <- st_read(shp_path, quiet = TRUE)
  shape$zone <- toupper(gsub("[^A-Z0-9]", "", shape$ADM2_PCODE))
  shape <- shape[match(zone_order, shape$zone), ]
  stopifnot(!anyNA(shape$zone))
  for (v in intersect(c("burden_100k", "eff_occ"), vars)) {
    m  <- cbind(shape, loc[variable == v, .(lisa, gi)])
    p1 <- ggplot(m) + geom_sf(aes(fill = lisa), colour = "white", linewidth = 0.1) +
      scale_fill_manual(values = c("High-High" = "#B2182B", "Low-Low" = "#2166AC", "High-Low" = "#EF8A62",
                                   "Low-High" = "#67A9CF", "Not significant" = "grey88"), name = "LISA (FDR<0.05)") +
      labs(title = paste("Local Moran's I:", v)) + theme_void()
    p2 <- ggplot(m) + geom_sf(aes(fill = gi), colour = "white", linewidth = 0.1) +
      scale_fill_manual(values = c("Hot (higher-value cluster)" = "#B2182B", "Cold (lower-value cluster)" = "#2166AC",
                                   "Not significant" = "grey88"), name = "Gi* (FDR<0.05)") +
      labs(title = paste("Getis-Ord Gi*:", v)) + theme_void()
    print(p1); print(p2)                                       # shown in the Plots pane (use the arrows to page back)
    if (SAVE_FIGS) {
      ggsave(file.path(output_path, paste0("Fig_S_LISA_", v, ".png")), p1, width = 7, height = 6, dpi = 300, bg = "white")
      ggsave(file.path(output_path, paste0("Fig_S_Gi_",   v, ".png")), p2, width = 7, height = 6, dpi = 300, bg = "white")
    }
  }
  cat("\nSaved tables and maps in", output_path, "\n")
  cat("READ: significant global Moran's I = neighbouring zones have similar values. Local clusters are\n",
      "descriptive (areas of higher/lower estimated risk), FDR-adjusted; they are not transmission\n",
      "hotspots and not a substitute for the BYM2 model, which already estimates the spatial structure.\n")
  
  
  ################################################################################
  # STEP 3 — PRE-DEPLOYMENT CHECKS
  ################################################################################
  # ==============================================================================
  # STAGE 1f — PRE-DEPLOYMENT MODEL ASSESSMENT (run in the SAME R session, right after
  #            Stage1_ALL_corrected.R has finished: uses dat, Q, fit_o, fit_i, forms_occ,
  #            forms_int, PI/MU/B draws, nat, joint, output_path, ...)
  #
  #  1  INLA fit health (convergence flag, CPO failures)
  #  2  Prior-data conflict in the variance hyperparameters
  #  3  Zero-truncation of the intensity likelihood (is an untruncated NB acceptable?)
  #  4  Posterior predictive checks: national months, zone totals, zeros, extremes
  #  5  Residual autocorrelation: within zone (time), across zones (space), national series
  #  6  Out-of-time validation: train on early years, predict later years
  #  7  Stability of the zone ranking (probability of being in the top 10)
  #  8  Reproducibility record
  # Every check ends in PASS / FAIL / INFO in Table_predeployment_checks.csv. A FAIL is not
  # a crash: it tells you which statement the paper must NOT make, or what to fix.
  # Each section is wrapped so that one failing check does not stop the others.
  # ==============================================================================
  library(INLA); library(data.table); library(pROC); library(spdep)
  stopifnot(all(sapply(c("dat", "Q", "fit_o", "fit_i", "forms_occ", "forms_int", "output_path"), exists)))
  n_area <- length(unique(dat$ID_area)); n_time <- max(dat$ID_time); NN <- nrow(dat)
  
  chk <- data.table(check = character(), result = character(), status = character())
  add <- function(check, result, ok) {
    chk <<- rbind(chk, data.table(check = check, result = result,
                                  status = if (is.na(ok)) "INFO" else if (ok) "PASS" else "FAIL"))
    cat(sprintf("  -> [%s] %s: %s\n", chk[.N, status], check, result)); invisible(NULL)
  }
  step <- function(name, expr) tryCatch(expr, error = function(e) {
    cat("  [", name, "] could not run:", conditionMessage(e), "\n"); add(name, paste("error:", conditionMessage(e)), NA) })
  size_of <- function(f) f$summary.hyperpar[grep("^size", rownames(f$summary.hyperpar))[1], "mean"]
  lw <- nb2listw(mat2listw(as.matrix(Q))$neighbours, style = "W", zero.policy = TRUE)
  
  # ------------------------------------------------------------------------------
  cat("\n=== 1. INLA FIT HEALTH ===\n")
  step("1 INLA fit health", {
    fl <- list(occurrence = fit_o, intensity = fit_i)
    if (exists("nat") && !is.null(nat)) fl <- c(fl, list(national_occ = nat$occ, national_int = nat$int))
    if (exists("joint") && !is.null(joint)) fl$joint <- joint$fit
    H <- rbindlist(lapply(names(fl), function(n) { f <- fl[[n]]
    data.table(model = n, mode_status = f$mode$mode.status,
               cpo_failures = sum(f$cpo$failure > 0, na.rm = TRUE), n_obs_cpo = sum(!is.na(f$cpo$cpo)))}))
    H[, share_failed := round(cpo_failures / pmax(n_obs_cpo, 1), 4)]
    print(H)
    add("1 INLA mode.status = 0 (converged) for every model", paste(H$model[H$mode_status != 0], collapse = ", "), all(H$mode_status == 0))
    add("1 CPO failures < 2% (occurrence) and < 5% (intensity)",
        sprintf("occ %.2f%% | int %.2f%%", 100 * H[model == "occurrence", share_failed], 100 * H[model == "intensity", share_failed]),
        H[model == "occurrence", share_failed] < 0.02 && H[model == "intensity", share_failed] < 0.05)
  })
  
  # ------------------------------------------------------------------------------
  cat("\n=== 2. PRIOR-DATA CONFLICT (is the prior overriding, or being overridden by, the data?) ===\n")
  step("2 prior-data conflict", {
    uo <- U_RW_OCC_EFF; ui <- U_RW_INT_EFF
    ub_o <- if (exists("sel_occ")) sel_occ$U_bym else 1; uj_o <- if (exists("sel_occ")) sel_occ$U_iid else 0.5
    ub_i <- if (exists("sel_int")) sel_int$U_bym else 1; uj_i <- if (exists("sel_int")) sel_int$U_iid else 0.5
    terms <- rbindlist(list(
      data.table(part = "occurrence", pat = c("ID_area_bym", "ID_time_rw", "ID_time_iid"), U = c(ub_o, uo, uj_o), fit = list(fit_o)),
      data.table(part = "intensity",  pat = c("ID_area_bym", "ID_time_rw", "ID_time_iid"), U = c(ub_i, ui, uj_i), fit = list(fit_i))))
    out <- rbindlist(lapply(seq_len(nrow(terms)), function(i) {
      h <- terms$fit[[i]]$summary.hyperpar; r <- h[grep(paste0("^Precision for ", terms$pat[i]), rownames(h))[1], ]
      data.table(part = terms$part[i], term = terms$pat[i], sd_median = round(1 / sqrt(r[["0.5quant"]]), 3),
                 sd_lo = round(1 / sqrt(r[["0.975quant"]]), 3), sd_hi = round(1 / sqrt(r[["0.025quant"]]), 3),
                 prior_U = terms$U[i], sd_lo_above_U = (1 / sqrt(r[["0.975quant"]])) > terms$U[i])
    }))
    print(out)
    add("2 no posterior sd entirely above its prior bound U (prior-data conflict)",
        paste(out[sd_lo_above_U == TRUE, paste(part, term)], collapse = "; "), !any(out$sd_lo_above_U))
    cat("  A conflict means the prior says 'sd above U is very unlikely' while the data say the sd is above U.\n",
        "  Fix: raise U (U_RW_OCC / U_RW_INT in the main script) and refit; report the sensitivity.\n")
  })
  
  # ------------------------------------------------------------------------------
  cat("\n=== 3. ZERO-TRUNCATION OF THE INTENSITY LIKELIHOOD ===\n")
  step("3 zero-truncation", {
    pos <- which(!is.na(dat$O_int)); sz <- size_of(fit_i); mu <- mu_nb(fit_i, pos)
    p0 <- (sz / (sz + mu))^sz
    cat(sprintf("  NB probability of a zero at the positive months: median %.3f | mean %.3f | max %.3f\n",
                median(p0), mean(p0), max(p0)))
    cat(sprintf("  Share of positive months with P(0) > 5%%: %.1f%% | mean/truncated-mean ratio (median) %.3f\n",
                100 * mean(p0 > 0.05), median(1 - p0)))
    if (is_zt_fit(fit_i)) {
      add("3 intensity: zero-truncated NB likelihood (zero is excluded by construction)",
          sprintf("truncated NB; the untruncated P(0) (median %.3f) is removed from the likelihood", median(p0)), TRUE)
    } else {
      add("3 intensity NB: median P(0) at positive months < 5% (untruncated NB acceptable)",
          sprintf("median %.3f, %.0f%% of months > 5%%", median(p0), 100 * mean(p0 > 0.05)), median(p0) < 0.05)
    }
    cat("  If this FAILS: the fitted mean mu is the mean of an NB that puts mass on zero, but the data are positive.\n",
        "  The conditional expected caseload is mu / (1 - P0), and the PIT histogram lacks small values.\n",
        "  Options: report B = pi * mu/(1 - P0); or use a zero-truncated likelihood if your INLA build has one\n",
        "  (check names(inla.models()$likelihood) for 'nz'/'trunc'); or state it as a limitation.\n")
  })
  
  # ------------------------------------------------------------------------------
  cat("\n=== 4. POSTERIOR PREDICTIVE CHECKS ===\n")
  Tm_nat <- NULL
  step("4 posterior predictive", {
    # the predictive draws need the NB MEAN PARAMETER (MU_NB, before truncation); MU is E(Y | Y > 0) when the intensity is a truncated NB
    if (exists("PI") && exists("MU_NB")) { Pm <- PI[, seq_len(ncol(MU_NB)), drop = FALSE]; Mm <- MU_NB } else if (exists("PI") && exists("MU") && !INT_ZT) { Pm <- PI; Mm <- MU } else {
      if (INT_ZT) stop("MU_NB is not in memory: run Step 1 in this R session (the saved file keeps only E(Y | Y > 0))")
      SS <- readRDS(OUT_RDS); Pm <- SS$draws_for_stage2$pi; Mm <- SS$draws_for_stage2$mu }
    R <- min(500, ncol(Pm), ncol(Mm)); set.seed(99); sz <- size_of(fit_i)
    Yrep <- matrix(0, NN, R)
    for (b in 1:R) {
      mu <- Mm[, b]; p0 <- (sz / (sz + mu))^sz
      y  <- qnbinom(pmin(p0 + runif(NN) * (1 - p0), 1 - 1e-12), size = sz, mu = mu)   # zero-truncated NB draw
      Yrep[, b] <- rbinom(NN, 1, Pm[, b]) * y
    }
    q <- function(M) t(apply(M, 1, quantile, c(.025, .975)))
    tm <- rowsum(Yrep, dat$ID_time); obs_t <- dat[, .(c = sum(O)), by = ID_time][order(ID_time)]$c; qt <- q(tm)
    zt <- rowsum(Yrep, dat$ID_area); obs_z <- dat[, .(c = sum(O)), by = ID_area][order(ID_area)]$c; qz <- q(zt)
    Tm_nat <<- tm
    cov_t <- mean(obs_t >= qt[, 1] & obs_t <= qt[, 2]); cov_z <- mean(obs_z >= qz[, 1] & obs_z <= qz[, 2])
    zs <- colMeans(Yrep == 0); zq <- quantile(zs, c(.025, .975)); mx <- apply(Yrep, 2, max); mq <- quantile(mx, c(.025, .975))
    tot <- colSums(Yrep); tq <- quantile(tot, c(.025, .975))
    add("4 national monthly cases: observed within the 95% predictive interval in >= 90% of months",
        sprintf("%.1f%% of %d months", 100 * cov_t, n_time), cov_t >= 0.90)
    add("4 zone total cases: observed within the 95% predictive interval in >= 90% of zones",
        sprintf("%.1f%% of %d zones", 100 * cov_z, n_area), cov_z >= 0.90)
    add("4 share of zero zone-months: observed inside the 95% predictive interval",
        sprintf("observed %.3f | predicted %.3f-%.3f", mean(dat$O == 0), zq[1], zq[2]), mean(dat$O == 0) >= zq[1] && mean(dat$O == 0) <= zq[2])
    add("4 total cases: observed inside the 95% predictive interval",
        sprintf("observed %d | predicted %.0f-%.0f", sum(dat$O), tq[1], tq[2]), sum(dat$O) >= tq[1] && sum(dat$O) <= tq[2])
    add("4 largest zone-month count: observed inside the 95% predictive interval",
        sprintf("observed %d | predicted %.0f-%.0f", max(dat$O), mq[1], mq[2]), max(dat$O) >= mq[1] && max(dat$O) <= mq[2])
    draw_ppc <- function() {
      matplot(1:n_time, log1p(cbind(qt, obs_t)), type = "n", xlab = "Month index (1 = Jan 2015)", ylab = "log(1 + cases)",
              main = "National monthly cases: observed vs 95% predictive interval")
      polygon(c(1:n_time, n_time:1), c(log1p(qt[, 1]), rev(log1p(qt[, 2]))), col = "grey85", border = NA)
      lines(1:n_time, log1p(obs_t), lwd = 2) }
    draw_ppc()                                                 # shown in the Plots pane
    if (SAVE_FIGS) { png(file.path(output_path, "Fig_S_posterior_predictive_national.png"), width = 1800, height = 900, res = 200); draw_ppc(); dev.off() }
  })
  
  # ------------------------------------------------------------------------------
  cat("\n=== 5. RESIDUAL AUTOCORRELATION (what the model has NOT captured) ===\n")
  step("5 residuals: occurrence, time within zone", {
    p <- fit_o$summary.fitted.values$mean[1:NN]
    dd <- data.table(a = dat$ID_area, t = dat$ID_time, r = (dat$Y_bin - p) / sqrt(p * (1 - p))); setorder(dd, a, t)
    dd[, r_prev := shift(r), by = a]; rho <- cor(dd$r, dd$r_prev, use = "complete.obs")
    add("5 occurrence: lag-1 residual correlation within zones |r| < 0.10", sprintf("r = %.3f", rho), abs(rho) < 0.10)
    cat("  A clearly positive value = outbreaks persist from one month to the next beyond what the Baseline explains\n",
        "  (the persistence extension in Part E addresses this). The Baseline then describes risk but forecasts poorly.\n")
  })
  step("5 residuals: occurrence, space", {
    p <- fit_o$summary.fitted.values$mean[1:NN]
    zr <- tapply(dat$Y_bin - p, dat$ID_area, sum) / sqrt(tapply(p * (1 - p), dat$ID_area, sum))
    mt <- moran.test(as.numeric(zr), lw, zero.policy = TRUE)
    add("5 occurrence: zone-level residuals not spatially autocorrelated (Moran p > 0.05)",
        sprintf("I = %.3f, p = %.3f", mt$estimate[["Moran I statistic"]], mt$p.value), mt$p.value > 0.05)
  })
  step("5 residuals: intensity, space", {
    r <- ifelse(is.na(dat$O_int), NA_real_, pearson_int(fit_i, dat$O_int, 1:NN))      # Pearson residual under the fitted (truncated) likelihood
    zr <- tapply(r, dat$ID_area, function(x) { x <- x[!is.na(x)]; if (length(x) >= 3) sum(x) / sqrt(length(x)) else NA_real_ })   # zones with >= 3 outbreak months
    keep <- !is.na(zr)
    lw_k <- spdep::subset.listw(lw, keep, zero.policy = TRUE)
    mt <- moran.test(as.numeric(zr[keep]), lw_k, zero.policy = TRUE)
    add("5 intensity: zone-level residuals not spatially autocorrelated (Moran p > 0.05)",
        sprintf("I = %.3f, p = %.3f (%d zones with >= 3 outbreak months)", mt$estimate[["Moran I statistic"]], mt$p.value, sum(keep)), mt$p.value > 0.05)
  })
  step("5 residuals: intensity, time within zone", {
    r <- ifelse(is.na(dat$O_int), NA_real_, pearson_int(fit_i, dat$O_int, 1:NN))
    dd <- data.table(a = dat$ID_area, t = dat$ID_time, r = r); setorder(dd, a, t); dd[, r_prev := shift(r), by = a]
    ok <- complete.cases(dd[, .(r, r_prev)]); rho <- cor(dd$r[ok], dd$r_prev[ok])
    add("5 intensity: lag-1 residual correlation between consecutive outbreak months |r| < 0.10",
        sprintf("r = %.3f (%d pairs)", rho, sum(ok)), abs(rho) < 0.10)
  })
  step("5 residuals: national series", {
    stopifnot(!is.null(Tm_nat))
    obs_t <- dat[, .(c = sum(O)), by = ID_time][order(ID_time)]$c
    z <- (log1p(obs_t) - rowMeans(log1p(Tm_nat))) / apply(log1p(Tm_nat), 1, sd)
    a1 <- cor(z[-1], z[-length(z)]); bt <- Box.test(z, lag = 12, type = "Ljung-Box")
    add("5 national monthly residuals: no autocorrelation (Ljung-Box p > 0.05)",
        sprintf("lag-1 r = %.2f, p = %.3g", a1, bt$p.value), bt$p.value > 0.05)
  })
  
  # ------------------------------------------------------------------------------
  cat("\n=== 6. OUT-OF-TIME VALIDATION (train on early years, predict later years) ===\n")
  lag_dat <- function(d) {            # lags built from OBSERVED earlier months: one-step-ahead forecasting
    d <- copy(d); setorder(d, ID_area, ID_time)
    d[, prev_cases := shift(O, 1L, fill = 0L), by = ID_area]
    d[, `:=`(lag_own = as.numeric(prev_cases > 0), prev_occ = as.numeric(prev_cases > 0), prev_size = log1p(prev_cases))]
    Ymat <- matrix(as.numeric(d$O > 0), nrow = n_time, ncol = n_area)
    Ybef <- rbind(0, Ymat[-n_time, , drop = FALSE])
    d[, lag_nb := as.vector(as.matrix(Ybef %*% t(as.matrix(Q))) / rep(rowSums(as.matrix(Q)), each = n_time))]
    setorder(d, ID_time, ID_area); d[]
  }
  run_f <- function(f, dd, fam) inla(f, family = fam, data = dd, control.family = if (is_int_fam(fam)) INT_CF else list(),
                                     control.predictor = list(compute = TRUE, link = 1),
                                     control.compute = list(waic = FALSE, cpo = FALSE),
                                     control.inla = list(strategy = "simplified.laplace"))
  fwd <- function(train_end) {
    d <- lag_dat(dat); test <- d$year > train_end
    y <- d$Y_bin[test]; prev <- mean(d$Y_bin[!test])
    d_o <- copy(d); d_o[test, Y_bin := NA_integer_]
    f1 <- run_f(forms_occ[["Baseline"]], d_o, "binomial")
    f2 <- run_f(update(forms_occ[["Baseline"]], . ~ . + lag_own + lag_nb), d_o, "binomial")
    sc <- function(p) c(AUC = as.numeric(roc(y, p, quiet = TRUE)$auc), Brier = mean((y - p)^2),
                        LogScore = -mean(y * log(p) + (1 - y) * log(1 - p)))
    p1 <- f1$summary.fitted.values$mean[test]; p2 <- f2$summary.fitted.values$mean[test]
    null <- c(AUC = 0.5, Brier = mean((y - prev)^2), LogScore = -mean(y * log(prev) + (1 - y) * log(1 - prev)))
    occ <- rbind(data.table(window = sprintf("train<=%d, test %d-%d", train_end, train_end + 1, max(d$year)), component = "occurrence",
                            model = c("Climatology (training frequency)", "Baseline", "Baseline + persistence"),
                            rbind(null, sc(p1), sc(p2))))
    # intensity: positive months in the test period
    d_i <- copy(d); d_i[test, O_int := NA_integer_]
    g1 <- run_f(forms_int[["Baseline"]], d_i, INT_FAM)
    g2 <- run_f(update(forms_int[["Baseline"]], . ~ . + prev_occ + prev_size), d_i, INT_FAM)
    tp <- test & !is.na(d$O_int); yi <- d$O_int[tp]
    isc <- function(g) { sz <- size_of(g); mu <- mu_nb(g, which(tp)); p0 <- (sz / (sz + mu))^sz
    c(LogScore = -mean(dnbinom(yi, size = sz, mu = mu, log = TRUE) - log1p(-p0)),
      MedianAbsErr_log = median(abs(log(yi) - log(mu / (1 - p0))))) }
    int <- rbind(data.table(window = occ$window[1], component = "intensity", model = c("Baseline", "Baseline + last-month terms"),
                            rbind(isc(g1), isc(g2)), n_test = length(yi)), fill = TRUE)
    rbind(occ, int, fill = TRUE)
  }
  step("6 out-of-time validation", {
    fv <- rbindlist(list(fwd(2018), fwd(2021)), fill = TRUE)
    print(fv)
    fwrite(fv, file.path(output_path, "Table_S_out_of_time_validation.csv"))
    b <- fv[component == "occurrence" & model == "Baseline"]; n0 <- fv[component == "occurrence" & grepl("Climatology", model)]
    add("6 occurrence Baseline: out-of-time AUC >= 0.75 in every window", paste(round(b$AUC, 3), collapse = " / "), all(b$AUC >= 0.75))
    add("6 occurrence Baseline: out-of-time Brier better than climatology in every window",
        paste(sprintf("%.4f vs %.4f", b$Brier, n0$Brier), collapse = " / "), all(b$Brier < n0$Brier))
    pp <- fv[component == "occurrence" & model == "Baseline + persistence"]
    add("6 persistence terms improve out-of-time occurrence prediction (lower log score)",
        paste(sprintf("%.4f vs %.4f", pp$LogScore, b$LogScore), collapse = " / "), all(pp$LogScore < b$LogScore))
    ii <- fv[component == "intensity"]
    cat("  Intensity out-of-time results are INFO: compare Baseline with last-month terms in the table above.\n")
  })
  
  # ------------------------------------------------------------------------------
  cat("\n=== 7. STABILITY OF THE ZONE RANKING ===\n")
  step("7 rank stability", {
    Bm <- if (exists("B")) B else readRDS(OUT_RDS)$draws_for_stage2$burden
    Bz <- rowsum(Bm * (12 * 1e5 / dat$monthly_population), dat$ID_area) / n_time          # zone burden per 100,000 per year, per draw
    top <- apply(Bz, 2, function(x) order(-x)[1:10])
    pr  <- table(factor(as.vector(top), levels = seq_len(n_area))) / ncol(Bz)
    nm  <- if (exists("zone")) zone$zone[order(zone$ID_area)] else paste0("zone", seq_len(n_area))
    tab <- data.table(zone = nm, P_in_top10 = round(as.numeric(pr), 3))[order(-P_in_top10)][1:15]
    print(tab); fwrite(tab, file.path(output_path, "Table_S_rank_stability_top10.csv"))
    add("7 number of zones with >= 80% posterior probability of being in the top 10 by burden",
        sprintf("%d of 10 (>=50%%: %d)", sum(pr >= 0.8), sum(pr >= 0.5)), NA)
    cat("  Report the ranking with these probabilities; do not present a point ranking as fact.\n")
  })
  
  # ------------------------------------------------------------------------------
  cat("\n=== 8. REPRODUCIBILITY RECORD ===\n")
  step("8 run info", {
    info <- c(sprintf("Run: %s", format(Sys.time(), "%Y-%m-%d %H:%M")), R.version.string,
              sprintf("INLA %s | spdep %s | data.table %s", as.character(packageVersion("INLA")),
                      as.character(packageVersion("spdep")), as.character(packageVersion("data.table"))),
              sprintf("zones %d | months %d | zone-months %d | cases %d", n_area, n_time, NN, sum(dat$O)),
              sprintf("OCC_POP=%s INT_POP=%s U_RW_OCC=%s U_RW_INT=%s NDRAW=%s",
                      if (exists("OCC_POP")) OCC_POP else NA, if (exists("INT_POP")) INT_POP else NA,
                      if (exists("U_RW_OCC")) U_RW_OCC else NA, if (exists("U_RW_INT")) U_RW_INT else NA,
                      if (exists("NDRAW")) NDRAW else NA),
              "Seeds: draws 2024/2025, stage-2 subsample 7, joint 42, predictive checks 99",
              "", capture.output(sessionInfo()))
    writeLines(info, file.path(output_path, "Stage1_run_info.txt"))
    add("8 run information saved (R/INLA versions, settings, seeds)", "Stage1_run_info.txt", TRUE)
  })
  
  # ------------------------------------------------------------------------------
  cat("\n================ PRE-DEPLOYMENT SUMMARY ================\n")
  print(chk, row.names = FALSE)
  fwrite(chk, file.path(output_path, "Table_predeployment_checks.csv"))
  cat(sprintf("\nPASS %d | FAIL %d | INFO %d\n", sum(chk$status == "PASS"), sum(chk$status == "FAIL"), sum(chk$status == "INFO")))
  cat("Deploy the model only for statements whose checks PASS. For every FAIL, either fix and refit, or write the\n",
      "limitation into the paper (see the notes printed under each section).\n")
  
  
  
  ################################################################################
  # STEP 3B — OVERFITTING CHECK FOR ALL FIVE CANDIDATE MODELS (Baseline, Type I-IV)
  #   (a) complexity: effective number of parameters (p_WAIC) against the number of observations
  #   (b) random-effect SDs: is an interaction term soaking up the variance?
  #   (c) HELD-OUT validation. Two splits, each model refitted with the test responses set to NA:
  #         random20 : 20% of zone-months hidden at random (optimistic: neighbours/months stay in training)
  #         forward  : train 2015-2021, predict 2022-2024 (the honest test for a surveillance model)
  #       In-sample vs held-out AUC, log loss and Brier (occurrence); truncated log score (intensity).
  #       A large drop from in-sample to held-out = overfitting. Fast INLA setting (Gaussian approximation).
  #   The PIT histogram (Part C) is the calibration check; this step is the generalisation check.
  ################################################################################
  RUN_OVERFIT    <- TRUE
  OVERFIT_MODELS <- names(forms_occ)        # drop slow ones if needed, e.g. c("Baseline", "Type_I", "Type_II")
  if (RUN_OVERFIT) {
    cat("\n================ STEP 3B: overfitting checks (Baseline vs interaction models) ================\n")
    # ---- (a) complexity
    cx <- rbindlist(list(
      data.table(part = "occurrence", model = names(fits_occ), p_WAIC = sapply(fits_occ, function(f) f$waic$p.eff), n_obs = sum(!is.na(dat$Y_bin))),
      data.table(part = "intensity",  model = names(fits_int), p_WAIC = sapply(fits_int, function(f) f$waic$p.eff), n_obs = sum(!is.na(dat$O_int)))))
    cx[, `:=`(p_WAIC = round(p_WAIC, 1), p_WAIC_share_of_n = round(p_WAIC / n_obs, 3))]
    cat("\n(a) Effective number of parameters (p_WAIC) relative to the number of observations:\n"); print(cx)
    fwrite(cx, file.path(output_path, "Table_S_overfit_complexity.csv"))
    # ---- (b) random-effect standard deviations
    sd_tab <- rbindlist(lapply(c("occurrence", "intensity"), function(part) {
      fl <- if (part == "occurrence") fits_occ else fits_int
      rbindlist(lapply(names(fl), function(m) { h <- fl[[m]]$summary.hyperpar; i <- grep("^Precision for", rownames(h))
      if (!length(i)) return(NULL)
      data.table(part = part, model = m, term = sub("^Precision for ", "", rownames(h)[i]), sd_median = round(1 / sqrt(h[i, "0.5quant"]), 3)) })) }))
    cat("\n(b) Posterior median SD of each random effect (large interaction SD = it is absorbing variance):\n")
    print(dcast(sd_tab, part + term ~ model, value.var = "sd_median"))
    fwrite(sd_tab, file.path(output_path, "Table_S_overfit_random_effect_sd.csv"))
    # ---- (c) held-out validation
    split_defs <- list(random20 = function(d) { set.seed(5); sample.int(nrow(d), round(0.2 * nrow(d))) },
                       forward  = function(d) which(d$year >= 2022))
    ho_fit <- function(form, family, d) inla(form, family = family, data = d, control.family = if (is_int_fam(family)) INT_CF else list(),
                                             control.predictor = list(compute = TRUE, link = 1),
                                             control.compute = list(waic = FALSE, cpo = FALSE, dic = FALSE),
                                             control.inla = list(strategy = "gaussian", int.strategy = "eb"))
    NR <- nrow(dat); y_all <- dat$Y_bin; ho <- list()
    for (sp in names(split_defs)) {
      te <- split_defs[[sp]](dat); tr <- setdiff(seq_len(NR), te)
      cat(sprintf("\nSplit %s: %d training and %d test zone-months (%d test positive months)\n", sp, length(tr), length(te), sum(!is.na(dat$O_int[te]))))
      for (m in OVERFIT_MODELS) {
        t0 <- Sys.time()
        # occurrence
        d <- copy(dat); d[te, Y_bin := NA_integer_]
        f <- tryCatch(ho_fit(forms_occ[[m]], "binomial", d), error = function(e) { message("  occ ", m, " failed: ", conditionMessage(e)); NULL })
        if (!is.null(f)) {
          p <- pmin(pmax(f$summary.fitted.values$mean[seq_len(NR)], 1e-9), 1 - 1e-9)
          auc <- function(i) as.numeric(roc(y_all[i], p[i], quiet = TRUE)$auc)
          ll  <- function(i) -mean(y_all[i] * log(p[i]) + (1 - y_all[i]) * log(1 - p[i]))
          prev <- mean(y_all[tr])
          ho[[length(ho) + 1]] <- data.table(split = sp, component = "occurrence", model = m, in_sample = auc(tr), held_out = auc(te),
                                             in_sample_loss = ll(tr), held_out_loss = ll(te), held_out_Brier = mean((y_all[te] - p[te])^2),
                                             climatology_Brier = mean((y_all[te] - prev)^2), climatology_loss = -mean(y_all[te] * log(prev) + (1 - y_all[te]) * log(1 - prev)), n_test = length(te))
        }
        # intensity
        d <- copy(dat); d[te, O_int := NA_integer_]
        g <- tryCatch(ho_fit(forms_int[[m]], INT_FAM, d), error = function(e) { message("  int ", m, " failed: ", conditionMessage(e)); NULL })
        if (!is.null(g)) {
          mu <- mu_nb(g, seq_len(NR)); sz <- size_of(g)
          tp <- te[!is.na(dat$O_int[te])]; trp <- tr[!is.na(dat$O_int[tr])]
          isc <- function(i) { yi <- dat$O_int[i]; mm <- mu[i]; p0 <- (sz / (sz + mm))^sz
          c(ls = -mean(dnbinom(yi, size = sz, mu = mm, log = TRUE) - log1p(-p0)), mae = median(abs(log(yi) - log(mm / (1 - p0))))) }
          a <- isc(trp); b <- isc(tp)
          ho[[length(ho) + 1]] <- data.table(split = sp, component = "intensity", model = m, in_sample = a[["ls"]], held_out = b[["ls"]],
                                             in_sample_loss = a[["mae"]], held_out_loss = b[["mae"]], n_test = length(tp))
        }
        cat(sprintf("  %-8s done in %.1f min\n", m, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
      }
    }
    HO <- rbindlist(ho, fill = TRUE)
    fwrite(HO, file.path(output_path, "Table_S_overfit_heldout_validation.csv"))
    cat("\n--- (c) OCCURRENCE: AUC in-sample vs held-out (loss = log loss, lower is better) ---\n")
    print(HO[component == "occurrence", .(split, model, AUC_in = round(in_sample, 3), AUC_out = round(held_out, 3), AUC_drop = round(in_sample - held_out, 3),
                                          logloss_in = round(in_sample_loss, 3), logloss_out = round(held_out_loss, 3), Brier_out = round(held_out_Brier, 4),
                                          Brier_climatology = round(climatology_Brier, 4), n_test)][order(split, -AUC_out)])
    cat("\n--- (c) INTENSITY: truncated NB log score on positive months (lower is better); MAE = median abs error on the log scale ---\n")
    print(HO[component == "intensity", .(split, model, logscore_in = round(in_sample, 2), logscore_out = round(held_out, 2), change = round(held_out - in_sample, 2),
                                         MAE_log_in = round(in_sample_loss, 2), MAE_log_out = round(held_out_loss, 2), n_test)][order(split, logscore_out)])
    print(ggplot(HO[component == "occurrence"], aes(model)) +
            geom_col(aes(y = in_sample), fill = "grey75", width = 0.4, position = position_nudge(x = -0.2)) +
            geom_col(aes(y = held_out), fill = "#B2182B", width = 0.4, position = position_nudge(x = 0.2)) +
            facet_wrap(~split) + coord_cartesian(ylim = c(0.5, 1)) +
            labs(x = NULL, y = "AUC", title = "Occurrence: in-sample (grey) vs held-out (red) AUC") + theme_minimal())
    cat("\nREAD: a model that overfits shows a high in-sample AUC (grey) with a much lower held-out AUC (red), higher held-out\n",
        "log loss than the Baseline, and a large p_WAIC share of n in (a). The forward split is the stricter test: if the\n",
        "Baseline has the best or equal held-out scores there, selecting it is supported without relying on WAIC/DIC.\n")
  } else cat("RUN_OVERFIT is FALSE: Step 3B skipped.\n")
  
  
  
  ################################################################################
  # STEP 3C — PIT (CALIBRATION) FOR ALL FIVE CANDIDATE MODELS, BOTH PARTS
  #   PIT = where the observed value falls in the model's leave-one-out predictive distribution.
  #   For a well-calibrated model the PIT values are uniform (flat histogram).
  #   Counts are discrete, so the plain INLA PIT piles up on a few values; here it is RANDOMISED
  #   (a random point inside the jump of the predictive CDF), which gives a proper flat histogram when calibrated.
  #   Occurrence: leave-one-out P(outbreak) = CPO if y = 1, 1 - CPO if y = 0.
  #   Intensity : CDF jump = [PIT - CPO, PIT]. Fitted to positive months (untruncated NB): a deficit in the
  #               lowest bin is the signature of the missing zero-truncation (see Step 3, check 3).
  #   Shapes: U-shaped = predictive distribution too narrow (overconfident); hump = too wide;
  #           left-heavy = predicts too high; right-heavy = predicts too low (under-predicts large outbreaks).
  ################################################################################
  RUN_PIT <- TRUE
  if (RUN_PIT) {
    cat("\n================ STEP 3C: PIT for all candidate models ================\n")
    cpo_ok <- function(f) !is.na(f$cpo$cpo) & (is.na(f$cpo$failure) | f$cpo$failure == 0)
    pit_occ <- function(f) { ok <- cpo_ok(f); y <- dat$Y_bin[ok]; cp <- f$cpo$cpo[ok]; set.seed(11)
    runif(length(y), ifelse(y == 1, 1 - cp, 0), ifelse(y == 1, 1, cp)) }
    pit_int <- function(f) { ok <- cpo_ok(f) & !is.na(f$cpo$pit) & !is.na(dat$O_int); pt <- f$cpo$pit[ok]; cp <- f$cpo$cpo[ok]; set.seed(12)
    runif(length(pt), pmax(pt - cp, 0), pmin(pt, 1)) }
    PL <- rbindlist(c(
      lapply(names(fits_occ), function(m) data.table(part = "Occurrence", model = m, u = pit_occ(fits_occ[[m]]))),
      lapply(names(fits_int), function(m) data.table(part = "Intensity",  model = m, u = pit_int(fits_int[[m]])))))
    PL[, model := factor(model, levels = names(fits_occ))]
    sm <- PL[, .(n = .N, lowest_10pct = round(mean(u < 0.1), 3), highest_10pct = round(mean(u > 0.9), 3),
                 KS_distance_from_uniform = round(unname(ks.test(u, "punif")$statistic), 3)), by = .(part, model)]
    cat("\nShare of PIT values in the lowest / highest 10% (0.10 expected each) and KS distance from uniform (0 = perfect):\n")
    print(sm[order(part, KS_distance_from_uniform)])
    fwrite(sm, file.path(output_path, "Table_S_PIT_all_models.csv"))
    print(ggplot(PL, aes(u)) + geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 1, 0.1), fill = "grey75", colour = "white") +
            geom_hline(yintercept = 1, linetype = 2, colour = "#B2182B") + facet_grid(part ~ model, scales = "free_y") +
            labs(x = "Randomised leave-one-out PIT", y = "Density", title = "PIT: flat = calibrated (dashed line = uniform)") + theme_minimal())
    cat("\nREAD: compare each model with the dashed line. Look at the Baseline row first. Interaction models can look\n",
        "calibrated in-sample yet predict worse (Table 2, Step 3B); PIT does not replace the held-out validation.\n")
  } else cat("RUN_PIT is FALSE: Step 3C skipped.\n")
  
  
  ################################################################################
  # STEP 5 — COMPLEMENTARY DESCRIPTIVE ANALYSES (adapted from a published spatiotemporal INLA study)
  #   5a  STL decomposition of the national monthly incidence (trend, season, remainder)
  #   5b  Month-of-year effects: negative-binomial mixed model, offset log(population per 100,000),
  #       random intercept for year (which months have higher incidence, given year-to-year differences)
  #   5c  Getis-Ord Gi* on the 10-year mean incidence (queen contiguity, row-standardised) with
  #       90 / 95 / 99% confidence classes (the classical hotspot-analysis convention)
  #   5d  Model fit to the data: observed vs model-predicted zone incidence (RMSE, MAE, correlation)
  #       and uncertainty map (relative width of the 95% credible interval)
  # Not used here: environmental covariates, VIF, GAM and lags (your data have no covariates;
  # add them only if you bring climate/WASH data, and describe them as associations).
  # Wording: higher/lower-incidence clusters, not transmission hotspots. Figures: Plots pane only.
  ################################################################################
  library(ggplot2)
  n_area <- length(unique(dat$ID_area)); n_time <- max(dat$ID_time)
  
  # ---------------- 5a. STL decomposition of national monthly incidence per 100,000 ----------------
  nat_m <- dat[, .(cases = sum(O), pop = sum(monthly_population)), by = .(ID_time, year, month)][order(ID_time)]
  nat_m[, rate100k := cases / pop * 1e5]
  ts_rate <- ts(nat_m$rate100k, start = c(min(dat$year), 1), frequency = 12)
  stl_fit <- stl(ts_rate, s.window = "periodic", robust = TRUE)
  plot(stl_fit, main = "STL decomposition: national monthly incidence per 100,000")        # Plots pane
  cmp <- as.data.table(stl_fit$time.series)
  cat(sprintf("\n5a. STL: share of variance - trend %.1f%% | seasonal %.1f%% | remainder %.1f%%\n",
              100 * var(cmp$trend) / var(nat_m$rate100k), 100 * var(cmp$seasonal) / var(nat_m$rate100k),
              100 * var(cmp$remainder) / var(nat_m$rate100k)))
  seas <- data.table(month = 1:12, seasonal = tapply(cmp$seasonal, nat_m$month, mean))
  cat("Peak month of the seasonal component:", month.abb[seas$month[which.max(seas$seasonal)]],
      "| lowest:", month.abb[seas$month[which.min(seas$seasonal)]], "\n")
  fwrite(cbind(nat_m, cmp), file.path(output_path, "Table_S_STL_national_monthly.csv"))
  cat("READ: a large trend/remainder share and a small seasonal share means incidence is driven by epidemic\n",
      "waves, not a regular yearly cycle.\n")
  
  # ---------------- 5b. Month-of-year effects (negative binomial, year random intercept) ----------------
  if (requireNamespace("lme4", quietly = TRUE)) {
    tryCatch({
      nat_m[, `:=`(month_f = relevel(factor(month), ref = "2"), year_f = factor(year), off = log(pop / 1e5))]
      mm <- lme4::glmer.nb(cases ~ month_f + offset(off) + (1 | year_f), data = nat_m)
      co <- summary(mm)$coefficients; est <- co[grep("^month_f", rownames(co)), , drop = FALSE]
      irr <- data.table(month = as.integer(sub("month_f", "", rownames(est))),
                        IRR = exp(est[, 1]), lo = exp(est[, 1] - 1.96 * est[, 2]), hi = exp(est[, 1] + 1.96 * est[, 2]),
                        p = est[, 4])
      irr <- rbind(data.table(month = 2L, IRR = 1, lo = NA_real_, hi = NA_real_, p = NA_real_), irr)[order(month)]
      irr[, month_name := month.abb[month]]
      cat("\n5b. Month effects, incidence rate ratio vs February (95% Wald CI):\n"); print(irr[, .(month_name, IRR = round(IRR, 2), lo = round(lo, 2), hi = round(hi, 2), p = signif(p, 3))])
      fwrite(irr, file.path(output_path, "Table_S_month_effects_NB_mixed.csv"))
      print(ggplot(irr, aes(factor(month_name, levels = month.abb), IRR)) +
              geom_hline(yintercept = 1, linetype = 2, colour = "grey50") + geom_pointrange(aes(ymin = lo, ymax = hi)) +
              scale_y_log10() + labs(x = NULL, y = "Incidence rate ratio vs February (log scale)",
                                     title = "Month-of-year effects (NB mixed model, random year intercept)") + theme_minimal())
      cat("READ: months with the interval above 1 have higher incidence than February in the same year, on average.\n",
          "Only 120 national monthly values and a few large waves: treat as descriptive. Year-level waves are not seasonality.\n")
    }, error = function(e) cat("5b could not run:", conditionMessage(e), "\n"))
  } else cat("\n5b skipped: install.packages('lme4') to run the month-effects model.\n")
  
  # ---------------- 5c. Gi* with 90/95/99% classes on the 10-year mean incidence ----------------
  nbS  <- mat2listw(as.matrix(Q), style = "B")$neighbours
  lwS  <- nb2listw(include.self(nbS), style = "W", zero.policy = TRUE)          # row-standardised, self included (Gi*)
  crude5 <- dat[, .(cases = sum(O), pm = sum(monthly_population)), by = ID_area][order(ID_area)]
  inc10 <- crude5$cases / crude5$pm * 12 * 1e5                                  # mean annual incidence per 100,000, 10 years
  gz5 <- as.numeric(localG(inc10, lwS, zero.policy = TRUE))
  cls5 <- cut(gz5, breaks = c(-Inf, -2.576, -1.960, -1.645, 1.645, 1.960, 2.576, Inf),
              labels = c("Cold 99%", "Cold 95%", "Cold 90%", "Not significant", "Hot 90%", "Hot 95%", "Hot 99%"))
  cat("\n5c. Gi* classes on the 10-year mean incidence (no multiple-testing correction; see Step 2 for FDR/permutation):\n")
  print(table(cls5))
  shp5 <- st_read(shp_path, quiet = TRUE); shp5$zone <- toupper(gsub("[^A-Z0-9]", "", shp5$ADM2_PCODE))
  shp5 <- shp5[match(zone_order, shp5$zone), ]; stopifnot(!anyNA(shp5$zone))
  shp5$gi_class <- cls5; shp5$inc10 <- inc10
  print(ggplot(shp5) + geom_sf(aes(fill = gi_class), colour = "white", linewidth = 0.1) +
          scale_fill_manual(values = c("Cold 99%" = "#08519C", "Cold 95%" = "#4292C6", "Cold 90%" = "#9ECAE1", "Not significant" = "grey88",
                                       "Hot 90%" = "#FCAE91", "Hot 95%" = "#EF3B2C", "Hot 99%" = "#99000D"), name = "Gi* class", drop = FALSE) +
          labs(title = "Getis-Ord Gi*: 10-year mean incidence per 100,000 per year") + theme_void())
  fwrite(data.table(ID_area = 1:n_area, zone = zone_order, incidence_100k_year = inc10, Gi_star_z = gz5, class = as.character(cls5)),
         file.path(output_path, "Table_S_Gi_star_10yr_mean_incidence.csv"))
  cat("READ: hot/cold = clusters of zones with higher/lower observed incidence than the national pattern. With 92 zones\n",
      "and no correction, 90% classes are exploratory; rely on 95/99% and on the FDR version in Step 2.\n")
  
  # ---------------- 5d. Observed vs predicted, and uncertainty map ----------------
  zz <- merge(data.table(ID_area = 1:n_area, observed = inc10), zone[, .(ID_area, predicted = burden_100k, lo = burden_lo, hi = burden_hi)], by = "ID_area")
  zz[, rel_width := (hi - lo) / pmax(predicted, 1e-6)]
  fit_stats <- zz[, .(RMSE = sqrt(mean((observed - predicted)^2)), MAE = mean(abs(observed - predicted)),
                      spearman = cor(observed, predicted, method = "spearman"),
                      pearson_log = cor(log1p(observed), log1p(predicted)),
                      inside_95CrI = mean(observed >= lo & observed <= hi))]
  cat("\n5d. Zone-level observed vs model-predicted annual incidence per 100,000 (IN-SAMPLE; the out-of-time check is in Step 3):\n")
  print(round(fit_stats, 3))
  cat(sprintf("Median relative width of the 95%% CrI: %.2f | zones with width > 2x the estimate: %d\n",
              median(zz$rel_width), sum(zz$rel_width > 2)))
  print(ggplot(zz, aes(observed + 0.1, predicted + 0.1)) + geom_abline(linetype = 2, colour = "grey50") +
          geom_point(alpha = 0.7) + scale_x_log10() + scale_y_log10() +
          labs(x = "Observed incidence per 100,000 per year (+0.1)", y = "Model-predicted expected caseload per 100,000 per year (+0.1)",
               title = "Observed vs predicted, 92 zones") + theme_minimal())
  shp5$rel_width <- zz$rel_width[match(1:n_area, zz$ID_area)]
  print(ggplot(shp5) + geom_sf(aes(fill = rel_width), colour = "white", linewidth = 0.1) +
          scale_fill_viridis_c(option = "magma", direction = -1, trans = "sqrt", name = "Relative CrI width") +
          labs(title = "Uncertainty of the expected caseload: (upper - lower) / estimate") + theme_void())
  fwrite(zz, file.path(output_path, "Table_S_observed_vs_predicted_uncertainty.csv"))
  cat("READ: zones with no observed cases have a smoothed, non-zero prediction by design (borrowing from neighbours),\n",
      "so observed and predicted differ there. A wide interval means the estimate is uncertain, not that risk is high.\n")
  if (SAVE_FIGS) cat("Note: Step 5 figures are shown in the Plots pane; use the Export button there, or add ggsave() calls.\n")
  
  
  
  # (The prior grid search now runs INSIDE Step 1, before the final models are fitted: see "6A. PRIOR GRID".)
  
  ################################################################################
  # STEP 7 — TWO MODEL-STRUCTURE CHECKS (adapted from a zero-inflated space-time study)
  #   7A  Hurdle versus zero-inflated models: is the hurdle a better description of the zeros?
  #       Same Baseline structure; the zero part of INLA's zero-inflated families has a CONSTANT
  #       probability, while the hurdle models it with BYM2 + time effects. Compared by log-CPO and WAIC.
  #   7B  Population adjustment and spatial confounding: the zone effect (BYM2) and population are both
  #       spatial, so adding or removing log(population) can move the zone effects. The Baseline is refitted
  #       with the opposite population choice and the zone effects are compared (this is also the
  #       OCC_POP sensitivity for Reviewer comment 2).
  #   Not added: the two-stage covariate / cluster-specific-coefficient / NIMBLE-MCMC machinery. It serves
  #   covariate effects (temperature, etc.); your paper estimates no covariate effects.
  ################################################################################
  RUN_STRUCT <- TRUE
  if (RUN_STRUCT) {
    cat("\n================ STEP 7: model-structure checks ================\n")
    mk_formula <- function(lhs, U_rw, extra, int = FALSE) {     # same priors as the main model of that part (selected by the grid)
      pri_bym7 <- if (int) prior_bym2_int else prior_bym2
      pri_iid7 <- if (int) prior_iid_int else prior_iid
      pri_rw7 <- list(theta = list(prior = "pc.prec", param = c(U_rw, 0.01)))
      as.formula(paste(lhs, paste(SEAS, collapse = " + "),
                       "+ f(ID_area_bym, model='bym2', graph=Q, scale.model=TRUE, constr=TRUE, hyper=pri_bym7)",
                       "+ f(ID_time_rw, model='rw2', scale.model=TRUE, constr=TRUE, hyper=pri_rw7)",
                       "+ f(ID_time_iid, model='iid', hyper=pri_iid7)", extra), env = environment())
    }
    run7 <- function(f, family, ...) inla(f, family = family, data = dat, control.predictor = list(compute = TRUE, link = 1),
                                          control.compute = list(dic = FALSE, waic = TRUE, cpo = TRUE),
                                          control.inla = list(strategy = "simplified.laplace"), ...)
    cpo_mean <- function(f, idx = NULL) { cp <- f$cpo$cpo; fl <- f$cpo$failure; if (!is.null(idx)) { cp <- cp[idx]; fl <- fl[idx] }
    ok <- !is.na(cp) & (is.na(fl) | fl == 0); mean(log(cp[ok])) }
    phi7 <- function(f) f$summary.hyperpar[grep("^Phi", rownames(f$summary.hyperpar))[1], "mean"]
    
    # ---------------- 7A. hurdle versus zero-inflated ----------------
    NN7 <- nrow(dat)
    lp_o <- sum(log(fit_o$cpo$cpo[1:NN7]), na.rm = TRUE)                       # hurdle: occurrence part, all zone-months
    pos7 <- which(!is.na(dat$O_int))
    lp_i <- sum(log(fit_i$cpo$cpo[pos7]), na.rm = TRUE)                        # hurdle: intensity part, positive months only
    hurdle <- data.table(model = "Hurdle (Baseline occurrence + Baseline intensity)", family = if (INT_ZT) "bernoulli + zero-truncated NB" else "bernoulli + NB",
                         WAIC = fit_o$waic$waic + fit_i$waic$waic, sum_logCPO = lp_o + lp_i, mean_logCPO_per_zone_month = (lp_o + lp_i) / NN7)
    cat(sprintf("\n7A. Hurdle: WAIC %.1f | mean log-CPO per zone-month %.4f\n", hurdle$WAIC, hurdle$mean_logCPO_per_zone_month))
    cmpf <- list(c("NB, no zero-inflation", "nbinomial"), c("Zero-inflated Poisson (constant p)", "zeroinflatedpoisson1"),
                 c("Zero-inflated NB (constant p)", "zeroinflatednbinomial1"))
    zi_rows <- lapply(cmpf, function(a) tryCatch({
      cat("  fitting", a[1], "...\n")
      f <- run7(mk_formula("O ~ offset(logE) +", U_RW_INT_EFF, int_extra, int = TRUE), a[2])
      zp <- f$summary.hyperpar[grep("zero|Zero", rownames(f$summary.hyperpar)), "mean"]
      data.table(model = a[1], family = a[2], WAIC = f$waic$waic, sum_logCPO = sum(log(f$cpo$cpo), na.rm = TRUE),
                 mean_logCPO_per_zone_month = cpo_mean(f), CPO_failures = sum(f$cpo$failure > 0, na.rm = TRUE),
                 zero_inflation_p = if (length(zp)) round(zp[1], 3) else NA_real_)
    }, error = function(e) { cat("  ", a[1], "failed:", conditionMessage(e), "\n"); NULL }))
    cmp7a <- rbindlist(c(list(hurdle), zi_rows), fill = TRUE)
    cmp7a[, `:=`(WAIC = round(WAIC, 1), sum_logCPO = round(sum_logCPO, 1), mean_logCPO_per_zone_month = round(mean_logCPO_per_zone_month, 4))]
    cat("\n--- 7A. Hurdle vs zero-inflated (same data, Baseline structure) ---\n"); print(cmp7a)
    fwrite(cmp7a, file.path(output_path, "Table_S_hurdle_vs_zeroinflated.csv"))
    cat("READ: lower WAIC and a higher (less negative) sum_logCPO = better. Do NOT compare mean_logCPO_per_zone_month across rows: for models with CPO failures it averages over the non-failed observations only, so it is not on the same footing as the hurdle. The hurdle's intensity part is the main intensity likelihood (zero-truncated NB when INT_ZT is TRUE).\n",
        " If a zero-inflated model were clearly better,\n",
        "the zeros would look like a constant 'extra zeros' process; a hurdle win supports modelling occurrence as a\n",
        "structured space-time process. Note: the zero-inflated probability here is a single national value.\n")
    
    # ---------------- 7B. population adjustment / spatial confounding ----------------
    ze <- function(f) f$summary.random$ID_area_bym$mean[1:n_area]
    alt_occ_extra <- if (OCC_POP) "" else "+ log_pop_s"
    alt_int_extra <- if (INT_POP) "" else "+ log_pop_zone_s"
    cat("\n  refitting occurrence with", if (OCC_POP) "NO population term" else "log(population)", "...\n")
    f_oalt <- run7(mk_formula("Y_bin ~", U_RW_OCC_EFF, alt_occ_extra), "binomial")
    cat("  refitting intensity with", if (INT_POP) "NO population term" else "log(zone population)", "...\n")
    f_ialt <- run7(mk_formula("O_int ~ offset(logE) +", U_RW_INT_EFF, alt_int_extra, int = TRUE), INT_FAM, control.family = INT_CF)
    cmp7b <- rbindlist(lapply(list(
      list("occurrence", fit_o, f_oalt, if (OCC_POP) "with log(pop)" else "without log(pop)", if (OCC_POP) "without log(pop)" else "with log(pop)"),
      list("intensity",  fit_i, f_ialt, if (INT_POP) "with log(pop)" else "without log(pop)", if (INT_POP) "without log(pop)" else "with log(pop)")),
      function(a) data.table(part = a[[1]], main = a[[4]], alternative = a[[5]],
                             zone_effect_spearman = round(cor(ze(a[[2]]), ze(a[[3]]), method = "spearman"), 3),
                             top10_shared = length(intersect(order(-ze(a[[2]]))[1:10], order(-ze(a[[3]]))[1:10])),
                             phi_main = round(phi7(a[[2]]), 3), phi_alt = round(phi7(a[[3]]), 3),
                             logCPO_main = round(cpo_mean(a[[2]]), 4), logCPO_alt = round(cpo_mean(a[[3]]), 4),
                             WAIC_main = round(a[[2]]$waic$waic, 1), WAIC_alt = round(a[[3]]$waic$waic, 1))))
    cat("\n--- 7B. Zone effects with and without the population term ---\n"); print(cmp7b)
    popco <- function(f, nm) { r <- f$summary.fixed[nm, c("mean", "0.025quant", "0.975quant")]; round(unlist(r), 3) }
    if ("log_pop_s" %in% rownames(f_oalt$summary.fixed) || "log_pop_s" %in% rownames(fit_o$summary.fixed)) {
      fo_pop <- if (OCC_POP) fit_o else f_oalt
      cat("Occurrence: coefficient of scaled log(population) [mean, 95% CrI]:", popco(fo_pop, "log_pop_s"), "\n")
    }
    if ("log_pop_zone_s" %in% rownames(f_ialt$summary.fixed) || "log_pop_zone_s" %in% rownames(fit_i$summary.fixed)) {
      fi_pop <- if (INT_POP) fit_i else f_ialt
      cat("Intensity: coefficient of scaled log(zone population) [mean, 95% CrI]:", popco(fi_pop, "log_pop_zone_s"), "\n")
    }
    fwrite(cmp7b, file.path(output_path, "Table_S_population_adjustment_spatial_confounding.csv"))
    zz7 <- data.table(occ_main = ze(fit_o), occ_alt = ze(f_oalt), int_main = ze(fit_i), int_alt = ze(f_ialt))
    print(ggplot(rbind(zz7[, .(part = "Occurrence zone effect", main = occ_main, alt = occ_alt)],
                       zz7[, .(part = "Intensity zone effect", main = int_main, alt = int_alt)]), aes(main, alt)) +
            geom_abline(linetype = 2, colour = "grey50") + geom_point(alpha = 0.7) + facet_wrap(~part, scales = "free") +
            labs(x = "Main model", y = "With the opposite population choice", title = "Zone effects: main model vs population-choice alternative") +
            theme_minimal())
    cat("READ: Spearman near 1 and 8-10 shared top-10 zones = the zone results do not depend on the population choice\n",
        "(no strong spatial confounding between population and the zone effect). A clear drop means part of the zone\n",
        "effect reflects population size: report the main model and this contrast, and word the zone effect as\n",
        "'after adjustment for population', not as an absolute property of the zone.\n")
    
    # ---------------- 7C. intensity likelihood: zero-truncated NB versus ordinary NB ----------------
    # The intensity part is fitted to positive months only. MAIN (INT_ZT = TRUE): the zero-truncated NB (INLA 'zeroinflatednbinomial0' with the zero
    # probability fixed near 0), the correct likelihood for counts that cannot be 0. The ordinary NB (the earlier main model) is refitted here as the
    # sensitivity analysis; with INT_ZT = FALSE the two roles are swapped. Needs an INLA version that has the family (names(inla.models()$likelihood)).
    RUN_TRUNC <- TRUE
    if (RUN_TRUNC && "zeroinflatednbinomial0" %in% names(INLA::inla.models()$likelihood)) {
      alt_fam <- if (INT_ZT) "nbinomial" else "zeroinflatednbinomial0"
      sz_hyper <- if (INT_ZT) INT_CF$hyper$theta1 else INT_CF$hyper$theta               # keep the size prior of the main fit
      cf_alt <- mk_int_cf(sz_hyper, zt = !INT_ZT)
      cat("\n  7C. refitting the intensity part with the", if (INT_ZT) "ordinary NB" else "zero-truncated NB", "(sensitivity analysis)...\n")
      f_alt <- tryCatch(run7(mk_formula("O_int ~ offset(logE) +", U_RW_INT_EFF, if (INT_POP) "+ log_pop_zone_s" else "", int = TRUE),
                             alt_fam, control.family = cf_alt), error = function(e) { message("  7C failed: ", e$message); NULL })
      if (!is.null(f_alt)) {
        f_tr <- if (INT_ZT) fit_i else f_alt; f_nb <- if (INT_ZT) f_alt else fit_i          # f_tr = zero-truncated NB, f_nb = ordinary NB
        pos <- which(!is.na(dat$O_int)); allr <- seq_len(nrow(dat))
        sz_nb <- size_of(f_nb); sz_tr <- size_of(f_tr)
        mu_trv <- mu_nb(f_tr, pos); p0_tr <- (sz_tr / (sz_tr + mu_trv))^sz_tr
        exp_tr  <- function(f) if ("log_pop_zone_s" %in% rownames(f$summary.fixed))
          round(1 + unlist(f$summary.fixed["log_pop_zone_s", c("mean", "0.025quant", "0.975quant")]) / POP_SD, 2) else rep(NA_real_, 3)
        ord <- if (INT_ZT) list(f_tr, f_nb) else list(f_nb, f_tr)                           # main model first
        lab <- if (INT_ZT) c("zero-truncated NB (main)", "ordinary NB") else c("ordinary NB (main)", "zero-truncated NB")
        szs <- c(size_of(ord[[1]]), size_of(ord[[2]]))
        cmp7c <- data.table(
          model = lab,
          size = round(szs, 3), alpha = round(1 / szs, 2),
          pop_exponent = c(exp_tr(ord[[1]])[1], exp_tr(ord[[2]])[1]), pop_lo = c(exp_tr(ord[[1]])[2], exp_tr(ord[[2]])[2]), pop_hi = c(exp_tr(ord[[1]])[3], exp_tr(ord[[2]])[3]),
          phi = round(c(phi7(ord[[1]]), phi7(ord[[2]])), 3), WAIC = round(c(ord[[1]]$waic$waic, ord[[2]]$waic$waic), 1),
          mean_logCPO = round(c(cpo_mean(ord[[1]], pos), cpo_mean(ord[[2]], pos)), 4))
        cat("\n--- 7C. Intensity: zero-truncated NB vs ordinary NB (positive months only) ---\n"); print(cmp7c)
        sp_ze <- cor(ze(f_tr), ze(f_nb), method = "spearman"); t10_ze <- length(intersect(order(-ze(f_tr))[1:10], order(-ze(f_nb))[1:10]))
        cat(sprintf("Zone effects: Spearman = %.3f | top-10 shared = %d of 10\n", sp_ze, t10_ze))
        zsize <- function(v) tapply(v, dat$ID_area, mean)                                   # zone mean of the expected size per outbreak month
        zs_tr <- zsize(zt_mean(mu_nb(f_tr, allr), sz_tr)); zs_nb <- zsize(mu_nb(f_nb, allr))
        sp_sz <- cor(zs_tr, zs_nb, method = "spearman"); t10_sz <- length(intersect(order(-zs_tr)[1:10], order(-zs_nb)[1:10]))
        cat(sprintf("Zone-level expected size per outbreak month E(Y | Y > 0): Spearman = %.3f | top-10 shared = %d of 10\n", sp_sz, t10_sz))
        cat(sprintf("Under the truncated fit: median P(0) of the underlying NB at positive months = %.3f | share of months with P(0) > 5%% = %.0f%%\n",
                    median(p0_tr), 100 * mean(p0_tr > 0.05)))
        cat(sprintf("Mean expected size at positive months: zero-truncated NB %.1f | ordinary NB %.1f | observed mean %.1f\n",
                    mean(zt_mean(mu_trv, sz_tr)), mean(mu_nb(f_nb, pos)), mean(dat$O_int[pos])))
        fwrite(cmp7c, file.path(output_path, "Table_S_zero_truncated_intensity.csv"))
        fwrite(data.table(comparison = c("zone effects (BYM2)", "zone mean expected size per outbreak month"),
                          spearman = round(c(sp_ze, sp_sz), 3), top10_shared = c(t10_ze, t10_sz)),
               file.path(output_path, "Table_S_zero_truncated_ranking.csv"))
        cat("READ: log-CPO and WAIC of the two likelihoods are comparable on the same positive months (both are likelihoods of the observed counts); the truncated\n",
            "NB does not waste probability on 0. Spearman near 1 and a similar population exponent mean that the choice of likelihood did not change\n",
            "the zone results; alpha differs because the two families share the variance differently. Report the main model and this contrast.\n")
      }
    } else cat("7C skipped: RUN_TRUNC is FALSE or this INLA has no 'zeroinflatednbinomial0' family.\n")
  } else cat("RUN_STRUCT is FALSE: Step 7 skipped.\n")
  
  
  ################################################################################
  # STEP 8 — PERSISTENCE EXTENSION (the predictive model) AND ITS CHECKS
  #   The Baseline has no memory: a zone's risk this month does not depend on last month. The data say it should
  #   (occurrence odds about 35-fold after an outbreak month; larger counts follow larger counts). This step fits
  #   Baseline + persistence terms and compares it with the Baseline on:
  #     8a  forward one-step-ahead validation: train 2015-2021, predict each month of 2022-2024 knowing the PREVIOUS month
  #     8b  residual autocorrelation (the Step 3 check 5 that the Baseline failed)
  #     8c  what happens to the zone effects and phi (persistence takes over part of the spatial structure)
  #   Wording: persistence = association between consecutive months (continuing outbreaks, continuing surveillance,
  #   shared conditions), not transmission. The Baseline stays the DESCRIPTIVE model; this is the PREDICTIVE extension.
  ################################################################################
  RUN_PERSIST_MODEL <- TRUE
  if (RUN_PERSIST_MODEL) {
    cat("\n================ STEP 8: persistence extension ================\n")
    dP <- copy(dat); setorder(dP, ID_area, ID_time)
    dP[, lag_own := as.numeric(shift(Y_bin, 1L, fill = 0L)), by = ID_area]
    Ym <- matrix(dP$Y_bin, nrow = n_time, ncol = n_area); Yb <- rbind(0, Ym[-n_time, , drop = FALSE])
    dP[, lag_nb := as.vector(as.matrix(Yb %*% t(as.matrix(Q))) / rep(rowSums(as.matrix(Q)), each = n_time))]
    dP[, lag_nb := as.numeric(scale(lag_nb))]
    dP[, prev_cases := shift(O, 1L, fill = 0L), by = ID_area]
    dP[, `:=`(prev_occ = as.numeric(prev_cases > 0), prev_size = as.numeric(scale(log1p(prev_cases))))]
    setorder(dP, ID_time, ID_area); stopifnot(identical(dP$ID_area, dat$ID_area), identical(dP$ID_time, dat$ID_time))
    f_occ_P <- update(forms_occ[["Baseline"]], . ~ . + lag_own + lag_nb)
    f_int_P <- update(forms_int[["Baseline"]], . ~ . + prev_occ + prev_size)
    full_fit <- function(form, family, d) inla(form, family = family, data = d, control.family = if (is_int_fam(family)) INT_CF else list(),
                                               control.predictor = list(compute = TRUE, link = 1),
                                               control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                                               control.inla = list(strategy = "simplified.laplace"))
    cat("Fitting the persistence models on all data (a few minutes)...\n")
    fP_o <- full_fit(f_occ_P, "binomial", dP); fP_i <- full_fit(f_int_P, INT_FAM, dP)
    
    # ---- 8a. forward one-step-ahead validation
    te <- which(dat$year >= 2022); tr <- setdiff(seq_len(nrow(dat)), te); NR <- nrow(dat); yb <- dat$Y_bin
    cat(sprintf("\n8a. Forward validation: %d training zone-months (2015-2021), %d test (2022-2024)\n", length(tr), length(te)))
    fw_occ <- function(form, d0) { d <- copy(d0); d[te, Y_bin := NA_integer_]; f <- ho_fit(form, "binomial", d)
    p <- pmin(pmax(f$summary.fitted.values$mean[seq_len(NR)], 1e-9), 1 - 1e-9); prev <- mean(yb[tr])
    data.table(AUC = as.numeric(roc(yb[te], p[te], quiet = TRUE)$auc), logloss = -mean(yb[te] * log(p[te]) + (1 - yb[te]) * log(1 - p[te])),
               Brier = mean((yb[te] - p[te])^2), Brier_climatology = mean((yb[te] - prev)^2)) }
    fw_int <- function(form, d0) { d <- copy(d0); d[te, O_int := NA_integer_]; g <- ho_fit(form, INT_FAM, d)
    mu <- mu_nb(g, seq_len(NR)); sz <- size_of(g); tp <- te[!is.na(dat$O_int[te])]
    yi <- dat$O_int[tp]; mm <- mu[tp]; p0 <- (sz / (sz + mm))^sz
    data.table(logscore_truncNB = -mean(dnbinom(yi, size = sz, mu = mm, log = TRUE) - log1p(-p0)),
               MAE_log = median(abs(log(yi) - log(mm / (1 - p0)))), n_test_positive = length(tp)) }
    VA <- rbind(cbind(part = "occurrence", model = "Baseline",            fw_occ(forms_occ[["Baseline"]], dP)),
                cbind(part = "occurrence", model = "Baseline + persistence", fw_occ(f_occ_P, dP)), fill = TRUE)
    VI <- rbind(cbind(part = "intensity", model = "Baseline",              fw_int(forms_int[["Baseline"]], dP)),
                cbind(part = "intensity", model = "Baseline + persistence", fw_int(f_int_P, dP)), fill = TRUE)
    cat("\nOccurrence (higher AUC and lower log loss / Brier are better):\n"); print(VA[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)])
    cat("\nIntensity (lower is better):\n"); print(VI[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
    fwrite(rbind(VA, VI, fill = TRUE), file.path(output_path, "Table_S_persistence_forward_validation.csv"))
    
    # ---- 8b. residual autocorrelation
    res_occ <- function(f) { p <- pmin(pmax(f$summary.fitted.values$mean[seq_len(NR)], 1e-9), 1 - 1e-9); (yb - p) / sqrt(p * (1 - p)) }
    lag1_within <- function(r_all) { d <- data.table(a = dat$ID_area, t = dat$ID_time, r = r_all)[order(a, t)]; d[, rl := shift(r), by = a]; cor(d$r, d$rl, use = "complete.obs") }
    res_int <- function(f) ifelse(is.na(dat$O_int), NA_real_, pearson_int(f, dat$O_int, seq_len(NR)))
    lag1_int <- function(f) { r <- res_int(f); d <- data.table(a = dat$ID_area, t = dat$ID_time, r = r)[order(a, t)]
    d[, `:=`(rl = shift(r), tl = shift(t)), by = a]; d <- d[!is.na(r) & !is.na(rl) & t - tl == 1]; c(r = cor(d$r, d$rl), pairs = nrow(d)) }
    RS <- data.table(part = c("occurrence", "occurrence", "intensity", "intensity"), model = rep(c("Baseline", "Baseline + persistence"), 2),
                     lag1_residual_r = round(c(lag1_within(res_occ(fit_o)), lag1_within(res_occ(fP_o)), lag1_int(fit_i)[["r"]], lag1_int(fP_i)[["r"]]), 3),
                     target = "|r| < 0.10")
    cat("\n8b. Lag-1 residual correlation within zones (consecutive months):\n"); print(RS)
    fwrite(RS, file.path(output_path, "Table_S_persistence_residuals.csv"))
    
    # ---- 8c. effect on zone results
    ze8 <- function(f) f$summary.random$ID_area_bym$mean[1:n_area]
    cmp8 <- data.table(part = c("occurrence", "intensity"),
                       zone_effect_spearman = round(c(cor(ze8(fit_o), ze8(fP_o), method = "spearman"), cor(ze8(fit_i), ze8(fP_i), method = "spearman")), 3),
                       top10_shared = c(length(intersect(order(-ze8(fit_o))[1:10], order(-ze8(fP_o))[1:10])), length(intersect(order(-ze8(fit_i))[1:10], order(-ze8(fP_i))[1:10]))),
                       phi_baseline = round(c(phi7(fit_o), phi7(fit_i)), 3), phi_persistence = round(c(phi7(fP_o), phi7(fP_i)), 3),
                       WAIC_baseline = round(c(fit_o$waic$waic, fit_i$waic$waic), 1), WAIC_persistence = round(c(fP_o$waic$waic, fP_i$waic$waic), 1))
    cat("\n8c. Zone effects and phi, Baseline vs persistence:\n"); print(cmp8)
    cof <- rbind(cbind(part = "occurrence", term = c("lag_own", "lag_nb"), round(fP_o$summary.fixed[c("lag_own", "lag_nb"), c("mean", "0.025quant", "0.975quant")], 3)),
                 cbind(part = "intensity",  term = c("prev_occ", "prev_size"), round(fP_i$summary.fixed[c("prev_occ", "prev_size"), c("mean", "0.025quant", "0.975quant")], 3)))
    cat("\nPersistence coefficients (log scale):\n"); print(cof)
    fwrite(cmp8, file.path(output_path, "Table_S_persistence_zone_effects.csv"))
    cat("\nREAD: if forward AUC rises above 0.75, Brier beats climatology and lag-1 residual r falls below 0.10, report the extension as the\n",
        "predictive model. Spearman near 1 and a similar top-10 mean the zone ranking does not depend on persistence; a lower phi means part\n",
        "of the spatial structure was persistence. Keep the Baseline for description. Persistence = association between consecutive months.\n")
  } else cat("RUN_PERSIST_MODEL is FALSE: Step 8 skipped.\n")
  
  
  ################################################################################
  # STEP 9 — MAPS WITH HONEST HANDLING OF ZONES WITHOUT DATA
  #   Intensity is estimated from positive months only. The zones with no positive month have no intensity data (their value
  #   comes only from neighbours), so they are drawn GREY; zones with 1-2 positive months are drawn but outlined in black.
  #   Occurrence and expected caseload are shown for all 92 zones (a zero month is information). Plots pane only.
  ################################################################################
  RUN_MAPS <- FALSE      # quick-look maps; SUPERSEDED by the publication maps in STEP 11 (set TRUE only if you want these too)
  if (RUN_MAPS) {
    cat("\n================ STEP 9: maps ================\n")
    shp9 <- st_read(shp_path, quiet = TRUE); shp9$zone <- toupper(gsub("[^A-Z0-9]", "", shp9$ADM2_PCODE))
    shp9 <- shp9[match(zone_order, shp9$zone), ]; stopifnot(!anyNA(shp9$zone))
    z9 <- zone[order(ID_area)]; stopifnot(nrow(z9) == nrow(shp9))
    shp9$pi <- z9$pi; shp9$mu100k <- z9$mu_100k; shp9$burden <- z9$burden_100k; shp9$npos <- z9$N_outbreak_months
    shp9$mu100k_plot <- ifelse(shp9$npos == 0, NA_real_, shp9$mu100k)
    shp9$few <- shp9$npos %in% 1:2
    shp9$mu_width <- ifelse(shp9$npos == 0, NA_real_, (z9$mu_100k_hi - z9$mu_100k_lo) / pmax(z9$mu_100k, 1e-6))
    cat(sprintf("Zones with no positive month (grey on the intensity map): %d | with 1-2 positive months (outlined): %d\n", sum(shp9$npos == 0), sum(shp9$few)))
    mp <- function(var, ttl, nm, opt = "viridis", tr = "identity", dir = 1) {
      g <- ggplot(shp9) + geom_sf(aes(fill = .data[[var]]), colour = "white", linewidth = 0.1) +
        scale_fill_viridis_c(option = opt, trans = tr, direction = dir, name = nm, na.value = "grey80") + labs(title = ttl) + theme_void()
      if (any(shp9$few) && var %in% c("mu100k_plot", "mu_width")) g <- g + geom_sf(data = shp9[shp9$few, ], fill = NA, colour = "black", linewidth = 0.6)
      g }
    print(mp("pi", "Outbreak probability per zone-month (all 92 zones)", "Probability", "viridis", "sqrt"))
    print(mp("mu100k_plot", "Cases per 100,000 in an outbreak month (grey = no positive month; black outline = 1-2 positive months)", "Per 100,000", "magma", "sqrt", -1))
    print(mp("burden", "Expected annual caseload per 100,000 (probability x size)", "Per 100,000", "inferno", "sqrt", -1))
    print(mp("mu_width", "Uncertainty of intensity: relative 95% CrI width (grey = no positive month)", "Relative width", "mako", "sqrt", -1))
    if (SAVE_FIGS) cat("SAVE_FIGS is TRUE: add ggsave() calls or use the Plots-pane Export.\n")
    cat("READ: grey zones have no intensity estimate of their own. The probability and caseload maps are estimates for all zones; for the\n",
        "never-reporting zones they are low and rely on neighbours and the national pattern.\n")
  } else cat("RUN_MAPS is FALSE: Step 9 skipped.\n")
  
  
  ################################################################################
  # STEP 10 — MANUSCRIPT TABLES WITH gt (docx).  Needs install.packages("gt"); the code is also available as Stage_tables_gt.R
  ################################################################################
  RUN_TABLES_GT <- TRUE
  SAVE_TABLES   <- TRUE
  if (RUN_TABLES_GT && requireNamespace("gt", quietly = TRUE)) {
    # ==============================================================================
    # STAGE: MANUSCRIPT TABLES WITH gt  (each table is shown in the RStudio Viewer and saved as .docx)
    # Reads the CSV files written by the pipeline (so it can be run on its own, after the pipeline).
    # Needs: install.packages(c("gt", "data.table")). Word export needs a recent gt (>= 0.10): update.packages("gt") if gtsave fails.
    # Wording: associations only; phi = structured share of the spatial variance; no transmission / hub / hotspot language.
    # ==============================================================================
    suppressPackageStartupMessages({ library(data.table); library(gt) })
    if (!exists("output_path")) output_path <- getOption("cholera.proj", ".")
    if (!exists("SUFFIX"))      SUFFIX      <- "_OCCPOP"          # population in both parts is the main run
    if (!exists("SAVE_TABLES")) SAVE_TABLES <- TRUE               # TRUE = also save each table as .docx (gtsave)
    TAB_DIR <- file.path(output_path, "tables_docx"); if (SAVE_TABLES) dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
    if (utils::packageVersion("gt") < "0.10.0") message("NOTE: gt ", utils::packageVersion("gt"), " may not export .docx; update.packages('gt') if saving fails.")
    
    rd <- function(f) { fp <- file.path(output_path, f)
    if (!file.exists(fp)) { message("  (missing: ", f, ")"); return(NULL) }
    if (exists("RUN_START") && file.info(fp)$mtime < RUN_START) { message("  (SKIPPED, older than this run: ", f, " - written ", format(file.info(fp)$mtime), ")"); return(NULL) }
    fread(fp) }
    safe <- function(label, expr) tryCatch(expr, error = function(e) message("  TABLE FAILED [", label, "]: ", conditionMessage(e)))
    finish <- function(tb, name) {
      print(tb)                                                        # Viewer
      if (SAVE_TABLES) {
        fp <- file.path(TAB_DIR, paste0(name, ".docx"))
        ok <- tryCatch({ gtsave(tb, fp); TRUE }, error = function(e) { message("  docx failed for ", name, ": ", conditionMessage(e)); FALSE })
        if (ok) cat("  saved:", fp, "\n") else tryCatch(gtsave(tb, file.path(TAB_DIR, paste0(name, ".html"))), error = function(e) NULL)
      }
      invisible(tb)
    }
    pub <- function(tb, title, subtitle = NULL, notes = NULL) {
      tb <- tab_header(tb, title = title, subtitle = subtitle)
      for (n in notes) tb <- tab_source_note(tb, source_note = n)
      tryCatch(tab_options(tb, table.font.size = px(11), table.font.names = "Arial", heading.align = "left"), error = function(e) tb)
    }
    ci <- function(m, lo, hi, d = 2, k = 1) sprintf(paste0("%.", d, "f (%.", d, "f to %.", d, "f)"), m * k, lo * k, hi * k)
    model_lab <- c(Baseline = "Baseline (no interaction)", Type_I = "Type I (iid space x iid time)", Type_II = "Type II (iid space x RW2 time)",
                   Type_III = "Type III (BYM2 space x iid time)", Type_IV = "Type IV (BYM2 space x RW2 time)")
    # generic: any CSV -> gt (numeric columns rounded; whole numbers kept whole)
    gen <- function(df, title, name, notes = NULL, dec = 3, subtitle = NULL) {
      df <- as.data.frame(df); numc <- names(df)[vapply(df, is.numeric, logical(1))]
      whole <- numc[vapply(df[numc], function(x) all(is.na(x) | abs(x - round(x)) < 1e-9), logical(1))]
      frac <- setdiff(numc, whole)
      tb <- gt(df)
      if (length(whole)) tb <- fmt_number(tb, columns = tidyselect::all_of(whole), decimals = 0, use_seps = FALSE)
      if (length(frac))  tb <- fmt_number(tb, columns = tidyselect::all_of(frac),  decimals = dec)
      tb <- sub_missing(tb, missing_text = "")
      finish(pub(tb, title, subtitle, notes), name)
    }
    
    cat("\n================ gt TABLES ================\n")
    zs <- rd(paste0("TableS3_zone_summary", SUFFIX, ".csv"))     # no fallback to the old unsuffixed file
    zname <- if (!is.null(zs)) setNames(zs$zone, zs$admin2_clean) else NULL
    
    # ---------------- Table 2: model comparison ----------------
    t2 <- function(f, part, name) { d <- rd(f); if (is.null(d)) return(invisible()); d <- as.data.frame(d)
    d$Model <- ifelse(d$Model %in% names(model_lab), model_lab[d$Model], d$Model)
    d <- d[order(match(names(model_lab)[match(d$Model, model_lab)], names(model_lab))), ]
    keep <- intersect(c("Model", "WAIC", "DIC", "p_WAIC", "mean_logCPO", "CPO_failures", "mean_logCPO_common", "AUC_insample", "selected"), names(d)); d <- d[, keep, drop = FALSE]
    if ("selected" %in% names(d)) d$selected <- ifelse(d$selected %in% c(TRUE, "TRUE"), "Yes", "")
    tb <- gt(d)
    tb <- fmt_number(tb, columns = tidyselect::any_of(c("WAIC", "DIC", "p_WAIC")), decimals = 1)
    tb <- fmt_number(tb, columns = tidyselect::any_of(c("mean_logCPO", "mean_logCPO_common")), decimals = 3)
    tb <- fmt_number(tb, columns = tidyselect::any_of("AUC_insample"), decimals = 3)
    tb <- cols_label(tb, .list = as.list(c(Model = "Model", WAIC = "WAIC", DIC = "DIC", p_WAIC = "p(WAIC)", mean_logCPO = "Mean log-CPO", CPO_failures = "CPO failures",
                                           mean_logCPO_common = "Mean log-CPO (common set)", AUC_insample = "AUC (in-sample)", selected = "Selected")[keep]))
    finish(pub(tb, paste0("Table 2. Model comparison: ", part), notes = c(
      "Selected by leave-one-out log-CPO (ties go to the simpler model). p(WAIC) > 0.4 for individual observations flags WAIC as unreliable.",
      "In-sample AUC of the interaction models reflects over-fitting; see the held-out validation table.")), name) }
    safe("Table 2 occurrence", t2(paste0("Table2_occurrence_metrics", SUFFIX, ".csv"), "occurrence", "Table2a_model_comparison_occurrence"))
    safe("Table 2 intensity",  t2(paste0("Table2_intensity_metrics", SUFFIX, ".csv"),  "intensity",  "Table2b_model_comparison_intensity"))
    
    # ---------------- Table 3: parameters of the selected models ----------------
    safe("Table 3", { d <- rd(paste0("Table3_parameters", SUFFIX, ".csv")); if (!is.null(d)) {
      setnames(d, c("0.025quant", "0.975quant"), c("lo", "hi"), skip_absent = TRUE)
      lab <- c("(Intercept)" = "Intercept", sin_month = "sin(2*pi*month/12)", cos_month = "cos(2*pi*month/12)", sin2_month = "sin(4*pi*month/12)",
               cos2_month = "cos(4*pi*month/12)", log_pop_s = "log population (scaled)", log_pop_zone_s = "log zone population (scaled)",
               "Precision for ID_area_bym" = "Spatial SD (BYM2)", "Phi for ID_area_bym" = "Spatial mixing parameter (phi)",
               "Precision for ID_time_rw" = "National trend SD (RW2)", "Precision for ID_time_iid" = "Month effect SD (iid)")
      rows <- lapply(seq_len(nrow(d)), function(i) { t <- d$term[i]; m <- d$mean[i]; l <- d$lo[i]; h <- d$hi[i]
      nm <- if (t %in% names(lab)) lab[[t]] else if (grepl("^size", t)) "NB size" else t
      est <- if (grepl("^Precision for", t)) { if (m > 1e6) sprintf("approx. 0 (95%% upper bound %.2f)", 1 / sqrt(l)) else ci(1 / sqrt(m), 1 / sqrt(h), 1 / sqrt(l), 2) }
      else ci(m, l, h, 3)
      out <- data.table(Component = d$component[i], Parameter = nm, Estimate = est)
      if (grepl("^size", t)) out <- rbind(out, data.table(Component = d$component[i], Parameter = "Overdispersion alpha = 1/size", Estimate = ci(1 / m, 1 / h, 1 / l, 2)))
      out })
      d3 <- rbindlist(rows); d3[, Component := ifelse(Component == "occurrence", "Occurrence (outbreak yes/no)", "Intensity (size of outbreak)")]
      tb <- gt(d3, groupname_col = "Component"); tb <- cols_label(tb, Parameter = "Parameter", Estimate = "Posterior mean (95% CrI)")
      finish(pub(tb, "Table 3. Parameters of the selected hurdle models", notes = c(
        "Fixed effects on the log-odds (occurrence) and log-mean (intensity) scales; population coefficients are per SD of log population.",
        "Random-effect standard deviations are computed from the posterior mean precision; intervals use the precision quantiles.",
        "phi = structured share of the spatial variance, not a transmission measure. For intensity, phi is weakly identified.")), "Table3_parameters") } })
    
    # ---------------- Table 4: zone results ----------------
    if (!is.null(zs)) {
      zt <- copy(zs); zt[, `:=`(pos = N_outbreak_months,
                                P = ci(pi, pi_lo, pi_hi, 1, 100), M = ci(mu, mu_lo, mu_hi, 1), M100 = ci(mu_100k, mu_100k_lo, mu_100k_hi, 1), B = ci(burden_100k, burden_lo, burden_hi, 1))]
      top <- function(ord, cols, title, name, note, filt = rep(TRUE, nrow(zt)), n = 10) { d <- zt[filt][order(ord[filt])][1:n]
      tb <- gt(d[, c("zone", cols, "pos", "total_cases"), with = FALSE])
      lbl <- list(zone = "Zone", P = "Outbreak probability, % (95% CrI)", M = "Expected cases per outbreak month (95% CrI)",
                  M100 = "Cases per 100,000 per outbreak month (95% CrI)", B = "Expected cases per 100,000 per year (95% CrI)",
                  pos = "Months with cases", total_cases = "Total cases")
      tb <- cols_label(tb, .list = lbl[c("zone", cols, "pos", "total_cases")]); tb <- fmt_number(tb, columns = tidyselect::all_of("total_cases"), decimals = 0)
      finish(pub(tb, title, notes = note), name) }
      safe("Table 4a", top(-zt$pi, "P", "Table 4a. Ten zones with the highest estimated probability of an outbreak per zone-month", "Table4a_top10_probability",
                           "Probability for a zone-month, averaged over the study period. Intervals are 95% credible intervals; ranks are uncertain, so report them with intervals."))
      safe("Table 4b", top(-zt$mu, "M", "Table 4b. Ten zones with the largest expected cases per outbreak month (zones with at least 3 months with cases)", "Table4b_top10_size",
                           c("Absolute counts are larger in larger zones; outbreak size increases only weakly with population (exponent about 0.1).", "Zones with fewer than 3 positive months are excluded because their estimate relies mainly on neighbouring zones."), filt = zt$pos >= 3))
      safe("Table 4c", top(-zt$burden_100k, "B", "Table 4c. Ten zones with the highest expected annual caseload per 100,000", "Table4c_top10_burden",
                           "Expected caseload = probability of an outbreak x expected cases per outbreak month, per 100,000 population per year."))
      safe("Table S3", { d <- zt[order(-burden_100k)][, flag := fifelse(pos == 0, "none (estimate from neighbours)", fifelse(pos <= 2, "1-2 months (mainly neighbours)", "3 or more months"))]
      tb <- gt(d[, .(zone, P, M, M100, B, pos, total_cases, flag)])
      tb <- cols_label(tb, zone = "Zone", P = "Outbreak probability, % (95% CrI)", M = "Cases per outbreak month (95% CrI)", M100 = "Cases per 100,000 per outbreak month (95% CrI)",
                       B = "Expected cases per 100,000 per year (95% CrI)", pos = "Months with cases", total_cases = "Total cases", flag = "Data behind the intensity estimate")
      finish(pub(tb, "Table S3. Zone-level estimates for all 92 zones", notes = "Sorted by expected annual caseload per 100,000."), "TableS3_all_zones") })
    }
    
    # ---------------- Waves, seasonality ----------------
    safe("Waves caseload", { d <- rd("Table_waves_national_caseload.csv"); if (!is.null(d)) gen(d, "Table 5a. Waves in national expected caseload (occurrence x intensity)", "Table5a_waves_caseload",
                                                                                                notes = c("A wave is a run of months in which the model's expected national caseload is credibly above the average month; gaps of up to 4 months are merged.",
                                                                                                          "Observed cases are the reported cases in the wave months.")) })
    safe("Waves intensity", { d <- rd("Table_waves.csv"); if (!is.null(d)) gen(d, "Table 5b. Periods of larger outbreak size (national intensity term)", "Table5b_waves_intensity",
                                                                               notes = "Months in which the national intensity trend is credibly above the period mean. This describes outbreak size, not how many zones were affected.") })
    safe("Season window", { d <- rd(paste0("Table_S_seasonal_risk_window", SUFFIX, ".csv")); if (!is.null(d)) gen(d, "Table 6a. Seasonal peak and risk window", "Table6a_seasonal_window",
                                                                                                                  notes = "Window = consecutive months in which the ratio to the average month is above 1 (or credibly above 1), wrapping across December-January.", dec = 2) })
    safe("Season curve", { d <- rd(paste0("Table_S_seasonal_curve_by_month", SUFFIX, ".csv")); if (!is.null(d)) {
      d[, est := ci(ratio, lo, hi, 2)]; d[, month_f := factor(month_name, levels = month.abb)]; w <- dcast(d, month_f ~ part, value.var = "est"); setnames(w, 1, "Month")
      finish(pub(gt(w), "Table 6b. Seasonal curve by month (ratio to the average month, 95% CrI)", notes = "Occurrence: odds ratio; intensity: ratio of expected size."), "Table6b_seasonal_curve") } })
    
    # ---------------- Association between occurrence and size ----------------
    safe("D7", { d <- rd("Table_S_decoupling_population_adjustment.csv"); if (!is.null(d)) {
      o <- data.table(Occurrence = ifelse(d$occurrence_population, "with population", "without"), Intensity = ifelse(d$intensity_population, "with population", "without"),
                      beta = ci(d$beta, d$beta_lo, d$beta_hi, 2), cor = ci(d$cor, d$cor_lo, d$cor_hi, 2), exponent = ifelse(is.na(d$population_exponent_of_size), "", sprintf("%.2f", d$population_exponent_of_size)), WAIC = sprintf("%.1f", d$WAIC))
      tb <- cols_label(gt(o), Occurrence = "Population term in occurrence", Intensity = "Population term in intensity", beta = "Shared-component beta (95% CrI)",
                       cor = "Spatial correlation (95% CrI)", exponent = "Population exponent of size", WAIC = "WAIC")
      finish(pub(tb, "Table 7. Association between zone-level outbreak occurrence and outbreak size, by population adjustment", notes = c(
        "The row with population in both parts is the primary analysis. Without population adjustment the two opposite population effects cancel.",
        "Association between zone-level effects only; not a causal or transmission measure.")), "Table7_association_population") } })
    safe("D1-D3", { d <- rd("Table_S_decoupling_correlations.csv"); if (!is.null(d)) gen(d, "Table S7a. Posterior correlation between how often and how large (zone level)", "TableS7a_correlations", dec = 3) })
    safe("D4", { d <- rd("Table_S_decoupling_rank_overlap.csv"); if (!is.null(d)) gen(d, "Table S7b. Overlap of the top zones for occurrence and size", "TableS7b_rank_overlap", dec = 3,
                                                                                      notes = "Chance overlap is shown as expected_by_chance. A moderate correlation (about 0.5) still gives mostly different top-10 lists.") })
    safe("D5", { d <- rd("Table_S_decoupling_beta_prior.csv"); if (!is.null(d)) gen(d, "Table S7c. Joint-model beta under three priors", "TableS7c_beta_prior", dec = 3) })
    
    # ---------------- Persistence extension ----------------
    safe("Persist forward", { d <- rd("Table_S_persistence_forward_validation.csv"); if (!is.null(d)) {
      o <- d[part == "occurrence", .(model, AUC, logloss, Brier, Brier_climatology)]; i <- d[part == "intensity", .(model, logscore_truncNB, MAE_log, n_test_positive)]
      gen(o, "Table 8a. Forward validation 2022-2024, occurrence (trained on 2015-2021)", "Table8a_persistence_forward_occurrence", dec = 3,
          notes = c("One-step-ahead: the previous month's status is observed. Brier_climatology uses the training-period frequency.", "Persistence = association between consecutive months, not transmission."))
      gen(i, "Table 8b. Forward validation 2022-2024, intensity", "Table8b_persistence_forward_intensity", dec = 2,
          notes = "Forward prediction of outbreak size remains poor for all models because the national trend cannot be extrapolated into a new wave.") } })
    safe("Persist resid", { d <- rd("Table_S_persistence_residuals.csv"); if (!is.null(d)) gen(d, "Table 8c. Lag-1 residual autocorrelation within zones", "Table8c_persistence_residuals", dec = 3) })
    safe("Persist zones", { d <- rd("Table_S_persistence_zone_effects.csv"); if (!is.null(d)) gen(d, "Table 8d. Effect of persistence terms on zone effects and phi", "Table8d_persistence_zone_effects", dec = 3) })
    
    # ---------------- Model checks and sensitivity ----------------
    safe("Predeployment", { d <- rd("Table_predeployment_checks.csv"); if (!is.null(d)) {
      tb <- gt(as.data.frame(d)); tb <- tryCatch(tab_style(tb, style = cell_fill(color = "#D9EAD3"), locations = cells_body(columns = tidyselect::all_of("status"), rows = status == "PASS")), error = function(e) tb)
      tb <- tryCatch(tab_style(tb, style = cell_fill(color = "#F4CCCC"), locations = cells_body(columns = tidyselect::all_of("status"), rows = status == "FAIL")), error = function(e) tb)
      finish(pub(tb, "Table S4. Pre-deployment checks of the selected models", notes = "Failed checks are reported as limitations (no month-to-month persistence in the Baseline; weak forward prediction)."), "TableS4_predeployment_checks") } })
    safe("Prior sens", { d <- rd("Table_S_prior_sensitivity.csv"); if (!is.null(d)) gen(d, "Table S5. Sensitivity of the zone effects to the PC priors", "TableS5_prior_sensitivity", dec = 3) })
    safe("Sel priors", { d <- rd("Selected_priors.csv"); if (!is.null(d)) gen(d, "Table S5b. Priors selected by the grid search", "TableS5b_selected_priors", dec = 3) })
    safe("Hurdle", { d <- rd("Table_S_hurdle_vs_zeroinflated.csv"); if (!is.null(d)) gen(d, "Table S6a. Hurdle versus zero-inflated models", "TableS6a_hurdle_vs_zeroinflated", dec = 1,
                                                                                         notes = "Compare WAIC and the summed log-CPO; mean log-CPO is not comparable across models with CPO failures.") })
    safe("Pop adj", { d <- rd("Table_S_population_adjustment_spatial_confounding.csv"); if (!is.null(d)) gen(d, "Table S6b. Zone effects with and without the population term", "TableS6b_population_adjustment", dec = 3) })
    safe("Trunc", { d <- rd("Table_S_zero_truncated_intensity.csv"); if (!is.null(d)) gen(d, "Table S6c. Intensity likelihood: zero-truncated versus ordinary negative binomial (main model first)", "TableS6c_zero_truncated", dec = 3) })
    safe("Trunc rank", { d <- rd("Table_S_zero_truncated_ranking.csv"); if (!is.null(d)) gen(d, "Table S6d. Agreement of zone results between the two intensity likelihoods", "TableS6d_zero_truncated_ranking", dec = 3) })
    safe("Overfit cx", { d <- rd("Table_S_overfit_complexity.csv"); if (!is.null(d)) gen(d, "Table S8a. Model complexity (effective parameters)", "TableS8a_overfit_complexity", dec = 3) })
    safe("Overfit ho", { d <- rd("Table_S_overfit_heldout_validation.csv"); if (!is.null(d)) gen(d, "Table S8b. Held-out validation of the five candidate models", "TableS8b_heldout_validation", dec = 3,
                                                                                                 notes = "random20: 20% of zone-months hidden at random; forward: trained on 2015-2021, tested on 2022-2024.") })
    safe("Overfit sd", { d <- rd("Table_S_overfit_random_effect_sd.csv"); if (!is.null(d)) gen(d, "Table S8c. Standard deviations of the random effects, by model", "TableS8c_random_effect_sd", dec = 3) })
    safe("PIT", { d <- rd("Table_S_PIT_all_models.csv"); if (!is.null(d)) gen(d, "Table S9. Randomised PIT summary for all candidate models", "TableS9_PIT", dec = 3) })
    safe("Silent", { d <- rd("Table_S_silent_months_sensitivity.csv"); if (!is.null(d)) gen(d, "Table S10. Sensitivity to the 17 months without cases", "TableS10_silent_months", dec = 3) })
    safe("Moran", { d <- rd("Table_S_global_morans_I.csv"); if (!is.null(d)) gen(d, "Table S11a. Global Moran's I", "TableS11a_morans_I", dec = 3) })
    safe("Rank", { d <- rd("Table_S_rank_stability_top10.csv"); if (!is.null(d)) gen(d, "Table S11b. Stability of the top-10 zones", "TableS11b_rank_stability", dec = 3) })
    safe("Obs pred", { d <- rd("Table_S_observed_vs_predicted_uncertainty.csv"); if (!is.null(d)) { if (!is.null(zname) && "ID_area" %in% names(d)) d[, zone := zs$zone[match(ID_area, zs$ID_area)]]
      gen(d, "Table S12. Observed versus model-predicted annual incidence per 100,000, by zone", "TableS12_observed_vs_predicted", dec = 2) } })
    
    cat("\nTables done. Word files are in:", TAB_DIR, "\n")
    
  } else cat("STEP 10 skipped: RUN_TABLES_GT is FALSE or the gt package is not installed (install.packages('gt')).\n")
  
  ################################################################################
  # STEP 11 — PUBLICATION MAPS. Look at them in the Plots pane first; set SAVE_MAPS <- TRUE inside this block (or in Stage_maps_publication.R) to save.
  ################################################################################
  RUN_MAPS_PUB <- TRUE
  if (RUN_MAPS_PUB) {
    # ==============================================================================
    # STAGE: PUBLICATION MAPS  (shown in the RStudio Plots pane; saved only when SAVE_MAPS <- TRUE)
    # Reads the pipeline's CSV files and the shapefile, so it can be rerun on its own to adjust the look.
    # Needs: sf, ggplot2, dplyr, data.table. Optional (used when installed): patchwork, ggrepel, ggspatial.
    # Rules: zones with no positive month are GREY on the intensity map (no estimate of their own); zones with 1-2 positive months get a black outline.
    # Wording: "estimated risk", "clusters of higher/lower estimated risk"; no hotspot / transmission / hub language.
    # ==============================================================================
    suppressPackageStartupMessages({ library(sf); library(ggplot2); library(dplyr); library(data.table) })
    if (!exists("output_path")) output_path <- getOption("cholera.proj", ".")
    if (!exists("shp_path"))    shp_path    <- "F:/Cholera/shape_22/Adm_2/eth_admbnda_adm2_csa_bofedb_2021.shp"
    if (!exists("SUFFIX"))      SUFFIX      <- "_OCCPOP"
    SAVE_MAPS      <- FALSE      # FALSE = only look at the maps in the Plots pane; TRUE = also save PNG (600 dpi) and PDF
    REGION_BORDERS <- TRUE       # draw region (ADM1) outlines if the shapefile has an ADM1 column
    LABEL_TOP      <- 5          # label the top zones on the main maps (needs ggrepel); 0 = no labels
    FIG_DIR <- file.path(output_path, "figures_publication"); if (SAVE_MAPS) dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
    has <- function(p) requireNamespace(p, quietly = TRUE)
    
    # ---------------- data ----------------
    zs <- fread(file.path(output_path, paste0("TableS3_zone_summary", SUFFIX, ".csv")))
    shp <- st_make_valid(st_read(shp_path, quiet = TRUE))
    shp$admin2_clean <- toupper(gsub("[^A-Z0-9]", "", shp$ADM2_PCODE))
    shp <- left_join(shp, as.data.frame(zs), by = "admin2_clean")
    stopifnot("zone summary does not match the shapefile" = sum(!is.na(shp$pi)) > 80)
    shp$pi_pct   <- 100 * shp$pi
    shp$has_int  <- shp$N_outbreak_months > 0
    shp$few      <- shp$N_outbreak_months %in% 1:2
    shp$mu_plot  <- ifelse(shp$has_int, shp$mu, NA_real_)
    shp$mu100_plot <- ifelse(shp$has_int, shp$mu_100k, NA_real_)
    shp$burden_width <- (shp$burden_hi - shp$burden_lo) / pmax(shp$burden_100k, 1e-6)
    shp$mu_width <- ifelse(shp$has_int, (shp$mu_hi - shp$mu_lo) / pmax(shp$mu, 1e-6), NA_real_)
    reg <- NULL
    if (REGION_BORDERS) reg <- tryCatch({ col <- intersect(c("ADM1_EN", "ADM1_PCODE"), names(shp))[1]
    if (is.na(col)) NULL else shp %>% group_by(.data[[col]]) %>% summarise(.groups = "drop") }, error = function(e) NULL)
    cat(sprintf("Zones: %d | no positive month (grey on intensity maps): %d | 1-2 positive months (outlined): %d\n", nrow(shp), sum(!shp$has_int), sum(shp$few)))
    
    # ---------------- look ----------------
    tr_arg <- if (utils::packageVersion("ggplot2") >= "3.5.0") "transform" else "trans"
    theme_map <- function(b = 11) theme_void(base_size = b) +
      theme(plot.title = element_text(face = "bold", size = b + 1, hjust = 0, margin = margin(b = 4)),
            legend.position = "bottom", legend.key.width = unit(1.8, "cm"), legend.key.height = unit(0.28, "cm"),
            legend.title = element_text(size = b - 1, vjust = 0.9), legend.text = element_text(size = b - 2),
            plot.caption = element_text(size = b - 3, hjust = 0, colour = "grey30", margin = margin(t = 4)),
            plot.margin = margin(4, 4, 4, 4))
    decor <- function(p) { if (!is.null(reg)) p <- p + geom_sf(data = reg, fill = NA, colour = "grey20", linewidth = 0.3)
    if (has("ggspatial")) p <- p + ggspatial::annotation_scale(location = "bl", width_hint = 0.25, height = unit(0.12, "cm"), text_cex = 0.6) +
        ggspatial::annotation_north_arrow(location = "tr", which_north = "true", height = unit(0.8, "cm"), width = unit(0.6, "cm"), style = ggspatial::north_arrow_minimal(text_size = 7))
    p }
    cmap <- function(var, title, name, option = "mako", dir = -1, trans = "identity", outline_few = FALSE, label_top = 0, caption = NULL, limits = NULL) {
      sc <- do.call(scale_fill_viridis_c, c(list(option = option, direction = dir, name = name, na.value = "grey85", limits = limits,
                                                 guide = guide_colourbar(title.position = "top", ticks.colour = "white")), setNames(list(trans), tr_arg)))
      p <- ggplot(shp) + geom_sf(aes(fill = .data[[var]]), colour = "white", linewidth = 0.15) + sc
      if (outline_few && any(shp$few)) p <- p + geom_sf(data = shp[shp$few, ], fill = NA, colour = "black", linewidth = 0.5)
      p <- decor(p)
      if (label_top > 0 && has("ggrepel")) { v <- shp[[var]]; idx <- order(-v, na.last = NA)[seq_len(label_top)]
      p <- p + ggrepel::geom_text_repel(data = shp[idx, ], aes(label = zone, geometry = geometry), stat = "sf_coordinates", size = 2.6, fontface = "italic",
                                        min.segment.length = 0, segment.size = 0.25, segment.colour = "grey30", max.overlaps = Inf, seed = 1, box.padding = 0.35, bg.color = "white", bg.r = 0.12) }
      p + labs(title = title, caption = caption) + theme_map() }
    save_map <- function(p, name, w = 10, h = 6.2) {
      print(p)
      if (SAVE_MAPS) { ggsave(file.path(FIG_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 600, bg = "white")
        tryCatch(ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p, width = w, height = h, device = PDF_DEVICE, bg = "white"), error = function(e) ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p, width = w, height = h))
        cat("  saved:", file.path(FIG_DIR, name), "(.png, .pdf)\n") } }
    panels <- function(plist, name, w = 10, h = 6.2, ncol = 2) {
      if (has("patchwork")) { p <- patchwork::wrap_plots(plist, ncol = ncol) + patchwork::plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") &
        theme(plot.tag = element_text(face = "bold", size = 13)); save_map(p, name, w, h)
      } else { cat("patchwork not installed: panels are drawn and saved one by one.\n"); for (i in seq_along(plist)) save_map(plist[[i]], paste0(name, "_", letters[i]), w / ncol, h) } }
    cap_int <- "Grey: no month with reported cases (outcome size not estimable). Black outline: 1-2 months with cases (estimate relies mainly on neighbouring zones)."
    
    # ---------------- Figure 4: occurrence and size ----------------
    f4a <- cmap("pi_pct", "Probability of an outbreak per zone-month", "Posterior mean probability (%)", "mako", -1, "sqrt", label_top = LABEL_TOP)
    f4b <- cmap("mu_plot", "Expected cases per outbreak month", "Posterior mean cases", "rocket", -1, "sqrt", outline_few = TRUE, label_top = LABEL_TOP, caption = cap_int)
    panels(list(f4a, f4b), "Fig4_occurrence_and_size", 11, 6.4)
    
    # ---------------- Figure S4: per-population size and expected caseload ----------------
    f5a <- cmap("mu100_plot", "Cases per 100,000 per outbreak month", "Posterior mean per 100,000", "rocket", -1, "sqrt", outline_few = TRUE, caption = cap_int)
    f5b <- cmap("burden_100k", "Expected cases per 100,000 per year", "Probability x size, per 100,000", "inferno", -1, "sqrt", label_top = LABEL_TOP)
    panels(list(f5a, f5b), "FigS4_size_per100k_and_caseload", 11, 6.4)
    
    # ---------------- Figure S5: uncertainty ----------------
    u1 <- cmap("burden_width", "Uncertainty of expected annual caseload", "Relative 95% CrI width", "magma", -1, "sqrt")
    u2 <- cmap("mu_width", "Uncertainty of outbreak size", "Relative 95% CrI width", "magma", -1, "sqrt", outline_few = TRUE, caption = cap_int)
    panels(list(u1, u2), "FigS5_uncertainty", 11, 6.4)
    
    # ---------------- Figure S6: how much data sit behind each zone ----------------
    shp$data_class <- cut(shp$N_outbreak_months, c(-Inf, 0, 2, 11, 35, Inf), labels = c("0 months", "1-2 months", "3-11 months", "12-35 months", "36+ months"))
    f6 <- ggplot(shp) + geom_sf(aes(fill = data_class), colour = "white", linewidth = 0.15) +
      scale_fill_manual(values = c("0 months" = "grey85", "1-2 months" = "#FDE0C5", "3-11 months" = "#F59E5B", "12-35 months" = "#D6562B", "36+ months" = "#8C1C13"),
                        name = "Months with reported cases (of 120)", drop = FALSE, guide = guide_legend(nrow = 1, title.position = "top"))
    f6 <- decor(f6) + labs(title = "Months with reported cases, by zone, 2015-2024") + theme_map()
    save_map(f6, "FigS6_data_availability", 7, 6.4)
    
    # ---------------- Figure S7: Getis-Ord Gi* classes on 10-year mean incidence ----------------
    gi <- fread(file.path(output_path, "Table_S_Gi_star_10yr_mean_incidence.csv"))
    if (nrow(gi)) { shp$gi_class <- factor(gi$class[match(shp$admin2_clean, gi$zone)], levels = c("Cold 99%", "Cold 95%", "Cold 90%", "Not significant", "Hot 90%", "Hot 95%", "Hot 99%"))
    shp$inc10 <- gi$incidence_100k_year[match(shp$admin2_clean, gi$zone)]
    g1 <- cmap("inc10", "Mean annual incidence per 100,000", "Cases per 100,000 per year", "viridis", -1, "sqrt")
    g2 <- ggplot(shp) + geom_sf(aes(fill = gi_class), colour = "white", linewidth = 0.15) +
      scale_fill_manual(values = c("Cold 99%" = "#08519C", "Cold 95%" = "#4292C6", "Cold 90%" = "#9ECAE1", "Not significant" = "grey88", "Hot 90%" = "#FCAE91", "Hot 95%" = "#EF3B2C", "Hot 99%" = "#99000D"),
                        name = "Gi* class", drop = FALSE, guide = guide_legend(nrow = 2, title.position = "top"))
    g2 <- decor(g2) + labs(title = "Clusters of higher and lower observed incidence (Getis-Ord Gi*)", caption = "No multiple-testing correction; 90% classes are exploratory.") + theme_map()
    panels(list(g1, g2), "FigS7_incidence_and_Gi_star", 11, 6.4) }
    
    cat("\nMaps done. SAVE_MAPS is", SAVE_MAPS, "- set it to TRUE and rerun this script to write PNG (600 dpi) and PDF files to", FIG_DIR, "\n")
    
  } else cat("STEP 11 skipped: RUN_MAPS_PUB is FALSE.\n")
  
}   # end PART 1
RUN <- function(k) isTRUE(getOption(paste0("cholera.run_", k), FALSE))   # re-defined: STEP 1 starts with rm(list = ls()), which removes it

# ##############################################################################
# PART 2 - EXTRAS (use the objects of PART 1 that are in memory; each block is isolated in local())
# ##############################################################################
if (RUN("part2")) {
  need <- c("dat", "fit_o", "fit_i", "NN", "shape", "output_path")
  if (!all(sapply(need, exists))) stop("PART 2 needs these objects from PART 1 in memory: ", paste(need[!sapply(need, exists)], collapse = ", "), ". Run PART 1 in this session first.")
  cat("\n######## PART 2a: residual semivariograms ########\n");           local({
    # ==============================================================================
    # Stage_variogram_residuals.R  -  empirical semivariograms of the model residuals (spatial check, complements Moran's I)
    # Run in the SAME R session after Step 3 of the main pipeline (needs: dat, fit_o, fit_i, NN, shape, pearson_int, output_path; zone order = shape rows).
    # Reading: if the model captured the spatial structure, the residual semivariogram is flat and stays inside the grey envelope
    # (envelope = 199 random re-arrangements of the same residuals over the zones, i.e. no spatial dependence).
    # Install once: install.packages(c("gstat", "sf"))
    # Wording: "residual spatial dependence"; no transmission language.
    # ==============================================================================
    suppressPackageStartupMessages({ library(gstat); library(sf); library(data.table); library(ggplot2) })
    stopifnot(all(sapply(c("dat", "fit_o", "fit_i", "NN", "shape", "output_path"), exists)))
    if (!exists("SAVE_VARIO")) SAVE_VARIO <- FALSE
    set.seed(2026)
    
    # ---- zone residuals (same definitions as Step 3, check 5)
    p  <- fit_o$summary.fitted.values$mean[1:NN]
    r_occ <- as.numeric(tapply(dat$Y_bin - p, dat$ID_area, sum) / sqrt(tapply(p * (1 - p), dat$ID_area, sum)))
    r_int_m <- ifelse(is.na(dat$O_int), NA_real_, pearson_int(fit_i, dat$O_int, 1:NN))
    r_int <- as.numeric(tapply(r_int_m, dat$ID_area, function(x) { x <- x[!is.na(x)]; if (length(x) >= 3) sum(x) / sqrt(length(x)) else NA_real_ }))
    crude <- as.numeric(tapply(dat$O, dat$ID_area, sum) / tapply(dat$monthly_population, dat$ID_area, sum) * 1e5 * 12)   # cases per 100,000 per year
    zone_df <- data.frame(ID_area = 1:length(r_occ), r_occ = r_occ, r_int = r_int, log_crude = log1p(crude))
    
    # ---- zone centroids in metres (UTM 37N), same order as the model zones
    pts <- suppressWarnings(st_centroid(st_transform(st_geometry(shape), 20137)))
    pts <- st_sf(zone_df, geometry = pts[match(zone_df$ID_area, seq_along(pts))])
    mx <- max(dist(st_coordinates(pts))); cutoff <- mx / 3; width <- cutoff / 10
    cat(sprintf("Maximum distance between centroids: %.0f km | variogram cutoff %.0f km | bin width %.0f km\n", mx / 1000, cutoff / 1000, width / 1000))
    
    vg <- function(var) {
      d <- pts[!is.na(pts[[var]]), ]; f <- as.formula(paste(var, "~ 1"))
      v <- variogram(f, d, cutoff = cutoff, width = width)
      env <- replicate(199, { dp <- d; dp[[var]] <- sample(dp[[var]]); variogram(f, dp, cutoff = cutoff, width = width)$gamma })
      data.table(variable = var, dist_km = v$dist / 1000, gamma = v$gamma, np = v$np,
                 lo = apply(env, 1, quantile, 0.025), hi = apply(env, 1, quantile, 0.975), n_zones = nrow(d))
    }
    V <- rbindlist(lapply(c("r_occ", "r_int", "log_crude"), vg))
    V[, label := factor(variable, c("r_occ", "r_int", "log_crude"),
                        c("Occurrence residuals", "Intensity residuals (zones with >= 3 outbreak months)", "Crude log incidence (data, no model)"))]
    V[, outside := gamma < lo | gamma > hi]
    fwrite(V[, -"label"], file.path(output_path, "Table_S_variogram_residuals.csv"))
    cat("\n--- Semivariogram bins outside the no-dependence envelope ---\n")
    print(V[, .(bins = .N, outside_envelope = sum(outside), first_bins_outside_km = paste(round(head(dist_km[outside], 3)), collapse = ", ")), by = variable])
    
    g <- ggplot(V, aes(dist_km)) + geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey80", alpha = .8) +
      geom_line(aes(y = gamma), colour = "#D55E00", linewidth = .9) + geom_point(aes(y = gamma, shape = outside), colour = "#D55E00", size = 2) +
      scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), guide = "none") +
      facet_wrap(~ label, scales = "free_y", ncol = 1) +
      labs(x = "Distance between zone centroids (km)", y = "Semivariance",
           title = "Empirical semivariograms", subtitle = "Grey band = range expected with no spatial dependence (199 permutations); filled points fall outside it") +
      theme_minimal(base_size = 10)
    print(g)
    if (SAVE_VARIO) { ggsave(file.path(output_path, "Fig_S_variogram_residuals.png"), g, width = 6.5, height = 8, dpi = 600); ggsave(file.path(output_path, "Fig_S_variogram_residuals.pdf"), g, width = 6.5, height = 8) }
    cat("\nREAD: residual semivariograms inside the band = the BYM2 effect captured the spatial dependence (agrees with Moran's I on the residuals).\n",
        " The crude-data variogram typically rises with distance (nearby zones are more alike) and shows what the model had to explain.\n")
    
  })
  cat("\n######## PART 2b: seasonal pattern, absolute and by year ########\n"); local({
    # ==============================================================================
    # Stage_seasonal_absolute_and_by_year.R  -  absolute seasonal probabilities and the year-by-year pattern
    # Reads only zone_month_results_OCCPOP.csv (written by the main pipeline). Run on its own, after the pipeline.
    # Caution 1: seasonal ratios compare months with each other; this script gives the ABSOLUTE chance of an outbreak per zone-month and expected size.
    # Caution 2: the seasonal curve is the average year; this script shows each year separately (the epidemic waves).
    # Wording: association with time of year; no causal statement about rainfall.
    # ==============================================================================
    suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(gt) })
    if (!exists("output_path")) output_path <- getOption("cholera.proj", ".")
    if (!exists("SAVE_SEAS")) SAVE_SEAS <- FALSE
    out_s <- file.path(output_path, "seasonal"); if (SAVE_SEAS) dir.create(out_s, showWarnings = FALSE, recursive = TRUE)
    if (!exists("SUFFIX")) SUFFIX <- "_OCCPOP"
    zf <- file.path(output_path, paste0("zone_month_results", SUFFIX, ".csv"))
    if (!file.exists(zf)) stop("Not found: ", zf, "\nThe main pipeline has not written it yet (it is written after the models and draws). Files present: ",
                               paste(list.files(output_path, pattern = "zone_month_results"), collapse = ", "))
    cat("Using:", zf, "| last written:", format(file.info(zf)$mtime), "(check this is from your latest run)\n")
    zm <- fread(zf)
    stopifnot("pi" %in% names(zm))
    
    # ---- 1. Absolute seasonal table (average year): per calendar month
    by_m <- zm[, .(model_prob_pct = 100 * mean(pi), observed_pct = 100 * mean(O > 0),
                   expected_size_if_outbreak = mean(mu[O > 0]), observed_mean_size = mean(O[O > 0]),
                   zone_months_with_cases = sum(O > 0)), by = month][order(month)]
    by_m[, month_name := month.abb[month]]
    cat("\n--- Absolute seasonal pattern (average over 2015-2024) ---\n")
    print(by_m[, .(month_name, model_prob_pct = round(model_prob_pct, 1), observed_pct = round(observed_pct, 1),
                   expected_size_if_outbreak = round(expected_size_if_outbreak, 1), observed_mean_size = round(observed_mean_size, 1))])
    cat(sprintf("\nAll months: model %.1f%% | observed %.1f%%\n", 100 * mean(zm$pi), 100 * mean(zm$O > 0)))
    fwrite(by_m, file.path(output_path, "Table_S_seasonal_absolute.csv"))
    
    # ---- 2. Year-by-month: share of zones with cases (observed) and model probability
    ym <- zm[, .(observed_pct = 100 * mean(O > 0), model_pct = 100 * mean(pi), cases = sum(O)), by = .(year, month)]
    fwrite(ym, file.path(output_path, "Table_S_seasonal_by_year.csv"))
    cat("\n--- Peak month in each year (observed share of zones with cases) ---\n")
    pk <- ym[, .SD[which.max(observed_pct)], by = year][, .(year, peak_month = month.abb[month], observed_pct = round(observed_pct, 1), cases)]
    print(pk); fwrite(pk, file.path(output_path, "Table_S_seasonal_peak_by_year.csv"))
    
    g1 <- ggplot(ym, aes(factor(month, labels = month.abb), factor(year), fill = observed_pct)) + geom_tile(colour = "white") +
      scale_fill_viridis_c(option = "C", direction = -1, name = "% of zones\nwith cases") + labs(x = NULL, y = NULL, title = "Share of zones with cholera cases, by year and month") +
      theme_minimal(base_size = 10) + theme(panel.grid = element_blank())
    g2 <- ggplot(by_m, aes(factor(month, labels = month.abb))) + geom_col(aes(y = observed_pct), fill = "grey80") +
      geom_line(aes(y = model_prob_pct, group = 1), colour = "#D55E00", linewidth = 1) + geom_point(aes(y = model_prob_pct), colour = "#D55E00") +
      labs(x = NULL, y = "% of zone-months with an outbreak", title = "Average year: observed (bars) and model probability (line)") + theme_minimal(base_size = 10)
    print(g2); print(g1)
    if (SAVE_SEAS) { ggsave(file.path(out_s, "seasonal_average_year.png"), g2, width = 7, height = 4, dpi = 600); ggsave(file.path(out_s, "seasonal_by_year.png"), g1, width = 7, height = 5, dpi = 600) }
    
    tb <- gt(by_m[, .(month_name, model_prob_pct, observed_pct, expected_size_if_outbreak, observed_mean_size)]) |>
      fmt_number(everything(), decimals = 1) |>
      cols_label(month_name = "Month", model_prob_pct = "Model probability of an outbreak (%)", observed_pct = "Observed zone-months with cases (%)",
                 expected_size_if_outbreak = "Model expected cases per outbreak month", observed_mean_size = "Observed mean cases per outbreak month") |>
      tab_header(title = "Table S. Seasonal pattern in absolute terms, Ethiopia, 2015-2024") |>
      tab_footnote("Average over all zones and years. Association with time of year, not a causal effect of season.")
    print(tb); if (SAVE_SEAS) tryCatch(gtsave(tb, file.path(out_s, "Table_seasonal_absolute.docx")), error = function(e) message("docx failed: ", conditionMessage(e)))
    
  })
  cat("\n######## PART 2c: descriptive Figure 1, Figure 2, Table 1 ########\n"); local({
    # ==============================================================================
    # Stage_descriptive_Fig1_Fig2_Table1.R  -  study-area map (Figure 1), epidemic curve (Figure 2), yearly descriptive table (Table 1, gt)
    # Needs only the saved state file from the main pipeline and the shapefile. Looked at in the Plots pane / Viewer; saved only if SAVE_DESC <- TRUE.
    # Table 1 here covers the zone-month data (cases, zones with cases, zone-months with cases, incidence). Age / sex / case-type columns are NOT in the
    # zone-month data; they must come from the line list (see Handoff data items) - add them to the table when that file is available.
    # ==============================================================================
    suppressPackageStartupMessages({ library(data.table); library(sf); library(ggplot2); library(gt) })
    if (!exists("output_path")) output_path <- getOption("cholera.proj", ".")
    if (!exists("shp_path"))    shp_path    <- "F:/Cholera/shape_22/Adm_2/eth_admbnda_adm2_csa_bofedb_2021.shp"
    if (!exists("SAVE_DESC"))   SAVE_DESC   <- FALSE
    out_d <- file.path(output_path, "descriptive"); if (SAVE_DESC) dir.create(out_d, showWarnings = FALSE, recursive = TRUE)
    S <- readRDS(file.path(output_path, "INLA_Rebuilt_From_Zero_OCCPOP.rds")); dat <- as.data.table(S$data)
    stopifnot(all(c("O", "year", "month", "monthly_population", "ID_area") %in% names(dat)))
    dat[, date := as.Date(sprintf("%d-%02d-01", year, month))]
    
    # ---- Figure 1: study area
    shp <- st_make_valid(st_read(shp_path, quiet = TRUE)); shp$zone <- toupper(gsub("[^A-Z0-9]", "", shp$ADM2_PCODE))
    shp$in_analysis <- shp$zone %in% S$zone_order
    f1 <- ggplot(shp) + geom_sf(aes(fill = in_analysis), colour = "white", linewidth = .1) +
      scale_fill_manual(values = c(`TRUE` = "#4C78A8", `FALSE` = "grey85"), labels = c(`TRUE` = "Included (n = 92)", `FALSE` = "Not included"), name = NULL) +
      labs(title = "Study area: second-level administrative zones, Ethiopia") + theme_void(base_size = 10) + theme(legend.position = "bottom")
    print(f1)
    
    # ---- Figure 2: epidemic curve (monthly national cases) with seasons shaded (Bega Dec-Feb, Belg Mar-May, Kiremt Jun-Sep, Bega Oct-Nov)
    nat <- dat[, .(cases = sum(O), zones_with_cases = sum(O > 0)), by = date]
    season <- function(m) fifelse(m %in% 6:9, "Kiremt (Jun-Sep)", fifelse(m %in% 3:5, "Belg (Mar-May)", "Bega (Oct-Feb)"))
    nat[, season := season(as.integer(format(date, "%m")))]
    f2 <- ggplot(nat, aes(date, cases, fill = season)) + geom_col(width = 28) +
      scale_fill_manual(values = c(`Kiremt (Jun-Sep)` = "#4C78A8", `Belg (Mar-May)` = "#F58518", `Bega (Oct-Feb)` = "grey60"), name = NULL) +
      scale_x_date(date_breaks = "1 year", date_labels = "%Y") + labs(x = NULL, y = "Reported cases per month", title = "Monthly reported cholera cases, Ethiopia, 2015-2024") +
      theme_minimal(base_size = 10) + theme(legend.position = "bottom")
    print(f2)
    if (SAVE_DESC) { ggsave(file.path(out_d, "Fig1_study_area.png"), f1, width = 6, height = 6, dpi = 600); ggsave(file.path(out_d, "Fig2_epidemic_curve.png"), f2, width = 9, height = 4.5, dpi = 600) }
    
    # ---- Table 1 (gt): by year
    yr <- dat[, .(cases = sum(O), zones_with_cases = uniqueN(ID_area[O > 0]), zone_months_with_cases = sum(O > 0),
                  person_years = sum(monthly_population) / 12), by = year]
    yr[, `:=`(share_pct = 100 * cases / sum(cases), incidence_100k = 1e5 * cases / person_years)]
    tot <- data.table(year = "All years", cases = sum(dat$O), zones_with_cases = uniqueN(dat$ID_area[dat$O > 0]), zone_months_with_cases = sum(dat$O > 0),
                      person_years = sum(dat$monthly_population) / 12, share_pct = 100, incidence_100k = 1e5 * sum(dat$O) / (sum(dat$monthly_population) / 12))
    t1 <- rbind(yr[, year := as.character(year)], tot)[order(year == "All years", year)]
    tb <- gt(t1[, .(year, cases, share_pct, zones_with_cases, zone_months_with_cases, incidence_100k)]) |>
      fmt_number(c(cases, zones_with_cases, zone_months_with_cases), decimals = 0) |> fmt_number(c(share_pct, incidence_100k), decimals = 1) |>
      cols_label(year = "Year", cases = "Reported cases", share_pct = "% of all cases", zones_with_cases = "Zones with cases", zone_months_with_cases = "Zone-months with cases", incidence_100k = "Cases per 100,000 person-years") |>
      tab_header(title = "Table 1. Reported cholera cases by year, Ethiopia, 2015-2024") |>
      tab_footnote("Zone-month surveillance data (92 zones x 120 months). Age, sex and case-type breakdowns come from the line list.")
    print(tb); fwrite(t1, file.path(output_path, "Table1_by_year.csv"))
    if (SAVE_DESC) tryCatch(gtsave(tb, file.path(out_d, "Table1_by_year.docx")), error = function(e) message("docx failed: ", conditionMessage(e)))
    
  })
  cat("\n######## PART 2d: numbers for the yellow [XX] placeholders ########\n"); local({
    # ==============================================================================
    # Check_yellow_items.R  -  prints every number that is still yellow ([XX]) in the Methods
    # Run in the SAME R session, after the full pipeline has finished (needs dat, fit_o, fit_i,
    # output_path, SUFFIX). Nothing is refitted except one intercept-only NB (seconds).
    # ==============================================================================
    library(data.table)
    if (!exists("SUFFIX")) SUFFIX <- "_OCCPOP"
    stopifnot(exists("dat"), exists("output_path"))
    rd <- function(f) { p <- file.path(output_path, f); if (file.exists(p)) fread(p) else { cat("   (file not found:", f, ")\n"); NULL } }
    hd <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
    
    # ---- 1. Exploratory data distribution (the five numbers in the first yellow note) -----------------
    hd("1. DATA DISTRIBUTION  (Methods: 'Exploratory analysis of data distribution')")
    O <- dat$O
    cat(sprintf("Zone-months: %d | zones: %d | months: %d | total cases: %s\n",
                nrow(dat), uniqueN(dat$ID_area), uniqueN(dat$ID_time), format(sum(O), big.mark = ",")))
    cat(sprintf("Mean = %.1f | variance = %.1f | dispersion ratio (var/mean) = %.1f\n", mean(O), var(O), var(O) / mean(O)))
    cat(sprintf("Zeros: %d of %d (%.1f%%)\n", sum(O == 0), length(O), 100 * mean(O == 0)))
    if (requireNamespace("MASS", quietly = TRUE)) {
      nb0 <- tryCatch(MASS::glm.nb(O ~ 1), error = function(e) NULL)
      if (!is.null(nb0)) cat(sprintf("Intercept-only negative binomial: theta = %.3f | expected zeros = %.0f (observed %d)\n",
                                     nb0$theta, sum(dnbinom(0, mu = fitted(nb0), size = nb0$theta)), sum(O == 0)))
    } else cat("Install MASS for the intercept-only NB: install.packages('MASS')\n")
    cat(sprintf("Minimum expected count E: %.3f (Methods say 0.416) | any E = 0: %s\n", min(dat$E), any(dat$E == 0)))
    
    # ---- 2. Which priors were actually used -----------------------------------------------------------------------
    hd("2. PRIORS IN FORCE  (Methods: 'Prior distributions')")
    for (v in c("U_RW_OCC_EFF", "U_RW_INT_EFF", "U_INT", "OCC_POP", "INT_POP"))
      if (exists(v)) cat(sprintf("%-14s = %s\n", v, format(get(v)))) else cat(v, "not found in this session\n")
    cat("Methods say: RW2 5 (occurrence) and 10 (intensity); interaction 0.5; BYM2 sd Pr(sigma>1)=0.01; phi median 0.5.\n")
    
    # ---- 3. Prior grid numbers for the sentence 'Spearman >= 0.98' ------------------------------------------------
    hd("3. PRIOR GRID  (Methods: 'Spearman correlation >= 0.98', Supplementary Table S[XX])")
    for (part in c("occurrence", "intensity")) {
      g <- rd(paste0("Table_S_prior_grid_", part, ".csv"))
      if (!is.null(g)) {
        sp <- grep("spearman", names(g), value = TRUE, ignore.case = TRUE)[1]
        t10 <- grep("top10", names(g), value = TRUE, ignore.case = TRUE)[1]
        ph <- grep("^phi", names(g), value = TRUE, ignore.case = TRUE)[1]
        cat(sprintf("%-10s n = %d | Spearman min %.3f (median %.3f)", part, nrow(g), min(g[[sp]], na.rm = TRUE), median(g[[sp]], na.rm = TRUE)))
        if (!is.na(t10)) cat(sprintf(" | top-10 shared %d-%d", min(g[[t10]], na.rm = TRUE), max(g[[t10]], na.rm = TRUE)))
        if (!is.na(ph)) cat(sprintf(" | phi median %.2f (range %.2f-%.2f)", median(g[[ph]], na.rm = TRUE), min(g[[ph]], na.rm = TRUE), max(g[[ph]], na.rm = TRUE)))
        cat("\n")
      }
    }
    cat("The sentence is correct only if both Spearman minima are >= 0.98.\n")
    
    # ---- 4. Second harmonic (the Belg/Kiremt sentence) ----------------------------------------------------------
    hd("4. SECOND HARMONIC  (only claim a bimodal season if a CrI excludes 0)")
    for (nm in c("occurrence", "intensity")) {
      f <- if (nm == "occurrence") fit_o$summary.fixed else fit_i$summary.fixed
      h <- f[c("sin_month", "cos_month", "sin2_month", "cos2_month"), c("mean", "0.025quant", "0.975quant")]
      cat("\n", nm, ":\n", sep = ""); print(round(h, 3))
      s2 <- f[c("sin2_month", "cos2_month"), ]
      cat("  second harmonic:", if (any(s2[, "0.025quant"] > 0 | s2[, "0.975quant"] < 0)) "a CrI excludes 0 -> a second peak is supported" else "both CrIs include 0 -> do NOT write 'bimodal'", "\n")
    }
    
    # ---- 5. Zero-truncated NB sensitivity (rho and the Supplementary table) ----------------------------------
    hd("5. ZERO-TRUNCATED vs ORDINARY NB  (Methods: sensitivity analysis, Supplementary Tables S6c and S6d)")
    z <- rd(paste0("Table_S_zero_truncated_intensity.csv")); if (!is.null(z)) print(z)
    zr <- rd("Table_S_zero_truncated_ranking.csv"); if (!is.null(zr)) { cat("\nAgreement between the two intensity likelihoods (zone effects and zone-level expected size):\n"); print(zr) }
    cat("Use these two rows for the sentence on the ordinary-NB sensitivity analysis ('Spearman rho = [XX]; [XX] of the top 10 zones shared').\n",
        "Step 7C also prints them in the console.\n")
    
    # ---- 6. Association table (Supplementary table for the joint model) ----------------------------------------
    hd("6. ASSOCIATION BY POPULATION ADJUSTMENT  (Methods: joint model, Supplementary Table S[XX])")
    d7 <- rd("Table_S_decoupling_population_adjustment.csv"); if (!is.null(d7)) print(d7)
    cat("Main specification = population in both parts. Write the correlation and 95% CrI from that row.\n")
    
    # ---- 7. Held-out validation and calibration (the model evaluation sentence) ---------------------------------
    hd("7. HELD-OUT VALIDATION AND CALIBRATION  (Methods: 'Model evaluation')")
    ho <- rd("Table_S_overfit_heldout_validation.csv"); if (!is.null(ho)) print(ho)
    pt <- rd("Table_S_PIT_all_models.csv"); if (!is.null(pt)) { cat("\nPIT summary:\n"); print(pt) }
    cx <- rd("Table_S_overfit_complexity.csv"); if (!is.null(cx)) { cat("\nComplexity (p_WAIC) by model:\n"); print(cx) }
    cat("The sentence 'the interaction models did not predict better in the forward period' is true only for Types II and III.\n",
        "Check the Type IV forward AUC against the Baseline before you keep it.\n")
    
    # ---- 8. Moran's I on residuals (the Methods say queen contiguity with island links) -------------------
    hd("8. RESIDUAL MORAN'S I  (Methods: 'Model evaluation', last sentences)")
    cat("These come from the Step 3 pre-deployment table (checks 5 occurrence and intensity). Read them from the console or from\n",
        "Table_S4_predeployment_checks*.csv:\n")
    pc <- list.files(output_path, pattern = "predeploy|pre_deploy|Table_S4", full.names = TRUE)
    if (length(pc)) { x <- fread(pc[1]); print(x[grepl("Moran|autocorr", paste(x[[1]], x[[2]]), ignore.case = TRUE)]) } else cat("(no pre-deployment table found)\n")
    
    cat("\nDone. Every yellow [XX] in the Word file maps to one block above:\n",
        " 1 = data numbers | 2-3 = priors and grid | 4 = bimodal sentence | 5 = rho and its table | 6 = association table | 7 = model evaluation | 8 = Moran's I\n")
    
  })
}   # end PART 2
