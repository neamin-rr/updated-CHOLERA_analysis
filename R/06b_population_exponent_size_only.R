## ===========================================================================
## Population exponent of outbreak size  (1 = proportional to population)
## Run in RStudio on your PC.  Only edit the path on the next line.
## ===========================================================================

RDS <- STAGE1_RDS   # <-- your file

obj <- readRDS(RDS)

## ---- 1. find the intensity fit and the data inside the object -------------
find_it <- function(x, want, depth = 0) {            # recursive search
  if (depth > 4 || is.null(x)) return(NULL)
  if (want == "fit"  && inherits(x, "inla") &&
      "log_pop_zone_s" %in% rownames(x$summary.fixed)) return(x)
  if (want == "data" && is.data.frame(x) &&
      all(c("monthly_population", "ID_area") %in% names(x))) return(x)
  if (is.list(x)) for (e in x) { r <- find_it(e, want, depth + 1); if (!is.null(r)) return(r) }
  NULL
}

fit <- find_it(obj, "fit")
dat <- find_it(obj, "data")

if (is.null(fit)) stop("No intensity fit with 'log_pop_zone_s' found. Run str(obj, max.level = 2) and tell me the structure.")
if (is.null(dat)) stop("No data frame with monthly_population + ID_area found. Run str(obj, max.level = 2).")

## ---- 2. POP_SD, exactly as the main script computed it --------------------
## settings may already hold it; otherwise rebuild from the data
POP_SD <- NULL
if (is.list(obj) && !is.null(obj$settings$POP_SD)) POP_SD <- obj$settings$POP_SD
if (is.null(POP_SD)) {
  lp     <- tapply(dat$monthly_population, dat$ID_area, function(v) log(mean(v)))
  POP_SD <- sd(lp)
}

## ---- 3. the exponent ------------------------------------------------------
bp <- fit$summary.fixed["log_pop_zone_s", ]
e  <- 1 + c(mean = bp[["mean"]],
            lo   = bp[["0.025quant"]],
            hi   = bp[["0.975quant"]]) / POP_SD

cat("\n================= POPULATION SCALING OF OUTBREAK SIZE =================\n")
cat(sprintf("beta (log_pop_zone_s) : %.3f  (95%% CrI %.3f to %.3f)\n",
            bp[["mean"]], bp[["0.025quant"]], bp[["0.975quant"]]))
cat(sprintf("size ratio per SD     : %.2f  (95%% CrI %.2f to %.2f)\n",
            exp(bp[["mean"]]), exp(bp[["0.025quant"]]), exp(bp[["0.975quant"]])))
cat(sprintf("POP_SD (sd of log pop): %.3f   -> 1 SD = %.2f x the population\n",
            POP_SD, exp(POP_SD)))
cat(sprintf("\nPOPULATION EXPONENT   : %.2f  (95%% CrI %.2f to %.2f)   [1 = proportional]\n",
            e[1], e[2], e[3]))

## ---- 4. what it means in practice ----------------------------------------
cat("\n--- expected cases and cases per person, relative to a reference zone ---\n")
for (f in c(2, 3, 5, 10))
  cat(sprintf("  %2dx the people  ->  %.2fx the cases   (%.2fx the cases per person)\n",
              f, f^e[1], f^e[1] / f))

cat("\nSentence for the paper:\n")
cat(sprintf("  Outbreak size scaled almost flatly with population (exponent %.2f,\n  95%% CrI %.2f to %.2f, where 1 denotes proportional scaling).\n",
            e[1], e[2], e[3]))
cat("=======================================================================\n")
