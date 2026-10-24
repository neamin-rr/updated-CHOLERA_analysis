## ===========================================================================
## Population scaling for BOTH parts of the hurdle model
##   occurrence : odds of an outbreak      (no offset -> odds scale)
##   intensity  : outbreak size            (offset logE -> 1 = proportional)
## Run in RStudio. Edit only the path below.
## ===========================================================================

RDS <- STAGE1_RDS
obj <- readRDS(RDS)

find_it <- function(x, test, depth = 0) {
  if (depth > 4 || is.null(x)) return(NULL)
  if (test(x)) return(x)
  if (is.list(x)) for (e in x) { r <- find_it(e, test, depth + 1); if (!is.null(r)) return(r) }
  NULL
}
is_fit <- function(v) function(x) inherits(x, "inla") && v %in% rownames(x$summary.fixed)
fit_o <- find_it(obj, is_fit("log_pop_s"))        # occurrence
fit_i <- find_it(obj, is_fit("log_pop_zone_s"))   # intensity
dat   <- find_it(obj, function(x) is.data.frame(x) &&
                 all(c("monthly_population", "ID_area") %in% names(x)))
if (is.null(dat)) stop("data frame not found; run str(obj, max.level = 2)")

SD_occ <- sd(log(dat$monthly_population))                                   # all zone-months
SD_int <- sd(tapply(dat$monthly_population, dat$ID_area, function(v) log(mean(v))))

show <- function(b, SD, kind) {
  k <- b / SD                      # per one unit of log population
  cat(sprintf("\n--- %s ---\n", kind))
  cat(sprintf("beta per SD      : %+.3f (95%% CrI %+.3f to %+.3f)\n", b[1], b[2], b[3]))
  cat(sprintf("SD of log pop    : %.3f  (1 SD = %.2fx the population)\n", SD, exp(SD)))
  if (kind == "INTENSITY (outbreak size)") {
    e <- 1 + k
    cat(sprintf("EXPONENT         : %.2f (95%% CrI %.2f to %.2f)  [1 = proportional]\n", e[1], e[2], e[3]))
    for (f in c(2, 5, 10)) cat(sprintf("   %2dx people -> %.2fx cases\n", f, f^e[1]))
  } else {
    cat(sprintf("ODDS EXPONENT    : %.2f (95%% CrI %.2f to %.2f)  [odds ~ population^k]\n", k[1], k[2], k[3]))
    for (f in c(2, 5, 10)) cat(sprintf("   %2dx people -> odds x%.2f\n", f, exp(k[1] * log(f))))
  }
}

g <- function(fit, v) unlist(fit$summary.fixed[v, c("mean", "0.025quant", "0.975quant")])
cat("================ POPULATION SCALING, BOTH PARTS ================\n")
if (!is.null(fit_o)) show(g(fit_o, "log_pop_s"),      SD_occ, "OCCURRENCE (odds of an outbreak)")
if (!is.null(fit_i)) show(g(fit_i, "log_pop_zone_s"), SD_int, "INTENSITY (outbreak size)")
cat("\n================================================================\n")
