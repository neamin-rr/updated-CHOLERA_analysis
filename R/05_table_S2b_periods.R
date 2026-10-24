## ===========================================================================
## Supplementary Table S2b: periods of larger estimated outbreak size
## Builds the table from Table_temporalRR_monthly.csv (national model).
## Edit the two paths, then run.
## ===========================================================================
suppressPackageStartupMessages(library(data.table))

IN  <- file.path(PROJ)   # <-- your file
OUT <- file.path(PROJ)                                # <-- folder for the table
GAP <- 4L     # elevated runs separated by <= GAP months are one period (MERGE_GAP in the main script)

d <- fread(IN)
setnames(d, gsub("^RR \\(national\\)$", "RR", names(d)))
setnames(d, gsub("^95% CrI \\(national\\)$", "CrI", names(d)))
setnames(d, gsub("^Signal \\(national\\)$", "Signal", names(d)))
stopifnot(all(c("Year","Month","RR","CrI","Signal") %in% names(d)))

d[, RR := as.numeric(RR)]
d[, mn := match(Month, month.abb)]
d[, t  := (Year - min(Year)) * 12 + mn]
setorder(d, t)
d[, lo := as.numeric(sub("-.*", "", CrI))]
d[, hi := as.numeric(sub(".*-", "", CrI))]

el <- d[Signal == "credibly elevated"]
if (!nrow(el)) stop("no credibly elevated months found - check the Signal column")
el[, grp := cumsum(c(1L, diff(t)) > GAP)]

tab <- el[, {
  k <- which.max(RR)
  .(From            = sprintf("%s %d", Month[1], Year[1]),
    To              = sprintf("%s %d", Month[.N], Year[.N]),
    `Elevated months` = .N,
    `Span (months)`   = t[.N] - t[1] + 1L,
    `Peak month`      = sprintf("%s %d", Month[k], Year[k]),
    `Peak RR`         = sprintf("%.2f", RR[k]),
    `Peak 95% CrI`    = sprintf("%.2f-%.2f", lo[k], hi[k]),
    `Median RR`       = sprintf("%.2f", median(RR)))
}, by = grp][, grp := NULL][]
tab[, Period := .I]
setcolorder(tab, "Period")

print(tab, row.names = FALSE)
fwrite(tab, file.path(OUT, "Table_S2b_periods_of_larger_outbreak_size.csv"), bom = TRUE)

cat("\n--- sentence for the supplement ---\n")
cat(sprintf("Outbreak size was credibly larger than the 2015-2024 average month in %d separate periods: %s.\n",
            nrow(tab),
            paste(sprintf("%s to %s (%d elevated months; peak %s, RR %s, 95%% CrI %s)",
                          tab$From, tab$To, tab$`Elevated months`, tab$`Peak month`,
                          tab$`Peak RR`, tab$`Peak 95% CrI`), collapse = "; ")))
cat(sprintf("\nMonths credibly elevated in total: %d of %d.\n", nrow(el), nrow(d)))
