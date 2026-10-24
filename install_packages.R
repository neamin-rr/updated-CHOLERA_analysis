## Run once on a new Windows machine. INLA and JAGS need special handling.
pk <- c("data.table","ggplot2","sf","mclust","igraph","coda","rjags","loo",
        "officer","flextable","patchwork","scales","viridisLite","RColorBrewer",
        "readxl","writexl","spdep","units")
new <- pk[!pk %in% rownames(installed.packages())]
if (length(new)) install.packages(new, repos = "https://cloud.r-project.org")

if (!requireNamespace("INLA", quietly = TRUE))
  install.packages("INLA",
    repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"),
    dep = TRUE)

cat("\n--- versions used for the published analysis ---\n")
cat("R 4.6.1 | INLA 26.8.7 | rjags 4.15 | mclust 6.0.1 | igraph 2.3.4 | sf 1.1.3 | ggplot2 4.0.3\n")
cat("\n--- this machine ---\n")
cat(R.version.string, "\n")
for (p in c("INLA","rjags","mclust","igraph","sf","ggplot2","officer"))
  cat(sprintf("%-10s %s\n", p,
      tryCatch(as.character(packageVersion(p)), error = function(e) "NOT INSTALLED")))
cat("\nJAGS executable found: ",
    nzchar(Sys.which("jags")) || file.exists("C:/Program Files/JAGS/JAGS-4.3.2/x64/bin/jags-terminal.exe"), "\n")
