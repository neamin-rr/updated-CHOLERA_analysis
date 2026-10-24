## ===========================================================================
##  config.R - the only file you need to edit
##  Every script reads its paths from here. Works on Windows, macOS and Linux.
## ===========================================================================

## 1. Project folder (where data/, outputs/ and R/ live). "." = this repository.
PROJ <- normalizePath(getOption("cholera.proj", "."), winslash = "/", mustWork = FALSE)

## 2. Input data (not committed to GitHub - see README)
DATA_RAW <- file.path(PROJ, "data", "raw")
CASES_CSV <- file.path(DATA_RAW, "cholera_zone_month_2015_2024.csv")   # zone, year, month, cases, population
SHAPEFILE <- file.path(DATA_RAW, "Adm_2", "eth_admbnda_adm2_csa_bofedb_2021.shp")

## 3. Outputs
OUT        <- file.path(PROJ, "outputs")
FIG_DIR    <- file.path(OUT, "figures")
TAB_DIR    <- file.path(OUT, "tables")
MODEL_DIR  <- file.path(OUT, "models")
SUPP_DIR   <- file.path(OUT, "supplement")
for (p in c(OUT, FIG_DIR, TAB_DIR, MODEL_DIR, SUPP_DIR))
  dir.create(p, showWarnings = FALSE, recursive = TRUE)

## 4. Saved model objects produced by earlier stages
STAGE1_RDS <- file.path(MODEL_DIR, "INLA_Rebuilt_From_Zero_OCCPOP.rds")
JOINT_RDS  <- file.path(MODEL_DIR, "joint_fields_20.rds")
DRAWS_RDS  <- file.path(MODEL_DIR, "Draws_updated_indicators_tau2_2_hi170_rmin0.5_hm_FINAL.rds")
REFITS_RDS <- file.path(MODEL_DIR, "Draws_joint20_long_indicators.rds")

## 5. Compute. Windows has no fork(), so every script uses PSOCK clusters.
##    Leave a few cores free; JAGS chains are one per worker.
CORES <- as.integer(getOption("cholera.cores",
           Sys.getenv("CORES", max(1L, parallel::detectCores() - 2L))))

## 6. PDF device. cairo_pdf is often unavailable on Windows; fall back silently.
PDF_DEVICE <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf

save_fig <- function(plot, name, width, height, dpi = 600, png = TRUE) {
  ggplot2::ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), plot,
                  width = width, height = height, device = PDF_DEVICE, bg = "white")
  if (png) ggplot2::ggsave(file.path(FIG_DIR, paste0(name, ".png")), plot,
                  width = width, height = height, dpi = dpi, bg = "white")
  invisible(file.path(FIG_DIR, name))
}

## 7. Options the individual scripts read
options(cholera.proj   = PROJ,
        cholera.cores  = CORES,
        cholera.draws  = DRAWS_RDS,
        cholera.refits = REFITS_RDS,
        cholera.stage1 = STAGE1_RDS)

message(sprintf("config.R loaded | project: %s | cores: %d | cairo: %s",
                PROJ, CORES, capabilities("cairo")))

## 8. Where the analysis outputs actually live.
##    If you keep the repository separate from your working analysis folder, point
##    this at that folder; find_output() then searches it recursively by file name.
##      options(cholera.analysis = "F:/Ub/Paper_one_revision"); source("config.R")
ANALYSIS_DIR <- normalizePath(getOption("cholera.analysis", PROJ),
                              winslash = "/", mustWork = FALSE)

.find_cache <- new.env(parent = emptyenv())

find_output <- function(name, kind = c("table", "figure"), pattern = NULL) {
  kind <- match.arg(kind)
  home <- if (kind == "figure") FIG_DIR else TAB_DIR
  direct <- file.path(home, name)
  if (file.exists(direct)) return(direct)
  ## published tables and figures shipped with the repository
  bundled <- file.path(PROJ, "inst",
                       if (kind == "figure") "figures" else "supplement_tables", name)
  if (file.exists(bundled)) return(bundled)
  if (!is.null(.find_cache$index)) idx <- .find_cache$index else {
    idx <- list.files(ANALYSIS_DIR, recursive = TRUE, full.names = TRUE,
                      all.files = FALSE, no.. = TRUE)
    .find_cache$index <- idx
  }
  hit <- idx[basename(idx) == name]
  if (length(hit)) return(hit[1])
  if (!is.null(pattern)) {                      # fall back to a name pattern
    hit <- idx[grepl(pattern, basename(idx))]
    hit <- hit[!grepl("_OLD_|backup|_old", hit, ignore.case = TRUE)]
    if (length(hit)) return(hit[order(nchar(hit))][1])
  }
  NA_character_
}

message(sprintf("analysis folder: %s", ANALYSIS_DIR))
