## ===========================================================================
## 16_build_supplement.R - build the whole Word supplement in R
## Reads the CSV/XLSX tables produced by stages 1-3 and the figure PNGs,
## and writes outputs/supplement/Supplementary_Material.docx
##   source("config.R"); source("R/16_build_supplement.R")
## Needs: officer, flextable, data.table, readxl (for .xlsx tables)
## ===========================================================================
suppressPackageStartupMessages({
  library(officer); library(flextable); library(data.table)
})
has_xl <- requireNamespace("readxl", quietly = TRUE)

TITLE <- paste("Disentangling cholera dynamics in Ethiopia: outbreak occurrence,",
               "size and persistence, and spatial risk communities")

## --- the supplement, in order -------------------------------------------
## type: "sec" section heading | "tab" table | "fig" figure
M <- list(
 list(type="sec", t="Stage 1: hurdle model"),
 list(type="tab", id="Table S1",  f="Table1_by_year.csv",
   cap="Reported cholera cases, zones with cases and national annual incidence, Ethiopia, 2015-2024.",
   note="National incidence is total reported cases divided by the mean national population of the year."),
 list(type="fig", id="Figure S1", f="Fig_S1_annual_incidence.png", pat="annual_incidence.*\\.png$", w=6.5, h=6.5,
   cap="Annual reported cholera incidence per 100,000 population by zone, Ethiopia, 2015-2024.",
   note="Each panel shows one year, with the national annual incidence in the panel title."),
 list(type="tab", id="Table S2",  f="Table_S_seasonal_curve_by_month_OCCPOP.csv",
   cap="Seasonal pattern of outbreak occurrence and outbreak size by calendar month, Baseline hurdle model.",
   note="OR, odds ratio of an outbreak relative to the average calendar month; size ratio, expected outbreak size relative to the average month, given an outbreak. * 95% CrI excludes 1."),
 list(type="fig", id="Figure S2", f="Fig_S2_seasonality.png", pat="seasonality.*\\.png$|Fig_S2.*\\.png$", w=6.5, h=3.2,
   cap="Seasonal pattern of outbreak occurrence and size by calendar month, relative to the average month (log scale).",
   note="Points and bars show posterior means and 95% credible intervals; blue, interval excludes 1. Shaded area: June-September (Kiremt)."),
 list(type="tab", id="Table S2b", f="Table_S2b_periods_of_larger_outbreak_size.csv",
   cap="Periods of larger estimated outbreak size, Baseline hurdle model.",
   note="Rate ratio relative to the average month of 2015-2024, after the seasonal terms and conditional on an outbreak. Periods are runs of consecutive months whose 95% credible interval excluded 1, merged when four months or fewer apart. Peak RR is the highest monthly value in the period."),
 list(type="tab", id="Table S3",  f="TableS3_zone_summary_OCCPOP.csv", land=TRUE,
   cap="Zone-level probability of an outbreak, expected outbreak size and expected burden.",
   note="Posterior means and 95% credible intervals, averaged over the 120 months within each posterior draw."),
 list(type="fig", id="Figure S3", f="Fig_S3_residual_spatial.png", pat="residual_spatial.*\\.png$", w=6.5, h=6.0,
   cap="Residual spatial pattern of the Baseline hurdle model.",
   note="(A, B) zone-level standardized residuals of the occurrence and size components. (C, D) Moran scatterplots against the mean residual of neighbouring zones."),
 list(type="tab", id="Table S4",  f="Table_S_zero_truncated_vs_NB.xlsx", pat="hurdle_vs_zeroinflated|zero_truncated_vs_NB",
   cap="Hurdle model compared with single-likelihood alternatives.",
   note="Lower WAIC and higher sum of log-CPO indicate better fit."),
 list(type="tab", id="Table S5",  f="Table_S_likelihood_and_population_sensitivity.xlsx", pat="population_adjustment_spatial_confounding",
   cap="Sensitivity of the zone effects to the population term.",
   note="The alternative removes the standardized log-population covariate; the intensity offset ln(E) is kept in both models."),
 list(type="tab", id="Table S6",  f="Table_S_occurrence_vs_size_association.xlsx", pat="occurrence_vs_size_association",
   cap="Association between the spatial patterns of outbreak occurrence and outbreak size.",
   note="exp(beta) is the ratio of expected outbreak size for a one-unit higher occurrence spatial effect; 1 = no association. Associations only."),
 list(type="tab", id="Table S7",  f="Table_S_persistence_extension.xlsx", pat="persistence_model|persistence_extension",
   cap="Baseline compared with the persistence extension.",
   note="Associations between consecutive months, not transmission."),

 list(type="sec", t="Stage 2: hidden Markov model"),
 list(type="tab", id="Table S8",  f="Table_S_HMM_convergence.csv", pat="^Table_S_HMM_convergence.*\\.csv$",
   cap="Convergence of the hidden Markov model.",
   note="R-hat, potential scale reduction factor; ESS, effective sample size."),
 list(type="tab", id="Table S9",  f="Table_S_stage2_checklist.csv", pat="stage2_checklist",
   cap="Checks of the hidden Markov model.",
   note="Each row is a pre-specified check with its result."),
 list(type="tab", id="Table S10", f="Table_S_state_means.csv", pat="state_means|HMM_state",
   cap="State means, expected stay and population coefficient.",
   note="Expected stay = 1/(1 - p) from the national transition matrix in each posterior draw."),
 list(type="tab", id="Table S11", f="Table_S_HMM_insample_error.csv", pat="^Table_S_HMM_insample_error.*\\.csv$",
   cap="In-sample prediction error of the hidden Markov model.", note="MAE and RMSE against simple methods."),
 list(type="tab", id="Table S12", f="Table_S_HMM_holdout_prediction.csv", pat="^Table_S_HMM_holdout.*\\.csv$",
   cap="Hold-out prediction of 2024.", note="CRPS and 95% interval coverage against simple baselines."),
 list(type="tab", id="Table S13", f="Table_S_HMM_leftover_associations.csv", pat="^Table_S_HMM_(leftover|indicators_by_zone).*\\.csv$",
   cap="Zone indicators and associations left in the residuals.", note="Associations only."),
 list(type="tab", id="Table S14", f="Table_S_prior_sensitivity_refits_hmm.csv", pat="prior_sensitivity",
   cap="Sensitivity of the hidden Markov model to the high-intensity prior and to population.",
   note="Compared with the main model by rank correlation of the zone indicators and the ten highest-ranked zones."),
 list(type="tab", id="Table S15", f="Table_S_stability_refits.csv", pat="stability_refits",
   cap="Refits on 20 joint draws of the first-stage spatial fields.",
   note="Each refit used two chains with its own convergence check."),
 list(type="tab", id="Table S16", f="Table_S_number_of_states.csv", pat="number_of_states",
   cap="Comparison of models with two, three and four hidden states.",
   note="One-step-ahead WAIC from the forward algorithm, hold-out 2024 CRPS and coverage, and identifiability of the top state."),

 list(type="sec", t="Stage 3: latent profiles and spatial risk communities"),
 list(type="tab", id="Table S17", f="Table_S_LPA_selection.csv", pat="LPA_selection",
   cap="Selection of the number of latent profiles.",
   note="BIC in the mclust sign convention (higher is better), entropy and the smallest profile size."),
 list(type="tab", id="Table S18", f="Table_S_LPA_profile_means.csv", pat="LPA_profile_means",
   cap="Mean indicators by latent profile.", note="Posterior means of the zone indicators."),
 list(type="fig", id="Figure S4", f="Fig_S_profile_indicators.png", pat="profile_indicators.*\\.png$", w=6.5, h=2.6,
   cap="Hidden Markov model indicators by latent profile.",
   note="Each dot is one zone (posterior mean); boxes show the median and interquartile range."),
 list(type="tab", id="Table S19", f="Table_S_community_indicators.csv", pat="community_indicators|indicators_by_community",
   cap="Indicators by spatial risk community and class.",
   note="Medians and ranges of zone posterior means."),
 list(type="tab", id="Table S20", f="Table_S_stability_summary.csv", pat="stability_summary",
   cap="Community structure and stability of profiles, communities and classes.",
   note="ARI, adjusted Rand index (1 = identical grouping). Profiles were held fixed and zones assigned with the fitted mixture model."),
 list(type="tab", id="Table S21", f="Table_S_zone_profile_community_class.csv", pat="zone_profile_community_class", land=TRUE,
   cap="Profile, community, class and stability for each zone.",
   note="Shares are over 500 posterior draws; uncertain = share below 0.70 for profile, community or class. Centrality is eigenvector centrality scaled 0-1 and describes position in the adjacency network only.")
)

## Read a table. Excel sheets written for print often carry a title row above the
## real header, which makes readxl invent "...2", "...3" names; find the header row.
read_any <- function(path) {
  if (!grepl("\\.xlsx?$", path, ignore.case = TRUE))
    return(as.data.frame(fread(path)))
  if (!has_xl) return(NULL)
  raw <- suppressMessages(as.data.frame(
           readxl::read_excel(path, col_names = FALSE, .name_repair = "minimal")))
  if (!nrow(raw)) return(NULL)
  filled <- apply(raw, 1, function(r) {
    v <- trimws(as.character(r)); sum(!is.na(v) & nzchar(v) & v != "NA") })
  hdr <- which(filled >= 2)[1]; if (is.na(hdr)) hdr <- 1L
  nm <- trimws(as.character(unlist(raw[hdr, ], use.names = FALSE)))
  blank <- is.na(nm) | !nzchar(nm) | nm == "NA"
  nm[blank] <- paste0("V", seq_along(nm)[blank])
  df <- raw[-seq_len(hdr), , drop = FALSE]
  names(df) <- make.unique(nm)
  keep <- apply(df, 1, function(r) {
    v <- trimws(as.character(r)); any(!is.na(v) & nzchar(v) & v != "NA") })
  df <- df[keep, , drop = FALSE]
  rownames(df) <- NULL
  df[] <- lapply(df, function(x) { x <- as.character(x); x[is.na(x)] <- ""; x })
  df
}

mk <- function(df) {
  ft <- flextable(df)
  ft <- fontsize(ft, size = 8, part = "all")
  ft <- font(ft, fontname = "Times New Roman", part = "all")
  ft <- bold(ft, part = "header"); ft <- theme_booktabs(ft)
  autofit(ft)
}

doc <- read_docx()
doc <- body_add_par(doc, "Supplementary Material", style = "heading 1")
doc <- body_add_par(doc, TITLE, style = "Normal")
doc <- body_add_par(doc, "", style = "Normal")
doc <- body_add_par(doc, "Contents", style = "heading 2")
for (m in M) {
  if (m$type == "sec") doc <- body_add_par(doc, m$t, style = "Normal")
  else doc <- body_add_par(doc, paste0("   ", m$id, ". ", m$cap), style = "Normal")
}
doc <- body_add_break(doc)

missing <- character()
for (m in M) {
  if (m$type == "sec") { doc <- body_add_par(doc, m$t, style = "heading 2"); next }
  src <- find_output(m$f, if (m$type == "fig") "figure" else "table", m$pat)
  doc <- body_add_par(doc, paste0(m$id, ": ", m$cap), style = "heading 3")
  ok <- !is.na(src) && file.exists(src)
  if (ok && m$type == "tab") {
    df <- read_any(src)
    if (is.null(df)) { ok <- FALSE
    } else {
      if (isTRUE(m$land)) doc <- body_end_section_continuous(doc)
      doc <- body_add_flextable(doc, mk(df))
      if (isTRUE(m$land)) doc <- body_end_section_landscape(doc)
    }
  } else if (ok) {
    doc <- body_add_img(doc, src, width = m$w, height = m$h)
  }
  if (!ok) {
    missing <- c(missing, m$f)
    doc <- body_add_par(doc, paste0("[not found: ", m$f, " - run the stage that produces it]"),
                        style = "Normal")
  }
  doc <- body_add_par(doc, paste("Notes:", m$note), style = "Normal")
  doc <- body_add_par(doc, "", style = "Normal")
}

out <- file.path(SUPP_DIR, "Supplementary_Material.docx")
print(doc, target = out)
cat("Written:", out, "\n")
if (length(missing)) {
  cat("\nMissing sources (placeholders inserted) -", length(missing), "of", sum(sapply(M, function(x) x$type != "sec")), ":\n")
  cat(paste0("  - ", missing, collapse = "\n"), "\n")
  if (!has_xl) cat("\n  (install 'readxl' to include the .xlsx tables)\n")
}
