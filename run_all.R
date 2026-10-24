## ===========================================================================
##  run_all.R - end-to-end pipeline, Ethiopia cholera 2015-2024
##  Windows / macOS / Linux.  Open the project in RStudio and run this file.
##
##    source("run_all.R")                 # everything
##    source("run_all.R"); run(1:6)       # Stage 1 only
##    source("run_all.R"); run(13:15)     # Stage 3 only
##
##  Stages 7-12 fit JAGS models and take hours. They are skipped by default
##  if their saved output already exists (set FORCE = TRUE to refit).
## ===========================================================================

source("config.R")
FORCE <- getOption("cholera.force", FALSE)

STEPS <- read.csv(text = "
n,file,what,heavy
1,01_stage1_hurdle_inla.R,Stage 1: hurdle model (INLA) - all tables and Figures 1-2,TRUE
2,02_fig_S3_residuals.R,Figure S3: residual spatial pattern,FALSE
3,03_fig3_temporal_rr.R,Figure 3: temporal rate ratio,FALSE
4,04_fig4_occurrence_size_maps.R,Figure 4: occurrence and size maps,FALSE
4.1,04b_fig4_maps_v1.R,Figure 4 (earlier layout; kept for reference),FALSE
5,05_table_S2b_periods.R,Table S2b: periods of larger outbreak size,FALSE
6,06_population_scaling.R,Population scaling exponents (occurrence and size),FALSE
6.1,06b_population_exponent_size_only.R,Population exponent - size component only,FALSE
7,07_joint_fields_draws.R,Stage 2 input: 20 joint draws of the Stage 1 fields,TRUE
7.1,07b_draw_fields_small_machine.R,Draw the fields on a small machine (laptop),FALSE
8,08_hmm_main_and_joint_refits.R,Stage 2: HMM main fit and 20 joint refits [HOURS],TRUE
8.1,08b_hmm_joint_refits_population.R,Stage 2: joint refits with the population term [HOURS],TRUE
9,09_hmm_sensitivity.R,Stage 2: sensitivity analyses [HOURS],TRUE
10,10_hmm_number_of_states.R,Stage 2: two three and four states compared [HOURS],TRUE
11,11_hmm_checks.R,Stage 2: convergence and model checks,FALSE
11.1,11b_check_persistence_leakage.R,Check: forward prediction leakage in the persistence model,FALSE
11.2,11c_hmm_error_by_state.R,Check: prediction error by state,FALSE
12,12_hmm_variants.R,Stage 2: prior and likelihood variants [HOURS],TRUE
13,13_stage3_lpa_louvain.R,Stage 3: latent profiles Louvain communities stability,FALSE
13.1,13b_stage3_stability_options.R,Stage 3: alternative stability settings,FALSE
13.2,13c_stage3_options_fixed.R,Stage 3: fixed-option variant,FALSE
13.3,13d_class_ari_interval.R,Stage 3: interval for the class ARI (run after step 13),FALSE
14,14_fig5_fig6_stage3.R,Figures 5 and 6: profiles communities central zones,FALSE
14.1,14b_fig_stage3_alt.R,Figure 5 (earlier layout; kept for reference),FALSE
15,15_fig_S4_profile_indicators.R,Figure S4: HMM indicators by profile,FALSE
16,16_build_supplement.R,Assemble the Word supplement from the output tables,FALSE
17,17_pipeline_extra.R,Extra descriptive outputs,FALSE
18,18_run_on_small_machine.R,Reduced settings for a laptop,FALSE
90,90_poster_fig1.R,Poster version of Figure 1,FALSE
91,91_poster_fig2_fig3.R,Poster versions of Figures 2 and 3,FALSE
92,92_conference_figures.R,Conference-sized figures,FALSE
", stringsAsFactors = FALSE)

run <- function(which = MAIN) {
  sel <- STEPS[STEPS$n %in% which, ]
  for (i in seq_len(nrow(sel))) {
    s <- sel[i, ]
    cat(sprintf("\n=== [%s] %s\n    %s\n", s$n, s$file, s$what))
    if (s$heavy && !FORCE) cat("    (heavy step - set options(cholera.force = TRUE) to refit)\n")
    t0 <- Sys.time()
    ok <- tryCatch({ source(file.path("R", s$file), local = new.env()); TRUE },
                   error = function(e) { message("    FAILED: ", conditionMessage(e)); FALSE })
    cat(sprintf("    %s in %.1f min\n", if (ok) "done" else "stopped",
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  invisible(NULL)
}

MAIN <- c(1:16)
cat("\nSteps available:\n"); print(STEPS[, c("n", "what")], row.names = FALSE)
cat("\nRun everything with:  run(STEPS$n)\n",
    "Light steps only:     run(c(2:6, 11, 13:16))\n",
    "One step:             run(13)\n\n")
