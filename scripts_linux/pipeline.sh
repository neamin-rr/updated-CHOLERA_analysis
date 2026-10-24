#!/bin/bash
# ONE PIPELINE for the new main HMM (with population). Each step is skipped if already finished, waited for if already running, started otherwise.
#   start:  bash ~/Documents/spatial_anlysis/update_HMM/sh/pipeline.sh       (it detaches itself; safe to close the terminal)
#   watch:  bash ~/Documents/spatial_anlysis/update_HMM/sh/dash_all.sh
#   stop:   bash ~/Documents/spatial_anlysis/update_HMM/sh/pipeline.sh stop  (stops the pipeline only; running R jobs keep going)
B=~/Documents/spatial_anlysis/update_HMM; P=$B/pipeline; mkdir -p $P; cd $B
if [ "$1" = stop ]; then pkill -f "sh/pipeline.sh run"; echo "pipeline stopped (R jobs keep running; start again to continue)"; exit 0; fi
if [ "$1" != run ]; then pgrep -f "sh/pipeline.sh run" >/dev/null && { echo "pipeline already running - watch with: bash sh/dash_all.sh"; exit 0; }
  setsid nohup bash $B/sh/pipeline.sh run > $P/pipeline_stdout.txt 2>&1 < /dev/null & echo "pipeline started - watch with: bash $B/sh/dash_all.sh"; exit 0; fi
. $B/sh/common.sh
T0=$(date +%s); log() { echo "$(date '+%d %b %H:%M:%S') $*" >> $P/pipeline_log.txt; }
state() { echo "$1|$2|$(date '+%d %b %H:%M')" > $P/state_$1; }       # step|status|time
hrs() { local s=$(( $(date +%s) - $1 )); printf "%dh %02dm" $((s/3600)) $(( (s%3600)/60 )); }
TAG=tau2_2_hi170_pop_rmin0.5_hm
OPT="cholera.tau23 = c(2, 2), cholera.anchhi = 170, cholera.hiermean = TRUE, cholera.rmin = 0.5, cholera.pop = TRUE"
F1=$B/final_pop/improve_hmm/Methods_Results_sentences_${TAG}_FINAL.txt
F2=$B/test1_population/improve_hmm/Table_S_HMM_holdout_prediction_${TAG}.csv
F3=$B/joint_20draws_long_pop/Joint20_long_sentences.txt
F4=$B/final_pop/improve_hmm/Table_S_HMM_error_by_state_pop_FINAL.csv
F5=$B/test3_stage3_pop/stage2_final/Stage2_sentences.txt
nch3() { ls $B/joint_20draws_long_pop/draw*_chain*.rds 2>/dev/null | wc -l; }
log "===== PIPELINE START ====="
for f in HMM_check_ALL.R INLA_Rebuilt_From_Zero_OCCPOP.rds joint_fields_20.rds Joint_20draws_long_pop.R Stage2_LPA_Louvain_final.R test1_population/improve_hmm/cache_${TAG}_nooff_n1091_b2000.rds; do
  [ -e $B/$f ] || { log "MISSING FILE: $f - pipeline not started"; state 0 "FAILED: missing $f"; exit 1; }; done
grep -q keeping_degseq $B/Stage2_LPA_Louvain_final.R || { log "Stage2_LPA_Louvain_final.R is the OLD version (download the updated one) - pipeline not started"; state 0 "FAILED: old Stage 3 script"; exit 1; }
for p in rjags coda data.table mclust igraph sf ggplot2; do Rscript -e "q(status = !requireNamespace('$p', quietly = TRUE))" >/dev/null 2>&1 || { log "R package missing: $p"; state 0 "FAILED: R package $p"; exit 1; }; done
state 0 "checks passed"
pgrep -f "sh/guard_temp.sh" >/dev/null || { setsid nohup bash $B/sh/guard_temp.sh > $B/sh/guard_temp.log 2>&1 < /dev/null & log "temperature guard started (pauses the refits at 90 C, resumes at 80 C)"; }
echo 1 | sudo -n tee /sys/devices/system/cpu/intel_pstate/no_turbo >/dev/null 2>&1 && log "turbo off"

## step 1: final population HMM (4 chains x 300,000) ---------------------------------------------------------------
s1() { [ -f $F1 ] && { state 1 done; return; }
  if ! pgrep -f "cholera.final = TRUE" >/dev/null; then mkdir -p final_pop/improve_hmm
    ln -sf $B/INLA_Rebuilt_From_Zero_OCCPOP.rds $B/HMM_check_ALL.R final_pop/
    ln -sf $B/test1_population/improve_hmm/cache_${TAG}_nooff_n1091_b2000.rds final_pop/improve_hmm/
    (cd final_pop && setsid nohup nice -n 19 env R_DATATABLE_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript -e "options(cholera.proj = getwd(), $OPT, cholera.cores = 4, cholera.nch = 4, cholera.holdout = FALSE, cholera.final = TRUE, cholera.finaltest = FALSE); source('HMM_check_ALL.R')" >> final_log.txt 2>&1 < /dev/null &)
    log "step 1 started (final population HMM)"; sleep 60; fi
  state 1 running; }
## step 2: 2024 hold-out of the population model (short run, 11 chains) --------------------------------------------
s2() { [ -f $F2 ] && { state 2 done; return; }
  pgrep -f "cholera.holdout = TRUE" >/dev/null || { (cd test1_population && setsid nohup nice -n 19 Rscript -e "options(cholera.proj = getwd(), $OPT, cholera.nch = 11, cholera.cores = 11, cholera.holdout = TRUE, cholera.final = FALSE); source('HMM_check_ALL.R')" >> holdout_log.txt 2>&1 < /dev/null &); log "step 2 started (hold-out)"; sleep 30; }
  state 2 running; }
## step 3a: 20 refits with population (2 chains each), can run while step 1 runs ----------------------------------
s3a() { [ "$(nch3)" -ge 40 ] && { state 3a done; return; }
  pgrep -f "Joint_20draws_long_pop" >/dev/null || { mkdir -p joint_20draws_long_pop; CORES=${CORES3:-20} setsid nohup nice -n 19 Rscript Joint_20draws_long_pop.R >> joint_20draws_long_pop/log.txt 2>&1 < /dev/null & log "step 3a started (20 refits, ${CORES3:-20} workers)"; sleep 30; }
  state 3a "running ($(nch3)/40 chains)"; }

S1T=$(date +%s); fails2=0; fails3=0
while true; do
  s1; s2; s3a
  [ -f $F1 ] && [ "$(cat $P/state_1 2>/dev/null | cut -d'|' -f2)" != done ] && { state 1 done; log "step 1 FINISHED"; }
  [ -f $F2 ] && [ -z "$done2" ] && { done2=1; log "step 2 FINISHED"; }
  [ "$(nch3)" -ge 40 ] && [ -z "$done3" ] && { done3=1; log "step 3a FINISHED (40 chains)"; }
  # restart a step that stopped without finishing (at most twice)
  if [ ! -f $F2 ] && ! pgrep -f "cholera.holdout = TRUE" >/dev/null; then fails2=$((fails2+1)); log "step 2 stopped early (attempt $fails2): see test1_population/holdout_log.txt"; [ $fails2 -ge 3 ] && { state 2 FAILED; touch $F2.skip; F2=$F2.skip; }; fi
  if [ "$(nch3)" -lt 40 ] && ! pgrep -f "Joint_20draws_long_pop" >/dev/null; then fails3=$((fails3+1)); log "step 3a stopped early (attempt $fails3): see joint_20draws_long_pop/log.txt"; [ $fails3 -ge 3 ] && { state 3a FAILED; log "PIPELINE STOPPED at step 3a"; exit 1; }; fi
  if [ ! -f $F1 ] && ! pgrep -f "cholera.final = TRUE" >/dev/null && [ -f $B/final_pop/final_log.txt ] && grep -q "Execution halted" $B/final_pop/final_log.txt; then
    K=tau2_2_hi170_pop_rmin0.5_hm_main_n100000_b200000; D=$B/final_pop/improve_hmm
    if [ $(ls $D/chain_${K}_ch[1-4].rds 2>/dev/null | wc -l) -eq 4 ] && [ ! -f $D/cache_$K.rds ] && [ -z "$combined" ]; then combined=1
      log "step 1 stopped after the 4 chains finished (as on 13 Oct): combining the chains and restarting the checks"
      (cd $D && Rscript -e "library(coda); f <- sprintf('chain_${K}_ch%d.rds', 1:4); x <- lapply(f, readRDS); saveRDS(mcmc.list(lapply(x, as.mcmc)), 'cache_${K}.rds.tmp'); file.rename('cache_${K}.rds.tmp', 'cache_${K}.rds')" >> $P/combine_log.txt 2>&1)
      mv $B/final_pop/final_log.txt $B/final_pop/final_log_first_attempt.txt; continue
    fi
    state 1 FAILED; log "STEP 1 FAILED: see final_pop/final_log.txt (tail -40). Pipeline stopped."; exit 1; fi
  [ -f $F1 ] && [ -f $F2 ] && [ "$(nch3)" -ge 40 ] && break
  sleep 120
done
log "steps 1, 2 and 3a finished after $(hrs $T0)"

## step 3b: comparison of the refits with the final run (reloads the 40 chains) -------------------------------------
if [ ! -f $F3 ]; then state 3b running; log "step 3b started (comparison + main draws file for Stage 3)"
  nice -n 19 Rscript Joint_20draws_long_pop.R >> joint_20draws_long_pop/log.txt 2>&1 < /dev/null
  [ -f $F3 ] && { state 3b done; log "step 3b FINISHED"; } || { state 3b FAILED; log "STEP 3b FAILED: tail -40 joint_20draws_long_pop/log.txt. Pipeline stopped."; exit 1; }
else state 3b done; fi
## step 4: Table 4, state means, stays, error by state (population model) -------------------------------------------
if [ ! -f $F4 ]; then state 4 running; log "step 4 started (Table 4, stays, error by state)"
  nice -n 19 Rscript $B/sh/pipeline_extra.R > $P/step4_output.txt 2>&1
  [ -f $F4 ] && { state 4 done; log "step 4 FINISHED"; } || { state 4 FAILED; log "step 4 FAILED: see pipeline/step4_output.txt (continuing)"; }
else state 4 done; fi
## step 5: Stage 3 (LPA, Louvain, stability over 500 draws and 20 refits) ------------------------------------------
if [ ! -f $F5 ]; then state 5 running; log "step 5 started (Stage 3)"; mkdir -p test3_stage3_pop
  nice -n 19 Rscript -e "options(cholera.out = '$B/test3_stage3_pop', cholera.draws = '$B/final_pop/improve_hmm/Draws_updated_indicators_${TAG}_FINAL.rds', cholera.refits = '$B/joint_20draws_long_pop/Draws_joint20_long_indicators.rds', cholera.ndraws = 500); source('Stage2_LPA_Louvain_final.R')" > test3_stage3_pop/stage3_log.txt 2>&1 < /dev/null
  [ -f $F5 ] && { state 5 done; log "step 5 FINISHED"; } || { state 5 FAILED; log "STEP 5 FAILED: tail -40 test3_stage3_pop/stage3_log.txt"; }
else state 5 done; fi

## summary of everything to send ----------------------------------------------------------------------------------
S=$P/RESULTS_SUMMARY.txt
{ echo "RESULTS SUMMARY  $(date '+%d %b %Y %H:%M')   (pipeline run time $(hrs $T0))"
  for f in $F1 $B/final_pop/improve_hmm/HMM_check_summary_${TAG}_FINAL.txt $P/step4_output.txt; do echo; echo "=== $(basename $f)"; cat $f 2>/dev/null; done
  for f in Table_S_HMM_small_check Table_S_HMM_convergence Table_S_HMM_insample_error Table_S_HMM_updated_indicator_summary Table_S_HMM_leftover_associations; do
    echo; echo "=== ${f}_FINAL"; cat $B/final_pop/improve_hmm/${f}_${TAG}_FINAL.csv 2>/dev/null; done
  echo; echo "=== hold-out (population model)"; cat $B/test1_population/improve_hmm/Table_S_HMM_holdout_prediction_${TAG}.csv 2>/dev/null
  grep -E "^(9|10|11)\." $B/test1_population/improve_hmm/Table_S_HMM_small_check_${TAG}.csv 2>/dev/null
  echo; echo "=== 20 refits"; cat $F3 $B/joint_20draws_long_pop/Table_S_joint20_long_comparison.csv 2>/dev/null
  echo; echo "=== Stage 3"; cat $F5 2>/dev/null; for f in Table_S_LPA_selection Table_S_LPA_profile_means Table_S_stability_summary Table_S_stability_split_half Table_S_stability_refits Table_S_stage2_checklist; do echo "--- $f"; cat $B/test3_stage3_pop/stage2_final/$f.csv 2>/dev/null; done
  echo; echo "=== pipeline log"; cat $P/pipeline_log.txt; } > $S
log "===== PIPELINE DONE after $(hrs $T0). Send: cat $S ====="
