#!/bin/bash
# FINAL HMM run (version F settings): 4 chains x 300,000 iterations (200,000 burn-in, thinning 50), memory capped at 40 GB.
#   bash start_final.sh test    -> 15-minute end-to-end test of the long-run machinery first (recommended)
#   bash start_final.sh         -> the real run (about 21 hours; if the PC stops, run the same command again and it resumes)
MODE=${1:-final}; MEMGB=${MEMGB:-40}
cd ~/Documents/spatial_anlysis/update_HMM || { echo "folder not found"; exit 1; }
for f in HMM_check_ALL.R guard.sh dash_final.sh INLA_Rebuilt_From_Zero_OCCPOP.rds improve_hmm/cache_tau2_2_hi170_rmin0.5_hm_nooff_n1091_b2000.rds; do [ -f "$f" ] || { echo "missing file: $f"; exit 1; }; done
grep -q "cholera.final" HMM_check_ALL.R || { echo "HMM_check_ALL.R is the old version (no final mode). Copy the new one into this folder."; exit 1; }
if [ "$(ps aux | grep -c '[w]orkRSOCK\|[R]un_all_variants')" != "0" ]; then echo "A run is already going. Stop it first (pkill -u $USER -f workRSOCK; pkill -u $USER -f Rscript)."; exit 1; fi
hot=$(sensors | grep -E "^Package id" | awk '{v=$4; gsub(/[+°C]/,"",v); print int(v)}' | sort -n | tail -1)
[ "$hot" -ge ${MAXT:-66} ] && { echo "Processor is at ${hot}°C. Wait until both packages are below 65°C, then run again."; exit 1; }
echo 1 | sudo tee /sys/devices/system/cpu/intel_pstate/no_turbo >/dev/null && echo "turbo boost off"
FT=FALSE; [ "$MODE" = "test" ] && FT=TRUE
OPTS="cholera.tau23 = c(2, 2), cholera.anchhi = 170, cholera.hiermean = TRUE, cholera.rmin = 0.5, cholera.cores = 4, cholera.nch = 4, cholera.holdout = FALSE, cholera.final = TRUE, cholera.finaltest = $FT"
if systemd-run --user --scope -p MemoryMax=${MEMGB}G true >/dev/null 2>&1; then CAP="systemd-run --user --scope -p MemoryMax=${MEMGB}G"; else CAP=""; echo "(memory cap not available here: running without it)"; fi
[ "$CPUSET" = "auto" ] && CPUSET=$(lscpu -p=CPU,CORE,SOCKET | awk -F, '!/^#/ && $3==0 && !seen[$2]++ {print $1}' | head -8 | paste -sd,)   # CPUSET=auto: one thread on each of 8 physical cores of the first processor
PIN=""; [ -n "$CPUSET" ] && { PIN="taskset -c $CPUSET"; echo "Pinned to CPUs $CPUSET"; }   # e.g. CPUSET=0-15 keeps the run on the first processor, which runs cooler
ls improve_hmm/chain_*_state.rds >/dev/null 2>&1 && echo "Saved chain states found: the run RESUMES where it stopped."
echo "Starting the $MODE run: 4 chains, memory cap ${MEMGB} GB. Log: improve_hmm/final_log.txt"
setsid nohup nice -n 19 $PIN $CAP env R_DATATABLE_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript -e "options($OPTS); source('HMM_check_ALL.R')" > improve_hmm/final_log.txt 2>&1 &
pgrep -f guard.sh >/dev/null || setsid nohup bash guard.sh > guard_log.txt 2>&1 &
sleep 5
[ "$BG" = "1" ] && { echo "Running in the background (safe to close this terminal or VS Code). Look at it any time with: bash dash_final.sh"; exit 0; }
exec bash dash_final.sh
