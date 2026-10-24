#!/bin/bash
# ONE command to start everything:   bash start_all.sh
# 3 versions at a time x 11 chains = 33 cores, memory capped at 50 GB, turbo off, temperature guard on, live dashboard.
# Change with:  CHAINS=5 NCONC=1 MEMGB=50 ONLY=D,F,A bash start_all.sh
CHAINS=${CHAINS:-11}; NCONC=${NCONC:-3}; MEMGB=${MEMGB:-50}; ONLY=${ONLY:-D,F,A}
cd ~/Documents/spatial_anlysis/update_HMM || { echo "folder not found"; exit 1; }
for f in HMM_check_ALL.R Run_all_variants.R dash.sh guard.sh INLA_Rebuilt_From_Zero_OCCPOP.rds; do [ -f $f ] || { echo "missing file: $f"; exit 1; }; done
if [ "$(ps aux | grep -c '[w]orkRSOCK\|[R]un_all_variants')" != "0" ]; then echo "A run is already going. Stop it first (pkill -u $USER -f Run_all_variants; pkill -u $USER -f workRSOCK)."; exit 1; fi
hot=$(sensors | grep -E "^Package id" | awk '{v=$4; gsub(/[+°C]/,"",v); print int(v)}' | sort -n | tail -1)
[ "$hot" -ge 66 ] && { echo "Processor is at ${hot}°C. Wait until both packages are below 65°C, then run again."; exit 1; }
echo 1 | sudo tee /sys/devices/system/cpu/intel_pstate/no_turbo >/dev/null && echo "turbo boost off"
rm -f improve_hmm/cache_*.tmp improve_hmm/chain_*.tmp
export ONLY SKIPTEST=1 CHAINS NCONC CORES_TOTAL=$((CHAINS*NCONC))
echo "Starting: $NCONC versions at a time x $CHAINS chains = $((CHAINS*NCONC)) cores, memory cap ${MEMGB} GB"
if systemd-run --user --scope -p MemoryMax=${MEMGB}G true >/dev/null 2>&1; then CAP="systemd-run --user --scope -p MemoryMax=${MEMGB}G"; else CAP=""; echo "(memory cap not available on this system: running without it; memory use was only 10 GB before)"; fi
nohup nice -n 19 $CAP Rscript Run_all_variants.R > variants_log.txt 2>&1 &
nohup bash guard.sh > guard_log.txt 2>&1 &
sleep 5; echo "Started. Guard pauses the jobs at 80°C and resumes at 70°C (guard_log.txt)."; sleep 3
exec bash dash.sh
