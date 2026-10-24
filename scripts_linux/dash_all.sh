#!/bin/bash
# Dashboard for all runs:  bash ~/Documents/spatial_anlysis/update_HMM/sh/dash_all.sh   (Ctrl+C leaves; runs keep going)
B=~/Documents/spatial_anlysis/update_HMM; . $B/sh/common.sh
ago() { [ -f "$1" ] && echo "$(( ($(date +%s) - $(stat -c %Y "$1")) / 60 )) min ago" || echo "-"; }
while true; do clear; echo "==================== ALL RUNS   $(date '+%d %b %H:%M:%S') ===================="
  echo; echo "--- PIPELINE (new main model with population) ---   $(pgrep -f 'sh/pipeline.sh run' >/dev/null && echo RUNNING || echo 'not running')"
  for k in 0 1 2 3a 3b 4 5; do [ -f $B/pipeline/state_$k ] && IFS='|' read a b c < $B/pipeline/state_$k && printf "  step %-3s %-28s (%s)\n" "$a" "$b" "$c"; done
  [ -d $B/final_pop/improve_hmm ] && for f in $B/final_pop/improve_hmm/chain_*_progress.txt; do [ -f "$f" ] && echo "  final chain $(basename $f | grep -o 'ch[0-9]*'): $(tail -1 $f | cut -c1-60)"; done
  tail -2 $B/pipeline/pipeline_log.txt 2>/dev/null | sed 's/^/  log: /' | cut -c1-130
  echo; echo "--- workers by run: population runs $(workers_of pop | wc -l) | refits $(workers_of long | wc -l) ---"
  echo; echo "--- (earlier) population short test ---"
  f=$B/test1_population/pop_log.txt; if [ -f $f ]; then echo "log changed $(ago $f)"; grep -E "^Fitting|Loaded|population \(per|Calibration" $f | tail -3 | cut -c1-130
    [ -f $B/test1_population/improve_hmm/Table_S_HMM_small_check_tau2_2_hi170_pop_rmin0.5_hm.csv ] && echo ">>> FINISHED: check table written"; else echo "not started"; fi
  echo; echo "--- Step 2: long 20 refits (joint_20draws_long) ---   workers: $(workers_of long | wc -l)"
  echo "chains finished: $(ls $B/joint_20draws_long/draw*_chain*.rds 2>/dev/null | wc -l) of 40"; tail -1 $B/joint_20draws_long/log.txt 2>/dev/null | cut -c1-130
  [ -f $B/sh/guard_state ] && { read gs gt gtime < $B/sh/guard_state; echo "temperature guard: $gs (last check ${gt} C, $(( $(date +%s) - gtime )) s ago)"; } || echo "temperature guard: not running"
  echo; echo "--- Step 4: Stage 3 (test3_stage3) ---"; [ -f $B/test3_stage3/stage3_log.txt ] && tail -2 $B/test3_stage3/stage3_log.txt | cut -c1-130 || echo "not started"
  echo; echo "--- machine ---"; o=$(( $(workers_of orphan | wc -l) + $(workers_of smoke | wc -l) )); echo "R workers total: $(pgrep -fc workRSOCK)   leftover/smoke: $o"
  uptime | sed 's/.*load/load/'; free -g | awk 'NR==2{print "memory used "$3" of "$2" GB (available "$7")"}'
  sensors 2>/dev/null | grep -E "^Package id" | head -2; tail -1 $B/sh/guard_temp.log 2>/dev/null
  sleep 20; done
