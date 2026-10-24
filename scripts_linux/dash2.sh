#!/bin/bash
# Detailed dashboard:  bash ~/Documents/spatial_anlysis/update_HMM/sh/dash2.sh     (Ctrl+C leaves; runs keep going)
B=~/Documents/spatial_anlysis/update_HMM; P=$B/pipeline; . $B/sh/common.sh; mkdir -p $P
TOT=300000; TAG=tau2_2_hi170_pop_rmin0.5_hm; D1=$B/final_pop/improve_hmm; D3=$B/joint_20draws_long_pop
bar() { local p=$1 n=$(( $1 / 5 )); printf "["; printf "%${n}s" | tr ' ' '#'; printf "%$((20-n))s" | tr ' ' '.'; printf "] %3d%%" $p; }
clk() { date -d @$1 '+%a %H:%M' 2>/dev/null; }
dur() { local s=$1; [ $s -lt 0 ] && s=0; printf "%dh %02dm" $((s/3600)) $(((s%3600)/60)); }
while true; do NOW=$(date +%s); W=""
  clear; echo "======================= CHOLERA PIPELINE   $(date '+%a %d %b %H:%M:%S') ======================="
  echo "pipeline: $(pgrep -f 'sh/pipeline.sh run' >/dev/null && echo RUNNING || echo 'NOT running')"
  printf "  steps: "; for k in 0 1 2 3a 3b 4 5; do [ -f $P/state_$k ] && { s=$(cut -d'|' -f2 $P/state_$k | cut -c1-12); printf "%s=%s  " $k "$s"; echo "$s" | grep -q FAIL && W="$W step $k FAILED;"; }; done; echo
  # ---------------- step 1: final HMM -----------------------------------------------------------
  echo; echo "--- STEP 1  final population HMM (4 chains x 300,000) ---"
  st=$(stat -c %Y $D1/hmm_check_model_${TAG}*.jags 2>/dev/null | sort -n | head -1); mn=$TOT; any=0
  for c in 1 2 3 4; do f=$(ls $D1/chain_*ch${c}.rds_progress.txt 2>/dev/null | head -1); [ -f "$f" ] || continue; any=1
    it=$(tail -1 $f | grep -o 'iteration [0-9e+.]*' | awk '{printf "%d", $2}'); [ -z "$it" ] && it=0; [ $it -lt $mn ] && mn=$it
    age=$(( (NOW - $(stat -c %Y $f)) / 60 )); [ $age -gt 60 ] && W="$W chain $c no progress for ${age} min;"
    printf "  chain %d  %7d / %d  %s   (last update %d min ago)\n" $c $it $TOT "$(bar $((100*it/TOT)))" $age; done
  if [ -f $B/final_pop/improve_hmm/Methods_Results_sentences_${TAG}_FINAL.txt ]; then echo "  >>> STEP 1 FINISHED"; E1=$NOW
  elif [ $any = 1 ]; then
    H=$P/rate_history; echo "$NOW $mn" >> $H; old=$(awk -v t=$((NOW-3600)) '$1>=t' $H | head -1)
    [ -n "$old" ] && { t0=${old% *}; i0=${old#* }; } ; rr=0; [ -n "$old" ] && [ $((NOW-t0)) -gt 600 ] && rr=$(( (mn-i0)*3600/(NOW-t0) ))
    ro=0; [ -n "$st" ] && [ $((NOW-st)) -gt 0 ] && ro=$(( mn*3600/(NOW-st) )); r=$rr; [ $r -le 0 ] && r=$ro
    if [ $r -gt 0 ]; then E1=$(( NOW + (TOT-mn)*3600/r )); echo "  speed: last hour $rr/h | overall $ro/h   running for $(dur $((NOW-st)))   ->  expected finish $(clk $E1) (in $(dur $((E1-NOW))))"; fi
    tail -200 $H > $H.tmp && mv $H.tmp $H
  else echo "  not started"; fi
  # ---------------- step 3a: refits -----------------------------------------------------------------
  echo; echo "--- STEP 3a  20 refits with population (40 chains, 2 per refit) ---"
  n=$(ls $D3/draw*_chain*.rds 2>/dev/null | wc -l); s3=$(stat -c %Y $D3/hmm_model_joint.jags 2>/dev/null)
  printf "  chains finished %2d / 40  %s   workers %d\n" $n "$(bar $((100*n/40)))" $(workers_of long | wc -l)
  if [ $n -ge 40 ]; then echo "  >>> REFITS FINISHED ($(clk $(ls -t $D3/draw*_chain*.rds | head -1 | xargs stat -c %Y)))"
  elif [ -n "$s3" ] && [ $n -gt 0 ]; then f1=$(ls -tr $D3/draw*_chain*.rds | head -1 | xargs stat -c %Y); rd=$((f1-s3))
    echo "  first round took $(dur $rd) (finished $(clk $f1))  ->  remaining chains expected about $(clk $((f1+rd)))"
  elif [ -n "$s3" ]; then echo "  first round running for $(dur $((NOW-s3)))"; fi
  [ -f $D3/Joint20_long_sentences.txt ] && echo "  >>> COMPARISON WITH MAIN RUN DONE (step 3b)"
  # ---------------- later steps ---------------------------------------------------------------------
  echo; echo "--- LATER STEPS ---"
  if [ -n "$E1" ]; then echo "  3b comparison  ~$(clk $((E1+2400)))   |   4 tables  ~$(clk $((E1+3600)))   |   5 Stage 3  ~$(clk $((E1+9000)))   |   ALL DONE ~$(clk $((E1+9000)))"; fi
  [ -f $P/RESULTS_SUMMARY.txt ] && echo "  >>> PIPELINE DONE: send pipeline/RESULTS_SUMMARY.txt"
  # ---------------- machine -----------------------------------------------------------------------
  echo; echo "--- MACHINE ---"
  po=$(workers_of pop | wc -l); lo=$(workers_of long | wc -l); ot=$(( $(workers_of orphan | wc -l) + $(workers_of smoke | wc -l) ))
  echo "  R workers: final/hold-out $po | refits $lo | leftover $ot   (cores 64)    load $(cut -d' ' -f1-3 /proc/loadavg)"
  [ $ot -gt 0 ] && W="$W $ot leftover workers;"
  free -g | awk 'NR==2{printf "  memory: %d of %d GB used (%d available)", $3, $2, $7}'; echo "    disk free: $(df -h $B | awk 'NR==2{print $4}')"
  sensors 2>/dev/null | grep -E "^Package id" | head -2 | while read a b c t rest; do v=${t//[+°C]/}; v=${v%.*}; printf "  CPU %s %3d C  " "$c" $v; bar $((v>100?100:v)); [ $v -ge 90 ] && printf "  HOT"; echo; done
  mt=$(maxtemp); [ "${mt:-0}" -ge 90 ] && W="$W CPU ${mt} C;"
  if [ -f $B/sh/guard_state ]; then read gs gt gtime < $B/sh/guard_state; np=$(grep -c PAUSED $B/sh/guard_temp.log 2>/dev/null)
    echo "  temperature guard: $gs (checked $(( NOW - gtime )) s ago; paused $np times today; pauses refits at 90 C, resumes at 80 C)"
    [ $((NOW-gtime)) -gt 300 ] && W="$W guard not running;"; else echo "  temperature guard: not running"; W="$W guard not running;"; fi
  # ---------------- log + warnings -----------------------------------------------------------------
  echo; echo "--- LAST PIPELINE EVENTS ---"; tail -4 $P/pipeline_log.txt 2>/dev/null | sed 's/^/  /' | cut -c1-120
  echo; [ -n "$W" ] && echo "!!! ATTENTION:$W" || echo "OK: no problems detected"
  sleep 30; done
