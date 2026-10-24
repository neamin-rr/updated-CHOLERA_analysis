#!/bin/bash
# Live dashboard. Run:  bash dash.sh      (Ctrl+C leaves it; the analysis keeps running)
cd ~/Documents/spatial_anlysis/update_HMM
while true; do
  clear
  echo "=== HMM CHECK DASHBOARD   $(date +%T) ==="
  echo; echo "--- launcher ---"; tail -4 variants_log.txt 2>/dev/null
  echo; echo "--- versions ---"
  for f in improve_hmm/run_log_[A-O]_*.txt; do
    [ -f "$f" ] || continue
    n=$(basename $f .txt | sed 's/run_log_//')
    if   grep -q "^Error" $f;            then s="ERROR (see: tail $f)"
    elif grep -q " SUMMARY " $f;         then s="DONE"
    elif grep -q "THE FOUR INDICATORS" $f; then s="step 4/5: indicators + summary"
    elif grep -q "HOLD-OUT" $f;          then s="step 3/5: hold-out done, indicators next"
    elif grep -q "Fitting: without" $f;  then s="step 2/5: first fit done, fitting without offset"
    elif grep -q "Fitting: with" $f;     then s="step 1/5: first fit running"
    else s="starting"; fi
    printf "%-34s %s\n" "$n" "$s"
  done
  echo; echo "--- saved fits ---"; ls -l --time-style=+%H:%M improve_hmm/cache_*.rds improve_hmm/chain_*.rds 2>/dev/null | awk '{print $6, $5/1e6" MB", $7}'
  echo; echo "--- machine ---"
  echo "R processes: $(ps aux | grep -c '[R]script\|[/]exec/R')"; uptime | sed 's/.*load/load/'
  free -g | awk 'NR==2{print "memory used "$3" of "$2" GB"}'
  sensors 2>/dev/null | grep -E "^Package id" | head -2
  sleep 15
done
