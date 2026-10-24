#!/bin/bash
# Live dashboard for the FINAL run.  bash dash_final.sh   (Ctrl+C leaves it; the run keeps going)
cd ~/Documents/spatial_anlysis/update_HMM
while true; do
  clear; echo "=== FINAL HMM RUN   $(date '+%d %b %H:%M:%S') ==="
  echo; echo "--- the four chains (iteration done of total) ---"
  for f in improve_hmm/chain_*_progress.txt; do [ -f "$f" ] || continue
    l=$(tail -1 "$f"); d=$(echo "$l" | awk '{print $(NF-2)}'); t=$(echo "$l" | awk '{print $NF}')
    p=$(awk -v d="$d" -v t="$t" 'BEGIN{ if (t>0) printf "%.1f", 100*d/t }')
    echo "$(echo $f | sed 's/.*_\(ch[0-9]*\).rds_progress.txt/\1/')  $l   (${p}%)"; done
  first=$(ls -t improve_hmm/chain_*_progress.txt 2>/dev/null | tail -1)
  if [ -n "$first" ]; then s=$(head -1 "$first" | awk '{print $1,$2,$3}'); echo; echo "started: $s   now: $(date '+%d %b %H:%M:%S')"; fi
  echo; echo "--- log ---"; tail -3 improve_hmm/final_log.txt 2>/dev/null | cut -c1-150
  echo; echo "--- machine ---"
  echo "R workers: $(ps aux | grep -c '[w]orkRSOCK')"; uptime | sed 's/.*load/load/'
  free -g | awk 'NR==2{print "memory used "$3" of "$2" GB"}'
  sensors 2>/dev/null | grep -E "^Package id" | head -2
  echo; tail -1 guard_log.txt 2>/dev/null
  sleep 20
done
