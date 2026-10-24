#!/bin/bash
# Temperature guard. Run in its own terminal:  bash guard.sh
# If the hottest processor reaches PAUSE_AT, it PAUSES the R jobs (they keep their progress), and resumes them when it cools to RESUME_AT.
# If it reaches STOP_AT (should not happen) it stops the R jobs.
PAUSE_AT=${PAUSE_AT:-80}; RESUME_AT=${RESUME_AT:-70}; STOP_AT=${STOP_AT:-90}
state=running
while true; do
  t=$(sensors 2>/dev/null | grep -E "^Package id" | awk '{v=$4; gsub(/[+°C]/,"",v); print int(v)}' | sort -n | tail -1)
  [ -z "$t" ] && { echo "$(date +%T) cannot read temperature (is lm-sensors installed?)"; sleep 10; continue; }
  load=$(awk '{print $1}' /proc/loadavg); mem=$(free -g | awk 'NR==2{print $3"/"$2" GB"}')
  if [ "$t" -ge "$STOP_AT" ]; then pkill -u "$USER" -f Run_all_variants; pkill -u "$USER" -f Rscript; pkill -u "$USER" -x R; state=STOPPED
  elif [ "$t" -ge "$PAUSE_AT" ] && [ "$state" = running ]; then pkill -STOP -u "$USER" -f Rscript; pkill -STOP -u "$USER" -x R; state=PAUSED
  elif [ "$t" -le "$RESUME_AT" ] && [ "$state" = PAUSED ]; then pkill -CONT -u "$USER" -f Rscript; pkill -CONT -u "$USER" -x R; state=running; fi
  echo "$(date +%T)  hottest processor ${t}°C | load $load | memory $mem | jobs: $state"
  sleep 5
done
