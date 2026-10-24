# Linux helper scripts (optional)

Used to run the long JAGS stages on an EPHI Linux server and watch them from a terminal.
**Not needed on Windows** — `run_all.R` does the same work, and the heavy steps already
save each chain as it finishes, so they are resumable.

| Script | Purpose |
|---|---|
| `pipeline.sh`, `start_all.sh`, `start_final.sh` | launch the long runs detached (`setsid nohup`) |
| `guard.sh` | restart a run if it dies |
| `dash.sh`, `dash2.sh`, `dash_all.sh`, `dash_final.sh` | live progress dashboards from the chain logs |
| `mk.js` | small helper for assembling log output |

The Windows equivalent is simply:

```r
source("run_all.R"); run(c(7:10, 12))
```
