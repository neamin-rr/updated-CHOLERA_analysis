# Published tables and figures

`supplement_tables/` holds the final tables of the published supplement, and `figures/`
the final figure images. They are included so the supplement can be rebuilt from a fresh
clone with no access to the model objects:

```r
source("config.R"); source("run_all.R"); run(16)
```

`config.R`'s `find_output()` looks here after `outputs/` and your own analysis folder, so
a freshly produced table always takes precedence over the bundled copy.

**Provenance.** Stage 1 and Stage 3 tables were regenerated from
`INLA_Rebuilt_From_Zero_OCCPOP.rds` and the Stage 3 indicator file by `R/01` – `R/15`
(see `build_tables.R` in the release notes). Stage 2 tables (S8–S16) are the published
values from the JAGS runs described in `R/08` – `R/12`; refitting them takes several hours
and requires the HMM draws object, which is not distributed here.

These are **results**, not input data. The underlying surveillance data are not included —
see `data/raw/README.md`.
