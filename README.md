# Cholera dynamics in Ethiopia, 2015–2024

Code for *Disentangling cholera dynamics in Ethiopia: outbreak occurrence, size and persistence, and spatial risk communities.*

National surveillance data for **92 zones over 120 months** (112,906 reported cases; 89.5% of zone-months with none) are analysed in three stages:

1. **Stage 1** — Bayesian spatiotemporal **hurdle model** (INLA): where outbreaks occur and how large they are, with population as a covariate in both components and as the offset of the intensity component, plus a joint extension for the spatial association between them.
2. **Stage 2** — three-state first-order Bayesian **hidden Markov model** (JAGS): states of estimated activity, with two- and four-state models compared.
3. **Stage 3** — **latent profile analysis** on burden and peak reached, and **Louvain community detection** on the zone adjacency network, with stability assessed across 500 posterior draws.

---

## Quick start (Windows)

```r
# 1. open cholera-ethiopia.Rproj in RStudio
source("install_packages.R")   # once
source("config.R")             # paths; edit only this file
source("run_all.R")
run(c(2:6, 11, 13:15))         # light steps, minutes
run(STEPS$n)                   # everything, hours
```

`run_all.R` prints the full step list with what each one does.

## What you need installed

| | Version used | Notes |
|---|---|---|
| R | 4.6.1 | |
| **Rtools** | matching your R | Windows only, needed to build packages — https://cran.r-project.org/bin/windows/Rtools/ |
| **JAGS** | 4.3.2 | install the **Windows binary first**, then `rjags` — https://sourceforge.net/projects/mcmc-jags/ |
| **INLA** | 26.8.7 | not on CRAN; `install_packages.R` adds the INLA repository |
| sf | 1.1.3 | ships with GDAL/GEOS/PROJ on Windows, no separate install |
| mclust | 6.0.1 | |
| igraph | 2.3.4 | |
| ggplot2 | 4.0.3 | |

`install_packages.R` reports what is missing and whether JAGS is on the path.

## Data

Not in this repository. See `data/raw/README.md` for the two required inputs and their expected columns. Boundaries are public (OCHA Ethiopia admin-2); the surveillance data are held by the Ethiopian Public Health Institute and available on reasonable request.

## Layout

```
config.R                  all paths, core count, PDF device   <- the only file to edit
run_all.R                 master runner; prints all 31 steps
install_packages.R        one-time setup

R/01 … R/06b              Stage 1: hurdle model, Figures 3, 4, S3, Table S2b, scaling exponents
R/07 … R/12               Stage 2: JAGS — joint draws, main fit, sensitivity, 2/3/4 states, checks
R/13 … R/15               Stage 3: LPA, Louvain, stability, Figures 5, 6, S4
R/16_build_supplement.R   builds the whole Word supplement (21 tables + 4 figures)
R/17, R/18                extra descriptive outputs; reduced settings for a laptop
R/90 … R/92               poster and conference figure versions

scripts_linux/            optional helpers for running the long stages on a Linux server

data/raw/                 inputs (not committed) — see data/raw/README.md
outputs/figures|tables|models|supplement
```

Scripts with a letter suffix (`04b`, `11c`, `13d`) are variants, checks and earlier layouts,
kept so every result in the paper and supplement can be reproduced. `run_all.R` lists them
all with a one-line description; `run(MAIN)` runs only the main path.

## Runtimes

Measured on a 40-core Linux server; a typical Windows laptop with 8 cores will be several times slower.

| Step | Time |
|---|---|
| Stage 1 (INLA, all models and sensitivity analyses) | ~2 h |
| 07 joint field draws | ~20 min |
| 08 HMM main fit + 20 joint refits | ~6 h |
| 09 sensitivity analyses | ~2 h |
| 10 two/three/four states | ~4 h |
| 12 prior and likelihood variants | ~3 h |
| Stage 3 + all figures | ~10 min |

Heavy steps save their chains as they finish and are **resumable** — rerunning skips completed chains. Each JAGS chain runs in its own worker. `config.R` sets `CORES = detectCores() - 2`; override with `options(cholera.cores = 8)`.

## Windows notes

* Every parallel step uses **PSOCK** clusters (`makeCluster`), which work on Windows. No `fork()`, no `mclapply`.
* `cairo_pdf` is often unavailable on Windows; `config.R` detects this and falls back to the standard `pdf` device. Use `save_fig()` for new figures.
* Paths are built with `file.path()` throughout. Keep the project out of OneDrive-synced folders — file locking during long runs can cause write failures.
* JAGS must be installed **before** `rjags`, or the package will not find the library.

## Model specification

**Stage 1.** Hurdle: Bernoulli occurrence + zero-truncated negative binomial intensity. BYM2 spatial effect, RW2 national trend, iid month effect, Fourier seasonality. Occurrence carries `log_pop_s` (standardised log monthly population); intensity carries `log_pop_zone_s` (standardised log zone population) with `offset(log E)`, where E = population × the pooled national rate.

**Stage 2.** λ = μ · exp(δ) · E/Ē · exp(β_pop x_z); negative binomial with r ~ Gamma(2, 0.5)T(0.5,); Dirichlet transition rows, concentration 10, data-informed means; positive-gap ordering with η23 ~ N(log(170 − 28), 1); 4 chains × 300,000, burn-in 200,000, thinning 50.

**Stage 3.** mclust VVI on burden and peak reached (persistence and entry excluded for collinearity); profiles held fixed via `predict()`; Louvain on the BYM2 adjacency with edge weights = mean of rescaled burden and peak; 1,000 degree-preserving null networks; stability over 500 draws and 20 joint refits; zones below 70% flagged uncertain.

## Rebuilding the supplement from a fresh clone

The repository ships the published tables and figures in `inst/`, so the full supplement
builds with nothing else installed beyond `officer` and `flextable`:

```r
source("config.R")
source("run_all.R"); run(16)
# -> outputs/supplement/Supplementary_Material.docx
```

Freshly produced tables always win: `find_output()` looks in `outputs/`, then your own
analysis folder, then the bundled copies in `inst/`.

## If your outputs live elsewhere

The repository ships with empty `outputs/` folders. If your analysis outputs sit in a
separate working folder, point `cholera.analysis` at it and every script will find its
inputs by file name, wherever they are nested:

```r
options(cholera.analysis = "F:/Ub/Paper_one_revision")
source("config.R")
source("run_all.R"); run(16)
```

## Everything is R

There are no Python or shell dependencies in the analysis path. `R/16_build_supplement.R`
builds the complete Word supplement with `officer` and `flextable`, reading the CSV and
XLSX tables and the figure PNGs that the earlier stages write. Tables that have not been
produced yet appear as labelled placeholders, so the document always builds.

`scripts_linux/` holds optional helpers for running the long JAGS stages detached on a
server; they are not needed on Windows.

## Citing

Tesfaye Gebreyesus N, et al. *Disentangling cholera dynamics in Ethiopia: outbreak occurrence, size and persistence, and spatial risk communities.* (under review)

## Licence

Code released under the MIT Licence (`LICENSE`). The surveillance data are not covered by this licence.
