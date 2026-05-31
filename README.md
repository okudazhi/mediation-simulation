# mediation-package

Last updated: 31May2026

Source-based R code for comparing causal mediation estimators with
time-varying mediators.  The codebase has two main workflows:

- Simulation comparison under controlled longitudinal data-generating mechanisms.
- Real-data empirical comparison (draft).

This directory is not organized as an installable R package with a
`DESCRIPTION` file.  The entry-point scripts source the required files under
`R/` directly and are intended to be run from, or with `chdir = TRUE` into, the
package root.

## Current Structure

- `run_analysis.R`: real-data arm-comparison entry point.
- `run_simulation.R`: serial simulation entry point.
- `run_simulation_parallel.R`: replication-sharded parallel simulation wrapper
  using the settings in `run_simulation.R`.
- `R/`: core estimands, DGPs, nuisance fitting, estimators, real-data
  preparation, simulation orchestration, summaries, and figure builders.
- `scripts/`: focused ltmle_exact diagnostic, acceptance, static-audit, and
  aggregation scripts.
- `tests/`: `testthat` regression and acceptance tests.
- `data/`: input-data placement area and dummy structure-check CSV.
- `output/`: generated outputs and diagnostics.
- `document/`: manuscript, appendix, and abstract artifacts.

The real-data layer is randomized treatment-arm comparison only.  Legacy
targeted real-data APIs have been removed.  The simulation and real-data TMLE
paths both use the full exact-EIF targeted minimum loss-based estimator for the
nine fixed-horizon component means.

The manuscript comparison estimators are:

- MSM/IPW.
- g-computation.
- `ltmle_exact`.

Yamamuro-type fixed-horizon effects and Tai-type joint-draw supplement effects
are both reported from the shared nine-component mean registry.

## Dependencies

Core code uses base R plus optional packages where needed:

- `haven`: required to read the production SAS real-data input.
- `SuperLearner`: required only when the Super Learner nuisance backend is
  requested.
- `dplyr`, `tidyr`, `ggplot2`: required for simulation figure generation.
- `statmod`: optional; used for Gaussian quadrature when available.
- `testthat`: required to run the test suite.
- `rstudioapi`: optional convenience for locating the package root from
  RStudio.

## Real-Data Input

Place the authoritative analysis SAS data set here:

```text
data/medidata.sas7bdat
```

An optional SAS catalog can be placed here:

```text
data/medidata.sas7bcat
```

`run_analysis.R` requires `data/medidata.sas7bdat` and does not automatically
fall back to CSV.  `data/medidata_dummy.csv` is only for explicit structure
checks and dry-run development.

The expected real-data layout is a wide subject-level data set with one row per
subject, visit-specific columns with suffixes `_0`, `_1`, ..., final
outcome , and censoring indicators `censor_0`, ... .
See `data/README_DATA.txt` for the current real-data variable map.

## Running the Real-Data Empirical Comparison

```r
source("/full/path/to/mediation-package/run_analysis.R", chdir = TRUE)
```

Current defaults in `run_analysis.R`:

```r
T <- 13L
estimands <- c("yamamuro", "tai")
methods <- c("msm_ipw", "gcomp", "ltmle_exact")
nuisance_engines <- c("glm")
```

The default output directory is a timestamped run folder:

```text
output/analysis/<timestamp>/
```

The default method folders are:

```text
arm_comparison/msm_ipw/
arm_comparison/gcomp/
arm_comparison/ltmle_glm/
```

If `superlearner` is added to `nuisance_engines`, the ltmle_exact output folder
is `arm_comparison/ltmle_sl/`.  The current default `default_method_specs()`
keeps MSM/IPW and g-computation on the GLM backend and varies only the
ltmle_exact nuisance backend.

Each real-data estimator folder contains:

- `<estimator>_result.rds`.
- `<estimator>_component_means.csv`.
- `<estimator>_week52_fixed.csv` and `.tex`.
- `<estimator>_joint_draw_supplement.csv` and `.tex`.
- `<estimator>_time_resolved.csv`, when time-resolved rows are available.
- estimator-specific diagnostics such as MC diagnostics, targeting
  diagnostics, weight summaries, EIF diagnostics, and generated figures.

The run root also stores `run_config.rds`, `sessionInfo.txt`, and
`realdata_results.rds`.  Within `arm_comparison/`, the prepared analysis
specification and metadata are stored as `analysis_spec.rds` and
`analysis_meta.rds`.

### Real-Data Bootstrap

Subject-level bootstrap is configured in `run_analysis.R`:

```r
bootstrap_cfg <- default_esax_bootstrap_cfg(
  enabled = TRUE,
  B = 1000L,
  conf_level = 0.95,
  seed = 20260531L
)
```

Set `B = 0L` for a quick dry run.  Bootstrap checkpoints are written under each
estimator folder's `bootstrap/` directory and resume by default.

### Real-Data ltmle_exact Adapter

The real-data ltmle_exact adapter calls the same full engine as the simulation
layer through an explicit `node_spec`.  It requires prepared wide data with the
resolved baseline, treatment, mediator, post-mediator covariate, outcome, and
censoring columns.  It fails closed on missing columns and does not fabricate
intermediate outcomes or collapse multiple post-mediator covariate blocks into a
single scalar L.

## Running the Simulation

```r
source("/full/path/to/mediation-package/run_simulation.R", chdir = TRUE)
```

Current serial simulation defaults in `run_simulation.R`:

```r
seed_base <- 20260531
output_dir <- file.path("output", "simulation_8scenarios")
R_reps <- 1000L
B_mc <- 2000L
B_boot <- 0L
T <- 5
B_truth <- 5000000L
truth_n_batches <- 25L
n_dgp <- c(500L)
DGP_treat_mech <- "baseline_rct"
p_rct <- 0.5
```

The default scenario subset contains 8 scenarios:

- pathway magnitude: `high`, `low`.
- residual mediator dependence: `resid-corr-present`.
- mediator interdependence structure: `none`, `full`.
- Q model: `correct`, `wrong`.

The default treatment regimens are always treated versus never treated:

```r
reg_a <- rep(1, T)
reg_as <- rep(0, T)
```

The serial simulation registry currently evaluates:

- `msm_ipw_quantile_truncated`.
- `msm_ipw_untruncated`.
- `gcomp`.
- `ltmle_exact_quantile_truncated`.
- `ltmle_exact_untruncated`.

Core simulation outputs are written directly under `output_dir`:

```text
run_config.rds
sessionInfo.txt
scenario_manifest.csv
truth_worldmeans.csv
truth_effects.csv
truth_figures.pdf
estimates_long.csv
worldmeans_estimates_long.csv
performance_summary.csv
performance_figures.pdf
```

Figure PDFs are controlled by `make_truth_figures` and `make_figures`.

## Running the Parallel Simulation

```r
source("/full/path/to/mediation-package/run_simulation_parallel.R", chdir = TRUE)
```

The parallel wrapper reads and evaluates the USER SETTINGS block from
`run_simulation.R`, then splits replication indices across shards.

Current parallel defaults:

```r
parallel_splits <- 5L
parallel_workers <- NULL
```

When `parallel_workers` is `NULL`, it uses one worker per shard, capped by the
number of shards.  On Windows, the wrapper falls back to one worker because it
uses `parallel::mclapply()`.

Parallel shard artifacts are written under:

```text
output_dir/parallel/shard_manifest.csv
output_dir/parallel/shards/shard_<id>_result.rds
output_dir/parallel/integrated_result.rds
output_dir/diagnostics/parallel_shards/shard_<id>/
```

After all shards complete, the parent process recomputes full truth and
integrated performance so that final CSV outputs match the serial wrapper's
output structure.

## Simulation Diagnostics

Production simulation diagnostics are lightweight by default and are written
under:

```text
output_dir/diagnostics/
```

Possible diagnostics include:

```text
worldmean_bias_by_estimator.csv
worldmean_bias_summary.csv
ltmle_exact_score_equations.csv
ltmle_exact_run_summary.csv
ltmle_exact_fold_bounds.csv
ltmle_exact_component_law_registry.csv
ltmle_exact_factor_tasks.csv
ltmle_exact_component_eif_summary.csv
ltmle_exact_component_eif_terms.csv
truncation_diagnostics.csv
truncation_diagnostics_summary.csv
msm_weight_diagnostics.csv
gcomp_mcse.csv
estimator_attempt_status.csv
estimator_runtime_summary.csv
estimator_failures.csv
full_rds_index.csv
```

`estimator_failures.csv` is written only when estimator failures occur.
`ltmle_exact_score_equations.csv` is written for outlier ltmle_exact runs by
default, or for all ltmle_exact runs if `save_score_equations = "all"`.
Full per-run RDS diagnostics are disabled by default (`save_full = "none"`).

The separate diagnostic and acceptance scripts under `scripts/` are intended
for implementation validation, root/source-graph debugging, and static audits.
They are not part of the default production simulation run.

## Truncation and Bounds

The common truncation framework is quantile-based:

- MSM/IPW quantile-truncates final cumulative weights at the 1st and 99th
  percentiles.
- ltmle_exact quantile-truncates the targeting clever covariate H at the 1st and
  99th percentiles.
- g-computation has no weight-like truncation target and is reported as not
  applicable.

No fixed-bound density-ratio truncation is used.  The legacy
`density_ratio_bounds` argument is rejected by the ltmle_exact public APIs.
Probability bounds, when used, are numerical probability bounds and are not
counted as truncation.  The simulation ltmle_exact default uses training-fold
outcome bounds because the Gaussian DGP has unbounded outcomes.

## Model Specification

Simulation uses `Q_model = "correct"` or `"wrong"` through the scenario
manifest.  The DGP is fixed by the scenario, while the Q-model label controls
which nuisance-regression terms are included.

Real-data analyses use named model sets such as `main` and
`sensitivity_restricted_covariates`; real-data code intentionally avoids
calling model sets correct or wrong.

The code keeps model specification separate from the learner backend:

```text
model_set / Q_model = variables and terms included
nuisance_engine     = fitting backend, e.g. glm or SuperLearner
```

## Tests

Run tests from the package root:

```sh
cd /full/path/to/mediation-package
Rscript tests/testthat.R
```

The tests cover core estimand definitions, real-data preparation behavior,
bootstrap checkpointing, estimator dispatch, truncation policy, ltmle_exact
task-graph semantics, source-row handling, dynamic acceptance gates, and
production diagnostic output schemas.

## Notes

Runtime validation should be performed in the target R environment with the
required packages installed.  Large defaults such as `B_truth = 5000000L` and
`B = 1000L` are intended for manuscript-grade runs, not quick smoke checks.
