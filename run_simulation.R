################################################################################
# run_simulation.R
#
# MAIN ENTRY POINT:
# 1) Edit the USER PARAMETERS section.
# 2) Safest run: source("/full/path/to/run_simulation.R", chdir = TRUE)
#
# Outputs will be written under the output_dir set in USER SETTINGS.
################################################################################

# ---- Robust project root -----------------------------------------------------
# We resolve the package root from multiple sources so the script works when
# called via Rscript, source(), or the RStudio Source button.
.find_pkg_dir <- function(required_rel_paths,
                          script_name = "run_simulation.R",
                          max_up = 4L) {
  .is_scalar_string <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  }
  .norm <- function(p) normalizePath(p, winslash = "/", mustWork = FALSE)
  .parent_candidates <- function(path) {
    out <- character()
    if (!.is_scalar_string(path)) return(out)
    cur <- .norm(path)
    if (file.exists(cur) && !dir.exists(cur)) cur <- dirname(cur)
    out <- c(out, cur)
    for (k in seq_len(max_up)) {
      parent <- dirname(cur)
      if (identical(parent, cur)) break
      out <- c(out, parent)
      cur <- parent
    }
    unique(out[nzchar(out)])
  }
  .is_pkg_root <- function(path) {
    if (!.is_scalar_string(path)) return(FALSE)
    all(vapply(required_rel_paths, function(rel) {
      file.exists(file.path(path, rel))
    }, logical(1)))
  }

  cands <- character()

  ca <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", ca, value = TRUE)
  if (length(file_arg)) {
    cands <- c(cands, .parent_candidates(sub("^--file=", "", file_arg[1L])))
  }

  frames <- sys.frames()
  ofiles <- unique(vapply(frames, function(fr) {
    val <- tryCatch(fr$ofile, error = function(e) "")
    if (.is_scalar_string(val)) val else ""
  }, character(1)))
  ofiles <- ofiles[nzchar(ofiles)]
  if (length(ofiles)) {
    for (p in ofiles) cands <- c(cands, .parent_candidates(p))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    active_path <- tryCatch(rstudioapi::getSourceEditorContext()$path,
                            error = function(e) "")
    if (.is_scalar_string(active_path)) {
      cands <- c(cands, .parent_candidates(active_path))
    }
  }

  wd <- .norm(getwd())
  cands <- c(cands, .parent_candidates(wd))
  cands <- c(cands, .parent_candidates(file.path(wd, script_name)))

  cands <- unique(cands[nzchar(cands)])
  hits <- cands[vapply(cands, .is_pkg_root, logical(1))]
  if (length(hits)) return(hits[1L])

  stop(
    paste0(
      "Could not locate the package root. Expected to find: ",
      paste(required_rel_paths, collapse = ", "),
      ".\n\n",
      "Safest invocation:\n",
      "  source('/full/path/to/", script_name, "', chdir = TRUE)\n",
      "or setwd('/full/path/to/package_root') before sourcing the script."
    ),
    call. = FALSE
  )
}

pkg_dir <- .find_pkg_dir(
  required_rel_paths = c("R/core_utils.R", "run_simulation.R"),
  script_name = "run_simulation.R"
)
setwd(pkg_dir)

# ---- Source package code -----------------------------------------------------
src_files <- c(
  "core_utils.R",
  "core_output.R",
  "core_targets.R",
  "nuisance_learners.R",
  "model_specs.R",
  "dgp.R",
  "truth.R",
  "ltmle_branch_registry.R",
  "estimator_msm_ipw.R",
  "estimator_gcomp.R",
  "estimator_ltmle_exact.R",
  "simulation_runner.R",
  "summarize.R",
  "figures.R"
)
for (f in src_files) {
  source(file.path(pkg_dir, "R", f), chdir = TRUE)
}


# ==============================================================================
# USER SETTINGS (edit this block)
# ==============================================================================

# ---- Reproducibility ---------------------------------------------------------
seed_base <- 20260531

# ---- Output directory --------------------------------------------------------
output_dir <- file.path("output", "simulation_8scenarios")

# ---- Simulation controls -----------------------------------------------------
R_reps <- 10L
B_mc   <- 2000L
B_boot <- 0L        # bootstrap replicates per dataset (0 = off)
conf_level <- 0.95

# ---- Core time axis ----------------------------------------------------------
# T defines the final outcome time Y_T.
# Package-specific timing convention:
#   - M1_0 and M2_0 are pre-treatment baseline mediators generated before A_0.
#   - At baseline visit (covariate time 0): M1_0, M2_0 -> A_0 -> L_0 -> Y_1.
#   - For follow-up visits t=1,..,T-1: A_t -> M1_t -> M2_t -> L_t -> Y_{t+1}.
T <- 5

# ---- Truth cohort + DGP sample sizes -----------------------------------------
B_truth <- 5000000L
truth_n_batches <- 25L
n_dgp   <- c(500L)

# ---- Treatment assignment in the DGP -----------------------------------------
# Supported: baseline_rct / sequential_rct / observational / fixed
DGP_treat_mech <- "baseline_rct"
p_rct <- 0.5

# ---- DGP coefficients --------------------------------------------------------
params0 <- default_dgp_params(T)
params0$treatment$p_rct <- p_rct

# Manual overrides example (keep commented):
# params0$Y$A <- 0.2
# params0$M1$A <- 0.3
# params0$M2$A <- 0.2

# ---- Scenario grid: pathway magnitude ----------------------------------------
# Manuscript mapping:
#   c1: A -> M1
#   c2: A -> M2
#   c3: all M1 -> Y terms incl. A:M1
#   c4: all M2 -> Y terms incl. A:M2 and delta
PM_MAP <- list(
  high = list(c1 = 0.50, c2 = 0.50, c3 = 0.50, c4 = 0.50),
  low  = list(c1 = 0.25, c2 = 0.25, c3 = 0.25, c4 = 0.25),
  `near-null` = list(c1 = 0.05, c2 = 0.05, c3 = 0.05, c4 = 0.05)
)

# ---- Scenario grid: residual mediator dependence (rho0, rho1) ---------------
# rho_mode controls how resid-corr-present residual correlation is assigned:
#   "rho1_only" : (0,rho_high)
#   "rho0_rho1" : (rho_low,rho_high)
rho_mode <- "rho0_rho1"
rho_low <- 0.1
rho_high <- 0.1


# ---- Scenario grid: mediator interdependence modes ---------------------------
# create_scenario_grid() includes:
#   none        : gamma0=0, gammaA=0, delta=0
#   gamma0-only : gamma0 kept at the base value below, gammaA=0, delta=0
#   full        : gamma0 and gammaA kept at the base values below, delta kept
#                 after PM$c4 scaling
#
# Base values used before PM scaling and scenario-specific zeroing.  Change
# these to tune structural M1 -> M2 dependence and the M1*M2 outcome
# interaction without editing R/dgp.R.
gamma0_base <- 0.1
gammaA_base <- 0.1
delta_base  <- 0.05
params0$gamma0 <- gamma0_base
params0$gammaA <- gammaA_base
params0$delta  <- delta_base


# ---- Regimens defining the estimands -----------------------------------------
# Active regimen a_bar and reference regimen a*_bar (length-1 recycled ok).
# NOTE: reg_a[1] corresponds to A_0.
reg_a  <- rep(1, T)  # e.g., always treated
reg_as <- rep(0, T)  # e.g., never treated
validate_baseline_rct_regimens(reg_a, reg_as, DGP_treat_mech)

# ---- Full ltmle_exact controls ----------------------------------------------
# The implementation solves the full exact-EIF score equations for all nine
# component means. Training-fold outcome bounds are used because Y is unbounded
# under the Gaussian DGP; fixed deterministic bounds may be supplied by setting
# ltmle_exact_y_bounds_mode <- "fixed" and providing ltmle_exact_y_bounds.

msm_ipw_trunc <- c(0.01, 0.99)
msm_ipw_truncation_enabled <- TRUE
msm_ipw_truncation_policy <- "quantile"
msm_ipw_truncation_quantile_lower <- 0.01
msm_ipw_truncation_quantile_upper <- 0.99
msm_ipw_truncation_target <- "final_cumulative_weight"
mediator_density_mc_n <- B_mc
ltmle_exact_learner <- "glm"
ltmle_exact_density_ratio_mc_n <- 2000L
ltmle_exact_verbose <- identical(tolower(Sys.getenv("LTMLE_EXACT_VERBOSE", unset = "true")), "true")
ltmle_exact_diagnostics_level <- "summary"
ltmle_exact_probability_bounds <- c(0.01, 0.99)
ltmle_exact_truncation_enabled <- TRUE
ltmle_exact_truncation_policy <- "quantile"
ltmle_exact_truncation_quantile_lower <- 0.01
ltmle_exact_truncation_quantile_upper <- 0.99
ltmle_exact_truncation_target <- "clever_covariate_H"
ltmle_exact_y_bounds_mode <- "train_fold"
ltmle_exact_y_bounds <- NULL
ltmle_exact_score_tolerance <- 1e-4
ltmle_exact_component_tolerance <- 1e-5
ltmle_exact_scaled_z_tolerance <- 3
finite_estimate_rule <- "finite numeric estimate; failed estimator attempts contribute NA rows to primary effects"
stop_on_estimator_failure <- identical(
  tolower(Sys.getenv("SIM_STOP_ON_ESTIMATOR_FAILURE", unset = "false")),
  "true"
)


# ---- Scenario subset for this run -------------------------------------------
scenario_pathway_levels <- c("high", "low")
scenario_rho_setting_levels <- c("resid-corr-present")
scenario_structure_levels <- c("none", "full")
scenario_Q_model_levels <- c("correct", "wrong")
include_near_null_scenarios <- FALSE
include_sensitivity_scenarios <- FALSE


# ---- Figure output -------------------------------------------
# If TRUE, write PDFs under output_dir.
make_figures <- TRUE
make_truth_figures <- TRUE


# ==============================================================================
# END USER PARAMETERS
# ==============================================================================


# ---- Build scenario grid -----------------------------------------------------
scenarios <- create_scenario_grid(
  n_vec = n_dgp,
  R = R_reps,
  B_truth_init = B_truth,
  B_MC_init = B_mc,
  truncation_rule = "T1",
  include_near_null = include_near_null_scenarios,
  include_sensitivity = include_sensitivity_scenarios,
  rho_mode = rho_mode,
  rho_low = rho_low,
  rho_high = rho_high
)

scenarios <- scenarios[
  scenarios$n %in% n_dgp &
    scenarios$pathway_setting %in% scenario_pathway_levels &
    scenarios$rho_setting %in% scenario_rho_setting_levels &
    scenarios$structure_setting %in% scenario_structure_levels &
    scenarios$Q_model %in% scenario_Q_model_levels,
  ,
  drop = FALSE
]

if (nrow(scenarios) == 0L) {
  stop("Scenario subset is empty. Check scenario_*_levels in USER SETTINGS.", call. = FALSE)
}

scenarios$scenario_label <- scenarios$scenario_id
scenarios$R_reps <- R_reps
scenarios$B_truth <- B_truth
scenarios$B_mc <- B_mc
scenarios$truth_n_batches <- truth_n_batches
scenarios$msm_ipw_trunc <- paste(msm_ipw_trunc, collapse = ",")
scenarios$msm_ipw_truncation_enabled <- msm_ipw_truncation_enabled
scenarios$msm_ipw_truncation_policy <- msm_ipw_truncation_policy
scenarios$msm_ipw_truncation_quantile_lower <- msm_ipw_truncation_quantile_lower
scenarios$msm_ipw_truncation_quantile_upper <- msm_ipw_truncation_quantile_upper
scenarios$msm_ipw_truncation_target <- msm_ipw_truncation_target
scenarios$mediator_density_mc_n <- mediator_density_mc_n
scenarios$ltmle_exact_density_ratio_mc_n <- ltmle_exact_density_ratio_mc_n
scenarios$ltmle_exact_verbose <- ltmle_exact_verbose
scenarios$ltmle_exact_diagnostics_level <- ltmle_exact_diagnostics_level
scenarios$ltmle_exact_probability_bounds <- paste(ltmle_exact_probability_bounds, collapse = ",")
scenarios$ltmle_exact_truncation_enabled <- ltmle_exact_truncation_enabled
scenarios$ltmle_exact_truncation_policy <- ltmle_exact_truncation_policy
scenarios$ltmle_exact_truncation_quantile_lower <- ltmle_exact_truncation_quantile_lower
scenarios$ltmle_exact_truncation_quantile_upper <- ltmle_exact_truncation_quantile_upper
scenarios$ltmle_exact_truncation_target <- ltmle_exact_truncation_target
scenarios$ltmle_exact_score_tolerance <- ltmle_exact_score_tolerance
scenarios$ltmle_exact_component_tolerance <- ltmle_exact_component_tolerance
scenarios$finite_estimate_rule <- finite_estimate_rule

# ---- Estimator registry ------------------------------------------------------
# Manuscript outputs distinguish the predeclared truncation variants.
estimators <- list()

estimators$msm_ipw_quantile_truncated <- function(dat, cfg) {
  msm_ipw_estimate_worldmeans(
    dat = dat,
    T = cfg$T,
    reg_a = cfg$reg_a,
    reg_as = cfg$reg_as,
    treat_mech = cfg$DGP_treat_mech,
    p_rct = cfg$p_rct,
    trunc = cfg$msm_ipw_trunc,
    truncation_enabled = TRUE,
    truncation_policy = "quantile",
    truncation_quantile_lower = cfg$msm_ipw_truncation_quantile_lower,
    truncation_quantile_upper = cfg$msm_ipw_truncation_quantile_upper,
    truncation_target = cfg$msm_ipw_truncation_target,
    Q_model = cfg$Q_model,
    B_mc = cfg$B_mc,
    mediator_density_mc_n = cfg$mediator_density_mc_n,
    seed = cfg$seed
  )
}

estimators$msm_ipw_untruncated <- function(dat, cfg) {
  msm_ipw_estimate_worldmeans(
    dat = dat,
    T = cfg$T,
    reg_a = cfg$reg_a,
    reg_as = cfg$reg_as,
    treat_mech = cfg$DGP_treat_mech,
    p_rct = cfg$p_rct,
    truncation_enabled = FALSE,
    truncation_policy = "none",
    truncation_target = cfg$msm_ipw_truncation_target,
    Q_model = cfg$Q_model,
    B_mc = cfg$B_mc,
    mediator_density_mc_n = cfg$mediator_density_mc_n,
    seed = cfg$seed
  )
}

estimators$gcomp <- function(dat, cfg) {
  gcomp_estimate_worldmeans(
    dat = dat,
    T = cfg$T,
    reg_a = cfg$reg_a,
    reg_as = cfg$reg_as,
    B_mc = cfg$B_mc,
    Q_model = cfg$Q_model,
    seed = cfg$seed,
    progress = FALSE,
    treat_mech = cfg$DGP_treat_mech
  )
}

estimators$ltmle_exact_quantile_truncated <- function(dat, cfg) {
  ltmle_exact_estimate_worldmeans(
    dat = dat, T = cfg$T,
    reg_a = cfg$reg_a, reg_as = cfg$reg_as,
    treat_mech = cfg$DGP_treat_mech,
    p_rct = cfg$p_rct,
    learner = cfg$ltmle_exact_learner,
    Q_model = cfg$Q_model,
    probability_bounds = cfg$ltmle_exact_probability_bounds,
    truncation_enabled = cfg$ltmle_exact_truncation_enabled,
    truncation_policy = cfg$ltmle_exact_truncation_policy,
    truncation_quantile_lower = cfg$ltmle_exact_truncation_quantile_lower,
    truncation_quantile_upper = cfg$ltmle_exact_truncation_quantile_upper,
    truncation_target = cfg$ltmle_exact_truncation_target,
    ltmle_exact_density_ratio_mc_n = cfg$ltmle_exact_density_ratio_mc_n,
    diagnostics_level = cfg$ltmle_exact_diagnostics_level,
    y_bounds_mode = cfg$ltmle_exact_y_bounds_mode,
    y_bounds = cfg$ltmle_exact_y_bounds,
    score_tolerance = cfg$ltmle_exact_score_tolerance,
    component_tolerance = cfg$ltmle_exact_component_tolerance,
    scaled_z_tolerance = cfg$ltmle_exact_scaled_z_tolerance,
    seed = cfg$seed,
    verbose = cfg$ltmle_exact_verbose
  )
}

estimators$ltmle_exact_untruncated <- function(dat, cfg) {
  ltmle_exact_estimate_worldmeans(
    dat = dat, T = cfg$T,
    reg_a = cfg$reg_a, reg_as = cfg$reg_as,
    treat_mech = cfg$DGP_treat_mech,
    p_rct = cfg$p_rct,
    learner = cfg$ltmle_exact_learner,
    Q_model = cfg$Q_model,
    probability_bounds = cfg$ltmle_exact_probability_bounds,
    truncation_enabled = FALSE,
    truncation_policy = "none",
    truncation_target = cfg$ltmle_exact_truncation_target,
    ltmle_exact_density_ratio_mc_n = cfg$ltmle_exact_density_ratio_mc_n,
    diagnostics_level = cfg$ltmle_exact_diagnostics_level,
    y_bounds_mode = cfg$ltmle_exact_y_bounds_mode,
    y_bounds = cfg$ltmle_exact_y_bounds,
    score_tolerance = cfg$ltmle_exact_score_tolerance,
    component_tolerance = cfg$ltmle_exact_component_tolerance,
    scaled_z_tolerance = cfg$ltmle_exact_scaled_z_tolerance,
    seed = cfg$seed,
    verbose = cfg$ltmle_exact_verbose
  )
}

# ---- Master configuration passed to runner -----------------------------------
run_cfg <- list(
  output_dir = output_dir,
  T = T,
  reg_a = reg_a,
  reg_as = reg_as,
  R_reps = R_reps,
  B_truth = B_truth,
  truth_n_batches = truth_n_batches,
  B_mc = B_mc,
  mediator_density_mc_n = mediator_density_mc_n,
  seed = seed_base,
  DGP_treat_mech = DGP_treat_mech,
  p_rct = p_rct,
  PM_MAP = PM_MAP,
  rho_mode = rho_mode,
  rho_low = rho_low,
  rho_high = rho_high,
  gamma0_base = gamma0_base,
  gammaA_base = gammaA_base,
  delta_base = delta_base,
  B_boot = B_boot,
  conf_level = conf_level,
  msm_ipw_trunc = msm_ipw_trunc,
  msm_ipw_truncation_enabled = msm_ipw_truncation_enabled,
  msm_ipw_truncation_policy = msm_ipw_truncation_policy,
  msm_ipw_truncation_quantile_lower = msm_ipw_truncation_quantile_lower,
  msm_ipw_truncation_quantile_upper = msm_ipw_truncation_quantile_upper,
  msm_ipw_truncation_target = msm_ipw_truncation_target,
  ltmle_exact_learner = ltmle_exact_learner,
  ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
  ltmle_exact_verbose = ltmle_exact_verbose,
  ltmle_exact_diagnostics_level = ltmle_exact_diagnostics_level,
  ltmle_exact_probability_bounds = ltmle_exact_probability_bounds,
  ltmle_exact_truncation_enabled = ltmle_exact_truncation_enabled,
  ltmle_exact_truncation_policy = ltmle_exact_truncation_policy,
  ltmle_exact_truncation_quantile_lower = ltmle_exact_truncation_quantile_lower,
  ltmle_exact_truncation_quantile_upper = ltmle_exact_truncation_quantile_upper,
  ltmle_exact_truncation_target = ltmle_exact_truncation_target,
  ltmle_exact_y_bounds_mode = ltmle_exact_y_bounds_mode,
  ltmle_exact_y_bounds = ltmle_exact_y_bounds,
  ltmle_exact_score_tolerance = ltmle_exact_score_tolerance,
  ltmle_exact_component_tolerance = ltmle_exact_component_tolerance,
  ltmle_exact_scaled_z_tolerance = ltmle_exact_scaled_z_tolerance,
  finite_estimate_rule = finite_estimate_rule,
  stop_on_estimator_failure = stop_on_estimator_failure
)


# ---- Diagnostics settings (optional but recommended) -------------------------
# These settings collect lightweight per-rep diagnostics and save full RDS dumps
# only for outlier replications. All diagnostics go under output_dir/diagnostics.

run_cfg$diagnostics <- list(
  enabled = TRUE,
  dir = file.path(output_dir, "diagnostics"),
  save_full = "none",        # none/outlier/all
  save_score_equations = "outlier",    # none/outlier/all
  thr = list(
    max_abs_effect = 1e6,    # outlier if any effect estimate exceeds this
    max_H = 500,             # outlier if max(H) exceeds this
    min_ess_frac = 0.05,     # outlier if ESS/N is below this
    max_abs_eps = 25         # outlier if |eps| exceeds this
  ),
  outlier_triggers = c(
    "nonfinite estimate",
    "estimator failure",
    "absolute effect estimate exceeds diagnostics$thr$max_abs_effect",
    "ltmle_exact clever covariate H exceeds diagnostics$thr$max_H",
    "ltmle_exact ESS fraction below diagnostics$thr$min_ess_frac",
    "absolute targeting epsilon exceeds diagnostics$thr$max_abs_eps"
  )
)

# ---- Consistency check: p_rct used in DGP vs estimators ----------------------
# DGP uses params0$treatment$p_rct; IPW/ltmle_exact use run_cfg$p_rct.
# If these diverge, estimates can be catastrophically wrong.
if (!is.null(params0$treatment$p_rct) && !is.null(run_cfg$p_rct)) {
  if (!isTRUE(all.equal(params0$treatment$p_rct, run_cfg$p_rct))) {
    stop("p_rct mismatch: params0$treatment$p_rct != run_cfg$p_rct")
  }
}

# ---- Run simulation ----------------------------------------------------------
# Persist run configuration and session information for auditability.
# These files allow reviewers to verify that the reported settings (rep count,
# MC sizes, truncation, and ltmle_exact controls) match the produced outputs.
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(run_cfg, file.path(output_dir, "run_config.rds"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

res <- run_simulation(
  scenarios = scenarios,
  params0 = params0,
  run_cfg = run_cfg,
  estimators = estimators,
  R_reps = R_reps,
  progress_every_dgp = 1L,
  progress_every_analysis = 1L
)

# ---- Guardrail: enforce presence of requested estimators ---------------------
# Prevents "silent disappearance" of an estimator when it fails in all reps.
requested <- names(estimators)
present <- if (!is.null(res$estimates_long) && nrow(res$estimates_long) > 0) {
  unique(res$estimates_long$estimator)
} else {
  character(0)
}

missing <- setdiff(requested, present)
if (length(missing) > 0) {
  stop(
    "Missing estimator outputs (no rows in estimates_long): ", paste(missing, collapse = ", "),
    ". Check ", file.path(output_dir, "diagnostics", "estimator_failures.csv"),
    " and package availability.",
    call. = FALSE
  )
}


# ---- Write outputs -----------------------------------------------------------
truth_effects_out <- order_effects_long(res$truth_effects)
estimates_long_out <- order_effects_long(.add_effect_truth_bias(res$estimates_long, res$truth_effects))
performance_out <- order_effects_long(res$performance)

# Replicate-level world-mean estimates (mu_*) for diagnosing which world broke
worldmeans_long_out <- .add_worldmean_truth_bias(res$worldmeans_estimates_long, res$truth_worldmeans)
if (!is.null(worldmeans_long_out) && nrow(worldmeans_long_out) > 0) {
  worldmeans_long_out <- worldmeans_long_out[order(worldmeans_long_out$scenario_id,
                                                  worldmeans_long_out$n,
                                                  worldmeans_long_out$rep,
                                                  worldmeans_long_out$estimator,
                                                  worldmeans_long_out$world), , drop = FALSE]
}

write.csv(scenarios,            file = file.path(output_dir, "scenario_manifest.csv"), row.names = FALSE)
write.csv(res$truth_worldmeans, file = file.path(output_dir, "truth_worldmeans.csv"), row.names = FALSE)
write.csv(truth_effects_out,    file = file.path(output_dir, "truth_effects.csv"),    row.names = FALSE)

if (isTRUE(make_truth_figures)) {
  truth_figures_written <- FALSE
  tryCatch({
    create_truth_figures(truth_df = truth_effects_out,
                         out_pdf = file.path(output_dir, "truth_figures.pdf"))
    message("Wrote ", file.path(output_dir, "truth_figures.pdf"))
    truth_figures_written <- TRUE
  }, error = function(e) {
    warning("figure_generation_failed: truth_figures: ", conditionMessage(e), call. = FALSE)
  })
} else {
  truth_figures_written <- FALSE
}
write.csv(estimates_long_out,   file = file.path(output_dir, "estimates_long.csv"),   row.names = FALSE)
write.csv(worldmeans_long_out,  file = file.path(output_dir, "worldmeans_estimates_long.csv"), row.names = FALSE)
write.csv(performance_out,      file = file.path(output_dir, "performance_summary.csv"), row.names = FALSE)

cat("\nDone. Output files written to:\n")
cat("  ", normalizePath(file.path(output_dir, "scenario_manifest.csv")), "\n")
cat("  ", normalizePath(file.path(output_dir, "truth_worldmeans.csv")), "\n")
cat("  ", normalizePath(file.path(output_dir, "truth_effects.csv")), "\n")
if (isTRUE(truth_figures_written)) {
  cat("  ", normalizePath(file.path(output_dir, "truth_figures.pdf")), "\n")
}
cat("  ", normalizePath(file.path(output_dir, "estimates_long.csv")), "\n")
cat("  ", normalizePath(file.path(output_dir, "worldmeans_estimates_long.csv")), "\n")
cat("  ", normalizePath(file.path(output_dir, "performance_summary.csv")), "\n")

# ---- Write diagnostics (CSV) -------------------------------------------------
if (!is.null(res$diagnostics)) {
  diag_dir <- file.path(output_dir, "diagnostics")
  dir.create(diag_dir, showWarnings = FALSE, recursive = TRUE)

  if (!is.null(worldmeans_long_out) && nrow(worldmeans_long_out) > 0 &&
      !is.null(res$truth_worldmeans) && nrow(res$truth_worldmeans) > 0) {
    wm_bias <- summarize_worldmean_bias(worldmeans_long_out, res$truth_worldmeans)
    write.csv(wm_bias$detail,
              file = file.path(diag_dir, "worldmean_bias_by_estimator.csv"),
              row.names = FALSE)
    write.csv(wm_bias$summary,
              file = file.path(diag_dir, "worldmean_bias_summary.csv"),
              row.names = FALSE)
  }

  if (!is.null(res$diagnostics$ltmle_exact_score_equations)) {
    write.csv(res$diagnostics$ltmle_exact_score_equations,
              file = file.path(diag_dir, "ltmle_exact_score_equations.csv"),
              row.names = FALSE)
  }
  if (!is.null(res$diagnostics$ltmle_exact_run)) {
    write.csv(res$diagnostics$ltmle_exact_run,
              file = file.path(diag_dir, "ltmle_exact_run_summary.csv"),
              row.names = FALSE)
  }
  if (!is.null(res$diagnostics$ltmle_exact_fold)) {
    write.csv(res$diagnostics$ltmle_exact_fold,
              file = file.path(diag_dir, "ltmle_exact_fold_bounds.csv"),
              row.names = FALSE)
  }
  if (!is.null(res$diagnostics$ltmle_exact_component_law_registry)) {
    write.csv(res$diagnostics$ltmle_exact_component_law_registry,
              file = file.path(diag_dir, "ltmle_exact_component_law_registry.csv"),
              row.names = FALSE)
  }
  if (!is.null(res$diagnostics$ltmle_exact_factor_tasks)) {
    write.csv(res$diagnostics$ltmle_exact_factor_tasks,
              file = file.path(diag_dir, "ltmle_exact_factor_tasks.csv"),
              row.names = FALSE)
  }
  if (!is.null(res$diagnostics$ltmle_exact_component_eif_summary)) {
    write.csv(res$diagnostics$ltmle_exact_component_eif_summary,
              file = file.path(diag_dir, "ltmle_exact_component_eif_summary.csv"),
              row.names = FALSE)
  }
  if (!is.null(res$diagnostics$ltmle_exact_component_eif_terms)) {
    write.csv(res$diagnostics$ltmle_exact_component_eif_terms,
              file = file.path(diag_dir, "ltmle_exact_component_eif_terms.csv"),
              row.names = FALSE)
  }
  if (!is.null(res$diagnostics$estimator_attempt_status)) {
    write.csv(res$diagnostics$estimator_attempt_status,
              file = file.path(diag_dir, "estimator_attempt_status.csv"),
              row.names = FALSE)
  }
  if (!is.null(res$diagnostics$estimator_runtime_summary)) {
    write.csv(res$diagnostics$estimator_runtime_summary,
              file = file.path(diag_dir, "estimator_runtime_summary.csv"),
              row.names = FALSE)
  }
  if (!is.null(res$diagnostics$truncation_diagnostics)) {
    write.csv(res$diagnostics$truncation_diagnostics,
              file = file.path(diag_dir, "truncation_diagnostics.csv"),
              row.names = FALSE)
  }
  if (!is.null(res$diagnostics$truncation_diagnostics_summary)) {
    write.csv(res$diagnostics$truncation_diagnostics_summary,
              file = file.path(diag_dir, "truncation_diagnostics_summary.csv"),
              row.names = FALSE)
  }
  if (!is.null(res$diagnostics$msm)) {
    write.csv(res$diagnostics$msm,
              file = file.path(diag_dir, "msm_weight_diagnostics.csv"),
              row.names = FALSE)
  }
  if (!is.null(res$diagnostics$gcomp)) {
    write.csv(res$diagnostics$gcomp,
              file = file.path(diag_dir, "gcomp_mcse.csv"),
              row.names = FALSE)
  }
  if (!is.null(res$diagnostics$failures)) {
    write.csv(res$diagnostics$failures,
              file = file.path(diag_dir, "estimator_failures.csv"),
              row.names = FALSE)
  }
  if (!is.null(res$diagnostics$full_files)) {
    write.csv(data.frame(path = res$diagnostics$full_files, stringsAsFactors = FALSE),
              file = file.path(diag_dir, "full_rds_index.csv"),
              row.names = FALSE)
  }
}


# ---- Create Figures ----------------------------------------------------------
if (make_figures) {
  tryCatch({
    create_figures(performance_out, out_pdf = file.path(output_dir, "performance_figures.pdf"))
  }, error = function(e) {
    warning("figure_generation_failed: performance_figures: ", conditionMessage(e), call. = FALSE)
  })
}
