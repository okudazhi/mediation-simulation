################################################################################
# run_simulation_parallel.R
#
# PARALLEL ENTRY POINT:
# 1) Edit USER SETTINGS in run_simulation.R.
# 2) Edit PARALLEL SETTINGS below if needed.
# 3) Safest run: source("/full/path/to/run_simulation_parallel.R", chdir = TRUE)
#
# Final integrated outputs are written to the output_dir set in run_simulation.R.
#
# Replication-shard intermediate RDS files are written under:
# <package_root>/<output_dir>/parallel/shards/
################################################################################

# ---- PARALLEL SETTINGS -------------------------------------------------------
parallel_splits <- 5L
parallel_workers <- NULL

# ---- Load serial wrapper setup without running it ----------------------------
.parallel_read_serial_wrapper <- function(path = "run_simulation.R") {
  if (!file.exists(path)) {
    stop("Could not find ", path, ". Run from the package root.", call. = FALSE)
  }
  readLines(path, warn = FALSE)
}

.parallel_eval_until <- function(lines, end_pattern, env = parent.frame()) {
  end_idx <- grep(end_pattern, lines, fixed = TRUE)
  if (length(end_idx) == 0L) {
    stop("Could not find marker in run_simulation.R: ", end_pattern, call. = FALSE)
  }
  eval(parse(text = lines[seq_len(end_idx[1L])]), envir = env)
}

.parallel_eval_between <- function(lines, start_pattern, end_pattern, env = parent.frame()) {
  start_idx <- grep(start_pattern, lines, fixed = TRUE)
  end_idx <- grep(end_pattern, lines, fixed = TRUE)
  if (length(start_idx) == 0L) {
    stop("Could not find marker in run_simulation.R: ", start_pattern, call. = FALSE)
  }
  if (length(end_idx) == 0L) {
    stop("Could not find marker in run_simulation.R: ", end_pattern, call. = FALSE)
  }
  start_idx <- start_idx[1L]
  end_idx <- end_idx[end_idx > start_idx][1L] - 1L
  if (is.na(end_idx) || end_idx < start_idx) {
    stop("Invalid marker order in run_simulation.R.", call. = FALSE)
  }
  eval(parse(text = lines[start_idx:end_idx]), envir = env)
}

serial_lines <- .parallel_read_serial_wrapper("run_simulation.R")
.parallel_eval_until(serial_lines, "# END USER PARAMETERS", env = .GlobalEnv)
.parallel_eval_between(
  serial_lines,
  "# ---- Build scenario grid",
  "# ---- Run simulation",
  env = .GlobalEnv
)


# ---- Helpers -----------------------------------------------------------------
.fmt_hms_parallel <- function(seconds) {
  if (is.null(seconds) || !is.finite(seconds) || seconds < 0) return(NA_character_)
  s <- as.integer(round(seconds))
  h <- s %/% 3600L
  m <- (s %% 3600L) %/% 60L
  ss <- s %% 60L
  sprintf("%02d:%02d:%02d", h, m, ss)
}

.parallel_split_indices <- function(n, k) {
  if (!is.numeric(n) || length(n) != 1L || n <= 0L) {
    stop("n must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(k) || length(k) != 1L || k <= 0L) {
    stop("parallel_splits must be a positive integer.", call. = FALSE)
  }
  k <- min(as.integer(k), as.integer(n))
  base <- as.integer(n) %/% k
  extra <- as.integer(n) %% k
  sizes <- c(rep(base + 1L, extra), rep(base, k - extra))
  sizes <- sizes[sizes > 0L]
  ends <- cumsum(sizes)
  starts <- ends - sizes + 1L
  Map(seq.int, starts, ends)
}

.parallel_rbind <- function(xs) {
  xs <- Filter(function(x) !is.null(x) && is.data.frame(x) && nrow(x) > 0L, xs)
  if (length(xs) == 0L) return(NULL)
  do.call(rbind, xs)
}

.parallel_combine_diagnostics <- function(res_list) {
  diag_list <- lapply(res_list, function(x) x$diagnostics)
  diag_list <- Filter(Negate(is.null), diag_list)
  if (length(diag_list) == 0L) return(NULL)

  .get <- function(name) .parallel_rbind(lapply(diag_list, function(d) d[[name]]))
  full_files <- unlist(lapply(diag_list, function(d) d$full_files), use.names = FALSE)
  if (length(full_files) == 0L) full_files <- NULL

  list(
    ltmle_exact_score_equations = .get("ltmle_exact_score_equations"),
    ltmle_exact_run = .get("ltmle_exact_run"),
    ltmle_exact_fold = .get("ltmle_exact_fold"),
    ltmle_exact_component_law_registry = .get("ltmle_exact_component_law_registry"),
    ltmle_exact_factor_tasks = .get("ltmle_exact_factor_tasks"),
    ltmle_exact_component_eif_summary = .get("ltmle_exact_component_eif_summary"),
    ltmle_exact_component_eif_terms = .get("ltmle_exact_component_eif_terms"),
    estimator_attempt_status = .get("estimator_attempt_status"),
    estimator_runtime_summary = .get("estimator_runtime_summary"),
    truncation_diagnostics = .get("truncation_diagnostics"),
    truncation_diagnostics_summary = .get("truncation_diagnostics_summary"),
    msm = .get("msm"),
    gcomp = .get("gcomp"),
    failures = .get("failures"),
    full_files = full_files
  )
}

.parallel_add_performance_descriptors <- function(perf, scenarios) {
  scenario_desc_cols <- intersect(
    c("scenario_id", "scenario_label", "analysis_tier", "pathway_setting", "rho_setting", "structure_setting",
      "Q_model", "fold_count", "R", "B_truth_init", "B_MC_init", "truncation_rule",
      "R_reps", "B_truth", "B_mc", "truth_n_batches", "msm_ipw_trunc",
      "msm_ipw_truncation_enabled", "msm_ipw_truncation_policy",
      "msm_ipw_truncation_quantile_lower", "msm_ipw_truncation_quantile_upper",
      "msm_ipw_truncation_target", "mediator_density_mc_n", "ltmle_exact_density_ratio_mc_n",
      "ltmle_exact_probability_bounds", "ltmle_exact_truncation_enabled",
      "ltmle_exact_truncation_policy", "ltmle_exact_truncation_quantile_lower",
      "ltmle_exact_truncation_quantile_upper", "ltmle_exact_truncation_target",
      "ltmle_exact_score_tolerance",
      "ltmle_exact_component_tolerance", "finite_estimate_rule",
      "PM", "rho", "rho0", "rho1", "MI_mode"),
    names(scenarios)
  )
  perf <- merge(
    perf,
    unique(scenarios[, scenario_desc_cols, drop = FALSE]),
    by = "scenario_id",
    all.x = TRUE,
    sort = FALSE
  )

  core_cols <- c("scenario_id", "scenario_label", "analysis_tier", "pathway_setting", "rho_setting", "structure_setting",
                 "Q_model", "fold_count", "PM", "rho", "rho0", "rho1", "MI_mode", "n",
                 "estimator", "estimand", "effect")
  other_cols <- setdiff(names(perf), core_cols)
  perf[, c(core_cols, other_cols), drop = FALSE]
}

.parallel_write_diagnostics <- function(diagnostics, diag_dir = file.path(output_dir, "diagnostics")) {
  if (is.null(diagnostics)) return(invisible(NULL))

  dir.create(diag_dir, showWarnings = FALSE, recursive = TRUE)

  if (!is.null(diagnostics$ltmle_exact_score_equations)) {
    write.csv(diagnostics$ltmle_exact_score_equations,
              file = file.path(diag_dir, "ltmle_exact_score_equations.csv"),
              row.names = FALSE)
  }
  if (!is.null(diagnostics$ltmle_exact_run)) {
    write.csv(diagnostics$ltmle_exact_run,
              file = file.path(diag_dir, "ltmle_exact_run_summary.csv"),
              row.names = FALSE)
  }
  if (!is.null(diagnostics$ltmle_exact_fold)) {
    write.csv(diagnostics$ltmle_exact_fold,
              file = file.path(diag_dir, "ltmle_exact_fold_bounds.csv"),
              row.names = FALSE)
  }
  if (!is.null(diagnostics$ltmle_exact_component_law_registry)) {
    write.csv(diagnostics$ltmle_exact_component_law_registry,
              file = file.path(diag_dir, "ltmle_exact_component_law_registry.csv"),
              row.names = FALSE)
  }
  if (!is.null(diagnostics$ltmle_exact_factor_tasks)) {
    write.csv(diagnostics$ltmle_exact_factor_tasks,
              file = file.path(diag_dir, "ltmle_exact_factor_tasks.csv"),
              row.names = FALSE)
  }
  if (!is.null(diagnostics$ltmle_exact_component_eif_summary)) {
    write.csv(diagnostics$ltmle_exact_component_eif_summary,
              file = file.path(diag_dir, "ltmle_exact_component_eif_summary.csv"),
              row.names = FALSE)
  }
  if (!is.null(diagnostics$ltmle_exact_component_eif_terms)) {
    write.csv(diagnostics$ltmle_exact_component_eif_terms,
              file = file.path(diag_dir, "ltmle_exact_component_eif_terms.csv"),
              row.names = FALSE)
  }
  if (!is.null(diagnostics$estimator_attempt_status)) {
    write.csv(diagnostics$estimator_attempt_status,
              file = file.path(diag_dir, "estimator_attempt_status.csv"),
              row.names = FALSE)
  }
  if (!is.null(diagnostics$estimator_runtime_summary)) {
    write.csv(diagnostics$estimator_runtime_summary,
              file = file.path(diag_dir, "estimator_runtime_summary.csv"),
              row.names = FALSE)
  }
  if (!is.null(diagnostics$truncation_diagnostics)) {
    write.csv(diagnostics$truncation_diagnostics,
              file = file.path(diag_dir, "truncation_diagnostics.csv"),
              row.names = FALSE)
  }
  if (!is.null(diagnostics$truncation_diagnostics_summary)) {
    write.csv(diagnostics$truncation_diagnostics_summary,
              file = file.path(diag_dir, "truncation_diagnostics_summary.csv"),
              row.names = FALSE)
  }
  if (!is.null(diagnostics$msm)) {
    write.csv(diagnostics$msm,
              file = file.path(diag_dir, "msm_weight_diagnostics.csv"),
              row.names = FALSE)
  }
  if (!is.null(diagnostics$gcomp)) {
    write.csv(diagnostics$gcomp,
              file = file.path(diag_dir, "gcomp_mcse.csv"),
              row.names = FALSE)
  }
  if (!is.null(diagnostics$failures)) {
    write.csv(diagnostics$failures,
              file = file.path(diag_dir, "estimator_failures.csv"),
              row.names = FALSE)
  }
  if (!is.null(diagnostics$full_files)) {
    write.csv(data.frame(path = diagnostics$full_files, stringsAsFactors = FALSE),
              file = file.path(diag_dir, "full_rds_index.csv"),
              row.names = FALSE)
  }

  invisible(NULL)
}

.parallel_run_shard <- function(job) {
  shard_id <- job$shard_id
  rep_idx <- job$rep_indices
  shard_dir <- file.path(output_dir, "parallel", "shards")
  diag_dir <- file.path(output_dir, "diagnostics", "parallel_shards",
                        sprintf("shard_%02d", shard_id))
  started_at <- Sys.time()
  worker_pid <- Sys.getpid()

  run_cfg_shard <- run_cfg
  run_cfg_shard$rep_indices <- rep_idx
  run_cfg_shard$skip_truth <- TRUE
  run_cfg_shard$skip_performance <- TRUE
  if (!is.null(run_cfg_shard$diagnostics)) {
    run_cfg_shard$diagnostics$dir <- diag_dir
  }

  cat(sprintf(
    "[Parallel] Shard %02d started | reps=%s | scenario rows=%d | seed=%d\n",
    shard_id, paste(rep_idx, collapse = ","), nrow(scenarios), as.integer(run_cfg_shard$seed)
  ))

  t0 <- proc.time()[["elapsed"]]
  shard_error <- NULL
  res <- tryCatch(
    run_simulation(
      scenarios = scenarios,
      params0 = params0,
      run_cfg = run_cfg_shard,
      estimators = estimators,
      R_reps = R_reps,
      show_progress = TRUE,
      progress_every_dgp = 1L,
      progress_every_analysis = 1L
    ),
    error = function(e) {
      shard_error <<- e
      NULL
    }
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  completed_at <- Sys.time()

  out_path <- file.path(shard_dir, sprintf("shard_%02d_result.rds", shard_id))
  if (!is.null(res)) {
    saveRDS(
      list(
        shard_id = shard_id,
        rep_indices = rep_idx,
        seed = run_cfg_shard$seed,
        elapsed_sec = elapsed,
        result = res
      ),
      out_path
    )
  } else {
    out_path <- NA_character_
  }

  status <- if (is.null(shard_error)) "completed" else "failed"
  cat(sprintf("[Parallel] Shard %02d %s | elapsed=%s\n",
              shard_id, status, .fmt_hms_parallel(elapsed)))

  list(
    shard_id = shard_id,
    rep_indices = rep_idx,
    seed = run_cfg_shard$seed,
    status = status,
    started_at = format(started_at, "%Y-%m-%d %H:%M:%S %Z"),
    completed_at = format(completed_at, "%Y-%m-%d %H:%M:%S %Z"),
    elapsed_sec = elapsed,
    worker_pid = worker_pid,
    result_path = out_path,
    error_message = if (!is.null(shard_error)) conditionMessage(shard_error) else NA_character_,
    n_failed_estimator_attempts = if (!is.null(res$diagnostics$estimator_attempt_status)) {
      sum(as.logical(res$diagnostics$estimator_attempt_status$failed), na.rm = TRUE)
    } else NA_integer_,
    result = res
  )
}


# ---- Validate and prepare shards --------------------------------------------
parallel_splits <- as.integer(parallel_splits)
if (length(parallel_splits) != 1L || is.na(parallel_splits) || parallel_splits <= 0L) {
  stop("parallel_splits must be a positive integer.", call. = FALSE)
}

if (is.null(parallel_workers)) {
  parallel_workers <- parallel_splits
}
parallel_workers <- as.integer(parallel_workers)
if (length(parallel_workers) != 1L || is.na(parallel_workers) || parallel_workers <= 0L) {
  stop("parallel_workers must be NULL or a positive integer.", call. = FALSE)
}

shard_rep_indices <- .parallel_split_indices(as.integer(R_reps), parallel_splits)
parallel_workers <- min(parallel_workers, length(shard_rep_indices))

if (.Platform$OS.type == "windows" && parallel_workers > 1L) {
  warning("parallel::mclapply is not parallel on Windows; using one worker.", call. = FALSE)
  parallel_workers <- 1L
}

dir.create(file.path(output_dir, "parallel", "shards"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "diagnostics"), showWarnings = FALSE, recursive = TRUE)

shard_manifest <- do.call(rbind, lapply(seq_along(shard_rep_indices), function(i) {
  rep_idx <- shard_rep_indices[[i]]
  data.frame(
    shard_id = i,
    rep_start = min(rep_idx),
    rep_end = max(rep_idx),
    n_reps = length(rep_idx),
    rep_count = length(rep_idx),
    rep_indices = paste(rep_idx, collapse = ","),
    scenario_rows = nrow(scenarios),
    seed = as.integer(run_cfg$seed),
    status = "pending",
    started_at = NA_character_,
    completed_at = NA_character_,
    elapsed_seconds = NA_real_,
    worker_pid = NA_integer_,
    output_file = NA_character_,
    error_message = NA_character_,
    n_failed_estimator_attempts = NA_integer_,
    stringsAsFactors = FALSE
  )
}))
write.csv(shard_manifest, file = file.path(output_dir, "parallel", "shard_manifest.csv"),
          row.names = FALSE)

run_cfg_parallel <- run_cfg
run_cfg_parallel$parallel <- list(
  enabled = TRUE,
  splits_requested = parallel_splits,
  split_unit = "replication",
  shards = length(shard_rep_indices),
  workers = parallel_workers,
  shard_manifest = shard_manifest
)
saveRDS(run_cfg_parallel, file.path(output_dir, "run_config.rds"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

cat(sprintf(
  "[Parallel] Prepared %d replications across %d scenario rows into %d replication shards | workers=%d\n",
  as.integer(R_reps), nrow(scenarios), length(shard_rep_indices), parallel_workers
))
cat("[Parallel] Shard manifest: ", normalizePath(file.path(output_dir, "parallel", "shard_manifest.csv")), "\n", sep = "")

# ---- Run shards --------------------------------------------------------------
t_parallel0 <- proc.time()[["elapsed"]]
jobs <- lapply(seq_along(shard_rep_indices), function(i) {
  list(shard_id = i, rep_indices = shard_rep_indices[[i]])
})

if (parallel_workers > 1L) {
  shard_results <- parallel::mclapply(
    jobs,
    .parallel_run_shard,
    mc.cores = parallel_workers,
    mc.preschedule = FALSE
  )
} else {
  shard_results <- lapply(jobs, .parallel_run_shard)
}

parallel_elapsed <- proc.time()[["elapsed"]] - t_parallel0

shard_manifest_final <- do.call(rbind, lapply(shard_results, function(x) {
  data.frame(
    shard_id = x$shard_id,
    rep_start = min(x$rep_indices),
    rep_end = max(x$rep_indices),
    n_reps = length(x$rep_indices),
    rep_count = length(x$rep_indices),
    rep_indices = paste(x$rep_indices, collapse = ","),
    scenario_rows = nrow(scenarios),
    seed = as.integer(x$seed),
    status = x$status,
    started_at = x$started_at,
    completed_at = x$completed_at,
    elapsed_seconds = as.numeric(x$elapsed_sec),
    worker_pid = as.integer(x$worker_pid),
    output_file = x$result_path,
    error_message = x$error_message,
    n_failed_estimator_attempts = as.integer(x$n_failed_estimator_attempts),
    stringsAsFactors = FALSE
  )
}))
write.csv(shard_manifest_final, file = file.path(output_dir, "parallel", "shard_manifest.csv"),
          row.names = FALSE)

failed_shards <- shard_manifest_final[shard_manifest_final$status != "completed", , drop = FALSE]
if (nrow(failed_shards)) {
  stop(
    "One or more parallel shards failed: ",
    paste(sprintf("shard_%02d: %s", failed_shards$shard_id, failed_shards$error_message),
          collapse = " | "),
    call. = FALSE
  )
}

res_list <- lapply(shard_results, function(x) x$result)
estimates_long <- .parallel_rbind(lapply(res_list, function(x) x$estimates_long))
worldmeans_estimates_long <- .parallel_rbind(lapply(res_list, function(x) x$worldmeans_estimates_long))
diagnostics <- .parallel_combine_diagnostics(res_list)
if (!is.null(diagnostics)) {
  diagnostics$estimator_runtime_summary <-
    .summarize_estimator_runtime(diagnostics$estimator_attempt_status)
  diagnostics$truncation_diagnostics_summary <-
    .summarize_truncation_diagnostics(
      diagnostics$truncation_diagnostics,
      attempt_status = diagnostics$estimator_attempt_status
    )
}


# ---- Recompute full truth and integrated performance -------------------------
# Shards only compute assigned replications.  To keep integrated truth and
# performance identical to the serial wrapper, compute truth once here with the
# original full scenario grid and original seed.
cat("[Parallel] Recomputing full truth and integrated performance.\n")
t_truth0 <- proc.time()[["elapsed"]]
truth_obj <- .compute_truth_for_unique_scenarios(unique(scenarios), params0, run_cfg)
truth_elapsed <- proc.time()[["elapsed"]] - t_truth0

performance <- summarize_performance(estimates_long, truth_obj$truth_effects)
performance <- .augment_performance_summary(
  performance,
  attempt_status = if (!is.null(diagnostics)) diagnostics$estimator_attempt_status else NULL,
  truncation_summary = if (!is.null(diagnostics)) diagnostics$truncation_diagnostics_summary else NULL
)
performance <- .parallel_add_performance_descriptors(performance, scenarios)

res <- list(
  truth_worldmeans = truth_obj$truth_worldmeans,
  truth_effects = truth_obj$truth_effects,
  estimates_long = estimates_long,
  worldmeans_estimates_long = worldmeans_estimates_long,
  performance = performance,
  diagnostics = diagnostics,
  shard_results = lapply(shard_results, function(x) {
    x$result <- NULL
    x
  })
)


# ---- Guardrail: enforce presence of requested estimators ---------------------
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


# ---- Write integrated outputs ------------------------------------------------
truth_effects_out <- order_effects_long(res$truth_effects)
estimates_long_out <- order_effects_long(.add_effect_truth_bias(res$estimates_long, res$truth_effects))
performance_out <- order_effects_long(res$performance)

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
write.csv(estimates_long_out,   file = file.path(output_dir, "estimates_long.csv"),   row.names = FALSE)
write.csv(worldmeans_long_out,  file = file.path(output_dir, "worldmeans_estimates_long.csv"), row.names = FALSE)
write.csv(performance_out,      file = file.path(output_dir, "performance_summary.csv"), row.names = FALSE)

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

.parallel_write_diagnostics(res$diagnostics, diag_dir = file.path(output_dir, "diagnostics"))

if (isTRUE(make_figures)) {
  tryCatch({
    create_figures(performance_out, out_pdf = file.path(output_dir, "performance_figures.pdf"))
  }, error = function(e) {
    warning("figure_generation_failed: performance_figures: ", conditionMessage(e), call. = FALSE)
  })
}

saveRDS(res, file.path(output_dir, "parallel", "integrated_result.rds"))

cat("\nDone. Integrated parallel output files written to:\n")
cat("  ", normalizePath(file.path(output_dir, "scenario_manifest.csv")), "\n")
cat("  ", normalizePath(file.path(output_dir, "truth_worldmeans.csv")), "\n")
cat("  ", normalizePath(file.path(output_dir, "truth_effects.csv")), "\n")
if (isTRUE(truth_figures_written)) {
  cat("  ", normalizePath(file.path(output_dir, "truth_figures.pdf")), "\n")
}
cat("  ", normalizePath(file.path(output_dir, "estimates_long.csv")), "\n")
cat("  ", normalizePath(file.path(output_dir, "worldmeans_estimates_long.csv")), "\n")
cat("  ", normalizePath(file.path(output_dir, "performance_summary.csv")), "\n")
cat("  ", normalizePath(file.path(output_dir, "parallel", "integrated_result.rds")), "\n")
cat(sprintf("[Parallel] Shard elapsed=%s | Full truth recompute=%s | Total=%s\n",
            .fmt_hms_parallel(parallel_elapsed),
            .fmt_hms_parallel(truth_elapsed),
            .fmt_hms_parallel(parallel_elapsed + truth_elapsed)))
