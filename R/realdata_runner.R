################################################################################
# realdata_runner.R
#
# Orchestration helpers for the real-data empirical comparison.
################################################################################

.write_realdata_figures <- function(result, outdir, prefix) {
  fig_dir <- file.path(outdir, "figures")
  ensure_dir(fig_dir)

  if (exists("plot_weight_diagnostics", mode = "function")) {
    plot_weight_diagnostics(
      result,
      file.path(fig_dir, paste0(prefix, "_weight_diagnostics.pdf"))
    )
  }
  if (exists("plot_mc_diagnostics", mode = "function")) {
    plot_mc_diagnostics(
      result,
      file.path(fig_dir, paste0(prefix, "_mc_diagnostics.pdf"))
    )
  }
  if (exists("plot_time_resolved_decomposition", mode = "function")) {
    plot_time_resolved_decomposition(
      result,
      file.path(fig_dir, paste0(prefix, "_time_resolved_decomposition.pdf"))
    )
  }

  invisible(fig_dir)
}

.run_single_realdata_estimator <- function(prepared,
                                           estimator_name,
                                           estimator_fn,
                                           control,
                                           bootstrap_cfg = NULL,
                                           outdir = NULL,
                                           seed = NULL,
                                           save_diagnostics = TRUE,
                                           progress = TRUE) {
  if (isTRUE(progress)) {
    message("[", prepared$spec$analysis_name, "] fitting ", estimator_name, " ...")
  }

  res <- estimator_fn(prepared, control = control, seed = seed)

  if (!is.null(bootstrap_cfg) && isTRUE(bootstrap_cfg$enabled %||% TRUE) &&
      isTRUE((bootstrap_cfg$B %||% 0L) > 0L)) {
    boot <- bootstrap_realdata_estimator(
      prepared = prepared,
      estimator_fn = estimator_fn,
      control = control,
      B = bootstrap_cfg$B,
      conf_level = bootstrap_cfg$conf_level %||% 0.95,
      seed_base = bootstrap_cfg$seed_base %||% 20260424L,
      checkpoint_dir = file.path(outdir, "bootstrap"),
      prefix = estimator_name,
      save_every = bootstrap_cfg$save_every %||%
        default_esax_bootstrap_cfg(enabled = TRUE)$save_every,
      resume = isTRUE(bootstrap_cfg$resume %||% TRUE),
      progress = progress
    )
    res$bootstrap <- boot
  }
  res$output_name <- estimator_name
  res$internal_estimator <- res$estimator %||% estimator_name

  if (!is.null(outdir)) {
    ensure_dir(outdir)
    save_rds_safe(res, file.path(outdir, paste0(estimator_name, "_result.rds")))
    write_analysis_tables(res, outdir = outdir,
                          prefix = estimator_name,
                          digits = 3,
                          write_tex = TRUE)

    if (isTRUE(save_diagnostics)) {
      if (!is.null(res$mc_diagnostics)) {
        write_csv_safe(res$mc_diagnostics, file.path(outdir, paste0(estimator_name, "_mc_diagnostics.csv")))
      }
      if (!is.null(res$targeting_diagnostics)) {
        write_csv_safe(res$targeting_diagnostics, file.path(outdir, paste0(estimator_name, "_targeting_diagnostics.csv")))
      }
      if (!identical(res$estimator, "ltmle_exact") &&
          !is.null(res$weight_models) &&
          !is.null(res$weight_models$subject_weights)) {
        write_csv_safe(res$weight_models$subject_weights, file.path(outdir, paste0(estimator_name, "_subject_weights.csv")))
      }
      if (!identical(res$estimator, "ltmle_exact") &&
          !is.null(res$weight_models) &&
          !is.null(res$weight_models$long)) {
        write_csv_safe(res$weight_models$long, file.path(outdir, paste0(estimator_name, "_long_weights.csv")))
      }
      if (!is.null(res$diagnostics$weight_summary)) {
        write_csv_safe(res$diagnostics$weight_summary,
                       file.path(outdir, paste0(estimator_name, "_weight_summary.csv")))
      }
      if (!is.null(res$diagnostics$visit_weight_summary)) {
        write_csv_safe(res$diagnostics$visit_weight_summary,
                       file.path(outdir, paste0(estimator_name, "_visit_weight_summary.csv")))
      }
      if (!is.null(res$diagnostics$eif)) {
        write_csv_safe(res$diagnostics$eif, file.path(outdir, paste0(estimator_name, "_eif_diagnostics.csv")))
      }
    }
    .write_realdata_figures(res, outdir = outdir, prefix = estimator_name)
  }

  res
}

.realdata_control_for_spec <- function(spec, ltmle_control, msm_control = NULL, gcomp_control = NULL) {
  method <- spec$method
  eng <- spec$nuisance_engine %||% "glm"
  eng <- if (eng %in% c("sl", "superlearner")) "superlearner" else "glm"

  if (identical(method, "ltmle_exact")) {
    ctl <- ltmle_control %||% default_esax_ltmle_exact_cfg(nuisance_engine = eng)
    ctl$nuisance_engine <- eng
    ctl$learner <- if (identical(eng, "superlearner")) "sl" else "glm"
    return(ctl)
  }

  if (identical(method, "msm_ipw")) {
    ctl <- msm_control %||% default_esax_msm_cfg(nuisance_engine = eng)
    ctl$nuisance_engine <- eng
    return(ctl)
  }

  if (identical(method, "gcomp")) {
    ctl <- gcomp_control %||% default_esax_gcomp_cfg(nuisance_engine = eng)
    ctl$nuisance_engine <- eng
    ctl$learner <- if (identical(eng, "superlearner")) "sl" else "glm"
    return(ctl)
  }

  .stop("Unknown real-data method: ", method)
}

.realdata_estimator_for_spec <- function(spec) {
  if (identical(spec$method, "ltmle_exact")) return(realdata_ltmle_exact_estimate)
  if (identical(spec$method, "msm_ipw")) return(realdata_msm_ipw_estimate)
  if (identical(spec$method, "gcomp")) return(realdata_gcomp_estimate)
  .stop("Unknown real-data method: ", spec$method)
}

run_realdata_pipeline <- function(dat_raw,
                                  analyses,
                                  var_spec,
                                  T,
                                  output_dir,
                                  method_specs = default_method_specs(),
                                  ltmle_control = default_esax_ltmle_exact_cfg(),
                                  msm_control = default_esax_msm_cfg(),
                                  gcomp_control = default_esax_gcomp_cfg(),
                                  bootstrap_cfg = NULL,
                                  save_diagnostics = TRUE,
                                  progress = TRUE) {
  ensure_dir(output_dir)
  results <- list()

  for (nm in names(analyses)) {
    cfg <- analyses[[nm]]
    if (!identical(cfg$analysis_type %||% "arm", "arm")) {
      .stop("Only arm_comparison is supported.")
    }

    prepared <- prepare_realdata_analysis(
      dat_raw = dat_raw,
      var_spec = var_spec,
      analysis_cfg = cfg,
      T = T,
      analysis_name = cfg$name %||% nm
    )
    analysis_dir <- file.path(output_dir, prepared$spec$analysis_name)
    ensure_dir(analysis_dir)
    save_rds_safe(prepared$spec, file.path(analysis_dir, "analysis_spec.rds"))
    save_rds_safe(prepared$meta, file.path(analysis_dir, "analysis_meta.rds"))

    res_an <- list(spec = prepared$spec, meta = prepared$meta)

    for (ms_nm in names(method_specs)) {
      ms <- method_specs[[ms_nm]]
      ctl <- .realdata_control_for_spec(ms, ltmle_control, msm_control, gcomp_control)
      est_fn <- .realdata_estimator_for_spec(ms)
      estimator_name <- ms$name %||% paste(ms$method, ms$nuisance_engine, sep = "_")
      res_an[[estimator_name]] <- .run_single_realdata_estimator(
        prepared = prepared,
        estimator_name = estimator_name,
        estimator_fn = est_fn,
        control = ctl,
        bootstrap_cfg = bootstrap_cfg,
        outdir = file.path(analysis_dir, estimator_name),
        seed = 20260424L + length(res_an),
        save_diagnostics = save_diagnostics,
        progress = progress
      )
    }

    results[[prepared$spec$analysis_name]] <- res_an
  }

  save_rds_safe(results, file.path(output_dir, "realdata_results.rds"))
  results
}
