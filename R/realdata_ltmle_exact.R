################################################################################
# realdata_ltmle_exact.R
#
# Real-data adapter for the same full exact-EIF ltmle_exact engine used by the
# simulation layer.
################################################################################

.rd_numeric <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  as.numeric(factor(x))
}

.rd_ltmle_exact_input_from_prepared <- function(prepared) {
  spec <- prepared$spec
  wide <- prepared$wide

  baseline_vars <- spec$baseline_vars
  treatment_vars <- paste0("A", 0:(spec$T - 1L))
  mediator1_vars <- paste0("M1_", 0:(spec$T - 1L))
  mediator2_vars <- paste0("M2_", 0:(spec$T - 1L))
  L_blocks <- stats::setNames(
    lapply(spec$L_names, function(nm) paste0("L_", nm, "_", 0:(spec$T - 1L))),
    spec$L_names
  )
  L_lag_init <- spec$L_lag_init %||% list()
  L_lag_init_cols <- unlist(
    L_lag_init[vapply(L_lag_init, is.character, logical(1))],
    use.names = FALSE
  )
  outcome_vars <- spec$Y_cols %||% if ("Y_final" %in% names(wide)) "Y_final" else spec$Y_final_col %||% "Y_final"
  censoring_vars <- grep("^R_(visit|final)", names(wide), value = TRUE)

  required <- unique(c(
    baseline_vars,
    treatment_vars,
    mediator1_vars,
    mediator2_vars,
    unlist(L_blocks, use.names = FALSE),
    L_lag_init_cols,
    outcome_vars,
    censoring_vars
  ))
  assert_cols(wide, required, "prepared real-data input for ltmle_exact")

  list(
    data = wide,
    node_spec = list(
      baseline_vars = baseline_vars,
      treatment_vars = treatment_vars,
      mediator1_vars = mediator1_vars,
      mediator2_vars = mediator2_vars,
      L_blocks = L_blocks,
      L_order = paste0("L_", spec$L_names),
      L_lag_init = L_lag_init,
      outcome_vars = outcome_vars,
      censoring_vars = censoring_vars,
      T = spec$T
    )
  )
}

estimate_realdata_ltmle_exact <- function(prepared,
                                          control = default_esax_ltmle_exact_cfg(),
                                          seed = NULL) {
  input <- .rd_ltmle_exact_input_from_prepared(prepared)
  spec <- prepared$spec
  learner <- control$learner %||% if (identical(control$nuisance_engine %||% "glm", "superlearner")) "sl" else "glm"
  fit <- fit_ltmle_exact(
    dat = input$data,
    node_spec = input$node_spec,
    T = input$node_spec$T,
    reg_a = spec$reg_a,
    reg_as = spec$reg_as,
    learner = learner,
    sl_library = control$sl_library %||% c("SL.glm", "SL.mean"),
    Q_model = "correct",
    seed = seed %||% 20240505L,
    probability_bounds = control$probability_bounds %||% c(0.01, 0.99),
    truncation_enabled = control$truncation_enabled %||% TRUE,
    truncation_policy = control$truncation_policy %||% "quantile",
    truncation_quantile_lower = control$truncation_quantile_lower %||% 0.01,
    truncation_quantile_upper = control$truncation_quantile_upper %||% 0.99,
    truncation_target = control$truncation_target %||% "clever_covariate_H",
    ltmle_exact_density_ratio_mc_n = control$ltmle_exact_density_ratio_mc_n %||% 2000L,
    y_bounds_mode = "train_fold",
    treat_mech = spec$treat_mech %||% "baseline_rct",
    p_rct = spec$p_rct %||% 0.5,
    censoring_vars = input$node_spec$censoring_vars,
    censoring_mech = if (length(input$node_spec$censoring_vars)) "estimated" else "none",
    score_tolerance = control$score_tolerance %||% 1e-4,
    component_tolerance = control$component_tolerance %||% 1e-5,
    scaled_z_tolerance = control$scaled_z_tolerance %||% 3,
    V = control$fold_count %||% 5L,
    diagnostics_level = control$diagnostics_level %||% "summary",
    verbose = isTRUE(control$verbose %||% FALSE)
  )

  fixed <- compute_fixed_horizon_main_effects_from_means(fit$means)
  joint <- compute_joint_draw_effects_from_means(fit$means)
  effects <- list(
    fixed = fixed,
    joint = joint,
    time_resolved = data.frame()
  )

  out <- list(
    estimator = "ltmle_exact",
    method = "ltmle_exact",
    estimator_class = .ltmle_exact_class(),
    learner = learner,
    analysis_name = prepared$spec$analysis_name,
    spec = prepared$spec,
    means = fit$means,
    registry_means = fit$means,
    effects = effects,
    diagnostics = c(
      fit$diagnostics,
      list(
        node_spec = input$node_spec,
        eif = fit$diagnostics$component_eif_summary
      )
    ),
    targeting_diagnostics = fit$diagnostics$score_diagnostics,
    mc_diagnostics = fit$diagnostics$mc_integration_diagnostics,
    metadata = fit$metadata
  )
  class(out) <- c("realdata_ltmle_exact_result", "list")
  out
}

realdata_ltmle_exact_estimate <- estimate_realdata_ltmle_exact
