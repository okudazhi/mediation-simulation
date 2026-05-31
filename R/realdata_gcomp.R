################################################################################
# realdata_gcomp.R
#
# Real-data Monte Carlo g-computation for the fixed-horizon component means.
################################################################################

.realdata_gcomp_mc_diag <- function(component, y, n_subjects, mc_n) {
  data.frame(
    component = component,
    n_subjects = as.integer(n_subjects),
    mc_n = as.integer(mc_n),
    mean = mean(y, na.rm = TRUE),
    sd = stats::sd(y, na.rm = TRUE),
    mcse = stats::sd(y, na.rm = TRUE) / sqrt(length(y)),
    stringsAsFactors = FALSE
  )
}

fit_realdata_gcomp <- function(prepared,
                               control = default_esax_gcomp_cfg(),
                               seed = NULL) {
  input <- .rd_ltmle_exact_input_from_prepared(prepared)
  spec <- prepared$spec
  T <- input$node_spec$T
  learner <- control$learner %||% if (identical(control$nuisance_engine %||% "glm", "superlearner")) "sl" else "glm"
  mc_n <- as.integer(control$mc_B_main %||% 2000L)
  if (!is.finite(mc_n) || mc_n < 1L) .stop("g-computation mc_B_main must be a positive integer.")

  dat_wide <- .ltmle_exact_canonicalize_node_spec(input$data, input$node_spec, T)
  node_spec <- attr(dat_wide, "ltmle_exact_node_spec")
  censoring_vars <- .ltmle_exact_normalize_censoring_vars(
    node_spec$censoring_vars %||% input$node_spec$censoring_vars %||% NULL,
    T
  )
  censoring_adjust <- isTRUE(control$censoring_adjust %||% TRUE) && !is.null(censoring_vars)
  long <- .ltmle_exact_node_spec_to_long(dat_wide, T, node_spec)
  for (nm in as.character(node_spec$baseline_vars %||% character(0))) {
    if (nm %in% names(dat_wide) && !nm %in% names(long)) {
      long[[nm]] <- as.numeric(dat_wide[[nm]])[long$id]
    }
  }
  if (isTRUE(censoring_adjust)) {
    long <- .ltmle_exact_attach_censoring_to_long(
      long = long,
      dat_wide = dat_wide,
      T = T,
      censoring_vars = censoring_vars
    )
  }

  models <- .ltmle_exact_fit_node_models(
    long = long,
    T = T,
    learner = learner,
    sl_library = control$sl_library %||% c("SL.glm", "SL.mean"),
    Q_model = "correct",
    node_spec = node_spec,
    censoring_vars = censoring_vars,
    censoring_adjust = censoring_adjust
  )

  world_specs <- .ltmle_exact_world_spec(spec$reg_a, spec$reg_as, T)
  registry <- .ltmle_exact_component_registry(spec$reg_a, spec$reg_as, T)
  means <- setNames(rep(NA_real_, nrow(registry)), registry$component)
  mc_rows <- list()

  for (ii in seq_len(nrow(registry))) {
    comp <- registry$component[ii]
    sp <- world_specs[world_specs$component == comp, , drop = FALSE]
    histories <- .ltmle_exact_draw_component_histories(
      dat_wide = dat_wide,
      T = T,
      spec = sp,
      models = models,
      mc_n = mc_n,
      seed = (seed %||% 20240505L) + ii * 1009L,
      node_spec = node_spec
    )
    y <- histories$outcome_process[[T]]$row$Y
    means[comp] <- mean(y, na.rm = TRUE)
    mc_rows[[length(mc_rows) + 1L]] <- .realdata_gcomp_mc_diag(comp, y, nrow(dat_wide), mc_n)
  }

  fixed <- compute_fixed_horizon_main_effects_from_means(as.list(means))
  joint <- compute_joint_draw_effects_from_means(as.list(means))
  effects <- list(
    fixed = fixed,
    joint = joint,
    time_resolved = data.frame()
  )

  out <- list(
    estimator = "gcomp",
    analysis_name = spec$analysis_name,
    spec = spec,
    means = as.list(means),
    registry_means = as.list(means),
    effects = effects,
    mc_diagnostics = do.call(rbind, mc_rows),
    targeting_diagnostics = NULL,
    diagnostics = list(
      component_registry = registry,
      mediator_density_engine = models$mediator_density_engine %||% "gaussian_location_scale",
      censoring_adjusted = isTRUE(censoring_adjust),
      censoring_vars = censoring_vars
    ),
    model_fits = models,
    control = control,
    registry = list(fixed = registry, time = list())
  )
  class(out) <- c("realdata_gcomp_result", "list")
  out
}

realdata_gcomp_estimate <- function(prepared,
                                    control = default_esax_gcomp_cfg(),
                                    seed = NULL) {
  fit_realdata_gcomp(prepared = prepared, control = control, seed = seed)
}
