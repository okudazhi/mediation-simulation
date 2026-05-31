################################################################################
# estimator_ltmle_exact.R
#
# Full exact-EIF targeted minimum loss-based estimator for the nine fixed-horizon
# component means used in this project.  The public entry point intentionally has
# one production path only: fit the component-specific sequential regressions,
# target the outcome, post-mediator covariate transition, and stochastic mediator
# law continuation equations, then return the plug-in means and diagnostics.
################################################################################

.ltmle_exact_class <- function() {
  "exact_eif_targeted_minimum_loss_estimator"
}

.ltmle_exact_log <- function(verbose, ...) {
  if (isTRUE(verbose)) {
    message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  }
}

.ltmle_exact_clamp01 <- function(x, eps = 1e-6) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  pmin(pmax(x, eps), 1 - eps)
}

.ltmle_exact_bounds <- function(y, y_bounds_mode = c("train_fold", "fixed"),
                                y_bounds = NULL) {
  y_bounds_mode <- match.arg(y_bounds_mode)
  yy <- as.numeric(y)
  if (y_bounds_mode == "fixed") {
    if (is.null(y_bounds) || length(y_bounds) != 2L) {
      .stop("y_bounds must be length 2 when y_bounds_mode='fixed'.")
    }
    a <- as.numeric(y_bounds[1L])
    b <- as.numeric(y_bounds[2L])
  } else {
    a <- min(yy, na.rm = TRUE)
    b <- max(yy, na.rm = TRUE)
  }
  if (!is.finite(a) || !is.finite(b)) {
    .stop("Outcome bounds are not finite. Check pseudo-outcome values or specify fixed y_bounds.")
  }
  if ((b - a) <= 1e-8) {
    center <- mean(yy, na.rm = TRUE)
    if (!is.finite(center)) center <- 0
    half_width <- max(1, abs(center) * 0.01, 1e-4)
    a <- center - half_width
    b <- center + half_width
  }
  list(a = a, b = b, mode = y_bounds_mode)
}

.ltmle_exact_scale01 <- function(y, bounds, eps = 1e-6) {
  .ltmle_exact_clamp01((as.numeric(y) - bounds$a) / (bounds$b - bounds$a), eps = eps)
}

.ltmle_exact_unscale01 <- function(y01, bounds) {
  bounds$a + (bounds$b - bounds$a) * as.numeric(y01)
}

.ltmle_exact_formula <- function(lhs, rhs) {
  stats::as.formula(paste(lhs, "~", paste(unique(rhs), collapse = " + ")))
}

.ltmle_exact_model_data <- function(df, cols) {
  cols <- intersect(unique(cols), names(df))
  out <- df[, cols, drop = FALSE]
  for (nm in names(out)) out[[nm]] <- as.numeric(out[[nm]])
  out
}

.ltmle_exact_q_covs_y <- function(Q_model = c("correct", "wrong")) {
  Q_model <- match.arg(Q_model)
  x <- c("W1", "W2", "Y0", "Y_lag", "L_lag", "M1_lag", "M2_lag",
         "A", "M1", "M2", "L")
  if (Q_model == "correct") {
    x <- c(x, "AM1", "AM2", "MM12")
  }
  unique(x)
}

.ltmle_exact_q_covs_l <- function(Q_model = c("correct", "wrong")) {
  Q_model <- match.arg(Q_model)
  x <- c("W1", "W2", "Y0", "Y_lag", "L_lag", "M1_lag", "M2_lag",
         "A", "M1", "M2")
  if (Q_model == "correct") {
    x <- c(x, "AM1", "AM2", "MM12")
  }
  unique(x)
}

.ltmle_exact_add_terms <- function(df) {
  if (all(c("A", "M1") %in% names(df))) df$AM1 <- df$A * df$M1
  if (all(c("A", "M2") %in% names(df))) df$AM2 <- df$A * df$M2
  if (all(c("M1", "M2") %in% names(df))) df$MM12 <- df$M1 * df$M2
  df
}

.ltmle_exact_parent_marker_name <- function(task) {
  task <- as.list(task)
  raw <- paste(
    task$component %||% "component",
    task$task_id %||% "task",
    task$t %||% "t",
    task$node %||% "node",
    sep = "_"
  )
  paste0(".ltmle_parent_", gsub("[^A-Za-z0-9_]+", "_", raw))
}

.ltmle_exact_particle_marker_cols <- function(...) {
  frames <- list(...)
  frames <- frames[vapply(frames, function(x) is.data.frame(x) && nrow(x) >= 0L, logical(1))]
  if (!length(frames)) return(character())
  common <- Reduce(intersect, lapply(frames, names))
  grep("^\\.ltmle_parent_", common, value = TRUE)
}

.ltmle_exact_branch_group_key <- function(particles, marker_cols = character()) {
  if (!nrow(particles)) return(character())
  marker_cols <- intersect(marker_cols, names(particles))
  pieces <- c(
    list(as.character(as.integer(particles$id0))),
    lapply(marker_cols, function(nm) as.character(particles[[nm]]))
  )
  do.call(paste, c(pieces, sep = "\r"))
}

.ltmle_exact_unique_particle_rows <- function(particles, exclude = character()) {
  if (!nrow(particles)) return(particles)
  cols <- setdiff(names(particles), unique(c("id", ".particle_id", ".orig_index", exclude)))
  if (!length(cols)) return(particles[!duplicated(seq_len(nrow(particles))), , drop = FALSE])
  key <- do.call(paste, c(lapply(cols, function(nm) as.character(particles[[nm]])), sep = "\r"))
  particles[!duplicated(key), , drop = FALSE]
}

.ltmle_exact_set_parent_marker <- function(branch_state, marker) {
  for (key in .ltmle_exact_branch_keys()) {
    branch_state[[key]]$particles[[marker]] <- seq_len(nrow(branch_state[[key]]$particles))
  }
  branch_state
}

.ltmle_exact_integrate_child_particle_q <- function(child,
                                                    marker,
                                                    n0,
                                                    ng,
                                                    grid_w,
                                                    task_id) {
  expected_n <- as.integer(n0) * as.integer(ng)
  q <- as.numeric(child$particle_q)
  frame <- child$particle_frame
  cell_q <- NULL
  if (is.data.frame(frame) && marker %in% names(frame) && length(q) == nrow(frame)) {
    marker_id <- as.integer(frame[[marker]])
    keep <- is.finite(marker_id) & marker_id >= 1L & marker_id <= expected_n & is.finite(q)
    if (!all(keep)) {
      .stop("Child branch evaluation returned non-finite or unmapped particle values for task=", task_id)
    }
    w <- as.numeric(frame$.branch_weight %||% rep(1, nrow(frame)))
    w[!is.finite(w) | w < 0] <- 0
    num <- rowsum(q * w, group = marker_id, reorder = TRUE)
    den <- rowsum(w, group = marker_id, reorder = TRUE)
    idx <- as.integer(rownames(num))
    cell_q <- rep(NA_real_, expected_n)
    den_v <- as.numeric(den[, 1L])
    if (any(!is.finite(den_v) | den_v <= 0)) {
      .stop("Child branch evaluation has zero product weight for task=", task_id)
    }
    cell_q[idx] <- as.numeric(num[, 1L]) / den_v
  } else if (length(q) == expected_n) {
    cell_q <- q
  } else {
    .stop(
      "Child branch evaluation particle count mismatch for task=", task_id,
      ": got ", length(q), ", expected ", expected_n,
      ". Product-join children must carry parent-particle markers."
    )
  }
  if (length(cell_q) != expected_n || any(!is.finite(cell_q))) {
    .stop("Child branch evaluation did not collapse to finite integration cells for task=", task_id)
  }
  as.numeric(matrix(cell_q, nrow = n0, ncol = ng) %*% grid_w)
}

.ltmle_exact_set_A <- function(df, A) {
  df$A <- as.numeric(A)
  .ltmle_exact_add_terms(df)
}

.ltmle_exact_component_registry <- function(reg_a, reg_as, T) {
  if (!exists("ltmle_exact_component_registry", mode = "function")) {
    .stop("ltmle_exact_component_registry() must be available before estimator_ltmle_exact.R is used.")
  }
  reg <- ltmle_exact_component_registry(reg_a, reg_as, T)
  reg$outer_label <- reg$outcome_regimen
  reg$m1_label <- ifelse(reg$first_mediator_regimen == "outcome",
                         reg$outer_label, reg$first_mediator_regimen)
  reg$m2_label <- ifelse(reg$second_mediator_regimen == "outcome",
                         reg$outer_label, reg$second_mediator_regimen)
  reg
}

.ltmle_exact_world_spec <- function(reg_a, reg_as, T) {
  reg_a <- normalize_regimen(reg_a, T, "reg_a")
  reg_as <- normalize_regimen(reg_as, T, "reg_as")
  reg_by_label <- function(lbl, outer) {
    if (identical(lbl, "a")) return(reg_a)
    if (identical(lbl, "as")) return(reg_as)
    if (identical(lbl, "outer")) return(outer)
    .stop("Unknown regimen label: ", lbl)
  }
  reg <- .ltmle_exact_component_registry(reg_a, reg_as, T)
  out <- vector("list", nrow(reg))
  for (ii in seq_len(nrow(reg))) {
    outer <- reg_by_label(reg$outer_label[ii], reg_a)
    m1 <- reg_by_label(reg$m1_label[ii], outer)
    m2 <- reg_by_label(reg$m2_label[ii], outer)
    out[[ii]] <- data.frame(
      component_id = reg$component_id[ii],
      component = reg$component[ii],
      world_type = reg$world_type[ii],
      t = seq_len(T),
      outer_A = outer,
      m1_A = m1,
      m2_A = m2,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

.ltmle_exact_factor_tasks <- function(reg_a, reg_as, T, node_spec = NULL) {
  if (!exists("ltmle_exact_factor_task_registry", mode = "function")) {
    .stop(
      "ltmle_exact_factor_task_registry() must be available before ",
      "estimator_ltmle_exact.R is used. The full exact efficient influence ",
      "function-based TMLE requires the explicit factor-task registry."
    )
  }

  out <- ltmle_exact_factor_task_registry(reg_a = reg_a, reg_as = reg_as, T = T, node_spec = node_spec)

  required <- c(
    "component", "t", "node", "process_type", "task_id",
    "observed_pseudooutcome_source_task_id",
    "assigned_treatment_regimen_label",
    "conditioning_history_type",
    "uses_current_mediator_density_ratio",
    "uses_current_censoring_factor",
    "mediator_role"
  )
  missing <- setdiff(required, names(out))
  if (length(missing)) {
    .stop(
      "ltmle_exact_factor_task_registry() returned an incomplete registry. Missing: ",
      paste(missing, collapse = ", ")
    )
  }

  out
}

.ltmle_exact_L_nodes <- function(node_spec = NULL) {
  if (is.null(node_spec)) return("L")
  L_order <- as.character(node_spec$L_order %||% character(0))
  if (length(L_order)) return(L_order)
  L_blocks <- node_spec$L_blocks %||% NULL
  if (is.null(L_blocks) || !length(L_blocks)) return("L")
  nm <- names(L_blocks)
  if (is.null(nm) || any(!nzchar(nm))) nm <- paste0("L", seq_along(L_blocks))
  ifelse(startsWith(nm, "L_") | identical(nm, "L"), nm, paste0("L_", nm))
}

.ltmle_exact_multi_L <- function(node_spec = NULL) {
  nodes <- .ltmle_exact_L_nodes(node_spec)
  length(nodes) > 1L || !identical(nodes, "L")
}

.ltmle_exact_lag_name <- function(node) {
  if (identical(node, "L")) return("L_lag")
  paste0(node, "_lag")
}

.ltmle_exact_outcome_type <- function(node_spec = NULL) {
  node_spec$outcome_type %||% "longitudinal"
}

.ltmle_exact_prior_history_covariates <- function(node_spec = NULL) {
  out <- c(vapply(.ltmle_exact_L_nodes(node_spec), .ltmle_exact_lag_name, character(1)),
           "M1_lag", "M2_lag")
  if (!identical(.ltmle_exact_outcome_type(node_spec), "terminal_only")) out <- c("Y_lag", out)
  unique(out)
}

.ltmle_exact_visit_censoring_covariates <- function(node_spec = NULL) {
  unique(c("W1", "W2", "Y0",
           as.character(node_spec$baseline_vars %||% character(0)),
           .ltmle_exact_prior_history_covariates(node_spec)))
}

.ltmle_exact_history_covs_for_L_node <- function(L_node, Q_model = c("correct", "wrong"),
                                                 node_spec = NULL) {
  Q_model <- match.arg(Q_model)
  L_nodes <- .ltmle_exact_L_nodes(node_spec)
  pos <- match(L_node, L_nodes)
  preceding_L <- if (is.na(pos) || pos <= 1L) character(0) else L_nodes[seq_len(pos - 1L)]
  baseline_vars <- as.character(node_spec$baseline_vars %||% character(0))
  base <- unique(c("W1", "W2", "Y0", baseline_vars,
                   .ltmle_exact_prior_history_covariates(node_spec), "A", "M1", "M2",
                   preceding_L))
  unique(base)
}

.ltmle_exact_history_covs <- function(node, Q_model = c("correct", "wrong"), node_spec = NULL) {
  Q_model <- match.arg(Q_model)
  if (node %in% .ltmle_exact_L_nodes(node_spec)) {
    return(.ltmle_exact_history_covs_for_L_node(node, Q_model, node_spec))
  }
  baseline_vars <- as.character(node_spec$baseline_vars %||% character(0))
  base <- unique(c("W1", "W2", "Y0", baseline_vars,
                   .ltmle_exact_prior_history_covariates(node_spec), "A"))
  if (node == "M1") return(base)
  if (node == "M2") {
    if (Q_model == "correct") return(unique(c(base, "M1", "AM1")))
    return(unique(base))
  }
  if (node == "L") {
    return(.ltmle_exact_history_covs_for_L_node("L", Q_model, node_spec))
  }
  if (node == "Y") {
    L_nodes <- .ltmle_exact_L_nodes(node_spec)
    out <- c(base, "M1", "M2", L_nodes)
    if (Q_model == "correct") {
      out <- c(out, "AM1", "AM2", "MM12")
    }
    return(unique(out))
  }
  .stop("Unknown node: ", node)
}

.ltmle_exact_continuation_Q_interaction_covs <- function(Q_model = c("correct", "wrong")) {
  Q_model <- match.arg(Q_model)
  if (identical(Q_model, "correct")) c("AM1", "AM2", "MM12") else character(0)
}

.ltmle_exact_L_transition_continuation_Q_covs <- function(L_node, Q_model = c("correct", "wrong"),
                                                          node_spec = NULL) {
  Q_model <- match.arg(Q_model)
  unique(c(
    .ltmle_exact_history_covs_for_L_node(L_node, Q_model, node_spec),
    .ltmle_exact_continuation_Q_interaction_covs(Q_model)
  ))
}

.ltmle_exact_virtual_mixed_L_continuation_Q_covs <- function(L_node, Q_model = c("correct", "wrong"),
                                                             node_spec = NULL) {
  .ltmle_exact_L_transition_continuation_Q_covs(L_node, Q_model, node_spec)
}

.ltmle_exact_continuation_Q_covs <- function(node, Q_model = c("correct", "wrong"),
                                             node_spec = NULL) {
  Q_model <- match.arg(Q_model)
  if (node %in% .ltmle_exact_L_nodes(node_spec)) {
    return(.ltmle_exact_L_transition_continuation_Q_covs(node, Q_model, node_spec))
  }
  covs <- .ltmle_exact_history_covs(node, Q_model, node_spec)
  if (identical(node, "Y")) {
    covs <- c(covs, .ltmle_exact_continuation_Q_interaction_covs(Q_model))
  }
  unique(covs)
}

.ltmle_exact_virtual_mixed_covariate_set_role <- function(task, node_spec = NULL) {
  task <- as.list(task)
  explicit <- as.character(task$virtual_mixed_covariate_set_role %||% NA_character_)
  allowed <- c(
    "post_mediator_L_transition",
    "pre_mediator_continuation",
    "standard_node_covariates"
  )
  if (!is.na(explicit) && nzchar(explicit)) {
    if (!explicit %in% allowed) {
      .stop("Unknown virtual_mixed_covariate_set_role: ", explicit)
    }
    return(explicit)
  }
  if (!.ltmle_exact_is_virtual_mixed_task(task)) {
    return("standard_node_covariates")
  }

  target_label <- as.character(
    task$virtual_mixed_continuation_target_label %||% NA_character_
  )
  task_node <- as.character(task$node %||% NA_character_)
  if (identical(target_label, "mixed_outer_LY_history_with_mediator_law_path") &&
      task_node %in% .ltmle_exact_L_nodes(node_spec)) {
    return("post_mediator_L_transition")
  }
  if (identical(target_label, "mixed_outer_LY_history_with_mediator_law_path")) {
    return("pre_mediator_continuation")
  }
  "standard_node_covariates"
}

.ltmle_exact_virtual_mixed_needs_pre_mediator_covariates <- function(task, node_spec = NULL) {
  identical(
    .ltmle_exact_virtual_mixed_covariate_set_role(task, node_spec),
    "pre_mediator_continuation"
  )
}

.ltmle_exact_virtual_mixed_can_use_downstream_cached_source <- function(task) {
  task <- as.list(task)
  if (!.ltmle_exact_is_virtual_mixed_task(task) ||
      !identical(as.character(task$world_type %||% NA_character_), "separate")) {
    return(FALSE)
  }
  outer_regimen <- as.character(task$source_boundary_outer_regimen %||% NA_character_)
  m1_regimen <- as.character(task$source_boundary_m1_regimen %||% NA_character_)
  m2_regimen <- as.character(task$source_boundary_m2_regimen %||% NA_character_)
  if (!is.na(outer_regimen) && nzchar(outer_regimen) &&
      !is.na(m1_regimen) && nzchar(m1_regimen) &&
      !is.na(m2_regimen) && nzchar(m2_regimen)) {
    return(identical(outer_regimen, m1_regimen) && identical(m1_regimen, m2_regimen))
  }
  FALSE
}

.ltmle_exact_virtual_mixed_direct_continuation_status <- function(task) {
  task <- as.list(task)
  can_represent <- isTRUE(.ltmle_exact_virtual_mixed_can_use_downstream_cached_source(task))
  list(
    direct_continuation_allowed = can_represent,
    downstream_source_can_represent_virtual_target = can_represent,
    dedicated_virtual_Q_target_used = .ltmle_exact_is_virtual_mixed_task(task) && !can_represent
  )
}

.ltmle_exact_task_history_covs <- function(task, Q_model = c("correct", "wrong"), node_spec = NULL) {
  Q_model <- match.arg(Q_model)
  task <- as.list(task)
  virtual_mixed_covariate_set_role <-
    .ltmle_exact_virtual_mixed_covariate_set_role(task, node_spec)
  if (identical(virtual_mixed_covariate_set_role, "post_mediator_L_transition")) {
    return(.ltmle_exact_virtual_mixed_L_continuation_Q_covs(
      as.character(task$node %||% .ltmle_exact_L_nodes(node_spec)[1L]),
      Q_model,
      node_spec
    ))
  }
  if (identical(virtual_mixed_covariate_set_role, "pre_mediator_continuation")) {
    return(.ltmle_exact_history_covs("M1", Q_model, node_spec))
  }
  .ltmle_exact_continuation_Q_covs(task$node, Q_model, node_spec)
}

.ltmle_exact_process_type <- function(task) {
  task$process_type
}

.ltmle_exact_is_virtual_mixed_task <- function(task) {
  identical(as.character(.ltmle_exact_process_type(as.list(task))), "virtual_mixed_continuation_task")
}

.ltmle_exact_continuation_target_label <- function(task) {
  task <- as.list(task)
  if (.ltmle_exact_is_virtual_mixed_task(task)) {
    return(as.character(task$virtual_mixed_continuation_target_label %||%
                          "mixed_outer_LY_history_with_mediator_law_path"))
  }
  as.character(task$conditioning_history_type %||% "standard_factor_task")
}

.ltmle_exact_semantic_cache_key <- function(task) {
  task <- as.list(task)
  fields <- c(
    component = as.character(task$component %||% NA_character_),
    task_id = as.character(task$task_id %||% NA_character_),
    continuation_target_label = .ltmle_exact_continuation_target_label(task),
    outer_regimen = as.character(task$source_boundary_outer_regimen %||% NA_character_),
    m1_regimen = as.character(task$source_boundary_m1_regimen %||% NA_character_),
    m2_regimen = as.character(task$source_boundary_m2_regimen %||% NA_character_),
    outcome_history_state =
      as.character(task$source_boundary_outcome_history_state %||% NA_character_),
    m1_history_state =
      as.character(task$source_boundary_m1_history_state %||% NA_character_),
    m2_history_state =
      as.character(task$source_boundary_m2_history_state %||% NA_character_),
    auxiliary_mediator_history_state =
      as.character(task$source_boundary_auxiliary_mediator_history_state %||% NA_character_)
  )
  fields[is.na(fields) | !nzchar(fields)] <- "not_available"
  paste(paste(names(fields), fields, sep = "="), collapse = "|")
}

.ltmle_exact_is_mixed_boundary_type <- function(boundary_type) {
  as.character(boundary_type) %in% c("cross_regimen_mixed", "virtual_mixed_continuation")
}

.ltmle_exact_fit_nuis <- function(formula, data, family, learner, sl_library,
                                  component, weights = NULL, fold_id = NULL) {
  engine <- if (learner == "sl") "superlearner" else "glm"
  fit <- fit_nuisance_model(
    formula = formula,
    data = data,
    family = family,
    engine = engine,
    sl_library = sl_library,
    component = component,
    weights = weights
  )
  if (identical(learner, "sl") && !is.null(fold_id)) {
    fold_id <- as.integer(fold_id)
    if (length(fold_id) == nrow(data) && length(unique(fold_id)) >= 2L) {
      cf_fits <- list()
      for (ff in sort(unique(fold_id))) {
        train <- fold_id != ff
        if (!any(train)) next
        cf_fits[[as.character(ff)]] <- fit_nuisance_model(
          formula = formula,
          data = data[train, , drop = FALSE],
          family = family,
          engine = engine,
          sl_library = sl_library,
          component = paste0(component, " fold=", ff),
          weights = if (is.null(weights)) NULL else as.numeric(weights)[train]
        )
      }
      if (length(cf_fits) >= 2L) {
        fit$crossfit_fits <- cf_fits
        fit$crossfit_used <- TRUE
      }
    }
  }
  fit
}

.ltmle_exact_predict_chunk_size <- function() {
  100000L
}

.ltmle_exact_normalize_density_ratio_mc_n <- function(ltmle_exact_density_ratio_mc_n) {
  ltmle_exact_density_ratio_mc_n <- as.integer(ltmle_exact_density_ratio_mc_n)
  if (!is.finite(ltmle_exact_density_ratio_mc_n) || ltmle_exact_density_ratio_mc_n < 2L) {
    .stop("ltmle_exact_density_ratio_mc_n must be an integer >= 2.")
  }
  ltmle_exact_density_ratio_mc_n
}

.ltmle_exact_prediction_data <- function(newdata, required, keep_fold_id = FALSE) {
  cols <- unique(c(required, if (keep_fold_id && ".fold_id" %in% names(newdata)) ".fold_id"))
  out <- newdata[, cols, drop = FALSE]
  attr(out, "ltmle_exact_long") <- NULL
  attr(out, "ltmle_exact_row_source") <- NULL
  out
}

.ltmle_exact_can_fast_predict_glm <- function(fit, pred_data, required) {
  if (!identical(fit$type, "glm")) return(FALSE)
  if (!length(required)) return(FALSE)
  if (!all(required %in% names(pred_data))) return(FALSE)
  if (!all(vapply(pred_data[, required, drop = FALSE], is.numeric, logical(1)))) return(FALSE)

  terms_obj <- stats::terms(fit$fit)
  if (attr(terms_obj, "offset") %||% FALSE) return(FALSE)
  term_labels <- attr(terms_obj, "term.labels") %||% character(0)
  if (!setequal(term_labels, required)) return(FALSE)

  TRUE
}

.ltmle_exact_predict_glm_fast <- function(fit, pred_data, type = "response", eps = 1e-8) {
  beta <- stats::coef(fit$fit)
  beta[!is.finite(beta)] <- 0

  intercept <- if ("(Intercept)" %in% names(beta)) as.numeric(beta["(Intercept)"]) else 0
  vars <- setdiff(names(beta), "(Intercept)")

  missing <- setdiff(vars, names(pred_data))
  if (length(missing)) {
    .stop("Fast glm prediction missing columns: ", paste(missing, collapse = ", "))
  }

  eta <- rep(intercept, nrow(pred_data))
  for (v in vars) {
    eta <- eta + as.numeric(beta[v]) * as.numeric(pred_data[[v]])
  }

  fam <- .nuisance_family_name(fit$family)
  if (identical(fam, "binomial")) {
    return(.ltmle_exact_clamp01(stats::plogis(eta), eps = eps))
  }
  if (identical(type, "response")) {
    return(as.numeric(eta))
  }
  as.numeric(eta)
}

.ltmle_exact_predict_nuis <- function(fit, newdata, type = "response", eps = 1e-8) {
  n <- nrow(newdata)
  if (!n) return(numeric(0))

  required <- fit$rhs_names %||% character(0)
  missing <- setdiff(required, names(newdata))
  if (length(missing)) {
    .stop(
      "ltmle_exact strict prediction failed for ", fit$component,
      ": missing columns: ", paste(missing, collapse = ", ")
    )
  }

  if (!is.null(fit$crossfit_fits) && !(".fold_id" %in% names(newdata))) {
    .stop("ltmle_exact cross-fit prediction requires .fold_id for ", fit$component, ".")
  }
  if (!is.null(fit$crossfit_fits) && ".fold_id" %in% names(newdata)) {
    out <- rep(NA_real_, n)
    fold_data <- .ltmle_exact_prediction_data(newdata, required, keep_fold_id = TRUE)
    for (ff in names(fit$crossfit_fits)) {
      idx <- as.character(fold_data$.fold_id) == ff
      if (any(idx)) {
        out[idx] <- .ltmle_exact_predict_nuis(fit$crossfit_fits[[ff]], fold_data[idx, , drop = FALSE],
                                              type = type, eps = eps)
      }
    }
    if (length(out) == n && all(is.finite(out))) return(out)
    .stop("ltmle_exact cross-fit prediction failed for ", fit$component, ".")
  }

  pred_data <- .ltmle_exact_prediction_data(newdata, required, keep_fold_id = FALSE)

  if (length(required)) {
    cc <- stats::complete.cases(pred_data[, required, drop = FALSE])
    if (!all(cc)) {
      .stop(
        "ltmle_exact strict prediction failed for ", fit$component,
        ": incomplete prediction rows detected."
      )
    }
  }

  if (.ltmle_exact_can_fast_predict_glm(fit, pred_data, required)) {
    pred <- .ltmle_exact_predict_glm_fast(fit, pred_data, type = type, eps = eps)
    if (length(pred) != n || any(!is.finite(pred))) {
      .stop("ltmle_exact fast glm prediction returned non-finite or wrong-length predictions for ",
            fit$component, ".")
    }
    return(pred)
  }

  predict_one <- function(nd) {
    if (identical(fit$type, "superlearner")) {
      return(tryCatch(
        as.numeric(SuperLearner::predict.SuperLearner(fit$fit, newdata = nd)$pred),
        error = function(e) e
      ))
    }
    fam <- .nuisance_family_name(fit$family)
    ptype <- if (identical(fam, "binomial") || identical(type, "response")) "response" else "response"
    tryCatch(
      as.numeric(stats::predict(fit$fit, newdata = nd, type = ptype)),
      error = function(e) e
    )
  }

  chunk_size <- .ltmle_exact_predict_chunk_size()
  if (n > chunk_size) {
    pred <- rep(NA_real_, n)
    starts <- seq.int(1L, n, by = chunk_size)
    for (start in starts) {
      idx <- seq.int(start, min(n, start + chunk_size - 1L))
      chunk_pred <- predict_one(pred_data[idx, , drop = FALSE])
      if (inherits(chunk_pred, "error")) {
        pred <- chunk_pred
        break
      }
      if (length(chunk_pred) != length(idx)) {
        pred <- structure(
          list(message = paste0("wrong-length chunk prediction for ", fit$component)),
          class = c("simpleError", "error", "condition")
        )
        break
      }
      pred[idx] <- chunk_pred
    }
  } else {
    pred <- predict_one(pred_data)
  }

  if (inherits(pred, "error")) {
    .stop("ltmle_exact strict prediction failed for ", fit$component, ": ", pred$message)
  }
  if (length(pred) != n || any(!is.finite(pred))) {
    .stop("ltmle_exact strict prediction returned non-finite or wrong-length predictions for ", fit$component, ".")
  }

  fam <- .nuisance_family_name(fit$family)
  if (identical(fam, "binomial") || identical(type, "response")) {
    pred <- .ltmle_exact_clamp01(pred, eps = eps)
  }
  pred
}

.ltmle_exact_sigma <- function(fit, fallback = 1) {
  sig <- nuisance_sigma_safe(fit, fallback = fallback)
  if (!is.finite(sig) || sig <= 0) sig <- fallback
  sig
}

.ltmle_exact_fit_treatment_models <- function(long, T, learner, sl_library, treat_mech, p_rct,
                                              node_spec = NULL, fold_id = NULL) {
  treat_mech <- match.arg(treat_mech, c("observational", "baseline_rct"))

  if (treat_mech == "baseline_rct") {
    return(list(type = "baseline_rct", p_rct = as.numeric(p_rct), T = T))
  }

  out <- vector("list", T)
  for (tt in seq_len(T)) {
    dt <- long[long$t == tt, , drop = FALSE]
    rhs <- unique(c("W1", "W2", "Y0", as.character(node_spec$baseline_vars %||% character(0)),
                    .ltmle_exact_prior_history_covariates(node_spec)))
    rhs <- intersect(rhs, names(dt))
    out[[tt]] <- .ltmle_exact_fit_nuis(
      .ltmle_exact_formula("A", rhs),
      dt,
      stats::binomial(),
      learner,
      sl_library,
      paste0("A model t=", tt),
      fold_id = if (".fold_id" %in% names(dt)) dt$.fold_id else NULL
    )
  }
  list(type = "observational", fits = out, T = T)
}

.ltmle_exact_treatment_ratio_matrix <- function(long, treatment_models, regimen,
                                                treat_mech, p_rct,
                                                probability_bounds = c(0.01, 0.99),
                                                eps = 1e-8,
                                                max_t = NULL) {
  T <- treatment_models$T
  max_t <- if (is.null(max_t)) T else min(as.integer(max_t), T)
  if (!is.finite(max_t) || max_t < 1L) max_t <- 1L
  ids <- sort(unique(long$id))
  n <- length(ids)
  regimen <- normalize_regimen(regimen, T, "assigned_regimen")
  out <- matrix(1, nrow = n, ncol = T)

  if (treat_mech == "baseline_rct") {
    if (any(regimen != regimen[1L])) {
      .stop("baseline_rct in ltmle_exact requires a constant treatment regimen across visits.")
    }
    dt1 <- long[long$t == 1L, , drop = FALSE]
    dt1 <- dt1[match(ids, dt1$id), , drop = FALSE]
    p <- if (regimen[1L] == 1) p_rct else 1 - p_rct
    out[, 1L] <- as.numeric(dt1$A == regimen[1L]) / pmax(p, eps)
    if (T >= 2L) out[, 2:T] <- 1
    return(out)
  }

  for (tt in seq_len(max_t)) {
    dt <- long[long$t == tt, , drop = FALSE]
    dt <- dt[match(ids, dt$id), , drop = FALSE]
    p1 <- .ltmle_exact_predict_nuis(treatment_models$fits[[tt]], dt, type = "response", eps = eps)
    p1 <- .ltmle_exact_bound_probability(p1, probability_bounds)
    p_reg <- ifelse(regimen[tt] == 1, p1, 1 - p1)
    p_reg <- .ltmle_exact_bound_probability(p_reg, probability_bounds)
    out[, tt] <- as.numeric(dt$A == regimen[tt]) / pmax(p_reg, eps)
  }
  out
}

.ltmle_exact_fit_censoring_models <- function(long, T, censoring_vars, learner, sl_library,
                                             node_spec = NULL, fold_id = NULL) {
  if (is.null(censoring_vars) || !length(censoring_vars)) return(list(type = "none", T = T))

  visit_vars <- as.character(censoring_vars$visit %||% character(0))
  final_var <- censoring_vars$final %||% NA_character_
  if (length(visit_vars) != T) {
    .stop("visit censoring vars must have length T for ltmle_exact. Got ", length(visit_vars), ".")
  }

  out <- vector("list", T)
  for (tt in seq_len(T)) {
    R_col <- visit_vars[tt]
    if (!R_col %in% names(long)) {
      .stop("Missing censoring variable for t=", tt, ": ", R_col)
    }
    dt <- long[long$t == tt, , drop = FALSE]
    rhs <- .ltmle_exact_visit_censoring_covariates(node_spec)
    rhs <- intersect(rhs, names(dt))
    out[[tt]] <- .ltmle_exact_fit_nuis(
      .ltmle_exact_formula(R_col, rhs),
      dt,
      stats::binomial(),
      learner,
      sl_library,
      paste0("Censoring model t=", tt),
      fold_id = if (".fold_id" %in% names(dt)) dt$.fold_id else NULL
    )
  }

  final_fit <- NULL
  if (is.character(final_var) && length(final_var) == 1L && !is.na(final_var)) {
    if (!final_var %in% names(long)) {
      .stop("Missing final censoring variable: ", final_var)
    }
    dt <- long[long$t == T, , drop = FALSE]
    rhs <- unique(c("W1", "W2", "Y0", as.character(node_spec$baseline_vars %||% character(0)),
                    .ltmle_exact_prior_history_covariates(node_spec), "A", "M1", "M2",
                    .ltmle_exact_L_nodes(node_spec)))
    rhs <- intersect(rhs, names(dt))
    final_fit <- .ltmle_exact_fit_nuis(
      .ltmle_exact_formula(final_var, rhs),
      dt,
      stats::binomial(),
      learner,
      sl_library,
      "Final censoring model",
      fold_id = if (".fold_id" %in% names(dt)) dt$.fold_id else NULL
    )
  }

  list(type = "estimated", fits = out, visit_fits = out, final_fit = final_fit,
       vars = list(visit = visit_vars, final = final_var), T = T)
}

.ltmle_exact_uncensored_node_rows <- function(dt, tt, T, censoring_vars = NULL,
                                             require_final = FALSE) {
  censoring_vars <- .ltmle_exact_normalize_censoring_vars(censoring_vars, T)
  if (is.null(censoring_vars)) return(rep(TRUE, nrow(dt)))

  keep <- rep(TRUE, nrow(dt))

  visit_var <- censoring_vars$visit[as.integer(tt)] %||% NA_character_
  if (is.character(visit_var) && length(visit_var) == 1L && !is.na(visit_var)) {
    if (!visit_var %in% names(dt)) {
      .stop("Missing visit retention variable in long data: ", visit_var)
    }
    keep <- keep & is.finite(dt[[visit_var]]) & as.numeric(dt[[visit_var]]) == 1
  }

  final_var <- censoring_vars$final %||% NA_character_
  if (isTRUE(require_final) &&
      is.character(final_var) &&
      length(final_var) == 1L &&
      !is.na(final_var)) {
    if (!final_var %in% names(dt)) {
      .stop("Missing final retention variable in long data: ", final_var)
    }
    keep <- keep & is.finite(dt[[final_var]]) & as.numeric(dt[[final_var]]) == 1
  }

  keep[is.na(keep)] <- FALSE
  keep
}

.ltmle_exact_require_fit_rows <- function(dt, component) {
  if (!nrow(dt)) {
    .stop(
      "No uncensored rows available for fitting ", component,
      ". Check censoring variables and visit indexing."
    )
  }
  dt
}

.ltmle_exact_fit_node_models <- function(long, T, learner, sl_library, Q_model, node_spec = NULL,
                                         fold_id = NULL,
                                         censoring_vars = NULL,
                                         censoring_adjust = FALSE) {
  long <- .ltmle_exact_add_terms(long)
  censoring_vars <- if (isTRUE(censoring_adjust)) {
    .ltmle_exact_normalize_censoring_vars(censoring_vars, T)
  } else {
    NULL
  }
  L_nodes <- .ltmle_exact_L_nodes(node_spec)
  out <- list(M1 = vector("list", T), M2 = vector("list", T),
              L = stats::setNames(vector("list", length(L_nodes)), L_nodes),
              Y = vector("list", T))
  for (L_node in L_nodes) out$L[[L_node]] <- vector("list", T)
  for (tt in seq_len(T)) {
    dt_all <- long[long$t == tt, , drop = FALSE]
    keep_current <- .ltmle_exact_uncensored_node_rows(
      dt_all,
      tt = tt,
      T = T,
      censoring_vars = censoring_vars,
      require_final = FALSE
    )
    dt <- .ltmle_exact_require_fit_rows(
      dt_all[keep_current, , drop = FALSE],
      paste0("time-varying nodes at t=", tt)
    )
    if (tt >= 2L) {
      out$M1[[tt]] <- .ltmle_exact_fit_mediator_density(
        node = "M1",
        tt = tt,
        data = dt,
        history_covariates = .ltmle_exact_history_covs("M1", Q_model, node_spec),
        learner = learner,
        density_library = sl_library
      )
      out$M2[[tt]] <- .ltmle_exact_fit_mediator_density(
        node = "M2",
        tt = tt,
        data = dt,
        history_covariates = .ltmle_exact_history_covs("M2", Q_model, node_spec),
        learner = learner,
        density_library = sl_library
      )
    }
    for (L_node in L_nodes) {
      out$L[[L_node]][[tt]] <- .ltmle_exact_fit_nuis(
        .ltmle_exact_formula(L_node, .ltmle_exact_history_covs_for_L_node(L_node, Q_model, node_spec)), dt,
        stats::gaussian(), learner, sl_library, paste0(L_node, " transition t=", tt),
        fold_id = if (".fold_id" %in% names(dt)) dt$.fold_id else NULL
      )
    }
    if (identical(.ltmle_exact_outcome_type(node_spec), "longitudinal") || tt == T) {
      require_final <- identical(.ltmle_exact_outcome_type(node_spec), "terminal_only") && tt == T
      keep_y <- .ltmle_exact_uncensored_node_rows(
        dt_all,
        tt = tt,
        T = T,
        censoring_vars = censoring_vars,
        require_final = require_final
      )
      dt_y <- .ltmle_exact_require_fit_rows(
        dt_all[keep_y, , drop = FALSE],
        paste0("Y model at t=", tt)
      )
      out$Y[[tt]] <- .ltmle_exact_fit_nuis(
        .ltmle_exact_formula("Y", .ltmle_exact_history_covs("Y", Q_model, node_spec)), dt_y,
        stats::gaussian(), learner, sl_library, paste0("Y continuation t=", tt),
        fold_id = if (".fold_id" %in% names(dt_y)) dt_y$.fold_id else NULL
      )
    } else {
      out$Y[[tt]] <- NULL
    }
  }
  out$mediator_density_engine <- "gaussian_location_scale"
  out$mediator_density_models <- .ltmle_exact_mediator_density_summary(out)
  out
}

.ltmle_exact_fit_mediator_density <- function(node, tt, data, history_covariates,
                                              learner, density_library) {
  if (!node %in% c("M1", "M2")) {
    .stop("Unsupported mediator density node: ", node)
  }
  mean_fit <- .ltmle_exact_fit_nuis(
    .ltmle_exact_formula(node, history_covariates),
    data,
    stats::gaussian(),
    learner,
    density_library,
    paste0(node, " gaussian_location_scale mediator density t=", tt),
    fold_id = if (".fold_id" %in% names(data)) data$.fold_id else NULL
  )
  structure(
    list(
      engine = "gaussian_location_scale",
      node = node,
      t = as.integer(tt),
      mean_fit = mean_fit,
      sigma = .ltmle_exact_sigma(mean_fit),
      crossfit_used = isTRUE(mean_fit$crossfit_used),
      density_library = density_library
    ),
    class = c("ltmle_exact_mediator_density", "list")
  )
}

.ltmle_exact_is_mediator_density <- function(fit) {
  inherits(fit, "ltmle_exact_mediator_density")
}

.ltmle_exact_mediator_mean_fit <- function(fit) {
  if (.ltmle_exact_is_mediator_density(fit)) return(fit$mean_fit)
  fit
}

.ltmle_exact_predict_m1 <- function(fit, newdata) {
  .ltmle_exact_predict_nuis(.ltmle_exact_mediator_mean_fit(fit), newdata, "numeric")
}

.ltmle_exact_predict_m2 <- function(fit, newdata) {
  .ltmle_exact_predict_nuis(.ltmle_exact_mediator_mean_fit(fit), newdata, "numeric")
}

.ltmle_exact_dnorm_safe <- function(x, mu, sigma, eps = 1e-12) {
  pmax(stats::dnorm(as.numeric(x), mean = as.numeric(mu), sd = as.numeric(sigma)), eps)
}

.ltmle_exact_predict_mediator_sigma <- function(fit, history) {
  n <- nrow(history)
  if (!n) return(numeric(0))
  if (!.ltmle_exact_is_mediator_density(fit)) {
    return(rep(.ltmle_exact_sigma(fit), n))
  }
  mean_fit <- fit$mean_fit
  if (!is.null(mean_fit$crossfit_fits)) {
    if (!".fold_id" %in% names(history)) {
      .stop("ltmle_exact cross-fit mediator density prediction requires .fold_id for ", mean_fit$component, ".")
    }
    out <- rep(NA_real_, n)
    for (ff in names(mean_fit$crossfit_fits)) {
      idx <- as.character(history$.fold_id) == ff
      if (any(idx)) {
        out[idx] <- .ltmle_exact_sigma(mean_fit$crossfit_fits[[ff]])
      }
    }
    if (all(is.finite(out) & out > 0)) return(out)
    .stop("ltmle_exact cross-fit mediator density sigma lookup failed for ", mean_fit$component, ".")
  }
  rep(fit$sigma %||% .ltmle_exact_sigma(mean_fit), n)
}

.ltmle_exact_predict_mediator_density <- function(fit, value, history, eps = 1e-12) {
  mu <- .ltmle_exact_predict_nuis(.ltmle_exact_mediator_mean_fit(fit), history, "numeric")
  sigma <- .ltmle_exact_predict_mediator_sigma(fit, history)
  .ltmle_exact_dnorm_safe(value, mu, sigma, eps)
}

.ltmle_exact_draw_mediator <- function(fit, history) {
  mu <- .ltmle_exact_predict_nuis(.ltmle_exact_mediator_mean_fit(fit), history, "numeric")
  sigma <- .ltmle_exact_predict_mediator_sigma(fit, history)
  stats::rnorm(nrow(history), mu, sigma)
}

.ltmle_exact_mediator_density_summary <- function(models) {
  m1 <- models$M1 %||% list()
  m2 <- models$M2 %||% list()
  fits <- c(m1, m2)
  fits <- Filter(.ltmle_exact_is_mediator_density, fits)
  crossfit_used <- length(fits) > 0L && all(vapply(fits, function(fit) {
    isTRUE(fit$crossfit_used)
  }, logical(1)))
  list(
    engine = models$mediator_density_engine %||% "gaussian_location_scale",
    crossfit_used = crossfit_used
  )
}

.ltmle_exact_bound_probability <- function(x, bounds) {
  x <- as.numeric(x)
  bounds <- as.numeric(bounds)
  if (length(bounds) != 2L || any(!is.finite(bounds)) ||
      bounds[1L] <= 0 || bounds[2L] >= 1 || bounds[1L] >= bounds[2L]) {
    .stop("probability_bounds must satisfy 0 < lower < upper < 1.")
  }
  pmin(pmax(x, bounds[1L]), bounds[2L])
}

.ltmle_exact_removed_density_ratio_arg <- function() {
  paste0("density", "_ratio", "_bounds")
}

.ltmle_exact_removed_density_ratio_failure_class <- function() {
  paste0("removed_", .ltmle_exact_removed_density_ratio_arg(), "_argument_used")
}

.ltmle_exact_legacy_truncation_arg_names <- function() {
  c(
    .ltmle_exact_removed_density_ratio_arg(),
    paste0("enable_", "density", "_ratio", "_truncation"),
    paste0("density", "_ratio", "_lower_bound"),
    paste0("density", "_ratio", "_upper_bound"),
    paste0("density", "_ratio", "_bound_lower"),
    paste0("density", "_ratio", "_bound_upper")
  )
}

.ltmle_exact_stop_removed_density_ratio_arg <- function() {
  arg <- .ltmle_exact_removed_density_ratio_arg()
  cnd <- structure(
    list(
      message = paste0(
        arg,
        " has been removed from ltmle_exact. Fixed-bound density-ratio clipping is no longer supported. ",
        "Use truncation_policy='quantile' for clever_covariate_H quantile truncation, ",
        "or truncation_policy='none' for the untruncated sensitivity variant."
      ),
      call = NULL,
      failure_class = .ltmle_exact_removed_density_ratio_failure_class()
    ),
    class = c(.ltmle_exact_removed_density_ratio_failure_class(), "error", "condition")
  )
  stop(cnd)
}

.ltmle_exact_assert_no_legacy_truncation_args <- function(x) {
  nms <- names(x %||% list())
  if (is.null(nms)) return(invisible(TRUE))
  if (any(nms %in% .ltmle_exact_legacy_truncation_arg_names())) {
    .ltmle_exact_stop_removed_density_ratio_arg()
  }
  invisible(TRUE)
}

.ltmle_exact_normalize_truncation <- function(truncation_enabled = TRUE,
                                             truncation_policy = "quantile",
                                             truncation_quantile_lower = 0.01,
                                             truncation_quantile_upper = 0.99,
                                             truncation_target = "clever_covariate_H") {
  truncation_policy <- as.character(truncation_policy %||% "quantile")[1L]
  if (!truncation_policy %in% c("quantile", "none")) {
    .stop("truncation_policy must be either 'quantile' or 'none' for ltmle_exact.")
  }
  truncation_enabled <- isTRUE(truncation_enabled) && !identical(truncation_policy, "none")
  if (!isTRUE(truncation_enabled)) truncation_policy <- "none"
  truncation_target <- as.character(truncation_target %||% "clever_covariate_H")[1L]
  if (!identical(truncation_target, "clever_covariate_H")) {
    .stop("ltmle_exact truncation_target must be 'clever_covariate_H'.")
  }
  q_lower <- as.numeric(truncation_quantile_lower)[1L]
  q_upper <- as.numeric(truncation_quantile_upper)[1L]
  if (identical(truncation_policy, "quantile")) {
    if (!is.finite(q_lower) || !is.finite(q_upper) ||
        q_lower < 0 || q_upper > 1 || q_lower >= q_upper) {
      .stop("truncation quantiles must satisfy 0 <= lower < upper <= 1.")
    }
  } else {
    q_lower <- NA_real_
    q_upper <- NA_real_
  }
  list(
    enabled = truncation_enabled,
    policy = truncation_policy,
    target = truncation_target,
    quantile_lower = q_lower,
    quantile_upper = q_upper
  )
}

.ltmle_exact_quantile_value <- function(x, p) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x) || !is.finite(p)) return(NA_real_)
  as.numeric(stats::quantile(x, probs = p, na.rm = TRUE, names = FALSE, type = 7))
}

.ltmle_exact_truncation_summary_row <- function(raw,
                                                truncated,
                                                task = NULL,
                                                truncation,
                                                estimator_variant,
                                                row_evaluation_context = NA_character_,
                                                q_lower_value = NA_real_,
                                                q_upper_value = NA_real_,
                                                n_lower = 0L,
                                                n_upper = 0L) {
  task <- as.list(task %||% list())
  raw <- as.numeric(raw)
  truncated <- as.numeric(truncated)
  good_raw <- is.finite(raw)
  good_trunc <- is.finite(truncated)
  raw_finite <- raw[good_raw]
  trunc_finite <- truncated[good_trunc]
  if (!length(raw_finite)) raw_finite <- NA_real_
  if (!length(trunc_finite)) trunc_finite <- NA_real_
  n_values <- sum(good_raw)
  data.frame(
    estimator = "ltmle_exact",
    estimator_variant = estimator_variant,
    component = task$component %||% NA_character_,
    task_id = task$task_id %||% NA_character_,
    t = as.integer(task$t %||% NA_integer_),
    node = task$node %||% NA_character_,
    process_type = .ltmle_exact_process_type(task),
    row_evaluation_context = row_evaluation_context,
    truncation_enabled = isTRUE(truncation$enabled),
    truncation_policy = truncation$policy,
    truncation_target = truncation$target,
    truncation_rule = if (isTRUE(truncation$enabled) && identical(truncation$policy, "quantile")) {
      "sample_quantile"
    } else {
      "none"
    },
    requested_quantile_lower = truncation$quantile_lower,
    requested_quantile_upper = truncation$quantile_upper,
    effective_quantile_lower = if (isTRUE(truncation$enabled)) truncation$quantile_lower else NA_real_,
    effective_quantile_upper = if (isTRUE(truncation$enabled)) truncation$quantile_upper else NA_real_,
    fixed_bound_truncation_enabled = FALSE,
    fixed_bound_truncation_used = FALSE,
    density_ratio_factor_truncation_used = FALSE,
    n_values = as.integer(n_values),
    n_values_truncated_lower = as.integer(n_lower),
    n_values_truncated_upper = as.integer(n_upper),
    fraction_truncated_lower = if (n_values > 0L) n_lower / n_values else 0,
    fraction_truncated_upper = if (n_values > 0L) n_upper / n_values else 0,
    fraction_truncated_total = if (n_values > 0L) (n_lower + n_upper) / n_values else 0,
    min_raw_value = min(raw_finite, na.rm = TRUE),
    p01_raw_value = .ltmle_exact_quantile_value(raw, 0.01),
    median_raw_value = stats::median(raw_finite, na.rm = TRUE),
    p95_raw_value = .ltmle_exact_quantile_value(raw, 0.95),
    p99_raw_value = .ltmle_exact_quantile_value(raw, 0.99),
    max_raw_value = max(raw_finite, na.rm = TRUE),
    min_truncated_value = min(trunc_finite, na.rm = TRUE),
    p95_truncated_value = .ltmle_exact_quantile_value(truncated, 0.95),
    p99_truncated_value = .ltmle_exact_quantile_value(truncated, 0.99),
    max_truncated_value = max(trunc_finite, na.rm = TRUE),
    mean_raw_value = mean(raw_finite, na.rm = TRUE),
    sd_raw_value = stats::sd(raw_finite, na.rm = TRUE),
    mean_truncated_value = mean(trunc_finite, na.rm = TRUE),
    sd_truncated_value = stats::sd(trunc_finite, na.rm = TRUE),
    ess_raw = .ltmle_exact_ess(abs(raw)),
    ess_truncated = .ltmle_exact_ess(abs(truncated)),
    H_quantile_lower_value = q_lower_value,
    H_quantile_upper_value = q_upper_value,
    H_truncation_applied = (n_lower + n_upper) > 0L,
    n_H_truncated_lower = as.integer(n_lower),
    n_H_truncated_upper = as.integer(n_upper),
    fraction_H_truncated_total = if (n_values > 0L) (n_lower + n_upper) / n_values else 0,
    passed = TRUE,
    failure_class = "no_failure",
    stringsAsFactors = FALSE
  )
}

.ltmle_exact_decomposition_summary_row <- function(H_raw,
                                                   H_truncated,
                                                   task = NULL,
                                                   truncation,
                                                   estimator_variant,
                                                   row_evaluation_context = NA_character_) {
  task <- as.list(task %||% list())
  attrs <- attributes(H_raw)
  A_part <- as.numeric(attrs$A_part %||% NA_real_)
  M_part <- as.numeric(attrs$M_part %||% NA_real_)
  C_part <- as.numeric(attrs$C_part %||% NA_real_)
  finite_max <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    max(abs(x), na.rm = TRUE)
  }
  finite_mean <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    mean(x, na.rm = TRUE)
  }
  data.frame(
    estimator = "ltmle_exact",
    estimator_variant = estimator_variant,
    component = task$component %||% NA_character_,
    task_id = task$task_id %||% NA_character_,
    t = as.integer(task$t %||% NA_integer_),
    node = task$node %||% NA_character_,
    process_type = .ltmle_exact_process_type(task),
    row_evaluation_context = row_evaluation_context,
    A_part = finite_mean(A_part),
    M_part = finite_mean(M_part),
    C_part = finite_mean(C_part),
    H_raw = finite_mean(H_raw),
    H_truncated = finite_mean(H_truncated),
    A_part_max = finite_max(A_part),
    M_part_max = finite_max(M_part),
    C_part_max = finite_max(C_part),
    H_raw_max = finite_max(H_raw),
    H_truncated_max = finite_max(H_truncated),
    clever_covariate_truncation_applied =
      any(abs(as.numeric(H_raw) - as.numeric(H_truncated)) > 0, na.rm = TRUE),
    clever_covariate_truncation_policy = truncation$policy,
    clever_covariate_quantile_lower = truncation$quantile_lower,
    clever_covariate_quantile_upper = truncation$quantile_upper,
    density_ratio_factor_truncation_used = FALSE,
    stringsAsFactors = FALSE
  )
}

.ltmle_exact_truncate_clever_covariate <- function(H,
                                                   truncation,
                                                   task = NULL,
                                                   estimator_variant = "ltmle_exact_quantile_truncated",
                                                   row_evaluation_context = NA_character_) {
  attrs <- attributes(H)
  raw <- as.numeric(H)
  out <- raw
  q_lower_value <- NA_real_
  q_upper_value <- NA_real_
  n_lower <- 0L
  n_upper <- 0L
  if (isTRUE(truncation$enabled) && identical(truncation$policy, "quantile")) {
    q_lower_value <- .ltmle_exact_quantile_value(raw, truncation$quantile_lower)
    q_upper_value <- .ltmle_exact_quantile_value(raw, truncation$quantile_upper)
    if (is.finite(q_lower_value) && is.finite(q_upper_value) && q_lower_value <= q_upper_value) {
      n_lower <- sum(is.finite(raw) & raw < q_lower_value)
      n_upper <- sum(is.finite(raw) & raw > q_upper_value)
      out <- pmin(pmax(raw, q_lower_value), q_upper_value)
    }
  }
  attrs$n_clever_covariate_values_truncated <- as.integer(n_lower + n_upper)
  attrs$fraction_clever_covariate_values_truncated <- if (sum(is.finite(raw)) > 0L) {
    (n_lower + n_upper) / sum(is.finite(raw))
  } else {
    0
  }
  attrs$H_raw <- raw
  attrs$H_truncated <- out
  attrs$H_quantile_lower_value <- q_lower_value
  attrs$H_quantile_upper_value <- q_upper_value
  attrs$n_H_truncated_lower <- as.integer(n_lower)
  attrs$n_H_truncated_upper <- as.integer(n_upper)
  attrs$truncation_diagnostics <- .ltmle_exact_truncation_summary_row(
    raw = raw,
    truncated = out,
    task = task,
    truncation = truncation,
    estimator_variant = estimator_variant,
    row_evaluation_context = row_evaluation_context,
    q_lower_value = q_lower_value,
    q_upper_value = q_upper_value,
    n_lower = n_lower,
    n_upper = n_upper
  )
  attrs$clever_covariate_decomposition_diagnostics <- .ltmle_exact_decomposition_summary_row(
    H_raw = H,
    H_truncated = out,
    task = task,
    truncation = truncation,
    estimator_variant = estimator_variant,
    row_evaluation_context = row_evaluation_context
  )
  attributes(out) <- attrs
  out
}

.ltmle_exact_normal_quadrature <- function(n) {
  if (requireNamespace("statmod", quietly = TRUE)) {
    gh <- statmod::gauss.quad.prob(n, dist = "normal")
    return(list(nodes = gh$nodes, weights = gh$weights))
  }
  p <- (seq_len(n) - 0.5) / n
  list(nodes = stats::qnorm(p), weights = rep(1 / n, n))
}

.ltmle_exact_second_mediator_law_density_ratio <- function(observed_long,
                                                          models,
                                                          T,
                                                          second_regimen,
                                                          n_particles,
                                                          seed,
                                                          eps = 1e-12,
                                                          node_spec = NULL,
                                                          max_t = NULL) {
  if (!is.null(seed)) set.seed(seed)
  ids <- sort(unique(observed_long$id))
  n <- length(ids)
  second_regimen <- normalize_regimen(second_regimen, T, "second_regimen")
  n_particles <- as.integer(n_particles)
  if (!is.finite(n_particles) || n_particles < 2L) n_particles <- 2000L

  out <- matrix(1, nrow = n, ncol = T)
  max_t <- if (is.null(max_t)) T else min(as.integer(max_t), T)
  if (!is.finite(max_t) || max_t < 1L) max_t <- 1L
  if (max_t < 2L) {
    return(out)
  }

  base <- observed_long[observed_long$t == 1L, , drop = FALSE]
  base <- base[match(ids, base$id), , drop = FALSE]
  history <- base[rep(seq_len(nrow(base)), each = n_particles), , drop = FALSE]
  history$id0 <- rep(seq_len(n), each = n_particles)
  history$id <- seq_len(nrow(history))
  history$.w <- rep(1 / n_particles, nrow(history))

  for (tt in seq_len(max_t)) {
    obs_tt <- observed_long[observed_long$t == tt, , drop = FALSE]
    obs_tt <- obs_tt[match(ids, obs_tt$id), , drop = FALSE]

    row <- data.frame(
      id0 = history$id0,
      id = seq_len(nrow(history)),
      W1 = history$W1,
      W2 = history$W2,
      Y0 = history$Y0,
      M1_0 = history$M1_0,
      M2_0 = history$M2_0,
      M1_lag = history$M1_lag,
      M2_lag = history$M2_lag,
      A = second_regimen[tt],
      stringsAsFactors = FALSE
    )
    if ("Y_lag" %in% names(history)) row$Y_lag <- history$Y_lag
    for (L_node in .ltmle_exact_L_nodes(node_spec)) {
      lag_nm <- .ltmle_exact_lag_name(L_node)
      if (lag_nm %in% names(history)) row[[lag_nm]] <- history[[lag_nm]]
    }
    for (nm in as.character(node_spec$baseline_vars %||% character(0))) {
      if (nm %in% names(history) && !nm %in% names(row)) row[[nm]] <- history[[nm]]
    }
    if (".fold_id" %in% names(history)) row$.fold_id <- history$.fold_id

    if (tt == 1L) {
      row$M1 <- row$M1_0
      row$M2 <- row$M2_0
    } else {
      m1_hist <- .ltmle_exact_set_A(row, second_regimen[tt])
      row$M1 <- .ltmle_exact_draw_mediator(models$M1[[tt]], m1_hist)

      m2_hist <- .ltmle_exact_set_A(row, second_regimen[tt])
      m2_hist$M1 <- row$M1
      m2_hist <- .ltmle_exact_add_terms(m2_hist)

      target_m2 <- obs_tt$M2[row$id0]
      inc_num <- .ltmle_exact_predict_mediator_density(models$M2[[tt]], target_m2, m2_hist, eps)
      num_weighted <- rowsum(history$.w * inc_num, group = row$id0, reorder = TRUE)[, 1L]
      weight_sum <- rowsum(history$.w, group = row$id0, reorder = TRUE)[, 1L]
      incremental_numerator <- as.numeric(num_weighted / pmax(weight_sum, eps))

      den_hist <- .ltmle_exact_set_A(obs_tt, obs_tt$A)
      denominator <- .ltmle_exact_predict_mediator_density(models$M2[[tt]], obs_tt$M2, den_hist, eps)
      out[, tt] <- incremental_numerator / denominator

      history$.w <- history$.w * inc_num
      weight_sum_particle <- rowsum(history$.w, group = row$id0, reorder = TRUE)
      weight_sum_particle <- weight_sum_particle[
        match(row$id0, as.integer(rownames(weight_sum_particle))), 1L
      ]
      history$.w <- history$.w / pmax(weight_sum_particle, eps)

      row$M2 <- target_m2
    }

    row <- .ltmle_exact_add_terms(row)
    row <- .ltmle_exact_draw_L_nodes(row, models, tt, node_spec)
    row <- .ltmle_exact_add_terms(row)
    draw_y <- identical(.ltmle_exact_outcome_type(node_spec), "longitudinal") || tt == T
    if (isTRUE(draw_y) && !is.null(models$Y[[tt]])) {
      row$Y <- .ltmle_exact_draw_y_history(row, models$Y[[tt]])
    }

    history$M1_lag <- row$M1
    history$M2_lag <- row$M2
    for (L_node in .ltmle_exact_L_nodes(node_spec)) {
      history[[.ltmle_exact_lag_name(L_node)]] <- row[[L_node]]
    }
    if (!identical(.ltmle_exact_outcome_type(node_spec), "terminal_only") && "Y" %in% names(row)) {
      history$Y_lag <- row$Y
    }
  }

  if (any(!is.finite(out))) .stop("Non-finite second-mediator trajectory density ratio.")
  out
}

.ltmle_exact_generated_second_mediator_conditional_density_ratio <- function(generated_row,
                                                                             models,
                                                                             tt,
                                                                             second_regimen_value,
                                                                             eps = 1e-12) {
  if (tt < 2L) return(rep(1, nrow(generated_row)))
  num_hist <- .ltmle_exact_set_A(generated_row, second_regimen_value)
  num_hist$M1 <- generated_row$M1
  num_hist <- .ltmle_exact_add_terms(num_hist)
  den_hist <- .ltmle_exact_set_A(generated_row, generated_row$A)
  ratio <- .ltmle_exact_predict_mediator_density(models$M2[[tt]], generated_row$M2, num_hist, eps) /
    .ltmle_exact_predict_mediator_density(models$M2[[tt]], generated_row$M2, den_hist, eps)
  if (any(!is.finite(ratio))) {
    .stop("Non-finite generated second-mediator density ratio at t=", tt)
  }
  ratio
}

.ltmle_exact_generated_second_mediator_marginal_density_ratio_for_outcome_process <- function(generated_long,
                                                                                              models,
                                                                                              T,
                                                                                              second_regimen,
                                                                                              n_particles,
                                                                                              seed,
                                                                                              eps = 1e-12,
                                                                                              node_spec = NULL,
                                                                                              max_t = NULL) {
  .ltmle_exact_second_mediator_law_density_ratio(
    observed_long = generated_long,
    models = models,
    T = T,
    second_regimen = second_regimen,
    n_particles = n_particles,
    seed = seed,
    eps = eps,
    node_spec = node_spec,
    max_t = max_t
  )
}

.ltmle_exact_independent_second_M2_marginal_reference <- function(long,
                                                                  models,
                                                                  T,
                                                                  second_regimen,
                                                                  n_particles,
                                                                  seed,
                                                                  eps = 1e-12,
                                                                  node_spec = NULL,
                                                                  max_t = NULL) {
  if (!is.null(seed)) set.seed(seed)
  ids <- sort(unique(long$id))
  n <- length(ids)
  second_regimen <- normalize_regimen(second_regimen, T, "independent_second_M2_regimen")
  n_particles <- as.integer(n_particles)
  if (!is.finite(n_particles) || n_particles < 2L) n_particles <- 2000L

  out <- matrix(1, nrow = n, ncol = T)
  max_t <- if (is.null(max_t)) T else min(as.integer(max_t), T)
  if (!is.finite(max_t) || max_t < 1L) max_t <- 1L
  if (max_t < 2L) {
    return(out)
  }

  base <- long[long$t == 1L, , drop = FALSE]
  base <- base[match(ids, base$id), , drop = FALSE]
  particle_history <- base[rep(seq_len(nrow(base)), each = n_particles), , drop = FALSE]
  particle_history$id0 <- rep(seq_len(n), each = n_particles)
  particle_history$id <- seq_len(nrow(particle_history))
  particle_history$.reference_weight <- rep(1 / n_particles, nrow(particle_history))

  for (tt in seq_len(max_t)) {
    observed_t <- long[long$t == tt, , drop = FALSE]
    observed_t <- observed_t[match(ids, observed_t$id), , drop = FALSE]

    candidate <- data.frame(
      id0 = particle_history$id0,
      id = seq_len(nrow(particle_history)),
      W1 = particle_history$W1,
      W2 = particle_history$W2,
      Y0 = particle_history$Y0,
      M1_0 = particle_history$M1_0,
      M2_0 = particle_history$M2_0,
      M1_lag = particle_history$M1_lag,
      M2_lag = particle_history$M2_lag,
      A = second_regimen[tt],
      stringsAsFactors = FALSE
    )
    if ("Y_lag" %in% names(particle_history)) candidate$Y_lag <- particle_history$Y_lag
    for (L_node in .ltmle_exact_L_nodes(node_spec)) {
      lag_nm <- .ltmle_exact_lag_name(L_node)
      if (lag_nm %in% names(particle_history)) candidate[[lag_nm]] <- particle_history[[lag_nm]]
    }
    for (nm in as.character(node_spec$baseline_vars %||% character(0))) {
      if (nm %in% names(particle_history) && !nm %in% names(candidate)) {
        candidate[[nm]] <- particle_history[[nm]]
      }
    }
    if (".fold_id" %in% names(particle_history)) {
      candidate$.fold_id <- particle_history$.fold_id
    }

    if (tt == 1L) {
      candidate$M1 <- candidate$M1_0
      candidate$M2 <- candidate$M2_0
    } else {
      m1_reference_history <- .ltmle_exact_set_A(candidate, second_regimen[tt])
      candidate$M1 <- .ltmle_exact_draw_mediator(models$M1[[tt]], m1_reference_history)

      m2_reference_history <- .ltmle_exact_set_A(candidate, second_regimen[tt])
      m2_reference_history$M1 <- candidate$M1
      m2_reference_history <- .ltmle_exact_add_terms(m2_reference_history)

      target_m2 <- observed_t$M2[candidate$id0]
      numerator_density <- .ltmle_exact_predict_mediator_density(
        models$M2[[tt]],
        target_m2,
        m2_reference_history,
        eps
      )
      numerator_weighted <- rowsum(
        particle_history$.reference_weight * numerator_density,
        group = candidate$id0,
        reorder = TRUE
      )[, 1L]
      denominator_weight_sum <- rowsum(
        particle_history$.reference_weight,
        group = candidate$id0,
        reorder = TRUE
      )[, 1L]
      numerator <- as.numeric(numerator_weighted / pmax(denominator_weight_sum, eps))

      observed_denominator_history <- .ltmle_exact_set_A(observed_t, observed_t$A)
      denominator <- .ltmle_exact_predict_mediator_density(
        models$M2[[tt]],
        observed_t$M2,
        observed_denominator_history,
        eps
      )
      out[, tt] <- numerator / denominator

      particle_history$.reference_weight <- particle_history$.reference_weight * numerator_density
      particle_weight_sum <- rowsum(
        particle_history$.reference_weight,
        group = candidate$id0,
        reorder = TRUE
      )
      particle_weight_sum <- particle_weight_sum[
        match(candidate$id0, as.integer(rownames(particle_weight_sum))), 1L
      ]
      particle_history$.reference_weight <- particle_history$.reference_weight /
        pmax(particle_weight_sum, eps)

      candidate$M2 <- target_m2
    }

    candidate <- .ltmle_exact_add_terms(candidate)
    candidate <- .ltmle_exact_draw_L_nodes(candidate, models, tt, node_spec)
    candidate <- .ltmle_exact_add_terms(candidate)
    draw_y <- identical(.ltmle_exact_outcome_type(node_spec), "longitudinal") || tt == T
    if (isTRUE(draw_y) && !is.null(models$Y[[tt]])) {
      candidate$Y <- .ltmle_exact_draw_y_history(candidate, models$Y[[tt]])
    }

    particle_history$M1_lag <- candidate$M1
    particle_history$M2_lag <- candidate$M2
    for (L_node in .ltmle_exact_L_nodes(node_spec)) {
      particle_history[[.ltmle_exact_lag_name(L_node)]] <- candidate[[L_node]]
    }
    if (!identical(.ltmle_exact_outcome_type(node_spec), "terminal_only") &&
        "Y" %in% names(candidate)) {
      particle_history$Y_lag <- candidate$Y
    }
  }

  if (any(!is.finite(out))) {
    .stop("Non-finite independent second-M2 marginal reference.")
  }
  out
}

.ltmle_exact_second_M2_marginal_reference_check <- function(component,
                                                           component_tasks,
                                                           long,
                                                           models,
                                                           T,
                                                           spec,
                                                           production_mediator_cache,
                                                           n_particles,
                                                           seed = 202405L,
                                                           eps = 1e-12,
                                                           node_spec = NULL) {
  if (!identical(as.character(spec$world_type[1L]), "separate")) {
    return(data.frame())
  }
  if (is.null(production_mediator_cache$second_M2_marginal_for_outcome_process)) {
    .stop("Missing production second_M2_marginal_for_outcome_process for component: ", component)
  }
  reference <- .ltmle_exact_independent_second_M2_marginal_reference(
    long = long,
    models = models,
    T = T,
    second_regimen = as.numeric(spec$m2_A),
    n_particles = n_particles,
    seed = seed,
    eps = eps,
    node_spec = node_spec,
    max_t = T
  )
  production <- production_mediator_cache$second_M2_marginal_for_outcome_process
  ids <- sort(unique(long$id))
  expected <- .ltmle_exact_expected_separate_labels(component)
  outer_regimen <- .ltmle_exact_regimen_label(spec$outer_A)
  m1_regimen <- .ltmle_exact_regimen_label(spec$m1_A)
  m2_regimen <- .ltmle_exact_regimen_label(spec$m2_A)
  factor_ok <- identical(unname(c(outer_regimen, m1_regimen, m2_regimen)), unname(expected))

  relevant_process <- c(
    "outcome_process",
    "post_mediator_covariate_transition",
    "virtual_mixed_continuation_task"
  )
  tasks <- component_tasks[
    as.character(component_tasks$component) == as.character(component) &
      as.character(component_tasks$process_type) %in% relevant_process &
      as.integer(component_tasks$t) >= 2L,
    ,
    drop = FALSE
  ]
  if (!nrow(tasks)) {
    tasks <- data.frame(
      component = component,
      task_id = paste0(component, "::second_M2_marginal_for_outcome_process"),
      t = seq.int(2L, T),
      node = "outcome",
      process_type = "outcome_process",
      stringsAsFactors = FALSE
    )
  }

  out <- vector("list", nrow(tasks))
  for (ii in seq_len(nrow(tasks))) {
    task <- tasks[ii, , drop = FALSE]
    tt <- as.integer(task$t[1L])
    observed_t <- long[long$t == tt, , drop = FALSE]
    observed_t <- observed_t[match(ids, observed_t$id), , drop = FALSE]
    prod_val <- as.numeric(production[, tt])
    ref_val <- as.numeric(reference[, tt])
    tol <- pmax(1e-8, 1e-7 * pmax(1, abs(prod_val), abs(ref_val)))
    abs_diff <- abs(prod_val - ref_val)
    rel_diff <- abs_diff / pmax(1, abs(ref_val))
    identity_match <- is.finite(max(abs_diff, na.rm = TRUE)) && all(abs_diff <= tol, na.rm = TRUE)
    conditioning_ok <- identical(as.character(task$process_type[1L]), "virtual_mixed_continuation_task") ||
      as.character(task$process_type[1L]) %in% c("outcome_process", "post_mediator_covariate_transition")
    passed <- isTRUE(identity_match) && isTRUE(factor_ok) && isTRUE(conditioning_ok)
    out[[ii]] <- data.frame(
      component = component,
      task_id = as.character(task$task_id[1L]),
      t = tt,
      node = as.character(task$node[1L]),
      process_type = as.character(task$process_type[1L]),
      outer_regimen = outer_regimen,
      m1_regimen = m1_regimen,
      m2_regimen = m2_regimen,
      source_boundary_type = "observed_targeting_rows",
      outcome_history_state = "outcome",
      m1_history_state = "first_law",
      m2_history_state = "second_law",
      M2_value = as.numeric(observed_t$M2),
      conditioning_history_label = "second_M2_marginal_for_outcome_process",
      conditioning_history_source = "independent_second_law_reference_mc",
      production_second_M2_marginal = prod_val,
      reference_second_M2_marginal = ref_val,
      abs_difference = abs_diff,
      relative_difference = rel_diff,
      tolerance = tol,
      reference_source = "independent_bruteforce_or_reference_mc",
      reference_uses_ltmle_second_mediator_law_density_ratio = FALSE,
      reference_uses_ltmle_mediator_ratio_helper = FALSE,
      reference_uses_ratio_cache = FALSE,
      reference_uses_targeting_H_attributes = FALSE,
      reference_uses_production_H_helper = FALSE,
      independent_reference_used = TRUE,
      marginal_identity_matches = abs_diff <= tol,
      factor_regimen_assignment_correct = factor_ok,
      conditioning_history_correct = conditioning_ok,
      passed = passed & (abs_diff <= tol),
      failure_class = if (passed) {
        "no_failure"
      } else if (!identity_match) {
        "second_M2_marginal_identity_mismatch"
      } else if (!factor_ok) {
        "second_M2_marginal_wrong_regimen_assignment"
      } else if (!conditioning_ok) {
        "second_M2_marginal_wrong_conditioning_history"
      } else {
        "second_M2_reference_check_failed"
      },
      stringsAsFactors = FALSE
    )
  }
  .ltmle_exact_rbind_fill(out)
}

.ltmle_exact_build_mediator_density_ratio_cache <- function(long, models, T, spec,
                                                           eps = 1e-12,
                                                           history_source = c("observed", "generated"),
                                                           node_spec = NULL,
                                                           n_particles = 2000L,
                                                           seed = 202405L,
                                                           max_t = NULL) {
  history_source <- match.arg(history_source)
  max_t <- if (is.null(max_t)) T else min(as.integer(max_t), T)
  if (!is.finite(max_t) || max_t < 1L) max_t <- 1L
  ids <- sort(unique(long$id))
  n <- length(ids)

  out <- list(
    joint_M1 = matrix(1, n, T),
    joint_M2 = matrix(1, n, T),
    first_M1 = matrix(1, n, T),
    first_M2_conditional_for_first_mediator_law = matrix(1, n, T),
    second_M1 = matrix(1, n, T),
    second_M2_marginal_for_outcome_process = matrix(1, n, T),
    second_M2_conditional_for_second_mediator_law = matrix(1, n, T),
    natural_M1 = matrix(1, n, T),
    natural_M2 = matrix(1, n, T)
  )
  second_M2_marginal <- matrix(1, nrow = n, ncol = T)
  if (identical(spec$world_type[1L], "separate") && max_t >= 2L) {
    second_M2_marginal <- if (identical(history_source, "observed")) {
      .ltmle_exact_second_mediator_law_density_ratio(
        observed_long = long,
        models = models,
        T = T,
        second_regimen = as.numeric(spec$m2_A),
        n_particles = n_particles,
        seed = seed,
        eps = eps,
        node_spec = node_spec,
        max_t = max_t
      )
    } else {
      .ltmle_exact_generated_second_mediator_marginal_density_ratio_for_outcome_process(
        generated_long = long,
        models = models,
        T = T,
        second_regimen = as.numeric(spec$m2_A),
        n_particles = n_particles,
        seed = seed,
        eps = eps,
        node_spec = node_spec,
        max_t = max_t
      )
    }
  }

  for (tt in seq_len(max_t)) {
    if (tt < 2L) next
    dt <- long[long$t == tt, , drop = FALSE]
    dt <- dt[match(ids, dt$id), , drop = FALSE]
    dt <- .ltmle_exact_add_terms(dt)

    sp <- spec[spec$t == tt, , drop = FALSE]

    den_M1_hist <- .ltmle_exact_set_A(dt, dt$A)
    den_M2_hist <- .ltmle_exact_set_A(dt, dt$A)

    den_M1_density <- .ltmle_exact_predict_mediator_density(models$M1[[tt]], dt$M1, den_M1_hist, eps)
    den_M2_density <- .ltmle_exact_predict_mediator_density(models$M2[[tt]], dt$M2, den_M2_hist, eps)

    if (identical(sp$world_type[1L], "joint")) {
      joint_M1_hist <- .ltmle_exact_set_A(dt, sp$m1_A[1L])
      joint_M2_hist <- .ltmle_exact_set_A(dt, sp$m2_A[1L])
      joint_M2_hist$M1 <- dt$M1
      joint_M2_hist <- .ltmle_exact_add_terms(joint_M2_hist)

      out$joint_M1[, tt] <- .ltmle_exact_predict_mediator_density(models$M1[[tt]], dt$M1, joint_M1_hist, eps) /
        den_M1_density
      out$joint_M2[, tt] <- .ltmle_exact_predict_mediator_density(models$M2[[tt]], dt$M2, joint_M2_hist, eps) /
        den_M2_density
    }

    if (identical(sp$world_type[1L], "separate")) {
      first_M1_hist <- .ltmle_exact_set_A(dt, sp$m1_A[1L])
      out$first_M1[, tt] <- .ltmle_exact_predict_mediator_density(models$M1[[tt]], dt$M1, first_M1_hist, eps) /
        den_M1_density

      first_M2_conditional_hist <- .ltmle_exact_set_A(dt, sp$m1_A[1L])
      first_M2_conditional_hist$M1 <- dt$M1
      first_M2_conditional_hist <- .ltmle_exact_add_terms(first_M2_conditional_hist)
      if (identical(history_source, "generated")) {
        out$first_M2_conditional_for_first_mediator_law[, tt] <-
          .ltmle_exact_generated_second_mediator_conditional_density_ratio(
            generated_row = dt,
            models = models,
            tt = tt,
            second_regimen_value = sp$m1_A[1L],
            eps = eps
          )
      } else {
        out$first_M2_conditional_for_first_mediator_law[, tt] <-
          .ltmle_exact_predict_mediator_density(models$M2[[tt]], dt$M2, first_M2_conditional_hist, eps) /
          den_M2_density
      }

      second_M1_hist <- .ltmle_exact_set_A(dt, sp$m2_A[1L])
      out$second_M1[, tt] <- .ltmle_exact_predict_mediator_density(models$M1[[tt]], dt$M1, second_M1_hist, eps) /
        den_M1_density

      second_M2_conditional_hist <- .ltmle_exact_set_A(dt, sp$m2_A[1L])
      second_M2_conditional_hist$M1 <- dt$M1
      second_M2_conditional_hist <- .ltmle_exact_add_terms(second_M2_conditional_hist)
      if (identical(history_source, "generated")) {
        out$second_M2_conditional_for_second_mediator_law[, tt] <-
          .ltmle_exact_generated_second_mediator_conditional_density_ratio(
            generated_row = dt,
            models = models,
            tt = tt,
            second_regimen_value = sp$m2_A[1L],
            eps = eps
          )
      } else {
        out$second_M2_conditional_for_second_mediator_law[, tt] <-
          .ltmle_exact_predict_mediator_density(models$M2[[tt]], dt$M2, second_M2_conditional_hist, eps) /
          den_M2_density
      }
      out$second_M2_marginal_for_outcome_process[, tt] <- second_M2_marginal[, tt]
    }
  }

  for (nm in names(out)) {
    if (any(!is.finite(out[[nm]]))) .stop("Non-finite mediator density ratio in ", nm, ".")
  }
  out
}

.ltmle_exact_censoring_ratio_matrix <- function(long, censoring_models,
                                                probability_bounds = c(0.01, 0.99),
                                                eps = 1e-8,
                                                max_t = NULL) {
  T <- censoring_models$T
  max_t <- if (is.null(max_t)) T else min(as.integer(max_t), T)
  if (!is.finite(max_t) || max_t < 1L) max_t <- 1L
  ids <- sort(unique(long$id))
  n <- length(ids)
  out <- matrix(1, nrow = n, ncol = T)
  final <- rep(1, n)
  if (is.null(censoring_models) || identical(censoring_models$type, "none")) {
    return(list(visit = out, final = final, final_var = NA_character_))
  }

  for (tt in seq_len(max_t)) {
    dt <- long[long$t == tt, , drop = FALSE]
    dt <- dt[match(ids, dt$id), , drop = FALSE]
    R_col <- censoring_models$vars$visit[tt]
    if (!R_col %in% names(dt)) dt[[R_col]] <- 1
    dt[[R_col]][!is.finite(as.numeric(dt[[R_col]]))] <- 1
    p_uncensored <- .ltmle_exact_predict_nuis(censoring_models$visit_fits[[tt]] %||% censoring_models$fits[[tt]],
                                              dt, type = "response", eps = eps)
    p_uncensored <- .ltmle_exact_bound_probability(p_uncensored, probability_bounds)
    out[, tt] <- as.numeric(dt[[R_col]] == 1) / pmax(p_uncensored, eps)
  }
  final_var <- censoring_models$vars$final %||% NA_character_
  if (max_t >= T &&
      !is.null(censoring_models$final_fit) &&
      is.character(final_var) && length(final_var) == 1L && !is.na(final_var)) {
    dt <- long[long$t == T, , drop = FALSE]
    dt <- dt[match(ids, dt$id), , drop = FALSE]
    if (!final_var %in% names(dt)) dt[[final_var]] <- 1
    dt[[final_var]][!is.finite(as.numeric(dt[[final_var]]))] <- 1
    p_final <- .ltmle_exact_predict_nuis(censoring_models$final_fit, dt, type = "response", eps = eps)
    p_final <- .ltmle_exact_bound_probability(p_final, probability_bounds)
    final <- as.numeric(dt[[final_var]] == 1) / pmax(p_final, eps)
  }
  list(visit = out, final = final, final_var = final_var)
}

.ltmle_exact_build_ratio_cache <- function(long, models, treatment_models, T, spec,
                                           probability_bounds,
                                           treat_mech, p_rct, censoring_models = NULL,
                                           history_source = c("observed", "generated"),
                                           node_spec = NULL,
                                           ltmle_exact_density_ratio_mc_n = 2000L,
                                           density_ratio_seed = 202405L,
                                           max_t = NULL) {
  history_source <- match.arg(history_source)
  max_t <- if (is.null(max_t)) T else min(as.integer(max_t), T)
  if (!is.finite(max_t) || max_t < 1L) max_t <- 1L
  outer_regimen <- as.numeric(spec$outer_A)
  first_regimen <- as.numeric(spec$m1_A)
  second_regimen <- as.numeric(spec$m2_A)

  A_outer <- .ltmle_exact_treatment_ratio_matrix(long, treatment_models, outer_regimen, treat_mech, p_rct,
                                                 probability_bounds = probability_bounds,
                                                 max_t = max_t)
  A_first <- .ltmle_exact_treatment_ratio_matrix(long, treatment_models, first_regimen, treat_mech, p_rct,
                                                 probability_bounds = probability_bounds,
                                                 max_t = max_t)
  A_second <- .ltmle_exact_treatment_ratio_matrix(long, treatment_models, second_regimen, treat_mech, p_rct,
                                                  probability_bounds = probability_bounds,
                                                  max_t = max_t)

  if (identical(spec$world_type[1L], "joint")) {
    A_joint <- .ltmle_exact_treatment_ratio_matrix(long, treatment_models, first_regimen, treat_mech, p_rct,
                                                   probability_bounds = probability_bounds,
                                                   max_t = max_t)
  } else {
    A_joint <- matrix(1, nrow = nrow(A_outer), ncol = T)
  }

  mediator_ratios <- .ltmle_exact_build_mediator_density_ratio_cache(
    long = long,
    models = models,
    T = T,
    spec = spec,
    history_source = history_source,
    node_spec = node_spec,
    n_particles = ltmle_exact_density_ratio_mc_n,
    seed = density_ratio_seed,
    max_t = max_t
  )

  C <- .ltmle_exact_censoring_ratio_matrix(long, censoring_models %||% list(type = "none", T = T),
                                           probability_bounds = probability_bounds,
                                           max_t = max_t)

  out <- list(
    A_outer = A_outer,
    A_first_mediator_law = A_first,
    A_second_mediator_law = A_second,
    A_joint_mediator_law = A_joint,
    M = mediator_ratios,
    C = C$visit,
    C_final = C$final,
    C_final_var = C$final_var,
    T = T
  )
  mats <- c(out[grepl("^A_", names(out))], list(C = out$C, C_final = matrix(out$C_final, ncol = 1L)), out$M)
  for (nm in names(mats)) {
    if (any(!is.finite(mats[[nm]]))) .stop("Non-finite likelihood ratio in ", nm, ".")
  }
  out
}

.ltmle_exact_cumprod_to <- function(mat, tt) {
  if (tt <= 0L) return(rep(1, nrow(mat)))
  apply(mat[, seq_len(tt), drop = FALSE], 1L, prod)
}

.ltmle_exact_mediator_ratio_parts_for_task <- function(task, mediator_cache) {
  node <- task$node
  if (startsWith(node, "L_")) node <- "L"
  tt <- task$t
  process_type <- .ltmle_exact_process_type(task)
  world_type <- task$world_type %||% NA_character_
  n <- nrow(mediator_cache$joint_M1)

  prior <- function(mat) {
    if (tt <= 1L) return(rep(1, n))
    .ltmle_exact_cumprod_to(mat, tt - 1L)
  }
  through <- function(mat) {
    if (tt <= 1L) return(rep(1, n))
    .ltmle_exact_cumprod_to(mat, tt)
  }

  if (process_type %in% c(
    "outcome_process",
    "post_mediator_covariate_transition",
    "virtual_mixed_continuation_task"
  )) {
    if (identical(world_type, "joint")) {
      M1_part <- through(mediator_cache$joint_M1)
      M2_part <- through(mediator_cache$joint_M2)
      return(list(M1_part = M1_part, M2_part = M2_part, M_part = M1_part * M2_part))
    }
    if (identical(world_type, "separate")) {
      M1_part <- through(mediator_cache$first_M1)
      M2_part <- through(mediator_cache$second_M2_marginal_for_outcome_process)
      return(list(M1_part = M1_part, M2_part = M2_part, M_part = M1_part * M2_part))
    }
    return(list(M1_part = rep(1, n), M2_part = rep(1, n), M_part = rep(1, n)))
  }
  if (identical(process_type, "observed_mediator_process")) {
    return(list(M1_part = rep(1, n), M2_part = rep(1, n), M_part = rep(1, n)))
  }
  if (identical(process_type, "joint_stochastic_mediator_intervention_law")) {
    if (identical(node, "M1")) {
      M1_part <- prior(mediator_cache$joint_M1)
      M2_part <- prior(mediator_cache$joint_M2)
      return(list(M1_part = M1_part, M2_part = M2_part, M_part = M1_part * M2_part))
    }
    if (identical(node, "M2")) {
      M1_part <- through(mediator_cache$joint_M1)
      M2_part <- prior(mediator_cache$joint_M2)
      return(list(M1_part = M1_part, M2_part = M2_part, M_part = M1_part * M2_part))
    }
    if (identical(node, "L") || identical(node, "Y")) {
      M1_part <- through(mediator_cache$joint_M1)
      M2_part <- through(mediator_cache$joint_M2)
      return(list(M1_part = M1_part, M2_part = M2_part, M_part = M1_part * M2_part))
    }
  }
  if (identical(process_type, "first_mediator_stochastic_intervention_law")) {
    if (identical(node, "M1")) {
      M1_part <- prior(mediator_cache$first_M1)
      M2_part <- prior(mediator_cache$first_M2_conditional_for_first_mediator_law)
      return(list(M1_part = M1_part, M2_part = M2_part, M_part = M1_part * M2_part))
    }
    if (identical(node, "M2")) {
      M1_part <- through(mediator_cache$first_M1)
      M2_part <- prior(mediator_cache$first_M2_conditional_for_first_mediator_law)
      return(list(M1_part = M1_part, M2_part = M2_part, M_part = M1_part * M2_part))
    }
    if (identical(node, "L") || identical(node, "Y")) {
      M1_part <- through(mediator_cache$first_M1)
      M2_part <- through(mediator_cache$first_M2_conditional_for_first_mediator_law)
      return(list(M1_part = M1_part, M2_part = M2_part, M_part = M1_part * M2_part))
    }
  }
  if (identical(process_type, "second_mediator_stochastic_intervention_law")) {
    if (identical(node, "M1")) {
      M1_part <- prior(mediator_cache$second_M1)
      M2_part <- prior(mediator_cache$second_M2_conditional_for_second_mediator_law)
      return(list(M1_part = M1_part, M2_part = M2_part, M_part = M1_part * M2_part))
    }
    if (identical(node, "M2")) {
      M1_part <- through(mediator_cache$second_M1)
      M2_part <- prior(mediator_cache$second_M2_conditional_for_second_mediator_law)
      return(list(M1_part = M1_part, M2_part = M2_part, M_part = M1_part * M2_part))
    }
    if (identical(node, "L") || identical(node, "Y")) {
      M1_part <- through(mediator_cache$second_M1)
      M2_part <- through(mediator_cache$second_M2_conditional_for_second_mediator_law)
      return(list(M1_part = M1_part, M2_part = M2_part, M_part = M1_part * M2_part))
    }
  }
  list(M1_part = rep(1, n), M2_part = rep(1, n), M_part = rep(1, n))
}

.ltmle_exact_mediator_ratio_for_task <- function(task, mediator_cache) {
  .ltmle_exact_mediator_ratio_parts_for_task(task, mediator_cache)$M_part
}

.ltmle_exact_censoring_ratio_for_task <- function(task, ratio_cache) {
  out <- .ltmle_exact_cumprod_to(ratio_cache$C, task$t)
  node <- task$node
  if (startsWith(node, "L_")) node <- "L"
  if (identical(node, "Y") && task$t == ratio_cache$T && !is.null(ratio_cache$C_final)) {
    out <- out * as.numeric(ratio_cache$C_final)
  }
  out
}

.ltmle_exact_describe_ratio_factors <- function(task, T) {
  process_type <- .ltmle_exact_process_type(task)
  final_c <- if (identical(task$node, "Y") && task$t == T) {
    ";C_final:terminal-if-present"
  } else {
    ""
  }
  paste0(
    "A:", process_type, ":1-", task$t,
    ";M:", task$world_type, ":", task$node, ":1-", task$t,
    ";C_visit:1-", task$t, final_c
  )
}

.ltmle_exact_clever_covariate <- function(task, ratio_cache) {
  node <- task$node
  if (startsWith(node, "L_")) node <- "L"
  tt <- task$t
  process_type <- .ltmle_exact_process_type(task)

  A_mat <- switch(process_type,
    outcome_process = ratio_cache$A_outer,
    post_mediator_covariate_transition = ratio_cache$A_outer,
    virtual_mixed_continuation_task = ratio_cache$A_outer,
    observed_mediator_process = ratio_cache$A_outer,
    joint_stochastic_mediator_intervention_law = ratio_cache$A_joint_mediator_law,
    first_mediator_stochastic_intervention_law = ratio_cache$A_first_mediator_law,
    second_mediator_stochastic_intervention_law = ratio_cache$A_second_mediator_law,
    .stop("Unknown process_type in ltmle_exact clever covariate: ", process_type)
  )

  A_part <- .ltmle_exact_cumprod_to(A_mat, tt)
  M_parts <- .ltmle_exact_mediator_ratio_parts_for_task(task, ratio_cache$M)
  M_part <- M_parts$M_part
  C_part <- .ltmle_exact_censoring_ratio_for_task(task, ratio_cache)

  H <- as.numeric(A_part * M_part * C_part)
  if (any(!is.finite(H))) {
    .stop("Non-finite clever covariate for ", task$component, " node=", node, " t=", tt)
  }
  attr(H, "A_part") <- A_part
  attr(H, "M1_part") <- M_parts$M1_part
  attr(H, "M2_part") <- M_parts$M2_part
  attr(H, "M_part") <- M_part
  attr(H, "C_part") <- C_part
  H
}

.ltmle_exact_identity_prod_to <- function(mat, tt) {
  mat <- as.matrix(mat)
  n <- nrow(mat)
  tt <- as.integer(tt)[1L]
  if (!is.finite(tt) || tt <= 0L) return(rep(1, n))
  tt <- min(tt, ncol(mat))
  apply(mat[, seq_len(tt), drop = FALSE], 1L, prod)
}

.ltmle_exact_identity_reference_parts <- function(task, ratio_cache) {
  task <- as.list(task)
  node <- as.character(task$node %||% NA_character_)[1L]
  if (startsWith(node, "L_")) node <- "L"
  tt <- as.integer(task$t %||% NA_integer_)[1L]
  process_type <- .ltmle_exact_process_type(task)
  world_type <- as.character(task$world_type %||% NA_character_)[1L]
  n <- nrow(ratio_cache$A_outer)
  one <- rep(1, n)
  prod_to <- .ltmle_exact_identity_prod_to

  A_mat <- switch(process_type,
    outcome_process = ratio_cache$A_outer,
    post_mediator_covariate_transition = ratio_cache$A_outer,
    virtual_mixed_continuation_task = ratio_cache$A_outer,
    observed_mediator_process = ratio_cache$A_outer,
    joint_stochastic_mediator_intervention_law = ratio_cache$A_joint_mediator_law,
    first_mediator_stochastic_intervention_law = ratio_cache$A_first_mediator_law,
    second_mediator_stochastic_intervention_law = ratio_cache$A_second_mediator_law,
    ratio_cache$A_outer
  )
  A_part <- prod_to(A_mat, tt)
  M1_part <- one
  M2_part <- one
  med <- ratio_cache$M
  prior <- function(mat) prod_to(mat, tt - 1L)
  through <- function(mat) prod_to(mat, tt)

  if (process_type %in% c(
    "outcome_process",
    "post_mediator_covariate_transition",
    "virtual_mixed_continuation_task"
  )) {
    if (identical(world_type, "joint")) {
      M1_part <- through(med$joint_M1)
      M2_part <- through(med$joint_M2)
    } else if (identical(world_type, "separate")) {
      M1_part <- through(med$first_M1)
      M2_part <- through(med$second_M2_marginal_for_outcome_process)
    }
  } else if (identical(process_type, "joint_stochastic_mediator_intervention_law")) {
    if (identical(node, "M1")) {
      M1_part <- prior(med$joint_M1)
      M2_part <- prior(med$joint_M2)
    } else if (identical(node, "M2")) {
      M1_part <- through(med$joint_M1)
      M2_part <- prior(med$joint_M2)
    } else if (identical(node, "L") || identical(node, "Y")) {
      M1_part <- through(med$joint_M1)
      M2_part <- through(med$joint_M2)
    }
  } else if (identical(process_type, "first_mediator_stochastic_intervention_law")) {
    if (identical(node, "M1")) {
      M1_part <- prior(med$first_M1)
      M2_part <- prior(med$first_M2_conditional_for_first_mediator_law)
    } else if (identical(node, "M2")) {
      M1_part <- through(med$first_M1)
      M2_part <- prior(med$first_M2_conditional_for_first_mediator_law)
    } else if (identical(node, "L") || identical(node, "Y")) {
      M1_part <- through(med$first_M1)
      M2_part <- through(med$first_M2_conditional_for_first_mediator_law)
    }
  } else if (identical(process_type, "second_mediator_stochastic_intervention_law")) {
    if (identical(node, "M1")) {
      M1_part <- prior(med$second_M1)
      M2_part <- prior(med$second_M2_conditional_for_second_mediator_law)
    } else if (identical(node, "M2")) {
      M1_part <- through(med$second_M1)
      M2_part <- prior(med$second_M2_conditional_for_second_mediator_law)
    } else if (identical(node, "L") || identical(node, "Y")) {
      M1_part <- through(med$second_M1)
      M2_part <- through(med$second_M2_conditional_for_second_mediator_law)
    }
  }

  C_part <- prod_to(ratio_cache$C, tt)
  if (identical(node, "Y") && identical(tt, as.integer(ratio_cache$T)) &&
      !is.null(ratio_cache$C_final)) {
    C_part <- C_part * as.numeric(ratio_cache$C_final)
  }
  H <- as.numeric(A_part * M1_part * M2_part * C_part)
  list(
    A_part = A_part,
    M1_part = M1_part,
    M2_part = M2_part,
    C_part = C_part,
    H = H
  )
}

.ltmle_exact_build_independent_separate_h_reference <- function(task,
                                                                long,
                                                                models,
                                                                treatment_models,
                                                                censoring_models,
                                                                T,
                                                                spec,
                                                                probability_bounds,
                                                                treat_mech,
                                                                p_rct,
                                                                node_spec = NULL,
                                                                n_particles = 2000L,
                                                                seed = 202405L,
                                                                eps = 1e-12,
                                                                max_t = NULL) {
  task <- as.list(task)
  node <- as.character(task$node %||% NA_character_)[1L]
  if (startsWith(node, "L_")) node <- "L"
  tt <- as.integer(task$t %||% NA_integer_)[1L]
  if (!is.finite(tt) || tt < 1L) tt <- 1L
  max_t <- if (is.null(max_t)) tt else min(as.integer(max_t), T)
  if (!is.finite(max_t) || max_t < 1L) max_t <- 1L
  process_type <- .ltmle_exact_process_type(task)
  world_type <- as.character(task$world_type %||% spec$world_type[1L] %||% NA_character_)[1L]
  ids <- sort(unique(long$id))
  n <- length(ids)
  one <- rep(1, n)
  prod_to <- function(mat, upto) .ltmle_exact_identity_prod_to(mat, upto)
  prior <- function(mat) prod_to(mat, tt - 1L)
  through <- function(mat) prod_to(mat, tt)

  direct_treatment_ratio_matrix <- function(regimen) {
    regimen <- normalize_regimen(regimen, T, "independent_reference_regimen")
    out <- matrix(1, nrow = n, ncol = T)
    if (identical(treat_mech, "baseline_rct")) {
      if (any(regimen != regimen[1L])) {
        .stop("baseline_rct in independent H reference requires a constant treatment regimen.")
      }
      dt1 <- long[long$t == 1L, , drop = FALSE]
      dt1 <- dt1[match(ids, dt1$id), , drop = FALSE]
      p <- if (regimen[1L] == 1) p_rct else 1 - p_rct
      out[, 1L] <- as.numeric(dt1$A == regimen[1L]) / pmax(p, 1e-8)
      if (T >= 2L) out[, 2:T] <- 1
      return(out)
    }
    for (jj in seq_len(max_t)) {
      dt <- long[long$t == jj, , drop = FALSE]
      dt <- dt[match(ids, dt$id), , drop = FALSE]
      p1 <- .ltmle_exact_predict_nuis(treatment_models$fits[[jj]], dt, type = "response", eps = 1e-8)
      p1 <- .ltmle_exact_bound_probability(p1, probability_bounds)
      p_reg <- ifelse(regimen[jj] == 1, p1, 1 - p1)
      p_reg <- .ltmle_exact_bound_probability(p_reg, probability_bounds)
      out[, jj] <- as.numeric(dt$A == regimen[jj]) / pmax(p_reg, 1e-8)
    }
    out
  }

  direct_mediator_parts <- function() {
    out <- list(
      joint_M1 = matrix(1, n, T),
      joint_M2 = matrix(1, n, T),
      first_M1 = matrix(1, n, T),
      first_M2_conditional_for_first_mediator_law = matrix(1, n, T),
      second_M1 = matrix(1, n, T),
      second_M2_marginal_for_outcome_process = matrix(1, n, T),
      second_M2_conditional_for_second_mediator_law = matrix(1, n, T)
    )
    second_marginal <- matrix(1, nrow = n, ncol = T)
    if (identical(world_type, "separate") && max_t >= 2L) {
      second_marginal <- .ltmle_exact_second_mediator_law_density_ratio(
        observed_long = long,
        models = models,
        T = T,
        second_regimen = as.numeric(spec$m2_A),
        n_particles = n_particles,
        seed = seed,
        eps = eps,
        node_spec = node_spec,
        max_t = max_t
      )
    }
    for (jj in seq_len(max_t)) {
      if (jj < 2L) next
      dt <- long[long$t == jj, , drop = FALSE]
      dt <- dt[match(ids, dt$id), , drop = FALSE]
      dt <- .ltmle_exact_add_terms(dt)
      sp <- spec[spec$t == jj, , drop = FALSE]
      if (!nrow(sp)) sp <- spec[1L, , drop = FALSE]

      den_M1_hist <- .ltmle_exact_set_A(dt, dt$A)
      den_M2_hist <- .ltmle_exact_set_A(dt, dt$A)
      den_M1_density <- .ltmle_exact_predict_mediator_density(models$M1[[jj]], dt$M1, den_M1_hist, eps)
      den_M2_density <- .ltmle_exact_predict_mediator_density(models$M2[[jj]], dt$M2, den_M2_hist, eps)

      if (identical(world_type, "joint")) {
        joint_M1_hist <- .ltmle_exact_set_A(dt, sp$m1_A[1L])
        joint_M2_hist <- .ltmle_exact_set_A(dt, sp$m2_A[1L])
        joint_M2_hist$M1 <- dt$M1
        joint_M2_hist <- .ltmle_exact_add_terms(joint_M2_hist)
        out$joint_M1[, jj] <-
          .ltmle_exact_predict_mediator_density(models$M1[[jj]], dt$M1, joint_M1_hist, eps) /
          den_M1_density
        out$joint_M2[, jj] <-
          .ltmle_exact_predict_mediator_density(models$M2[[jj]], dt$M2, joint_M2_hist, eps) /
          den_M2_density
      }

      if (identical(world_type, "separate")) {
        first_M1_hist <- .ltmle_exact_set_A(dt, sp$m1_A[1L])
        out$first_M1[, jj] <-
          .ltmle_exact_predict_mediator_density(models$M1[[jj]], dt$M1, first_M1_hist, eps) /
          den_M1_density

        first_M2_hist <- .ltmle_exact_set_A(dt, sp$m1_A[1L])
        first_M2_hist$M1 <- dt$M1
        first_M2_hist <- .ltmle_exact_add_terms(first_M2_hist)
        out$first_M2_conditional_for_first_mediator_law[, jj] <-
          .ltmle_exact_predict_mediator_density(models$M2[[jj]], dt$M2, first_M2_hist, eps) /
          den_M2_density

        second_M1_hist <- .ltmle_exact_set_A(dt, sp$m2_A[1L])
        out$second_M1[, jj] <-
          .ltmle_exact_predict_mediator_density(models$M1[[jj]], dt$M1, second_M1_hist, eps) /
          den_M1_density

        second_M2_hist <- .ltmle_exact_set_A(dt, sp$m2_A[1L])
        second_M2_hist$M1 <- dt$M1
        second_M2_hist <- .ltmle_exact_add_terms(second_M2_hist)
        out$second_M2_conditional_for_second_mediator_law[, jj] <-
          .ltmle_exact_predict_mediator_density(models$M2[[jj]], dt$M2, second_M2_hist, eps) /
          den_M2_density
        out$second_M2_marginal_for_outcome_process[, jj] <- second_marginal[, jj]
      }
    }
    out
  }

  A_outer <- direct_treatment_ratio_matrix(as.numeric(spec$outer_A))
  A_first <- direct_treatment_ratio_matrix(as.numeric(spec$m1_A))
  A_second <- direct_treatment_ratio_matrix(as.numeric(spec$m2_A))
  A_joint <- if (identical(world_type, "joint")) A_first else matrix(1, nrow = n, ncol = T)
  A_mat <- switch(process_type,
    outcome_process = A_outer,
    post_mediator_covariate_transition = A_outer,
    virtual_mixed_continuation_task = A_outer,
    observed_mediator_process = A_outer,
    joint_stochastic_mediator_intervention_law = A_joint,
    first_mediator_stochastic_intervention_law = A_first,
    second_mediator_stochastic_intervention_law = A_second,
    A_outer
  )
  A_part <- prod_to(A_mat, tt)

  med <- direct_mediator_parts()
  M1_part <- one
  M2_part <- one
  if (process_type %in% c(
    "outcome_process",
    "post_mediator_covariate_transition",
    "virtual_mixed_continuation_task"
  )) {
    if (identical(world_type, "joint")) {
      M1_part <- through(med$joint_M1)
      M2_part <- through(med$joint_M2)
    } else if (identical(world_type, "separate")) {
      M1_part <- through(med$first_M1)
      M2_part <- through(med$second_M2_marginal_for_outcome_process)
    }
  } else if (identical(process_type, "joint_stochastic_mediator_intervention_law")) {
    if (identical(node, "M1")) {
      M1_part <- prior(med$joint_M1)
      M2_part <- prior(med$joint_M2)
    } else if (identical(node, "M2")) {
      M1_part <- through(med$joint_M1)
      M2_part <- prior(med$joint_M2)
    } else if (identical(node, "L") || identical(node, "Y")) {
      M1_part <- through(med$joint_M1)
      M2_part <- through(med$joint_M2)
    }
  } else if (identical(process_type, "first_mediator_stochastic_intervention_law")) {
    if (identical(node, "M1")) {
      M1_part <- prior(med$first_M1)
      M2_part <- prior(med$first_M2_conditional_for_first_mediator_law)
    } else if (identical(node, "M2")) {
      M1_part <- through(med$first_M1)
      M2_part <- prior(med$first_M2_conditional_for_first_mediator_law)
    } else if (identical(node, "L") || identical(node, "Y")) {
      M1_part <- through(med$first_M1)
      M2_part <- through(med$first_M2_conditional_for_first_mediator_law)
    }
  } else if (identical(process_type, "second_mediator_stochastic_intervention_law")) {
    if (identical(node, "M1")) {
      M1_part <- prior(med$second_M1)
      M2_part <- prior(med$second_M2_conditional_for_second_mediator_law)
    } else if (identical(node, "M2")) {
      M1_part <- through(med$second_M1)
      M2_part <- prior(med$second_M2_conditional_for_second_mediator_law)
    } else if (identical(node, "L") || identical(node, "Y")) {
      M1_part <- through(med$second_M1)
      M2_part <- through(med$second_M2_conditional_for_second_mediator_law)
    }
  }

  C_direct <- .ltmle_exact_censoring_ratio_matrix(
    long,
    censoring_models %||% list(type = "none", T = T),
    probability_bounds = probability_bounds,
    max_t = max_t
  )
  C_part <- prod_to(C_direct$visit, tt)
  if (identical(node, "Y") && identical(tt, as.integer(T)) &&
      !is.null(C_direct$final)) {
    C_part <- C_part * as.numeric(C_direct$final)
  }

  H <- as.numeric(A_part * M1_part * M2_part * C_part)
  list(
    A_part = A_part,
    M1_part = M1_part,
    M2_part = M2_part,
    C_part = C_part,
    H = H,
    reference_source = "independent_separate_eif_reference",
    reference_builder = ".ltmle_exact_build_independent_separate_h_reference",
    reference_uses_ratio_cache = FALSE,
    reference_uses_targeting_H_attributes = FALSE,
    reference_uses_ltmle_mediator_ratio_helper = FALSE,
    reference_uses_ltmle_density_ratio_cache = FALSE
  )
}

.ltmle_exact_regimen_label <- function(x) {
  x <- as.numeric(x)
  if (!length(x)) return("unknown")
  if (all(x == 1, na.rm = TRUE)) return("a")
  if (all(x == 0, na.rm = TRUE)) return("as")
  paste(signif(x, 8), collapse = "|")
}

.ltmle_exact_expected_separate_labels <- function(component) {
  switch(as.character(component)[1L],
    mu_sep_aaa = c(outer = "a", m1 = "a", m2 = "a"),
    mu_sep_asas_asas = c(outer = "as", m1 = "as", m2 = "as"),
    mu_sep_a_asas = c(outer = "a", m1 = "as", m2 = "as"),
    mu_sep_a_aas = c(outer = "a", m1 = "a", m2 = "as"),
    c(outer = NA_character_, m1 = NA_character_, m2 = NA_character_)
  )
}

.ltmle_exact_theoretical_eif_term_label <- function(task) {
  task <- as.list(task)
  process_type <- .ltmle_exact_process_type(task)
  world_type <- as.character(task$world_type %||% NA_character_)[1L]
  node <- as.character(task$node %||% NA_character_)[1L]
  if (startsWith(node, "L_")) node <- "L"
  if (identical(world_type, "separate") &&
      process_type %in% c("outcome_process", "post_mediator_covariate_transition",
                          "virtual_mixed_continuation_task")) {
    return("A_outer * M1_first_law_through_t * M2_second_law_marginal_for_outcome_through_t * C")
  }
  paste("A", process_type, node, "mediator_factorization", "C", sep = " * ")
}

.ltmle_exact_clever_covariate_identity_row <- function(task,
                                                       H_used,
                                                       H_truncated_used = NULL,
                                                       ratio_cache = NULL,
                                                       reference_parts = NULL,
                                                       T,
                                                       spec,
                                                       row_evaluation_context,
                                                       truncation = NULL,
                                                       estimator_variant = "ltmle_exact_quantile_truncated") {
  task <- as.list(task)
  used_A <- as.numeric(attr(H_used, "A_part") %||% rep(NA_real_, length(H_used)))
  used_M <- as.numeric(attr(H_used, "M_part") %||% rep(NA_real_, length(H_used)))
  used_M1 <- as.numeric(attr(H_used, "M1_part") %||% rep(NA_real_, length(H_used)))
  used_M2 <- as.numeric(attr(H_used, "M2_part") %||% rep(NA_real_, length(H_used)))
  used_C <- as.numeric(attr(H_used, "C_part") %||% rep(NA_real_, length(H_used)))
  ref <- reference_parts %||% NULL
  H_used <- as.numeric(H_used)
  H_truncated_used <- as.numeric(H_truncated_used %||% H_used)
  if (is.null(ref)) {
    ref <- .ltmle_exact_identity_reference_parts(task, ratio_cache)
    ref$reference_source <- "ratio_cache"
    ref$reference_builder <- ".ltmle_exact_identity_reference_parts"
    ref$reference_uses_ratio_cache <- TRUE
    ref$reference_uses_targeting_H_attributes <- FALSE
    ref$reference_uses_ltmle_mediator_ratio_helper <- TRUE
    ref$reference_uses_ltmle_density_ratio_cache <- TRUE
  }
  reference_source <- as.character(ref$reference_source %||% "unknown_reference_source")[1L]
  reference_uses_ratio_cache <- isTRUE(ref$reference_uses_ratio_cache)
  reference_uses_targeting_H_attributes <- isTRUE(ref$reference_uses_targeting_H_attributes)
  reference_uses_ltmle_mediator_ratio_helper <- isTRUE(ref$reference_uses_ltmle_mediator_ratio_helper)
  reference_uses_ltmle_density_ratio_cache <- isTRUE(ref$reference_uses_ltmle_density_ratio_cache)
  independent_reference_used <- identical(reference_source, "independent_separate_eif_reference") &&
    !reference_uses_ratio_cache &&
    !reference_uses_targeting_H_attributes &&
    !reference_uses_ltmle_mediator_ratio_helper &&
    !reference_uses_ltmle_density_ratio_cache
  truncate_reference <- function(x) {
    x <- as.numeric(x)
    if (isTRUE(truncation$enabled) && identical(truncation$policy, "quantile")) {
      q_lower <- .ltmle_exact_quantile_value(x, truncation$quantile_lower)
      q_upper <- .ltmle_exact_quantile_value(x, truncation$quantile_upper)
      if (is.finite(q_lower) && is.finite(q_upper) && q_lower <= q_upper) {
        return(pmin(pmax(x, q_lower), q_upper))
      }
    }
    x
  }
  H_reference_truncated <- truncate_reference(ref$H)
  tol <- max(1e-8, 1e-7 * max(1, max(abs(ref$H), abs(H_used), na.rm = TRUE)))
  raw_abs_diff <- max(abs(H_used - ref$H), na.rm = TRUE)
  trunc_tol <- max(
    1e-8,
    1e-7 * max(1, max(abs(H_reference_truncated), abs(H_truncated_used), na.rm = TRUE))
  )
  trunc_abs_diff <- max(abs(H_truncated_used - H_reference_truncated), na.rm = TRUE)
  abs_diff <- raw_abs_diff
  rel_diff <- abs_diff / max(1, max(abs(ref$H), na.rm = TRUE))
  trunc_rel_diff <- trunc_abs_diff / max(1, max(abs(H_reference_truncated), na.rm = TRUE))
  A_diff <- max(abs(used_A - ref$A_part), na.rm = TRUE)
  M_diff <- max(abs(used_M - ref$M1_part * ref$M2_part), na.rm = TRUE)
  M1_diff <- max(abs(used_M1 - ref$M1_part), na.rm = TRUE)
  M2_diff <- max(abs(used_M2 - ref$M2_part), na.rm = TRUE)
  C_diff <- max(abs(used_C - ref$C_part), na.rm = TRUE)
  expected <- .ltmle_exact_expected_separate_labels(task$component)
  outer_regimen <- .ltmle_exact_regimen_label(spec$outer_A)
  m1_regimen <- .ltmle_exact_regimen_label(spec$m1_A)
  m2_regimen <- .ltmle_exact_regimen_label(spec$m2_A)
  is_separate <- identical(as.character(task$world_type %||% NA_character_)[1L], "separate")
  factor_regimen_assignment_correct <- if (isTRUE(is_separate)) {
    identical(unname(c(outer_regimen, m1_regimen, m2_regimen)), unname(expected))
  } else {
    TRUE
  }
  H_match <- is.finite(raw_abs_diff) && raw_abs_diff <= tol
  H_truncated_match <- is.finite(trunc_abs_diff) && trunc_abs_diff <= trunc_tol
  A_match <- is.finite(A_diff) && A_diff <= tol
  M_match <- is.finite(M_diff) && M_diff <= tol
  M1_match <- is.finite(M1_diff) && M1_diff <= tol
  M2_match <- is.finite(M2_diff) && M2_diff <= tol
  C_match <- is.finite(C_diff) && C_diff <= tol
  data.frame(
    component = task$component %||% NA_character_,
    task_id = task$task_id %||% NA_character_,
    t = as.integer(task$t %||% NA_integer_),
    node = task$node %||% NA_character_,
    process_type = .ltmle_exact_process_type(task),
    row_evaluation_context = row_evaluation_context,
    clever_covariate_term_label = .ltmle_exact_describe_ratio_factors(task, T),
    theoretical_eif_term_label = .ltmle_exact_theoretical_eif_term_label(task),
    reference_source = reference_source,
    reference_builder = as.character(ref$reference_builder %||% NA_character_)[1L],
    reference_uses_ratio_cache = reference_uses_ratio_cache,
    reference_uses_targeting_H_attributes = reference_uses_targeting_H_attributes,
    reference_uses_ltmle_mediator_ratio_helper = reference_uses_ltmle_mediator_ratio_helper,
    reference_uses_ltmle_density_ratio_cache = reference_uses_ltmle_density_ratio_cache,
    H_raw_used_before_truncation = mean(H_used, na.rm = TRUE),
    H_reference_independent_before_truncation = mean(ref$H, na.rm = TRUE),
    H_raw_abs_difference = raw_abs_diff,
    H_raw_relative_difference = rel_diff,
    H_raw_identity_matches = H_match,
    H_truncated_used_in_targeting = mean(H_truncated_used, na.rm = TRUE),
    H_truncated_reference_after_quantile_rule = mean(H_reference_truncated, na.rm = TRUE),
    H_truncated_abs_difference = trunc_abs_diff,
    H_truncated_relative_difference = trunc_rel_diff,
    H_truncated_identity_matches = H_truncated_match,
    H_used_in_targeting = mean(H_used, na.rm = TRUE),
    H_reference_independent = mean(ref$H, na.rm = TRUE),
    H_abs_difference = abs_diff,
    H_relative_difference = rel_diff,
    H_identity_matches = H_match,
    A_part_used = mean(used_A, na.rm = TRUE),
    M1_part_used = mean(used_M1, na.rm = TRUE),
    M2_part_used = mean(used_M2, na.rm = TRUE),
    C_part_used = mean(used_C, na.rm = TRUE),
    A_part_reference = mean(ref$A_part, na.rm = TRUE),
    M1_part_reference = mean(ref$M1_part, na.rm = TRUE),
    M2_part_reference = mean(ref$M2_part, na.rm = TRUE),
    C_part_reference = mean(ref$C_part, na.rm = TRUE),
    A_part_matches = A_match,
    M_part_matches = M_match,
    M1_part_matches = M1_match,
    M2_part_matches = M2_match,
    C_part_matches = C_match,
    first_law_id = paste0(task$component %||% NA_character_, "::first_law"),
    second_law_id = paste0(task$component %||% NA_character_, "::second_law"),
    product_join_event_id = NA_integer_,
    first_particle_id = NA_integer_,
    second_particle_id = NA_integer_,
    first_weight_raw = NA_real_,
    second_weight_raw = NA_real_,
    joined_weight_raw = NA_real_,
    joined_weight_normalized = NA_real_,
    H_factor_uses_joined_row_ownership = NA,
    outer_regimen = outer_regimen,
    m1_regimen = m1_regimen,
    m2_regimen = m2_regimen,
    factor_regimen_assignment_correct = factor_regimen_assignment_correct,
    product_join_weight_consistent_with_eif = TRUE,
    independent_reference_used = independent_reference_used,
    tolerance = tol,
    passed = independent_reference_used && H_match && H_truncated_match &&
      A_match && M_match && M1_match && M2_match && C_match &&
      isTRUE(factor_regimen_assignment_correct),
    failure_class = if (independent_reference_used && H_match && H_truncated_match &&
                        A_match && M_match && M1_match && M2_match && C_match &&
                        isTRUE(factor_regimen_assignment_correct)) {
      "no_failure"
    } else if (!independent_reference_used) {
      "independent_reference_not_used"
    } else if (!H_match) {
      "clever_covariate_identity_mismatch"
    } else if (!H_truncated_match) {
      "clever_covariate_truncation_identity_mismatch"
    } else if (!isTRUE(factor_regimen_assignment_correct)) {
      "factor_regimen_assignment_mismatch"
    } else {
      "clever_covariate_factor_mismatch"
    },
    stringsAsFactors = FALSE
  )
}

.ltmle_exact_ratio_stats <- function(x, task, ratio_name) {
  x <- as.numeric(x)
  good <- is.finite(x)
  xx <- x[good]
  if (!length(xx)) xx <- NA_real_
  data.frame(
    component = task$component,
    task_id = task$task_id %||% NA_character_,
    t = task$t,
    node = task$node,
    process_type = .ltmle_exact_process_type(task),
    ratio_name = ratio_name,
    min_ratio = min(xx, na.rm = TRUE),
    p01_ratio = as.numeric(stats::quantile(xx, 0.01, na.rm = TRUE, names = FALSE)),
    median_ratio = stats::median(xx, na.rm = TRUE),
    p95_ratio = as.numeric(stats::quantile(xx, 0.95, na.rm = TRUE, names = FALSE)),
    p99_ratio = as.numeric(stats::quantile(xx, 0.99, na.rm = TRUE, names = FALSE)),
    max_ratio = max(xx, na.rm = TRUE),
    mean_ratio = mean(xx, na.rm = TRUE),
    sd_ratio = stats::sd(xx, na.rm = TRUE),
    effective_sample_size = .ltmle_exact_ess(abs(xx)),
    fixed_bound_truncation_enabled = FALSE,
    fixed_bound_truncation_used = FALSE,
    density_ratio_factor_truncation_used = FALSE,
    ratio_was_bounded = FALSE,
    n_bounded_lower = 0L,
    n_bounded_upper = 0L,
    stringsAsFactors = FALSE
  )
}

.ltmle_exact_density_ratio_diagnostics_for_task <- function(task, ratio_cache) {
  process_type <- .ltmle_exact_process_type(task)
  A_mat <- switch(process_type,
    outcome_process = ratio_cache$A_outer,
    post_mediator_covariate_transition = ratio_cache$A_outer,
    virtual_mixed_continuation_task = ratio_cache$A_outer,
    observed_mediator_process = ratio_cache$A_outer,
    joint_stochastic_mediator_intervention_law = ratio_cache$A_joint_mediator_law,
    first_mediator_stochastic_intervention_law = ratio_cache$A_first_mediator_law,
    second_mediator_stochastic_intervention_law = ratio_cache$A_second_mediator_law,
    ratio_cache$A_outer
  )
  tt <- task$t
  rows <- list(
    .ltmle_exact_ratio_stats(.ltmle_exact_cumprod_to(A_mat, tt), task,
                             "treatment likelihood factor"),
    .ltmle_exact_ratio_stats(.ltmle_exact_censoring_ratio_for_task(task, ratio_cache), task,
                             "censoring factor"),
    .ltmle_exact_ratio_stats(.ltmle_exact_mediator_ratio_for_task(task, ratio_cache$M), task,
                             "mediator density ratio factor")
  )
  if (tt >= 2L) {
    rows[[length(rows) + 1L]] <- .ltmle_exact_ratio_stats(ratio_cache$M$joint_M1[, tt], task,
                                                          "joint mediator M1 density ratio")
    rows[[length(rows) + 1L]] <- .ltmle_exact_ratio_stats(ratio_cache$M$joint_M2[, tt], task,
                                                          "joint mediator M2 density ratio")
    rows[[length(rows) + 1L]] <- .ltmle_exact_ratio_stats(ratio_cache$M$first_M1[, tt], task,
                                                          "first-mediator stochastic intervention law M1 density ratio")
    rows[[length(rows) + 1L]] <- .ltmle_exact_ratio_stats(
      ratio_cache$M$first_M2_conditional_for_first_mediator_law[, tt], task,
      "first-mediator stochastic intervention law auxiliary M2 conditional density ratio")
    rows[[length(rows) + 1L]] <- .ltmle_exact_ratio_stats(
      ratio_cache$M$second_M2_marginal_for_outcome_process[, tt], task,
      "second-mediator stochastic intervention law M2 marginal density ratio for outcome process")
    rows[[length(rows) + 1L]] <- .ltmle_exact_ratio_stats(
      ratio_cache$M$second_M2_conditional_for_second_mediator_law[, tt], task,
      "second-mediator stochastic intervention law M2 conditional density ratio")
  }
  do.call(rbind, rows)
}

.ltmle_exact_clever_covariate_for_rows <- function(task,
                                                   rows,
                                                   models,
                                                   treatment_models,
                                                   censoring_models,
                                                   T,
                                                   spec,
                                                   probability_bounds,
                                                   treat_mech,
                                                   p_rct,
                                                   node_spec = NULL,
                                                   row_long = NULL,
                                                   row_source = NULL,
                                                   ltmle_exact_density_ratio_mc_n = 2000L,
                                                   max_t = NULL) {
  task <- as.list(task)
  row_long <- row_long %||% attr(rows, "ltmle_exact_long")
  if (is.null(row_long)) {
    row_long <- rows
    row_long$id <- seq_len(nrow(row_long))
    row_long$t <- task$t
  }
  row_source <- row_source %||% attr(rows, "ltmle_exact_row_source") %||% "observed"
  row_source <- match.arg(row_source, c("observed", "generated"))
  ratio_cache <- .ltmle_exact_build_ratio_cache(
    long = row_long,
    models = models,
    treatment_models = treatment_models,
    T = T,
    spec = spec,
    probability_bounds = probability_bounds,
    treat_mech = treat_mech,
    p_rct = p_rct,
    censoring_models = censoring_models,
    history_source = row_source,
    node_spec = node_spec,
    ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
    max_t = max_t %||% task$t
  )
  H <- .ltmle_exact_clever_covariate(task, ratio_cache)
  diagnostics <- .ltmle_exact_density_ratio_diagnostics_for_task(task, ratio_cache)
  attr(H, "density_ratio_diagnostics") <- diagnostics
  H
}

.ltmle_exact_draw_L_nodes <- function(row, models, tt, node_spec = NULL) {
  for (L_node in .ltmle_exact_L_nodes(node_spec)) {
    fit <- models$L[[L_node]][[tt]]
    mu <- .ltmle_exact_predict_nuis(fit, row, type = "numeric")
    sig <- .ltmle_exact_sigma(fit)
    row[[L_node]] <- stats::rnorm(nrow(row), mu, sig)
    if (identical(L_node, "L")) row$L <- row[[L_node]]
    row <- .ltmle_exact_add_terms(row)
  }
  row
}

.ltmle_exact_draw_y_history <- function(row, fit_Y) {
  mu <- .ltmle_exact_predict_nuis(fit_Y, row, type = "numeric")
  sig <- .ltmle_exact_sigma(fit_Y)
  stats::rnorm(nrow(row), mu, sig)
}

.ltmle_exact_generate_one_law_history <- function(history, tt, A_value_M1, A_value_M2, models,
                                                 law_type, M1_override = NULL,
                                                 node_spec = NULL) {
  row <- data.frame(
    id0 = history$id0,
    W1 = history$W1,
    W2 = history$W2,
    Y0 = history$Y0,
    M1_0 = history$M1_0,
    M2_0 = history$M2_0,
    M1_lag = history$M1_lag,
    M2_lag = history$M2_lag,
    A = A_value_M1,
    stringsAsFactors = FALSE
  )
  if ("Y_lag" %in% names(history)) row$Y_lag <- history$Y_lag
  for (L_node in .ltmle_exact_L_nodes(node_spec)) {
    lag_nm <- .ltmle_exact_lag_name(L_node)
    if (lag_nm %in% names(history)) row[[lag_nm]] <- history[[lag_nm]]
  }
  extra_baseline <- setdiff(names(history), c(names(row), "id0", "M1_0", "M2_0",
                                             "Y_lag", "L_lag", "M1_lag", "M2_lag"))
  for (nm in extra_baseline) row[[nm]] <- history[[nm]]

  if (tt == 1L) {
    row$M1 <- row$M1_0
    row$M2 <- row$M2_0
  } else {
    if (!is.null(M1_override)) {
      row$M1 <- as.numeric(M1_override)
    } else {
      M1_hist <- .ltmle_exact_set_A(row, A_value_M1)
      row$M1 <- .ltmle_exact_draw_mediator(models$M1[[tt]], M1_hist)
    }

    M2_hist <- .ltmle_exact_set_A(row, A_value_M2)
    M2_hist$M1 <- row$M1
    M2_hist <- .ltmle_exact_add_terms(M2_hist)
    row[["M2"]] <- .ltmle_exact_draw_mediator(models$M2[[tt]], M2_hist)
  }

  row$A <- A_value_M2
  row <- .ltmle_exact_add_terms(row)
  row <- .ltmle_exact_draw_L_nodes(row, models, tt, node_spec)
  row <- .ltmle_exact_add_terms(row)
  draw_y <- identical(.ltmle_exact_outcome_type(node_spec), "longitudinal") && tt < length(models$Y)
  if (isTRUE(draw_y) && !is.null(models$Y[[tt]])) {
    y_next <- .ltmle_exact_draw_y_history(row, models$Y[[tt]])
    row$Y <- y_next
  } else {
    y_next <- NULL
    row$Y <- NA_real_
  }

  next_history <- history
  next_history$M1_lag <- row$M1
  next_history$M2_lag <- row$M2
  for (L_node in .ltmle_exact_L_nodes(node_spec)) {
    next_history[[.ltmle_exact_lag_name(L_node)]] <- row[[L_node]]
  }
  if (!identical(.ltmle_exact_outcome_type(node_spec), "terminal_only") && !is.null(y_next)) {
    next_history$Y_lag <- y_next
  }

  list(row = row, next_history = next_history)
}

.ltmle_exact_assemble_outcome_process_row <- function(base_history,
                                                     tt,
                                                     spec_row,
                                                     joint_row,
                                                     first_row,
                                                     second_row,
                                                     models,
                                                     node_spec = NULL) {
  row <- data.frame(
    id0 = base_history$id0,
    W1 = base_history$W1,
    W2 = base_history$W2,
    Y0 = base_history$Y0,
    M1_0 = base_history$M1_0,
    M2_0 = base_history$M2_0,
    M1_lag = base_history$M1_lag,
    M2_lag = base_history$M2_lag,
    A = spec_row$outer_A[1L],
    stringsAsFactors = FALSE
  )
  if ("Y_lag" %in% names(base_history)) row$Y_lag <- base_history$Y_lag
  for (L_node in .ltmle_exact_L_nodes(node_spec)) {
    lag_nm <- .ltmle_exact_lag_name(L_node)
    if (lag_nm %in% names(base_history)) row[[lag_nm]] <- base_history[[lag_nm]]
  }
  extra_baseline <- setdiff(names(base_history), c(names(row), "id0", "M1_0", "M2_0",
                                                   "Y_lag", "L_lag", "M1_lag", "M2_lag"))
  for (nm in extra_baseline) row[[nm]] <- base_history[[nm]]

  if (tt == 1L) {
    row$M1 <- row$M1_0
    row$M2 <- row$M2_0
  } else if (identical(spec_row$world_type[1L], "joint")) {
    row$M1 <- joint_row$M1
    row$M2 <- joint_row$M2
  } else if (identical(spec_row$world_type[1L], "separate")) {
    row$M1 <- first_row$M1
    row$M2 <- second_row$M2
  } else if (identical(spec_row$world_type[1L], "natural")) {
    natural_row <- .ltmle_exact_generate_one_law_history(
      history = base_history,
      tt = tt,
      A_value_M1 = spec_row$outer_A[1L],
      A_value_M2 = spec_row$outer_A[1L],
      models = models,
      law_type = "observed_mediator_process",
      node_spec = node_spec
    )$row
    row$M1 <- natural_row$M1
    row$M2 <- natural_row$M2
  } else {
    .stop("Unknown world_type in outcome process assembly: ", spec_row$world_type[1L])
  }

  row <- .ltmle_exact_add_terms(row)
  row <- .ltmle_exact_draw_L_nodes(row, models, tt, node_spec)
  row <- .ltmle_exact_add_terms(row)
  draw_y <- identical(.ltmle_exact_outcome_type(node_spec), "longitudinal") || tt == length(models$Y)
  if (isTRUE(draw_y) && !is.null(models$Y[[tt]])) {
    y_next <- .ltmle_exact_draw_y_history(row, models$Y[[tt]])
    row$Y <- y_next
  } else {
    y_next <- NULL
    row$Y <- NA_real_
  }

  next_history <- base_history
  next_history$M1_lag <- row$M1
  next_history$M2_lag <- row$M2
  for (L_node in .ltmle_exact_L_nodes(node_spec)) {
    next_history[[.ltmle_exact_lag_name(L_node)]] <- row[[L_node]]
  }
  if (!identical(.ltmle_exact_outcome_type(node_spec), "terminal_only") && !is.null(y_next)) {
    next_history$Y_lag <- y_next
  }

  list(row = row, next_history = next_history)
}

.ltmle_exact_init_mc_history <- function(base, node_spec = NULL) {
  out <- data.frame(
    id0 = base$id0,
    W1 = base$W1,
    W2 = base$W2,
    Y0 = base$Y0,
    M1_0 = base$M1_0,
    M2_0 = base$M2_0,
    M1_lag = base$M1_0,
    M2_lag = base$M2_0,
    stringsAsFactors = FALSE
  )
  if (!identical(.ltmle_exact_outcome_type(node_spec), "terminal_only")) out$Y_lag <- base$Y0
  for (L_node in .ltmle_exact_L_nodes(node_spec)) {
    init_spec <- node_spec$L_lag_init[[sub("^L_", "", L_node)]] %||% 0
    lag_nm <- .ltmle_exact_lag_name(L_node)
    if (is.character(init_spec) && length(init_spec) == 1L) {
      if (!init_spec %in% names(base)) .stop("Missing L_lag_init column for ", L_node, ": ", init_spec)
      out[[lag_nm]] <- as.numeric(base[[init_spec]])
    } else if (is.numeric(init_spec) && length(init_spec) == 1L) {
      out[[lag_nm]] <- rep(as.numeric(init_spec), nrow(base))
    } else {
      .stop("Invalid L_lag_init for ", L_node)
    }
  }
  for (nm in as.character(node_spec$baseline_vars %||% character(0))) {
    if (nm %in% names(base) && !nm %in% names(out)) out[[nm]] <- as.numeric(base[[nm]])
  }
  if (".fold_id" %in% names(base)) out$.fold_id <- base$.fold_id
  out
}

.ltmle_exact_draw_component_histories <- function(dat_wide, T, spec, models, mc_n, seed,
                                                 node_spec = NULL) {
  if (!is.null(seed)) set.seed(seed)
  base <- dat_wide[rep(seq_len(nrow(dat_wide)), each = mc_n), , drop = FALSE]
  base$id0 <- rep(seq_len(nrow(dat_wide)), each = mc_n)

  histories <- list(
    outcome_process = .ltmle_exact_init_mc_history(base, node_spec),
    joint_stochastic_mediator_intervention_law = .ltmle_exact_init_mc_history(base, node_spec),
    first_mediator_stochastic_intervention_law = .ltmle_exact_init_mc_history(base, node_spec),
    second_mediator_stochastic_intervention_law = .ltmle_exact_init_mc_history(base, node_spec)
  )

  rows <- list(
    outcome_process = vector("list", T),
    joint_stochastic_mediator_intervention_law = vector("list", T),
    first_mediator_stochastic_intervention_law = vector("list", T),
    second_mediator_stochastic_intervention_law = vector("list", T)
  )

  for (tt in seq_len(T)) {
    sp <- spec[spec$t == tt, , drop = FALSE]

    joint_row <- .ltmle_exact_generate_one_law_history(
      history = histories$joint_stochastic_mediator_intervention_law,
      tt = tt,
      A_value_M1 = sp$m1_A[1L],
      A_value_M2 = sp$m2_A[1L],
      models = models,
      law_type = "joint_stochastic_mediator_intervention_law",
      node_spec = node_spec
    )
    histories$joint_stochastic_mediator_intervention_law <- joint_row$next_history
    joint_row$next_history <- NULL
    rows$joint_stochastic_mediator_intervention_law[[tt]] <- joint_row

    first_row <- .ltmle_exact_generate_one_law_history(
      history = histories$first_mediator_stochastic_intervention_law,
      tt = tt,
      A_value_M1 = sp$m1_A[1L],
      A_value_M2 = sp$m1_A[1L],
      models = models,
      law_type = "first_mediator_stochastic_intervention_law",
      node_spec = node_spec
    )
    histories$first_mediator_stochastic_intervention_law <- first_row$next_history
    first_row$next_history <- NULL
    rows$first_mediator_stochastic_intervention_law[[tt]] <- first_row

    second_row <- .ltmle_exact_generate_one_law_history(
      history = histories$second_mediator_stochastic_intervention_law,
      tt = tt,
      A_value_M1 = sp$m2_A[1L],
      A_value_M2 = sp$m2_A[1L],
      models = models,
      law_type = "second_mediator_stochastic_intervention_law",
      M1_override = NULL,
      node_spec = node_spec
    )
    histories$second_mediator_stochastic_intervention_law <- second_row$next_history
    second_row$next_history <- NULL
    rows$second_mediator_stochastic_intervention_law[[tt]] <- second_row

    outcome_row <- .ltmle_exact_assemble_outcome_process_row(
      base_history = histories$outcome_process,
      tt = tt,
      spec_row = sp,
      joint_row = rows$joint_stochastic_mediator_intervention_law[[tt]]$row,
      first_row = rows$first_mediator_stochastic_intervention_law[[tt]]$row,
      second_row = rows$second_mediator_stochastic_intervention_law[[tt]]$row,
      models = models,
      node_spec = node_spec
    )
    histories$outcome_process <- outcome_row$next_history
    outcome_row$next_history <- NULL
    rows$outcome_process[[tt]] <- outcome_row
  }

  rows
}

.ltmle_exact_history_name_for_process <- function(process_type) {
  switch(process_type,
    outcome_process = "outcome_process",
    post_mediator_covariate_transition = "outcome_process",
    virtual_mixed_continuation_task = "outcome_process",
    observed_mediator_process = "outcome_process",
    joint_stochastic_mediator_intervention_law = "joint_stochastic_mediator_intervention_law",
    first_mediator_stochastic_intervention_law = "first_mediator_stochastic_intervention_law",
    second_mediator_stochastic_intervention_law = "second_mediator_stochastic_intervention_law",
    .stop("Unknown process_type for generated history: ", process_type)
  )
}

.ltmle_exact_generated_row_to_long_block <- function(x, tt, node_spec = NULL) {
  out <- data.frame(
    id = seq_len(nrow(x)),
    t = tt,
    W1 = x$W1,
    W2 = x$W2,
    Y0 = x$Y0,
    M1_0 = x$M1_0,
    M2_0 = x$M2_0,
    A = x$A,
    M1 = x$M1,
    M2 = x$M2,
    Y = if ("Y" %in% names(x)) x$Y else NA_real_,
    M1_lag = x$M1_lag,
    M2_lag = x$M2_lag,
    stringsAsFactors = FALSE
  )
  if ("Y_lag" %in% names(x)) out$Y_lag <- x$Y_lag
  for (L_node in .ltmle_exact_L_nodes(node_spec)) {
    out[[L_node]] <- x[[L_node]]
    lag_nm <- .ltmle_exact_lag_name(L_node)
    if (lag_nm %in% names(x)) out[[lag_nm]] <- x[[lag_nm]]
  }
  if ("L" %in% .ltmle_exact_L_nodes(node_spec) && "L" %in% names(out)) {
    out$L_lag <- x$L_lag %||% out$L_lag
  }
  extra_baseline <- setdiff(names(x), c(names(out), "id0", "AM1", "AM2", "MM12",
                                        "M1_sq", "M2_sq", "L_sq"))
  for (nm in extra_baseline) out[[nm]] <- x[[nm]]
  out
}

.ltmle_exact_generated_history_to_long <- function(generated_histories, process_type, spec,
                                                  node_spec = NULL) {
  history_name <- .ltmle_exact_history_name_for_process(process_type)
  T <- length(generated_histories[[history_name]])
  rows <- vector("list", T)
  for (tt in seq_len(T)) {
    x <- generated_histories[[history_name]][[tt]]$row
    rows[[tt]] <- .ltmle_exact_generated_row_to_long_block(x, tt, node_spec)
  }
  .ltmle_exact_add_terms(do.call(rbind, rows))
}

.ltmle_exact_fit_continuation <- function(y, xdf, covs, learner, sl_library, component,
                                          weights = NULL) {
  missing <- setdiff(covs, names(xdf))
  if (length(missing)) {
    .stop("Continuation fit data missing columns for ", component, ": ", paste(missing, collapse = ", "))
  }
  dat <- .ltmle_exact_model_data(xdf, covs)
  dat$.Y <- as.numeric(y)
  .ltmle_exact_fit_nuis(
    .ltmle_exact_formula(".Y", names(dat)[names(dat) != ".Y"]),
    dat,
    stats::gaussian(),
    learner,
    sl_library,
    component,
    weights = weights,
    fold_id = if (".fold_id" %in% names(xdf)) xdf$.fold_id else NULL
  )
}

.ltmle_exact_predict_continuation <- function(fit, newdata) {
  .ltmle_exact_predict_nuis(fit, newdata, "numeric")
}

.ltmle_exact_design_columns <- function(fit, newdata = NULL) {
  rhs <- as.character(fit$rhs_names %||% character(0))
  if (identical(as.character(fit$type %||% NA_character_), "glm") &&
      !is.null(fit$fit)) {
    cols <- tryCatch(
      {
        if (is.null(newdata)) {
          colnames(stats::model.matrix(fit$fit))
        } else {
          colnames(stats::model.matrix(stats::delete.response(stats::terms(fit$fit)), newdata))
        }
      },
      error = function(e) character(0)
    )
    cols <- setdiff(as.character(cols), "(Intercept)")
    if (length(cols)) return(unique(cols))
  }
  unique(rhs)
}

.ltmle_exact_prediction_model_frame <- function(fit, newdata) {
  required <- as.character(fit$rhs_names %||% character(0))
  pred_data <- .ltmle_exact_prediction_data(newdata, required, keep_fold_id = FALSE)
  if (identical(as.character(fit$type %||% NA_character_), "glm") &&
      !is.null(fit$fit)) {
    out <- tryCatch(
      stats::model.frame(
        stats::delete.response(stats::terms(fit$fit)),
        data = pred_data,
        na.action = stats::na.pass,
        xlev = fit$fit$xlevels %||% NULL
      ),
      error = function(e) NULL
    )
    if (is.data.frame(out) && nrow(out) == nrow(pred_data)) return(out)
  }
  pred_data
}

.ltmle_exact_prediction_design_matrix <- function(fit, newdata) {
  mf <- .ltmle_exact_prediction_model_frame(fit, newdata)
  if (identical(as.character(fit$type %||% NA_character_), "glm") &&
      !is.null(fit$fit)) {
    out <- tryCatch(
      stats::model.matrix(stats::delete.response(stats::terms(fit$fit)), data = mf),
      error = function(e) NULL
    )
    if (is.matrix(out) && nrow(out) == nrow(mf)) return(out)
  }
  if (!ncol(mf)) {
    return(matrix(numeric(0), nrow = nrow(mf), ncol = 0L))
  }
  out <- data.matrix(mf)
  colnames(out) <- names(mf)
  out
}

.ltmle_exact_prediction_linear_predictor <- function(fit, newdata) {
  X <- .ltmle_exact_prediction_design_matrix(fit, newdata)
  pred <- .ltmle_exact_predict_continuation(fit, newdata)
  if (!nrow(X)) return(pred)
  beta <- if (identical(as.character(fit$type %||% NA_character_), "glm") &&
              !is.null(fit$fit)) {
    stats::coef(fit$fit)
  } else {
    numeric(0)
  }
  if (!length(beta)) return(pred)
  common <- intersect(colnames(X), names(beta))
  if (!length(common)) return(pred)
  beta <- as.numeric(beta[common])
  beta[!is.finite(beta)] <- 0
  as.numeric(X[, common, drop = FALSE] %*% beta)
}

.ltmle_exact_serialize_named_numeric_row <- function(x) {
  if (!length(x)) return("")
  paste(names(x), signif(as.numeric(x), 12), sep = "=", collapse = "|")
}

.ltmle_exact_serialize_named_numeric_matrix <- function(x) {
  if (!is.matrix(x) && !is.data.frame(x)) return(character(0))
  if (!nrow(x)) return(character(0))
  mat <- as.matrix(x)
  apply(mat, 1L, .ltmle_exact_serialize_named_numeric_row)
}

.ltmle_exact_numeric_vector_hash <- function(x) {
  x <- as.numeric(x)
  ok <- is.finite(x)
  xx <- x[ok]
  if (!length(xx)) {
    return("n=0|finite_n=0|mean=NA|checksum=NA")
  }
  idx <- seq_along(xx)
  paste(
    "n", length(x),
    "finite_n", length(xx),
    "mean", signif(mean(xx), 14),
    "sd", signif(stats::sd(xx), 14),
    "min", signif(min(xx), 14),
    "max", signif(max(xx), 14),
    "checksum", signif(sum(xx * idx), 14),
    "checksum2", signif(sum((xx^2) * idx), 14),
    sep = "="
  )
}

.ltmle_exact_diag_col <- function(x, nm, numeric = TRUE) {
  if (!is.data.frame(x) || !(nm %in% names(x))) {
    if (isTRUE(numeric)) return(rep(NA_real_, nrow(x)))
    return(rep(NA_character_, nrow(x)))
  }
  if (isTRUE(numeric)) suppressWarnings(as.numeric(x[[nm]])) else as.character(x[[nm]])
}

.ltmle_exact_dedicated_L_transition_instrumentation_enabled <- function(virtual_task,
                                                                        source_task) {
  virtual_task <- as.list(virtual_task)
  source_task <- as.list(source_task)
  component <- as.character(virtual_task$component %||% NA_character_)
  virtual_task_id <- as.character(virtual_task$task_id %||% NA_character_)
  component %in% c("mu_joint_aas", "mu_sep_a_asas") &&
    grepl("::1::virtual_mixed_continuation_after_Y::to_outer_L_2$", virtual_task_id) &&
    identical(as.integer(source_task$t %||% NA_integer_), 2L) &&
    identical(as.character(source_task$node %||% NA_character_), "L") &&
    identical(
      as.character(.ltmle_exact_process_type(source_task)),
      "post_mediator_covariate_transition"
    )
}

.ltmle_exact_dedicated_L_training_row_trace_enabled <- function(task) {
  task <- as.list(task)
  component <- as.character(task$component %||% NA_character_)
  component %in% c("mu_joint_aas", "mu_sep_a_asas") &&
    identical(as.integer(task$t %||% NA_integer_), 2L) &&
    identical(as.character(task$node %||% NA_character_), "L") &&
    identical(
      as.character(.ltmle_exact_process_type(task)),
      "post_mediator_covariate_transition"
    )
}

.ltmle_exact_dedicated_L_fit_training_row_trace <- function(task,
                                                            cont_fit,
                                                            rows,
                                                            response,
                                                            q0_train) {
  task <- as.list(task)
  component <- as.character(task$component %||% NA_character_)
  mf <- .ltmle_exact_prediction_model_frame(cont_fit, rows)
  X <- .ltmle_exact_prediction_design_matrix(cont_fit, rows)
  lp <- .ltmle_exact_prediction_linear_predictor(cont_fit, rows)
  raw_pred <- .ltmle_exact_predict_continuation(cont_fit, rows)
  n <- nrow(rows)
  design_cols <- colnames(X) %||% character(0)
  model_cols <- names(mf) %||% character(0)
  design_values <- .ltmle_exact_serialize_named_numeric_matrix(X)
  model_values <- .ltmle_exact_serialize_named_numeric_matrix(mf)
  design_m1_cols <- grep("(^|[^A-Za-z0-9_])M1($|[^A-Za-z0-9_])|^M1$", design_cols, value = TRUE)
  design_m2_cols <- grep("(^|[^A-Za-z0-9_])M2($|[^A-Za-z0-9_])|^M2$", design_cols, value = TRUE)
  design_value_for <- function(cols) {
    if (!length(cols) || !nrow(X)) return(rep(NA_character_, n))
    apply(X[, cols, drop = FALSE], 1L, .ltmle_exact_serialize_named_numeric_row)
  }
  response <- as.numeric(response)
  q0_train <- as.numeric(q0_train)
  residual <- response - q0_train
  passed <- length(q0_train) == n &&
    length(response) == n &&
    length(lp) == n &&
    all(is.finite(q0_train)) &&
    all(is.finite(response)) &&
    all(is.finite(lp))
  data.frame(
    component = component,
    virtual_task_id = paste0(component, "::1::virtual_mixed_continuation_after_Y::to_outer_L_2"),
    source_fit_task_id = as.character(task$task_id %||% NA_character_),
    task_id = as.character(task$task_id %||% NA_character_),
    t = as.integer(task$t %||% NA_integer_),
    node = as.character(task$node %||% NA_character_),
    process_type = .ltmle_exact_process_type(task),
    row_id = seq_len(n),
    subject_id = as.integer(.ltmle_exact_diag_col(rows, "id0")),
    raw_M1 = .ltmle_exact_diag_col(rows, "M1"),
    raw_M2 = .ltmle_exact_diag_col(rows, "M2"),
    raw_M1_lag = .ltmle_exact_diag_col(rows, "M1_lag"),
    raw_M2_lag = .ltmle_exact_diag_col(rows, "M2_lag"),
    raw_A = .ltmle_exact_diag_col(rows, "A"),
    raw_L_lag = .ltmle_exact_diag_col(rows, "L_lag"),
    raw_Y_lag = .ltmle_exact_diag_col(rows, "Y_lag"),
    model_frame_M1 = .ltmle_exact_diag_col(mf, "M1"),
    model_frame_M2 = .ltmle_exact_diag_col(mf, "M2"),
    model_frame_M1_lag = .ltmle_exact_diag_col(mf, "M1_lag"),
    model_frame_M2_lag = .ltmle_exact_diag_col(mf, "M2_lag"),
    model_frame_A = .ltmle_exact_diag_col(mf, "A"),
    model_frame_L_lag = .ltmle_exact_diag_col(mf, "L_lag"),
    model_frame_Y_lag = .ltmle_exact_diag_col(mf, "Y_lag"),
    design_M1_values = design_value_for(design_m1_cols),
    design_M2_values = design_value_for(design_m2_cols),
    design_M1_columns_present = length(design_m1_cols) > 0L,
    design_M2_columns_present = length(design_m2_cols) > 0L,
    model_frame_columns = paste(model_cols, collapse = "|"),
    design_matrix_columns = paste(design_cols, collapse = "|"),
    model_frame_values = model_values,
    design_matrix_values = design_values,
    training_response = response,
    training_prediction = q0_train,
    raw_training_prediction = as.numeric(raw_pred),
    training_residual = residual,
    linear_predictor = as.numeric(lp),
    passed = passed,
    failure_class = if (passed) {
      "no_failure"
    } else {
      "L_fit_training_row_trace_nonfinite"
    },
    stringsAsFactors = FALSE
  )
}

.ltmle_exact_dedicated_L_prediction_row_trace <- function(virtual_task,
                                                          source_task,
                                                          source_fit,
                                                          rows,
                                                          q_initial,
                                                          metadata) {
  virtual_task <- as.list(virtual_task)
  source_task <- as.list(source_task)
  fit <- source_fit$cont_fit
  mf <- .ltmle_exact_prediction_model_frame(fit, rows)
  X <- .ltmle_exact_prediction_design_matrix(fit, rows)
  lp <- .ltmle_exact_prediction_linear_predictor(fit, rows)
  raw_pred <- .ltmle_exact_predict_continuation(fit, rows)
  n <- nrow(rows)
  design_cols <- colnames(X) %||% character(0)
  model_cols <- names(mf) %||% character(0)
  design_values <- .ltmle_exact_serialize_named_numeric_matrix(X)
  model_values <- .ltmle_exact_serialize_named_numeric_matrix(mf)
  design_m1_cols <- grep("(^|[^A-Za-z0-9_])M1($|[^A-Za-z0-9_])|^M1$", design_cols, value = TRUE)
  design_m2_cols <- grep("(^|[^A-Za-z0-9_])M2($|[^A-Za-z0-9_])|^M2$", design_cols, value = TRUE)
  design_value_for <- function(cols) {
    if (!length(cols) || !nrow(X)) return(rep(NA_character_, n))
    apply(X[, cols, drop = FALSE], 1L, .ltmle_exact_serialize_named_numeric_row)
  }
  passed <- length(design_m1_cols) > 0L &&
    length(design_m2_cols) > 0L &&
    length(lp) == n &&
    length(q_initial) == n &&
    all(is.finite(lp)) &&
    all(is.finite(q_initial))
  data.frame(
    component = as.character(virtual_task$component %||% NA_character_),
    virtual_task_id = as.character(virtual_task$task_id %||% NA_character_),
    task_id = as.character(virtual_task$task_id %||% NA_character_),
    source_task_id = as.character(source_task$task_id %||% NA_character_),
    t = as.integer(source_task$t %||% NA_integer_),
    node = as.character(source_task$node %||% NA_character_),
    process_type = .ltmle_exact_process_type(source_task),
    row_id = seq_len(n),
    subject_id = as.integer(.ltmle_exact_diag_col(rows, "id0")),
    event_id = if ("event_id" %in% names(rows)) {
      as.character(rows$event_id)
    } else {
      as.character(seq_len(n))
    },
    parent_particle_id = if ("parent_particle_id" %in% names(rows)) {
      as.character(rows$parent_particle_id)
    } else if (".parent_particle_id" %in% names(rows)) {
      as.character(rows$.parent_particle_id)
    } else {
      NA_character_
    },
    outer_regimen = as.character(virtual_task$source_boundary_outer_regimen %||%
                                   virtual_task$outer_regimen %||% NA_character_),
    m1_regimen = as.character(virtual_task$source_boundary_m1_regimen %||%
                                virtual_task$m1_regimen %||% NA_character_),
    m2_regimen = as.character(virtual_task$source_boundary_m2_regimen %||%
                                virtual_task$m2_regimen %||% NA_character_),
    outcome_history_state = as.character(metadata$source_boundary_outcome_history_state %||% NA_character_),
    m1_history_state = as.character(metadata$source_boundary_m1_history_state %||% NA_character_),
    m2_history_state = as.character(metadata$source_boundary_m2_history_state %||% NA_character_),
    auxiliary_mediator_history_state =
      as.character(metadata$source_boundary_auxiliary_mediator_history_state %||% NA_character_),
    raw_M1 = .ltmle_exact_diag_col(rows, "M1"),
    raw_M2 = .ltmle_exact_diag_col(rows, "M2"),
    raw_M1_lag = .ltmle_exact_diag_col(rows, "M1_lag"),
    raw_M2_lag = .ltmle_exact_diag_col(rows, "M2_lag"),
    raw_A = .ltmle_exact_diag_col(rows, "A"),
    raw_L_lag = .ltmle_exact_diag_col(rows, "L_lag"),
    raw_Y_lag = .ltmle_exact_diag_col(rows, "Y_lag"),
    model_frame_M1 = .ltmle_exact_diag_col(mf, "M1"),
    model_frame_M2 = .ltmle_exact_diag_col(mf, "M2"),
    model_frame_M1_lag = .ltmle_exact_diag_col(mf, "M1_lag"),
    model_frame_M2_lag = .ltmle_exact_diag_col(mf, "M2_lag"),
    model_frame_A = .ltmle_exact_diag_col(mf, "A"),
    model_frame_L_lag = .ltmle_exact_diag_col(mf, "L_lag"),
    model_frame_Y_lag = .ltmle_exact_diag_col(mf, "Y_lag"),
    design_M1_values = design_value_for(design_m1_cols),
    design_M2_values = design_value_for(design_m2_cols),
    design_M1_columns_present = length(design_m1_cols) > 0L,
    design_M2_columns_present = length(design_m2_cols) > 0L,
    model_frame_columns = paste(model_cols, collapse = "|"),
    design_matrix_columns = paste(design_cols, collapse = "|"),
    model_frame_values = model_values,
    design_matrix_values = design_values,
    linear_predictor = as.numeric(lp),
    raw_L_fit_prediction = as.numeric(raw_pred),
    predicted_L_value = as.numeric(q_initial),
    integrated_node_oracle_mean = NA_real_,
    prediction_minus_integrated_node_oracle = NA_real_,
    passed = passed,
    failure_class = if (passed) {
      "no_failure"
    } else if (!length(design_m1_cols) || !length(design_m2_cols)) {
      "L_fit_design_matrix_missing_current_mediator_columns"
    } else {
      "L_fit_prediction_row_trace_nonfinite"
    },
    stringsAsFactors = FALSE
  )
}

.ltmle_exact_dedicated_L_coefficient_contribution_check <- function(virtual_task,
                                                                    source_task,
                                                                    source_fit,
                                                                    rows) {
  virtual_task <- as.list(virtual_task)
  source_task <- as.list(source_task)
  fit <- source_fit$cont_fit
  if (!identical(as.character(fit$type %||% NA_character_), "glm") ||
      is.null(fit$fit)) {
    return(data.frame())
  }
  X <- .ltmle_exact_prediction_design_matrix(fit, rows)
  beta <- stats::coef(fit$fit)
  terms <- intersect(colnames(X), names(beta))
  if (!length(terms)) return(data.frame())
  row_for_term <- lapply(terms, function(term_name) {
    contribution <- as.numeric(X[, term_name]) * as.numeric(beta[[term_name]])
    is_m1 <- grepl("(^|[^A-Za-z0-9_])M1($|[^A-Za-z0-9_])|^M1$", term_name)
    is_m2 <- grepl("(^|[^A-Za-z0-9_])M2($|[^A-Za-z0-9_])|^M2$", term_name)
    data.frame(
      row_type = "term",
      component = as.character(virtual_task$component %||% NA_character_),
      virtual_task_id = as.character(virtual_task$task_id %||% NA_character_),
      L_fit_task_id = as.character(source_task$task_id %||% NA_character_),
      term_name = term_name,
      coefficient = as.numeric(beta[[term_name]]),
      mean_design_value = mean(as.numeric(X[, term_name]), na.rm = TRUE),
      mean_contribution = mean(contribution, na.rm = TRUE),
      sd_contribution = stats::sd(contribution, na.rm = TRUE),
      min_contribution = min(contribution, na.rm = TRUE),
      max_contribution = max(contribution, na.rm = TRUE),
      mean_abs_contribution = mean(abs(contribution), na.rm = TRUE),
      is_M1_term = is_m1,
      is_M2_term = is_m2,
      is_mediator_interaction_term = is_m1 && is_m2,
      coefficient_abs_sum_M1_terms = NA_real_,
      coefficient_abs_sum_M2_terms = NA_real_,
      coefficient_abs_sum_mediator_interaction_terms = NA_real_,
      mean_abs_contribution_M1_terms = NA_real_,
      mean_abs_contribution_M2_terms = NA_real_,
      mean_abs_contribution_mediator_interaction_terms = NA_real_,
      mediator_terms_have_nonzero_contribution = NA,
      passed = TRUE,
      failure_class = "no_failure",
      stringsAsFactors = FALSE
    )
  })
  term_rows <- .ltmle_exact_rbind_fill(row_for_term)
  m1_idx <- as.logical(term_rows$is_M1_term)
  m2_idx <- as.logical(term_rows$is_M2_term)
  int_idx <- as.logical(term_rows$is_mediator_interaction_term)
  mediator_nonzero <- any((m1_idx | m2_idx) &
                            abs(term_rows$coefficient) > 1e-12 &
                            term_rows$mean_abs_contribution > 1e-12,
                          na.rm = TRUE)
  summary_row <- data.frame(
    row_type = "summary",
    component = as.character(virtual_task$component %||% NA_character_),
    virtual_task_id = as.character(virtual_task$task_id %||% NA_character_),
    L_fit_task_id = as.character(source_task$task_id %||% NA_character_),
    term_name = "__summary__",
    coefficient = NA_real_,
    mean_design_value = NA_real_,
    mean_contribution = NA_real_,
    sd_contribution = NA_real_,
    min_contribution = NA_real_,
    max_contribution = NA_real_,
    mean_abs_contribution = NA_real_,
    is_M1_term = FALSE,
    is_M2_term = FALSE,
    is_mediator_interaction_term = FALSE,
    coefficient_abs_sum_M1_terms = sum(abs(term_rows$coefficient[m1_idx]), na.rm = TRUE),
    coefficient_abs_sum_M2_terms = sum(abs(term_rows$coefficient[m2_idx]), na.rm = TRUE),
    coefficient_abs_sum_mediator_interaction_terms =
      sum(abs(term_rows$coefficient[int_idx]), na.rm = TRUE),
    mean_abs_contribution_M1_terms = sum(term_rows$mean_abs_contribution[m1_idx], na.rm = TRUE),
    mean_abs_contribution_M2_terms = sum(term_rows$mean_abs_contribution[m2_idx], na.rm = TRUE),
    mean_abs_contribution_mediator_interaction_terms =
      sum(term_rows$mean_abs_contribution[int_idx], na.rm = TRUE),
    mediator_terms_have_nonzero_contribution = mediator_nonzero,
    passed = mediator_nonzero,
    failure_class = if (mediator_nonzero) {
      "no_failure"
    } else {
      "L_fit_mediator_terms_have_zero_prediction_contribution"
    },
    stringsAsFactors = FALSE
  )
  .ltmle_exact_rbind_fill(list(term_rows, summary_row))
}

.ltmle_exact_dedicated_L_prediction_weighting_collapse_check <- function(virtual_task,
                                                                         source_task,
                                                                         source_boundary,
                                                                         q_initial,
                                                                         subject_initial) {
  virtual_task <- as.list(virtual_task)
  source_task <- as.list(source_task)
  particles <- source_boundary$branch_state[[source_boundary$state_key]]$particles
  weights <- as.numeric(particles$.branch_weight %||% rep(1, nrow(particles)))
  weights[!is.finite(weights) | weights < 0] <- 0
  unweighted <- mean(as.numeric(q_initial), na.rm = TRUE)
  weighted <- if (sum(weights) > 0) {
    sum(as.numeric(q_initial) * weights, na.rm = TRUE) / sum(weights, na.rm = TRUE)
  } else {
    NA_real_
  }
  final <- mean(as.numeric(subject_initial), na.rm = TRUE)
  passed <- is.finite(unweighted) && is.finite(weighted) && is.finite(final) &&
    is.finite(sum(weights)) && sum(weights) > 0
  data.frame(
    component = as.character(virtual_task$component %||% NA_character_),
    virtual_task_id = as.character(virtual_task$task_id %||% NA_character_),
    task_id = as.character(virtual_task$task_id %||% NA_character_),
    source_task_id = as.character(source_task$task_id %||% NA_character_),
    t = as.integer(source_task$t %||% NA_integer_),
    node = as.character(source_task$node %||% NA_character_),
    process_type = .ltmle_exact_process_type(source_task),
    unweighted_prediction_mean = unweighted,
    mediator_weighted_prediction_mean = weighted,
    final_collapse_mean = final,
    integrated_node_oracle_mean = NA_real_,
    mediator_weights_sum = sum(weights, na.rm = TRUE),
    mediator_weights_min = min(weights, na.rm = TRUE),
    mediator_weights_max = max(weights, na.rm = TRUE),
    mediator_weights_ess = .ltmle_exact_ess(weights),
    event_weights_sum = sum(weights, na.rm = TRUE),
    event_weights_min = min(weights, na.rm = TRUE),
    event_weights_max = max(weights, na.rm = TRUE),
    event_weights_ess = .ltmle_exact_ess(weights),
    prediction_minus_oracle = NA_real_,
    weighted_prediction_minus_oracle = NA_real_,
    final_collapse_minus_oracle = NA_real_,
    weighting_changes_prediction = is.finite(unweighted) &&
      is.finite(weighted) &&
      abs(weighted - unweighted) > 1e-10,
    collapse_changes_weighted_prediction = is.finite(weighted) &&
      is.finite(final) &&
      abs(final - weighted) > 1e-10,
    passed = passed,
    failure_class = if (passed) {
      "no_failure"
    } else {
      "L_transition_prediction_weighting_or_collapse_nonfinite"
    },
    stringsAsFactors = FALSE
  )
}

.ltmle_exact_dedicated_L_mediator_path_weight_trace <- function(virtual_task,
                                                                source_task,
                                                                source_boundary,
                                                                q_initial) {
  virtual_task <- as.list(virtual_task)
  source_task <- as.list(source_task)
  particles <- source_boundary$branch_state[[source_boundary$state_key]]$particles
  n <- nrow(particles)
  q_initial <- as.numeric(q_initial)
  weight_col <- function(nm, default = 1) {
    x <- if (nm %in% names(particles)) {
      as.numeric(particles[[nm]])
    } else {
      rep(default, n)
    }
    x[!is.finite(x)] <- default
    x
  }
  value_col <- function(...) {
    candidates <- c(...)
    hit <- candidates[candidates %in% names(particles)]
    if (!length(hit)) return(rep(NA_real_, n))
    suppressWarnings(as.numeric(particles[[hit[1L]]]))
  }
  char_col <- function(...) {
    candidates <- c(...)
    hit <- candidates[candidates %in% names(particles)]
    if (!length(hit)) return(rep(NA_character_, n))
    as.character(particles[[hit[1L]]])
  }
  branch_weight <- weight_col(".branch_weight")
  outcome_weight <- weight_col(".outcome_weight")
  joint_law_weight <- weight_col(".joint_law_weight")
  first_law_weight <- weight_col(".first_law_weight")
  second_law_weight <- weight_col(".second_law_weight")
  product_join_weight <- weight_col(".product_join_weight")
  integration_weight <- weight_col(".integration_weight")
  component <- as.character(virtual_task$component %||% NA_character_)
  metadata <- source_boundary$provenance %||% source_boundary$metadata %||% list()
  m1_state <- as.character(
    metadata$source_boundary_m1_history_state %||%
      virtual_task$source_boundary_m1_history_state %||% NA_character_
  )
  m2_state <- as.character(
    metadata$source_boundary_m2_history_state %||%
      virtual_task$source_boundary_m2_history_state %||% NA_character_
  )
  expected_weight <- if (identical(m1_state, "joint_law") &&
                         identical(m2_state, "joint_law")) {
    joint_law_weight
  } else if (identical(m1_state, "first_law") &&
             identical(m2_state, "second_law")) {
    product_join_weight
  } else if (identical(m1_state, "first_law") &&
             identical(m2_state, "first_law")) {
    first_law_weight
  } else if (identical(m1_state, "second_law") &&
             identical(m2_state, "second_law")) {
    second_law_weight
  } else if (identical(m1_state, "outcome") &&
             identical(m2_state, "outcome")) {
    outcome_weight
  } else {
    branch_weight
  }
  actual_weight <- branch_weight
  expected_sum_by_id <- ave(expected_weight, particles$id0, FUN = sum)
  actual_sum_by_id <- ave(actual_weight, particles$id0, FUN = sum)
  expected_norm <- expected_weight / pmax(expected_sum_by_id, 1e-12)
  actual_norm <- actual_weight / pmax(actual_sum_by_id, 1e-12)
  weight_matches <- is.finite(expected_norm) & is.finite(actual_norm) &
    abs(expected_norm - actual_norm) <= 1e-10
  uniform_by_id <- ave(actual_norm, particles$id0, FUN = function(x) {
    max(abs(x - mean(x)), na.rm = TRUE)
  }) <= 1e-10
  passed <- length(q_initial) == n &&
    all(is.finite(q_initial)) &&
    all(is.finite(actual_weight)) &&
    all(is.finite(expected_weight)) &&
    all(weight_matches)
  data.frame(
    component = component,
    virtual_task_id = as.character(virtual_task$task_id %||% NA_character_),
    source_fit_task_id = as.character(source_task$task_id %||% NA_character_),
    task_id = as.character(virtual_task$task_id %||% NA_character_),
    source_task_id = as.character(source_task$task_id %||% NA_character_),
    row_id = seq_len(n),
    subject_id = as.integer(value_col("id0")),
    event_id = char_col("event_id", ".event_id"),
    parent_particle_id = char_col("parent_particle_id", ".parent_particle_id"),
    M1_value = value_col("M1"),
    M2_value = value_col("M2"),
    M1_grid_z = value_col("M1_grid_z", ".M1_grid_z", "M1_z", ".M1_z"),
    M2_grid_z = value_col("M2_grid_z", ".M2_grid_z", "M2_z", ".M2_z"),
    M1_grid_weight = value_col("M1_grid_weight", ".M1_grid_weight", "M1_weight", ".M1_weight"),
    M2_grid_weight = value_col("M2_grid_weight", ".M2_grid_weight", "M2_weight", ".M2_weight"),
    branch_weight = branch_weight,
    outcome_weight = outcome_weight,
    joint_law_weight = joint_law_weight,
    first_law_weight = first_law_weight,
    second_law_weight = second_law_weight,
    product_join_weight = product_join_weight,
    integration_weight = integration_weight,
    expected_component_law_weight = expected_weight,
    actual_aggregation_weight = actual_weight,
    expected_component_law_weight_normalized_by_subject = expected_norm,
    actual_aggregation_weight_normalized_by_subject = actual_norm,
    row_prediction = q_initial,
    weight_matches_expected_component_law = weight_matches,
    weight_is_uniform = uniform_by_id,
    weight_component_specific = !uniform_by_id,
    passed = passed,
    failure_class = if (passed) {
      "no_failure"
    } else if (length(q_initial) != n || any(!is.finite(q_initial))) {
      "mediator_path_weight_trace_missing"
    } else {
      "mediator_path_weight_not_used_in_aggregation"
    },
    stringsAsFactors = FALSE
  )
}

.ltmle_exact_apply_fluctuation <- function(q0_new, H_new, epsilon, bounds, eps = 1e-8) {
  q01_new <- .ltmle_exact_scale01(q0_new, bounds, eps)
  H_new <- as.numeric(H_new)

  qstar_new01 <- .ltmle_exact_clamp01(
    stats::plogis(stats::qlogis(q01_new) + epsilon * H_new),
    eps
  )

  .ltmle_exact_unscale01(qstar_new01, bounds)
}

.ltmle_exact_ess <- function(w) {
  w <- as.numeric(w)
  w <- w[is.finite(w) & w >= 0]
  if (!length(w) || sum(w) <= 0) return(NA_real_)
  sum(w)^2 / sum(w^2)
}

.ltmle_exact_rbind_fill <- function(rows) {
  rows <- rows[!vapply(rows, is.null, logical(1))]
  rows <- rows[vapply(rows, function(x) is.data.frame(x) && nrow(x) > 0L, logical(1))]
  if (!length(rows)) return(data.frame())
  cols <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(x) {
    miss <- setdiff(cols, names(x))
    for (cc in miss) x[[cc]] <- NA
    x[, cols, drop = FALSE]
  })
  do.call(rbind, rows)
}

.ltmle_exact_target_continuation <- function(y_train, q0_train,
                                             q0_new = NULL,
                                             H_obs,
                                             H_new = NULL,
                                             y_bounds_mode, y_bounds, learner, node,
                                             component, tt, process_type,
                                             score_tolerance, included_ratio_factors,
                                             task_id = NA_character_,
                                             source_task_id = NA_character_,
                                             parent_task_id = NA_character_,
                                             role = NA_character_,
                                             eps = 1e-8,
                                             force_epsilon_zero = FALSE,
                                             compute_new = !is.null(q0_new) && !is.null(H_new)) {
  bnd <- .ltmle_exact_bounds(y_train, y_bounds_mode, y_bounds)

  y01 <- .ltmle_exact_scale01(y_train, bnd, eps)
  q01 <- .ltmle_exact_scale01(q0_train, bnd, eps)
  density_diag <- attr(H_obs, "density_ratio_diagnostics")
  H_obs <- as.numeric(H_obs)

  keep <- is.finite(y01) & is.finite(q01) & is.finite(H_obs) & abs(H_obs) > 0
  if (!any(keep)) {
    .stop("No usable observations for exact targeting: ", component, " node=", node, " t=", tt)
  }

  dat <- data.frame(
    y = y01[keep],
    H = H_obs[keep],
    off = stats::qlogis(q01[keep])
  )

  raw_score_tolerance <- as.numeric(score_tolerance)[1L]
  scaled_score_tolerance <- 1e-3
  standardized_score_tolerance <- 0.1
  score_for_epsilon <- function(epsilon_value, clamp = FALSE) {
    q01_eval <- stats::plogis(dat$off + as.numeric(epsilon_value)[1L] * dat$H)
    if (isTRUE(clamp)) q01_eval <- .ltmle_exact_clamp01(q01_eval, eps)
    q_eval <- .ltmle_exact_unscale01(q01_eval, bnd)
    mean(H_obs[keep] * (as.numeric(y_train)[keep] - q_eval), na.rm = TRUE)
  }
  solve_by_bisection <- function(center, tolerance) {
    center <- as.numeric(center)[1L]
    if (!is.finite(center)) center <- 0
    f_center <- score_for_epsilon(center, clamp = FALSE)
    if (is.finite(f_center) && abs(f_center) <= tolerance) {
      return(list(
        epsilon = center,
        score = f_center,
        converged = TRUE,
        iterations = 0L,
        failure_class = "no_failure"
      ))
    }
    width <- 1
    lower <- center - width
    upper <- center + width
    f_lower <- score_for_epsilon(lower, clamp = FALSE)
    f_upper <- score_for_epsilon(upper, clamp = FALSE)
    iter <- 0L
    while ((is.na(f_lower) || is.na(f_upper) || f_lower * f_upper > 0) &&
           width < 1e6 && iter < 80L) {
      width <- width * 2
      lower <- center - width
      upper <- center + width
      f_lower <- score_for_epsilon(lower, clamp = FALSE)
      f_upper <- score_for_epsilon(upper, clamp = FALSE)
      iter <- iter + 1L
    }
    if (!is.finite(f_lower) || !is.finite(f_upper) || f_lower * f_upper > 0) {
      return(list(
        epsilon = center,
        score = f_center,
        converged = FALSE,
        iterations = iter,
        failure_class = "root_bracket_not_found"
      ))
    }
    lo <- lower
    hi <- upper
    f_lo <- f_lower
    mid <- center
    f_mid <- f_center
    for (jj in seq_len(200L)) {
      mid <- (lo + hi) / 2
      f_mid <- score_for_epsilon(mid, clamp = FALSE)
      iter <- iter + 1L
      if (is.finite(f_mid) && abs(f_mid) <= tolerance) break
      if (!is.finite(f_mid)) break
      if (f_lo * f_mid <= 0) {
        hi <- mid
      } else {
        lo <- mid
        f_lo <- f_mid
      }
    }
    list(
      epsilon = mid,
      score = f_mid,
      converged = is.finite(f_mid) && abs(f_mid) <= tolerance,
      iterations = iter,
      failure_class = if (is.finite(f_mid) && abs(f_mid) <= tolerance) {
        "no_failure"
      } else {
        "root_bisection_not_converged"
      }
    )
  }

  epsilon_forced_zero <- isTRUE(force_epsilon_zero)
  solver_method <- "glm_quasibinomial"
  solver_converged <- TRUE
  solver_iterations <- NA_integer_
  solver_failure_class <- "no_failure"
  if (epsilon_forced_zero) {
    epsilon <- 0
    solver_method <- "force_epsilon_zero"
    solver_converged <- NA
    solver_iterations <- 0L
    solver_failure_class <- "not_applicable_force_epsilon_zero"
  } else {
    fl <- tryCatch(
      stats::glm(
        y ~ 0 + H,
        family = stats::quasibinomial(),
        offset = off,
        data = dat
      ),
      error = function(e) e
    )
    epsilon_glm <- if (inherits(fl, "error")) {
      NA_real_
    } else {
      as.numeric(stats::coef(fl)[["H"]])
    }
    if (!is.finite(epsilon_glm)) epsilon_glm <- 0
    root <- solve_by_bisection(epsilon_glm, tolerance = min(raw_score_tolerance / 10, 1e-8))
    score_glm <- score_for_epsilon(epsilon_glm, clamp = FALSE)
    if (isTRUE(root$converged) ||
        (is.finite(root$score) && is.finite(score_glm) && abs(root$score) < abs(score_glm))) {
      epsilon <- root$epsilon
      solver_method <- if (inherits(fl, "error")) {
        "bracketed_bisection_after_glm_failure"
      } else {
        "glm_quasibinomial_plus_bracketed_bisection"
      }
      solver_converged <- isTRUE(root$converged)
      solver_iterations <- as.integer(root$iterations)
      solver_failure_class <- root$failure_class
    } else {
      epsilon <- epsilon_glm
      solver_method <- if (inherits(fl, "error")) "glm_quasibinomial_failed" else "glm_quasibinomial"
      solver_converged <- !inherits(fl, "error") && is.finite(score_glm) &&
        abs(score_glm) <= raw_score_tolerance
      solver_iterations <- NA_integer_
      solver_failure_class <- if (solver_converged) "no_failure" else "glm_score_not_converged"
    }
  }

  qstar_train01_unclamped <- stats::plogis(stats::qlogis(q01) + epsilon * H_obs)
  qstar_train01 <- .ltmle_exact_clamp01(qstar_train01_unclamped, eps)
  qstar_train <- .ltmle_exact_unscale01(qstar_train01, bnd)

  if (isTRUE(compute_new)) {
    qstar_new <- .ltmle_exact_apply_fluctuation(
      q0_new = q0_new,
      H_new = H_new,
      epsilon = epsilon,
      bounds = bnd,
      eps = eps
    )
    H_new <- as.numeric(H_new)
    H_new_mean <- mean(H_new[is.finite(H_new)], na.rm = TRUE)
    H_new_max <- max(abs(H_new[is.finite(H_new)]), na.rm = TRUE)
  } else {
    qstar_new <- NULL
    H_new_mean <- NA_real_
    H_new_max <- NA_real_
  }

  D_task <- H_obs * (as.numeric(y_train) - qstar_train)
  D_task[!is.finite(D_task)] <- NA_real_

  score_after <- mean(D_task[keep], na.rm = TRUE)
  qstar_train_unclamped <- .ltmle_exact_unscale01(qstar_train01_unclamped, bnd)
  score_after_before_clamp <- mean(
    H_obs[keep] * (as.numeric(y_train)[keep] - qstar_train_unclamped[keep]),
    na.rm = TRUE
  )
  score_after_after_clamp <- score_after
  score_before <- mean(
    H_obs[keep] * (as.numeric(y_train)[keep] - q0_train[keep]),
    na.rm = TRUE
  )
  score_contrib <- D_task[keep]
  score_contribution_mean_abs <- mean(abs(score_contrib), na.rm = TRUE)
  score_contribution_sd <- stats::sd(score_contrib, na.rm = TRUE)
  top_abs_fraction <- function(x, frac) {
    x <- abs(as.numeric(x))
    x <- x[is.finite(x)]
    if (!length(x) || sum(x) <= 0) return(NA_real_)
    k <- max(1L, ceiling(length(x) * frac))
    sum(sort(x, decreasing = TRUE)[seq_len(k)], na.rm = TRUE) / sum(x, na.rm = TRUE)
  }
  scaled_score <- abs(score_after) / max(score_contribution_mean_abs, .Machine$double.eps)
  standardized_score <- abs(score_after) /
    max(score_contribution_sd / sqrt(sum(keep)), .Machine$double.eps)
  passed_raw_score <- is.finite(score_after) && abs(score_after) <= raw_score_tolerance
  passed_scaled_score <- is.finite(scaled_score) && scaled_score <= scaled_score_tolerance
  passed_standardized_score <- is.finite(standardized_score) &&
    standardized_score <= standardized_score_tolerance
  empirical_equation_required <- !epsilon_forced_zero
  passed_empirical_equation <- if (empirical_equation_required) {
    passed_raw_score || (passed_scaled_score && passed_standardized_score)
  } else {
    TRUE
  }
  clamp_delta <- abs(qstar_train01 - qstar_train01_unclamped)
  clamp_applied <- any(clamp_delta > 0, na.rm = TRUE)
  clamp_fraction <- mean(clamp_delta[is.finite(clamp_delta)] > 0, na.rm = TRUE)
  if (!is.finite(clamp_fraction)) clamp_fraction <- 0
  density_ratio_min <- if (!is.null(density_diag) && nrow(density_diag) && "min_ratio" %in% names(density_diag)) {
    min(as.numeric(density_diag$min_ratio), na.rm = TRUE)
  } else {
    min(abs(H_obs[keep]), na.rm = TRUE)
  }
  density_ratio_mean <- if (!is.null(density_diag) && nrow(density_diag) && "mean_ratio" %in% names(density_diag)) {
    mean(as.numeric(density_diag$mean_ratio), na.rm = TRUE)
  } else {
    mean(abs(H_obs[keep]), na.rm = TRUE)
  }
  density_ratio_max <- if (!is.null(density_diag) && nrow(density_diag) && "max_ratio" %in% names(density_diag)) {
    max(as.numeric(density_diag$max_ratio), na.rm = TRUE)
  } else {
    max(abs(H_obs[keep]), na.rm = TRUE)
  }
  row_hash <- function(...) {
    vals <- unlist(list(...), use.names = FALSE)
    vals <- vals[is.finite(vals)]
    paste0(
      "n=", length(vals),
      ";sum=", signif(sum(vals), 12),
      ";mean=", signif(mean(vals), 12)
    )
  }
  training_rows_hash <- row_hash(y01[keep], q01[keep], H_obs[keep])
  score_rows_hash <- row_hash(score_contrib)
  H_keep <- H_obs[keep]
  q0_keep <- q0_train[keep]
  qstar_keep <- qstar_train[keep]
  failure_class <- if (!empirical_equation_required) {
    "not_applicable_force_epsilon_zero"
  } else if (passed_empirical_equation) {
    "no_failure"
  } else if (!isTRUE(solver_converged)) {
    solver_failure_class
  } else {
    "empirical_equation_not_solved"
  }

  diag <- data.frame(
    component = component,
    task_id = task_id,
    source_task_id = source_task_id,
    parent_task_id = parent_task_id,
    t = tt,
    node = node,
    process_type = process_type,
    role = role,
    fold = NA_character_,
    n_rows = sum(keep),
    n_nonzero_H = sum(keep),
    n_train = sum(keep),
    epsilon = epsilon,
    abs_epsilon = abs(epsilon),
    max_abs_epsilon = abs(epsilon),
    epsilon_forced_zero = epsilon_forced_zero,
    force_epsilon_zero = epsilon_forced_zero,
    targeting_disabled = epsilon_forced_zero,
    targeting_applied = !epsilon_forced_zero,
    score_before = score_before,
    score_after = score_after,
    abs_score_after = abs(score_after),
    score_abs_after = abs(score_after),
    raw_score_tolerance = raw_score_tolerance,
    score_tolerance = score_tolerance,
    passed_raw_score = passed_raw_score,
    score_contribution_mean_abs = score_contribution_mean_abs,
    score_contribution_sd = score_contribution_sd,
    top_1pct_abs_score_contribution_fraction = top_abs_fraction(score_contrib, 0.01),
    top_5pct_abs_score_contribution_fraction = top_abs_fraction(score_contrib, 0.05),
    top_10pct_abs_score_contribution_fraction = top_abs_fraction(score_contrib, 0.10),
    scaled_score = scaled_score,
    scaled_score_tolerance = scaled_score_tolerance,
    passed_scaled_score = passed_scaled_score,
    standardized_score = standardized_score,
    standardized_score_tolerance = standardized_score_tolerance,
    passed_standardized_score = passed_standardized_score,
    empirical_equation_required = empirical_equation_required,
    empirical_equation_check_status = if (empirical_equation_required) {
      "evaluated"
    } else {
      "not_applicable_force_epsilon_zero"
    },
    score_equation_solved = passed_empirical_equation,
    solved = passed_empirical_equation,
    solver_method = solver_method,
    solver_converged = solver_converged,
    solver_iterations = solver_iterations,
    solver_failure_class = solver_failure_class,
    H_mean = mean(H_keep, na.rm = TRUE),
    H_sd = stats::sd(H_keep, na.rm = TRUE),
    H_min = min(H_keep, na.rm = TRUE),
    H_max = max(H_keep, na.rm = TRUE),
    H_abs_max = max(abs(H_keep), na.rm = TRUE),
    H_ess = .ltmle_exact_ess(abs(H_obs[keep])),
    Q_initial_mean = mean(q0_keep, na.rm = TRUE),
    Q_targeted_mean = mean(qstar_keep, na.rm = TRUE),
    Q_targeted_minus_initial = mean(qstar_keep - q0_keep, na.rm = TRUE),
    Q_initial_min = min(q0_keep, na.rm = TRUE),
    Q_initial_max = max(q0_keep, na.rm = TRUE),
    Q_targeted_min = min(qstar_keep, na.rm = TRUE),
    Q_targeted_max = max(qstar_keep, na.rm = TRUE),
    clamp_applied = clamp_applied,
    clamp_fraction = clamp_fraction,
    score_after_before_clamp = score_after_before_clamp,
    score_after_after_clamp = score_after_after_clamp,
    probability_bound_lower = bnd$a,
    probability_bound_upper = bnd$b,
    density_ratio_min = density_ratio_min,
    density_ratio_mean = density_ratio_mean,
    density_ratio_max = density_ratio_max,
    training_rows_hash = training_rows_hash,
    score_rows_hash = score_rows_hash,
    id_alignment_checked = TRUE,
    n_rows_targeting = sum(keep),
    n_rows_score_check = sum(keep),
    passed = passed_empirical_equation,
    failure_class = failure_class,
    generated_evaluated = isTRUE(compute_new),
    H_new_mean = H_new_mean,
    H_new_max = H_new_max,
    H_source_mean = NA_real_,
    H_source_max = NA_real_,
    H_obs_source_max_abs_diff = NA_real_,
    source_rows_evaluated = FALSE,
    terminal_plugin_batches = NA_integer_,
    terminal_plugin_type = NA_character_,
    terminal_plugin_mc_active = NA,
    effective_mc_n = NA_integer_,
    included_ratio_factors = included_ratio_factors,
    y_bounds_lower = bnd$a,
    y_bounds_upper = bnd$b,
    stringsAsFactors = FALSE
  )

  list(
    train = qstar_train,
    new = qstar_new,
    diagnostics = diag,
    D_task = D_task,
    epsilon = epsilon,
    bounds = bnd
  )
}

.ltmle_exact_sort_tasks_for_graph <- function(component_tasks) {
  node <- as.character(component_tasks$node)
  node_group <- ifelse(node == "Y", 1L,
                       ifelse(node == "L" | startsWith(node, "L_"), 2L,
                              ifelse(node == "M2", 3L, 4L)))
  role_group <- ifelse(component_tasks$process_type == "second_mediator_stochastic_intervention_law" &
                         component_tasks$node == "M1", 1L, 2L)
  component_tasks[order(-as.integer(component_tasks$t), node_group, role_group), , drop = FALSE]
}

.ltmle_exact_topological_sort_from_sources <- function(component_tasks) {
  ids <- as.character(component_tasks$task_id)
  if (anyDuplicated(ids)) .stop("Duplicate ltmle_exact task_id values in factor task registry.")
  deps <- stats::setNames(vector("list", length(ids)), ids)
  for (ii in seq_along(ids)) {
    src <- as.character(component_tasks$observed_pseudooutcome_source_task_id[ii])
    deps[[ids[ii]]] <- intersect(src, ids)
  }

  done <- character(0)
  remaining <- ids
  while (length(remaining)) {
    ready <- remaining[vapply(remaining, function(id) all(deps[[id]] %in% done), logical(1))]
    if (!length(ready)) {
      blocked <- remaining[seq_len(min(5L, length(remaining)))]
      .stop("ltmle_exact task graph has cyclic or unsatisfied pseudo-outcome dependencies near: ",
            paste(blocked, collapse = ", "))
    }
    ready_rows <- component_tasks[match(ready, component_tasks$task_id), , drop = FALSE]
    ready_rows <- .ltmle_exact_sort_tasks_for_graph(ready_rows)
    ready <- as.character(ready_rows$task_id)
    done <- c(done, ready)
    remaining <- setdiff(remaining, ready)
  }
  done
}

.ltmle_exact_build_task_graph <- function(component_tasks, component_spec, node_spec, T) {
  required <- c(
    "component", "t", "node", "process_type", "task_id",
    "observed_pseudooutcome_source_task_id"
  )
  miss <- setdiff(required, names(component_tasks))
  if (length(miss)) {
    .stop("ltmle_exact factor task registry is missing task dependency columns: ",
          paste(miss, collapse = ", "))
  }
  if (any(!nzchar(as.character(component_tasks$task_id)))) {
    .stop("Every ltmle_exact factor task must have a non-empty task_id.")
  }
  allowed_sources <- c(
    as.character(component_tasks$task_id),
    "observed_terminal_outcome"
  )
  bad_obs <- setdiff(as.character(component_tasks$observed_pseudooutcome_source_task_id), allowed_sources)
  if (length(bad_obs)) {
    .stop("ltmle_exact task graph contains invalid pseudo-outcome source task_id.")
  }

  tasks <- list()
  for (ii in seq_len(nrow(component_tasks))) {
    task <- as.list(component_tasks[ii, , drop = FALSE])
    task$time <- task$t
    task$conditioning_history_columns <- .ltmle_exact_task_history_covs(task, "correct", node_spec)
    task$assigned_treatment_regimen <- task$assigned_treatment_regimen_label %||% NA_character_
    task$parent_task_id <- task$observed_pseudooutcome_source_task_id
    tasks[[task$task_id]] <- task
  }
  reverse_order <- .ltmle_exact_topological_sort_from_sources(component_tasks)
  L_nodes <- .ltmle_exact_L_nodes(node_spec)
  if (!length(L_nodes)) {
    .stop("ltmle_exact deterministic root plug-in requires at least one post-mediator L node.")
  }
  root <- component_tasks[
    component_tasks$t == 1L &
      component_tasks$process_type == "post_mediator_covariate_transition" &
      component_tasks$node == L_nodes[1L],
    ,
    drop = FALSE
  ]
  if (!nrow(root)) {
    .stop(
      "ltmle_exact requires deterministic root task at t=1, node=first L node, ",
      "process_type='post_mediator_covariate_transition'. MC fallback is not supported."
    )
  }
  terminal_task_id <- as.character(root$task_id[1L])
  list(
    tasks = tasks,
    reverse_topological_order = reverse_order,
    observed_values = list(),
    generated_values = list(),
    terminal_task_id = terminal_task_id,
    graph_used = TRUE,
    T = T
  )
}

.ltmle_exact_is_supported_deterministic_terminal_root <- function(task, node_spec = NULL) {
  task <- as.list(task)
  L_nodes <- .ltmle_exact_L_nodes(node_spec)

  length(L_nodes) >= 1L &&
    identical(as.integer(task$t), 1L) &&
    identical(as.character(task$process_type), "post_mediator_covariate_transition") &&
    identical(as.character(task$node), as.character(L_nodes[1L]))
}

.ltmle_exact_assert_supported_deterministic_terminal_root <- function(task, node_spec = NULL) {
  if (!.ltmle_exact_is_supported_deterministic_terminal_root(task, node_spec)) {
    .stop(
      "Unsupported ltmle_exact terminal/root task. ",
      "This implementation supports only deterministic root plug-in ",
      "with t=1, node=first L node, process_type='post_mediator_covariate_transition'. ",
      "Got t=", task$t,
      ", node=", task$node,
      ", process_type=", task$process_type,
      ". MC fallback is intentionally not supported."
    )
  }
  invisible(TRUE)
}

.ltmle_exact_make_deterministic_root_plugin_data <- function(dat_wide, spec, node_spec = NULL) {
  sp <- spec[spec$t == 1L, , drop = FALSE]
  if (!nrow(sp)) .stop("No t=1 spec row for deterministic root plug-in.")

  base <- dat_wide
  base$id0 <- seq_len(nrow(dat_wide))
  history <- .ltmle_exact_init_mc_history(base, node_spec)

  row <- data.frame(
    id0 = history$id0,
    W1 = history$W1,
    W2 = history$W2,
    Y0 = history$Y0,
    M1_0 = history$M1_0,
    M2_0 = history$M2_0,
    M1_lag = history$M1_lag,
    M2_lag = history$M2_lag,
    A = sp$outer_A[1L],
    stringsAsFactors = FALSE
  )

  if ("Y_lag" %in% names(history)) row$Y_lag <- history$Y_lag

  for (L_node in .ltmle_exact_L_nodes(node_spec)) {
    lag_nm <- .ltmle_exact_lag_name(L_node)
    if (lag_nm %in% names(history)) row[[lag_nm]] <- history[[lag_nm]]
  }

  for (nm in as.character(node_spec$baseline_vars %||% character(0))) {
    if (nm %in% names(history) && !nm %in% names(row)) row[[nm]] <- history[[nm]]
  }

  if (".fold_id" %in% names(history)) row$.fold_id <- history$.fold_id

  row$M1 <- row$M1_0
  row$M2 <- row$M2_0

  for (L_node in .ltmle_exact_L_nodes(node_spec)) {
    if (!L_node %in% names(row)) row[[L_node]] <- NA_real_
  }
  row$Y <- NA_real_

  row <- .ltmle_exact_add_terms(row)
  long <- .ltmle_exact_generated_row_to_long_block(row, tt = 1L, node_spec = node_spec)
  long <- .ltmle_exact_add_terms(long)

  list(row = row, long = long, id0 = row$id0)
}

.ltmle_exact_terminal_root_plugin_deterministic <- function(dat_wide, T, spec, models,
                                                           treatment_models, censoring_models,
                                                           terminal_task, cont_fit,
                                                           epsilon, bounds,
                                                           probability_bounds,
                                                           treat_mech, p_rct,
                                                           node_spec = NULL,
                                                           ltmle_exact_density_ratio_mc_n = 2000L,
                                                           verbose = FALSE) {
  .ltmle_exact_assert_supported_deterministic_terminal_root(terminal_task, node_spec)

  .ltmle_exact_log(
    verbose,
    "[ltmle_exact] deterministic terminal/root plug-in",
    " rows=", nrow(dat_wide),
    " component=", spec$component[1L]
  )

  gen <- .ltmle_exact_make_deterministic_root_plugin_data(dat_wide, spec, node_spec)
  q0_gen <- .ltmle_exact_predict_continuation(cont_fit, gen$row)

  H_gen <- .ltmle_exact_clever_covariate_for_rows(
    task = terminal_task,
    rows = gen$row,
    models = models,
    treatment_models = treatment_models,
    censoring_models = censoring_models,
    T = T,
    spec = spec,
    probability_bounds = probability_bounds,
    treat_mech = treat_mech,
    p_rct = p_rct,
    node_spec = node_spec,
    row_long = gen$long,
    row_source = "generated",
    ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
    max_t = terminal_task$t
  )

  qstar_gen <- .ltmle_exact_apply_fluctuation(
    q0_new = q0_gen,
    H_new = H_gen,
    epsilon = epsilon,
    bounds = bounds
  )

  if (length(qstar_gen) != nrow(dat_wide)) {
    .stop("Deterministic terminal/root plug-in did not return one value per subject.")
  }

  H_finite <- H_gen[is.finite(H_gen)]
  H_new_mean <- if (length(H_finite)) mean(H_finite, na.rm = TRUE) else NA_real_
  H_new_max <- if (length(H_finite)) max(abs(H_finite), na.rm = TRUE) else NA_real_

  density_diag <- attr(H_gen, "density_ratio_diagnostics")
  if (!is.null(density_diag) && nrow(density_diag)) {
    density_diag$plugin_batch <- NA_integer_
    density_diag$terminal_plugin_type <- "deterministic_root"
  }

  terminal_diag <- data.frame(
    component = spec$component[1L],
    terminal_plugin_type = "deterministic_root",
    terminal_plugin_mc_active = FALSE,
    effective_mc_n = 1L,
    n_batches = 1L,
    subject_rows = nrow(dat_wide),
    batch_mean_min = mean(qstar_gen, na.rm = TRUE),
    batch_mean_max = mean(qstar_gen, na.rm = TRUE),
    batch_mean_sd = 0,
    mc_standard_error = 0,
    stringsAsFactors = FALSE
  )

  list(
    subject_q = as.numeric(qstar_gen),
    subject_q_initial = as.numeric(q0_gen),
    initial_mean = mean(q0_gen, na.rm = TRUE),
    targeted_mean = mean(qstar_gen, na.rm = TRUE),
    mean_targeted_minus_initial = mean(qstar_gen - q0_gen, na.rm = TRUE),
    mc_integration_diagnostics = terminal_diag,
    density_ratio_diagnostics = density_diag,
    H_new_mean = H_new_mean,
    H_new_max = H_new_max,
    n_batches = 1L,
    terminal_plugin_type = "deterministic_root",
    terminal_plugin_mc_active = FALSE,
    effective_mc_n = 1L
  )
}

.ltmle_exact_get_observed_pseudooutcome <- function(task_graph, task_id, observed_rows) {
  task <- task_graph$tasks[[task_id]]
  source <- task$observed_pseudooutcome_source_task_id
  if (identical(source, "observed_terminal_outcome")) return(as.numeric(observed_rows$Y))
  .stop("nonterminal_pseudooutcome_legacy_lookup_forbidden: task=", task_id,
        ", source_task_id=", source)
  out <- task_graph$observed_values[[source]]
  if (is.null(out)) .stop("Missing observed pseudo-outcome source for task ", task_id, ": ", source)
  out
}

.ltmle_exact_store_targeted_continuation <- function(task_graph, task_id,
                                                    qstar_observed,
                                                    qstar_generated = NULL,
                                                    store_generated = TRUE) {
  if (any(!is.finite(qstar_observed))) {
    .stop("Task graph received non-finite observed targeted continuation for ", task_id)
  }

  task_graph$observed_values[[task_id]] <- as.numeric(qstar_observed)

  if (isTRUE(store_generated)) {
    if (is.null(qstar_generated) || any(!is.finite(qstar_generated))) {
      .stop("Task graph received non-finite generated targeted continuation for ", task_id)
    }
    task_graph$generated_values[[task_id]] <- as.numeric(qstar_generated)
  }

  task_graph
}

.ltmle_exact_generated_history_process_for_task <- function(task) {
  .ltmle_exact_process_type(task)
}

.ltmle_exact_generated_row_for_task <- function(task, generated_histories) {
  process_type <- .ltmle_exact_generated_history_process_for_task(task)
  if (process_type == "outcome_process") {
    return(generated_histories$outcome_process[[task$t]]$row)
  }
  if (process_type == "post_mediator_covariate_transition") {
    return(generated_histories$outcome_process[[task$t]]$row)
  }
  if (process_type == "virtual_mixed_continuation_task") {
    return(generated_histories$outcome_process[[task$t]]$row)
  }
  if (process_type == "joint_stochastic_mediator_intervention_law") {
    return(generated_histories$joint_stochastic_mediator_intervention_law[[task$t]]$row)
  }
  if (process_type == "first_mediator_stochastic_intervention_law") {
    return(generated_histories$first_mediator_stochastic_intervention_law[[task$t]]$row)
  }
  if (process_type == "second_mediator_stochastic_intervention_law") {
    return(generated_histories$second_mediator_stochastic_intervention_law[[task$t]]$row)
  }
  if (process_type == "observed_mediator_process") {
    return(generated_histories$outcome_process[[task$t]]$row)
  }
  .stop("Unknown process_type for generated row lookup: ", process_type)
}

.ltmle_exact_get_task_data <- function(task, observed_long, generated_histories,
                                       Q_model, spec, node_spec = NULL,
                                       U_observed = NULL,
                                       include_generated = TRUE) {
  task <- as.list(task)

  obs <- observed_long[observed_long$t == task$t, , drop = FALSE]
  obs <- .ltmle_exact_add_terms(obs[order(obs$id), , drop = FALSE])

  if (!nrow(obs)) {
    .stop("Empty observed task data for ", task$component, " node=", task$node, " t=", task$t)
  }

  U <- U_observed
  if (is.null(U) || length(U) != nrow(obs) || any(!is.finite(U))) {
    .stop("Invalid continuation outcome for ", task$component, " node=", task$node, " t=", task$t)
  }

  covs <- .ltmle_exact_task_history_covs(task, Q_model, node_spec)

  missing_obs <- setdiff(covs, names(obs))
  if (length(missing_obs)) {
    .stop("Observed task data missing columns: ", paste(missing_obs, collapse = ", "))
  }

  out <- list(
    U_observed = as.numeric(U),
    H_observed_data = obs,
    H_observed_long = observed_long,
    H_observed_source = "observed",
    H_generated_data = NULL,
    H_generated_long = NULL,
    H_generated_source = "generated",
    conditioning_covariates = covs,
    subject_id = obs$id
  )

  if (isTRUE(include_generated)) {
    if (is.null(generated_histories)) {
      .stop("include_generated=TRUE but generated_histories is NULL for ",
            task$component, " node=", task$node, " t=", task$t)
    }

    gen <- .ltmle_exact_generated_row_for_task(task, generated_histories)
    gen <- .ltmle_exact_add_terms(gen)

    if (!nrow(gen)) {
      .stop("Empty generated task data for ", task$component, " node=", task$node, " t=", task$t)
    }

    missing_gen <- setdiff(covs, names(gen))
    if (length(missing_gen)) {
      .stop("Generated task data missing columns: ", paste(missing_gen, collapse = ", "))
    }

    gen_long <- .ltmle_exact_generated_history_to_long(
      generated_histories = generated_histories,
      process_type = .ltmle_exact_generated_history_process_for_task(task),
      spec = spec,
      node_spec = node_spec
    )

    out$H_generated_data <- gen
    out$H_generated_long <- gen_long
  }

  out
}

.ltmle_exact_subject_sum <- function(x, subject_id) {
  subject_id <- as.integer(subject_id)
  n <- max(subject_id, na.rm = TRUE)
  out <- rowsum(as.numeric(x), group = subject_id, reorder = TRUE)
  ans <- rep(0, n)
  ans[as.integer(rownames(out))] <- as.numeric(out[, 1L])
  ans
}

.ltmle_exact_make_eif_term_rows <- function(task, task_data, D_task, H, U, Q_star) {
  data.frame(
    id = task_data$subject_id,
    component = task$component,
    task_id = task$task_id %||% NA_character_,
    parent_task_id = task$parent_task_id %||% NA_character_,
    t = task$t,
    time = task$t,
    node = task$node,
    process_type = .ltmle_exact_process_type(task),
    role = task$role %||% NA_character_,
    D_task = D_task,
    H = as.numeric(H),
    U = as.numeric(U),
    Q_star = as.numeric(Q_star),
    stringsAsFactors = FALSE
  )
}

.ltmle_exact_process_assigned_A <- function(process_type, sp) {
  process_type <- as.character(process_type)
  if (process_type %in% c(
    "outcome_process",
    "post_mediator_covariate_transition",
    "virtual_mixed_continuation_task",
    "observed_mediator_process"
  )) {
    return(as.numeric(sp$outer_A[1L]))
  }
  if (identical(process_type, "joint_stochastic_mediator_intervention_law")) {
    return(as.numeric(sp$m1_A[1L]))
  }
  if (identical(process_type, "first_mediator_stochastic_intervention_law")) {
    return(as.numeric(sp$m1_A[1L]))
  }
  if (identical(process_type, "second_mediator_stochastic_intervention_law")) {
    return(as.numeric(sp$m2_A[1L]))
  }
  .stop("Unknown process_type for assigned treatment lookup: ", process_type)
}

.ltmle_exact_assigned_label_A <- function(assigned_label, sp, process_type = NA_character_) {
  label <- as.character(assigned_label %||% NA_character_)
  if (length(label) != 1L || is.na(label) || !nzchar(label)) {
    return(.ltmle_exact_process_assigned_A(process_type, sp))
  }
  if (identical(label, "outcome_process_regimen")) {
    return(as.numeric(sp$outer_A[1L]))
  }
  if (identical(label, "joint_mediator_law_regimen")) {
    return(as.numeric(sp$m1_A[1L]))
  }
  if (identical(label, "first_mediator_law_regimen")) {
    return(as.numeric(sp$m1_A[1L]))
  }
  if (identical(label, "second_mediator_law_regimen")) {
    return(as.numeric(sp$m2_A[1L]))
  }
  .ltmle_exact_process_assigned_A(process_type, sp)
}

.ltmle_exact_task_assigned_A <- function(task, spec) {
  task <- as.list(task)
  sp <- spec[spec$t == as.integer(task$t), , drop = FALSE]
  if (!nrow(sp)) .stop("Missing component spec row for task t=", task$t)
  .ltmle_exact_assigned_label_A(
    task$assigned_treatment_regimen_label %||% task$assigned_treatment_regimen %||% NA_character_,
    sp,
    process_type = .ltmle_exact_process_type(task)
  )
}

.ltmle_exact_branch_keys <- function() {
  c("outcome", "joint_law", "first_law", "second_law")
}

.ltmle_exact_task_branch_key <- function(task) {
  task <- as.list(task)
  process_type <- .ltmle_exact_process_type(task)

  if (process_type %in% c(
    "outcome_process",
    "post_mediator_covariate_transition",
    "virtual_mixed_continuation_task",
    "observed_mediator_process"
  )) {
    return("outcome")
  }
  if (identical(process_type, "joint_stochastic_mediator_intervention_law")) {
    return("joint_law")
  }
  if (identical(process_type, "first_mediator_stochastic_intervention_law")) {
    return("first_law")
  }
  if (identical(process_type, "second_mediator_stochastic_intervention_law")) {
    return("second_law")
  }

  .stop("Unknown ltmle_exact process_type for branch-state lookup: ", process_type)
}

.ltmle_exact_task_mediator_role <- function(task) {
  task <- as.list(task)
  if ("mediator_role" %in% names(task) &&
      !is.na(task$mediator_role) &&
      nzchar(as.character(task$mediator_role))) {
    return(as.character(task$mediator_role))
  }

  process_type <- .ltmle_exact_process_type(task)
  node <- as.character(task$node)

  if (identical(process_type, "joint_stochastic_mediator_intervention_law") &&
      identical(node, "M1")) {
    return("joint_target_M1")
  }
  if (identical(process_type, "joint_stochastic_mediator_intervention_law") &&
      identical(node, "M2")) {
    return("joint_target_M2")
  }
  if (identical(process_type, "first_mediator_stochastic_intervention_law") &&
      identical(node, "M1")) {
    return("first_target_M1")
  }
  if (identical(process_type, "first_mediator_stochastic_intervention_law") &&
      identical(node, "M2")) {
    return("first_auxiliary_M2")
  }
  if (identical(process_type, "second_mediator_stochastic_intervention_law") &&
      identical(node, "M1")) {
    return("second_auxiliary_M1")
  }
  if (identical(process_type, "second_mediator_stochastic_intervention_law") &&
      identical(node, "M2")) {
    return("second_target_M2")
  }

  NA_character_
}

.ltmle_exact_branch_state_process_type <- function(state_key) {
  switch(state_key,
    outcome = "outcome_process",
    joint_law = "joint_stochastic_mediator_intervention_law",
    first_law = "first_mediator_stochastic_intervention_law",
    second_law = "second_mediator_stochastic_intervention_law",
    .stop("Unknown branch state key: ", state_key)
  )
}

.ltmle_exact_branch_state_A_label <- function(state_key) {
  switch(state_key,
    outcome = "outcome_process_regimen",
    joint_law = "joint_mediator_law_regimen",
    first_law = "first_mediator_law_regimen",
    second_law = "second_mediator_law_regimen",
    .stop("Unknown branch state key: ", state_key)
  )
}

.ltmle_exact_branch_state_A_at <- function(state_key, spec, tt) {
  sp <- spec[spec$t == as.integer(tt), , drop = FALSE]
  if (!nrow(sp)) .stop("Missing component spec row for branch-state t=", tt)
  .ltmle_exact_process_assigned_A(
    .ltmle_exact_branch_state_process_type(state_key),
    sp
  )
}

.ltmle_exact_branch_state_prepare_rows <- function(rows, state_key, spec, node_spec = NULL) {
  rows <- rows
  if (!"id" %in% names(rows)) rows$id <- seq_len(nrow(rows))
  if (!"id0" %in% names(rows)) rows$id0 <- rows$id
  tt <- if ("t" %in% names(rows)) {
    as.integer(rows$t[1L])
  } else if (nrow(spec)) {
    as.integer(spec$t[1L])
  } else {
    1L
  }
  rows$A <- as.numeric(.ltmle_exact_branch_state_A_at(state_key, spec, tt))
  if (!"M1" %in% names(rows) && "M1_lag" %in% names(rows)) rows$M1 <- rows$M1_lag
  if (!"M2" %in% names(rows) && "M2_lag" %in% names(rows)) rows$M2 <- rows$M2_lag
  for (L_node in .ltmle_exact_L_nodes(node_spec)) {
    if (!L_node %in% names(rows)) {
      lag_nm <- .ltmle_exact_lag_name(L_node)
      rows[[L_node]] <- if (lag_nm %in% names(rows)) rows[[lag_nm]] else NA_real_
    }
  }
  if (!"Y" %in% names(rows)) rows$Y <- NA_real_
  .ltmle_exact_add_terms(rows)
}

.ltmle_exact_branch_state_make_long <- function(rows,
                                                observed_long,
                                                state_key,
                                                spec,
                                                node_spec = NULL,
                                                tt = NULL) {
  rows <- .ltmle_exact_branch_state_prepare_rows(rows, state_key, spec, node_spec)
  tt <- tt %||% if ("t" %in% names(rows)) as.integer(rows$t[1L]) else max(as.integer(spec$t), na.rm = TRUE)
  ids <- as.integer(rows$id)
  long <- observed_long[
    observed_long$t <= tt & observed_long$id %in% ids,
    ,
    drop = FALSE
  ]
  if (!nrow(long)) {
    long <- .ltmle_exact_generated_row_to_long_block(rows, tt = as.integer(tt), node_spec = node_spec)
  } else {
    long <- long[order(match(long$id, ids), long$t), , drop = FALSE]
  }

  for (ss in seq_len(as.integer(tt))) {
    idx <- long$t == ss
    if (any(idx)) {
      long$A[idx] <- .ltmle_exact_branch_state_A_at(state_key, spec, ss)
    }
  }

  idx_cur <- long$t == as.integer(tt)
  if (any(idx_cur)) {
    rr <- rows[match(long$id[idx_cur], rows$id), , drop = FALSE]
    common <- setdiff(intersect(names(rr), names(long)), c("id", "id0", "t"))
    for (nm in common) long[idx_cur, nm] <- rr[[nm]]
  }
  .ltmle_exact_add_terms(long)
}

.ltmle_exact_init_branch_state <- function(base_rows,
                                           observed_long,
                                           spec,
                                           node_spec = NULL) {
  force(observed_long)
  .ltmle_exact_init_branch_state_from_root(
    root_rows = base_rows,
    spec = spec,
    T = max(as.integer(spec$t), na.rm = TRUE),
    node_spec = node_spec
  )
}

.ltmle_exact_new_history <- function(T) {
  vector("list", as.integer(T))
}

.ltmle_exact_new_L_history <- function(T, node_spec = NULL) {
  stats::setNames(
    lapply(.ltmle_exact_L_nodes(node_spec), function(...) .ltmle_exact_new_history(T)),
    .ltmle_exact_L_nodes(node_spec)
  )
}

.ltmle_exact_history_get <- function(hist, tt, default = NULL, node = NULL) {
  tt <- as.integer(tt)
  if (!is.null(node)) {
    hist <- hist[[node]] %||% list()
  }
  out <- hist[[tt]]
  if (is.null(out)) default else out
}

.ltmle_exact_history_set <- function(hist, tt, value, node = NULL) {
  tt <- as.integer(tt)
  value <- as.numeric(value)
  if (!is.null(node)) {
    if (is.null(hist[[node]])) hist[[node]] <- list()
    hist[[node]][[tt]] <- value
    return(hist)
  }
  hist[[tt]] <- value
  hist
}

.ltmle_exact_branch_state_subject_n <- function(branch_state) {
  n <- attr(branch_state, "n_subjects") %||% NA_integer_
  if (is.finite(n)) return(as.integer(n))
  ids <- branch_state$outcome$particles$id0 %||% integer(0)
  if (!length(ids)) return(0L)
  max(as.integer(ids), na.rm = TRUE)
}

.ltmle_exact_prepare_root_particles <- function(root_rows, state_key, spec, T, node_spec = NULL) {
  rows <- root_rows
  if (!"id0" %in% names(rows)) {
    rows$id0 <- if ("id" %in% names(rows)) as.integer(rows$id) else seq_len(nrow(rows))
  }
  rows$id0 <- as.integer(rows$id0)
  rows$id <- seq_len(nrow(rows))
  rows$.particle_id <- seq_len(nrow(rows))
  rows$.branch_weight <- rep(1, nrow(rows))
  rows$.outcome_weight <- rep(1, nrow(rows))
  rows$.joint_law_weight <- rep(1, nrow(rows))
  rows$.first_law_weight <- rep(1, nrow(rows))
  rows$.second_law_weight <- rep(1, nrow(rows))
  rows$.product_join_weight <- rep(1, nrow(rows))
  rows$.integration_weight <- rep(1, nrow(rows))
  tt <- if ("t" %in% names(rows)) as.integer(rows$t[1L]) else 1L
  if (!is.finite(tt) || tt < 1L) tt <- 1L
  rows$t <- tt
  rows$A <- as.numeric(.ltmle_exact_branch_state_A_at(state_key, spec, tt))

  if (!"M1_0" %in% names(rows) && "M1_lag" %in% names(rows)) rows$M1_0 <- rows$M1_lag
  if (!"M2_0" %in% names(rows) && "M2_lag" %in% names(rows)) rows$M2_0 <- rows$M2_lag
  if (!"M1_lag" %in% names(rows)) rows$M1_lag <- rows$M1_0 %||% 0
  if (!"M2_lag" %in% names(rows)) rows$M2_lag <- rows$M2_0 %||% 0
  if (!"M1" %in% names(rows)) rows$M1 <- rows$M1_lag
  if (!"M2" %in% names(rows)) rows$M2 <- rows$M2_lag
  if (!"Y_lag" %in% names(rows) && "Y0" %in% names(rows)) rows$Y_lag <- rows$Y0
  if (!"Y" %in% names(rows)) rows$Y <- NA_real_
  for (L_node in .ltmle_exact_L_nodes(node_spec)) {
    lag_nm <- .ltmle_exact_lag_name(L_node)
    if (!lag_nm %in% names(rows)) rows[[lag_nm]] <- 0
    if (!L_node %in% names(rows)) rows[[L_node]] <- NA_real_
  }
  .ltmle_exact_add_terms(rows)
}

.ltmle_exact_init_state_histories <- function(particles, T, node_spec = NULL) {
  tt <- as.integer(particles$t[1L])
  out <- list(
    A_hist = .ltmle_exact_new_history(T),
    L_hist = .ltmle_exact_new_L_history(T, node_spec),
    Y_hist = .ltmle_exact_new_history(T),
    M1_hist = .ltmle_exact_new_history(T),
    M2_hist = .ltmle_exact_new_history(T)
  )
  for (ss in seq_len(tt)) out$A_hist[[ss]] <- as.numeric(particles$A)
  if ("M1" %in% names(particles)) out$M1_hist[[tt]] <- as.numeric(particles$M1)
  if ("M2" %in% names(particles)) out$M2_hist[[tt]] <- as.numeric(particles$M2)
  if (tt > 1L && "M1_lag" %in% names(particles)) out$M1_hist[[tt - 1L]] <- as.numeric(particles$M1_lag)
  if (tt > 1L && "M2_lag" %in% names(particles)) out$M2_hist[[tt - 1L]] <- as.numeric(particles$M2_lag)
  for (L_node in .ltmle_exact_L_nodes(node_spec)) {
    lag_nm <- .ltmle_exact_lag_name(L_node)
    if (L_node %in% names(particles) && any(is.finite(particles[[L_node]]))) {
      out$L_hist <- .ltmle_exact_history_set(out$L_hist, tt, particles[[L_node]], node = L_node)
    }
    if (tt > 1L && lag_nm %in% names(particles)) {
      out$L_hist <- .ltmle_exact_history_set(out$L_hist, tt - 1L, particles[[lag_nm]], node = L_node)
    }
  }
  if ("Y" %in% names(particles) && any(is.finite(particles$Y))) out$Y_hist[[tt]] <- as.numeric(particles$Y)
  if (tt > 1L && "Y_lag" %in% names(particles)) out$Y_hist[[tt - 1L]] <- as.numeric(particles$Y_lag)
  out
}

.ltmle_exact_init_branch_state_from_root <- function(root_rows,
                                                     spec,
                                                     T,
                                                     node_spec = NULL) {
  T <- as.integer(T)
  if (!is.finite(T) || T < 1L) .stop("Invalid T for branch-state initialization.")
  if (!nrow(root_rows)) .stop("Cannot initialize branch_state from empty root_rows.")

  make_one <- function(state_key) {
    particles <- .ltmle_exact_prepare_root_particles(root_rows, state_key, spec, T, node_spec)
    hist <- .ltmle_exact_init_state_histories(particles, T, node_spec)
    state <- c(
      list(
        state_key = state_key,
        particles = particles,
        rows = particles,
        long = NULL,
        A_label = .ltmle_exact_branch_state_A_label(state_key)
      ),
      hist
    )
    state$M1_target_hist <- .ltmle_exact_new_history(T)
    state$M2_target_hist <- .ltmle_exact_new_history(T)
    state$M1_aux_hist <- .ltmle_exact_new_history(T)
    state$M2_aux_hist <- .ltmle_exact_new_history(T)
    state$M1_target <- if ("M1" %in% names(particles)) as.numeric(particles$M1) else NULL
    state$M2_target <- if ("M2" %in% names(particles)) as.numeric(particles$M2) else NULL
    state$M1_aux <- NULL
    state$M2_aux <- NULL
    state
  }

  out <- stats::setNames(lapply(.ltmle_exact_branch_keys(), make_one), .ltmle_exact_branch_keys())
  out$first_law$M2_aux <- if ("M2" %in% names(out$first_law$particles)) as.numeric(out$first_law$particles$M2) else NULL
  out$second_law$M1_aux <- if ("M1" %in% names(out$second_law$particles)) as.numeric(out$second_law$particles$M1) else NULL
  attr(out, "spec") <- spec
  attr(out, "T") <- T
  attr(out, "node_spec") <- node_spec
  attr(out, "n_subjects") <- length(unique(out$outcome$particles$id0))
  attr(out, "source_eval_legacy_used_for_training_pseudooutcome") <- FALSE
  for (key in .ltmle_exact_branch_keys()) {
    out[[key]]$rows <- .ltmle_exact_materialize_branch_rows(out, key, out[[key]]$particles$t[1L], node_spec)
    out[[key]]$long <- .ltmle_exact_materialize_branch_long(out, key, out[[key]]$particles$t[1L], node_spec)
  }
  out
}

.ltmle_exact_init_branch_state_from_observed_task <- function(observed_task_rows,
                                                              task,
                                                              spec,
                                                              node_spec = NULL) {
  task <- as.list(task)
  rows <- .ltmle_exact_add_terms(observed_task_rows)
  if (!"t" %in% names(rows)) rows$t <- as.integer(task$t)
  rows$t <- as.integer(task$t)
  if (!"id0" %in% names(rows)) {
    rows$id0 <- if ("id" %in% names(rows)) as.integer(rows$id) else seq_len(nrow(rows))
  }
  T <- max(as.integer(spec$t), as.integer(task$t), na.rm = TRUE)
  state <- .ltmle_exact_init_branch_state_from_root(
    root_rows = rows,
    spec = spec,
    T = T,
    node_spec = node_spec
  )
  if (identical(as.character(spec$world_type[1L] %||% NA_character_), "separate") &&
      identical(.ltmle_exact_task_mediator_role(task), "second_target_M2") &&
      "M1" %in% names(rows)) {
    tt <- as.integer(task$t)
    target_m1 <- as.numeric(rows$M1)
    state$first_law$M1_target <- target_m1
    state$first_law$M1_target_hist <- .ltmle_exact_history_set(
      state$first_law$M1_target_hist,
      tt,
      target_m1
    )
  }
  attr(state, "conditioning_boundary_type") <- "observed_task"
  attr(state, "observed_history_used_as_conditioning") <- TRUE
  attr(state, "source_eval_legacy_used_for_training_pseudooutcome") <- FALSE
  state
}

.ltmle_exact_materialize_root_eval <- function(branch_state,
                                               root_task,
                                               spec,
                                               node_spec = NULL) {
  root_task <- as.list(root_task)
  state_key <- .ltmle_exact_task_branch_key(root_task)
  list(
    rows = .ltmle_exact_branch_state_rows(branch_state, root_task, spec, node_spec),
    long = .ltmle_exact_branch_state_long(branch_state, root_task, spec, node_spec),
    state_key = state_key
  )
}

.ltmle_exact_branch_state_rows <- function(branch_state,
                                           task,
                                           spec,
                                           node_spec = NULL) {
  key <- .ltmle_exact_task_branch_key(task)
  metadata <- .ltmle_exact_task_evaluation_metadata(task, spec)
  if (.ltmle_exact_is_mixed_boundary_type(metadata$source_boundary_type)) {
    rows <- .ltmle_exact_materialize_mixed_boundary_rows_at(
      branch_state = branch_state,
      metadata = metadata,
      tt = as.integer(as.list(task)$t),
      node_spec = node_spec
    )
    rows$A <- as.numeric(.ltmle_exact_task_assigned_A(task, spec))
    attr(rows, "spec") <- spec
    return(rows)
  }
  rows <- .ltmle_exact_materialize_branch_rows(branch_state, key, as.integer(as.list(task)$t), node_spec)
  rows$A <- as.numeric(.ltmle_exact_task_assigned_A(task, spec))
  rows
}

.ltmle_exact_branch_state_long <- function(branch_state,
                                           task,
                                           spec,
                                           node_spec = NULL) {
  key <- .ltmle_exact_task_branch_key(task)
  metadata <- .ltmle_exact_task_evaluation_metadata(task, spec)
  if (.ltmle_exact_is_mixed_boundary_type(metadata$source_boundary_type)) {
    long <- .ltmle_exact_materialize_mixed_boundary_long(
      branch_state = branch_state,
      metadata = metadata,
      max_t = as.integer(as.list(task)$t),
      node_spec = node_spec
    )
    tt <- as.integer(as.list(task)$t)
    long$A[as.integer(long$t) == tt] <- as.numeric(.ltmle_exact_task_assigned_A(task, spec))
    attr(long, "spec") <- spec
    return(long)
  }
  long <- .ltmle_exact_materialize_branch_long(branch_state, key, as.integer(as.list(task)$t), node_spec)
  tt <- as.integer(as.list(task)$t)
  long$A[as.integer(long$t) == tt] <- as.numeric(.ltmle_exact_task_assigned_A(task, spec))
  long
}

.ltmle_exact_state_mediator_hist <- function(state, node, tt, default) {
  tt <- as.integer(tt)
  if (identical(node, "M1")) {
    return(.ltmle_exact_history_get(state$M1_hist, tt, default))
  }
  if (identical(node, "M2")) {
    return(.ltmle_exact_history_get(state$M2_hist, tt, default))
  }
  default
}

.ltmle_exact_materialize_branch_rows <- function(branch_state,
                                                 state_key,
                                                 tt,
                                                 node_spec = NULL) {
  spec <- attr(branch_state, "spec")
  if (is.null(spec)) .stop("branch_state is missing component spec metadata.")
  node_spec <- node_spec %||% attr(branch_state, "node_spec")
  state <- branch_state[[state_key]]
  if (is.null(state) || is.null(state$particles)) {
    .stop("Missing particle branch_state for state_key=", state_key)
  }
  tt <- as.integer(tt)
  p <- state$particles
  rows <- p
  rows$id <- seq_len(nrow(rows))
  rows$t <- tt
  rows$A <- as.numeric(.ltmle_exact_branch_state_A_at(state_key, spec, tt))
  state$A_hist[[tt]] <- rows$A

  rows$M1_lag <- if (tt == 1L) {
    rows$M1_0 %||% rows$M1_lag %||% 0
  } else {
    .ltmle_exact_state_mediator_hist(state, "M1", tt - 1L, rows$M1_lag %||% rows$M1 %||% 0)
  }
  rows$M2_lag <- if (tt == 1L) {
    rows$M2_0 %||% rows$M2_lag %||% 0
  } else {
    .ltmle_exact_state_mediator_hist(state, "M2", tt - 1L, rows$M2_lag %||% rows$M2 %||% 0)
  }
  rows$M1 <- .ltmle_exact_state_mediator_hist(state, "M1", tt, rows$M1 %||% rows$M1_lag)
  rows$M2 <- .ltmle_exact_state_mediator_hist(state, "M2", tt, rows$M2 %||% rows$M2_lag)

  for (L_node in .ltmle_exact_L_nodes(node_spec)) {
    lag_nm <- .ltmle_exact_lag_name(L_node)
    rows[[lag_nm]] <- if (tt == 1L) {
      rows[[lag_nm]] %||% 0
    } else {
      .ltmle_exact_history_get(state$L_hist, tt - 1L, rows[[lag_nm]] %||% rows[[L_node]] %||% 0, node = L_node)
    }
    rows[[L_node]] <- .ltmle_exact_history_get(
      state$L_hist,
      tt,
      rows[[L_node]] %||% rows[[lag_nm]] %||% 0,
      node = L_node
    )
  }
  if (!identical(.ltmle_exact_outcome_type(node_spec), "terminal_only")) {
    rows$Y_lag <- if (tt == 1L) {
      rows$Y0 %||% rows$Y_lag %||% 0
    } else {
      .ltmle_exact_history_get(state$Y_hist, tt - 1L, rows$Y_lag %||% rows$Y %||% 0)
    }
  }
  rows$Y <- .ltmle_exact_history_get(state$Y_hist, tt, rows$Y %||% NA_real_)
  .ltmle_exact_add_terms(rows)
}

.ltmle_exact_materialize_branch_long <- function(branch_state,
                                                 state_key,
                                                 max_t,
                                                 node_spec = NULL) {
  spec <- attr(branch_state, "spec")
  if (is.null(spec)) .stop("branch_state is missing component spec metadata.")
  node_spec <- node_spec %||% attr(branch_state, "node_spec")
  max_t <- as.integer(max_t)
  if (!is.finite(max_t) || max_t < 1L) max_t <- 1L
  rows <- vector("list", max_t)
  for (ss in seq_len(max_t)) {
    rows[[ss]] <- .ltmle_exact_materialize_branch_rows(branch_state, state_key, ss, node_spec)
  }
  out <- do.call(rbind, rows)
  out$id <- rep(seq_len(nrow(branch_state[[state_key]]$particles)), times = max_t)
  out <- out[order(out$id, out$t), , drop = FALSE]
  .ltmle_exact_add_terms(out)
}

.ltmle_exact_source_boundary_metadata <- function(consuming_task, source_task, spec) {
  consuming_task <- as.list(consuming_task)
  source_task <- as.list(source_task)
  eval_state <- as.character(
    consuming_task$source_boundary_eval_state %||%
      consuming_task$source_boundary_source_state %||%
      .ltmle_exact_task_branch_key(source_task)
  )
  outcome_state <- as.character(
    consuming_task$source_boundary_outcome_history_state %||% eval_state
  )
  m1_state <- as.character(
    consuming_task$source_boundary_m1_history_state %||% eval_state
  )
  m2_state <- as.character(
    consuming_task$source_boundary_m2_history_state %||% eval_state
  )
  aux_state <- as.character(
    consuming_task$source_boundary_auxiliary_mediator_history_state %||% "not_applicable"
  )
  direction <- as.character(
    consuming_task$source_boundary_direction %||%
      if (identical(.ltmle_exact_task_branch_key(consuming_task), eval_state)) {
        "within_state"
      } else {
        "branch_boundary"
      }
  )
  boundary_type <- as.character(consuming_task$source_boundary_type %||% NA_character_)
  if (is.na(boundary_type) || !nzchar(boundary_type)) {
    states <- unique(c(eval_state, outcome_state, m1_state, m2_state, aux_state))
    states <- states[!is.na(states) & states != "not_applicable"]
    boundary_type <- if (length(states) > 1L) "cross_regimen_mixed" else "pure"
  }
  if (.ltmle_exact_is_virtual_mixed_task(source_task)) {
    boundary_type <- "virtual_mixed_continuation"
    eval_state <- "outcome"
    outcome_state <- "outcome"
    world_type <- as.character(consuming_task$world_type %||% source_task$world_type %||% "natural")
    m1_state <- if (identical(world_type, "joint")) {
      "joint_law"
    } else if (identical(world_type, "separate")) {
      "first_law"
    } else {
      "outcome"
    }
    m2_state <- if (identical(world_type, "joint")) {
      "joint_law"
    } else if (identical(world_type, "separate")) {
      "second_law"
    } else {
      "outcome"
    }
    aux_state <- "not_applicable"
    direction <- "outcome_consumes_virtual_mixed_continuation"
  }

  list(
    source_boundary_type = boundary_type,
    source_boundary_eval_state = eval_state,
    source_boundary_outcome_history_state = outcome_state,
    source_boundary_m1_history_state = m1_state,
    source_boundary_m2_history_state = m2_state,
    source_boundary_auxiliary_mediator_history_state = aux_state,
    source_boundary_outer_regimen = as.character(consuming_task$source_boundary_outer_regimen %||% NA_character_),
    source_boundary_m1_regimen = as.character(consuming_task$source_boundary_m1_regimen %||% NA_character_),
    source_boundary_m2_regimen = as.character(consuming_task$source_boundary_m2_regimen %||% NA_character_),
    source_boundary_direction = direction,
    source_boundary_virtual_mixed_task_id = if (.ltmle_exact_is_virtual_mixed_task(source_task)) {
      as.character(source_task$task_id %||% NA_character_)
    } else {
      as.character(consuming_task$source_boundary_virtual_mixed_task_id %||% NA_character_)
    },
    source_boundary_provenance_rule = as.character(
      consuming_task$source_boundary_provenance_rule %||% paste0("pure_history_from:", eval_state)
    )
  )
}

.ltmle_exact_task_evaluation_metadata <- function(task, spec) {
  task <- as.list(task)
  eval_state <- .ltmle_exact_task_branch_key(task)
  world_type <- as.character(spec$world_type[1L] %||% "natural")
  process_type <- .ltmle_exact_process_type(task)
  target_m1_state <- if (identical(world_type, "joint")) {
    "joint_law"
  } else if (identical(world_type, "separate")) {
    "first_law"
  } else {
    "outcome"
  }
  target_m2_state <- if (identical(world_type, "joint")) {
    "joint_law"
  } else if (identical(world_type, "separate")) {
    "second_law"
  } else {
    "outcome"
  }
  if (identical(world_type, "natural")) {
    outcome_state <- eval_state
    m1_state <- eval_state
    m2_state <- eval_state
    aux_state <- if (eval_state %in% c("joint_law", "first_law", "second_law")) {
      eval_state
    } else {
      "not_applicable"
    }
  } else if (identical(world_type, "joint")) {
    outcome_state <- "outcome"
    m1_state <- "joint_law"
    m2_state <- "joint_law"
    aux_state <- if (identical(process_type, "joint_stochastic_mediator_intervention_law")) {
      "joint_law"
    } else {
      "not_applicable"
    }
  } else if (identical(world_type, "separate") &&
             identical(process_type, "first_mediator_stochastic_intervention_law")) {
    outcome_state <- "outcome"
    m1_state <- "first_law"
    m2_state <- "first_law"
    aux_state <- "first_law"
  } else if (identical(world_type, "separate") &&
             identical(process_type, "second_mediator_stochastic_intervention_law")) {
    outcome_state <- "outcome"
    m1_state <- "second_law"
    m2_state <- "second_law"
    aux_state <- "second_law"
  } else if (identical(world_type, "separate")) {
    outcome_state <- "outcome"
    m1_state <- target_m1_state
    m2_state <- target_m2_state
    aux_state <- "not_applicable"
  } else {
    outcome_state <- eval_state
    m1_state <- eval_state
    m2_state <- eval_state
    aux_state <- "not_applicable"
  }
  states <- unique(c(eval_state, outcome_state, m1_state, m2_state, aux_state))
  states <- states[!is.na(states) & nzchar(states) & states != "not_applicable"]
  boundary_type <- if (length(states) > 1L) "cross_regimen_mixed" else "pure"
  if (.ltmle_exact_is_virtual_mixed_task(task)) {
    outcome_state <- "outcome"
    m1_state <- target_m1_state
    m2_state <- target_m2_state
    aux_state <- "not_applicable"
    boundary_type <- "virtual_mixed_continuation"
  }
  direction <- if (identical(boundary_type, "pure")) {
    "within_state"
  } else if (identical(boundary_type, "virtual_mixed_continuation")) {
    "virtual_mixed_continuation_rows"
  } else if (identical(eval_state, "outcome")) {
    "outcome_task_uses_mediator_law_history"
  } else {
    "law_task_uses_outer_outcome_history"
  }
  list(
    source_boundary_type = boundary_type,
    source_boundary_eval_state = eval_state,
    source_boundary_outcome_history_state = outcome_state,
    source_boundary_m1_history_state = m1_state,
    source_boundary_m2_history_state = m2_state,
    source_boundary_auxiliary_mediator_history_state = aux_state,
    source_boundary_outer_regimen = as.character(ifelse(spec$outer_A[1L] == 1, "a", "as")),
    source_boundary_m1_regimen = as.character(ifelse(spec$m1_A[1L] == 1, "a", "as")),
    source_boundary_m2_regimen = as.character(ifelse(spec$m2_A[1L] == 1, "a", "as")),
    source_boundary_direction = direction,
    source_boundary_virtual_mixed_task_id = if (.ltmle_exact_is_virtual_mixed_task(task)) {
      as.character(task$task_id %||% NA_character_)
    } else {
      as.character(task$source_boundary_virtual_mixed_task_id %||% NA_character_)
    },
    source_boundary_provenance_rule = if (identical(boundary_type, "cross_regimen_mixed")) {
      paste0(
        "task_eval_mixed_history_from:",
        "eval=", eval_state,
        ";outcome=", outcome_state,
        ";m1=", m1_state,
        ";m2=", m2_state,
        ";aux=", aux_state
      )
    } else if (identical(boundary_type, "virtual_mixed_continuation")) {
      paste0(
        "task_eval_virtual_mixed_history_from:",
        "eval=", eval_state,
        ";outcome=", outcome_state,
        ";m1=", m1_state,
        ";m2=", m2_state,
        ";aux=", aux_state
      )
    } else {
      paste0("task_eval_pure_history_from:", eval_state)
    }
  )
}

.ltmle_exact_state_rows_checked <- function(branch_state, state_key, tt, node_spec = NULL) {
  if (is.na(state_key) || !nzchar(state_key) || identical(state_key, "not_applicable")) {
    return(NULL)
  }
  .ltmle_exact_materialize_branch_rows(
    branch_state = branch_state,
    state_key = state_key,
    tt = as.integer(tt),
    node_spec = node_spec
  )
}

.ltmle_exact_copy_existing_columns <- function(target, source, cols) {
  if (is.null(source)) return(target)
  cols <- intersect(cols, intersect(names(target), names(source)))
  if (!length(cols)) return(target)
  for (nm in cols) target[[nm]] <- source[[nm]]
  target
}

.ltmle_exact_materialize_mixed_boundary_rows_at <- function(branch_state,
                                                           metadata,
                                                           tt,
                                                           node_spec = NULL) {
  eval_state <- metadata$source_boundary_eval_state
  rows <- .ltmle_exact_state_rows_checked(branch_state, eval_state, tt, node_spec)
  outcome_rows <- .ltmle_exact_state_rows_checked(
    branch_state,
    metadata$source_boundary_outcome_history_state,
    tt,
    node_spec
  )
  m1_rows <- .ltmle_exact_state_rows_checked(
    branch_state,
    metadata$source_boundary_m1_history_state,
    tt,
    node_spec
  )
  m2_rows <- .ltmle_exact_state_rows_checked(
    branch_state,
    metadata$source_boundary_m2_history_state,
    tt,
    node_spec
  )
  compatible <- function(x) {
    !is.null(x) &&
      nrow(x) == nrow(rows) &&
      all(as.integer(x$id0 %||% x$id) == as.integer(rows$id0 %||% rows$id))
  }
  outcome_ok <- compatible(outcome_rows)
  m1_ok <- compatible(m1_rows)
  m2_ok <- compatible(m2_rows)
  if (outcome_ok) {
    rows <- .ltmle_exact_copy_existing_columns(
      rows,
      outcome_rows,
      c("Y", "Y_lag", .ltmle_exact_L_nodes(node_spec),
        paste0(.ltmle_exact_L_nodes(node_spec), "_lag"))
    )
  }
  if (m1_ok) {
    rows <- .ltmle_exact_copy_existing_columns(rows, m1_rows, c("M1", "M1_lag"))
  }
  if (m2_ok) {
    rows <- .ltmle_exact_copy_existing_columns(rows, m2_rows, c("M2", "M2_lag"))
  }
  rows$.outcome_history_alignment_ok <- outcome_ok
  rows$.m1_history_alignment_ok <- m1_ok
  rows$.m2_history_alignment_ok <- m2_ok
  rows$.source_boundary_type <- metadata$source_boundary_type
  rows$.virtual_mixed_task_id <- as.character(
    metadata$source_boundary_virtual_mixed_task_id %||% NA_character_
  )
  rows$.source_boundary_direction <- metadata$source_boundary_direction
  rows$.eval_state_key <- metadata$source_boundary_eval_state
  rows$.outcome_history_state_key <- metadata$source_boundary_outcome_history_state
  rows$.m1_history_state_key <- metadata$source_boundary_m1_history_state
  rows$.m2_history_state_key <- metadata$source_boundary_m2_history_state
  rows$.auxiliary_mediator_history_state_key <-
    metadata$source_boundary_auxiliary_mediator_history_state
  rows$.outer_regimen <- metadata$source_boundary_outer_regimen
  rows$.m1_regimen <- metadata$source_boundary_m1_regimen
  rows$.m2_regimen <- metadata$source_boundary_m2_regimen
  rows$.L_history_source <- metadata$source_boundary_outcome_history_state
  rows$.Y_history_source <- metadata$source_boundary_outcome_history_state
  rows$.M1_history_source <- metadata$source_boundary_m1_history_state
  rows$.M2_history_source <- metadata$source_boundary_m2_history_state
  rows$.auxiliary_mediator_history_source <-
    metadata$source_boundary_auxiliary_mediator_history_state
  if (identical(metadata$source_boundary_type, "virtual_mixed_continuation")) {
    rows$.source_boundary_type <- "virtual_mixed_continuation"
  }
  .ltmle_exact_add_terms(rows)
}

.ltmle_exact_materialize_mixed_boundary_long <- function(branch_state,
                                                        metadata,
                                                        max_t,
                                                        node_spec = NULL) {
  max_t <- as.integer(max_t)
  rows <- lapply(seq_len(max_t), function(tt) {
    .ltmle_exact_materialize_mixed_boundary_rows_at(
      branch_state = branch_state,
      metadata = metadata,
      tt = tt,
      node_spec = node_spec
    )
  })
  out <- do.call(rbind, rows)
  out$id <- rep(seq_len(nrow(rows[[1L]])), times = max_t)
  out <- out[order(out$id, out$t), , drop = FALSE]
  .ltmle_exact_add_terms(out)
}

.ltmle_exact_boundary_trace_row <- function(consuming_task,
                                            source_task,
                                            metadata,
                                            rows,
                                            source_fit = NULL,
                                            source_covariates = NULL,
                                            qsource_predicted = FALSE,
                                            cached_targeted_continuation_used = FALSE) {
  consuming_task <- as.list(consuming_task)
  source_task <- as.list(source_task)
  source_covariates <- source_covariates %||%
    as.character(source_fit$conditioning_covariates %||% character(0))
  missing_covariates <- setdiff(source_covariates, names(rows))
  covariate_na <- if (length(source_covariates)) {
    vapply(source_covariates, function(nm) {
      nm %in% names(rows) && any(!is.finite(as.numeric(rows[[nm]])))
    }, logical(1))
  } else {
    logical(0)
  }
  coverage_ok <- !length(missing_covariates) && !any(covariate_na)
  ranges <- source_fit$training_covariate_ranges %||% NULL
  extrap_fraction <- NA_real_
  if (is.data.frame(ranges) && nrow(ranges) && length(source_covariates)) {
    checks <- lapply(source_covariates, function(nm) {
      if (!nm %in% names(rows)) return(rep(NA, nrow(rows)))
      rr <- ranges[ranges$covariate == nm, , drop = FALSE]
      if (!nrow(rr)) return(rep(NA, nrow(rows)))
      val <- as.numeric(rows[[nm]])
      val < as.numeric(rr$min_value[1L]) | val > as.numeric(rr$max_value[1L])
    })
    all_checks <- unlist(checks, use.names = FALSE)
    all_checks <- all_checks[!is.na(all_checks)]
    if (length(all_checks)) extrap_fraction <- mean(all_checks)
  }
  A_expected <- .ltmle_exact_branch_state_A_at(
    metadata$source_boundary_eval_state,
    attr(rows, "spec") %||% data.frame(t = source_task$t),
    as.integer(source_task$t)
  )
  A_used <- if ("A" %in% names(rows) && nrow(rows)) unique(as.numeric(rows$A)) else NA_real_
  A_used_one <- if (length(A_used) == 1L) A_used[1L] else NA_real_
  mixed <- .ltmle_exact_is_mixed_boundary_type(metadata$source_boundary_type)
  virtual_mixed <- identical(metadata$source_boundary_type, "virtual_mixed_continuation")
  alignment_ok <- if (all(c(
    ".outcome_history_alignment_ok",
    ".m1_history_alignment_ok",
    ".m2_history_alignment_ok"
  ) %in% names(rows))) {
    all(as.logical(rows$.outcome_history_alignment_ok) %in% TRUE) &&
      all(as.logical(rows$.m1_history_alignment_ok) %in% TRUE) &&
      all(as.logical(rows$.m2_history_alignment_ok) %in% TRUE)
  } else {
    TRUE
  }
  passed <- isTRUE(coverage_ok) &&
    isTRUE(alignment_ok) &&
    (!isTRUE(mixed) || isTRUE(qsource_predicted)) &&
    (is.na(A_used_one) || abs(A_used_one - A_expected) <= 1e-12)
  data.frame(
    component = consuming_task$component %||% source_task$component %||% NA_character_,
    task_id = consuming_task$task_id %||% NA_character_,
    source_task_id = source_task$task_id %||% NA_character_,
    t = as.integer(consuming_task$t %||% NA_integer_),
    node = consuming_task$node %||% NA_character_,
    process_type = .ltmle_exact_process_type(consuming_task),
    source_boundary_type = metadata$source_boundary_type,
    source_boundary_direction = metadata$source_boundary_direction,
    eval_state_key = metadata$source_boundary_eval_state,
    outcome_history_state_key = metadata$source_boundary_outcome_history_state,
    m1_history_state_key = metadata$source_boundary_m1_history_state,
    m2_history_state_key = metadata$source_boundary_m2_history_state,
    auxiliary_mediator_history_state_key =
      metadata$source_boundary_auxiliary_mediator_history_state,
    outer_regimen = metadata$source_boundary_outer_regimen,
    m1_regimen = metadata$source_boundary_m1_regimen,
    m2_regimen = metadata$source_boundary_m2_regimen,
    A_value_used = A_used_one,
    A_value_expected_for_eval_state = A_expected,
    L_history_source = metadata$source_boundary_outcome_history_state,
    Y_history_source = metadata$source_boundary_outcome_history_state,
    M1_history_source = metadata$source_boundary_m1_history_state,
    M2_history_source = metadata$source_boundary_m2_history_state,
    auxiliary_mediator_history_source =
      metadata$source_boundary_auxiliary_mediator_history_state,
    virtual_mixed_task_used = isTRUE(virtual_mixed),
    virtual_mixed_task_id = as.character(
      metadata$source_boundary_virtual_mixed_task_id %||% NA_character_
    ),
    virtual_mixed_rows_materialized = isTRUE(virtual_mixed) && nrow(rows) > 0L,
    virtual_mixed_state_key = if (isTRUE(virtual_mixed)) {
      metadata$source_boundary_eval_state
    } else {
      NA_character_
    },
    virtual_task_matches_integrated_order = isTRUE(virtual_mixed),
    source_fit_covariates = paste(source_covariates, collapse = "|"),
    mixed_rows_covariate_coverage_ok = isTRUE(coverage_ok),
    outcome_history_alignment_ok = isTRUE(alignment_ok) ||
      !".outcome_history_alignment_ok" %in% names(rows) ||
      all(as.logical(rows$.outcome_history_alignment_ok) %in% TRUE),
    m1_history_alignment_ok = isTRUE(alignment_ok) ||
      !".m1_history_alignment_ok" %in% names(rows) ||
      all(as.logical(rows$.m1_history_alignment_ok) %in% TRUE),
    m2_history_alignment_ok = isTRUE(alignment_ok) ||
      !".m2_history_alignment_ok" %in% names(rows) ||
      all(as.logical(rows$.m2_history_alignment_ok) %in% TRUE),
    mixed_rows_extrapolation_fraction = extrap_fraction,
    mixed_boundary_rows_materialized = isTRUE(mixed),
    pure_law_boundary_rows_used = !isTRUE(mixed) &&
      metadata$source_boundary_eval_state %in% c("joint_law", "first_law", "second_law"),
    pure_outcome_boundary_rows_used = !isTRUE(mixed) &&
      identical(metadata$source_boundary_eval_state, "outcome"),
    qsource_predicted_on_mixed_boundary = isTRUE(mixed) && isTRUE(qsource_predicted),
    cached_targeted_continuation_used = isTRUE(cached_targeted_continuation_used),
    passed = isTRUE(passed),
    failure_class = if (isTRUE(passed)) {
      "no_failure"
    } else if (!isTRUE(alignment_ok)) {
      "mixed_boundary_row_alignment_failed"
    } else if (!isTRUE(coverage_ok)) {
      "mixed_boundary_covariate_coverage_failed"
    } else if (isTRUE(mixed) && !isTRUE(qsource_predicted)) {
      "mixed_boundary_not_used_for_qsource_prediction"
    } else {
      "mixed_boundary_trace_failed"
    },
    stringsAsFactors = FALSE
  )
}

.ltmle_exact_materialize_cross_regimen_source_boundary <- function(branch_state,
                                                                  consuming_task,
                                                                  source_task,
                                                                  source_boundary_metadata,
                                                                  spec,
                                                                  node_spec = NULL) {
  rows <- .ltmle_exact_materialize_mixed_boundary_rows_at(
    branch_state = branch_state,
    metadata = source_boundary_metadata,
    tt = as.integer(source_task$t),
    node_spec = node_spec
  )
  attr(rows, "spec") <- spec
  long <- .ltmle_exact_materialize_mixed_boundary_long(
    branch_state = branch_state,
    metadata = source_boundary_metadata,
    max_t = as.integer(source_task$t),
    node_spec = node_spec
  )
  attr(long, "spec") <- spec
  rows$.source_task_id <- source_task$task_id
  long$.source_task_id <- source_task$task_id
  list(
    rows = rows,
    long = long,
    state_key = source_boundary_metadata$source_boundary_eval_state,
    provenance = source_boundary_metadata,
    diagnostics = .ltmle_exact_boundary_trace_row(
      consuming_task = consuming_task,
      source_task = source_task,
      metadata = source_boundary_metadata,
      rows = rows
    )
  )
}

.ltmle_exact_integrate_virtual_mixed_mediator_path <- function(branch_state,
                                                              virtual_task,
                                                              task_graph,
                                                              models,
                                                              spec,
                                                              node_spec = NULL,
                                                              ltmle_exact_law_integration_n = 5L,
                                                              integration_context = "virtual_mixed_continuation") {
  virtual_task <- as.list(virtual_task)
  if (!.ltmle_exact_is_virtual_mixed_task(virtual_task)) {
    return(list(
      branch_state = branch_state,
      law_integration_diagnostics = data.frame(),
      handoff_trace = data.frame()
    ))
  }
  mediator_task_ids <- c(
    as.character(virtual_task$virtual_mixed_m1_task_id %||% NA_character_),
    as.character(virtual_task$virtual_mixed_m2_task_id %||% NA_character_)
  )
  mediator_task_ids <- mediator_task_ids[!is.na(mediator_task_ids) & nzchar(mediator_task_ids)]
  if (!length(mediator_task_ids)) {
    .stop("virtual_mixed_continuation_missing_mediator_tasks: ", virtual_task$task_id)
  }
  if (identical(as.character(virtual_task$world_type %||% NA_character_), "separate")) {
    m2_task_id <- as.character(virtual_task$virtual_mixed_m2_task_id %||% NA_character_)
    m1_task_id <- as.character(virtual_task$virtual_mixed_m1_task_id %||% NA_character_)
    m2_task <- if (!is.na(m2_task_id) && nzchar(m2_task_id)) task_graph$tasks[[m2_task_id]] else NULL
    if (!is.null(m2_task)) {
      m2_task <- as.list(m2_task)
      all_tasks <- task_graph$tasks
      aux_ids <- names(all_tasks)[vapply(all_tasks, function(x) {
        xx <- as.list(x)
        identical(as.character(xx$component %||% NA_character_), as.character(virtual_task$component)) &&
          identical(as.integer(xx$t %||% NA_integer_), as.integer(m2_task$t)) &&
          identical(as.character(xx$node %||% NA_character_), "M1") &&
          identical(
            as.character(.ltmle_exact_process_type(xx)),
            "second_mediator_stochastic_intervention_law"
          )
      }, logical(1))]
      aux_id <- if (length(aux_ids)) aux_ids[1L] else NA_character_
      mediator_task_ids <- c(aux_id, m2_task_id, m1_task_id)
      mediator_task_ids <- mediator_task_ids[!is.na(mediator_task_ids) & nzchar(mediator_task_ids)]
    }
  }

  state <- branch_state
  law_rows <- list()
  handoff_rows <- list()
  for (mediator_task_id in mediator_task_ids) {
    if (is.null(task_graph$tasks[[mediator_task_id]])) {
      .stop("virtual_mixed_continuation_mediator_task_missing: ", mediator_task_id)
    }
    mediator_task <- as.list(task_graph$tasks[[mediator_task_id]])
    tt <- as.integer(mediator_task$t)
    node <- as.character(mediator_task$node)
    if (!node %in% c("M1", "M2") || tt < 2L) {
      next
    }
    state_key <- .ltmle_exact_task_branch_key(mediator_task)
    rows_current <- .ltmle_exact_branch_state_rows(
      branch_state = state,
      task = mediator_task,
      spec = spec,
      node_spec = node_spec
    )
    n0 <- nrow(rows_current)
    requested_n <- as.integer(ltmle_exact_law_integration_n)
    if (!is.finite(requested_n) || requested_n < 3L) {
      .stop("ltmle_exact law integration grid size must be >= 3 for virtual mixed continuation.")
    }
    if (identical(node, "M1")) {
      mu <- .ltmle_exact_predict_m1(models$M1[[tt]], rows_current)
      sigma <- .ltmle_exact_predict_mediator_sigma(models$M1[[tt]], rows_current)
    } else {
      mu <- .ltmle_exact_predict_m2(models$M2[[tt]], rows_current)
      sigma <- .ltmle_exact_predict_mediator_sigma(models$M2[[tt]], rows_current)
    }
    before_n <- nrow(state[[state_key]]$particles)
    grid <- .ltmle_exact_normal_quantile_grid(requested_n)
    row_idx <- rep(seq_len(n0), times = length(grid$z))
    state_g <- .ltmle_exact_branch_state_expand(state, row_idx)
    values <- rep(mu, times = length(grid$z)) + rep(sigma, times = length(grid$z)) * rep(grid$z, each = n0)
    weights <- rep(grid$w, each = n0)
    state <- .ltmle_exact_branch_state_update_mediator(
      branch_state = state_g,
      task = mediator_task,
      node = node,
      value = values,
      weight_multiplier = weights,
      spec = spec,
      node_spec = node_spec
    )
    after_n <- nrow(state[[state_key]]$particles)
    law_rows[[length(law_rows) + 1L]] <- .ltmle_exact_integration_diag_row(
      mediator_task,
      "continuous_normal",
      paste0(integration_context, "_mediator_path_grid"),
      requested_n,
      requested_n,
      before_n,
      after_n,
      length(unique(as.integer(state[[state_key]]$particles$id0))),
      cap_applied = FALSE,
      pruning_applied = FALSE
    )
    handoff <- attr(state, "last_handoff") %||% data.frame()
    if (!is.null(handoff) && nrow(handoff)) {
      handoff$virtual_mixed_task_id <- as.character(virtual_task$task_id %||% NA_character_)
      handoff_rows[[length(handoff_rows) + 1L]] <- handoff
    }
  }
  list(
    branch_state = state,
    law_integration_diagnostics = .ltmle_exact_rbind_fill(law_rows),
    handoff_trace = .ltmle_exact_rbind_fill(handoff_rows)
  )
}

.ltmle_exact_covariate_ranges <- function(rows, covariates) {
  covariates <- intersect(as.character(covariates), names(rows))
  if (!length(covariates)) {
    return(data.frame(covariate = character(0), min_value = numeric(0), max_value = numeric(0)))
  }
  do.call(rbind, lapply(covariates, function(nm) {
    val <- as.numeric(rows[[nm]])
    data.frame(
      covariate = nm,
      min_value = suppressWarnings(min(val, na.rm = TRUE)),
      max_value = suppressWarnings(max(val, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  }))
}

.ltmle_exact_make_task_source_eval <- function(task,
                                               observed_rows,
                                               observed_long,
                                               spec,
                                               node_spec = NULL,
                                               branch_state = NULL) {
  if (is.null(branch_state)) {
    branch_state <- .ltmle_exact_init_branch_state(
      base_rows = observed_rows,
      observed_long = observed_long,
      spec = spec,
      node_spec = node_spec
    )
  }

  rows <- .ltmle_exact_branch_state_rows(
    branch_state = branch_state,
    task = task,
    spec = spec,
    node_spec = node_spec
  )
  long <- .ltmle_exact_branch_state_long(
    branch_state = branch_state,
    task = task,
    spec = spec,
    node_spec = node_spec
  )

  list(
    rows = rows,
    long = long,
    branch_state = branch_state,
    state_key = .ltmle_exact_task_branch_key(task)
  )
}

.ltmle_exact_branch_state_expand <- function(branch_state, row_idx) {
  row_idx <- as.integer(row_idx)
  out <- branch_state
  for (key in .ltmle_exact_branch_keys()) {
    p <- out[[key]]$particles
    p2 <- p[row_idx, , drop = FALSE]
    p2$id <- seq_along(row_idx)
    p2$.particle_id <- seq_along(row_idx)
    out[[key]]$particles <- .ltmle_exact_add_terms(p2)
    for (nm in c(
      "A_hist", "Y_hist", "M1_hist", "M2_hist",
      "M1_target_hist", "M2_target_hist", "M1_aux_hist", "M2_aux_hist"
    )) {
      if (is.null(out[[key]][[nm]])) next
      out[[key]][[nm]] <- lapply(out[[key]][[nm]], function(x) {
        if (is.null(x)) return(NULL)
        as.numeric(x)[row_idx]
      })
    }
    if (!is.null(out[[key]]$L_hist)) {
      out[[key]]$L_hist <- lapply(out[[key]]$L_hist, function(one_node) {
        lapply(one_node, function(x) {
          if (is.null(x)) return(NULL)
          as.numeric(x)[row_idx]
        })
      })
    }
    for (nm in c("M1_target", "M2_target", "M1_aux", "M2_aux")) {
      if (!is.null(out[[key]][[nm]])) out[[key]][[nm]] <- out[[key]][[nm]][row_idx]
    }
  }
  attr(out, "n_subjects") <- attr(branch_state, "n_subjects")
  out
}

.ltmle_exact_branch_weight_column <- function(state_key) {
  switch(state_key,
    outcome = ".outcome_weight",
    joint_law = ".joint_law_weight",
    first_law = ".first_law_weight",
    second_law = ".second_law_weight",
    .stop("Unknown branch state key for weight ownership: ", state_key)
  )
}

.ltmle_exact_branch_state_weight_sums <- function(branch_state) {
  out <- c(
    outcome_weight = NA_real_,
    joint_law_weight = NA_real_,
    first_law_weight = NA_real_,
    second_law_weight = NA_real_,
    product_join_weight = NA_real_,
    integration_weight = NA_real_
  )
  for (key in .ltmle_exact_branch_keys()) {
    particles <- branch_state[[key]]$particles
    col <- .ltmle_exact_branch_weight_column(key)
    if (col %in% names(particles)) {
      out[paste0(key, "_weight") %||% ""] <- sum(as.numeric(particles[[col]]), na.rm = TRUE)
    }
  }
  if (".product_join_weight" %in% names(branch_state$outcome$particles)) {
    out["product_join_weight"] <- sum(as.numeric(branch_state$outcome$particles$.product_join_weight), na.rm = TRUE)
  }
  if (".integration_weight" %in% names(branch_state$outcome$particles)) {
    out["integration_weight"] <- sum(as.numeric(branch_state$outcome$particles$.integration_weight), na.rm = TRUE)
  }
  out
}

.ltmle_exact_branch_state_multiply_weights <- function(branch_state,
                                                       state_key,
                                                       weight_multiplier = NULL) {
  if (is.null(weight_multiplier)) return(branch_state)
  weight_multiplier <- as.numeric(weight_multiplier)
  if (!state_key %in% .ltmle_exact_branch_keys()) {
    .stop("Unknown branch state key for weight update: ", state_key)
  }
  if (length(weight_multiplier) != nrow(branch_state[[state_key]]$particles)) {
    .stop("weight_multiplier length does not match active branch particle count.")
  }
  weight_col <- .ltmle_exact_branch_weight_column(state_key)
  if (!weight_col %in% names(branch_state[[state_key]]$particles)) {
    branch_state[[state_key]]$particles[[weight_col]] <- rep(1, nrow(branch_state[[state_key]]$particles))
  }
  if (!".integration_weight" %in% names(branch_state[[state_key]]$particles)) {
    branch_state[[state_key]]$particles$.integration_weight <- rep(1, nrow(branch_state[[state_key]]$particles))
  }
  branch_state[[state_key]]$particles[[weight_col]] <-
    as.numeric(branch_state[[state_key]]$particles[[weight_col]]) * weight_multiplier
  branch_state[[state_key]]$particles$.integration_weight <-
    as.numeric(branch_state[[state_key]]$particles$.integration_weight) * weight_multiplier
  branch_state[[state_key]]$particles$.branch_weight <-
    as.numeric(branch_state[[state_key]]$particles$.branch_weight) * weight_multiplier
  branch_state
}

.ltmle_exact_branch_state_refresh_aliases <- function(branch_state, tt, node_spec = NULL) {
  node_spec <- node_spec %||% attr(branch_state, "node_spec")
  for (key in .ltmle_exact_branch_keys()) {
    branch_state[[key]]$rows <- .ltmle_exact_materialize_branch_rows(branch_state, key, tt, node_spec)
    branch_state[[key]]$long <- .ltmle_exact_materialize_branch_long(branch_state, key, tt, node_spec)
  }
  branch_state
}

.ltmle_exact_branch_state_set_node <- function(branch_state,
                                               state_key,
                                               node,
                                               value,
                                               tt,
                                               weight_multiplier = NULL,
                                               node_spec = NULL) {
  value <- as.numeric(value)
  tt <- as.integer(tt)
  if (length(value) != nrow(branch_state[[state_key]]$particles)) {
    .stop("Branch-state update for ", node, " received ", length(value),
          " values for ", nrow(branch_state[[state_key]]$particles), " particles.")
  }
  before_weights <- .ltmle_exact_branch_state_weight_sums(branch_state)
  branch_state <- .ltmle_exact_branch_state_multiply_weights(branch_state, state_key, weight_multiplier)
  branch_state[[state_key]]$particles[[node]] <- value
  if (identical(node, "M1")) {
    branch_state[[state_key]]$M1_hist <- .ltmle_exact_history_set(branch_state[[state_key]]$M1_hist, tt, value)
  } else if (identical(node, "M2")) {
    branch_state[[state_key]]$M2_hist <- .ltmle_exact_history_set(branch_state[[state_key]]$M2_hist, tt, value)
  } else if (node %in% .ltmle_exact_L_nodes(node_spec %||% attr(branch_state, "node_spec"))) {
    branch_state[[state_key]]$L_hist <- .ltmle_exact_history_set(
      branch_state[[state_key]]$L_hist,
      tt,
      value,
      node = node
    )
  } else if (identical(node, "Y")) {
    branch_state[[state_key]]$Y_hist <- .ltmle_exact_history_set(branch_state[[state_key]]$Y_hist, tt, value)
  }
  branch_state <- .ltmle_exact_branch_state_refresh_aliases(branch_state, tt, node_spec)
  after_weights <- .ltmle_exact_branch_state_weight_sums(branch_state)
  changed <- abs(after_weights - before_weights) > 1e-12
  expected_names <- switch(state_key,
    outcome = "outcome_weight",
    joint_law = "joint_law_weight",
    first_law = "first_law_weight",
    second_law = "second_law_weight"
  )
  expected_changed <- names(changed) %in% c(expected_names, "integration_weight") &
    !is.null(weight_multiplier)
  if (!is.null(weight_multiplier)) {
    attr(branch_state, "last_weight_update") <- data.frame(
      state_key = state_key,
      weight_column_updated = .ltmle_exact_branch_weight_column(state_key),
      outcome_weight_changed = isTRUE(changed["outcome_weight"]),
      joint_law_weight_changed = isTRUE(changed["joint_law_weight"]),
      first_law_weight_changed = isTRUE(changed["first_law_weight"]),
      second_law_weight_changed = isTRUE(changed["second_law_weight"]),
      product_join_weight_changed = isTRUE(changed["product_join_weight"]),
      integration_weight_changed = isTRUE(changed["integration_weight"]),
      expected_weight_changed = any(expected_changed & changed),
      unexpected_weight_changed = any(changed & !expected_changed),
      stringsAsFactors = FALSE
    )
  }
  branch_state
}

.ltmle_exact_state_take_indices <- function(state, idx) {
  idx <- as.integer(idx)
  state$particles <- state$particles[idx, , drop = FALSE]
  state$particles$id <- seq_along(idx)
  state$particles$.particle_id <- seq_along(idx)
  if (!is.null(state$rows) && nrow(state$rows)) {
    state$rows <- state$rows[idx, , drop = FALSE]
    state$rows$id <- seq_along(idx)
  }
  for (nm in c(
    "A_hist", "Y_hist", "M1_hist", "M2_hist",
    "M1_target_hist", "M2_target_hist", "M1_aux_hist", "M2_aux_hist"
  )) {
    if (is.null(state[[nm]])) next
    state[[nm]] <- lapply(state[[nm]], function(x) {
      if (is.null(x)) return(NULL)
      as.numeric(x)[idx]
    })
  }
  if (!is.null(state$L_hist)) {
    state$L_hist <- lapply(state$L_hist, function(one_node) {
      lapply(one_node, function(x) {
        if (is.null(x)) return(NULL)
        as.numeric(x)[idx]
      })
    })
  }
  for (nm in c("M1_target", "M2_target", "M1_aux", "M2_aux")) {
    if (!is.null(state[[nm]])) state[[nm]] <- as.numeric(state[[nm]])[idx]
  }
  state
}

.ltmle_exact_unique_target_particles <- function(particles, value, marker_cols = character()) {
  value <- as.numeric(value)
  w <- as.numeric(particles$.branch_weight %||% rep(1, nrow(particles)))
  marker_cols <- intersect(marker_cols, names(particles))
  pieces <- c(
    list(as.character(as.integer(particles$id0))),
    lapply(marker_cols, function(nm) as.character(particles[[nm]])),
    list(as.character(signif(value, 12)), as.character(signif(w, 12)))
  )
  key <- do.call(paste, c(pieces, sep = "\r"))
  particles[!duplicated(key), , drop = FALSE]
}

.ltmle_exact_join_product_weights_for_test <- function(first_weights, second_weights) {
  first_weights <- as.numeric(first_weights)
  second_weights <- as.numeric(second_weights)
  cmb <- expand.grid(
    second_i = seq_along(second_weights),
    first_i = seq_along(first_weights)
  )
  data.frame(
    first_i = cmb$first_i,
    second_i = cmb$second_i,
    weight = first_weights[cmb$first_i] * second_weights[cmb$second_i],
    row_order_copy_used = FALSE,
    first_auxiliary_M2_copied_to_outcome = FALSE,
    second_auxiliary_M1_copied_to_outcome = FALSE,
    stringsAsFactors = FALSE
  )
}

.ltmle_exact_branch_state_product_join_separate_targets <- function(branch_state,
                                                                    tt,
                                                                    node_spec = NULL) {
  tt <- as.integer(tt)
  node_spec <- node_spec %||% attr(branch_state, "node_spec")
  first_value_all <- .ltmle_exact_history_get(branch_state$first_law$M1_target_hist, tt, NULL)
  second_value_all <- .ltmle_exact_history_get(branch_state$second_law$M2_target_hist, tt, NULL)
  if (is.null(first_value_all) || is.null(second_value_all)) {
    .stop("separate_world_product_join_requires_both_target_mediators")
  }
  first_all <- branch_state$first_law$particles
  second_all <- branch_state$second_law$particles
  outcome_all <- branch_state$outcome$particles
  first_all$.orig_index <- seq_len(nrow(first_all))
  second_all$.orig_index <- seq_len(nrow(second_all))
  outcome_all$.orig_index <- seq_len(nrow(outcome_all))
  marker_cols <- .ltmle_exact_particle_marker_cols(first_all, second_all, outcome_all)
  if (length(first_value_all) != nrow(first_all) ||
      length(second_value_all) != nrow(second_all) ||
      any(!is.finite(first_value_all)) ||
      any(!is.finite(second_value_all))) {
    .stop("separate_world_product_join_invalid_target_mediators")
  }
  first_all$.target_value <- as.numeric(first_value_all)
  second_all$.target_value <- as.numeric(second_value_all)
  first <- .ltmle_exact_unique_target_particles(first_all, first_all$.target_value, marker_cols)
  second <- .ltmle_exact_unique_target_particles(second_all, second_all$.target_value, marker_cols)
  outcome <- .ltmle_exact_unique_particle_rows(outcome_all, exclude = ".target_value")
  first$.group_key <- .ltmle_exact_branch_group_key(first, marker_cols)
  second$.group_key <- .ltmle_exact_branch_group_key(second, marker_cols)
  outcome$.group_key <- .ltmle_exact_branch_group_key(outcome, marker_cols)

  out_parts <- list()
  trace_rows <- list()
  first_idx_all <- integer()
  second_idx_all <- integer()
  outcome_idx_all <- integer()
  group_keys <- sort(unique(c(first$.group_key, second$.group_key, outcome$.group_key)))
  for (group_key in group_keys) {
    f <- first[first$.group_key == group_key, , drop = FALSE]
    s <- second[second$.group_key == group_key, , drop = FALSE]
    o <- outcome[outcome$.group_key == group_key, , drop = FALSE]
    if (!nrow(f) || !nrow(s) || !nrow(o)) next
    id0 <- as.integer(o$id0[1L])
    cmb <- expand.grid(
      second_i = seq_len(nrow(s)),
      first_i = seq_len(nrow(f))
    )
    base <- o[rep(1L, nrow(cmb)), , drop = FALSE]
    base$M1 <- f$.target_value[cmb$first_i]
    base$M2 <- s$.target_value[cmb$second_i]
    base$M1_lag <- base$M1_lag %||% base$M1
    base$M2_lag <- base$M2_lag %||% base$M2
    wf <- as.numeric(f$.branch_weight %||% rep(1, nrow(f)))[cmb$first_i]
    ws <- as.numeric(s$.branch_weight %||% rep(1, nrow(s)))[cmb$second_i]
    wo <- rep(1, nrow(cmb))
    product_w <- wf * ws
    base$.product_join_weight <- product_w
    base$.branch_weight <- product_w
    if (!".outcome_weight" %in% names(base)) base$.outcome_weight <- wo
    base$.outcome_weight <- wo
    base$.first_law_weight <- wf
    base$.second_law_weight <- ws
    base$.integration_weight <- wo * product_w
    out_parts[[length(out_parts) + 1L]] <- base
    first_idx_all <- c(first_idx_all, as.integer(f$.orig_index)[cmb$first_i])
    second_idx_all <- c(second_idx_all, as.integer(s$.orig_index)[cmb$second_i])
    outcome_idx_all <- c(outcome_idx_all, rep(as.integer(o$.orig_index)[1L], nrow(cmb)))
    trace_rows[[length(trace_rows) + 1L]] <- data.frame(
      world_type = "separate",
      handoff_type = "separate_world_product_join",
      id0 = id0,
      n_first_particles = nrow(f),
      n_second_particles = nrow(s),
      n_outcome_particles_before_join = nrow(o),
      n_outcome_particles_after_join = nrow(cmb),
      expected_n_outcome_particles = nrow(f) * nrow(s),
      first_particle_id = as.integer(f$.particle_id %||% seq_len(nrow(f)))[1L],
      second_particle_id = as.integer(s$.particle_id %||% seq_len(nrow(s)))[1L],
      outcome_particle_id = as.integer(o$.particle_id %||% seq_len(nrow(o)))[1L],
      weight_first = sum(as.numeric(f$.branch_weight %||% rep(1, nrow(f))), na.rm = TRUE),
      weight_second = sum(as.numeric(s$.branch_weight %||% rep(1, nrow(s))), na.rm = TRUE),
      weight_existing_outcome = 1,
      weight_product = sum(product_w, na.rm = TRUE),
      joined_weight = sum(product_w, na.rm = TRUE),
      expected_joined_weight = sum(wf * ws, na.rm = TRUE),
      row_order_copy_used = FALSE,
      first_auxiliary_M2_copied_to_outcome = FALSE,
      second_auxiliary_M1_copied_to_outcome = FALSE,
      target_M1_copied_to_outcome = TRUE,
      target_M2_copied_to_outcome = TRUE,
      passed = nrow(cmb) == nrow(f) * nrow(s),
      failure_class = if (nrow(cmb) == nrow(f) * nrow(s)) {
        "no_failure"
      } else {
        "separate_world_product_join_count_mismatch"
      },
      stringsAsFactors = FALSE
    )
  }
  if (!length(out_parts)) .stop("Separate-world product join produced no outcome particles.")
  outcome_new <- do.call(rbind, out_parts)
  product_sum_by_id0 <- ave(
    as.numeric(outcome_new$.product_join_weight),
    as.integer(outcome_new$id0),
    FUN = function(x) sum(x, na.rm = TRUE)
  )
  product_sum_by_id0[!is.finite(product_sum_by_id0) | product_sum_by_id0 <= 0] <- 1
  outcome_new$.product_join_weight <- as.numeric(outcome_new$.product_join_weight) / product_sum_by_id0
  outcome_new$.branch_weight <- as.numeric(outcome_new$.branch_weight) / product_sum_by_id0
  outcome_new$.integration_weight <- as.numeric(outcome_new$.integration_weight) / product_sum_by_id0
  outcome_new$.second_law_weight <- as.numeric(outcome_new$.second_law_weight) / product_sum_by_id0
  product_trace <- do.call(rbind, trace_rows)
  product_trace$joined_weight_raw <- as.numeric(product_trace$joined_weight)
  product_trace$expected_joined_weight_raw <- as.numeric(product_trace$expected_joined_weight)
  trace_sum_by_id0 <- ave(
    as.numeric(product_trace$joined_weight),
    as.integer(product_trace$id0),
    FUN = function(x) sum(x, na.rm = TRUE)
  )
  trace_sum_by_id0[!is.finite(trace_sum_by_id0) | trace_sum_by_id0 <= 0] <- 1
  product_trace$weight_second <- as.numeric(product_trace$weight_second) / trace_sum_by_id0
  product_trace$weight_product <- as.numeric(product_trace$weight_product) / trace_sum_by_id0
  product_trace$joined_weight <- as.numeric(product_trace$joined_weight) / trace_sum_by_id0
  product_trace$joined_weight_normalized <- as.numeric(product_trace$joined_weight)
  product_trace$expected_joined_weight <- as.numeric(product_trace$expected_joined_weight) / trace_sum_by_id0
  product_trace$sum_joined_weight_raw_by_event <- ave(
    as.numeric(product_trace$joined_weight_raw),
    as.integer(product_trace$id0),
    FUN = function(x) sum(x, na.rm = TRUE)
  )
  product_trace$sum_joined_weight_normalized_by_event <- ave(
    as.numeric(product_trace$joined_weight_normalized),
    as.integer(product_trace$id0),
    FUN = function(x) sum(x, na.rm = TRUE)
  )
  if (nrow(product_trace) > length(unique(as.integer(product_trace$id0)))) {
    product_trace <- do.call(rbind, lapply(unique(as.integer(product_trace$id0)), function(id) {
      rows <- product_trace[as.integer(product_trace$id0) == id, , drop = FALSE]
      out <- rows[1L, , drop = FALSE]
      out$n_first_particles <- max(as.integer(rows$n_first_particles), na.rm = TRUE)
      out$n_second_particles <- max(as.integer(rows$n_second_particles), na.rm = TRUE)
      out$n_outcome_particles_before_join <- sum(as.integer(rows$n_outcome_particles_before_join), na.rm = TRUE)
      out$n_outcome_particles_after_join <- sum(as.integer(rows$n_outcome_particles_after_join), na.rm = TRUE)
      out$expected_n_outcome_particles <- sum(as.integer(rows$expected_n_outcome_particles), na.rm = TRUE)
      out$weight_first <- sum(as.numeric(rows$weight_first), na.rm = TRUE)
      out$weight_second <- sum(as.numeric(rows$weight_second), na.rm = TRUE)
      out$weight_existing_outcome <- sum(as.numeric(rows$weight_existing_outcome), na.rm = TRUE)
      out$weight_product <- sum(as.numeric(rows$weight_product), na.rm = TRUE)
      out$joined_weight <- sum(as.numeric(rows$joined_weight), na.rm = TRUE)
      out$expected_joined_weight <- sum(as.numeric(rows$expected_joined_weight), na.rm = TRUE)
      out$joined_weight_raw <- sum(as.numeric(rows$joined_weight_raw), na.rm = TRUE)
      out$expected_joined_weight_raw <- sum(as.numeric(rows$expected_joined_weight_raw), na.rm = TRUE)
      out$joined_weight_normalized <- sum(as.numeric(rows$joined_weight_normalized), na.rm = TRUE)
      out$sum_joined_weight_raw_by_event <- sum(as.numeric(rows$joined_weight_raw), na.rm = TRUE)
      out$sum_joined_weight_normalized_by_event <- sum(as.numeric(rows$joined_weight_normalized), na.rm = TRUE)
      out$row_order_copy_used <- any(as.logical(rows$row_order_copy_used) %in% TRUE)
      out$first_auxiliary_M2_copied_to_outcome <- any(as.logical(rows$first_auxiliary_M2_copied_to_outcome) %in% TRUE)
      out$second_auxiliary_M1_copied_to_outcome <- any(as.logical(rows$second_auxiliary_M1_copied_to_outcome) %in% TRUE)
      out$target_M1_copied_to_outcome <- all(as.logical(rows$target_M1_copied_to_outcome) %in% TRUE)
      out$target_M2_copied_to_outcome <- all(as.logical(rows$target_M2_copied_to_outcome) %in% TRUE)
      out$passed <- all(as.logical(rows$passed) %in% TRUE) &&
        identical(as.integer(out$n_outcome_particles_after_join), as.integer(out$expected_n_outcome_particles))
      out$failure_class <- if (isTRUE(out$passed)) {
        "no_failure"
      } else {
        "separate_world_product_join_count_mismatch"
      }
      out
    }))
  }
  outcome_new$.target_value <- NULL
  outcome_new$.orig_index <- NULL
  outcome_new$.group_key <- NULL
  outcome_new$id <- seq_len(nrow(outcome_new))
  outcome_new$.particle_id <- seq_len(nrow(outcome_new))
  branch_state$outcome <- .ltmle_exact_state_take_indices(branch_state$outcome, outcome_idx_all)
  branch_state$outcome$particles <- .ltmle_exact_add_terms(outcome_new)
  branch_state$outcome$M1_target <- as.numeric(outcome_new$M1)
  branch_state$outcome$M2_target <- as.numeric(outcome_new$M2)
  branch_state$outcome$M1_hist <- .ltmle_exact_history_set(branch_state$outcome$M1_hist, tt, outcome_new$M1)
  branch_state$outcome$M2_hist <- .ltmle_exact_history_set(branch_state$outcome$M2_hist, tt, outcome_new$M2)
  branch_state$outcome$M1_target_hist <- .ltmle_exact_history_set(
    branch_state$outcome$M1_target_hist, tt, outcome_new$M1
  )
  branch_state$outcome$M2_target_hist <- .ltmle_exact_history_set(
    branch_state$outcome$M2_target_hist, tt, outcome_new$M2
  )
  attr(branch_state, "last_weight_update") <- data.frame(
    state_key = "outcome",
    weight_column_updated = ".product_join_weight",
    outcome_weight_changed = FALSE,
    joint_law_weight_changed = FALSE,
    first_law_weight_changed = FALSE,
    second_law_weight_changed = FALSE,
    product_join_weight_changed = TRUE,
    integration_weight_changed = TRUE,
    expected_weight_changed = TRUE,
    unexpected_weight_changed = FALSE,
    stringsAsFactors = FALSE
  )

  branch_state$first_law <- .ltmle_exact_state_take_indices(branch_state$first_law, first_idx_all)
  branch_state$second_law <- .ltmle_exact_state_take_indices(branch_state$second_law, second_idx_all)
  branch_state$joint_law <- .ltmle_exact_state_take_indices(branch_state$joint_law, outcome_idx_all)
  branch_state <- .ltmle_exact_branch_state_refresh_aliases(branch_state, tt, node_spec)
  attr(branch_state, "last_handoff_product_trace") <- product_trace
  branch_state
}

.ltmle_exact_branch_state_update_mediator <- function(branch_state,
                                                      task,
                                                      node,
                                                      value,
                                                      weight_multiplier = NULL,
                                                      spec,
                                                      node_spec = NULL) {
  if (missing(spec) && is.data.frame(weight_multiplier)) {
    spec <- weight_multiplier
    weight_multiplier <- NULL
  }
  task <- as.list(task)
  process_type <- .ltmle_exact_process_type(task)
  role <- .ltmle_exact_task_mediator_role(task)
  tt <- as.integer(task$t)
  key <- .ltmle_exact_task_branch_key(task)
  state <- branch_state
  value <- as.numeric(value)
  is_separate_world <- identical(spec$world_type[1L], "separate")
  separate_targets_ready <- function(st) {
    is_separate_world &&
      !is.null(.ltmle_exact_history_get(st$first_law$M1_target_hist, tt, NULL)) &&
      !is.null(.ltmle_exact_history_get(st$second_law$M2_target_hist, tt, NULL))
  }

  state <- .ltmle_exact_branch_state_set_node(
    state,
    key,
    node,
    value,
    tt,
    weight_multiplier = weight_multiplier,
    node_spec = node_spec
  )

  if (identical(role, "joint_target_M1")) {
    state$joint_law$M1_target <- value
    state$joint_law$M1_hist <- .ltmle_exact_history_set(state$joint_law$M1_hist, tt, value)
  } else if (identical(role, "joint_target_M2")) {
    state$joint_law$M2_target <- value
    state$joint_law$M2_hist <- .ltmle_exact_history_set(state$joint_law$M2_hist, tt, value)
    if (!is.null(state$joint_law$M1_target)) {
      state$outcome$M1_target <- state$joint_law$M1_target
      state$outcome$M1_target_hist <- .ltmle_exact_history_set(
        state$outcome$M1_target_hist, tt, state$joint_law$M1_target
      )
      state <- .ltmle_exact_branch_state_set_node(state, "outcome", "M1", state$joint_law$M1_target, tt, node_spec = node_spec)
    }
    state$outcome$M2_target <- value
    state$outcome$M2_target_hist <- .ltmle_exact_history_set(state$outcome$M2_target_hist, tt, value)
    state <- .ltmle_exact_branch_state_set_node(state, "outcome", "M2", value, tt, node_spec = node_spec)
  } else if (identical(role, "first_target_M1")) {
    state$first_law$M1_target <- value
    state$first_law$M1_target_hist <- .ltmle_exact_history_set(state$first_law$M1_target_hist, tt, value)
    if (!is_separate_world) {
      state$outcome$M1_target <- value
      state$outcome$M1_target_hist <- .ltmle_exact_history_set(state$outcome$M1_target_hist, tt, value)
      state <- .ltmle_exact_branch_state_set_node(state, "outcome", "M1", value, tt, node_spec = node_spec)
    } else {
      second_aux <- .ltmle_exact_history_get(state$second_law$M1_aux_hist, tt, NULL)
      if (is.null(second_aux)) {
        state$second_law$M1_aux <- value
        state$second_law$M1_aux_hist <- .ltmle_exact_history_set(state$second_law$M1_aux_hist, tt, value)
        state <- .ltmle_exact_branch_state_set_node(
          state,
          "second_law",
          "M1",
          value,
          tt,
          node_spec = node_spec
        )
      }
      if (separate_targets_ready(state)) {
        state <- .ltmle_exact_branch_state_product_join_separate_targets(
          branch_state = state,
          tt = tt,
          node_spec = node_spec
        )
      }
    }
  } else if (identical(role, "first_auxiliary_M2")) {
    state$first_law$M2_aux <- value
    state$first_law$M2_aux_hist <- .ltmle_exact_history_set(state$first_law$M2_aux_hist, tt, value)
  } else if (identical(role, "second_auxiliary_M1")) {
    state$second_law$M1_aux <- value
    state$second_law$M1_aux_hist <- .ltmle_exact_history_set(state$second_law$M1_aux_hist, tt, value)
  } else if (identical(role, "second_target_M2")) {
    state$second_law$M2_target <- value
    state$second_law$M2_target_hist <- .ltmle_exact_history_set(state$second_law$M2_target_hist, tt, value)
    if (is_separate_world && separate_targets_ready(state)) {
      state <- .ltmle_exact_branch_state_product_join_separate_targets(
        branch_state = state,
        tt = tt,
        node_spec = node_spec
      )
    } else {
      state$outcome$M2_target <- value
      state$outcome$M2_target_hist <- .ltmle_exact_history_set(state$outcome$M2_target_hist, tt, value)
      state <- .ltmle_exact_branch_state_set_node(state, "outcome", "M2", value, tt, node_spec = node_spec)
    }
  }

  handoff_type <- role
  if (is_separate_world &&
      role %in% c("first_target_M1", "second_target_M2") &&
      separate_targets_ready(state)) {
    handoff_type <- "separate_world_product_join"
  } else if (is_separate_world &&
             role %in% c("first_target_M1", "second_target_M2")) {
    handoff_type <- "separate_world_product_join_deferred_until_both_targets_ready"
  }
  copied <- role %in% c("joint_target_M2", "first_target_M1", "second_target_M2")
  if (is_separate_world && role %in% c("first_target_M1", "second_target_M2")) {
    copied <- identical(handoff_type, "separate_world_product_join")
  }
  product_trace <- attr(state, "last_handoff_product_trace")
  if (!is.null(product_trace) && nrow(product_trace) &&
      identical(handoff_type, "separate_world_product_join")) {
    product_summary <- product_trace[1L, , drop = FALSE]
    copied <- TRUE
  } else {
    product_summary <- data.frame(
      world_type = spec$world_type[1L],
      n_first_particles = NA_integer_,
      n_second_particles = NA_integer_,
      n_outcome_particles_before_join = NA_integer_,
      n_outcome_particles_after_join = NA_integer_,
      expected_n_outcome_particles = NA_integer_,
      first_particle_id = NA_integer_,
      second_particle_id = NA_integer_,
      outcome_particle_id = NA_integer_,
      weight_first = NA_real_,
      weight_second = NA_real_,
      weight_existing_outcome = NA_real_,
      weight_product = NA_real_,
      joined_weight = NA_real_,
      expected_joined_weight = NA_real_,
      row_order_copy_used = FALSE,
      first_auxiliary_M2_copied_to_outcome = FALSE,
      second_auxiliary_M1_copied_to_outcome = FALSE,
      target_M1_copied_to_outcome = FALSE,
      target_M2_copied_to_outcome = FALSE,
      stringsAsFactors = FALSE
    )
  }
  attr(state, "last_handoff") <- data.frame(
    world_type = product_summary$world_type[1L],
    component = task$component %||% NA_character_,
    task_id = task$task_id %||% NA_character_,
    t = tt,
    node = node,
    process_type = process_type,
    role = role,
    handoff_type = handoff_type,
    value_mean = mean(value, na.rm = TRUE),
    value_sd = stats::sd(value, na.rm = TRUE),
    from_state = key,
    to_state = "outcome",
    copied_to_outcome = copied,
    should_copy_to_outcome = copied,
    id0 = NA_integer_,
    n_first_particles = product_summary$n_first_particles[1L],
    n_second_particles = product_summary$n_second_particles[1L],
    n_outcome_particles_before_join = product_summary$n_outcome_particles_before_join[1L],
    n_outcome_particles_after_join = product_summary$n_outcome_particles_after_join[1L],
    expected_n_outcome_particles = product_summary$expected_n_outcome_particles[1L],
    first_particle_id = product_summary$first_particle_id[1L],
    second_particle_id = product_summary$second_particle_id[1L],
    outcome_particle_id = product_summary$outcome_particle_id[1L],
    weight_first = product_summary$weight_first[1L],
    weight_second = product_summary$weight_second[1L],
    weight_existing_outcome = product_summary$weight_existing_outcome[1L],
    weight_product = product_summary$weight_product[1L],
    joined_weight = product_summary$joined_weight[1L],
    expected_joined_weight = product_summary$expected_joined_weight[1L],
    row_order_copy_used = product_summary$row_order_copy_used[1L],
    first_auxiliary_M2_copied_to_outcome = product_summary$first_auxiliary_M2_copied_to_outcome[1L],
    second_auxiliary_M1_copied_to_outcome = product_summary$second_auxiliary_M1_copied_to_outcome[1L],
    target_M1_copied_to_outcome = product_summary$target_M1_copied_to_outcome[1L],
    target_M2_copied_to_outcome = product_summary$target_M2_copied_to_outcome[1L],
    passed = TRUE,
    failure_class = "no_failure",
    stringsAsFactors = FALSE
  )
  if (!is.null(product_trace) && nrow(product_trace) &&
      identical(handoff_type, "separate_world_product_join")) {
    repeated <- attr(state, "last_handoff")
    repeated <- repeated[rep(1L, nrow(product_trace)), , drop = FALSE]
    shared <- intersect(names(product_trace), names(repeated))
    for (nm in shared) repeated[[nm]] <- product_trace[[nm]]
    attr(state, "last_handoff") <- repeated
  }

  state
}

.ltmle_exact_branch_state_update_L <- function(branch_state,
                                               task,
                                               L_values,
                                               weight_multiplier = NULL,
                                               spec,
                                               node_spec = NULL) {
  task <- as.list(task)
  state_key <- .ltmle_exact_task_branch_key(task)
  .ltmle_exact_branch_state_set_node(
    branch_state,
    state_key,
    as.character(task$node),
    L_values,
    as.integer(task$t),
    weight_multiplier = weight_multiplier,
    node_spec = node_spec
  )
}

.ltmle_exact_branch_state_update_Y <- function(branch_state,
                                               task,
                                               Y_values,
                                               weight_multiplier = NULL,
                                               spec,
                                               node_spec = NULL) {
  task <- as.list(task)
  state_key <- .ltmle_exact_task_branch_key(task)
  .ltmle_exact_branch_state_set_node(
    branch_state,
    state_key,
    "Y",
    Y_values,
    as.integer(task$t),
    weight_multiplier = weight_multiplier,
    node_spec = node_spec
  )
}

.ltmle_exact_particle_history_hash <- function(rows) {
  cols <- intersect(
    c("A", "M1_lag", "M2_lag", "Y_lag", "M1", "M2", "Y", .ltmle_exact_L_nodes(attr(rows, "node_spec"))),
    names(rows)
  )
  if (!length(cols) || !nrow(rows)) return(NA_character_)
  vals <- vapply(cols, function(nm) mean(as.numeric(rows[[nm]]), na.rm = TRUE), numeric(1))
  paste(names(vals), signif(vals, 10), sep = "=", collapse = "|")
}

.ltmle_exact_branch_trace_row <- function(branch_state, task, state_key, rows) {
  node_spec <- attr(branch_state, "node_spec")
  attr(rows, "node_spec") <- node_spec
  tt <- as.integer(as.list(task)$t)
  state_hash <- function(key) {
    rr <- .ltmle_exact_materialize_branch_rows(branch_state, key, tt, node_spec)
    attr(rr, "node_spec") <- node_spec
    .ltmle_exact_particle_history_hash(rr)
  }
  expected_state_key <- .ltmle_exact_task_branch_key(task)
  data.frame(
    component = task$component %||% NA_character_,
    task_id = task$task_id %||% NA_character_,
    t = tt,
    node = task$node %||% NA_character_,
    process_type = .ltmle_exact_process_type(task),
    expected_state_key = expected_state_key,
    actual_state_key = state_key,
    active_A_mean = mean(as.numeric(rows$A), na.rm = TRUE),
    active_history_hash = .ltmle_exact_particle_history_hash(rows),
    outcome_history_hash = state_hash("outcome"),
    joint_law_history_hash = state_hash("joint_law"),
    first_law_history_hash = state_hash("first_law"),
    second_law_history_hash = state_hash("second_law"),
    node_value_source = "branch_state_node_history",
    history_value_source = "realized_or_integrated_node_value",
    qstar_used_as_history_value = FALSE,
    q0_used_as_history_value = FALSE,
    continuation_mean_used_as_history_value = FALSE,
    matches_expected_state = identical(expected_state_key, state_key),
    passed = identical(expected_state_key, state_key),
    failure_class = if (identical(expected_state_key, state_key)) "no_failure" else "branch_state_mismatch",
    stringsAsFactors = FALSE
  )
}

.ltmle_exact_collapse_particles_to_subjects <- function(q, branch_state, state_key) {
  q <- as.numeric(q)
  particles <- branch_state[[state_key]]$particles
  if (length(q) != nrow(particles)) {
    .stop("Cannot collapse branch particles: q length does not match active particle count.")
  }
  id0 <- as.integer(particles$id0)
  w <- as.numeric(particles$.branch_weight)
  w[!is.finite(w) | w < 0] <- 0
  n <- .ltmle_exact_branch_state_subject_n(branch_state)
  num <- rowsum(q * w, group = id0, reorder = TRUE)
  den <- rowsum(w, group = id0, reorder = TRUE)
  out <- rep(NA_real_, n)
  idx <- as.integer(rownames(num))
  out[idx] <- as.numeric(num[, 1L]) / pmax(as.numeric(den[, 1L]), 1e-12)
  out
}

.ltmle_exact_normal_quantile_grid <- function(n_grid) {
  n_grid <- as.integer(n_grid)
  if (!is.finite(n_grid) || n_grid < 3L) {
    .stop("ltmle_exact law integration grid size must be >= 3.")
  }
  p <- (seq_len(n_grid) - 0.5) / n_grid
  list(
    z = stats::qnorm(p),
    w = rep(1 / n_grid, n_grid),
    n_grid = n_grid
  )
}

.ltmle_exact_density_ratio_mc_cap <- function(requested_n,
                                             n_evaluation_rows,
                                             max_generated_particles = NULL) {
  requested_n <- as.integer(requested_n)
  if (!is.finite(requested_n) || requested_n < 2L) {
    .stop("ltmle_exact_density_ratio_mc_n must be an integer >= 2.")
  }
  n_evaluation_rows <- as.integer(n_evaluation_rows)
  if (!is.finite(n_evaluation_rows) || n_evaluation_rows <= 0L) {
    return(list(
      requested_n = requested_n,
      effective_n = requested_n,
      cap_applied = FALSE,
      max_generated_particles = 0L
    ))
  }
  if (is.null(max_generated_particles)) {
    opt <- getOption("ltmle_exact.density_ratio_max_generated_particles", NULL)
    env <- Sys.getenv("LTMLE_EXACT_DENSITY_RATIO_MAX_GENERATED_PARTICLES", unset = "")
    max_generated_particles <- if (!is.null(opt)) opt else if (nzchar(env)) env else 1000000L
  }
  max_generated_particles <- as.integer(max_generated_particles)
  if (!is.finite(max_generated_particles) || max_generated_particles < 2L) {
    max_generated_particles <- 1000000L
  }
  effective_n <- min(requested_n, max(2L, as.integer(floor(max_generated_particles / n_evaluation_rows))))
  list(
    requested_n = requested_n,
    effective_n = effective_n,
    cap_applied = !identical(effective_n, requested_n),
    max_generated_particles = as.integer(n_evaluation_rows * effective_n)
  )
}

.ltmle_exact_predict_task_qstar_on_rows <- function(task_fit,
                                                    rows,
                                                    row_long,
                                                    models,
                                                    treatment_models,
                                                    censoring_models,
                                                    T,
                                                    spec,
                                                    probability_bounds,
                                                    treat_mech,
	                                                    p_rct,
	                                                    node_spec = NULL,
	                                                    ltmle_exact_density_ratio_mc_n = 2000L) {
  task <- as.list(task_fit$task)
  rows <- .ltmle_exact_add_terms(rows)

  missing_covs <- setdiff(task_fit$conditioning_covariates, names(rows))
  if (length(missing_covs)) {
    .stop(
      "Rows for task Q* prediction are missing conditioning columns: ",
      paste(missing_covs, collapse = ", "),
      "; task=", task_fit$task_id
    )
  }

  q0 <- .ltmle_exact_predict_continuation(task_fit$cont_fit, rows)
  mc_plan <- .ltmle_exact_density_ratio_mc_cap(
    requested_n = ltmle_exact_density_ratio_mc_n,
    n_evaluation_rows = nrow(rows)
  )
  requested_density_ratio_mc_n <- mc_plan$requested_n
  effective_density_ratio_mc_n <- mc_plan$effective_n
  density_ratio_mc_cap_applied <- isTRUE(mc_plan$cap_applied)
  max_generated_particles <- mc_plan$max_generated_particles

  H <- .ltmle_exact_clever_covariate_for_rows(
    task = task,
    rows = rows,
    models = models,
    treatment_models = treatment_models,
    censoring_models = censoring_models,
    T = T,
    spec = spec,
    probability_bounds = probability_bounds,
    treat_mech = treat_mech,
    p_rct = p_rct,
    node_spec = node_spec,
    row_long = row_long,
    row_source = "generated",
    ltmle_exact_density_ratio_mc_n = effective_density_ratio_mc_n,
    max_t = task$t
  )
  H <- .ltmle_exact_truncate_clever_covariate(
    H,
    truncation = task_fit$truncation,
    task = task,
    estimator_variant = task_fit$estimator_variant %||% "ltmle_exact_quantile_truncated",
    row_evaluation_context = "branch_state"
  )

	  out <- .ltmle_exact_apply_fluctuation(
	    q0_new = q0,
	    H_new = H,
	    epsilon = task_fit$epsilon,
	    bounds = task_fit$bounds
	  )

  if (length(out) != nrow(rows) || any(!is.finite(out))) {
    .stop("Invalid Q* prediction for task=", task_fit$task_id)
  }

  density_diag <- attr(H, "density_ratio_diagnostics")
  if (!is.null(density_diag) && nrow(density_diag)) {
    density_diag$requested_density_ratio_mc_n <- requested_density_ratio_mc_n
    density_diag$effective_density_ratio_mc_n <- effective_density_ratio_mc_n
    density_diag$density_ratio_mc_cap_applied <- density_ratio_mc_cap_applied
    density_diag$requested_law_integration_n <- NA_integer_
    density_diag$effective_law_integration_n <- NA_integer_
    density_diag$law_integration_cap_applied <- FALSE
    density_diag$particle_pruning_applied <- FALSE
    density_diag$max_generated_particles <- max_generated_particles
    density_diag$run_mode <- NA_character_
    density_diag$is_acceptance_gate <- NA
    attr(out, "density_ratio_diagnostics") <- density_diag
  }
  attr(out, "H_new_mean") <- mean(as.numeric(H)[is.finite(H)], na.rm = TRUE)
  attr(out, "H_new_max") <- max(abs(as.numeric(H)[is.finite(H)]), na.rm = TRUE)
  attr(out, "truncation_diagnostics") <- attr(H, "truncation_diagnostics") %||% data.frame()
  attr(out, "clever_covariate_decomposition_diagnostics") <-
    attr(H, "clever_covariate_decomposition_diagnostics") %||% data.frame()
  out
}

.ltmle_exact_integration_diag_row <- function(task,
                                              node_family,
                                              integration_method,
                                              requested_n,
                                              effective_n,
                                              before_n,
                                              after_n,
                                              collapsed_n,
                                              cap_applied = FALSE,
                                              pruning_applied = FALSE,
                                              unsupported = FALSE) {
  task <- as.list(task)
  data.frame(
    component = task$component %||% NA_character_,
    task_id = task$task_id %||% NA_character_,
    t = as.integer(task$t %||% NA_integer_),
    node = task$node %||% NA_character_,
    process_type = .ltmle_exact_process_type(task),
    node_family = node_family,
    integration_method = integration_method,
    requested_density_ratio_mc_n = NA_integer_,
    effective_density_ratio_mc_n = NA_integer_,
    density_ratio_mc_cap_applied = FALSE,
    requested_law_integration_n = as.integer(requested_n),
    effective_law_integration_n = as.integer(effective_n),
    law_integration_cap_applied = isTRUE(cap_applied),
    particle_pruning_applied = isTRUE(pruning_applied),
    max_generated_particles = NA_integer_,
    run_mode = NA_character_,
    is_acceptance_gate = NA,
    particle_count_before_node = as.integer(before_n),
    particle_count_after_node = as.integer(after_n),
    particle_count_after_collapse = as.integer(collapsed_n),
    unsupported_family_fail_closed = isTRUE(unsupported),
    diagnostic_only = FALSE,
    passed = !isTRUE(cap_applied) && !isTRUE(pruning_applied) &&
      !isTRUE(unsupported) && identical(as.integer(effective_n), as.integer(requested_n)),
    failure_class = if (!isTRUE(cap_applied) && !isTRUE(pruning_applied) &&
      !isTRUE(unsupported) && identical(as.integer(effective_n), as.integer(requested_n))) {
      "no_failure"
    } else if (isTRUE(unsupported)) {
      "unsupported_node_family"
    } else {
      "law_integration_cap_or_pruning_applied"
    },
    stringsAsFactors = FALSE
  )
}

.ltmle_exact_combine_eval_diagnostics <- function(...) {
  parts <- list(...)
  out <- list(
    density_ratio_diagnostics = data.frame(),
    truncation_diagnostics = data.frame(),
    clever_covariate_decomposition_diagnostics = data.frame(),
    law_integration_diagnostics = data.frame(),
    H_new_mean = NA_real_,
    H_new_max = NA_real_
  )
  for (part in parts) {
    if (is.null(part)) next
    if (!is.null(part$density_ratio_diagnostics) && nrow(part$density_ratio_diagnostics)) {
      out$density_ratio_diagnostics <- rbind(out$density_ratio_diagnostics, part$density_ratio_diagnostics)
    }
    if (!is.null(part$truncation_diagnostics) && nrow(part$truncation_diagnostics)) {
      out$truncation_diagnostics <- rbind(out$truncation_diagnostics, part$truncation_diagnostics)
    }
    if (!is.null(part$clever_covariate_decomposition_diagnostics) &&
        nrow(part$clever_covariate_decomposition_diagnostics)) {
      out$clever_covariate_decomposition_diagnostics <- rbind(
        out$clever_covariate_decomposition_diagnostics,
        part$clever_covariate_decomposition_diagnostics
      )
    }
    if (!is.null(part$law_integration_diagnostics) && nrow(part$law_integration_diagnostics)) {
      out$law_integration_diagnostics <- rbind(out$law_integration_diagnostics, part$law_integration_diagnostics)
    }
    if (is.finite(part$H_new_mean %||% NA_real_)) out$H_new_mean <- part$H_new_mean
    if (is.finite(part$H_new_max %||% NA_real_)) out$H_new_max <- part$H_new_max
  }
  out
}

.ltmle_exact_eval_continuous_node_integral <- function(task,
                                                       branch_state,
                                                       child_task_id,
                                                       task_graph,
                                                       task_fit_cache,
                                                       models,
                                                       treatment_models,
                                                       censoring_models,
                                                       T,
                                                       spec,
                                                       probability_bounds,
                                                       treat_mech,
                                                       p_rct,
                                                       node_spec = NULL,
                                                       ltmle_exact_density_ratio_mc_n = 2000L,
                                                       ltmle_exact_law_integration_n = 5L) {
  task <- as.list(task)
  state_key <- .ltmle_exact_task_branch_key(task)
  tt <- as.integer(task$t)
  node <- as.character(task$node)
  rows <- .ltmle_exact_branch_state_rows(
    branch_state = branch_state,
    task = task,
    spec = spec,
    node_spec = node_spec
  )
  n0 <- nrow(rows)
  requested_n <- as.integer(ltmle_exact_law_integration_n)
  if (!is.finite(requested_n) || requested_n < 3L) {
    .stop("ltmle_exact law integration grid size must be >= 3 for acceptance recursion.")
  }
  effective_n <- requested_n
  cap_applied <- FALSE
  pruning_applied <- FALSE

  if (node %in% .ltmle_exact_L_nodes(node_spec)) {
    fit <- models$L[[node]][[tt]]
    mu <- .ltmle_exact_predict_nuis(fit, rows, type = "numeric")
    sigma <- rep(.ltmle_exact_sigma(fit), n0)
    update_fun <- .ltmle_exact_branch_state_update_L
  } else if (identical(node, "Y")) {
    fit <- models$Y[[tt]]
    if (is.null(fit)) {
      .stop("Missing Y model for branch-state recursion at t=", tt)
    }
    mu <- .ltmle_exact_predict_nuis(fit, rows, type = "numeric")
    sigma <- rep(.ltmle_exact_sigma(fit), n0)
    update_fun <- .ltmle_exact_branch_state_update_Y
  } else {
    .stop("Continuous branch integration called for unsupported node: ", node)
  }

  chunk_size <- as.integer(getOption("ltmle_exact_branch_chunk_size", 50L))
  if (!is.finite(chunk_size) || chunk_size < 1L) chunk_size <- 1L
  if (n0 > chunk_size) {
    acc <- numeric(n0)
    trace_rows <- list()
    handoff_rows <- list()
    diag_parts <- list()
    starts <- seq.int(1L, n0, by = chunk_size)
    for (start in starts) {
      idx <- seq.int(start, min(n0, start + chunk_size - 1L))
      chunk_state <- .ltmle_exact_branch_state_expand(branch_state, idx)
      chunk_eval <- .ltmle_exact_eval_continuous_node_integral(
        task = task,
        branch_state = chunk_state,
        child_task_id = child_task_id,
        task_graph = task_graph,
        task_fit_cache = task_fit_cache,
        models = models,
        treatment_models = treatment_models,
        censoring_models = censoring_models,
        T = T,
        spec = spec,
        probability_bounds = probability_bounds,
        treat_mech = treat_mech,
	          p_rct = p_rct,
	          node_spec = node_spec,
	          ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
		          ltmle_exact_law_integration_n = ltmle_exact_law_integration_n
		        )
      acc[idx] <- chunk_eval$particle_q
      if (!is.null(chunk_eval$branch_trace) && nrow(chunk_eval$branch_trace)) {
        trace_rows[[length(trace_rows) + 1L]] <- chunk_eval$branch_trace
      }
      if (!is.null(chunk_eval$handoff_trace) && nrow(chunk_eval$handoff_trace)) {
        handoff_rows[[length(handoff_rows) + 1L]] <- chunk_eval$handoff_trace
      }
      diag_parts[[length(diag_parts) + 1L]] <- chunk_eval$diagnostics
    }
    return(list(
      particle_q = acc,
      subject_q = .ltmle_exact_collapse_particles_to_subjects(acc, branch_state, state_key),
      particle_frame = branch_state[[state_key]]$particles,
      branch_trace = .ltmle_exact_rbind_fill(trace_rows),
      handoff_trace = .ltmle_exact_rbind_fill(handoff_rows),
      diagnostics = do.call(.ltmle_exact_combine_eval_diagnostics, diag_parts)
    ))
  }

  grid <- .ltmle_exact_normal_quantile_grid(effective_n)
  ng <- length(grid$z)
  row_idx <- rep(seq_len(n0), times = ng)
  grid_z <- rep(grid$z, each = n0)
  state_g <- .ltmle_exact_branch_state_expand(branch_state, row_idx)
  marker <- .ltmle_exact_parent_marker_name(task)
  state_g <- .ltmle_exact_set_parent_marker(state_g, marker)
  values <- rep(mu, times = ng) + rep(sigma, times = ng) * grid_z
  weights <- rep(grid$w, each = n0)
  if (identical(node, "Y")) {
    state_g <- .ltmle_exact_branch_state_update_Y(state_g, task, values, weights, spec, node_spec)
  } else {
    state_g <- .ltmle_exact_branch_state_update_L(state_g, task, values, weights, spec, node_spec)
  }
  child <- .ltmle_exact_eval_task_on_branch_state(
    task_id = child_task_id,
    branch_state = state_g,
    task_graph = task_graph,
    task_fit_cache = task_fit_cache,
    models = models,
    treatment_models = treatment_models,
    censoring_models = censoring_models,
    T = T,
    spec = spec,
    probability_bounds = probability_bounds,
    treat_mech = treat_mech,
    p_rct = p_rct,
    node_spec = node_spec,
    ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
    ltmle_exact_law_integration_n = ltmle_exact_law_integration_n
  )
  acc <- .ltmle_exact_integrate_child_particle_q(
    child = child,
    marker = marker,
    n0 = n0,
    ng = ng,
    grid_w = grid$w,
    task_id = task$task_id
  )
  if (any(!is.finite(acc))) .stop("Non-finite branch integrated value for task=", task$task_id)
  diag <- .ltmle_exact_integration_diag_row(
    task, "continuous_normal", "normal_quantile_grid",
    requested_n, effective_n, n0, nrow(state_g[[state_key]]$particles), n0,
    cap_applied = cap_applied, pruning_applied = pruning_applied
  )
  diagnostics <- .ltmle_exact_combine_eval_diagnostics(
    child$diagnostics,
    list(law_integration_diagnostics = diag)
  )
  list(
    particle_q = acc,
    subject_q = .ltmle_exact_collapse_particles_to_subjects(acc, branch_state, state_key),
    particle_frame = branch_state[[state_key]]$particles,
    branch_trace = rbind(.ltmle_exact_branch_trace_row(branch_state, task, state_key, rows), child$branch_trace),
    handoff_trace = child$handoff_trace,
    diagnostics = diagnostics
  )
}

.ltmle_exact_eval_mediator_node_integral <- function(task,
                                                     branch_state,
                                                     child_task_id,
                                                     task_graph,
                                                     task_fit_cache,
                                                     models,
                                                     treatment_models,
                                                     censoring_models,
                                                     T,
                                                     spec,
                                                     probability_bounds,
                                                     treat_mech,
                                                     p_rct,
                                                     node_spec = NULL,
                                                     ltmle_exact_density_ratio_mc_n = 2000L,
                                                     ltmle_exact_law_integration_n = 5L) {
  task <- as.list(task)
  state_key <- .ltmle_exact_task_branch_key(task)
  tt <- as.integer(task$t)
  node <- as.character(task$node)
  rows <- .ltmle_exact_branch_state_rows(
    branch_state = branch_state,
    task = task,
    spec = spec,
    node_spec = node_spec
  )
  n0 <- nrow(rows)
  requested_n <- as.integer(ltmle_exact_law_integration_n)
  if (!is.finite(requested_n) || requested_n < 3L) {
    .stop("ltmle_exact law integration grid size must be >= 3 for acceptance recursion.")
  }
  effective_n <- requested_n
  cap_applied <- FALSE
  pruning_applied <- FALSE

  if (tt < 2L) {
    child <- .ltmle_exact_eval_task_on_branch_state(
      task_id = child_task_id,
      branch_state = branch_state,
      task_graph = task_graph,
      task_fit_cache = task_fit_cache,
      models = models,
      treatment_models = treatment_models,
      censoring_models = censoring_models,
      T = T,
      spec = spec,
      probability_bounds = probability_bounds,
      treat_mech = treat_mech,
      p_rct = p_rct,
      node_spec = node_spec,
      ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
      ltmle_exact_law_integration_n = ltmle_exact_law_integration_n
    )
    return(child)
  }

  if (identical(node, "M1")) {
    mu <- .ltmle_exact_predict_m1(models$M1[[tt]], rows)
    sigma <- .ltmle_exact_predict_mediator_sigma(models$M1[[tt]], rows)
  } else if (identical(node, "M2")) {
    mu <- .ltmle_exact_predict_m2(models$M2[[tt]], rows)
    sigma <- .ltmle_exact_predict_mediator_sigma(models$M2[[tt]], rows)
  } else {
    .stop("Mediator branch integration called for unsupported node: ", node)
  }

  chunk_size <- as.integer(getOption("ltmle_exact_branch_chunk_size", 50L))
  if (!is.finite(chunk_size) || chunk_size < 1L) chunk_size <- 1L
  if (n0 > chunk_size) {
    acc <- numeric(n0)
    trace_rows <- list()
    handoff_rows <- list()
    diag_parts <- list()
    starts <- seq.int(1L, n0, by = chunk_size)
    for (start in starts) {
      idx <- seq.int(start, min(n0, start + chunk_size - 1L))
      chunk_state <- .ltmle_exact_branch_state_expand(branch_state, idx)
      chunk_eval <- .ltmle_exact_eval_mediator_node_integral(
        task = task,
        branch_state = chunk_state,
        child_task_id = child_task_id,
        task_graph = task_graph,
        task_fit_cache = task_fit_cache,
        models = models,
        treatment_models = treatment_models,
        censoring_models = censoring_models,
        T = T,
        spec = spec,
        probability_bounds = probability_bounds,
        treat_mech = treat_mech,
        p_rct = p_rct,
        node_spec = node_spec,
        ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
        ltmle_exact_law_integration_n = ltmle_exact_law_integration_n
      )
      acc[idx] <- chunk_eval$particle_q
      if (!is.null(chunk_eval$branch_trace) && nrow(chunk_eval$branch_trace)) {
        trace_rows[[length(trace_rows) + 1L]] <- chunk_eval$branch_trace
      }
      if (!is.null(chunk_eval$handoff_trace) && nrow(chunk_eval$handoff_trace)) {
        handoff_rows[[length(handoff_rows) + 1L]] <- chunk_eval$handoff_trace
      }
      diag_parts[[length(diag_parts) + 1L]] <- chunk_eval$diagnostics
    }
    return(list(
      particle_q = acc,
      subject_q = .ltmle_exact_collapse_particles_to_subjects(acc, branch_state, state_key),
      particle_frame = branch_state[[state_key]]$particles,
      branch_trace = .ltmle_exact_rbind_fill(trace_rows),
      handoff_trace = .ltmle_exact_rbind_fill(handoff_rows),
      diagnostics = do.call(.ltmle_exact_combine_eval_diagnostics, diag_parts)
    ))
  }

  grid <- .ltmle_exact_normal_quantile_grid(effective_n)
  ng <- length(grid$z)
  row_idx <- rep(seq_len(n0), times = ng)
  state_g <- .ltmle_exact_branch_state_expand(branch_state, row_idx)
  marker <- .ltmle_exact_parent_marker_name(task)
  state_g <- .ltmle_exact_set_parent_marker(state_g, marker)
  values <- rep(mu, times = ng) + rep(sigma, times = ng) * rep(grid$z, each = n0)
  state_g <- .ltmle_exact_branch_state_update_mediator(
    state_g,
    task,
    node,
    values,
    weight_multiplier = rep(grid$w, each = n0),
    spec = spec,
    node_spec = node_spec
  )
  child <- .ltmle_exact_eval_task_on_branch_state(
    task_id = child_task_id,
    branch_state = state_g,
    task_graph = task_graph,
    task_fit_cache = task_fit_cache,
    models = models,
    treatment_models = treatment_models,
    censoring_models = censoring_models,
    T = T,
    spec = spec,
    probability_bounds = probability_bounds,
    treat_mech = treat_mech,
    p_rct = p_rct,
    node_spec = node_spec,
    ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
    ltmle_exact_law_integration_n = ltmle_exact_law_integration_n
  )
  acc <- .ltmle_exact_integrate_child_particle_q(
    child = child,
    marker = marker,
    n0 = n0,
    ng = ng,
    grid_w = grid$w,
    task_id = task$task_id
  )
  if (any(!is.finite(acc))) .stop("Non-finite branch mediator integrated value for task=", task$task_id)
  diag <- .ltmle_exact_integration_diag_row(
    task, "continuous_normal", "normal_quantile_grid",
    requested_n, effective_n, n0, nrow(state_g[[state_key]]$particles), n0,
    cap_applied = cap_applied, pruning_applied = pruning_applied
  )
  diagnostics <- .ltmle_exact_combine_eval_diagnostics(
    child$diagnostics,
    list(law_integration_diagnostics = diag)
  )
  list(
    particle_q = acc,
    subject_q = .ltmle_exact_collapse_particles_to_subjects(acc, branch_state, state_key),
    particle_frame = branch_state[[state_key]]$particles,
    branch_trace = rbind(.ltmle_exact_branch_trace_row(branch_state, task, state_key, rows), child$branch_trace),
    handoff_trace = rbind(attr(state_g, "last_handoff"), child$handoff_trace),
    diagnostics = diagnostics
  )
}

.ltmle_exact_eval_task_on_branch_state <- function(task_id,
                                                   branch_state,
                                                   task_graph,
                                                   task_fit_cache,
                                                   models,
                                                   treatment_models,
                                                   censoring_models,
                                                   T,
                                                   spec,
                                                   probability_bounds,
                                                   treat_mech,
                                                   p_rct,
                                                   node_spec = NULL,
                                                   ltmle_exact_density_ratio_mc_n = 2000L,
                                                   ltmle_exact_law_integration_n = 5L) {
  if (identical(task_id, "observed_terminal_outcome")) {
    rows <- .ltmle_exact_materialize_branch_rows(branch_state, "outcome", T, node_spec)
    y <- as.numeric(rows$Y)
    if (any(!is.finite(y))) {
      .stop("Observed terminal outcome requested from branch_state before finite Y was generated.")
    }
    return(list(
      particle_q = y,
      subject_q = .ltmle_exact_collapse_particles_to_subjects(y, branch_state, "outcome"),
      particle_frame = branch_state$outcome$particles,
      branch_trace = data.frame(),
      handoff_trace = data.frame(),
      diagnostics = list(density_ratio_diagnostics = data.frame(), law_integration_diagnostics = data.frame())
    ))
  }

  task_fit <- task_fit_cache[[task_id]]
  if (is.null(task_fit)) .stop("Missing task fit for branch recursion task_id=", task_id)
  task <- as.list(task_fit$task)
  state_key <- .ltmle_exact_task_branch_key(task)
  node <- as.character(task$node)
  process_type <- .ltmle_exact_process_type(task)
	  child_id <- as.character(task$observed_pseudooutcome_source_task_id)
	  if (.ltmle_exact_is_virtual_mixed_task(task)) {
	    if (.ltmle_exact_virtual_mixed_can_use_downstream_cached_source(task)) {
	      return(.ltmle_exact_eval_virtual_mixed_cached_continuation(
	        virtual_task = task,
	        branch_state = branch_state,
	        task_graph = task_graph,
	        task_fit_cache = task_fit_cache,
	        models = models,
	        treatment_models = treatment_models,
	        censoring_models = censoring_models,
	        T = T,
	        spec = spec,
	        probability_bounds = probability_bounds,
	        treat_mech = treat_mech,
	        p_rct = p_rct,
	        node_spec = node_spec,
	        ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
	        ltmle_exact_law_integration_n = ltmle_exact_law_integration_n
	      ))
	    }
	    rows <- .ltmle_exact_branch_state_rows(
	      branch_state = branch_state,
	      task = task,
      spec = spec,
      node_spec = node_spec
    )
    long <- .ltmle_exact_branch_state_long(
      branch_state = branch_state,
      task = task,
      spec = spec,
      node_spec = node_spec
    )
    q_virtual_star <- .ltmle_exact_predict_task_qstar_on_rows(
      task_fit = task_fit,
      rows = rows,
      row_long = long,
      models = models,
      treatment_models = treatment_models,
      censoring_models = censoring_models,
      T = T,
      spec = spec,
      probability_bounds = probability_bounds,
	    treat_mech = treat_mech,
		    p_rct = p_rct,
		    node_spec = node_spec,
		    ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n
		  )
    virtual_downstream_id <- as.character(
      task$observed_pseudooutcome_source_task_id %||% NA_character_
    )
    virtual_downstream_fit <- if (!.ltmle_exact_missing_or_empty(virtual_downstream_id)) {
      task_fit_cache[[virtual_downstream_id]]
    } else {
      NULL
    }
    return(list(
      particle_q = as.numeric(q_virtual_star),
      subject_q = .ltmle_exact_collapse_particles_to_subjects(
        q_virtual_star,
        branch_state,
        state_key
      ),
      particle_frame = branch_state[[state_key]]$particles,
      branch_trace = .ltmle_exact_branch_trace_row(branch_state, task, state_key, rows),
      handoff_trace = data.frame(),
      diagnostics = list(
        density_ratio_diagnostics =
          attr(q_virtual_star, "density_ratio_diagnostics") %||% data.frame(),
        law_integration_diagnostics = data.frame(),
        H_new_mean = attr(q_virtual_star, "H_new_mean") %||% NA_real_,
        H_new_max = attr(q_virtual_star, "H_new_max") %||% NA_real_
      ),
      metadata = list(
        source_eval_mode = "branch_state_dp_continuation_value",
        cached_targeted_continuation_used = TRUE,
        continuation_task_id = task_id,
        continuation_fit_found = TRUE,
        continuation_fit_targeted = TRUE,
        terminal_outcome_base_case_used = FALSE,
        terminal_full_recursion_used = FALSE,
        local_qstar_source_prediction_used = FALSE,
        local_current_task_prediction_used = FALSE,
        virtual_mixed_direct_continuation_used = FALSE,
        direct_continuation_allowed = FALSE,
        downstream_source_can_represent_virtual_target = FALSE,
        dedicated_virtual_Q_target_used = TRUE,
        virtual_mixed_downstream_source_task_id = virtual_downstream_id,
        virtual_mixed_downstream_source_fit_found = !is.null(virtual_downstream_fit),
        virtual_mixed_downstream_source_fit_targeted =
          isTRUE(virtual_downstream_fit$targeted)
      )
    ))
  }
  rows <- .ltmle_exact_branch_state_rows(
    branch_state = branch_state,
    task = task,
    spec = spec,
    node_spec = node_spec
  )

  mediator_process <- process_type %in% c(
    "observed_mediator_process",
    "joint_stochastic_mediator_intervention_law",
    "first_mediator_stochastic_intervention_law",
    "second_mediator_stochastic_intervention_law"
  )
  if (isTRUE(mediator_process) && node %in% c("M1", "M2") && as.integer(task$t) >= 2L) {
    return(.ltmle_exact_eval_mediator_node_integral(
      task = task,
      branch_state = branch_state,
      child_task_id = child_id,
      task_graph = task_graph,
      task_fit_cache = task_fit_cache,
      models = models,
      treatment_models = treatment_models,
      censoring_models = censoring_models,
      T = T,
      spec = spec,
      probability_bounds = probability_bounds,
      treat_mech = treat_mech,
      p_rct = p_rct,
      node_spec = node_spec,
      ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
      ltmle_exact_law_integration_n = ltmle_exact_law_integration_n
    ))
  }

  if (node %in% .ltmle_exact_L_nodes(node_spec)) {
    return(.ltmle_exact_eval_continuous_node_integral(
      task = task,
      branch_state = branch_state,
      child_task_id = child_id,
      task_graph = task_graph,
      task_fit_cache = task_fit_cache,
      models = models,
      treatment_models = treatment_models,
      censoring_models = censoring_models,
      T = T,
      spec = spec,
      probability_bounds = probability_bounds,
      treat_mech = treat_mech,
      p_rct = p_rct,
      node_spec = node_spec,
      ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
      ltmle_exact_law_integration_n = ltmle_exact_law_integration_n
    ))
  }

  if (identical(node, "Y") && !identical(child_id, "observed_terminal_outcome")) {
    return(.ltmle_exact_eval_continuous_node_integral(
      task = task,
      branch_state = branch_state,
      child_task_id = child_id,
      task_graph = task_graph,
      task_fit_cache = task_fit_cache,
      models = models,
      treatment_models = treatment_models,
      censoring_models = censoring_models,
      T = T,
      spec = spec,
      probability_bounds = probability_bounds,
      treat_mech = treat_mech,
      p_rct = p_rct,
      node_spec = node_spec,
      ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
      ltmle_exact_law_integration_n = ltmle_exact_law_integration_n
    ))
  }

  long <- .ltmle_exact_branch_state_long(
    branch_state = branch_state,
    task = task,
    spec = spec,
    node_spec = node_spec
  )
  q <- .ltmle_exact_predict_task_qstar_on_rows(
    task_fit = task_fit,
    rows = rows,
    row_long = long,
    models = models,
    treatment_models = treatment_models,
    censoring_models = censoring_models,
    T = T,
    spec = spec,
    probability_bounds = probability_bounds,
    treat_mech = treat_mech,
    p_rct = p_rct,
    node_spec = node_spec,
    ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n
  )
  diagnostics <- list(
    density_ratio_diagnostics = attr(q, "density_ratio_diagnostics") %||% data.frame(),
    truncation_diagnostics = attr(q, "truncation_diagnostics") %||% data.frame(),
    clever_covariate_decomposition_diagnostics =
      attr(q, "clever_covariate_decomposition_diagnostics") %||% data.frame(),
    law_integration_diagnostics = data.frame(),
    H_new_mean = attr(q, "H_new_mean") %||% NA_real_,
    H_new_max = attr(q, "H_new_max") %||% NA_real_
  )
  list(
    particle_q = as.numeric(q),
    subject_q = .ltmle_exact_collapse_particles_to_subjects(q, branch_state, state_key),
    particle_frame = branch_state[[state_key]]$particles,
    branch_trace = .ltmle_exact_branch_trace_row(branch_state, task, state_key, rows),
    handoff_trace = data.frame(),
    diagnostics = diagnostics
  )
}

.ltmle_exact_eval_targeted_branch_recursion <- function(component,
                                                        root_task_id,
                                                        root_rows,
                                                        task_graph,
                                                        task_fit_cache,
                                                        models,
                                                        treatment_models,
                                                        censoring_models,
                                                        T,
                                                        spec,
                                                        probability_bounds,
                                                        treat_mech,
                                                        p_rct,
                                                        node_spec = NULL,
                                                        ltmle_exact_density_ratio_mc_n = 2000L,
                                                        ltmle_exact_law_integration_n = 5L,
                                                        verbose = FALSE) {
  branch_state <- .ltmle_exact_init_branch_state_from_root(
    root_rows = root_rows,
    spec = spec,
    T = T,
    node_spec = node_spec
  )
  .ltmle_exact_log(verbose, "[ltmle_exact] targeted branch recursion component=", component)
  eval <- .ltmle_exact_eval_task_on_branch_state(
    task_id = root_task_id,
    branch_state = branch_state,
    task_graph = task_graph,
    task_fit_cache = task_fit_cache,
    models = models,
    treatment_models = treatment_models,
    censoring_models = censoring_models,
    T = T,
    spec = spec,
    probability_bounds = probability_bounds,
    treat_mech = treat_mech,
    p_rct = p_rct,
    node_spec = node_spec,
    ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
    ltmle_exact_law_integration_n = ltmle_exact_law_integration_n
  )
  list(
    subject_q = eval$subject_q,
    mean = mean(eval$subject_q, na.rm = TRUE),
    branch_trace = eval$branch_trace,
    handoff_trace = eval$handoff_trace,
    diagnostics = eval$diagnostics
  )
}

.ltmle_exact_empty_allowed_local_source_prediction_map <- function() {
  data.frame(
    component = character(0),
    task_id = character(0),
    process_type = character(0),
    node = character(0),
    source_task_id = character(0),
    source_process_type = character(0),
    source_node = character(0),
    allowed_local_prediction = logical(0),
    allowed_reason_code = character(0),
    requires_same_observed_history = logical(0),
    requires_no_downstream_counterfactual_branch = logical(0),
    stringsAsFactors = FALSE
  )
}

.ltmle_exact_source_task_lookup <- function(task, task_graph) {
  task <- as.list(task)
  task_id <- as.character(task$task_id)
  source_task_id <- as.character(task$observed_pseudooutcome_source_task_id %||% NA_character_)
  source_task_is_terminal_outcome <- identical(source_task_id, "observed_terminal_outcome")

  if (is.na(source_task_id) || !nzchar(source_task_id)) {
    return(list(
      task_id = task_id,
      source_task_id = source_task_id,
      source_task = NULL,
      source_task_t = NA_integer_,
      source_task_node = NA_character_,
      source_task_process_type = NA_character_,
      source_task_is_terminal_outcome = FALSE,
      default_requires_branch = FALSE,
      requires_branch = FALSE,
      requires_branch_exception = TRUE,
      requires_branch_exception_reason = "missing_source_task_id"
    ))
  }

  if (source_task_is_terminal_outcome) {
    return(list(
      task_id = task_id,
      source_task_id = source_task_id,
      source_task = NULL,
      source_task_t = NA_integer_,
      source_task_node = "observed_terminal_outcome",
      source_task_process_type = "observed_terminal_outcome",
      source_task_is_terminal_outcome = TRUE,
      default_requires_branch = FALSE,
      requires_branch = FALSE,
      requires_branch_exception = TRUE,
      requires_branch_exception_reason = "terminal_source_task"
    ))
  }

  source_task <- NULL
  if (!is.null(task_graph$tasks) && source_task_id %in% names(task_graph$tasks)) {
    source_task <- as.list(task_graph$tasks[[source_task_id]])
  } else if (!is.null(task_graph$task_id)) {
    idx <- match(source_task_id, task_graph$task_id)
    if (!is.na(idx)) source_task <- as.list(task_graph[idx, , drop = FALSE])
  }
  if (is.null(source_task)) {
    .stop("source_task_id_not_in_task_graph: ", source_task_id, " for task=", task_id)
  }

  list(
    task_id = task_id,
    source_task_id = source_task_id,
    source_task = source_task,
    source_task_t = as.integer(source_task$t %||% NA_integer_),
    source_task_node = as.character(source_task$node %||% NA_character_),
    source_task_process_type = as.character(.ltmle_exact_process_type(source_task)),
    source_task_is_terminal_outcome = FALSE,
    default_requires_branch = TRUE,
    requires_branch = TRUE,
    requires_branch_exception = FALSE,
    requires_branch_exception_reason = NA_character_
  )
}

.ltmle_exact_match_allowed_local_source_prediction <- function(
  task,
  source_info,
  allowed_local_source_prediction_map
) {
  if (is.null(allowed_local_source_prediction_map) ||
      nrow(allowed_local_source_prediction_map) == 0L) {
    return(list(matched = FALSE, n_matches = 0L, allowed_local_prediction = FALSE,
                allowed_reason_code = NA_character_))
  }

  if (isTRUE(source_info$source_task_is_terminal_outcome)) {
    return(list(matched = FALSE, n_matches = 0L, allowed_local_prediction = FALSE,
                allowed_reason_code = NA_character_))
  }

  required_cols <- names(.ltmle_exact_empty_allowed_local_source_prediction_map())
  missing_cols <- setdiff(required_cols, names(allowed_local_source_prediction_map))
  if (length(missing_cols)) {
    .stop("allowed_local_source_prediction_map.csv is missing required columns: ",
          paste(missing_cols, collapse = ", "))
  }

  task <- as.list(task)
  pt <- .ltmle_exact_process_type(task)
  idx <- which(
    allowed_local_source_prediction_map$task_id == as.character(task$task_id) &
      allowed_local_source_prediction_map$process_type == as.character(pt) &
      allowed_local_source_prediction_map$node == as.character(task$node %||% NA_character_) &
      allowed_local_source_prediction_map$source_task_id == as.character(source_info$source_task_id) &
      allowed_local_source_prediction_map$source_process_type == as.character(source_info$source_task_process_type) &
      allowed_local_source_prediction_map$source_node == as.character(source_info$source_task_node)
  )

  if (length(idx) != 1L) {
    return(list(matched = FALSE, n_matches = length(idx), allowed_local_prediction = FALSE,
                allowed_reason_code = NA_character_))
  }

  row <- allowed_local_source_prediction_map[idx, , drop = FALSE]
  reason <- as.character(row$allowed_reason_code[[1]])
  allowed_reasons <- c(
    "source_task_same_observed_history_by_graph_definition",
    "no_downstream_counterfactual_branch_required"
  )

  if (!reason %in% allowed_reasons) {
    .stop("Invalid allowed_reason_code in allowed_local_source_prediction_map.csv for task=",
          as.character(task$task_id), ": ", reason)
  }

  list(
    matched = TRUE,
    n_matches = 1L,
    allowed_local_prediction = isTRUE(row$allowed_local_prediction[[1]]),
    allowed_reason_code = reason
  )
}

.ltmle_exact_missing_or_empty <- function(x) {
  is.null(x) || length(x) != 1L || is.na(x) || !nzchar(as.character(x))
}

.ltmle_exact_source_task_downstream <- function(task_id, source_task_id, task_graph) {
  task_id <- as.character(task_id)
  source_task_id <- as.character(source_task_id)
  if (.ltmle_exact_missing_or_empty(task_id) ||
      .ltmle_exact_missing_or_empty(source_task_id)) {
    return(FALSE)
  }
  if (identical(source_task_id, "observed_terminal_outcome")) return(TRUE)
  if (is.null(task_graph$tasks) || !task_id %in% names(task_graph$tasks)) return(FALSE)

  seen <- character(0)
  cursor <- task_id
  repeat {
    if (cursor %in% seen) return(FALSE)
    seen <- c(seen, cursor)
    task <- task_graph$tasks[[cursor]]
    if (is.null(task)) return(FALSE)
    next_id <- as.character(task$observed_pseudooutcome_source_task_id %||% NA_character_)
    if (.ltmle_exact_missing_or_empty(next_id)) return(FALSE)
    if (identical(next_id, source_task_id)) return(TRUE)
    if (identical(next_id, "observed_terminal_outcome")) return(FALSE)
    cursor <- next_id
  }
}

.ltmle_exact_validate_source_task_for_cached_continuation <- function(task,
                                                                     source_task_id,
                                                                     task_graph,
                                                                     task_fit_cache) {
  task <- as.list(task)
  task_id <- as.character(task$task_id)
  source_task_id <- as.character(source_task_id)

  if (.ltmle_exact_missing_or_empty(source_task_id)) {
    .stop("missing_source_task_id")
  }
  if (is.null(task_graph$tasks) || !source_task_id %in% names(task_graph$tasks)) {
    .stop("source_task_id_not_in_task_graph: ", source_task_id)
  }
  if (identical(source_task_id, task_id)) {
    .stop("source_task_self_reference: ", source_task_id)
  }
  if (!.ltmle_exact_source_task_downstream(task_id, source_task_id, task_graph)) {
    .stop("source_task_not_downstream: task=", task_id, ", source=", source_task_id)
  }

  source_fit <- task_fit_cache[[source_task_id]]
  if (is.null(source_fit)) {
    .stop("cached_targeted_continuation_missing: ", source_task_id)
  }
  if (!isTRUE(source_fit$targeted)) {
    .stop("cached_continuation_not_targeted: ", source_task_id)
  }
  invisible(TRUE)
}

.ltmle_exact_eval_terminal_source_base_case <- function(branch_state,
                                                        task,
                                                        spec,
                                                        node_spec = NULL) {
  force(task)
  rows <- .ltmle_exact_materialize_branch_rows(branch_state, "outcome", max(as.integer(spec$t)), node_spec)
  y <- as.numeric(rows$Y)
  if (any(!is.finite(y))) {
    .stop("Observed terminal source requested before finite Y was present in branch_state.")
  }
  value <- .ltmle_exact_collapse_particles_to_subjects(y, branch_state, "outcome")
  list(
    value = value,
    id = seq_along(value),
    state_key = rep("outcome", length(value)),
    particle_q = y,
    subject_q = value,
    particle_frame = branch_state$outcome$particles,
    branch_trace = data.frame(),
    handoff_trace = data.frame(),
    diagnostics = list(
      density_ratio_diagnostics = data.frame(),
      law_integration_diagnostics = data.frame(),
      H_new_mean = NA_real_,
      H_new_max = NA_real_
    ),
    metadata = list(
      source_eval_mode = "branch_state_terminal_outcome_base_case",
      branch_state_source_boundary_materialized = TRUE,
      source_boundary_reached = TRUE,
      source_boundary_task_id = "observed_terminal_outcome",
      evaluated_beyond_source_boundary = FALSE,
      terminal_task_evaluated_in_source_path = FALSE,
      max_task_depth_beyond_source_boundary = 0L,
      cached_targeted_continuation_used = FALSE,
      continuation_task_id = "not_applicable",
      continuation_fit_found = FALSE,
      continuation_fit_targeted = FALSE,
      branch_state_downstream_recursion_used = FALSE,
      terminal_full_recursion_used = FALSE,
      terminal_outcome_base_case_used = TRUE,
      local_qstar_source_prediction_used = FALSE,
      local_current_task_prediction_used = FALSE,
      nodes_updated = "observed_terminal_outcome",
      tasks_visited = "observed_terminal_outcome",
      used_observed_long_as_generated_history = FALSE,
      used_assigned_A_observed_long_mutation = FALSE,
      used_observed_long_as_source_history = FALSE,
      product_join_used = FALSE,
      n_source_boundary_rows = nrow(rows),
      n_source_boundary_state_keys = 1L,
      n_particles_before_collapse = length(y),
      n_particles_after_collapse = length(value),
      safe_collapse_applied = length(y) != length(value),
      safe_collapse_key = "id0|state_key|source_task_id|event_id|branch_weight_ownership|density_ratio_context",
      actual_helper_name = ".ltmle_exact_eval_terminal_source_base_case"
    )
  )
}

.ltmle_exact_materialize_source_task_boundary <- function(branch_state,
                                                          task = NULL,
                                                          source_task_id,
                                                          task_graph,
                                                          models,
                                                          treatment_models,
                                                          censoring_models,
                                                          T,
                                                          spec,
                                                          probability_bounds,
                                                          treat_mech,
                                                          p_rct,
                                                          node_spec = NULL,
                                                          ltmle_exact_density_ratio_mc_n = 2000L,
                                                          ltmle_exact_law_integration_n = 5L,
                                                          forbid_terminal_full_recursion = TRUE) {
  force(models)
  force(treatment_models)
  force(censoring_models)
  force(probability_bounds)
  force(treat_mech)
  force(p_rct)
  force(ltmle_exact_density_ratio_mc_n)
  force(ltmle_exact_law_integration_n)
  if (!isTRUE(forbid_terminal_full_recursion)) {
    .stop("terminal_full_recursion_used_in_source_path")
  }
  if (.ltmle_exact_missing_or_empty(source_task_id)) {
    .stop("missing_source_task_id")
  }
  if (identical(source_task_id, "observed_terminal_outcome")) {
    .stop("Terminal source boundary must use the explicit terminal base-case.")
  }
  if (is.null(task_graph$tasks) || !source_task_id %in% names(task_graph$tasks)) {
    .stop("source_task_id_not_in_task_graph: ", source_task_id)
  }

  transition_law_rows <- list()
  handoff_rows <- list()
  state <- branch_state
  current_task <- if (is.null(task)) NULL else as.list(task)

  integrate_transition_to_child_boundary <- function(state, current_task) {
    current_task <- as.list(current_task)
    state_key <- .ltmle_exact_task_branch_key(current_task)
    tt <- as.integer(current_task$t)
    node <- as.character(current_task$node)
    process_type <- .ltmle_exact_process_type(current_task)
    if (identical(process_type, "virtual_mixed_continuation_task")) {
      virtual_step <- .ltmle_exact_integrate_virtual_mixed_mediator_path(
        branch_state = state,
        virtual_task = current_task,
        task_graph = task_graph,
        models = models,
        spec = spec,
        node_spec = node_spec,
        ltmle_exact_law_integration_n = ltmle_exact_law_integration_n,
        integration_context = "virtual_mixed_source_boundary"
      )
      return(list(
        branch_state = virtual_step$branch_state,
        law_diag = virtual_step$law_integration_diagnostics %||% data.frame(),
        handoff_trace = virtual_step$handoff_trace %||% data.frame()
      ))
    }
    rows_current <- .ltmle_exact_branch_state_rows(
      branch_state = state,
      task = current_task,
      spec = spec,
      node_spec = node_spec
    )
    n0 <- nrow(rows_current)
    requested_n <- as.integer(ltmle_exact_law_integration_n)
    effective_n <- requested_n
    before_n <- nrow(state[[state_key]]$particles)

    if (node %in% .ltmle_exact_L_nodes(node_spec)) {
      fit <- models$L[[node]][[tt]]
      mu <- .ltmle_exact_predict_nuis(fit, rows_current, type = "numeric")
      sigma <- rep(.ltmle_exact_sigma(fit), n0)
      updater <- function(st, value, weight) {
        .ltmle_exact_branch_state_update_L(st, current_task, value, weight, spec, node_spec)
      }
    } else if (identical(node, "Y")) {
      fit <- models$Y[[tt]]
      if (is.null(fit)) .stop("Missing Y model for source-boundary materialization at t=", tt)
      mu <- .ltmle_exact_predict_nuis(fit, rows_current, type = "numeric")
      sigma <- rep(.ltmle_exact_sigma(fit), n0)
      updater <- function(st, value, weight) {
        .ltmle_exact_branch_state_update_Y(st, current_task, value, weight, spec, node_spec)
      }
    } else if (process_type %in% c(
      "observed_mediator_process",
      "joint_stochastic_mediator_intervention_law",
      "first_mediator_stochastic_intervention_law",
      "second_mediator_stochastic_intervention_law"
    ) && node %in% c("M1", "M2") && tt >= 2L) {
      if (identical(node, "M1")) {
        mu <- .ltmle_exact_predict_m1(models$M1[[tt]], rows_current)
        sigma <- .ltmle_exact_predict_mediator_sigma(models$M1[[tt]], rows_current)
      } else {
        mu <- .ltmle_exact_predict_m2(models$M2[[tt]], rows_current)
        sigma <- .ltmle_exact_predict_mediator_sigma(models$M2[[tt]], rows_current)
      }
      updater <- function(st, value, weight) {
        .ltmle_exact_branch_state_update_mediator(
          branch_state = st,
          task = current_task,
          node = node,
          value = value,
          weight_multiplier = weight,
          spec = spec,
          node_spec = node_spec
        )
      }
    } else {
      return(list(
        branch_state = state,
        law_diag = data.frame(),
        handoff_trace = data.frame()
      ))
    }

    grid <- .ltmle_exact_normal_quantile_grid(effective_n)
    ng <- length(grid$z)
    row_idx <- rep(seq_len(n0), times = ng)
    state_g <- .ltmle_exact_branch_state_expand(state, row_idx)
    values <- rep(mu, times = ng) + rep(sigma, times = ng) * rep(grid$z, each = n0)
    weights <- rep(grid$w, each = n0)
    state_g <- updater(state_g, values, weights)
    after_n <- nrow(state_g[[state_key]]$particles)
    collapse_n <- length(unique(as.integer(state_g[[state_key]]$particles$id0)))
    law_diag <- data.frame(
      component = current_task$component %||% NA_character_,
      task_id = current_task$task_id %||% NA_character_,
      source_task_id = source_task_id,
      t = tt,
      node = node,
      process_type = process_type,
      requested_density_ratio_mc_n = as.integer(ltmle_exact_density_ratio_mc_n),
      effective_density_ratio_mc_n = as.integer(ltmle_exact_density_ratio_mc_n),
      density_ratio_mc_cap_applied = FALSE,
      requested_law_integration_n = requested_n,
      effective_law_integration_n = effective_n,
      law_integration_cap_applied = FALSE,
      particle_pruning_applied = FALSE,
      max_generated_particles = as.integer(after_n * as.integer(ltmle_exact_density_ratio_mc_n)),
      integration_method = "source_boundary_transition_grid",
      source_eval_mode = "branch_state_dp_continuation_value",
      n_particles_before_integration = before_n,
      n_particles_after_integration = after_n,
      n_particles_after_collapse = collapse_n,
      particle_count_before_node = before_n,
      particle_count_after_node = after_n,
      particle_count_after_collapse = collapse_n,
      safe_collapse_applied = after_n != collapse_n,
      safe_collapse_key = "id0|state_key|source_task_id|event_id|branch_weight_ownership|density_ratio_context",
      run_mode = NA_character_,
      is_acceptance_gate = NA,
      passed = TRUE,
      failure_class = "no_failure",
      stringsAsFactors = FALSE
    )
    list(
      branch_state = state_g,
      law_diag = law_diag,
      handoff_trace = attr(state_g, "last_handoff") %||% data.frame()
    )
  }

  visited <- character(0)
  while (!is.null(current_task) &&
         !identical(as.character(current_task$task_id), as.character(source_task_id))) {
    current_id <- as.character(current_task$task_id)
    if (current_id %in% visited) .stop("branch_state_recursion_cycle_detected")
    visited <- c(visited, current_id)
    child_id <- as.character(current_task$observed_pseudooutcome_source_task_id %||% NA_character_)
    if (.ltmle_exact_missing_or_empty(child_id)) .stop("missing_source_task_id")
    if (identical(child_id, "observed_terminal_outcome")) {
      .stop("terminal_full_recursion_used_in_source_path")
    }
    step <- integrate_transition_to_child_boundary(state, current_task)
    state <- step$branch_state
    if (!is.null(step$law_diag) && nrow(step$law_diag)) {
      transition_law_rows[[length(transition_law_rows) + 1L]] <- step$law_diag
    }
    if (!is.null(step$handoff_trace) && nrow(step$handoff_trace)) {
      handoff_rows[[length(handoff_rows) + 1L]] <- step$handoff_trace
    }
    if (identical(child_id, source_task_id)) break
    if (is.null(task_graph$tasks[[child_id]])) {
      .stop("source_task_id_not_in_task_graph: ", child_id)
    }
    current_task <- as.list(task_graph$tasks[[child_id]])
  }

  source_task <- as.list(task_graph$tasks[[source_task_id]])
  source_boundary_metadata <- .ltmle_exact_source_boundary_metadata(
    consuming_task = current_task %||% task,
    source_task = source_task,
    spec = spec
  )
  if (.ltmle_exact_is_mixed_boundary_type(source_boundary_metadata$source_boundary_type)) {
    mixed_boundary <- .ltmle_exact_materialize_cross_regimen_source_boundary(
      branch_state = state,
      consuming_task = current_task %||% task,
      source_task = source_task,
      source_boundary_metadata = source_boundary_metadata,
      spec = spec,
      node_spec = node_spec
    )
    state_key <- mixed_boundary$state_key
    rows <- mixed_boundary$rows
    long <- mixed_boundary$long
    cross_regimen_source_boundary_trace <- mixed_boundary$diagnostics
  } else {
    state_key <- .ltmle_exact_task_branch_key(source_task)
    rows <- .ltmle_exact_branch_state_rows(
      branch_state = state,
      task = source_task,
      spec = spec,
      node_spec = node_spec
    )
    long <- .ltmle_exact_branch_state_long(
      branch_state = state,
      task = source_task,
      spec = spec,
      node_spec = node_spec
    )
    rows$.source_task_id <- source_task_id
    long$.source_task_id <- source_task_id
    attr(rows, "spec") <- spec
    cross_regimen_source_boundary_trace <- .ltmle_exact_boundary_trace_row(
      consuming_task = current_task %||% task,
      source_task = source_task,
      metadata = source_boundary_metadata,
      rows = rows
    )
  }
  particles <- state[[state_key]]$particles
  boundary_state_keys <- unique(c(
    source_boundary_metadata$source_boundary_eval_state,
    source_boundary_metadata$source_boundary_outcome_history_state,
    source_boundary_metadata$source_boundary_m1_history_state,
    source_boundary_metadata$source_boundary_m2_history_state,
    source_boundary_metadata$source_boundary_auxiliary_mediator_history_state
  ))
  boundary_state_keys <- boundary_state_keys[
    !is.na(boundary_state_keys) &
      nzchar(boundary_state_keys) &
      boundary_state_keys != "not_applicable"
  ]
  safe_collapse_key <- paste(
    c(
      "id0",
      "source_task_id",
      "source_task_design_history",
      "state_key",
      "event_id",
      "branch_law_id",
      "target_mediator_identity",
      "auxiliary_mediator_identity",
      "branch_weight_ownership",
      "source_boundary_time_node",
      "treatment_censoring_state",
      "density_ratio_context"
    ),
    collapse = "|"
  )
  boundary_law_diag <- data.frame(
    component = source_task$component %||% NA_character_,
    task_id = source_task_id,
    source_task_id = source_task_id,
    t = as.integer(source_task$t %||% NA_integer_),
    node = source_task$node %||% NA_character_,
    process_type = .ltmle_exact_process_type(source_task),
    requested_density_ratio_mc_n = as.integer(ltmle_exact_density_ratio_mc_n),
    effective_density_ratio_mc_n = as.integer(ltmle_exact_density_ratio_mc_n),
    density_ratio_mc_cap_applied = FALSE,
    requested_law_integration_n = as.integer(ltmle_exact_law_integration_n),
    effective_law_integration_n = as.integer(ltmle_exact_law_integration_n),
    law_integration_cap_applied = FALSE,
    particle_pruning_applied = FALSE,
    max_generated_particles = as.integer(nrow(rows) * as.integer(ltmle_exact_density_ratio_mc_n)),
    integration_method = "source_boundary_materialization_no_terminal_recursion",
    source_eval_mode = "branch_state_dp_continuation_value",
    n_particles_before_integration = nrow(particles),
    n_particles_after_integration = nrow(particles),
    n_particles_after_collapse = length(unique(as.integer(particles$id0))),
    particle_count_before_node = nrow(particles),
    particle_count_after_node = nrow(particles),
    particle_count_after_collapse = length(unique(as.integer(particles$id0))),
    safe_collapse_applied = nrow(particles) != length(unique(as.integer(particles$id0))),
    safe_collapse_key = safe_collapse_key,
    run_mode = NA_character_,
    is_acceptance_gate = NA,
    passed = TRUE,
    failure_class = "no_failure",
    stringsAsFactors = FALSE
  )
  list(
    rows = rows,
    long = long,
    branch_state = state,
    source_task = source_task,
    state_key = state_key,
    provenance = source_boundary_metadata,
    weights = as.numeric(particles$.branch_weight %||% rep(1, nrow(particles))),
    state_keys = rep(state_key, nrow(rows)),
    collapse_keys = rep(safe_collapse_key, nrow(rows)),
    law_integration_diagnostics = .ltmle_exact_rbind_fill(c(transition_law_rows, list(boundary_law_diag))),
    handoff_trace = .ltmle_exact_rbind_fill(handoff_rows),
    cross_regimen_source_boundary_trace = cross_regimen_source_boundary_trace,
    diagnostics = list(
      branch_state_source_boundary_materialized = TRUE,
      source_boundary_reached = TRUE,
      source_boundary_task_id = source_task_id,
      source_boundary_type = source_boundary_metadata$source_boundary_type,
      source_boundary_direction = source_boundary_metadata$source_boundary_direction,
      source_boundary_eval_state = source_boundary_metadata$source_boundary_eval_state,
      source_boundary_outcome_history_state = source_boundary_metadata$source_boundary_outcome_history_state,
      source_boundary_m1_history_state = source_boundary_metadata$source_boundary_m1_history_state,
      source_boundary_m2_history_state = source_boundary_metadata$source_boundary_m2_history_state,
      source_boundary_auxiliary_mediator_history_state =
        source_boundary_metadata$source_boundary_auxiliary_mediator_history_state,
      qsource_predicted_on_mixed_boundary =
        .ltmle_exact_is_mixed_boundary_type(source_boundary_metadata$source_boundary_type),
      evaluated_beyond_source_boundary = FALSE,
      terminal_task_evaluated_in_source_path = FALSE,
      max_task_depth_beyond_source_boundary = 0L,
      n_source_boundary_rows = nrow(rows),
      n_source_boundary_state_keys = length(boundary_state_keys),
      n_particles_before_collapse = nrow(particles),
      n_particles_after_collapse = length(unique(as.integer(particles$id0))),
      safe_collapse_applied = nrow(particles) != length(unique(as.integer(particles$id0))),
      safe_collapse_key = safe_collapse_key,
      used_observed_long_as_source_history = FALSE,
      observed_long_used_as_generated_history = FALSE,
      observed_rows_mutated_to_assigned_A = FALSE,
      terminal_full_recursion_used = FALSE,
      branch_state_downstream_recursion_used = FALSE
    )
  )
}

.ltmle_exact_audit_integrate_source_boundary_transition <- function(state,
                                                                    current_task,
                                                                    source_task_id,
                                                                    task_graph,
                                                                    models,
                                                                    T,
                                                                    spec,
                                                                    node_spec,
                                                                    ltmle_exact_density_ratio_mc_n,
                                                                    ltmle_exact_law_integration_n) {
  force(T)
  current_task <- as.list(current_task)
  state_key <- .ltmle_exact_task_branch_key(current_task)
  tt <- as.integer(current_task$t)
  node <- as.character(current_task$node)
  process_type <- .ltmle_exact_process_type(current_task)
  if (identical(process_type, "virtual_mixed_continuation_task")) {
    virtual_step <- .ltmle_exact_integrate_virtual_mixed_mediator_path(
      branch_state = state,
      virtual_task = current_task,
      task_graph = task_graph,
      models = models,
      spec = spec,
      node_spec = node_spec,
      ltmle_exact_law_integration_n = ltmle_exact_law_integration_n,
      integration_context = "audit_virtual_mixed_source_boundary"
    )
    return(list(
      branch_state = virtual_step$branch_state,
      law_diag = virtual_step$law_integration_diagnostics %||% data.frame(),
      handoff_trace = virtual_step$handoff_trace %||% data.frame()
    ))
  }
  rows_current <- .ltmle_exact_branch_state_rows(
    branch_state = state,
    task = current_task,
    spec = spec,
    node_spec = node_spec
  )
  n0 <- nrow(rows_current)
  requested_n <- as.integer(ltmle_exact_law_integration_n)
  effective_n <- requested_n
  before_n <- nrow(state[[state_key]]$particles)

  if (node %in% .ltmle_exact_L_nodes(node_spec)) {
    fit <- models$L[[node]][[tt]]
    mu <- .ltmle_exact_predict_nuis(fit, rows_current, type = "numeric")
    sigma <- rep(.ltmle_exact_sigma(fit), n0)
    updater <- function(st, value, weight) {
      .ltmle_exact_branch_state_update_L(st, current_task, value, weight, spec, node_spec)
    }
  } else if (identical(node, "Y")) {
    fit <- models$Y[[tt]]
    if (is.null(fit)) .stop("Missing Y model for audit source-boundary materialization at t=", tt)
    mu <- .ltmle_exact_predict_nuis(fit, rows_current, type = "numeric")
    sigma <- rep(.ltmle_exact_sigma(fit), n0)
    updater <- function(st, value, weight) {
      .ltmle_exact_branch_state_update_Y(st, current_task, value, weight, spec, node_spec)
    }
  } else if (process_type %in% c(
    "observed_mediator_process",
    "joint_stochastic_mediator_intervention_law",
    "first_mediator_stochastic_intervention_law",
    "second_mediator_stochastic_intervention_law"
  ) && node %in% c("M1", "M2") && tt >= 2L) {
    if (identical(node, "M1")) {
      mu <- .ltmle_exact_predict_m1(models$M1[[tt]], rows_current)
      sigma <- .ltmle_exact_predict_mediator_sigma(models$M1[[tt]], rows_current)
    } else {
      mu <- .ltmle_exact_predict_m2(models$M2[[tt]], rows_current)
      sigma <- .ltmle_exact_predict_mediator_sigma(models$M2[[tt]], rows_current)
    }
    updater <- function(st, value, weight) {
      .ltmle_exact_branch_state_update_mediator(
        branch_state = st,
        task = current_task,
        node = node,
        value = value,
        weight_multiplier = weight,
        spec = spec,
        node_spec = node_spec
      )
    }
  } else {
    return(list(
      branch_state = state,
      law_diag = data.frame(),
      handoff_trace = data.frame()
    ))
  }

  grid <- .ltmle_exact_normal_quantile_grid(effective_n)
  ng <- length(grid$z)
  row_idx <- rep(seq_len(n0), times = ng)
  state_g <- .ltmle_exact_branch_state_expand(state, row_idx)
  values <- rep(mu, times = ng) + rep(sigma, times = ng) * rep(grid$z, each = n0)
  weights <- rep(grid$w, each = n0)
  state_g <- updater(state_g, values, weights)
  after_n <- nrow(state_g[[state_key]]$particles)
  collapse_n <- length(unique(as.integer(state_g[[state_key]]$particles$id0)))
  law_diag <- data.frame(
    component = current_task$component %||% NA_character_,
    task_id = current_task$task_id %||% NA_character_,
    source_task_id = source_task_id,
    t = tt,
    node = node,
    process_type = process_type,
    requested_density_ratio_mc_n = as.integer(ltmle_exact_density_ratio_mc_n),
    effective_density_ratio_mc_n = as.integer(ltmle_exact_density_ratio_mc_n),
    density_ratio_mc_cap_applied = FALSE,
    requested_law_integration_n = requested_n,
    effective_law_integration_n = effective_n,
    law_integration_cap_applied = FALSE,
    particle_pruning_applied = FALSE,
    max_generated_particles = as.integer(after_n * as.integer(ltmle_exact_density_ratio_mc_n)),
    integration_method = "audit_source_boundary_transition_grid",
    source_eval_mode = "independent_branch_state_dp_continuation_value",
    n_particles_before_integration = before_n,
    n_particles_after_integration = after_n,
    n_particles_after_collapse = collapse_n,
    particle_count_before_node = before_n,
    particle_count_after_node = after_n,
    particle_count_after_collapse = collapse_n,
    safe_collapse_applied = after_n != collapse_n,
    safe_collapse_key = "id0|state_key|source_task_id|event_id|branch_weight_ownership|density_ratio_context",
    run_mode = NA_character_,
    is_acceptance_gate = NA,
    passed = TRUE,
    failure_class = "no_failure",
    stringsAsFactors = FALSE
  )
  list(
    branch_state = state_g,
    law_diag = law_diag,
    handoff_trace = attr(state_g, "last_handoff") %||% data.frame()
  )
}

.ltmle_exact_audit_materialize_source_task_boundary <- function(branch_state,
                                                                task = NULL,
                                                                source_task_id,
                                                                task_graph,
                                                                models,
                                                                treatment_models,
                                                                censoring_models,
                                                                T,
                                                                spec,
                                                                probability_bounds,
                                                                treat_mech,
                                                                p_rct,
                                                                node_spec = NULL,
                                                                ltmle_exact_density_ratio_mc_n = 2000L,
                                                                ltmle_exact_law_integration_n = 5L,
                                                                forbid_terminal_full_recursion = TRUE) {
  force(treatment_models)
  force(censoring_models)
  force(probability_bounds)
  force(treat_mech)
  force(p_rct)
  if (!isTRUE(forbid_terminal_full_recursion)) {
    .stop("terminal_full_recursion_used_in_independent_source_boundary_path")
  }
  if (.ltmle_exact_missing_or_empty(source_task_id)) {
    .stop("missing_source_task_id")
  }
  if (identical(source_task_id, "observed_terminal_outcome")) {
    .stop("Terminal source boundary must use the explicit terminal base-case.")
  }
  if (is.null(task_graph$tasks) || !source_task_id %in% names(task_graph$tasks)) {
    .stop("source_task_id_not_in_task_graph: ", source_task_id)
  }

  transition_law_rows <- list()
  handoff_rows <- list()
  state <- branch_state
  current_task <- if (is.null(task)) NULL else as.list(task)
  visited <- character(0)
  while (!is.null(current_task) &&
         !identical(as.character(current_task$task_id), as.character(source_task_id))) {
    current_id <- as.character(current_task$task_id)
    if (current_id %in% visited) .stop("audit_branch_state_recursion_cycle_detected")
    visited <- c(visited, current_id)
    child_id <- as.character(current_task$observed_pseudooutcome_source_task_id %||% NA_character_)
    if (.ltmle_exact_missing_or_empty(child_id)) .stop("missing_source_task_id")
    if (identical(child_id, "observed_terminal_outcome")) {
      .stop("terminal_full_recursion_used_in_independent_source_boundary_path")
    }

    step <- .ltmle_exact_audit_integrate_source_boundary_transition(
      state = state,
      current_task = current_task,
      source_task_id = source_task_id,
      task_graph = task_graph,
      models = models,
      T = T,
      spec = spec,
      node_spec = node_spec,
      ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
      ltmle_exact_law_integration_n = ltmle_exact_law_integration_n
    )
    state <- step$branch_state
    if (!is.null(step$law_diag) && nrow(step$law_diag)) {
      transition_law_rows[[length(transition_law_rows) + 1L]] <- step$law_diag
    }
    if (!is.null(step$handoff_trace) && nrow(step$handoff_trace)) {
      handoff_rows[[length(handoff_rows) + 1L]] <- step$handoff_trace
    }
    if (identical(child_id, source_task_id)) break
    if (is.null(task_graph$tasks[[child_id]])) {
      .stop("source_task_id_not_in_task_graph: ", child_id)
    }
    current_task <- as.list(task_graph$tasks[[child_id]])
  }

  source_task <- as.list(task_graph$tasks[[source_task_id]])
  source_boundary_metadata <- .ltmle_exact_source_boundary_metadata(
    consuming_task = current_task %||% task,
    source_task = source_task,
    spec = spec
  )
  if (.ltmle_exact_is_mixed_boundary_type(source_boundary_metadata$source_boundary_type)) {
    mixed_boundary <- .ltmle_exact_materialize_cross_regimen_source_boundary(
      branch_state = state,
      consuming_task = current_task %||% task,
      source_task = source_task,
      source_boundary_metadata = source_boundary_metadata,
      spec = spec,
      node_spec = node_spec
    )
    state_key <- mixed_boundary$state_key
    rows <- mixed_boundary$rows
    long <- mixed_boundary$long
    cross_regimen_source_boundary_trace <- mixed_boundary$diagnostics
  } else {
    state_key <- .ltmle_exact_task_branch_key(source_task)
    rows <- .ltmle_exact_branch_state_rows(
      branch_state = state,
      task = source_task,
      spec = spec,
      node_spec = node_spec
    )
    long <- .ltmle_exact_branch_state_long(
      branch_state = state,
      task = source_task,
      spec = spec,
      node_spec = node_spec
    )
    rows$.source_task_id <- source_task_id
    long$.source_task_id <- source_task_id
    attr(rows, "spec") <- spec
    cross_regimen_source_boundary_trace <- .ltmle_exact_boundary_trace_row(
      consuming_task = current_task %||% task,
      source_task = source_task,
      metadata = source_boundary_metadata,
      rows = rows
    )
  }
  particles <- state[[state_key]]$particles
  boundary_state_keys <- unique(c(
    source_boundary_metadata$source_boundary_eval_state,
    source_boundary_metadata$source_boundary_outcome_history_state,
    source_boundary_metadata$source_boundary_m1_history_state,
    source_boundary_metadata$source_boundary_m2_history_state,
    source_boundary_metadata$source_boundary_auxiliary_mediator_history_state
  ))
  boundary_state_keys <- boundary_state_keys[
    !is.na(boundary_state_keys) &
      nzchar(boundary_state_keys) &
      boundary_state_keys != "not_applicable"
  ]
  safe_collapse_key <- paste(
    c(
      "id0",
      "source_task_id",
      "source_task_design_history",
      "state_key",
      "event_id",
      "branch_law_id",
      "target_mediator_identity",
      "auxiliary_mediator_identity",
      "branch_weight_ownership",
      "source_boundary_time_node",
      "treatment_censoring_state",
      "density_ratio_context"
    ),
    collapse = "|"
  )
  boundary_law_diag <- data.frame(
    component = source_task$component %||% NA_character_,
    task_id = source_task_id,
    source_task_id = source_task_id,
    t = as.integer(source_task$t %||% NA_integer_),
    node = source_task$node %||% NA_character_,
    process_type = .ltmle_exact_process_type(source_task),
    requested_density_ratio_mc_n = as.integer(ltmle_exact_density_ratio_mc_n),
    effective_density_ratio_mc_n = as.integer(ltmle_exact_density_ratio_mc_n),
    density_ratio_mc_cap_applied = FALSE,
    requested_law_integration_n = as.integer(ltmle_exact_law_integration_n),
    effective_law_integration_n = as.integer(ltmle_exact_law_integration_n),
    law_integration_cap_applied = FALSE,
    particle_pruning_applied = FALSE,
    max_generated_particles = as.integer(nrow(rows) * as.integer(ltmle_exact_density_ratio_mc_n)),
    integration_method = "audit_source_boundary_materialization_no_terminal_recursion",
    source_eval_mode = "independent_branch_state_dp_continuation_value",
    n_particles_before_integration = nrow(particles),
    n_particles_after_integration = nrow(particles),
    n_particles_after_collapse = length(unique(as.integer(particles$id0))),
    particle_count_before_node = nrow(particles),
    particle_count_after_node = nrow(particles),
    particle_count_after_collapse = length(unique(as.integer(particles$id0))),
    safe_collapse_applied = nrow(particles) != length(unique(as.integer(particles$id0))),
    safe_collapse_key = safe_collapse_key,
    run_mode = NA_character_,
    is_acceptance_gate = NA,
    passed = TRUE,
    failure_class = "no_failure",
    stringsAsFactors = FALSE
  )
  list(
    rows = rows,
    long = long,
    branch_state = state,
    source_task = source_task,
    state_key = state_key,
    provenance = source_boundary_metadata,
    weights = as.numeric(particles$.branch_weight %||% rep(1, nrow(particles))),
    state_keys = rep(state_key, nrow(rows)),
    collapse_keys = rep(safe_collapse_key, nrow(rows)),
    law_integration_diagnostics = .ltmle_exact_rbind_fill(c(transition_law_rows, list(boundary_law_diag))),
    handoff_trace = .ltmle_exact_rbind_fill(handoff_rows),
    cross_regimen_source_boundary_trace = cross_regimen_source_boundary_trace,
    diagnostics = list(
      branch_state_source_boundary_materialized = TRUE,
      source_boundary_reached = TRUE,
      source_boundary_task_id = source_task_id,
      source_boundary_type = source_boundary_metadata$source_boundary_type,
      source_boundary_direction = source_boundary_metadata$source_boundary_direction,
      source_boundary_eval_state = source_boundary_metadata$source_boundary_eval_state,
      source_boundary_outcome_history_state = source_boundary_metadata$source_boundary_outcome_history_state,
      source_boundary_m1_history_state = source_boundary_metadata$source_boundary_m1_history_state,
      source_boundary_m2_history_state = source_boundary_metadata$source_boundary_m2_history_state,
      source_boundary_auxiliary_mediator_history_state =
        source_boundary_metadata$source_boundary_auxiliary_mediator_history_state,
      qsource_predicted_on_mixed_boundary =
        .ltmle_exact_is_mixed_boundary_type(source_boundary_metadata$source_boundary_type),
      evaluated_beyond_source_boundary = FALSE,
      terminal_task_evaluated_in_source_path = FALSE,
      max_task_depth_beyond_source_boundary = 0L,
      n_source_boundary_rows = nrow(rows),
      n_source_boundary_state_keys = length(boundary_state_keys),
      n_particles_before_collapse = nrow(particles),
      n_particles_after_collapse = length(unique(as.integer(particles$id0))),
      safe_collapse_applied = nrow(particles) != length(unique(as.integer(particles$id0))),
      safe_collapse_key = safe_collapse_key,
      used_observed_long_as_source_history = FALSE,
      observed_long_used_as_generated_history = FALSE,
      observed_rows_mutated_to_assigned_A = FALSE,
      terminal_full_recursion_used = FALSE,
      branch_state_downstream_recursion_used = FALSE,
      independent_source_boundary_materializer_completed = TRUE,
      independent_source_boundary_materializer_helper_name =
        ".ltmle_exact_audit_materialize_source_task_boundary",
      same_boundary_materializer_recomputation_check = FALSE
    )
  )
}

.ltmle_exact_collapse_source_qstar_to_subjects <- function(q_source_star,
                                                           source_boundary,
                                                           task,
                                                           source_task_id) {
  force(task)
  q_source_star <- as.numeric(q_source_star)
  state_key <- source_boundary$state_key
  branch_state <- source_boundary$branch_state
  value <- .ltmle_exact_collapse_particles_to_subjects(q_source_star, branch_state, state_key)
  if (any(!is.finite(value))) {
    .stop("Non-finite collapsed source Q* values for source_task_id=", source_task_id)
  }
  list(
    value = value,
    id = seq_along(value),
    state_key = rep(state_key, length(value)),
    particle_q = q_source_star,
    subject_q = value,
    particle_frame = branch_state[[state_key]]$particles,
    branch_trace = data.frame(),
    handoff_trace = source_boundary$handoff_trace %||% data.frame(),
    diagnostics = list(
      density_ratio_diagnostics = attr(q_source_star, "density_ratio_diagnostics") %||% data.frame(),
      law_integration_diagnostics = source_boundary$law_integration_diagnostics %||% data.frame(),
      cross_regimen_source_boundary_trace =
        source_boundary$cross_regimen_source_boundary_trace %||% data.frame(),
      H_new_mean = attr(q_source_star, "H_new_mean") %||% NA_real_,
      H_new_max = attr(q_source_star, "H_new_max") %||% NA_real_
    ),
    metadata = c(
      source_boundary$diagnostics,
      list(
        source_eval_mode = "branch_state_dp_continuation_value",
        cached_targeted_continuation_used = TRUE,
        continuation_task_id = source_task_id,
        continuation_fit_found = TRUE,
        continuation_fit_targeted = TRUE,
        terminal_outcome_base_case_used = FALSE,
        terminal_full_recursion_used = FALSE,
        local_qstar_source_prediction_used = FALSE,
        local_current_task_prediction_used = FALSE,
        actual_helper_name = ".ltmle_exact_eval_source_pseudooutcome_dp"
      )
    )
  )
}

.ltmle_exact_virtual_mixed_outer_L_source_task_id <- function(virtual_task,
                                                             task_graph) {
  virtual_task <- as.list(virtual_task)
  if (!.ltmle_exact_is_virtual_mixed_task(virtual_task)) {
    .stop("virtual_mixed_outer_L_source_requires_virtual_task")
  }
  if (is.null(task_graph$tasks) || !length(task_graph$tasks)) {
    .stop("virtual_mixed_outer_L_source_missing_task_graph")
  }
  component <- as.character(virtual_task$component %||% NA_character_)
  tt <- as.integer(virtual_task$t %||% NA_integer_)
  node <- as.character(virtual_task$node %||% NA_character_)
  ids <- names(task_graph$tasks)[vapply(task_graph$tasks, function(x) {
    xx <- as.list(x)
    identical(as.character(xx$component %||% NA_character_), component) &&
      identical(as.integer(xx$t %||% NA_integer_), tt) &&
      identical(as.character(xx$node %||% NA_character_), node) &&
      identical(
        as.character(.ltmle_exact_process_type(xx)),
        "post_mediator_covariate_transition"
      )
  }, logical(1))]
  if (length(ids) != 1L) {
    .stop("virtual_mixed_outer_L_source_not_unique: task=",
          as.character(virtual_task$task_id %||% NA_character_),
          ", n_matches=", length(ids))
  }
  ids[1L]
}

.ltmle_exact_eval_virtual_mixed_dedicated_L_continuation <- function(
  virtual_task,
  branch_state,
  task_graph,
  task_fit_cache,
  models,
  treatment_models,
  censoring_models,
  T,
  spec,
  probability_bounds,
  treat_mech,
  p_rct,
  node_spec = NULL,
  ltmle_exact_density_ratio_mc_n = 2000L,
  ltmle_exact_law_integration_n = 5L,
  actual_helper_name = ".ltmle_exact_eval_source_pseudooutcome_dp",
  integration_context = "virtual_mixed_dedicated_L_transition_target"
) {
  virtual_task <- as.list(virtual_task)
  if (!.ltmle_exact_is_virtual_mixed_task(virtual_task)) {
    .stop("dedicated_virtual_mixed_target_requires_virtual_task")
  }
  direct_status <- .ltmle_exact_virtual_mixed_direct_continuation_status(virtual_task)
  if (isTRUE(direct_status$downstream_source_can_represent_virtual_target)) {
    .stop("dedicated_virtual_mixed_target_requested_for_representable_downstream_source: ",
          as.character(virtual_task$task_id %||% NA_character_))
  }
  source_task_id <- .ltmle_exact_virtual_mixed_outer_L_source_task_id(
    virtual_task = virtual_task,
    task_graph = task_graph
  )
  source_fit <- task_fit_cache[[source_task_id]]
  if (is.null(source_fit)) {
    .stop("dedicated_virtual_mixed_L_transition_fit_missing: ", source_task_id)
  }
  source_task <- as.list(task_graph$tasks[[source_task_id]])
  virtual_step <- .ltmle_exact_integrate_virtual_mixed_mediator_path(
    branch_state = branch_state,
    virtual_task = virtual_task,
    task_graph = task_graph,
    models = models,
    spec = spec,
    node_spec = node_spec,
    ltmle_exact_law_integration_n = ltmle_exact_law_integration_n,
    integration_context = integration_context
  )
  metadata <- .ltmle_exact_task_evaluation_metadata(virtual_task, spec)
  rows <- .ltmle_exact_materialize_mixed_boundary_rows_at(
    branch_state = virtual_step$branch_state,
    metadata = metadata,
    tt = as.integer(source_task$t),
    node_spec = node_spec
  )
  rows$A <- as.numeric(.ltmle_exact_task_assigned_A(source_task, spec))
  rows$.source_task_id <- source_task_id
  attr(rows, "spec") <- spec
  long <- .ltmle_exact_materialize_mixed_boundary_long(
    branch_state = virtual_step$branch_state,
    metadata = metadata,
    max_t = as.integer(source_task$t),
    node_spec = node_spec
  )
  long$.source_task_id <- source_task_id
  attr(long, "spec") <- spec
	  q_source_star <- .ltmle_exact_predict_task_qstar_on_rows(
	    task_fit = source_fit,
	    rows = rows,
	    row_long = long,
    models = models,
    treatment_models = treatment_models,
    censoring_models = censoring_models,
    T = T,
    spec = spec,
    probability_bounds = probability_bounds,
	    treat_mech = treat_mech,
	    p_rct = p_rct,
	    node_spec = node_spec,
	    ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n
	  )
  q_source_initial <- .ltmle_exact_apply_fluctuation(
    q0_new = .ltmle_exact_predict_continuation(source_fit$cont_fit, rows),
    H_new = rep(0, nrow(rows)),
    epsilon = 0,
    bounds = source_fit$bounds
  )
  source_boundary <- list(
    rows = rows,
    long = long,
    branch_state = virtual_step$branch_state,
    source_task = source_task,
    state_key = metadata$source_boundary_eval_state,
    provenance = metadata,
    law_integration_diagnostics = virtual_step$law_integration_diagnostics %||% data.frame(),
    handoff_trace = virtual_step$handoff_trace %||% data.frame(),
    cross_regimen_source_boundary_trace = .ltmle_exact_boundary_trace_row(
      consuming_task = virtual_task,
      source_task = source_task,
      metadata = metadata,
      rows = rows,
      source_fit = source_fit,
      source_covariates = source_fit$conditioning_covariates %||% character(0),
      qsource_predicted = TRUE,
      cached_targeted_continuation_used = TRUE
    ),
    diagnostics = list(
      branch_state_source_boundary_materialized = TRUE,
      source_boundary_reached = TRUE,
      source_boundary_task_id = source_task_id,
      source_boundary_type = metadata$source_boundary_type,
      source_boundary_direction = metadata$source_boundary_direction,
      source_boundary_eval_state = metadata$source_boundary_eval_state,
      source_boundary_outcome_history_state = metadata$source_boundary_outcome_history_state,
      source_boundary_m1_history_state = metadata$source_boundary_m1_history_state,
      source_boundary_m2_history_state = metadata$source_boundary_m2_history_state,
      source_boundary_auxiliary_mediator_history_state =
        metadata$source_boundary_auxiliary_mediator_history_state,
      qsource_predicted_on_mixed_boundary = TRUE,
      evaluated_beyond_source_boundary = FALSE,
      terminal_task_evaluated_in_source_path = FALSE,
      max_task_depth_beyond_source_boundary = 0L,
      n_source_boundary_rows = nrow(rows),
      n_source_boundary_state_keys = length(unique(c(
        metadata$source_boundary_eval_state,
        metadata$source_boundary_outcome_history_state,
        metadata$source_boundary_m1_history_state,
        metadata$source_boundary_m2_history_state,
        metadata$source_boundary_auxiliary_mediator_history_state
      ))),
      n_particles_before_collapse = nrow(virtual_step$branch_state[[metadata$source_boundary_eval_state]]$particles),
      n_particles_after_collapse = length(unique(as.integer(
        virtual_step$branch_state[[metadata$source_boundary_eval_state]]$particles$id0
      ))),
      safe_collapse_applied = nrow(virtual_step$branch_state[[metadata$source_boundary_eval_state]]$particles) !=
        length(unique(as.integer(virtual_step$branch_state[[metadata$source_boundary_eval_state]]$particles$id0))),
      safe_collapse_key =
        "id0|source_task_id|source_task_design_history|state_key|event_id|branch_law_id|target_mediator_identity|auxiliary_mediator_identity|branch_weight_ownership|source_boundary_time_node|treatment_censoring_state|density_ratio_context",
      used_observed_long_as_source_history = FALSE,
      observed_long_used_as_generated_history = FALSE,
      observed_rows_mutated_to_assigned_A = FALSE,
      terminal_full_recursion_used = FALSE,
      branch_state_downstream_recursion_used = FALSE
    )
  )
  initial_value <- .ltmle_exact_collapse_particles_to_subjects(
    q_source_initial,
    source_boundary$branch_state,
    source_boundary$state_key
  )
  dedicated_L_prediction_row_trace <- data.frame()
  dedicated_L_coefficient_contribution_check <- data.frame()
  dedicated_L_prediction_weighting_collapse_check <- data.frame()
  dedicated_L_mediator_path_weight_trace <- data.frame()
  if (.ltmle_exact_dedicated_L_transition_instrumentation_enabled(
    virtual_task = virtual_task,
    source_task = source_task
  )) {
    dedicated_L_prediction_row_trace <-
      .ltmle_exact_dedicated_L_prediction_row_trace(
        virtual_task = virtual_task,
        source_task = source_task,
        source_fit = source_fit,
        rows = rows,
        q_initial = q_source_initial,
        metadata = metadata
      )
    dedicated_L_coefficient_contribution_check <-
      .ltmle_exact_dedicated_L_coefficient_contribution_check(
        virtual_task = virtual_task,
        source_task = source_task,
        source_fit = source_fit,
        rows = rows
      )
    dedicated_L_prediction_weighting_collapse_check <-
      .ltmle_exact_dedicated_L_prediction_weighting_collapse_check(
        virtual_task = virtual_task,
        source_task = source_task,
        source_boundary = source_boundary,
        q_initial = q_source_initial,
        subject_initial = initial_value
      )
    dedicated_L_mediator_path_weight_trace <-
      .ltmle_exact_dedicated_L_mediator_path_weight_trace(
        virtual_task = virtual_task,
        source_task = source_task,
        source_boundary = source_boundary,
        q_initial = q_source_initial
      )
  }
  out <- .ltmle_exact_collapse_source_qstar_to_subjects(
    q_source_star = q_source_star,
    source_boundary = source_boundary,
    task = virtual_task,
    source_task_id = source_task_id
  )
  out$initial_value <- initial_value
  out$diagnostics$dedicated_L_transition_L_fit_prediction_row_trace <-
    dedicated_L_prediction_row_trace
  out$diagnostics$dedicated_L_transition_L_fit_coefficient_contribution_check <-
    dedicated_L_coefficient_contribution_check
  out$diagnostics$dedicated_L_transition_prediction_weighting_collapse_check <-
    dedicated_L_prediction_weighting_collapse_check
  out$diagnostics$dedicated_L_transition_mediator_path_weight_trace <-
    dedicated_L_mediator_path_weight_trace
  out$diagnostics$cross_regimen_source_boundary_trace <-
    .ltmle_exact_rbind_fill(list(
      source_boundary$cross_regimen_source_boundary_trace %||% data.frame(),
      out$diagnostics$cross_regimen_source_boundary_trace %||% data.frame()
    ))
  out$diagnostics$law_integration_diagnostics <-
    .ltmle_exact_rbind_fill(list(
      source_boundary$law_integration_diagnostics %||% data.frame(),
      out$diagnostics$law_integration_diagnostics %||% data.frame()
    ))
  downstream_source_id <- as.character(
    virtual_task$observed_pseudooutcome_source_task_id %||% NA_character_
  )
  downstream_fit <- if (!.ltmle_exact_missing_or_empty(downstream_source_id)) {
    task_fit_cache[[downstream_source_id]]
  } else {
    NULL
  }
  out$metadata <- c(
    source_boundary$diagnostics,
    out$metadata %||% list(),
    list(
      source_eval_mode = "branch_state_dp_continuation_value",
      cached_targeted_continuation_used = TRUE,
      continuation_task_id = source_task_id,
      continuation_fit_found = !is.null(source_fit),
      continuation_fit_targeted = isTRUE(source_fit$targeted),
      terminal_outcome_base_case_used = FALSE,
      terminal_full_recursion_used = FALSE,
      local_qstar_source_prediction_used = FALSE,
      local_current_task_prediction_used = FALSE,
      actual_helper_name = actual_helper_name,
      virtual_mixed_direct_continuation_used = FALSE,
      direct_continuation_allowed = isTRUE(direct_status$direct_continuation_allowed),
      downstream_source_can_represent_virtual_target =
        isTRUE(direct_status$downstream_source_can_represent_virtual_target),
      dedicated_virtual_Q_target_used = TRUE,
      virtual_mixed_downstream_source_task_id = downstream_source_id,
      virtual_mixed_downstream_source_fit_found = !is.null(downstream_fit),
      virtual_mixed_downstream_source_fit_targeted = isTRUE(downstream_fit$targeted),
      dedicated_virtual_Q_source_task_id = source_task_id,
      dedicated_virtual_Q_source_process_type = .ltmle_exact_process_type(source_task),
      dedicated_virtual_Q_source_node = as.character(source_task$node %||% NA_character_),
      dedicated_virtual_Q_source_t = as.integer(source_task$t %||% NA_integer_),
      initial_dp_source_eval_available = length(initial_value) == length(out$value) &&
        all(is.finite(initial_value)),
      initial_dp_source_eval_mean = mean(initial_value, na.rm = TRUE),
      initial_dp_source_eval_source =
        "same_time_outer_L_transition_fit_epsilon_zero_on_component_specific_mixed_mediator_rows",
      targeted_dp_source_eval_mean = mean(out$value, na.rm = TRUE),
      targeted_minus_initial_dp_source_eval =
        mean(out$value, na.rm = TRUE) - mean(initial_value, na.rm = TRUE)
    )
  )
  out
}

.ltmle_exact_eval_virtual_mixed_task_fit_on_mixed_rows <- function(
  virtual_task,
  consuming_task,
  branch_state,
  task_graph,
  task_fit_cache,
  models,
  treatment_models,
  censoring_models,
  T,
  spec,
  probability_bounds,
  treat_mech,
  p_rct,
  node_spec = NULL,
  ltmle_exact_density_ratio_mc_n = 2000L,
  ltmle_exact_law_integration_n = 5L,
  actual_helper_name = ".ltmle_exact_eval_source_pseudooutcome_dp",
  integration_context = "virtual_mixed_task_fit_source_boundary"
) {
  virtual_task <- as.list(virtual_task)
  consuming_task <- as.list(consuming_task)
  if (!.ltmle_exact_is_virtual_mixed_task(virtual_task)) {
    .stop("virtual_mixed_task_fit_evaluation_requires_virtual_task")
  }
  source_task_id <- as.character(virtual_task$task_id %||% NA_character_)
  source_fit <- task_fit_cache[[source_task_id]]
  if (is.null(source_fit)) {
    .stop("virtual_mixed_task_fit_missing: ", source_task_id)
  }
  direct_status <- .ltmle_exact_virtual_mixed_direct_continuation_status(virtual_task)
  virtual_step <- .ltmle_exact_integrate_virtual_mixed_mediator_path(
    branch_state = branch_state,
    virtual_task = virtual_task,
    task_graph = task_graph,
    models = models,
    spec = spec,
    node_spec = node_spec,
    ltmle_exact_law_integration_n = ltmle_exact_law_integration_n,
    integration_context = integration_context
  )
  metadata <- .ltmle_exact_task_evaluation_metadata(virtual_task, spec)
  rows <- .ltmle_exact_materialize_mixed_boundary_rows_at(
    branch_state = virtual_step$branch_state,
    metadata = metadata,
    tt = as.integer(virtual_task$t),
    node_spec = node_spec
  )
  rows$A <- as.numeric(.ltmle_exact_task_assigned_A(virtual_task, spec))
  rows$.source_task_id <- source_task_id
  attr(rows, "spec") <- spec
  long <- .ltmle_exact_materialize_mixed_boundary_long(
    branch_state = virtual_step$branch_state,
    metadata = metadata,
    max_t = as.integer(virtual_task$t),
    node_spec = node_spec
  )
  long$A[as.integer(long$t) == as.integer(virtual_task$t)] <-
    as.numeric(.ltmle_exact_task_assigned_A(virtual_task, spec))
  long$.source_task_id <- source_task_id
  attr(long, "spec") <- spec
  q_source_star <- .ltmle_exact_predict_task_qstar_on_rows(
    task_fit = source_fit,
    rows = rows,
    row_long = long,
    models = models,
    treatment_models = treatment_models,
    censoring_models = censoring_models,
    T = T,
    spec = spec,
    probability_bounds = probability_bounds,
    treat_mech = treat_mech,
    p_rct = p_rct,
    node_spec = node_spec,
    ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n
  )
  source_boundary <- list(
    rows = rows,
    long = long,
    branch_state = virtual_step$branch_state,
    source_task = virtual_task,
    state_key = metadata$source_boundary_eval_state,
    provenance = metadata,
    law_integration_diagnostics = virtual_step$law_integration_diagnostics %||% data.frame(),
    handoff_trace = virtual_step$handoff_trace %||% data.frame(),
    cross_regimen_source_boundary_trace = .ltmle_exact_boundary_trace_row(
      consuming_task = consuming_task,
      source_task = virtual_task,
      metadata = metadata,
      rows = rows,
      source_fit = source_fit,
      source_covariates = source_fit$conditioning_covariates %||% character(0),
      qsource_predicted = TRUE,
      cached_targeted_continuation_used = TRUE
    ),
    diagnostics = list(
      branch_state_source_boundary_materialized = TRUE,
      source_boundary_reached = TRUE,
      source_boundary_task_id = source_task_id,
      source_boundary_type = metadata$source_boundary_type,
      source_boundary_direction = metadata$source_boundary_direction,
      source_boundary_eval_state = metadata$source_boundary_eval_state,
      source_boundary_outcome_history_state = metadata$source_boundary_outcome_history_state,
      source_boundary_m1_history_state = metadata$source_boundary_m1_history_state,
      source_boundary_m2_history_state = metadata$source_boundary_m2_history_state,
      source_boundary_auxiliary_mediator_history_state =
        metadata$source_boundary_auxiliary_mediator_history_state,
      qsource_predicted_on_mixed_boundary = TRUE,
      evaluated_beyond_source_boundary = FALSE,
      terminal_task_evaluated_in_source_path = FALSE,
      max_task_depth_beyond_source_boundary = 0L,
      n_source_boundary_rows = nrow(rows),
      n_source_boundary_state_keys = length(unique(c(
        metadata$source_boundary_eval_state,
        metadata$source_boundary_outcome_history_state,
        metadata$source_boundary_m1_history_state,
        metadata$source_boundary_m2_history_state,
        metadata$source_boundary_auxiliary_mediator_history_state
      ))),
      n_particles_before_collapse =
        nrow(virtual_step$branch_state[[metadata$source_boundary_eval_state]]$particles),
      n_particles_after_collapse = length(unique(as.integer(
        virtual_step$branch_state[[metadata$source_boundary_eval_state]]$particles$id0
      ))),
      safe_collapse_applied =
        nrow(virtual_step$branch_state[[metadata$source_boundary_eval_state]]$particles) !=
        length(unique(as.integer(
          virtual_step$branch_state[[metadata$source_boundary_eval_state]]$particles$id0
        ))),
      safe_collapse_key =
        "id0|source_task_id|source_task_design_history|state_key|event_id|branch_law_id|target_mediator_identity|auxiliary_mediator_identity|branch_weight_ownership|source_boundary_time_node|treatment_censoring_state|density_ratio_context",
      used_observed_long_as_source_history = FALSE,
      observed_long_used_as_generated_history = FALSE,
      observed_rows_mutated_to_assigned_A = FALSE,
      terminal_full_recursion_used = FALSE,
      branch_state_downstream_recursion_used = FALSE
    )
  )
  out <- .ltmle_exact_collapse_source_qstar_to_subjects(
    q_source_star = q_source_star,
    source_boundary = source_boundary,
    task = consuming_task,
    source_task_id = source_task_id
  )
  virtual_downstream_id <- as.character(
    virtual_task$observed_pseudooutcome_source_task_id %||% NA_character_
  )
  virtual_downstream_fit <- if (!.ltmle_exact_missing_or_empty(virtual_downstream_id)) {
    task_fit_cache[[virtual_downstream_id]]
  } else {
    NULL
  }
  out$metadata <- c(
    source_boundary$diagnostics,
    out$metadata %||% list(),
    list(
      source_eval_mode = "branch_state_dp_continuation_value",
      cached_targeted_continuation_used = TRUE,
      continuation_task_id = source_task_id,
      continuation_fit_found = !is.null(source_fit),
      continuation_fit_targeted = isTRUE(source_fit$targeted),
      terminal_outcome_base_case_used = FALSE,
      terminal_full_recursion_used = FALSE,
      local_qstar_source_prediction_used = FALSE,
      local_current_task_prediction_used = FALSE,
      actual_helper_name = actual_helper_name,
      virtual_mixed_direct_continuation_used = FALSE,
      direct_continuation_allowed = isTRUE(direct_status$direct_continuation_allowed),
      downstream_source_can_represent_virtual_target =
        isTRUE(direct_status$downstream_source_can_represent_virtual_target),
      dedicated_virtual_Q_target_used =
        isTRUE(direct_status$dedicated_virtual_Q_target_used),
      virtual_mixed_downstream_source_task_id = virtual_downstream_id,
      virtual_mixed_downstream_source_fit_found = !is.null(virtual_downstream_fit),
      virtual_mixed_downstream_source_fit_targeted =
        isTRUE(virtual_downstream_fit$targeted)
    )
  )
  out$diagnostics$cross_regimen_source_boundary_trace <-
    .ltmle_exact_rbind_fill(list(
      source_boundary$cross_regimen_source_boundary_trace %||% data.frame(),
      out$diagnostics$cross_regimen_source_boundary_trace %||% data.frame()
    ))
  out$diagnostics$law_integration_diagnostics <-
    .ltmle_exact_rbind_fill(list(
      source_boundary$law_integration_diagnostics %||% data.frame(),
      out$diagnostics$law_integration_diagnostics %||% data.frame()
    ))
  out$handoff_trace <- .ltmle_exact_rbind_fill(list(
    source_boundary$handoff_trace %||% data.frame(),
    out$handoff_trace %||% data.frame()
  ))
  out
}

.ltmle_exact_eval_virtual_mixed_task_fit_on_source_boundary <- function(
  virtual_task,
  consuming_task,
  source_boundary,
  task_fit_cache,
  models,
  treatment_models,
  censoring_models,
  T,
  spec,
  probability_bounds,
  treat_mech,
  p_rct,
  node_spec = NULL,
  ltmle_exact_density_ratio_mc_n = 2000L,
  actual_helper_name = ".ltmle_exact_eval_source_pseudooutcome_dp"
) {
  virtual_task <- as.list(virtual_task)
  consuming_task <- as.list(consuming_task)
  if (!.ltmle_exact_is_virtual_mixed_task(virtual_task)) {
    .stop("virtual_mixed_task_fit_source_boundary_requires_virtual_task")
  }
  source_task_id <- as.character(virtual_task$task_id %||% NA_character_)
  source_fit <- task_fit_cache[[source_task_id]]
  if (is.null(source_fit)) {
    .stop("virtual_mixed_task_fit_missing: ", source_task_id)
  }
  direct_status <- .ltmle_exact_virtual_mixed_direct_continuation_status(virtual_task)
  q_source_star <- .ltmle_exact_predict_task_qstar_on_rows(
    task_fit = source_fit,
    rows = source_boundary$rows,
    row_long = source_boundary$long,
    models = models,
    treatment_models = treatment_models,
    censoring_models = censoring_models,
    T = T,
    spec = spec,
    probability_bounds = probability_bounds,
    treat_mech = treat_mech,
    p_rct = p_rct,
    node_spec = node_spec,
    ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n
  )
  source_boundary$cross_regimen_source_boundary_trace <-
    .ltmle_exact_boundary_trace_row(
      consuming_task = consuming_task,
      source_task = virtual_task,
      metadata = source_boundary$provenance,
      rows = source_boundary$rows,
      source_fit = source_fit,
      source_covariates = source_fit$conditioning_covariates %||% character(0),
      qsource_predicted = TRUE,
      cached_targeted_continuation_used = TRUE
    )
  out <- .ltmle_exact_collapse_source_qstar_to_subjects(
    q_source_star = q_source_star,
    source_boundary = source_boundary,
    task = consuming_task,
    source_task_id = source_task_id
  )
  virtual_downstream_id <- as.character(
    virtual_task$observed_pseudooutcome_source_task_id %||% NA_character_
  )
  virtual_downstream_fit <- if (!.ltmle_exact_missing_or_empty(virtual_downstream_id)) {
    task_fit_cache[[virtual_downstream_id]]
  } else {
    NULL
  }
  out$metadata <- c(
    source_boundary$diagnostics,
    out$metadata %||% list(),
    list(
      source_eval_mode = "branch_state_dp_continuation_value",
      cached_targeted_continuation_used = TRUE,
      continuation_task_id = source_task_id,
      continuation_fit_found = !is.null(source_fit),
      continuation_fit_targeted = isTRUE(source_fit$targeted),
      terminal_outcome_base_case_used = FALSE,
      terminal_full_recursion_used = FALSE,
      local_qstar_source_prediction_used = FALSE,
      local_current_task_prediction_used = FALSE,
      actual_helper_name = actual_helper_name,
      virtual_mixed_direct_continuation_used = FALSE,
      direct_continuation_allowed = isTRUE(direct_status$direct_continuation_allowed),
      downstream_source_can_represent_virtual_target =
        isTRUE(direct_status$downstream_source_can_represent_virtual_target),
      dedicated_virtual_Q_target_used =
        isTRUE(direct_status$dedicated_virtual_Q_target_used),
      virtual_mixed_downstream_source_task_id = virtual_downstream_id,
      virtual_mixed_downstream_source_fit_found = !is.null(virtual_downstream_fit),
      virtual_mixed_downstream_source_fit_targeted =
        isTRUE(virtual_downstream_fit$targeted)
    )
  )
  out
}

.ltmle_exact_eval_virtual_mixed_cached_continuation <- function(
  virtual_task,
  branch_state,
  task_graph,
  task_fit_cache,
  models,
  treatment_models,
  censoring_models,
  T,
  spec,
  probability_bounds,
  treat_mech,
  p_rct,
  node_spec = NULL,
  ltmle_exact_density_ratio_mc_n = 2000L,
  ltmle_exact_law_integration_n = 5L
) {
  virtual_task <- as.list(virtual_task)
  if (!.ltmle_exact_is_virtual_mixed_task(virtual_task)) {
    .stop("virtual_mixed_cached_continuation_requires_virtual_task")
  }
  downstream_source_task_id <- as.character(
    virtual_task$observed_pseudooutcome_source_task_id %||% NA_character_
  )
  if (.ltmle_exact_missing_or_empty(downstream_source_task_id) ||
      identical(downstream_source_task_id, "observed_terminal_outcome")) {
    .stop("virtual_mixed_continuation_missing_downstream_cached_source: ",
          virtual_task$task_id)
  }
  source_fit <- task_fit_cache[[downstream_source_task_id]]
  if (is.null(source_fit)) {
    .stop("virtual_mixed_continuation_downstream_fit_missing: ",
          downstream_source_task_id)
  }
  source_boundary <- .ltmle_exact_materialize_source_task_boundary(
    branch_state = branch_state,
    task = virtual_task,
    source_task_id = downstream_source_task_id,
    task_graph = task_graph,
    models = models,
    treatment_models = treatment_models,
    censoring_models = censoring_models,
    T = T,
    spec = spec,
    probability_bounds = probability_bounds,
    treat_mech = treat_mech,
    p_rct = p_rct,
    node_spec = node_spec,
    ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
    ltmle_exact_law_integration_n = ltmle_exact_law_integration_n,
    forbid_terminal_full_recursion = TRUE
  )
  q_source_star <- .ltmle_exact_predict_task_qstar_on_rows(
    task_fit = source_fit,
    rows = source_boundary$rows,
    row_long = source_boundary$long,
    models = models,
    treatment_models = treatment_models,
    censoring_models = censoring_models,
    T = T,
    spec = spec,
    probability_bounds = probability_bounds,
    treat_mech = treat_mech,
    p_rct = p_rct,
    node_spec = node_spec,
    ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n
  )
  source_boundary$cross_regimen_source_boundary_trace <-
    .ltmle_exact_boundary_trace_row(
      consuming_task = virtual_task,
      source_task = source_boundary$source_task,
      metadata = source_boundary$provenance,
      rows = source_boundary$rows,
      source_fit = source_fit,
      source_covariates = source_fit$conditioning_covariates %||% character(0),
      qsource_predicted = TRUE,
      cached_targeted_continuation_used = TRUE
    )
  out <- .ltmle_exact_collapse_source_qstar_to_subjects(
    q_source_star = q_source_star,
    source_boundary = source_boundary,
    task = virtual_task,
    source_task_id = downstream_source_task_id
  )
  out$metadata <- c(
    out$metadata %||% list(),
    list(
      virtual_mixed_direct_continuation_used = TRUE,
      direct_continuation_allowed = TRUE,
      downstream_source_can_represent_virtual_target = TRUE,
      dedicated_virtual_Q_target_used = FALSE,
      virtual_mixed_downstream_source_task_id = downstream_source_task_id,
      virtual_mixed_downstream_source_fit_found = !is.null(source_fit),
      virtual_mixed_downstream_source_fit_targeted = isTRUE(source_fit$targeted)
    )
  )
  out
}

.ltmle_exact_independent_eval_virtual_mixed_cached_continuation <- function(
  virtual_task,
  branch_state,
  task_graph,
  task_fit_cache,
  models,
  treatment_models,
  censoring_models,
  T,
  spec,
  probability_bounds,
  treat_mech,
  p_rct,
  node_spec = NULL,
  ltmle_exact_density_ratio_mc_n = 2000L,
  ltmle_exact_law_integration_n = 5L
) {
  virtual_task <- as.list(virtual_task)
  if (!.ltmle_exact_is_virtual_mixed_task(virtual_task)) {
    .stop("virtual_mixed_cached_continuation_requires_virtual_task")
  }
  downstream_source_task_id <- as.character(
    virtual_task$observed_pseudooutcome_source_task_id %||% NA_character_
  )
  if (.ltmle_exact_missing_or_empty(downstream_source_task_id) ||
      identical(downstream_source_task_id, "observed_terminal_outcome")) {
    .stop("virtual_mixed_continuation_missing_downstream_cached_source: ",
          virtual_task$task_id)
  }
  source_fit <- task_fit_cache[[downstream_source_task_id]]
  if (is.null(source_fit)) {
    .stop("virtual_mixed_continuation_downstream_fit_missing: ",
          downstream_source_task_id)
  }
  source_boundary <- .ltmle_exact_audit_materialize_source_task_boundary(
    branch_state = branch_state,
    task = virtual_task,
    source_task_id = downstream_source_task_id,
    task_graph = task_graph,
    models = models,
    treatment_models = treatment_models,
    censoring_models = censoring_models,
    T = T,
    spec = spec,
    probability_bounds = probability_bounds,
    treat_mech = treat_mech,
    p_rct = p_rct,
    node_spec = node_spec,
    ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
    ltmle_exact_law_integration_n = ltmle_exact_law_integration_n,
    forbid_terminal_full_recursion = TRUE
  )
  q_source_star <- .ltmle_exact_predict_task_qstar_on_rows(
    task_fit = source_fit,
    rows = source_boundary$rows,
    row_long = source_boundary$long,
    models = models,
    treatment_models = treatment_models,
    censoring_models = censoring_models,
    T = T,
    spec = spec,
    probability_bounds = probability_bounds,
    treat_mech = treat_mech,
    p_rct = p_rct,
    node_spec = node_spec,
    ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n
  )
  state_key <- source_boundary$state_key
  value <- .ltmle_exact_collapse_particles_to_subjects(
    as.numeric(q_source_star),
    source_boundary$branch_state,
    state_key
  )
  list(
    value = value,
    id = seq_along(value),
    state_key = rep(state_key, length(value)),
    metadata = c(
      source_boundary$diagnostics,
      list(
        source_eval_mode = "branch_state_dp_continuation_value",
        cached_targeted_continuation_used = TRUE,
        continuation_task_id = downstream_source_task_id,
        continuation_fit_found = !is.null(source_fit),
        continuation_fit_targeted = isTRUE(source_fit$targeted),
        terminal_outcome_base_case_used = FALSE,
        terminal_full_recursion_used = FALSE,
        local_qstar_source_prediction_used = FALSE,
        local_current_task_prediction_used = FALSE,
        actual_helper_name =
          ".ltmle_exact_independent_source_boundary_dp_continuation",
        virtual_mixed_direct_continuation_used = TRUE,
        direct_continuation_allowed = TRUE,
        downstream_source_can_represent_virtual_target = TRUE,
        dedicated_virtual_Q_target_used = FALSE,
        virtual_mixed_downstream_source_task_id = downstream_source_task_id,
        virtual_mixed_downstream_source_fit_found = !is.null(source_fit),
        virtual_mixed_downstream_source_fit_targeted =
          isTRUE(source_fit$targeted)
      )
    )
  )
}

.ltmle_exact_eval_source_pseudooutcome_dp <- function(
  task,
  observed_training_rows,
  task_graph,
  task_fit_cache,
  models,
  treatment_models,
  censoring_models,
  T,
  spec,
  probability_bounds,
  treat_mech,
  p_rct,
  node_spec = NULL,
  ltmle_exact_density_ratio_mc_n = 2000L,
  ltmle_exact_law_integration_n = 5L,
  force_epsilon_zero = FALSE
) {
  task <- as.list(task)
  source_info <- .ltmle_exact_source_task_lookup(task, task_graph)
  source_task_id <- source_info$source_task_id
  source_task_is_terminal_outcome <- isTRUE(source_info$source_task_is_terminal_outcome)
  requires_branch <- isTRUE(source_info$requires_branch)

  if (.ltmle_exact_missing_or_empty(source_task_id)) {
    .stop("missing_source_task_id")
  }

  boundary_state <- .ltmle_exact_init_branch_state_from_observed_task(
    observed_task_rows = observed_training_rows,
    task = task,
    spec = spec,
    node_spec = node_spec
  )

  if (.ltmle_exact_is_virtual_mixed_task(task) &&
      !.ltmle_exact_virtual_mixed_can_use_downstream_cached_source(task)) {
    return(.ltmle_exact_eval_virtual_mixed_dedicated_L_continuation(
      virtual_task = task,
      branch_state = boundary_state,
      task_graph = task_graph,
      task_fit_cache = task_fit_cache,
      models = models,
      treatment_models = treatment_models,
      censoring_models = censoring_models,
      T = T,
      spec = spec,
      probability_bounds = probability_bounds,
      treat_mech = treat_mech,
      p_rct = p_rct,
      node_spec = node_spec,
      ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
      ltmle_exact_law_integration_n = ltmle_exact_law_integration_n,
      actual_helper_name = ".ltmle_exact_eval_source_pseudooutcome_dp",
      integration_context = "virtual_mixed_dedicated_L_transition_target"
    ))
  }

  if (source_task_is_terminal_outcome) {
    return(.ltmle_exact_eval_terminal_source_base_case(
      branch_state = boundary_state,
      task = task,
      spec = spec,
      node_spec = node_spec
    ))
  }

  if (requires_branch) {
    .ltmle_exact_validate_source_task_for_cached_continuation(
      task = task,
      source_task_id = source_task_id,
      task_graph = task_graph,
      task_fit_cache = task_fit_cache
    )

    source_boundary <- .ltmle_exact_materialize_source_task_boundary(
      branch_state = boundary_state,
      task = task,
      source_task_id = source_task_id,
      task_graph = task_graph,
      models = models,
      treatment_models = treatment_models,
      censoring_models = censoring_models,
      T = T,
      spec = spec,
      probability_bounds = probability_bounds,
      treat_mech = treat_mech,
      p_rct = p_rct,
      node_spec = node_spec,
      ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
      ltmle_exact_law_integration_n = ltmle_exact_law_integration_n,
      forbid_terminal_full_recursion = TRUE
    )

	    if (.ltmle_exact_is_virtual_mixed_task(source_boundary$source_task) &&
	        .ltmle_exact_virtual_mixed_can_use_downstream_cached_source(source_boundary$source_task)) {
	      virtual_eval <- .ltmle_exact_eval_virtual_mixed_cached_continuation(
	        virtual_task = source_boundary$source_task,
	        branch_state = source_boundary$branch_state,
	        task_graph = task_graph,
	        task_fit_cache = task_fit_cache,
	        models = models,
	        treatment_models = treatment_models,
	        censoring_models = censoring_models,
	        T = T,
	        spec = spec,
	        probability_bounds = probability_bounds,
	        treat_mech = treat_mech,
	        p_rct = p_rct,
	        node_spec = node_spec,
	        ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
	        ltmle_exact_law_integration_n = ltmle_exact_law_integration_n
	      )
	      virtual_source_trace <- source_boundary$cross_regimen_source_boundary_trace %||% data.frame()
	      if (nrow(virtual_source_trace)) {
	        virtual_rows <- as.character(virtual_source_trace$source_boundary_type) ==
	          "virtual_mixed_continuation"
	        virtual_source_trace$qsource_predicted_on_mixed_boundary[virtual_rows] <- TRUE
	        virtual_source_trace$cached_targeted_continuation_used[virtual_rows] <- TRUE
	        virtual_source_trace$passed[virtual_rows] <- TRUE
	        virtual_source_trace$failure_class[virtual_rows] <- "no_failure"
	      }
	      virtual_eval$diagnostics$cross_regimen_source_boundary_trace <-
	        .ltmle_exact_rbind_fill(list(
	          virtual_source_trace,
	          virtual_eval$diagnostics$cross_regimen_source_boundary_trace %||% data.frame()
	        ))
	      virtual_eval$diagnostics$law_integration_diagnostics <-
	        .ltmle_exact_rbind_fill(list(
	          source_boundary$law_integration_diagnostics %||% data.frame(),
	          virtual_eval$diagnostics$law_integration_diagnostics %||% data.frame()
	        ))
	      virtual_eval$handoff_trace <- .ltmle_exact_rbind_fill(list(
	        source_boundary$handoff_trace %||% data.frame(),
	        virtual_eval$handoff_trace %||% data.frame()
	      ))
	      virtual_eval$metadata <- c(
	        source_boundary$diagnostics,
	        virtual_eval$metadata %||% list(),
	        list(actual_helper_name = ".ltmle_exact_eval_source_pseudooutcome_dp")
	      )
		      return(virtual_eval)
		    }

		    if (.ltmle_exact_is_virtual_mixed_task(source_boundary$source_task)) {
	      virtual_eval <- .ltmle_exact_eval_virtual_mixed_task_fit_on_source_boundary(
	        virtual_task = source_boundary$source_task,
	        consuming_task = task,
	        source_boundary = source_boundary,
	        task_fit_cache = task_fit_cache,
	        models = models,
	        treatment_models = treatment_models,
	        censoring_models = censoring_models,
	        T = T,
	        spec = spec,
	        probability_bounds = probability_bounds,
	        treat_mech = treat_mech,
	        p_rct = p_rct,
	        node_spec = node_spec,
	        ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
	        actual_helper_name = ".ltmle_exact_eval_source_pseudooutcome_dp"
	      )
	      return(virtual_eval)
	    }

	    source_fit <- task_fit_cache[[source_task_id]]
	    q_source_star <- .ltmle_exact_predict_task_qstar_on_rows(
	      task_fit = source_fit,
	      rows = source_boundary$rows,
	      row_long = source_boundary$long,
	      models = models,
	      treatment_models = treatment_models,
	      censoring_models = censoring_models,
	      T = T,
	      spec = spec,
	      probability_bounds = probability_bounds,
	      treat_mech = treat_mech,
	      p_rct = p_rct,
	      node_spec = node_spec,
	      ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n
	    )
	    source_boundary$cross_regimen_source_boundary_trace <- .ltmle_exact_boundary_trace_row(
	      consuming_task = task,
	      source_task = source_boundary$source_task,
	      metadata = source_boundary$provenance %||% .ltmle_exact_source_boundary_metadata(
	        consuming_task = task,
	        source_task = source_boundary$source_task,
	        spec = spec
	      ),
	      rows = source_boundary$rows,
	      source_fit = source_fit,
	      source_covariates = source_fit$conditioning_covariates %||% character(0),
	      qsource_predicted = TRUE,
	      cached_targeted_continuation_used = TRUE
	    )
	    return(.ltmle_exact_collapse_source_qstar_to_subjects(
	      q_source_star = q_source_star,
	      source_boundary = source_boundary,
	      task = task,
	      source_task_id = source_task_id
	    ))
		  }
	  .stop("Local source prediction is not explicitly allowed for task=", task$task_id)
}

.ltmle_exact_eval_downstream_source_from_branch_state <- function(source_task_id,
                                                                 branch_state,
                                                                 task_graph,
                                                                 task_fit_cache,
                                                                 models,
                                                                 treatment_models,
                                                                 censoring_models,
                                                                 T,
                                                                 spec,
                                                                 probability_bounds,
                                                                 treat_mech,
                                                                 p_rct,
                                                                 node_spec = NULL,
                                                                 ltmle_exact_density_ratio_mc_n = 2000L,
                                                                 ltmle_exact_law_integration_n = 5L) {
  if (identical(source_task_id, "observed_terminal_outcome")) {
    rows <- .ltmle_exact_materialize_branch_rows(branch_state, "outcome", T, node_spec)
    y <- as.numeric(rows$Y)
    if (any(!is.finite(y))) {
      .stop("Observed terminal source requested before finite Y was present in branch_state.")
    }
    value <- .ltmle_exact_collapse_particles_to_subjects(y, branch_state, "outcome")
    return(list(
      value = value,
      id = seq_along(value),
      state_key = rep("outcome", length(value)),
      particle_q = y,
      subject_q = value,
      particle_frame = branch_state$outcome$particles,
      branch_trace = data.frame(),
      handoff_trace = data.frame(),
      diagnostics = list(
        density_ratio_diagnostics = data.frame(),
        law_integration_diagnostics = data.frame(),
        H_new_mean = NA_real_,
        H_new_max = NA_real_
      ),
      metadata = list(
        branch_state_downstream_recursion_used = FALSE,
        terminal_outcome_base_case_used = TRUE,
        local_qstar_source_prediction_used = FALSE,
        nodes_updated = "observed_terminal_outcome",
        tasks_visited = "observed_terminal_outcome",
        used_observed_long_as_generated_history = FALSE,
        used_assigned_A_observed_long_mutation = FALSE,
        product_join_used = FALSE,
        law_integration_n = as.integer(ltmle_exact_law_integration_n),
        actual_helper_name = ".ltmle_exact_eval_downstream_source_from_branch_state"
      )
    ))
  }

  source_task <- NULL
  if (!is.null(task_graph$tasks) && source_task_id %in% names(task_graph$tasks)) {
    source_task <- as.list(task_graph$tasks[[source_task_id]])
  }
  if (is.null(source_task)) {
    .stop("Source task not found for branch-state downstream recursion: ", source_task_id)
  }
  source_state_key <- .ltmle_exact_task_branch_key(source_task)
  eval <- .ltmle_exact_eval_task_on_branch_state(
    task_id = source_task_id,
    branch_state = branch_state,
    task_graph = task_graph,
    task_fit_cache = task_fit_cache,
    models = models,
    treatment_models = treatment_models,
    censoring_models = censoring_models,
    T = T,
    spec = spec,
    probability_bounds = probability_bounds,
    treat_mech = treat_mech,
    p_rct = p_rct,
    node_spec = node_spec,
    ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
    ltmle_exact_law_integration_n = ltmle_exact_law_integration_n
  )
  value <- as.numeric(eval$subject_q)
  if (length(value) != .ltmle_exact_branch_state_subject_n(branch_state) ||
      any(!is.finite(value))) {
    .stop("Branch-state source pseudooutcome returned invalid values for source_task_id=", source_task_id)
  }
  eval$value <- value
  eval$id <- seq_along(value)
  eval$state_key <- rep(source_state_key, length(value))
  eval$metadata <- eval$metadata %||% list()
  eval$metadata$branch_state_downstream_recursion_used <- TRUE
  eval$metadata$terminal_outcome_base_case_used <- FALSE
  eval$metadata$local_qstar_source_prediction_used <- FALSE
  eval$metadata$nodes_updated <- unique(as.character(eval$branch_trace$node %||% source_task$node %||% NA_character_))
  eval$metadata$tasks_visited <- unique(as.character(eval$branch_trace$task_id %||% source_task_id))
  eval$metadata$used_observed_long_as_generated_history <- FALSE
  eval$metadata$used_assigned_A_observed_long_mutation <- FALSE
  eval$metadata$product_join_used <- if (!is.null(eval$handoff_trace) && nrow(eval$handoff_trace) &&
      "handoff_type" %in% names(eval$handoff_trace)) {
    any(eval$handoff_trace$handoff_type == "separate_world_product_join", na.rm = TRUE)
  } else {
    FALSE
  }
  eval$metadata$law_integration_n <- as.integer(ltmle_exact_law_integration_n)
  eval$metadata$actual_helper_name <- ".ltmle_exact_eval_downstream_source_from_branch_state"
  eval
}

.ltmle_exact_source_eval_training_diag_row <- function(task,
                                                       source_task_id,
                                                       source_info = NULL,
                                                       source_eval_mode,
                                                       branch_state_source_boundary_materialized = FALSE,
                                                       source_boundary_reached = FALSE,
                                                       source_boundary_task_id = NA_character_,
                                                       source_boundary_type = "pure",
                                                       source_boundary_direction = "not_applicable",
                                                       source_boundary_eval_state = "not_applicable",
                                                       source_boundary_outcome_history_state = "not_applicable",
                                                       source_boundary_m1_history_state = "not_applicable",
                                                       source_boundary_m2_history_state = "not_applicable",
                                                       source_boundary_auxiliary_mediator_history_state = "not_applicable",
                                                       qsource_predicted_on_mixed_boundary = FALSE,
                                                       evaluated_beyond_source_boundary = FALSE,
                                                       terminal_task_evaluated_in_source_path = FALSE,
                                                       max_task_depth_beyond_source_boundary = 0L,
                                                       cached_targeted_continuation_used = FALSE,
                                                       continuation_task_id = NA_character_,
                                                       continuation_fit_found = FALSE,
                                                       continuation_fit_targeted = FALSE,
                                                       branch_state_downstream_recursion_used = FALSE,
                                                       terminal_full_recursion_used = FALSE,
                                                       terminal_outcome_base_case_used = FALSE,
                                                       local_qstar_source_prediction_used = FALSE,
                                                       local_current_task_prediction_used = FALSE,
                                                       source_eval_helper_called = TRUE,
                                                       helper_name = NA_character_,
                                                       actual_helper_name = helper_name,
                                                       requires_branch = FALSE,
                                                       requires_branch_exception = FALSE,
                                                       requires_branch_exception_reason = NA_character_,
                                                       justification_for_local_prediction = NA_character_,
                                                       matched_allowed_local_source_prediction_map = FALSE,
                                                       allowed_local_prediction = FALSE,
                                                       n_allowed_local_source_prediction_matches = 0L,
                                                       allowed_reason_code = NA_character_,
                                                       local_prediction_allowlist_audit_passed = TRUE,
                                                       observed_training_rows_used_for_training = TRUE,
                                                       observed_rows_mutated_to_assigned_A = FALSE,
                                                       observed_long_used_as_generated_history = FALSE,
                                                       used_observed_long_as_source_history = FALSE,
                                                       n_source_boundary_rows = 0L,
                                                       n_source_boundary_state_keys = 0L,
                                                       n_particles_before_collapse = 0L,
                                                       n_particles_after_collapse = 0L,
                                                       safe_collapse_applied = FALSE,
                                                       safe_collapse_key = "not_applicable",
                                                       used_as_training_outcome = FALSE,
                                                       training_outcome_source = "not_used",
                                                       training_outcome_mean = 0,
                                                       dp_source_eval_mean = 0,
                                                       initial_dp_source_eval_available = FALSE,
                                                       initial_dp_source_eval_mean = NA_real_,
                                                       initial_dp_source_eval_source = "not_available",
                                                       targeted_dp_source_eval_mean = NA_real_,
                                                       targeted_minus_initial_dp_source_eval = NA_real_,
                                                       legacy_graph_observed_value_used = FALSE,
                                                       legacy_graph_observed_value_mean = 0,
                                                       training_outcome_minus_dp_source_eval = 0,
                                                       stored_graph_value_quantity = "not_stored_yet",
                                                       stored_graph_value_mean = 0,
                                                       stored_graph_value_minus_current_qstar_observed = 0,
                                                       self_source_reference_allowed = FALSE,
                                                       self_source_reference_reason = NA_character_,
                                                       virtual_mixed_direct_continuation_used = FALSE,
                                                       direct_continuation_allowed = FALSE,
                                                       downstream_source_can_represent_virtual_target = FALSE,
                                                       dedicated_virtual_Q_target_used = FALSE,
                                                       virtual_mixed_downstream_source_task_id = NA_character_,
                                                       virtual_mixed_downstream_source_fit_found = FALSE,
                                                       virtual_mixed_downstream_source_fit_targeted = FALSE,
                                                       T = NA_integer_,
                                                       passed = TRUE,
                                                       failure_class = "no_failure") {
  task <- as.list(task)
  source_info <- source_info %||% list()
  source_task_id <- as.character(source_task_id %||% source_info$source_task_id %||% NA_character_)
  source_task_is_terminal_outcome <- isTRUE(source_info$source_task_is_terminal_outcome) ||
    identical(source_task_id, "observed_terminal_outcome")
  source_task_t <- as.integer(source_info$source_task_t %||% NA_integer_)
  if (source_task_is_terminal_outcome && (!is.finite(source_task_t) || is.na(source_task_t))) {
    source_task_t <- as.integer(T %||% task$t %||% NA_integer_)
  }
  source_task_node <- as.character(source_info$source_task_node %||% NA_character_)
  source_task_process_type <- as.character(source_info$source_task_process_type %||% NA_character_)
  if (source_task_is_terminal_outcome) {
    source_task_node <- "observed_terminal_outcome"
    source_task_process_type <- "observed_terminal_outcome"
  }

  helper_name <- as.character(helper_name %||% "not_applicable")
  actual_helper_name <- as.character(actual_helper_name %||% helper_name)
  source_boundary_task_id <- as.character(source_boundary_task_id %||% source_task_id)
  continuation_task_id <- as.character(continuation_task_id %||% "not_applicable")
  requires_branch_exception_reason <- as.character(requires_branch_exception_reason %||% NA_character_)
  justification_for_local_prediction <- as.character(justification_for_local_prediction %||% NA_character_)
  allowed_reason_code <- as.character(allowed_reason_code %||% NA_character_)

  if ((isTRUE(branch_state_downstream_recursion_used) || isTRUE(cached_targeted_continuation_used)) &&
      !isTRUE(local_qstar_source_prediction_used)) {
    if (is.na(requires_branch_exception_reason) || !nzchar(requires_branch_exception_reason)) {
      requires_branch_exception_reason <- "not_applicable"
    }
    if (is.na(justification_for_local_prediction) || !nzchar(justification_for_local_prediction)) {
      justification_for_local_prediction <- "not_applicable"
    }
    if (is.na(allowed_reason_code) || !nzchar(allowed_reason_code)) {
      allowed_reason_code <- "not_applicable"
    }
  }
  if (isTRUE(terminal_outcome_base_case_used)) {
    requires_branch_exception_reason <- "terminal_source_task"
    justification_for_local_prediction <- "terminal_source_task"
    allowed_reason_code <- "not_applicable"
  }
  if (isTRUE(local_qstar_source_prediction_used)) {
    if (is.na(justification_for_local_prediction) || !nzchar(justification_for_local_prediction)) {
      justification_for_local_prediction <- allowed_reason_code
    }
    if (is.na(requires_branch_exception_reason) || !nzchar(requires_branch_exception_reason)) {
      requires_branch_exception_reason <- allowed_reason_code
    }
  }
  if (is.na(self_source_reference_reason) || !nzchar(self_source_reference_reason)) {
    self_source_reference_reason <- if (isTRUE(self_source_reference_allowed)) {
      "self_source_reference_allowed_by_task_graph"
    } else {
      "not_applicable"
    }
  }

  reason_matches <- if (isTRUE(local_qstar_source_prediction_used)) {
    identical(requires_branch_exception_reason, allowed_reason_code) &&
      identical(justification_for_local_prediction, allowed_reason_code)
  } else {
    TRUE
  }
  local_prediction_allowlist_audit_passed <- if (isTRUE(local_qstar_source_prediction_used)) {
    isTRUE(matched_allowed_local_source_prediction_map) &&
      identical(as.integer(n_allowed_local_source_prediction_matches), 1L) &&
      isTRUE(allowed_local_prediction) &&
      isTRUE(reason_matches)
  } else {
    isTRUE(local_prediction_allowlist_audit_passed)
  }
	  virtual_dedicated_target_path <- .ltmle_exact_is_virtual_mixed_task(task) &&
	    isTRUE(dedicated_virtual_Q_target_used) &&
	    !isTRUE(virtual_mixed_direct_continuation_used) &&
	    !isTRUE(downstream_source_can_represent_virtual_target)
  allowed_path <- if (identical(source_task_id, "observed_terminal_outcome")) {
    isTRUE(terminal_outcome_base_case_used) &&
      !isTRUE(branch_state_downstream_recursion_used) &&
      !isTRUE(terminal_full_recursion_used) &&
      !isTRUE(local_qstar_source_prediction_used)
  } else if (identical(requires_branch, TRUE)) {
    common_branch_path <- isTRUE(branch_state_source_boundary_materialized) &&
      isTRUE(source_boundary_reached) &&
      !isTRUE(evaluated_beyond_source_boundary) &&
      !isTRUE(terminal_task_evaluated_in_source_path) &&
      identical(as.integer(max_task_depth_beyond_source_boundary), 0L) &&
      isTRUE(cached_targeted_continuation_used) &&
      isTRUE(continuation_fit_found) &&
      isTRUE(continuation_fit_targeted) &&
      !isTRUE(branch_state_downstream_recursion_used) &&
      !isTRUE(terminal_full_recursion_used) &&
      !isTRUE(terminal_outcome_base_case_used) &&
      !isTRUE(local_current_task_prediction_used) &&
      !isTRUE(local_qstar_source_prediction_used)
	    common_branch_path &&
	      if (isTRUE(virtual_dedicated_target_path)) {
	        !identical(as.character(source_boundary_task_id), as.character(source_task_id)) &&
	          identical(as.character(source_boundary_task_id), as.character(continuation_task_id))
	      } else {
	        identical(as.character(source_boundary_task_id), as.character(source_task_id)) &&
	          identical(as.character(continuation_task_id), as.character(source_task_id))
      }
  } else if (isTRUE(local_qstar_source_prediction_used)) {
    isTRUE(local_prediction_allowlist_audit_passed)
  } else {
    FALSE
  }
  no_required_source_missing <- !is.na(source_task_id) && nzchar(source_task_id)
  helper_ok <- identical(helper_name, actual_helper_name)
  passed <- isTRUE(passed) &&
    !isTRUE(observed_rows_mutated_to_assigned_A) &&
    !isTRUE(observed_long_used_as_generated_history) &&
    !isTRUE(used_observed_long_as_source_history) &&
    no_required_source_missing &&
    helper_ok &&
    isTRUE(local_prediction_allowlist_audit_passed) &&
    isTRUE(reason_matches) &&
    isTRUE(allowed_path) &&
    !identical(requires_branch_exception_reason, "missing_source_task_id")
  data.frame(
    component = task$component %||% NA_character_,
    task_id = task$task_id %||% NA_character_,
    t = as.integer(task$t %||% NA_integer_),
    node = task$node %||% NA_character_,
    process_type = .ltmle_exact_process_type(task),
    requires_branch = isTRUE(requires_branch),
    source_task_id = source_task_id,
    source_task_t = source_task_t,
    source_task_node = source_task_node,
    source_task_process_type = source_task_process_type,
    source_task_is_terminal_outcome = isTRUE(source_task_is_terminal_outcome),
    requires_branch_exception = isTRUE(requires_branch_exception),
    requires_branch_exception_reason = requires_branch_exception_reason,
    source_eval_mode = source_eval_mode,
    branch_state_source_boundary_materialized = isTRUE(branch_state_source_boundary_materialized),
    source_boundary_reached = isTRUE(source_boundary_reached),
    source_boundary_task_id = source_boundary_task_id,
    source_boundary_type = as.character(source_boundary_type %||% "pure"),
    source_boundary_direction = as.character(source_boundary_direction %||% "not_applicable"),
    source_boundary_eval_state = as.character(source_boundary_eval_state %||% "not_applicable"),
    source_boundary_outcome_history_state = as.character(source_boundary_outcome_history_state %||% "not_applicable"),
    source_boundary_m1_history_state = as.character(source_boundary_m1_history_state %||% "not_applicable"),
    source_boundary_m2_history_state = as.character(source_boundary_m2_history_state %||% "not_applicable"),
    source_boundary_auxiliary_mediator_history_state =
      as.character(source_boundary_auxiliary_mediator_history_state %||% "not_applicable"),
    qsource_predicted_on_mixed_boundary = isTRUE(qsource_predicted_on_mixed_boundary),
    evaluated_beyond_source_boundary = isTRUE(evaluated_beyond_source_boundary),
    terminal_task_evaluated_in_source_path = isTRUE(terminal_task_evaluated_in_source_path),
    max_task_depth_beyond_source_boundary = as.integer(max_task_depth_beyond_source_boundary),
    cached_targeted_continuation_used = isTRUE(cached_targeted_continuation_used),
    continuation_task_id = continuation_task_id,
    continuation_fit_found = isTRUE(continuation_fit_found),
    continuation_fit_targeted = isTRUE(continuation_fit_targeted),
    branch_state_downstream_recursion_used = isTRUE(branch_state_downstream_recursion_used),
    terminal_full_recursion_used = isTRUE(terminal_full_recursion_used),
    terminal_outcome_base_case_used = isTRUE(terminal_outcome_base_case_used),
    local_qstar_source_prediction_used = isTRUE(local_qstar_source_prediction_used),
    local_current_task_prediction_used = isTRUE(local_current_task_prediction_used),
    source_eval_helper_called = isTRUE(source_eval_helper_called),
    helper_name = helper_name,
    actual_helper_name = actual_helper_name,
    justification_for_local_prediction = justification_for_local_prediction,
    matched_allowed_local_source_prediction_map = isTRUE(matched_allowed_local_source_prediction_map),
    allowed_local_prediction = isTRUE(allowed_local_prediction),
    n_allowed_local_source_prediction_matches = as.integer(n_allowed_local_source_prediction_matches),
    allowed_reason_code = allowed_reason_code,
    requires_branch_exception_reason_matches_allowed_reason_code = isTRUE(reason_matches),
    local_prediction_allowlist_audit_passed = isTRUE(local_prediction_allowlist_audit_passed),
    observed_training_rows_used_for_training = isTRUE(observed_training_rows_used_for_training),
    observed_rows_mutated_to_assigned_A = isTRUE(observed_rows_mutated_to_assigned_A),
    observed_long_used_as_generated_history = isTRUE(observed_long_used_as_generated_history),
    used_observed_long_as_source_history = isTRUE(used_observed_long_as_source_history),
    n_source_boundary_rows = as.integer(n_source_boundary_rows),
    n_source_boundary_state_keys = as.integer(n_source_boundary_state_keys),
    n_particles_before_collapse = as.integer(n_particles_before_collapse),
    n_particles_after_collapse = as.integer(n_particles_after_collapse),
    safe_collapse_applied = isTRUE(safe_collapse_applied),
    safe_collapse_key = as.character(safe_collapse_key %||% "not_applicable"),
    used_as_training_outcome = isTRUE(used_as_training_outcome),
    training_outcome_source = as.character(training_outcome_source %||% "not_used"),
    training_outcome_mean = as.numeric(training_outcome_mean),
    dp_source_eval_mean = as.numeric(dp_source_eval_mean),
    initial_dp_source_eval_available = isTRUE(initial_dp_source_eval_available),
    initial_dp_source_eval_mean = as.numeric(initial_dp_source_eval_mean),
    initial_dp_source_eval_source = as.character(initial_dp_source_eval_source %||% "not_available"),
    targeted_dp_source_eval_mean = as.numeric(targeted_dp_source_eval_mean),
    targeted_minus_initial_dp_source_eval = as.numeric(targeted_minus_initial_dp_source_eval),
    legacy_graph_observed_value_used = isTRUE(legacy_graph_observed_value_used),
    legacy_graph_observed_value_mean = as.numeric(legacy_graph_observed_value_mean),
    training_outcome_minus_dp_source_eval = as.numeric(training_outcome_minus_dp_source_eval),
    stored_graph_value_quantity = as.character(stored_graph_value_quantity %||% "not_stored_yet"),
    stored_graph_value_mean = as.numeric(stored_graph_value_mean),
    stored_graph_value_minus_current_qstar_observed = as.numeric(stored_graph_value_minus_current_qstar_observed),
    self_source_reference_allowed = isTRUE(self_source_reference_allowed),
    self_source_reference_reason = self_source_reference_reason,
    virtual_mixed_direct_continuation_used = isTRUE(virtual_mixed_direct_continuation_used),
    direct_continuation_allowed = isTRUE(direct_continuation_allowed),
    downstream_source_can_represent_virtual_target =
      isTRUE(downstream_source_can_represent_virtual_target),
		    dedicated_virtual_Q_target_used = isTRUE(dedicated_virtual_Q_target_used),
		    virtual_mixed_downstream_source_task_id =
      as.character(virtual_mixed_downstream_source_task_id %||% NA_character_),
    virtual_mixed_downstream_source_fit_found =
      isTRUE(virtual_mixed_downstream_source_fit_found),
    virtual_mixed_downstream_source_fit_targeted =
      isTRUE(virtual_mixed_downstream_source_fit_targeted),
    passed = isTRUE(passed),
    failure_class = if (isTRUE(passed)) "no_failure" else failure_class,
    stringsAsFactors = FALSE
  )
}

.ltmle_exact_source_task_is_terminal_local <- function(task, task_graph) {
  task <- as.list(task)
  identical(as.character(task$observed_pseudooutcome_source_task_id %||% NA_character_),
            "observed_terminal_outcome")
}

.ltmle_exact_source_task_requires_branch_recursion <- function(task, task_graph, node_spec = NULL) {
  force(node_spec)
  isTRUE(.ltmle_exact_source_task_lookup(task, task_graph)$requires_branch)
}

.ltmle_exact_source_independent_check_required <- function(task, task_graph, node_spec = NULL) {
  .ltmle_exact_source_task_requires_branch_recursion(task, task_graph, node_spec)
}

.ltmle_exact_source_task_local_justification <- function(task, task_graph, node_spec = NULL) {
  force(task_graph)
  force(node_spec)
  task <- as.list(task)
  if (identical(as.character(task$observed_pseudooutcome_source_task_id %||% NA_character_),
                "observed_terminal_outcome")) {
    return("terminal_source_task")
  }
  NA_character_
}

.ltmle_exact_local_source_qstar_prediction <- function(task,
                                                       branch_state,
                                                       task_fit_cache,
                                                       models,
                                                       treatment_models,
                                                       censoring_models,
                                                       T,
                                                       spec,
                                                       probability_bounds,
                                                       treat_mech,
                                                       p_rct,
                                                       node_spec = NULL,
                                                       ltmle_exact_density_ratio_mc_n = 2000L) {
  task <- as.list(task)
  task_id <- as.character(task$task_id)
  task_fit <- task_fit_cache[[task_id]]
  if (is.null(task_fit)) .stop("Missing task fit for local source Q* prediction: ", task_id)
  source_eval <- .ltmle_exact_make_task_source_eval(
    task = task,
    observed_rows = branch_state[[.ltmle_exact_task_branch_key(task)]]$rows,
    observed_long = NULL,
    spec = spec,
    node_spec = node_spec,
    branch_state = branch_state
  )
  q <- .ltmle_exact_predict_task_qstar_on_rows(
    task_fit = task_fit,
    rows = source_eval$rows,
    row_long = source_eval$long,
    models = models,
    treatment_models = treatment_models,
    censoring_models = censoring_models,
    T = T,
    spec = spec,
    probability_bounds = probability_bounds,
    treat_mech = treat_mech,
    p_rct = p_rct,
    node_spec = node_spec,
    ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n
  )
  list(
    value = as.numeric(q),
    branch_trace = data.frame(),
    handoff_trace = data.frame(),
    diagnostics = list(
      density_ratio_diagnostics = attr(q, "density_ratio_diagnostics") %||% data.frame(),
      truncation_diagnostics = attr(q, "truncation_diagnostics") %||% data.frame(),
      clever_covariate_decomposition_diagnostics =
        attr(q, "clever_covariate_decomposition_diagnostics") %||% data.frame(),
      law_integration_diagnostics = data.frame(),
      H_new_mean = attr(q, "H_new_mean") %||% NA_real_,
      H_new_max = attr(q, "H_new_max") %||% NA_real_
    )
  )
}

.ltmle_exact_eval_allowed_local_source_prediction <- function(task,
                                                             branch_state,
                                                             task_fit_cache,
                                                             models,
                                                             treatment_models,
                                                             censoring_models,
                                                             T,
                                                             spec,
                                                             probability_bounds,
                                                             treat_mech,
                                                             p_rct,
                                                             node_spec = NULL,
                                                             ltmle_exact_density_ratio_mc_n = 2000L) {
  .ltmle_exact_local_source_qstar_prediction(
    task = task,
    branch_state = branch_state,
    task_fit_cache = task_fit_cache,
    models = models,
    treatment_models = treatment_models,
    censoring_models = censoring_models,
    T = T,
    spec = spec,
    probability_bounds = probability_bounds,
    treat_mech = treat_mech,
    p_rct = p_rct,
    node_spec = node_spec,
    ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n
  )
}

.ltmle_exact_eval_source_pseudooutcome_dispatch <- function(
  task,
  task_data,
  observed_training_rows,
  observed_long,
  task_graph,
  task_fit_cache,
  models,
  treatment_models,
  censoring_models,
  T,
  spec,
  probability_bounds,
  treat_mech,
  p_rct,
  node_spec = NULL,
  ltmle_exact_density_ratio_mc_n = 2000L,
  ltmle_exact_law_integration_n = 5L,
  allow_local_source_prediction = FALSE,
  allowed_local_source_prediction_map = NULL,
  force_epsilon_zero = FALSE
) {
  force(task_data)
  force(observed_long)
  task <- as.list(task)
  task_id <- as.character(task$task_id)
  source_info0 <- .ltmle_exact_source_task_lookup(task, task_graph)
  source_task_id <- source_info0$source_task_id
  source_task_is_terminal_outcome <- isTRUE(source_info0$source_task_is_terminal_outcome)
  requires_branch <- isTRUE(source_info0$requires_branch)

  if (.ltmle_exact_missing_or_empty(source_task_id)) {
    .stop("missing_source_task_id")
  }
  self_allowed <- isTRUE(task$source_self_reference_allowed %||% FALSE)
  if (identical(source_task_id, task_id) && !self_allowed) {
    .stop("source_task_self_reference: ", source_task_id)
  }

  allowed_local_source_prediction_map <- allowed_local_source_prediction_map %||%
    .ltmle_exact_empty_allowed_local_source_prediction_map()
  helper_name <- ".ltmle_exact_eval_source_pseudooutcome_dp"
  local_helper_name <- ".ltmle_exact_local_source_qstar_prediction"
  actual_helper_name <- helper_name
  terminal_used <- FALSE
  downstream_recursion_used <- FALSE
  local_used <- FALSE
  requires_branch_exception <- FALSE
  requires_branch_exception_reason <- NA_character_
  matched_allowed_local_source_prediction_map <- FALSE
  allowed_local_prediction <- FALSE
  n_allowed_local_source_prediction_matches <- 0L
  allowed_reason_code <- NA_character_
  justification <- NA_character_

  if (isTRUE(source_task_is_terminal_outcome)) {
    eval <- .ltmle_exact_eval_source_pseudooutcome_dp(
      task = task,
      observed_training_rows = observed_training_rows,
      task_graph = task_graph,
      task_fit_cache = task_fit_cache,
      models = models,
      treatment_models = treatment_models,
      censoring_models = censoring_models,
      T = T,
      spec = spec,
      probability_bounds = probability_bounds,
      treat_mech = treat_mech,
      p_rct = p_rct,
      node_spec = node_spec,
      ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
      ltmle_exact_law_integration_n = ltmle_exact_law_integration_n,
      force_epsilon_zero = force_epsilon_zero
    )
    source_eval_mode <- "branch_state_terminal_outcome_base_case"
    helper_name <- ".ltmle_exact_eval_terminal_source_base_case"
    actual_helper_name <- helper_name
    terminal_used <- TRUE
    requires_branch <- FALSE
    requires_branch_exception <- TRUE
    requires_branch_exception_reason <- "terminal_source_task"
    justification <- "terminal_source_task"
  } else {
    allow <- .ltmle_exact_match_allowed_local_source_prediction(
      task = task,
      source_info = source_info0,
      allowed_local_source_prediction_map = allowed_local_source_prediction_map
    )
    matched_allowed_local_source_prediction_map <- isTRUE(allow$matched)
    allowed_local_prediction <- isTRUE(allow$allowed_local_prediction)
    n_allowed_local_source_prediction_matches <- as.integer(allow$n_matches %||% 0L)
    allowed_reason_code <- as.character(allow$allowed_reason_code %||% NA_character_)

    requires_branch <- isTRUE(source_info0$default_requires_branch)
    requires_branch_exception <- FALSE
    requires_branch_exception_reason <- NA_character_
    if (isTRUE(allow$matched) && isTRUE(allow$allowed_local_prediction)) {
      requires_branch <- FALSE
      requires_branch_exception <- TRUE
      requires_branch_exception_reason <- as.character(allow$allowed_reason_code)
    }

    if (identical(requires_branch, TRUE)) {
      helper_name <- ".ltmle_exact_eval_source_pseudooutcome_dp"
      actual_helper_name <- helper_name
      eval <- .ltmle_exact_eval_source_pseudooutcome_dp(
        task = task,
        observed_training_rows = observed_training_rows,
        task_graph = task_graph,
        task_fit_cache = task_fit_cache,
        models = models,
        treatment_models = treatment_models,
        censoring_models = censoring_models,
        T = T,
        spec = spec,
        probability_bounds = probability_bounds,
        treat_mech = treat_mech,
        p_rct = p_rct,
        node_spec = node_spec,
        ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
        ltmle_exact_law_integration_n = ltmle_exact_law_integration_n,
        force_epsilon_zero = force_epsilon_zero
      )
      source_eval_mode <- "branch_state_dp_continuation_value"
      downstream_recursion_used <- FALSE
      local_used <- FALSE
      requires_branch_exception_reason <- "not_applicable"
      justification <- "not_applicable"
      allowed_reason_code <- "not_applicable"
    } else {
      if (!isTRUE(allow_local_source_prediction) ||
          !isTRUE(allow$matched) ||
          !isTRUE(allow$allowed_local_prediction)) {
        .stop("Nonterminal source task is not branch-recursive and does not match allowed_local_source_prediction_map: task=",
              task_id, ", source_task_id=", source_task_id)
      }
      boundary_state <- .ltmle_exact_init_branch_state_from_observed_task(
        observed_task_rows = observed_training_rows,
        task = task,
        spec = spec,
        node_spec = node_spec
      )
      helper_name <- local_helper_name
      actual_helper_name <- helper_name
      justification <- as.character(allow$allowed_reason_code)
      eval <- .ltmle_exact_eval_allowed_local_source_prediction(
        task = task,
        branch_state = boundary_state,
        task_fit_cache = task_fit_cache,
        models = models,
        treatment_models = treatment_models,
        censoring_models = censoring_models,
        T = T,
        spec = spec,
        probability_bounds = probability_bounds,
        treat_mech = treat_mech,
        p_rct = p_rct,
        node_spec = node_spec,
        ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n
      )
      eval$metadata <- eval$metadata %||% list()
      eval$metadata$source_eval_mode <- "local_qstar_source_prediction"
      eval$metadata$branch_state_source_boundary_materialized <- FALSE
      eval$metadata$source_boundary_reached <- FALSE
      eval$metadata$source_boundary_task_id <- "not_applicable"
      eval$metadata$evaluated_beyond_source_boundary <- FALSE
      eval$metadata$terminal_task_evaluated_in_source_path <- FALSE
      eval$metadata$max_task_depth_beyond_source_boundary <- 0L
      eval$metadata$cached_targeted_continuation_used <- FALSE
      eval$metadata$continuation_task_id <- "not_applicable"
      eval$metadata$continuation_fit_found <- FALSE
      eval$metadata$continuation_fit_targeted <- FALSE
      eval$metadata$branch_state_downstream_recursion_used <- FALSE
      eval$metadata$terminal_full_recursion_used <- FALSE
      eval$metadata$terminal_outcome_base_case_used <- FALSE
      eval$metadata$local_qstar_source_prediction_used <- TRUE
      eval$metadata$local_current_task_prediction_used <- TRUE
      eval$metadata$used_observed_long_as_source_history <- FALSE
      eval$metadata$n_source_boundary_rows <- 0L
      eval$metadata$n_source_boundary_state_keys <- 0L
      eval$metadata$n_particles_before_collapse <- 0L
      eval$metadata$n_particles_after_collapse <- 0L
      eval$metadata$safe_collapse_applied <- FALSE
      eval$metadata$safe_collapse_key <- "not_applicable"
      source_eval_mode <- "local_qstar_source_prediction"
      downstream_recursion_used <- FALSE
      local_used <- TRUE
    }
  }

  if (identical(requires_branch, TRUE) && identical(helper_name, local_helper_name)) {
    .stop("Internal error: branch-required source pseudooutcome used local Q* prediction for task=", task_id)
  }
  if (isTRUE(source_task_is_terminal_outcome) && isTRUE(local_used)) {
    .stop("Internal error: terminal source task was handled as local prediction for task=", task_id)
  }

  if (length(eval$value) != nrow(observed_training_rows) || any(!is.finite(eval$value))) {
    .stop("Source pseudooutcome dispatcher returned invalid values for task=", task_id)
  }
  meta <- eval$metadata %||% list()
  source_eval_mode <- as.character(meta$source_eval_mode %||% source_eval_mode)
  downstream_recursion_used <- isTRUE(meta$branch_state_downstream_recursion_used) ||
    isTRUE(downstream_recursion_used)
  terminal_used <- isTRUE(meta$terminal_outcome_base_case_used) || isTRUE(terminal_used)
  local_used <- isTRUE(meta$local_qstar_source_prediction_used) || isTRUE(local_used)
  path_row <- .ltmle_exact_source_eval_training_diag_row(
    task = task,
    source_info = source_info0,
    source_task_id = source_task_id,
    source_eval_mode = source_eval_mode,
    branch_state_source_boundary_materialized = isTRUE(meta$branch_state_source_boundary_materialized),
    source_boundary_reached = isTRUE(meta$source_boundary_reached),
    source_boundary_task_id = meta$source_boundary_task_id %||% source_task_id,
    source_boundary_type = meta$source_boundary_type %||% "pure",
    source_boundary_direction = meta$source_boundary_direction %||% "not_applicable",
    source_boundary_eval_state = meta$source_boundary_eval_state %||% "not_applicable",
    source_boundary_outcome_history_state = meta$source_boundary_outcome_history_state %||% "not_applicable",
    source_boundary_m1_history_state = meta$source_boundary_m1_history_state %||% "not_applicable",
    source_boundary_m2_history_state = meta$source_boundary_m2_history_state %||% "not_applicable",
    source_boundary_auxiliary_mediator_history_state =
      meta$source_boundary_auxiliary_mediator_history_state %||% "not_applicable",
    qsource_predicted_on_mixed_boundary =
      isTRUE(meta$qsource_predicted_on_mixed_boundary),
    evaluated_beyond_source_boundary = isTRUE(meta$evaluated_beyond_source_boundary),
    terminal_task_evaluated_in_source_path = isTRUE(meta$terminal_task_evaluated_in_source_path),
    max_task_depth_beyond_source_boundary = as.integer(meta$max_task_depth_beyond_source_boundary %||% 0L),
    cached_targeted_continuation_used = isTRUE(meta$cached_targeted_continuation_used),
    continuation_task_id = meta$continuation_task_id %||% "not_applicable",
    continuation_fit_found = isTRUE(meta$continuation_fit_found),
    continuation_fit_targeted = isTRUE(meta$continuation_fit_targeted),
    branch_state_downstream_recursion_used = downstream_recursion_used,
    terminal_full_recursion_used = isTRUE(meta$terminal_full_recursion_used),
    terminal_outcome_base_case_used = terminal_used,
    local_qstar_source_prediction_used = local_used,
    local_current_task_prediction_used = isTRUE(meta$local_current_task_prediction_used),
    source_eval_helper_called = TRUE,
    helper_name = helper_name,
    actual_helper_name = as.character(meta$actual_helper_name %||% actual_helper_name),
    requires_branch = requires_branch,
    requires_branch_exception = requires_branch_exception,
    requires_branch_exception_reason = requires_branch_exception_reason,
    justification_for_local_prediction = justification,
    matched_allowed_local_source_prediction_map = matched_allowed_local_source_prediction_map,
    allowed_local_prediction = allowed_local_prediction,
    n_allowed_local_source_prediction_matches = n_allowed_local_source_prediction_matches,
    allowed_reason_code = allowed_reason_code,
    local_prediction_allowlist_audit_passed = TRUE,
    observed_training_rows_used_for_training = TRUE,
    observed_rows_mutated_to_assigned_A = FALSE,
    observed_long_used_as_generated_history = FALSE,
    used_observed_long_as_source_history = isTRUE(meta$used_observed_long_as_source_history),
    n_source_boundary_rows = as.integer(meta$n_source_boundary_rows %||% 0L),
    n_source_boundary_state_keys = as.integer(meta$n_source_boundary_state_keys %||% 0L),
    n_particles_before_collapse = as.integer(meta$n_particles_before_collapse %||% 0L),
    n_particles_after_collapse = as.integer(meta$n_particles_after_collapse %||% 0L),
    safe_collapse_applied = isTRUE(meta$safe_collapse_applied),
    safe_collapse_key = as.character(meta$safe_collapse_key %||% "not_applicable"),
    initial_dp_source_eval_available = isTRUE(meta$initial_dp_source_eval_available),
    initial_dp_source_eval_mean = as.numeric(meta$initial_dp_source_eval_mean %||% NA_real_),
    initial_dp_source_eval_source =
      as.character(meta$initial_dp_source_eval_source %||% "not_available"),
    targeted_dp_source_eval_mean = as.numeric(meta$targeted_dp_source_eval_mean %||% NA_real_),
    targeted_minus_initial_dp_source_eval =
      as.numeric(meta$targeted_minus_initial_dp_source_eval %||% NA_real_),
	    self_source_reference_allowed = self_allowed,
	    self_source_reference_reason = task$self_source_reference_reason %||% NA_character_,
		    virtual_mixed_direct_continuation_used =
		      isTRUE(meta$virtual_mixed_direct_continuation_used),
		    direct_continuation_allowed =
		      isTRUE(meta$direct_continuation_allowed),
	    downstream_source_can_represent_virtual_target =
	      isTRUE(meta$downstream_source_can_represent_virtual_target),
	    dedicated_virtual_Q_target_used =
	      isTRUE(meta$dedicated_virtual_Q_target_used),
	    virtual_mixed_downstream_source_task_id =
      meta$virtual_mixed_downstream_source_task_id %||% NA_character_,
    virtual_mixed_downstream_source_fit_found =
      isTRUE(meta$virtual_mixed_downstream_source_fit_found),
    virtual_mixed_downstream_source_fit_targeted =
      isTRUE(meta$virtual_mixed_downstream_source_fit_targeted),
    T = T,
    passed = TRUE,
    failure_class = "no_failure"
  )
  eval$source_eval_mode <- source_eval_mode
  eval$branch_state_downstream_recursion_used <- downstream_recursion_used
  eval$cached_targeted_continuation_used <- isTRUE(meta$cached_targeted_continuation_used)
  eval$terminal_full_recursion_used <- isTRUE(meta$terminal_full_recursion_used)
	  eval$terminal_outcome_base_case_used <- terminal_used
	  eval$local_qstar_source_prediction_used <- local_used
	  eval$local_current_task_prediction_used <- isTRUE(meta$local_current_task_prediction_used)
		  eval$virtual_mixed_direct_continuation_used <- isTRUE(meta$virtual_mixed_direct_continuation_used)
		  eval$direct_continuation_allowed <- isTRUE(meta$direct_continuation_allowed)
	  eval$downstream_source_can_represent_virtual_target <-
	    isTRUE(meta$downstream_source_can_represent_virtual_target)
	  eval$dedicated_virtual_Q_target_used <- isTRUE(meta$dedicated_virtual_Q_target_used)
	  eval$helper_name <- helper_name
  eval$actual_helper_name <- as.character(meta$actual_helper_name %||% actual_helper_name)
  eval$justification_for_local_prediction <- justification
  eval$source_task_id <- source_task_id
  if (is.null(eval$id)) eval$id <- seq_along(eval$value)
  if (is.null(eval$state_key)) {
    state_key <- if (isTRUE(source_task_is_terminal_outcome)) {
      "outcome"
    } else {
      .ltmle_exact_task_branch_key(source_info0$source_task)
    }
    eval$state_key <- rep(state_key, length(eval$value))
  }
  eval$diagnostics <- eval$diagnostics %||% list()
  eval$diagnostics$source_eval_training_path_check <- path_row
  eval
}

.ltmle_exact_source_pseudooutcome_independent_check_row <- function(task,
                                                                    observed_training_rows,
                                                                    stored_qstar_source,
                                                                    stored_id = NULL,
                                                                    stored_state_key = NULL,
                                                                    task_graph,
                                                                    task_fit_cache,
                                                                    models,
                                                                    treatment_models,
                                                                    censoring_models,
                                                                    T,
                                                                    spec,
                                                                    probability_bounds,
                                                                    treat_mech,
                                                                    p_rct,
	                                                                    node_spec = NULL,
	                                                                    ltmle_exact_density_ratio_mc_n = 2000L,
	                                                                    ltmle_exact_law_integration_n = 5L,
	                                                                    force_epsilon_zero = FALSE,
	                                                                    tolerance = 1e-8) {
  task <- as.list(task)
  task_id <- as.character(task$task_id)
  source_info <- .ltmle_exact_source_task_lookup(task, task_graph)
  source_task_id <- source_info$source_task_id
  source_task_is_terminal_outcome <- isTRUE(source_info$source_task_is_terminal_outcome)
  requires_branch <- isTRUE(source_info$requires_branch)

  if (is.na(source_task_id) || !nzchar(source_task_id)) {
    .stop("Missing source_task_id in independent source pseudooutcome check for task=", task_id)
  }
  if (identical(source_task_id, task_id) &&
      !isTRUE(task$source_self_reference_allowed %||% FALSE)) {
    .stop("Unexpected self-source pseudooutcome recursion for task=", task_id,
          ". Check observed_pseudooutcome_source_task_id.")
  }

  boundary_state <- .ltmle_exact_init_branch_state_from_observed_task(
    observed_task_rows = observed_training_rows,
    task = task,
    spec = spec,
    node_spec = node_spec
  )

	  if (.ltmle_exact_is_virtual_mixed_task(task) &&
	      !.ltmle_exact_virtual_mixed_can_use_downstream_cached_source(task)) {
	    independent_eval <- .ltmle_exact_eval_virtual_mixed_dedicated_L_continuation(
	      virtual_task = task,
	      branch_state = boundary_state,
	      task_graph = task_graph,
	      task_fit_cache = task_fit_cache,
	      models = models,
	      treatment_models = treatment_models,
	      censoring_models = censoring_models,
	      T = T,
	      spec = spec,
	      probability_bounds = probability_bounds,
	      treat_mech = treat_mech,
	      p_rct = p_rct,
	      node_spec = node_spec,
	      ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
	      ltmle_exact_law_integration_n = ltmle_exact_law_integration_n,
	      actual_helper_name =
	        ".ltmle_exact_independent_source_boundary_dp_continuation",
	      integration_context = "independent_virtual_mixed_dedicated_L_transition_target"
	    )
	    independent_eval$metadata$independent_source_boundary_materializer_applicable <- TRUE
    independent_eval$metadata$independent_source_boundary_materializer_completed <- TRUE
    independent_eval$metadata$independent_source_boundary_materializer_helper_name <-
      ".ltmle_exact_audit_materialize_source_task_boundary"
    independent_eval$metadata$same_boundary_materializer_recomputation_check <- FALSE
  } else if (isTRUE(source_task_is_terminal_outcome)) {
    independent_eval <- .ltmle_exact_eval_terminal_source_base_case(
      branch_state = boundary_state,
      task = task,
      spec = spec,
      node_spec = node_spec
    )
    independent_eval$metadata$source_eval_mode <- "branch_state_terminal_outcome_base_case"
    independent_eval$metadata$actual_helper_name <- ".ltmle_exact_independent_terminal_source_base_case"
    independent_eval$metadata$independent_source_boundary_materializer_applicable <- FALSE
    independent_eval$metadata$independent_source_boundary_materializer_completed <- TRUE
    independent_eval$metadata$independent_source_boundary_materializer_helper_name <-
      ".ltmle_exact_independent_terminal_source_base_case"
    independent_eval$metadata$same_boundary_materializer_recomputation_check <- FALSE
  } else {
    .ltmle_exact_validate_source_task_for_cached_continuation(
      task = task,
      source_task_id = source_task_id,
      task_graph = task_graph,
      task_fit_cache = task_fit_cache
    )
    source_boundary <- .ltmle_exact_audit_materialize_source_task_boundary(
      branch_state = boundary_state,
      task = task,
      source_task_id = source_task_id,
      task_graph = task_graph,
      models = models,
      treatment_models = treatment_models,
      censoring_models = censoring_models,
      T = T,
      spec = spec,
      probability_bounds = probability_bounds,
      treat_mech = treat_mech,
      p_rct = p_rct,
      node_spec = node_spec,
      ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
      ltmle_exact_law_integration_n = ltmle_exact_law_integration_n,
      forbid_terminal_full_recursion = TRUE
    )
    if (.ltmle_exact_is_virtual_mixed_task(source_boundary$source_task) &&
        .ltmle_exact_virtual_mixed_can_use_downstream_cached_source(source_boundary$source_task)) {
      independent_eval <- .ltmle_exact_independent_eval_virtual_mixed_cached_continuation(
        virtual_task = source_boundary$source_task,
        branch_state = source_boundary$branch_state,
        task_graph = task_graph,
        task_fit_cache = task_fit_cache,
        models = models,
        treatment_models = treatment_models,
        censoring_models = censoring_models,
        T = T,
        spec = spec,
        probability_bounds = probability_bounds,
        treat_mech = treat_mech,
        p_rct = p_rct,
        node_spec = node_spec,
        ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
        ltmle_exact_law_integration_n = ltmle_exact_law_integration_n
      )
	      independent_eval$metadata <- c(
	        source_boundary$diagnostics,
	        independent_eval$metadata %||% list(),
	        list(actual_helper_name =
	               ".ltmle_exact_independent_source_boundary_dp_continuation")
	      )
		    } else if (.ltmle_exact_is_virtual_mixed_task(source_boundary$source_task)) {
	      independent_eval <- .ltmle_exact_eval_virtual_mixed_task_fit_on_source_boundary(
	        virtual_task = source_boundary$source_task,
	        consuming_task = task,
	        source_boundary = source_boundary,
	        task_fit_cache = task_fit_cache,
	        models = models,
	        treatment_models = treatment_models,
	        censoring_models = censoring_models,
	        T = T,
        spec = spec,
        probability_bounds = probability_bounds,
        treat_mech = treat_mech,
	        p_rct = p_rct,
	        node_spec = node_spec,
	        ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
	        actual_helper_name =
	          ".ltmle_exact_independent_source_boundary_dp_continuation"
	      )
    } else {
      q_source_star <- .ltmle_exact_predict_task_qstar_on_rows(
        task_fit = task_fit_cache[[source_task_id]],
        rows = source_boundary$rows,
        row_long = source_boundary$long,
        models = models,
        treatment_models = treatment_models,
        censoring_models = censoring_models,
        T = T,
        spec = spec,
        probability_bounds = probability_bounds,
        treat_mech = treat_mech,
        p_rct = p_rct,
        node_spec = node_spec,
        ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n
      )
      independent_eval <- .ltmle_exact_collapse_source_qstar_to_subjects(
        q_source_star = q_source_star,
        source_boundary = source_boundary,
        task = task,
        source_task_id = source_task_id
      )
    }
    independent_eval$metadata$actual_helper_name <-
      ".ltmle_exact_independent_source_boundary_dp_continuation"
    independent_eval$metadata$independent_source_boundary_materializer_applicable <- TRUE
    independent_eval$metadata$independent_source_boundary_materializer_completed <- TRUE
    independent_eval$metadata$independent_source_boundary_materializer_helper_name <-
      ".ltmle_exact_audit_materialize_source_task_boundary"
    independent_eval$metadata$same_boundary_materializer_recomputation_check <- FALSE
  }

  meta <- independent_eval$metadata %||% list()
  branch_used <- isTRUE(meta$branch_state_downstream_recursion_used)
  terminal_used <- isTRUE(meta$terminal_outcome_base_case_used)
  local_used <- isTRUE(meta$local_qstar_source_prediction_used)
  actual_helper_name <- as.character(meta$actual_helper_name %||%
    ".ltmle_exact_independent_source_boundary_dp_continuation")
  cached_continuation_used <- isTRUE(meta$cached_targeted_continuation_used)
  continuation_task_id <- as.character(meta$continuation_task_id %||% if (source_task_is_terminal_outcome) {
    "not_applicable"
  } else {
    source_task_id
  })
  continuation_fit_found <- isTRUE(meta$continuation_fit_found)
  continuation_fit_targeted <- isTRUE(meta$continuation_fit_targeted)
  virtual_mixed_direct_continuation_used <- isTRUE(
    meta$virtual_mixed_direct_continuation_used
  )
  direct_continuation_allowed <- isTRUE(meta$direct_continuation_allowed)
  downstream_source_can_represent_virtual_target <- isTRUE(
    meta$downstream_source_can_represent_virtual_target
  )
  dedicated_virtual_Q_target_used <- isTRUE(meta$dedicated_virtual_Q_target_used)
  virtual_mixed_downstream_source_task_id <- as.character(
    meta$virtual_mixed_downstream_source_task_id %||% NA_character_
  )
  virtual_mixed_downstream_source_fit_found <- isTRUE(
    meta$virtual_mixed_downstream_source_fit_found
  )
  virtual_mixed_downstream_source_fit_targeted <- isTRUE(
    meta$virtual_mixed_downstream_source_fit_targeted
  )
  terminal_full_recursion_used <- isTRUE(meta$terminal_full_recursion_used)
  independent_materializer_applicable <- isTRUE(
    meta$independent_source_boundary_materializer_applicable %||%
      !isTRUE(source_task_is_terminal_outcome)
  )
  independent_materializer_completed <- isTRUE(
    meta$independent_source_boundary_materializer_completed %||%
      isTRUE(source_task_is_terminal_outcome)
  )
  independent_materializer_helper_name <- as.character(
    meta$independent_source_boundary_materializer_helper_name %||%
      if (isTRUE(source_task_is_terminal_outcome)) {
        ".ltmle_exact_independent_terminal_source_base_case"
      } else {
        ".ltmle_exact_audit_materialize_source_task_boundary"
      }
  )
  same_boundary_materializer_recomputation_check <- isTRUE(
    meta$same_boundary_materializer_recomputation_check %||% FALSE
  )
  independent_eval_mode <- if (isTRUE(source_task_is_terminal_outcome)) {
    "independent_branch_state_terminal_outcome_base_case"
  } else {
    "independent_branch_state_dp_continuation_value"
  }

  allowed_path <- if (identical(source_task_id, "observed_terminal_outcome")) {
    terminal_used && !branch_used && !local_used && !terminal_full_recursion_used
  } else {
    cached_continuation_used &&
      continuation_fit_found &&
      continuation_fit_targeted &&
      independent_materializer_completed &&
      !same_boundary_materializer_recomputation_check &&
      !branch_used &&
      !terminal_used &&
      !local_used &&
      !terminal_full_recursion_used
  }

  base_fail_row <- function(failure_class) {
    data.frame(
      component = task$component %||% NA_character_,
      task_id = task_id,
      t = as.integer(task$t %||% NA_integer_),
      node = task$node %||% NA_character_,
      process_type = .ltmle_exact_process_type(task),
      requires_branch = requires_branch,
      source_task_id = source_task_id,
      continuation_task_id = continuation_task_id,
      source_eval_mode = if (source_task_is_terminal_outcome) {
        "branch_state_terminal_outcome_base_case"
      } else {
        "branch_state_dp_continuation_value"
      },
      independent_eval_mode = independent_eval_mode,
      independent_eval_helper_name = actual_helper_name,
      independent_source_boundary_materializer_applicable = independent_materializer_applicable,
      independent_source_boundary_materializer_completed = independent_materializer_completed,
      independent_source_boundary_materializer_helper_name = independent_materializer_helper_name,
      same_boundary_materializer_recomputation_check = same_boundary_materializer_recomputation_check,
      cached_targeted_continuation_used = cached_continuation_used,
      continuation_fit_found = continuation_fit_found,
      continuation_fit_targeted = continuation_fit_targeted,
      virtual_mixed_direct_continuation_used =
        virtual_mixed_direct_continuation_used,
      direct_continuation_allowed = direct_continuation_allowed,
      downstream_source_can_represent_virtual_target =
        downstream_source_can_represent_virtual_target,
      dedicated_virtual_Q_target_used = dedicated_virtual_Q_target_used,
      virtual_mixed_downstream_source_task_id =
        virtual_mixed_downstream_source_task_id,
      virtual_mixed_downstream_source_fit_found =
        virtual_mixed_downstream_source_fit_found,
      virtual_mixed_downstream_source_fit_targeted =
        virtual_mixed_downstream_source_fit_targeted,
      terminal_full_recursion_used = terminal_full_recursion_used,
      stored_qstar_source_mean = NA_real_,
      independent_branch_recursion_qstar_source_mean = NA_real_,
      difference = NA_real_,
      max_abs_subject_difference = NA_real_,
      mean_abs_subject_difference = NA_real_,
      n_subjects_compared = 0L,
      tolerance = tolerance,
      branch_state_downstream_recursion_used = branch_used,
      terminal_outcome_base_case_used = terminal_used,
      local_qstar_source_prediction_used = local_used,
      used_stored_qstar_as_reference = FALSE,
      used_subject_matrix_as_reference = FALSE,
      used_production_subject_matrix_as_reference = FALSE,
      used_returned_mean_as_reference = FALSE,
      used_same_dispatcher_path_as_reference = FALSE,
      id_alignment_checked = FALSE,
      state_key_alignment_checked = FALSE,
      n_ids_compared = 0L,
      n_state_keys_compared = 0L,
      passed = FALSE,
      failure_class = failure_class,
      stringsAsFactors = FALSE
    )
  }

  if (identical(requires_branch, TRUE) && !allowed_path) {
    return(base_fail_row("independent_source_pseudooutcome_recursion_not_implemented"))
  }

  stored_tbl <- data.frame(
    id = stored_id %||% observed_training_rows$id,
    state_key = stored_state_key %||% rep(if (source_task_is_terminal_outcome) {
      "outcome"
    } else {
      .ltmle_exact_task_branch_key(source_info$source_task)
    }, length(stored_qstar_source)),
    stored = as.numeric(stored_qstar_source),
    stringsAsFactors = FALSE
  )
  independent_tbl <- data.frame(
    id = independent_eval$id,
    state_key = independent_eval$state_key,
    independent = as.numeric(independent_eval$value),
    stringsAsFactors = FALSE
  )
  cmp <- merge(stored_tbl, independent_tbl, by = c("id", "state_key"), all = FALSE, sort = FALSE)
  id_alignment_checked <- nrow(cmp) == nrow(stored_tbl) && nrow(cmp) == nrow(independent_tbl)
  state_key_alignment_checked <- id_alignment_checked

  if (!id_alignment_checked || !state_key_alignment_checked) {
    out <- base_fail_row("source_pseudooutcome_subject_alignment_failed")
    out$stored_qstar_source_mean <- mean(stored_tbl$stored, na.rm = TRUE)
    out$independent_branch_recursion_qstar_source_mean <- mean(independent_tbl$independent, na.rm = TRUE)
    out$n_subjects_compared <- nrow(cmp)
    out$id_alignment_checked <- id_alignment_checked
    out$state_key_alignment_checked <- state_key_alignment_checked
    out$n_ids_compared <- length(unique(cmp$id))
    out$n_state_keys_compared <- length(unique(cmp$state_key))
    return(out)
  }

  diff_subject <- cmp$stored - cmp$independent
  max_abs <- max(abs(diff_subject), na.rm = TRUE)
  mean_abs <- mean(abs(diff_subject), na.rm = TRUE)
  mean_diff <- mean(cmp$stored, na.rm = TRUE) - mean(cmp$independent, na.rm = TRUE)

  passed <- is.finite(max_abs) && max_abs <= tolerance &&
    is.finite(mean_abs) && mean_abs <= tolerance &&
    is.finite(mean_diff) && abs(mean_diff) <= tolerance &&
    id_alignment_checked && state_key_alignment_checked &&
    allowed_path && !local_used

  failure_class <- if (passed) {
    "no_failure"
  } else {
    "source_pseudooutcome_independent_recursion_mismatch"
  }
  data.frame(
    component = task$component %||% NA_character_,
    task_id = task_id,
    t = as.integer(task$t %||% NA_integer_),
    node = task$node %||% NA_character_,
    process_type = .ltmle_exact_process_type(task),
    requires_branch = requires_branch,
    source_task_id = source_task_id,
    continuation_task_id = continuation_task_id,
    source_eval_mode = if (source_task_is_terminal_outcome) {
      "branch_state_terminal_outcome_base_case"
    } else {
      "branch_state_dp_continuation_value"
    },
    independent_eval_mode = independent_eval_mode,
    independent_eval_helper_name = actual_helper_name,
    independent_source_boundary_materializer_applicable = independent_materializer_applicable,
    independent_source_boundary_materializer_completed = independent_materializer_completed,
    independent_source_boundary_materializer_helper_name = independent_materializer_helper_name,
    same_boundary_materializer_recomputation_check = same_boundary_materializer_recomputation_check,
    cached_targeted_continuation_used = cached_continuation_used,
    continuation_fit_found = continuation_fit_found,
    continuation_fit_targeted = continuation_fit_targeted,
    virtual_mixed_direct_continuation_used =
      virtual_mixed_direct_continuation_used,
    direct_continuation_allowed = direct_continuation_allowed,
    downstream_source_can_represent_virtual_target =
      downstream_source_can_represent_virtual_target,
    dedicated_virtual_Q_target_used = dedicated_virtual_Q_target_used,
    virtual_mixed_downstream_source_task_id =
      virtual_mixed_downstream_source_task_id,
    virtual_mixed_downstream_source_fit_found =
      virtual_mixed_downstream_source_fit_found,
    virtual_mixed_downstream_source_fit_targeted =
      virtual_mixed_downstream_source_fit_targeted,
    terminal_full_recursion_used = terminal_full_recursion_used,
    stored_qstar_source_mean = mean(cmp$stored, na.rm = TRUE),
    independent_branch_recursion_qstar_source_mean = mean(cmp$independent, na.rm = TRUE),
    difference = mean_diff,
    max_abs_subject_difference = max_abs,
    mean_abs_subject_difference = mean_abs,
    n_subjects_compared = nrow(cmp),
    tolerance = tolerance,
    branch_state_downstream_recursion_used = branch_used,
    terminal_outcome_base_case_used = terminal_used,
    local_qstar_source_prediction_used = local_used,
    used_stored_qstar_as_reference = FALSE,
    used_subject_matrix_as_reference = FALSE,
    used_production_subject_matrix_as_reference = FALSE,
    used_returned_mean_as_reference = FALSE,
    used_same_dispatcher_path_as_reference = FALSE,
    id_alignment_checked = id_alignment_checked,
    state_key_alignment_checked = state_key_alignment_checked,
    n_ids_compared = length(unique(cmp$id)),
    n_state_keys_compared = length(unique(cmp$state_key)),
    passed = passed,
    failure_class = failure_class,
    stringsAsFactors = FALSE
  )
}

.ltmle_exact_observed_vs_source_row_role <- function(task,
                                                     row_role,
                                                     uses_task_data_H_observed_data,
                                                     uses_observed_long,
                                                     uses_branch_state,
                                                     uses_branch_state_long,
                                                     conditioning_boundary,
                                                     passed = TRUE,
                                                     failure_class = "no_failure") {
  task <- as.list(task)
  data.frame(
    component = task$component %||% NA_character_,
    task_id = task$task_id %||% NA_character_,
    t = as.integer(task$t %||% NA_integer_),
    node = task$node %||% NA_character_,
    process_type = .ltmle_exact_process_type(task),
    row_role = row_role,
    uses_task_data_H_observed_data = isTRUE(uses_task_data_H_observed_data),
    uses_observed_long = isTRUE(uses_observed_long),
    mutates_observed_A_to_assigned_A = FALSE,
    uses_branch_state = isTRUE(uses_branch_state),
    uses_branch_state_long = isTRUE(uses_branch_state_long),
    uses_observed_long_as_generated_history = FALSE,
    conditioning_boundary = conditioning_boundary,
    passed = isTRUE(passed),
    failure_class = if (isTRUE(passed)) "no_failure" else failure_class,
    stringsAsFactors = FALSE
  )
}

.ltmle_exact_eval_source_pseudooutcome_from_boundary <- function(task,
                                                                observed_training_rows,
                                                                task_graph,
                                                                task_fit_cache,
                                                                models,
                                                                treatment_models,
                                                                censoring_models,
                                                                T,
                                                                spec,
                                                                probability_bounds,
                                                                treat_mech,
                                                                p_rct,
                                                                node_spec = NULL,
                                                                ltmle_exact_density_ratio_mc_n = 2000L,
                                                                ltmle_exact_law_integration_n = 5L) {
  task <- as.list(task)
  source_task_id <- as.character(task$observed_pseudooutcome_source_task_id)
  source_info <- .ltmle_exact_source_task_lookup(task, task_graph)
  boundary_state <- .ltmle_exact_init_branch_state_from_observed_task(
    observed_task_rows = observed_training_rows,
    task = task,
    spec = spec,
    node_spec = node_spec
  )
  if (isTRUE(source_info$source_task_is_terminal_outcome)) {
    eval <- .ltmle_exact_eval_terminal_source_base_case(
      branch_state = boundary_state,
      task = task,
      spec = spec,
      node_spec = node_spec
    )
  } else {
    .ltmle_exact_validate_source_task_for_cached_continuation(
      task = task,
      source_task_id = source_task_id,
      task_graph = task_graph,
      task_fit_cache = task_fit_cache
    )
    source_boundary <- .ltmle_exact_materialize_source_task_boundary(
      branch_state = boundary_state,
      task = task,
      source_task_id = source_task_id,
      task_graph = task_graph,
      models = models,
      treatment_models = treatment_models,
      censoring_models = censoring_models,
      T = T,
      spec = spec,
      probability_bounds = probability_bounds,
      treat_mech = treat_mech,
      p_rct = p_rct,
      node_spec = node_spec,
      ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
      ltmle_exact_law_integration_n = ltmle_exact_law_integration_n,
      forbid_terminal_full_recursion = TRUE
    )
    q_source_star <- .ltmle_exact_predict_task_qstar_on_rows(
      task_fit = task_fit_cache[[source_task_id]],
      rows = source_boundary$rows,
      row_long = source_boundary$long,
      models = models,
      treatment_models = treatment_models,
      censoring_models = censoring_models,
      T = T,
      spec = spec,
      probability_bounds = probability_bounds,
      treat_mech = treat_mech,
      p_rct = p_rct,
      node_spec = node_spec,
      ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n
    )
    eval <- .ltmle_exact_collapse_source_qstar_to_subjects(
      q_source_star = q_source_star,
      source_boundary = source_boundary,
      task = task,
      source_task_id = source_task_id
    )
  }
  mode <- if (identical(source_task_id, "observed_terminal_outcome")) {
    "branch_state_terminal_outcome_base_case"
  } else {
    "branch_state_dp_continuation_value"
  }
  branch_used <- FALSE
  terminal_used <- identical(source_task_id, "observed_terminal_outcome")
  meta <- eval$metadata %||% list()
  list(
    value = as.numeric(eval$value),
    id = eval$id %||% seq_along(eval$value),
    state_key = eval$state_key %||% rep(if (terminal_used) "outcome" else .ltmle_exact_task_branch_key(source_info$source_task), length(eval$value)),
    source_task_id = source_task_id,
    boundary_state = boundary_state,
    source_eval_mode = mode,
    source_eval_training_path_check = .ltmle_exact_source_eval_training_diag_row(
      task = task,
      source_info = source_info,
      source_task_id = source_task_id,
      source_eval_mode = mode,
      branch_state_source_boundary_materialized = isTRUE(meta$branch_state_source_boundary_materialized),
      source_boundary_reached = isTRUE(meta$source_boundary_reached),
      source_boundary_task_id = meta$source_boundary_task_id %||% source_task_id,
      evaluated_beyond_source_boundary = isTRUE(meta$evaluated_beyond_source_boundary),
      terminal_task_evaluated_in_source_path = isTRUE(meta$terminal_task_evaluated_in_source_path),
      max_task_depth_beyond_source_boundary = as.integer(meta$max_task_depth_beyond_source_boundary %||% 0L),
      cached_targeted_continuation_used = isTRUE(meta$cached_targeted_continuation_used),
      continuation_task_id = meta$continuation_task_id %||% "not_applicable",
      continuation_fit_found = isTRUE(meta$continuation_fit_found),
      continuation_fit_targeted = isTRUE(meta$continuation_fit_targeted),
      branch_state_downstream_recursion_used = branch_used,
      terminal_full_recursion_used = isTRUE(meta$terminal_full_recursion_used),
      terminal_outcome_base_case_used = terminal_used,
      helper_name = if (terminal_used) ".ltmle_exact_eval_terminal_source_base_case" else ".ltmle_exact_eval_source_pseudooutcome_dp",
      actual_helper_name = meta$actual_helper_name %||% if (terminal_used) ".ltmle_exact_eval_terminal_source_base_case" else ".ltmle_exact_eval_source_pseudooutcome_dp",
      requires_branch = isTRUE(source_info$requires_branch),
      requires_branch_exception = isTRUE(source_info$requires_branch_exception),
      requires_branch_exception_reason = if (terminal_used) "terminal_source_task" else "not_applicable",
      justification_for_local_prediction = if (terminal_used) "terminal_source_task" else "not_applicable",
      allowed_reason_code = "not_applicable",
      used_observed_long_as_source_history = FALSE,
      n_source_boundary_rows = as.integer(meta$n_source_boundary_rows %||% 0L),
      n_source_boundary_state_keys = as.integer(meta$n_source_boundary_state_keys %||% 0L),
      n_particles_before_collapse = as.integer(meta$n_particles_before_collapse %||% 0L),
      n_particles_after_collapse = as.integer(meta$n_particles_after_collapse %||% 0L),
      safe_collapse_applied = isTRUE(meta$safe_collapse_applied),
      safe_collapse_key = meta$safe_collapse_key %||% "not_applicable",
      T = T
    ),
    observed_vs_source_row_role_check = .ltmle_exact_observed_vs_source_row_role(
      task = task,
      row_role = "source_evaluation",
      uses_task_data_H_observed_data = FALSE,
      uses_observed_long = FALSE,
      uses_branch_state = TRUE,
      uses_branch_state_long = TRUE,
      conditioning_boundary = "observed_task"
    ),
    branch_trace = eval$branch_trace %||% data.frame(),
    handoff_trace = eval$handoff_trace %||% data.frame(),
    diagnostics = eval$diagnostics %||% list(
      density_ratio_diagnostics = data.frame(),
      law_integration_diagnostics = data.frame()
    )
  )
}

.ltmle_exact_independent_Qroot_star_eval <- function(
  component,
  root_task_id,
  root_rows,
  task_graph,
  task_fit_cache,
  T,
  spec,
  treatment_models,
  censoring_models,
  probability_bounds,
  treat_mech,
  p_rct,
  node_spec = NULL,
  ltmle_exact_law_integration_n = 5L
) {
  force(task_graph)
  models <- attr(task_fit_cache, "models")
  if (is.null(models)) {
    .stop("Independent Qroot* eval requires nuisance models attached to task_fit_cache.")
  }
  density_ratio_mc_n <- attr(task_fit_cache, "ltmle_exact_density_ratio_mc_n") %||% 2000L
  independent <- .ltmle_exact_eval_targeted_branch_recursion(
    component = component,
    root_task_id = root_task_id,
    root_rows = root_rows,
    task_graph = task_graph,
    task_fit_cache = task_fit_cache,
    models = models,
    treatment_models = treatment_models,
    censoring_models = censoring_models,
    T = T,
    spec = spec,
    probability_bounds = probability_bounds,
    treat_mech = treat_mech,
    p_rct = p_rct,
    node_spec = node_spec,
    ltmle_exact_density_ratio_mc_n = density_ratio_mc_n,
    ltmle_exact_law_integration_n = ltmle_exact_law_integration_n
  )
  list(
    component = component,
    root_task_id = root_task_id,
    subject_q = independent$subject_q,
    mean = independent$mean,
    branch_trace = independent$branch_trace %||% data.frame(),
    handoff_trace = independent$handoff_trace %||% data.frame(),
    diagnostics = independent$diagnostics %||% list(),
    same_targeted_Q_state_used = TRUE,
    used_subject_matrix_as_reference = FALSE,
    used_root_plugin_diagnostics_as_reference = FALSE,
    used_returned_mean_as_reference = FALSE
  )
}

.ltmle_exact_terminal_root_plugin_branch_state <- function(dat_wide,
                                                           root_rows,
                                                           T,
                                                           spec,
                                                           models,
                                                           treatment_models,
                                                           censoring_models,
                                                           root_task,
                                                           root_task_id,
                                                           task_graph,
                                                           task_fit_cache,
                                                           probability_bounds,
                                                           treat_mech,
                                                           p_rct,
                                                           node_spec = NULL,
                                                           ltmle_exact_density_ratio_mc_n = 2000L,
                                                           ltmle_exact_law_integration_n = 5L,
                                                           verbose = FALSE) {
  force(dat_wide)
  force(ltmle_exact_law_integration_n)
  .ltmle_exact_assert_supported_deterministic_terminal_root(root_task, node_spec)
  root_fit <- task_fit_cache[[root_task_id]]
  if (is.null(root_fit)) .stop("Missing root task fit for P_n Q_root^*: ", root_task_id)
  root_state <- .ltmle_exact_init_branch_state_from_root(
    root_rows = root_rows,
    spec = spec,
    T = T,
    node_spec = node_spec
  )
  root_eval <- .ltmle_exact_materialize_root_eval(
    branch_state = root_state,
    root_task = root_task,
    spec = spec,
    node_spec = node_spec
  )
  q0_root <- .ltmle_exact_predict_continuation(root_fit$cont_fit, root_eval$rows)
  H_root <- .ltmle_exact_clever_covariate_for_rows(
    task = root_task,
    rows = root_eval$rows,
    models = models,
    treatment_models = treatment_models,
    censoring_models = censoring_models,
    T = T,
    spec = spec,
    probability_bounds = probability_bounds,
    treat_mech = treat_mech,
    p_rct = p_rct,
    node_spec = node_spec,
    row_long = root_eval$long,
    row_source = "generated",
    ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
    max_t = root_task$t
  )
  H_root <- .ltmle_exact_truncate_clever_covariate(
    H_root,
    truncation = root_fit$truncation,
    task = root_task,
    estimator_variant = root_fit$estimator_variant %||% "ltmle_exact_quantile_truncated",
    row_evaluation_context = "terminal_root"
  )
  qstar_reconstructed <- .ltmle_exact_apply_fluctuation(
    q0_new = q0_root,
    H_new = H_root,
    epsilon = root_fit$epsilon,
    bounds = root_fit$bounds
  )
  q_root_star <- qstar_reconstructed
  qstar_stored <- NULL
  if (!is.null(root_fit$predict_qstar) && is.function(root_fit$predict_qstar)) {
    qstar_stored <- root_fit$predict_qstar(root_eval$rows, root_eval$long)
    q_root_star <- as.numeric(qstar_stored)
  }
  if (length(q_root_star) != nrow(root_eval$rows) || any(!is.finite(q_root_star))) {
    .stop("P_n Q_root^* produced invalid root predictions for component: ", spec$component[1L])
  }
  if (length(q0_root) != nrow(root_eval$rows) || any(!is.finite(q0_root))) {
    .stop("P_n Q_root initial prediction produced invalid root predictions for component: ", spec$component[1L])
  }

  subject_q <- .ltmle_exact_collapse_particles_to_subjects(q_root_star, root_state, root_eval$state_key)
  subject_q_initial <- .ltmle_exact_collapse_particles_to_subjects(q0_root, root_state, root_eval$state_key)
  H_finite <- as.numeric(H_root)[is.finite(H_root)]
  density_diag <- attr(H_root, "density_ratio_diagnostics")
  if (!is.null(density_diag) && nrow(density_diag)) {
    density_diag$plugin_batch <- NA_integer_
    density_diag$terminal_plugin_type <- "deterministic_root"
    density_diag$row_evaluation_context <- "terminal_root"
  }
  max_abs_reconstructed_minus_stored <- if (!is.null(qstar_stored)) {
    max(abs(as.numeric(qstar_reconstructed) - as.numeric(qstar_stored)), na.rm = TRUE)
  } else {
    NA_real_
  }
  double_fluctuation_detected <- is.finite(max_abs_reconstructed_minus_stored) &&
    max_abs_reconstructed_minus_stored > 1e-8
  root_plugin_diag <- data.frame(
    component = spec$component[1L],
    root_plugin_evaluation_mode = "Pn_Q_root_star",
    root_task_id = root_task_id,
    root_fit_found = TRUE,
    root_rows_n = nrow(root_eval$rows),
    root_long_rows_n = nrow(root_eval$long),
    root_q0_mean = mean(q0_root, na.rm = TRUE),
    root_qstar_mean = mean(q_root_star, na.rm = TRUE),
    root_H_mean = if (length(H_finite)) mean(H_finite, na.rm = TRUE) else NA_real_,
    root_H_max_abs = if (length(H_finite)) max(abs(H_finite), na.rm = TRUE) else NA_real_,
    used_law_recursive_terminal_eval_for_final = FALSE,
    root_qstar_reconstructed_mean = mean(qstar_reconstructed, na.rm = TRUE),
    root_qstar_stored_fit_mean = if (!is.null(qstar_stored)) mean(qstar_stored, na.rm = TRUE) else NA_real_,
    max_abs_reconstructed_minus_stored = max_abs_reconstructed_minus_stored,
    double_fluctuation_detected = double_fluctuation_detected,
    passed = !double_fluctuation_detected,
    failure_class = if (!double_fluctuation_detected) "no_failure" else "double_fluctuation_detected",
    stringsAsFactors = FALSE
  )
  list(
    subject_q = subject_q,
    subject_q_initial = subject_q_initial,
    initial_mean = mean(subject_q_initial, na.rm = TRUE),
    targeted_mean = mean(subject_q, na.rm = TRUE),
    mean_targeted_minus_initial = mean(subject_q - subject_q_initial, na.rm = TRUE),
    terminal_plugin_type = "deterministic_root",
    terminal_plugin_mc_active = FALSE,
    effective_mc_n = 1L,
    n_batches = 1L,
    mc_integration_diagnostics = data.frame(
      component = spec$component[1L],
      terminal_plugin_type = "deterministic_root",
      terminal_plugin_mc_active = FALSE,
      effective_mc_n = 1L,
      n_batches = 1L,
      subject_rows = length(subject_q),
      batch_mean_min = mean(subject_q, na.rm = TRUE),
      batch_mean_max = mean(subject_q, na.rm = TRUE),
      batch_mean_sd = 0,
      mc_standard_error = 0,
      stringsAsFactors = FALSE
    ),
    density_ratio_diagnostics = density_diag %||% data.frame(),
    truncation_diagnostics = attr(H_root, "truncation_diagnostics") %||% data.frame(),
    clever_covariate_decomposition_diagnostics =
      attr(H_root, "clever_covariate_decomposition_diagnostics") %||% data.frame(),
    law_integration_diagnostics = data.frame(),
    branch_trace = data.frame(),
    handoff_trace = data.frame(),
    root_plugin_diagnostics = root_plugin_diag,
    observed_vs_source_row_role_check = .ltmle_exact_observed_vs_source_row_role(
      task = root_task,
      row_role = "root_plugin",
      uses_task_data_H_observed_data = FALSE,
      uses_observed_long = FALSE,
      uses_branch_state = TRUE,
      uses_branch_state_long = TRUE,
      conditioning_boundary = "root"
    ),
    H_new_mean = root_plugin_diag$root_H_mean[1L],
    H_new_max = root_plugin_diag$root_H_max_abs[1L]
  )
}

.ltmle_exact_integrate_M2_law_task <- function(task,
                                               source_branch_state_before_M2,
                                               child_task_fit,
                                               task_fit_cache,
                                               task_graph,
                                               models,
                                               treatment_models,
                                               censoring_models,
                                               T,
                                               spec,
                                               probability_bounds,
                                               treat_mech,
                                               p_rct,
                                               node_spec = NULL,
                                               ltmle_exact_density_ratio_mc_n = 2000L,
                                               ltmle_exact_law_integration_n = 5L) {
  task <- as.list(task)
  if (!identical(as.character(task$node), "M2")) {
    .stop("M2 law integration called for non-M2 task: ", task$task_id)
  }

  tt <- as.integer(task$t)
  if (tt < 2L) {
    rows0 <- .ltmle_exact_branch_state_rows(source_branch_state_before_M2, task, spec, node_spec)
    return(rep(1, nrow(rows0)))
  }

  hist <- .ltmle_exact_branch_state_rows(
    branch_state = source_branch_state_before_M2,
    task = task,
    spec = spec,
    node_spec = node_spec
  )

  mu <- .ltmle_exact_predict_m2(models$M2[[tt]], hist)
  sigma <- .ltmle_exact_predict_mediator_sigma(models$M2[[tt]], hist)

  grid <- .ltmle_exact_normal_quantile_grid(ltmle_exact_law_integration_n)
  child_task <- as.list(child_task_fit$task)
  n0 <- nrow(hist)
  ng <- length(grid$z)
  max_expanded_rows <- 200000L
  if (n0 * ng > max_expanded_rows) {
    chunk_size <- max(1L, floor(max_expanded_rows / ng))
    acc <- numeric(n0)
    starts <- seq.int(1L, n0, by = chunk_size)
    for (start in starts) {
      idx <- seq.int(start, min(n0, start + chunk_size - 1L))
      state_chunk <- source_branch_state_before_M2
      for (key in .ltmle_exact_branch_keys()) {
        state_chunk[[key]]$rows <- state_chunk[[key]]$rows[idx, , drop = FALSE]
        state_chunk[[key]]$rows$id <- seq_along(idx)
        state_chunk[[key]]$long <- state_chunk[[key]]$long[
          state_chunk[[key]]$long$id %in% source_branch_state_before_M2[[key]]$rows$id[idx],
          ,
          drop = FALSE
        ]
        if (nrow(state_chunk[[key]]$long)) {
          old_ids <- source_branch_state_before_M2[[key]]$rows$id[idx]
          state_chunk[[key]]$long$id <- match(state_chunk[[key]]$long$id, old_ids)
        }
      }
      acc[idx] <- .ltmle_exact_integrate_M2_law_task(
        task = task,
        source_branch_state_before_M2 = state_chunk,
        child_task_fit = child_task_fit,
        task_fit_cache = task_fit_cache,
        task_graph = task_graph,
        models = models,
        treatment_models = treatment_models,
        censoring_models = censoring_models,
        T = T,
        spec = spec,
        probability_bounds = probability_bounds,
        treat_mech = treat_mech,
        p_rct = p_rct,
        node_spec = node_spec,
        ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
        ltmle_exact_law_integration_n = ltmle_exact_law_integration_n
      )
    }
    return(acc)
  }

  row_idx <- rep(seq_len(n0), times = ng)
  grid_z <- rep(grid$z, each = n0)
  state_g <- .ltmle_exact_branch_state_expand(
    source_branch_state_before_M2,
    row_idx
  )
  state_g <- .ltmle_exact_branch_state_update_mediator(
    branch_state = state_g,
    task = task,
    node = "M2",
    value = rep(mu, times = ng) + rep(sigma, times = ng) * grid_z,
    spec = spec,
    node_spec = node_spec
  )
  child_eval <- .ltmle_exact_make_task_source_eval(
    task = child_task,
    observed_rows = state_g[[.ltmle_exact_task_branch_key(child_task)]]$rows,
    observed_long = state_g[[.ltmle_exact_task_branch_key(child_task)]]$long,
    spec = spec,
    node_spec = node_spec,
    branch_state = state_g
  )

  pred <- .ltmle_exact_predict_task_qstar_on_rows(
    task_fit = child_task_fit,
    rows = child_eval$rows,
    row_long = child_eval$long,
    models = models,
    treatment_models = treatment_models,
    censoring_models = censoring_models,
    T = T,
    spec = spec,
    probability_bounds = probability_bounds,
    treat_mech = treat_mech,
    p_rct = p_rct,
    node_spec = node_spec,
    ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n
  )

  acc <- as.numeric(matrix(pred, nrow = n0, ncol = ng) %*% grid$w)
  if (any(!is.finite(acc))) {
    .stop("Non-finite M2 law integrated pseudo-outcome for task=", task$task_id)
  }
  acc
}

.ltmle_exact_integrate_M1_law_task <- function(task,
                                               source_branch_state_before_M1,
                                               child_task_fit,
                                               task_fit_cache,
                                               task_graph,
                                               models,
                                               treatment_models,
                                               censoring_models,
                                               T,
                                               spec,
                                               probability_bounds,
                                               treat_mech,
                                               p_rct,
                                               node_spec = NULL,
                                               ltmle_exact_density_ratio_mc_n = 2000L,
                                               ltmle_exact_law_integration_n = 5L) {
  task <- as.list(task)
  if (!identical(as.character(task$node), "M1")) {
    .stop("M1 law integration called for non-M1 task: ", task$task_id)
  }

  tt <- as.integer(task$t)
  if (tt < 2L) {
    rows0 <- .ltmle_exact_branch_state_rows(source_branch_state_before_M1, task, spec, node_spec)
    return(rep(1, nrow(rows0)))
  }

  hist <- .ltmle_exact_branch_state_rows(
    branch_state = source_branch_state_before_M1,
    task = task,
    spec = spec,
    node_spec = node_spec
  )

  mu <- .ltmle_exact_predict_m1(models$M1[[tt]], hist)
  sigma <- .ltmle_exact_predict_mediator_sigma(models$M1[[tt]], hist)

  grid <- .ltmle_exact_normal_quantile_grid(ltmle_exact_law_integration_n)
  n0 <- nrow(hist)
  ng <- length(grid$z)
  child_task <- as.list(child_task_fit$task)
  row_idx <- rep(seq_len(n0), times = ng)
  grid_z <- rep(grid$z, each = n0)
  state_g <- .ltmle_exact_branch_state_expand(
    source_branch_state_before_M1,
    row_idx
  )
  state_g <- .ltmle_exact_branch_state_update_mediator(
    branch_state = state_g,
    task = task,
    node = "M1",
    value = rep(mu, times = ng) + rep(sigma, times = ng) * grid_z,
    spec = spec,
    node_spec = node_spec
  )

  if (identical(as.character(child_task$node), "M2") &&
      .ltmle_exact_process_type(child_task) %in% c(
        "joint_stochastic_mediator_intervention_law",
        "first_mediator_stochastic_intervention_law",
        "second_mediator_stochastic_intervention_law"
      )) {
    child_source_id <- child_task$observed_pseudooutcome_source_task_id
    child_child_fit <- task_fit_cache[[child_source_id]]
    if (is.null(child_child_fit)) {
      .stop("Missing child fit for nested M2 integration: ", child_source_id)
    }

    pred <- .ltmle_exact_integrate_M2_law_task(
      task = child_task,
      source_branch_state_before_M2 = state_g,
      child_task_fit = child_child_fit,
      task_fit_cache = task_fit_cache,
      task_graph = task_graph,
      models = models,
      treatment_models = treatment_models,
      censoring_models = censoring_models,
      T = T,
      spec = spec,
      probability_bounds = probability_bounds,
      treat_mech = treat_mech,
      p_rct = p_rct,
      node_spec = node_spec,
      ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
      ltmle_exact_law_integration_n = ltmle_exact_law_integration_n
    )
  } else {
    child_eval <- .ltmle_exact_make_task_source_eval(
      task = child_task,
      observed_rows = state_g[[.ltmle_exact_task_branch_key(child_task)]]$rows,
      observed_long = state_g[[.ltmle_exact_task_branch_key(child_task)]]$long,
      spec = spec,
      node_spec = node_spec,
      branch_state = state_g
    )

    pred <- .ltmle_exact_predict_task_qstar_on_rows(
      task_fit = child_task_fit,
      rows = child_eval$rows,
      row_long = child_eval$long,
      models = models,
      treatment_models = treatment_models,
      censoring_models = censoring_models,
      T = T,
      spec = spec,
      probability_bounds = probability_bounds,
      treat_mech = treat_mech,
      p_rct = p_rct,
      node_spec = node_spec,
      ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n
    )
  }

  acc <- as.numeric(matrix(pred, nrow = n0, ncol = ng) %*% grid$w)
  if (any(!is.finite(acc))) {
    .stop("Non-finite M1 law integrated pseudo-outcome for task=", task$task_id)
  }
  acc
}

.ltmle_exact_slim_task_graph <- function(task_graph) {
  task_graph$observed_values <- list()
  task_graph$generated_values <- list()
  task_graph
}

.ltmle_exact_baseline_frame <- function(dat_wide, node_spec = NULL) {
  cols <- unique(c("W1", "W2", "Y0", "M1_0", "M2_0", as.character(node_spec$baseline_vars %||% character(0))))
  cols <- intersect(cols, names(dat_wide))
  out <- dat_wide[, cols, drop = FALSE]
  for (nm in names(out)) out[[nm]] <- as.numeric(out[[nm]])
  out
}

.ltmle_exact_add_baseline_covariates <- function(row, baseline, node_spec = NULL) {
  base_cols <- unique(c("W1", "W2", "Y0", "M1_0", "M2_0", as.character(node_spec$baseline_vars %||% character(0))))
  for (nm in intersect(base_cols, names(baseline))) {
    if (!nm %in% names(row)) row[[nm]] <- as.numeric(baseline[[nm]])
  }
  row
}

.ltmle_exact_draw_L_for_recursion <- function(row, baseline, models, tt, node_spec = NULL) {
  row <- .ltmle_exact_add_baseline_covariates(row, baseline, node_spec)
  for (L_node in .ltmle_exact_L_nodes(node_spec)) {
    row <- .ltmle_exact_add_terms(row)
    mu <- .ltmle_exact_predict_nuis(models$L[[L_node]][[tt]], row, type = "numeric")
    sig <- .ltmle_exact_sigma(models$L[[L_node]][[tt]])
    row[[L_node]] <- mu + stats::rnorm(nrow(row), mean = 0, sd = sig)
    if (identical(L_node, "L")) row$L <- row[[L_node]]
  }
  .ltmle_exact_add_terms(row)
}

.ltmle_exact_simulate_full_paths_once <- function(baseline, models, regimen, node_spec = NULL) {
  n <- nrow(baseline)
  T <- length(regimen)
  L_nodes <- .ltmle_exact_L_nodes(node_spec)
  L_lag <- stats::setNames(vector("list", length(L_nodes)), L_nodes)
  for (L_node in L_nodes) L_lag[[L_node]] <- rep(0, n)
  Y_lag <- baseline$Y0
  M1_lag <- baseline$M1_0
  M2_lag <- baseline$M2_0
  M1_path <- matrix(NA_real_, n, T)
  M2_path <- matrix(NA_real_, n, T)
  L_path <- matrix(NA_real_, n, T)
  Y_path <- matrix(NA_real_, n, T)

  for (tt in seq_len(T)) {
    A_t <- rep(regimen[tt], n)
    if (tt == 1L) {
      M1_t <- baseline$M1_0
      M2_t <- baseline$M2_0
    } else {
      row_m1 <- data.frame(A = A_t, L_lag = L_lag[[L_nodes[1L]]], Y_lag = Y_lag,
                           M1_lag = M1_lag, M2_lag = M2_lag, stringsAsFactors = FALSE)
      row_m1 <- .ltmle_exact_add_baseline_covariates(row_m1, baseline, node_spec)
      for (L_node in L_nodes) row_m1[[.ltmle_exact_lag_name(L_node)]] <- L_lag[[L_node]]
      row_m1 <- .ltmle_exact_add_terms(row_m1)
      mu1 <- .ltmle_exact_predict_m1(models$M1[[tt]], row_m1)
      sig1 <- .ltmle_exact_predict_mediator_sigma(models$M1[[tt]], row_m1)
      M1_t <- mu1 + stats::rnorm(n, mean = 0, sd = sig1)

      row_m2 <- row_m1
      row_m2$M1 <- M1_t
      row_m2 <- .ltmle_exact_add_terms(row_m2)
      mu2 <- .ltmle_exact_predict_m2(models$M2[[tt]], row_m2)
      sig2 <- .ltmle_exact_predict_mediator_sigma(models$M2[[tt]], row_m2)
      M2_t <- mu2 + stats::rnorm(n, mean = 0, sd = sig2)
    }

    row_ly <- data.frame(A = A_t, L_lag = L_lag[[L_nodes[1L]]], Y_lag = Y_lag,
                         M1_lag = M1_lag, M2_lag = M2_lag,
                         M1 = M1_t, M2 = M2_t, stringsAsFactors = FALSE)
    for (L_node in L_nodes) row_ly[[.ltmle_exact_lag_name(L_node)]] <- L_lag[[L_node]]
    row_ly <- .ltmle_exact_draw_L_for_recursion(row_ly, baseline, models, tt, node_spec)
    muY <- .ltmle_exact_predict_nuis(models$Y[[tt]], row_ly, type = "numeric")
    sigY <- .ltmle_exact_sigma(models$Y[[tt]])
    Y_t <- muY + stats::rnorm(n, mean = 0, sd = sigY)

    M1_path[, tt] <- M1_t
    M2_path[, tt] <- M2_t
    L_path[, tt] <- row_ly[[L_nodes[1L]]]
    Y_path[, tt] <- Y_t
    M1_lag <- M1_t
    M2_lag <- M2_t
    for (L_node in L_nodes) L_lag[[L_node]] <- row_ly[[L_node]]
    Y_lag <- Y_t
  }
  list(M1 = M1_path, M2 = M2_path, L = L_path, Y = Y_path)
}

.ltmle_exact_simulate_LY_given_mediators_once <- function(baseline, models, outer_regimen,
                                                          M1_path, M2_path, node_spec = NULL) {
  n <- nrow(baseline)
  T <- length(outer_regimen)
  L_nodes <- .ltmle_exact_L_nodes(node_spec)
  L_lag <- stats::setNames(vector("list", length(L_nodes)), L_nodes)
  for (L_node in L_nodes) L_lag[[L_node]] <- rep(0, n)
  Y_lag <- baseline$Y0
  M1_lag <- baseline$M1_0
  M2_lag <- baseline$M2_0
  L_path <- matrix(NA_real_, n, T)
  Y_path <- matrix(NA_real_, n, T)

  for (tt in seq_len(T)) {
    row_ly <- data.frame(A = rep(outer_regimen[tt], n),
                         L_lag = L_lag[[L_nodes[1L]]],
                         Y_lag = Y_lag,
                         M1_lag = M1_lag,
                         M2_lag = M2_lag,
                         M1 = M1_path[, tt],
                         M2 = M2_path[, tt],
                         stringsAsFactors = FALSE)
    for (L_node in L_nodes) row_ly[[.ltmle_exact_lag_name(L_node)]] <- L_lag[[L_node]]
    row_ly <- .ltmle_exact_draw_L_for_recursion(row_ly, baseline, models, tt, node_spec)
    muY <- .ltmle_exact_predict_nuis(models$Y[[tt]], row_ly, type = "numeric")
    sigY <- .ltmle_exact_sigma(models$Y[[tt]])
    Y_t <- muY + stats::rnorm(n, mean = 0, sd = sigY)
    L_path[, tt] <- row_ly[[L_nodes[1L]]]
    Y_path[, tt] <- Y_t
    M1_lag <- M1_path[, tt]
    M2_lag <- M2_path[, tt]
    for (L_node in L_nodes) L_lag[[L_node]] <- row_ly[[L_node]]
    Y_lag <- Y_t
  }
  list(L = L_path, Y = Y_path)
}

.ltmle_exact_independent_L_training_row_oracle <- function(task,
                                                           rows,
                                                           training_response,
                                                           training_prediction,
                                                           models,
                                                           spec,
                                                           T,
                                                           node_spec = NULL,
                                                           B_mc = 50L,
                                                           seed = 202405L,
                                                           force_epsilon_zero = FALSE) {
  task <- as.list(task)
  if (!.ltmle_exact_dedicated_L_training_row_trace_enabled(task)) {
    return(data.frame())
  }
  rows <- as.data.frame(rows)
  n <- nrow(rows)
  if (!n) return(data.frame())
  B_mc <- as.integer(B_mc)
  if (!is.finite(B_mc) || B_mc < 1L) B_mc <- 50L
  start_t <- as.integer(task$t %||% NA_integer_)
  if (!is.finite(start_t) || start_t < 1L || start_t > T) return(data.frame())
  component <- as.character(task$component %||% NA_character_)
  source_fit_task_id <- as.character(task$task_id %||% NA_character_)
  virtual_task_id <- paste0(component, "::1::virtual_mixed_continuation_after_Y::to_outer_L_2")
  L_nodes <- .ltmle_exact_L_nodes(node_spec)
  first_L_node <- L_nodes[1L]
  baseline <- rows
  if (!"Y0" %in% names(baseline) && "Y_lag" %in% names(baseline)) {
    baseline$Y0 <- as.numeric(baseline$Y_lag)
  }
  if (!"M1_0" %in% names(baseline) && "M1_lag" %in% names(baseline)) {
    baseline$M1_0 <- as.numeric(baseline$M1_lag)
  }
  if (!"M2_0" %in% names(baseline) && "M2_lag" %in% names(baseline)) {
    baseline$M2_0 <- as.numeric(baseline$M2_lag)
  }
  numeric_col <- function(nm, default = 0) {
    out <- if (nm %in% names(rows)) {
      suppressWarnings(as.numeric(rows[[nm]]))
    } else {
      rep(default, n)
    }
    out[!is.finite(out)] <- default
    out
  }
  lag_for <- function(L_node) {
    lag_nm <- .ltmle_exact_lag_name(L_node)
    if (lag_nm %in% names(rows)) return(numeric_col(lag_nm))
    if (identical(L_node, "L") && "L_lag" %in% names(rows)) return(numeric_col("L_lag"))
    rep(0, n)
  }
  response <- as.numeric(training_response)
  prediction <- as.numeric(training_prediction)
  if (length(response) != n) response <- rep(NA_real_, n)
  if (length(prediction) != n) prediction <- rep(NA_real_, n)
  oracle_sum <- rep(0, n)
  oracle_sumsq <- rep(0, n)
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(as.integer(seed)[1L])
  for (bb in seq_len(B_mc)) {
    L_lag <- stats::setNames(vector("list", length(L_nodes)), L_nodes)
    for (L_node in L_nodes) L_lag[[L_node]] <- lag_for(L_node)
    Y_lag <- numeric_col("Y_lag")
    M1_lag <- numeric_col("M1_lag")
    M2_lag <- numeric_col("M2_lag")
    Y_t <- rep(NA_real_, n)
    for (tt in seq.int(start_t, T)) {
      if (identical(tt, start_t)) {
        A_for_LY <- numeric_col("A", default = spec$outer_A[tt])
        M1_t <- numeric_col("M1")
        M2_t <- numeric_col("M2")
      } else {
        row_m1 <- data.frame(
          A = rep(spec$m1_A[tt], n),
          L_lag = L_lag[[first_L_node]],
          Y_lag = Y_lag,
          M1_lag = M1_lag,
          M2_lag = M2_lag,
          stringsAsFactors = FALSE
        )
        for (L_node in L_nodes) row_m1[[.ltmle_exact_lag_name(L_node)]] <- L_lag[[L_node]]
        row_m1 <- .ltmle_exact_add_baseline_covariates(row_m1, baseline, node_spec)
        row_m1 <- .ltmle_exact_add_terms(row_m1)
        mu1 <- .ltmle_exact_predict_m1(models$M1[[tt]], row_m1)
        sig1 <- .ltmle_exact_predict_mediator_sigma(models$M1[[tt]], row_m1)
        M1_t <- mu1 + stats::rnorm(n, mean = 0, sd = sig1)

        row_m2 <- data.frame(
          A = rep(spec$m2_A[tt], n),
          L_lag = L_lag[[first_L_node]],
          Y_lag = Y_lag,
          M1_lag = M1_lag,
          M2_lag = M2_lag,
          M1 = M1_t,
          stringsAsFactors = FALSE
        )
        for (L_node in L_nodes) row_m2[[.ltmle_exact_lag_name(L_node)]] <- L_lag[[L_node]]
        row_m2 <- .ltmle_exact_add_baseline_covariates(row_m2, baseline, node_spec)
        row_m2 <- .ltmle_exact_add_terms(row_m2)
        mu2 <- .ltmle_exact_predict_m2(models$M2[[tt]], row_m2)
        sig2 <- .ltmle_exact_predict_mediator_sigma(models$M2[[tt]], row_m2)
        M2_t <- mu2 + stats::rnorm(n, mean = 0, sd = sig2)
        A_for_LY <- rep(spec$outer_A[tt], n)
      }
      row_ly <- data.frame(
        A = A_for_LY,
        L_lag = L_lag[[first_L_node]],
        Y_lag = Y_lag,
        M1_lag = M1_lag,
        M2_lag = M2_lag,
        M1 = M1_t,
        M2 = M2_t,
        stringsAsFactors = FALSE
      )
      for (L_node in L_nodes) row_ly[[.ltmle_exact_lag_name(L_node)]] <- L_lag[[L_node]]
      row_ly <- .ltmle_exact_draw_L_for_recursion(row_ly, baseline, models, tt, node_spec)
      muY <- .ltmle_exact_predict_nuis(models$Y[[tt]], row_ly, type = "numeric")
      sigY <- .ltmle_exact_sigma(models$Y[[tt]])
      Y_t <- muY + stats::rnorm(n, mean = 0, sd = sigY)
      M1_lag <- M1_t
      M2_lag <- M2_t
      for (L_node in L_nodes) L_lag[[L_node]] <- row_ly[[L_node]]
      Y_lag <- Y_t
    }
    oracle_sum <- oracle_sum + Y_t
    oracle_sumsq <- oracle_sumsq + Y_t^2
  }
  oracle <- oracle_sum / B_mc
  oracle_var <- if (B_mc > 1L) {
    pmax((oracle_sumsq - (oracle_sum^2 / B_mc)) / (B_mc - 1L), 0)
  } else {
    rep(NA_real_, n)
  }
  oracle_mcse <- sqrt(oracle_var / B_mc)
  response_diff <- response - oracle
  prediction_diff <- prediction - oracle
  tol <- pmax(0.02, 3 * oracle_mcse)
  response_passed <- is.finite(response_diff) & is.finite(tol) &
    abs(response_diff) <= tol
  prediction_passed <- is.finite(prediction_diff) & is.finite(tol) &
    abs(prediction_diff) <= tol
  data.frame(
    component = component,
    source_fit_task_id = source_fit_task_id,
    virtual_task_id = virtual_task_id,
    t = start_t,
    node = as.character(task$node %||% NA_character_),
    process_type = .ltmle_exact_process_type(task),
    training_row_id = seq_len(n),
    subject_id = if ("id0" %in% names(rows)) as.integer(rows$id0) else seq_len(n),
    event_id = paste0(source_fit_task_id, "::training_row_", seq_len(n)),
    parent_particle_id = NA_character_,
    training_response = response,
    training_prediction = prediction,
    independent_training_row_oracle = oracle,
    independent_training_row_oracle_mcse = oracle_mcse,
    training_response_minus_oracle = response_diff,
    training_prediction_minus_oracle = prediction_diff,
    abs_training_response_minus_oracle = abs(response_diff),
    abs_training_prediction_minus_oracle = abs(prediction_diff),
    oracle_reference_source =
      "integrated_gcomp_forward_simulation_from_training_L_transition_boundary_row",
    oracle_uses_ltmle_root_initial_recursion = FALSE,
    oracle_uses_task_fit_cache = FALSE,
    oracle_uses_ratio_cache = FALSE,
    oracle_uses_targeting_H_attributes = FALSE,
    independent_oracle_used = TRUE,
    force_epsilon_zero = isTRUE(force_epsilon_zero),
    response_oracle_alignment_passed = response_passed,
    prediction_oracle_alignment_passed = prediction_passed,
    passed = response_passed & prediction_passed,
    failure_class = ifelse(
      response_passed & prediction_passed,
      "no_failure",
      ifelse(
        !response_passed,
        "training_response_not_aligned_with_integrated_node_oracle",
        "training_prediction_not_aligned_with_integrated_node_oracle"
      )
    ),
    stringsAsFactors = FALSE
  )
}

.ltmle_exact_integrated_worldmeans_from_models <- function(dat_wide, models, reg_a, reg_as, T,
                                                           node_spec = NULL,
                                                           B_mc = 2000L,
                                                           seed = 202405L) {
  # Auxiliary fitted-nuisance g-computation diagnostic. This is not the
  # ltmle_exact targeted substitution estimator and must not overwrite means.
  B_mc <- as.integer(B_mc)
  if (!is.finite(B_mc) || B_mc < 1L) B_mc <- 2000L
  baseline <- .ltmle_exact_baseline_frame(dat_wide, node_spec)
  n_subjects <- nrow(baseline)
  components <- component_mean_keys()
  out <- matrix(NA_real_, nrow = B_mc, ncol = length(component_mean_keys()))
  colnames(out) <- component_mean_keys()
  L_transition_sum <- matrix(0, nrow = n_subjects, ncol = length(components))
  L_transition_sumsq <- matrix(0, nrow = n_subjects, ncol = length(components))
  colnames(L_transition_sum) <- components
  colnames(L_transition_sumsq) <- components
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)

  for (bb in seq_len(B_mc)) {
    nat_a <- .ltmle_exact_simulate_full_paths_once(baseline, models, reg_a, node_spec)
    nat_as <- .ltmle_exact_simulate_full_paths_once(baseline, models, reg_as, node_spec)
    joint_med_aa <- .ltmle_exact_simulate_full_paths_once(baseline, models, reg_a, node_spec)
    joint_med_asas <- .ltmle_exact_simulate_full_paths_once(baseline, models, reg_as, node_spec)
    joint_aa <- .ltmle_exact_simulate_LY_given_mediators_once(
      baseline, models, reg_a, joint_med_aa$M1, joint_med_aa$M2, node_spec)
    joint_asas <- .ltmle_exact_simulate_LY_given_mediators_once(
      baseline, models, reg_as, joint_med_asas$M1, joint_med_asas$M2, node_spec)
    joint_aas <- .ltmle_exact_simulate_LY_given_mediators_once(
      baseline, models, reg_a, joint_med_asas$M1, joint_med_asas$M2, node_spec)
    sep_M1_a <- .ltmle_exact_simulate_full_paths_once(baseline, models, reg_a, node_spec)
    sep_M1_as <- .ltmle_exact_simulate_full_paths_once(baseline, models, reg_as, node_spec)
    sep_M2_a <- .ltmle_exact_simulate_full_paths_once(baseline, models, reg_a, node_spec)
    sep_M2_as <- .ltmle_exact_simulate_full_paths_once(baseline, models, reg_as, node_spec)
    sep_aaa <- .ltmle_exact_simulate_LY_given_mediators_once(
      baseline, models, reg_a, sep_M1_a$M1, sep_M2_a$M2, node_spec)
    sep_asas_asas <- .ltmle_exact_simulate_LY_given_mediators_once(
      baseline, models, reg_as, sep_M1_as$M1, sep_M2_as$M2, node_spec)
    sep_a_asas <- .ltmle_exact_simulate_LY_given_mediators_once(
      baseline, models, reg_a, sep_M1_as$M1, sep_M2_as$M2, node_spec)
    sep_a_aas <- .ltmle_exact_simulate_LY_given_mediators_once(
      baseline, models, reg_a, sep_M1_a$M1, sep_M2_as$M2, node_spec)
    component_subject_y <- cbind(
      mu_nat_a = nat_a$Y[, T],
      mu_nat_as = nat_as$Y[, T],
      mu_joint_aa = joint_aa$Y[, T],
      mu_joint_asas = joint_asas$Y[, T],
      mu_joint_aas = joint_aas$Y[, T],
      mu_sep_aaa = sep_aaa$Y[, T],
      mu_sep_asas_asas = sep_asas_asas$Y[, T],
      mu_sep_a_asas = sep_a_asas$Y[, T],
      mu_sep_a_aas = sep_a_aas$Y[, T]
    )
    out[bb, ] <- colMeans(component_subject_y)
    # Boundary-row continuation values for the L1 transition oracle. These are
    # subject/boundary-level MC averages from the integrated forward path.
    L_transition_sum <- L_transition_sum + component_subject_y[, components, drop = FALSE]
    L_transition_sumsq <- L_transition_sumsq + component_subject_y[, components, drop = FALSE]^2
  }
  L_transition_mean <- L_transition_sum / B_mc
  L_transition_var <- if (B_mc > 1L) {
    pmax((L_transition_sumsq - (L_transition_sum^2 / B_mc)) / (B_mc - 1L), 0)
  } else {
    matrix(NA_real_, nrow = n_subjects, ncol = length(components),
           dimnames = dimnames(L_transition_mean))
  }
  L_transition_mcse <- sqrt(L_transition_var / B_mc)
  L_transition_rows <- list()
  for (component in components) {
    task_id <- paste0(component, "::1::L::post_mediator_covariate_transition::main")
    L_transition_rows[[component]] <- data.frame(
      component = component,
      task_id = task_id,
      t = 1L,
      node = "L",
      process_type = "post_mediator_covariate_transition",
      boundary_row_id = seq_len(n_subjects),
      subject_id = seq_len(n_subjects),
      event_id = paste0(task_id, "::subject_", seq_len(n_subjects)),
      parent_particle_id = NA_character_,
      branch_state_id = paste0("integrated_gcomp_forward_boundary::", component, "::", seq_len(n_subjects)),
      integrated_node_level_value = as.numeric(L_transition_mean[, component]),
      integrated_node_level_mcse = as.numeric(L_transition_mcse[, component]),
      n_downstream_continuation_samples = B_mc,
      node_level_value_computation_method =
        "integrated_gcomp_forward_simulation_boundary_row_downstream_mean",
      reference_source = "integrated_gcomp_forward_simulation",
      reference_uses_ltmle_root_initial_recursion = FALSE,
      reference_uses_ltmle_continuation_cache = FALSE,
      reference_uses_task_fit_cache = FALSE,
      reference_uses_ratio_cache = FALSE,
      reference_uses_targeting_H_attributes = FALSE,
      reference_uses_ltmle_source_graph = FALSE,
      reference_uses_ltmle_materialized_rows = FALSE,
      independent_oracle_used = TRUE,
      stringsAsFactors = FALSE
    )
  }
  L_transition_node_values <- do.call(rbind, L_transition_rows)
  list(
    means = colMeans(out),
    mcse = apply(out, 2L, stats::sd) / sqrt(B_mc),
    L_transition_node_values = L_transition_node_values,
    B_mc = B_mc
  )
}

.ltmle_exact_component_fit <- function(dat_wide, long, T, spec, models, treatment_models,
                                       learner, sl_library, Q_model,
                                       seed, probability_bounds,
                                       truncation,
                                       estimator_variant,
                                       y_bounds_mode, y_bounds,
                                       score_tolerance, treat_mech, p_rct,
                                       censoring_models = NULL,
                                       component_tasks = NULL,
                                       node_spec = NULL,
                                       ltmle_exact_density_ratio_mc_n = 2000L,
                                       ltmle_exact_law_integration_n = 5L,
                                       auxiliary_nuisance_gcomp_mc_n = 2000L,
                                       force_epsilon_zero = FALSE,
                                       diagnostics_level = c("summary", "full"),
                                       verbose = FALSE) {
  diagnostics_level <- match.arg(diagnostics_level)

  comp <- spec$component[1L]
  .ltmle_exact_log(verbose, "[ltmle_exact] start component=", comp,
                   " terminal_plugin=deterministic_root",
                   " n=", nrow(dat_wide))

  observed_cache <- .ltmle_exact_build_ratio_cache(
    long = long,
    models = models,
    treatment_models = treatment_models,
    T = T,
    spec = spec,
    probability_bounds = probability_bounds,
    treat_mech = treat_mech,
    p_rct = p_rct,
    censoring_models = censoring_models,
    node_spec = node_spec,
    ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n
  )

  if (is.null(component_tasks) || !nrow(component_tasks)) {
    .stop("No factor tasks registered for ltmle_exact component: ", comp)
  }
  task_graph <- .ltmle_exact_build_task_graph(
    component_tasks = component_tasks,
    component_spec = spec,
    node_spec = node_spec,
    T = T
  )

  score_rows <- list()
  eif_rows <- list()
  density_ratio_rows <- list()
  truncation_rows <- list()
  clever_decomposition_rows <- list()
  terminal_subject_q <- NULL
  terminal_subject_q_initial <- NULL
  terminal_mc_diag <- NULL
  terminal_plugin_density_diag <- NULL
  component_D_accum <- numeric(nrow(dat_wide))
  task_fit_cache <- list()
  law_task_integration_rows <- list()
  branch_trace_rows <- list()
  handoff_trace_rows <- list()
  branch_law_integration_rows <- list()
  source_eval_training_rows <- list()
  cross_regimen_source_boundary_rows <- list()
  source_pseudooutcome_independent_rows <- list()
  targeted_continuation_storage_rows <- list()
  observed_vs_source_rows <- list()
  dedicated_L_training_row_trace_rows <- list()
  dedicated_L_prediction_row_trace_rows <- list()
  dedicated_L_coefficient_contribution_rows <- list()
  dedicated_L_prediction_weighting_collapse_rows <- list()
  dedicated_L_mediator_path_weight_trace_rows <- list()
  dedicated_L_training_response_oracle_rows <- list()
  root_plugin_rows <- list()
  root_vs_branch_rows <- list()
  separate_clever_covariate_identity_rows <- list()
  second_M2_marginal_reference_rows <- list()
  targeted_storage_rows <- list()

  if (identical(as.character(spec$world_type[1L]), "separate")) {
    second_M2_marginal_reference_rows[[1L]] <-
      .ltmle_exact_second_M2_marginal_reference_check(
        component = comp,
        component_tasks = component_tasks,
        long = long,
        models = models,
        T = T,
        spec = spec,
        production_mediator_cache = observed_cache$M,
        n_particles = ltmle_exact_density_ratio_mc_n,
        seed = 202405L,
        node_spec = node_spec
      )
  }

  for (task_id in task_graph$reverse_topological_order) {
    task <- task_graph$tasks[[task_id]]
    is_terminal_task <- identical(task_id, task_graph$terminal_task_id)
    is_virtual_mixed_task <- .ltmle_exact_is_virtual_mixed_task(task)
    .ltmle_exact_log(verbose,
                     "[ltmle_exact] component=", comp,
                     " task=", task_id,
                     " t=", task$t,
                     " node=", task$node,
                     " process=", .ltmle_exact_process_type(task),
                     " terminal_generated=", is_terminal_task)

    obs_current <- long[long$t == task$t, , drop = FALSE]
    obs_current <- obs_current[order(obs_current$id), , drop = FALSE]
    obs_for_boundary <- .ltmle_exact_add_terms(obs_current)
    qstar_source_eval <- .ltmle_exact_eval_source_pseudooutcome_dispatch(
      task = task,
      task_data = NULL,
      observed_training_rows = obs_for_boundary,
      observed_long = long,
      task_graph = task_graph,
      task_fit_cache = task_fit_cache,
      models = models,
      treatment_models = treatment_models,
      censoring_models = censoring_models,
      T = T,
      spec = spec,
      probability_bounds = probability_bounds,
      treat_mech = treat_mech,
      p_rct = p_rct,
	      node_spec = node_spec,
	      ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
	      ltmle_exact_law_integration_n = ltmle_exact_law_integration_n,
	      allow_local_source_prediction = FALSE,
	      force_epsilon_zero = force_epsilon_zero
	    )
    U_task <- as.numeric(qstar_source_eval$value)
    if (length(U_task) != nrow(obs_for_boundary) || any(!is.finite(U_task))) {
      .stop("Invalid DP source pseudooutcome training value for task=", task_id)
    }
    source_info_for_training <- .ltmle_exact_source_task_lookup(task, task_graph)
    source_is_terminal_for_training <- isTRUE(source_info_for_training$source_task_is_terminal_outcome)
    source_requires_branch_for_training <- isTRUE(source_info_for_training$requires_branch)
    training_outcome_source <- if (isTRUE(source_is_terminal_for_training)) {
      "observed_terminal_outcome_base_case"
    } else if (isTRUE(source_requires_branch_for_training)) {
      "branch_state_dp_source_eval"
    } else if (isTRUE(qstar_source_eval$local_qstar_source_prediction_used)) {
      "allowed_local_source_prediction"
    } else {
      "branch_state_dp_source_eval"
    }
    if (!isTRUE(source_is_terminal_for_training) && isTRUE(source_requires_branch_for_training)) {
      if (!identical(as.character(qstar_source_eval$source_eval_mode), "branch_state_dp_continuation_value") ||
          !isTRUE(qstar_source_eval$cached_targeted_continuation_used) ||
          isTRUE(qstar_source_eval$terminal_full_recursion_used) ||
          isTRUE(qstar_source_eval$local_qstar_source_prediction_used)) {
        .stop("task_fitting_order_violation: task=", task_id,
              ", source_task_id=", source_info_for_training$source_task_id)
      }
    }
    observed_vs_source_rows[[length(observed_vs_source_rows) + 1L]] <-
      .ltmle_exact_observed_vs_source_row_role(
        task = task,
        row_role = "observed_training",
        uses_task_data_H_observed_data = TRUE,
        uses_observed_long = TRUE,
        uses_branch_state = FALSE,
        uses_branch_state_long = FALSE,
        conditioning_boundary = "observed_data"
      )
    task_data <- .ltmle_exact_get_task_data(
      task = task,
      observed_long = long,
      generated_histories = NULL,
      Q_model = Q_model,
      spec = spec,
      node_spec = node_spec,
      U_observed = U_task,
      include_generated = FALSE
    )
    task_covariate_set_role <- .ltmle_exact_virtual_mixed_covariate_set_role(task, node_spec)

    cont_fit <- .ltmle_exact_fit_continuation(
      y = task_data$U_observed,
      xdf = task_data$H_observed_data,
      covs = task_data$conditioning_covariates,
      learner = learner,
      sl_library = sl_library,
      component = paste0(comp, " ", task$node, " t=", task$t)
    )
    design_matrix_columns <- .ltmle_exact_design_columns(cont_fit, task_data$H_observed_data)

    q0_obs <- .ltmle_exact_predict_continuation(cont_fit, task_data$H_observed_data)
    if (.ltmle_exact_dedicated_L_training_row_trace_enabled(task)) {
      dedicated_L_training_row_trace_rows[[length(dedicated_L_training_row_trace_rows) + 1L]] <-
        .ltmle_exact_dedicated_L_fit_training_row_trace(
          task = task,
          cont_fit = cont_fit,
          rows = task_data$H_observed_data,
          response = task_data$U_observed,
          q0_train = q0_obs
        )
      dedicated_L_training_response_oracle_rows[[length(dedicated_L_training_response_oracle_rows) + 1L]] <-
        .ltmle_exact_independent_L_training_row_oracle(
          task = task,
          rows = task_data$H_observed_data,
          training_response = task_data$U_observed,
          training_prediction = q0_obs,
          models = models,
          spec = spec,
          T = T,
          node_spec = node_spec,
          B_mc = auxiliary_nuisance_gcomp_mc_n,
          seed = seed + as.integer(task$t %||% 0L) * 104729L,
          force_epsilon_zero = force_epsilon_zero
        )
    }
    source_boundary_state <- .ltmle_exact_init_branch_state_from_observed_task(
      observed_task_rows = task_data$H_observed_data,
      task = task,
      spec = spec,
      node_spec = node_spec
    )
    source_eval <- .ltmle_exact_make_task_source_eval(
      task = task,
      observed_rows = task_data$H_observed_data,
      observed_long = NULL,
      spec = spec,
      node_spec = node_spec,
      branch_state = source_boundary_state
    )
    source_rows <- source_eval$rows
    source_long <- source_eval$long
    source_branch_state <- source_eval$branch_state
    source_state_key <- source_eval$state_key
    q0_source <- .ltmle_exact_predict_continuation(cont_fit, source_rows)

    H_obs_raw <- .ltmle_exact_clever_covariate(task, observed_cache)
    attr(H_obs_raw, "density_ratio_diagnostics") <-
      .ltmle_exact_density_ratio_diagnostics_for_task(
        task,
        observed_cache
      )
    H_obs <- .ltmle_exact_truncate_clever_covariate(
      H_obs_raw,
      truncation = truncation,
      task = task,
      estimator_variant = estimator_variant,
      row_evaluation_context = "observed"
    )
    if (identical(as.character(task$world_type %||% NA_character_), "separate")) {
      independent_H_reference <- .ltmle_exact_build_independent_separate_h_reference(
        task = task,
        long = long,
        models = models,
        treatment_models = treatment_models,
        censoring_models = censoring_models,
        T = T,
        spec = spec,
        probability_bounds = probability_bounds,
        treat_mech = treat_mech,
        p_rct = p_rct,
        node_spec = node_spec,
        n_particles = ltmle_exact_density_ratio_mc_n,
        seed = 202405L,
        max_t = task$t
      )
      separate_clever_covariate_identity_rows[[length(separate_clever_covariate_identity_rows) + 1L]] <-
        .ltmle_exact_clever_covariate_identity_row(
          task = task,
          H_used = H_obs_raw,
          H_truncated_used = H_obs,
          ratio_cache = observed_cache,
          reference_parts = independent_H_reference,
          T = T,
          spec = spec,
          row_evaluation_context = "observed",
          truncation = truncation,
          estimator_variant = estimator_variant
        )
    }
    truncation_rows[[length(truncation_rows) + 1L]] <-
      attr(H_obs, "truncation_diagnostics") %||% data.frame()
    clever_decomposition_rows[[length(clever_decomposition_rows) + 1L]] <-
      attr(H_obs, "clever_covariate_decomposition_diagnostics") %||% data.frame()

    obs_density_diag <- attr(H_obs, "density_ratio_diagnostics")
    obs_density_diag$plugin_batch <- NA_integer_
    obs_density_diag$terminal_plugin_type <- NA_character_
    obs_density_diag$row_evaluation_context <- "observed"
    density_ratio_rows[[length(density_ratio_rows) + 1L]] <- obs_density_diag

    H_source <- .ltmle_exact_clever_covariate_for_rows(
      task = task,
      rows = source_rows,
      models = models,
      treatment_models = treatment_models,
      censoring_models = censoring_models,
      T = T,
      spec = spec,
      probability_bounds = probability_bounds,
      treat_mech = treat_mech,
      p_rct = p_rct,
      node_spec = node_spec,
      row_long = source_long,
      row_source = "generated",
      ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
      max_t = task$t
    )
    H_source <- .ltmle_exact_truncate_clever_covariate(
      H_source,
      truncation = truncation,
      task = task,
      estimator_variant = estimator_variant,
      row_evaluation_context = "source"
    )
    truncation_rows[[length(truncation_rows) + 1L]] <-
      attr(H_source, "truncation_diagnostics") %||% data.frame()
    clever_decomposition_rows[[length(clever_decomposition_rows) + 1L]] <-
      attr(H_source, "clever_covariate_decomposition_diagnostics") %||% data.frame()
    if (any(!is.finite(H_source))) {
      .stop("Non-finite H_source for task=", task_id)
    }
    source_density_diag <- NULL
    source_density_diag <- attr(H_source, "density_ratio_diagnostics")
    if (!is.null(source_density_diag) && nrow(source_density_diag)) {
      source_density_diag$plugin_batch <- NA_integer_
      source_density_diag$terminal_plugin_type <- NA_character_
      source_density_diag$row_evaluation_context <- "source"
      density_ratio_rows[[length(density_ratio_rows) + 1L]] <- source_density_diag
    }

    upd <- .ltmle_exact_target_continuation(
      y_train = task_data$U_observed,
      q0_train = q0_obs,
      q0_new = NULL,
      H_obs = H_obs,
      H_new = NULL,
      y_bounds_mode = y_bounds_mode,
      y_bounds = y_bounds,
      learner = learner,
      node = task$node,
      component = comp,
      tt = task$t,
      process_type = .ltmle_exact_process_type(task),
      score_tolerance = score_tolerance,
      included_ratio_factors = .ltmle_exact_describe_ratio_factors(task, T),
      task_id = task$task_id %||% NA_character_,
      source_task_id = as.character(source_info_for_training$source_task_id %||% NA_character_),
      parent_task_id = task$parent_task_id %||% NA_character_,
      role = task$role %||% NA_character_,
      force_epsilon_zero = isTRUE(force_epsilon_zero) || isTRUE(is_virtual_mixed_task),
      compute_new = FALSE
    )

    H_source_finite <- H_source[is.finite(H_source)]
    upd$diagnostics$H_source_mean <- mean(H_source_finite, na.rm = TRUE)
    upd$diagnostics$H_source_max <- max(abs(H_source_finite), na.rm = TRUE)
    upd$diagnostics$H_obs_source_max_abs_diff <- max(abs(H_obs - H_source), na.rm = TRUE)
    upd$diagnostics$source_rows_evaluated <- TRUE
    upd$diagnostics$clever_covariate_truncation_enabled <- isTRUE(truncation$enabled)
    upd$diagnostics$clever_covariate_truncation_policy <- truncation$policy
    upd$diagnostics$clever_covariate_truncation_target <- truncation$target
    upd$diagnostics$clever_covariate_quantile_lower <- truncation$quantile_lower
    upd$diagnostics$clever_covariate_quantile_upper <- truncation$quantile_upper
    upd$diagnostics$density_ratio_factor_truncation_used <- FALSE
    upd$diagnostics$n_clever_covariate_values_truncated <- as.integer(
      (attr(H_obs, "n_clever_covariate_values_truncated") %||% 0L) +
        (attr(H_source, "n_clever_covariate_values_truncated") %||% 0L)
    )

    qstar_source_diagnostic <- .ltmle_exact_apply_fluctuation(
      q0_new = q0_source,
      H_new = H_source,
      epsilon = upd$epsilon,
      bounds = upd$bounds
    )
    upd$diagnostics$mean_q0_source <- mean(q0_source, na.rm = TRUE)
    upd$diagnostics$mean_qstar_source <- mean(qstar_source_diagnostic, na.rm = TRUE)
    upd$diagnostics$mean_qstar_source_minus_q0_source <- mean(qstar_source_diagnostic - q0_source, na.rm = TRUE)
    upd$diagnostics$max_abs_qstar_source_minus_q0_source <- max(abs(qstar_source_diagnostic - q0_source), na.rm = TRUE)
    semantic_cache_key <- .ltmle_exact_semantic_cache_key(task)
    target_label <- .ltmle_exact_continuation_target_label(task)
    upd$diagnostics$expected_cache_key <- semantic_cache_key
    upd$diagnostics$actual_cache_key <- semantic_cache_key
    upd$diagnostics$cache_key_matches <- TRUE
    upd$diagnostics$cache_key_target_label <- target_label
    upd$diagnostics$cache_key_component <- comp
    upd$diagnostics$cache_key_outer_regimen <- as.character(
      task$source_boundary_outer_regimen %||% NA_character_
    )
    upd$diagnostics$cache_key_m1_regimen <- as.character(
      task$source_boundary_m1_regimen %||% NA_character_
    )
    upd$diagnostics$cache_key_m2_regimen <- as.character(
      task$source_boundary_m2_regimen %||% NA_character_
    )

    task_fit_cache[[task_id]] <- list(
      task = task,
      task_id = task_id,
      component = comp,
      cache_key = semantic_cache_key,
      cont_fit = cont_fit,
      epsilon = upd$epsilon,
      bounds = upd$bounds,
      conditioning_covariates = task_data$conditioning_covariates,
      virtual_mixed_covariate_set_role = task_covariate_set_role,
      design_matrix_columns = design_matrix_columns,
      training_covariate_ranges = .ltmle_exact_covariate_ranges(
        task_data$H_observed_data,
        task_data$conditioning_covariates
      ),
      targeted = TRUE,
      virtual_mixed_task = isTRUE(is_virtual_mixed_task),
      target_label = target_label,
      source_boundary_outer_regimen =
        as.character(task$source_boundary_outer_regimen %||% NA_character_),
      source_boundary_m1_regimen =
        as.character(task$source_boundary_m1_regimen %||% NA_character_),
      source_boundary_m2_regimen =
        as.character(task$source_boundary_m2_regimen %||% NA_character_),
      source_boundary_outcome_history_state =
        as.character(task$source_boundary_outcome_history_state %||% NA_character_),
      source_boundary_m1_history_state =
        as.character(task$source_boundary_m1_history_state %||% NA_character_),
      source_boundary_m2_history_state =
        as.character(task$source_boundary_m2_history_state %||% NA_character_),
      source_boundary_auxiliary_mediator_history_state =
        as.character(task$source_boundary_auxiliary_mediator_history_state %||% NA_character_),
      truncation = truncation,
      estimator_variant = estimator_variant
    )

    if (isTRUE(qstar_source_eval$cached_targeted_continuation_used) ||
        isTRUE(qstar_source_eval$branch_state_downstream_recursion_used) ||
        isTRUE(qstar_source_eval$terminal_outcome_base_case_used)) {
      source_pseudooutcome_independent_rows[[length(source_pseudooutcome_independent_rows) + 1L]] <-
        .ltmle_exact_source_pseudooutcome_independent_check_row(
          task = task,
          observed_training_rows = task_data$H_observed_data,
          stored_qstar_source = qstar_source_eval$value,
          stored_id = qstar_source_eval$id %||% task_data$H_observed_data$id,
          stored_state_key = qstar_source_eval$state_key,
          task_graph = task_graph,
          task_fit_cache = task_fit_cache,
          models = models,
          treatment_models = treatment_models,
          censoring_models = censoring_models,
          T = T,
          spec = spec,
          probability_bounds = probability_bounds,
          treat_mech = treat_mech,
          p_rct = p_rct,
          node_spec = node_spec,
          ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
          ltmle_exact_law_integration_n = ltmle_exact_law_integration_n
        )
      invisible(gc(FALSE))
    }
    if (!is.null(qstar_source_eval$branch_trace) && nrow(qstar_source_eval$branch_trace)) {
      branch_trace_rows[[length(branch_trace_rows) + 1L]] <- qstar_source_eval$branch_trace
    }
    if (!is.null(qstar_source_eval$handoff_trace) && nrow(qstar_source_eval$handoff_trace)) {
      handoff_trace_rows[[length(handoff_trace_rows) + 1L]] <- qstar_source_eval$handoff_trace
    }
    if (!is.null(qstar_source_eval$diagnostics$density_ratio_diagnostics) &&
        nrow(qstar_source_eval$diagnostics$density_ratio_diagnostics)) {
      src_dr <- qstar_source_eval$diagnostics$density_ratio_diagnostics
      src_dr$row_evaluation_context <- "source_recursion"
      density_ratio_rows[[length(density_ratio_rows) + 1L]] <- src_dr
    }
    if (!is.null(qstar_source_eval$diagnostics$truncation_diagnostics) &&
        nrow(qstar_source_eval$diagnostics$truncation_diagnostics)) {
      src_tr <- qstar_source_eval$diagnostics$truncation_diagnostics
      src_tr$row_evaluation_context <- "source_recursion"
      truncation_rows[[length(truncation_rows) + 1L]] <- src_tr
    }
    if (!is.null(qstar_source_eval$diagnostics$clever_covariate_decomposition_diagnostics) &&
        nrow(qstar_source_eval$diagnostics$clever_covariate_decomposition_diagnostics)) {
      src_dec <- qstar_source_eval$diagnostics$clever_covariate_decomposition_diagnostics
      src_dec$row_evaluation_context <- "source_recursion"
      clever_decomposition_rows[[length(clever_decomposition_rows) + 1L]] <- src_dec
    }
    if (!is.null(qstar_source_eval$diagnostics$law_integration_diagnostics) &&
        nrow(qstar_source_eval$diagnostics$law_integration_diagnostics)) {
      branch_law_integration_rows[[length(branch_law_integration_rows) + 1L]] <-
        qstar_source_eval$diagnostics$law_integration_diagnostics
    }
    if (!is.null(qstar_source_eval$diagnostics$cross_regimen_source_boundary_trace) &&
        nrow(qstar_source_eval$diagnostics$cross_regimen_source_boundary_trace)) {
      cross_regimen_source_boundary_rows[[length(cross_regimen_source_boundary_rows) + 1L]] <-
        qstar_source_eval$diagnostics$cross_regimen_source_boundary_trace
    }
    if (!is.null(qstar_source_eval$diagnostics$dedicated_L_transition_L_fit_prediction_row_trace) &&
        nrow(qstar_source_eval$diagnostics$dedicated_L_transition_L_fit_prediction_row_trace)) {
      dedicated_L_prediction_row_trace_rows[[length(dedicated_L_prediction_row_trace_rows) + 1L]] <-
        qstar_source_eval$diagnostics$dedicated_L_transition_L_fit_prediction_row_trace
    }
    if (!is.null(qstar_source_eval$diagnostics$dedicated_L_transition_L_fit_coefficient_contribution_check) &&
        nrow(qstar_source_eval$diagnostics$dedicated_L_transition_L_fit_coefficient_contribution_check)) {
      dedicated_L_coefficient_contribution_rows[[length(dedicated_L_coefficient_contribution_rows) + 1L]] <-
        qstar_source_eval$diagnostics$dedicated_L_transition_L_fit_coefficient_contribution_check
    }
    if (!is.null(qstar_source_eval$diagnostics$dedicated_L_transition_prediction_weighting_collapse_check) &&
        nrow(qstar_source_eval$diagnostics$dedicated_L_transition_prediction_weighting_collapse_check)) {
      dedicated_L_prediction_weighting_collapse_rows[[length(dedicated_L_prediction_weighting_collapse_rows) + 1L]] <-
        qstar_source_eval$diagnostics$dedicated_L_transition_prediction_weighting_collapse_check
    }
    if (!is.null(qstar_source_eval$diagnostics$dedicated_L_transition_mediator_path_weight_trace) &&
        nrow(qstar_source_eval$diagnostics$dedicated_L_transition_mediator_path_weight_trace)) {
      dedicated_L_mediator_path_weight_trace_rows[[length(dedicated_L_mediator_path_weight_trace_rows) + 1L]] <-
        qstar_source_eval$diagnostics$dedicated_L_transition_mediator_path_weight_trace
    }
    observed_vs_source_rows[[length(observed_vs_source_rows) + 1L]] <-
      .ltmle_exact_observed_vs_source_row_role(
        task = task,
        row_role = "source_evaluation",
        uses_task_data_H_observed_data = FALSE,
        uses_observed_long = FALSE,
        uses_branch_state = TRUE,
        uses_branch_state_long = TRUE,
        conditioning_boundary = "observed_task"
      )

    if (is_terminal_task) {
      .ltmle_exact_assert_supported_deterministic_terminal_root(task, node_spec)

      root_rows <- .ltmle_exact_make_deterministic_root_plugin_data(
        dat_wide,
        spec,
        node_spec
      )$row

      # legacy diagnostic-only fallback remains available as
      # .ltmle_exact_terminal_root_plugin_deterministic(
      plugin <- .ltmle_exact_terminal_root_plugin_branch_state(
        dat_wide = dat_wide,
        root_rows = root_rows,
        T = T,
        spec = spec,
        models = models,
        treatment_models = treatment_models,
        censoring_models = censoring_models,
        root_task = task,
        root_task_id = task_id,
        task_graph = task_graph,
        task_fit_cache = task_fit_cache,
        probability_bounds = probability_bounds,
        treat_mech = treat_mech,
        p_rct = p_rct,
        node_spec = node_spec,
        ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
        ltmle_exact_law_integration_n = ltmle_exact_law_integration_n,
        verbose = verbose
      )

      terminal_subject_q <- plugin$subject_q
      terminal_subject_q_initial <- plugin$subject_q_initial
      terminal_mc_diag <- plugin$mc_integration_diagnostics
      terminal_plugin_density_diag <- plugin$density_ratio_diagnostics

      upd$diagnostics$generated_evaluated <- TRUE
      upd$diagnostics$H_new_mean <- plugin$H_new_mean
      upd$diagnostics$H_new_max <- plugin$H_new_max
      upd$diagnostics$terminal_plugin_batches <- plugin$n_batches
      upd$diagnostics$terminal_plugin_type <- plugin$terminal_plugin_type
      upd$diagnostics$terminal_plugin_mc_active <- plugin$terminal_plugin_mc_active
      upd$diagnostics$effective_mc_n <- plugin$effective_mc_n

      if (!is.null(terminal_plugin_density_diag) && nrow(terminal_plugin_density_diag)) {
        terminal_plugin_density_diag$row_evaluation_context <- "terminal_root"
        density_ratio_rows[[length(density_ratio_rows) + 1L]] <- terminal_plugin_density_diag
      }
      if (!is.null(plugin$truncation_diagnostics) && nrow(plugin$truncation_diagnostics)) {
        truncation_rows[[length(truncation_rows) + 1L]] <- plugin$truncation_diagnostics
      }
      if (!is.null(plugin$clever_covariate_decomposition_diagnostics) &&
          nrow(plugin$clever_covariate_decomposition_diagnostics)) {
        clever_decomposition_rows[[length(clever_decomposition_rows) + 1L]] <-
          plugin$clever_covariate_decomposition_diagnostics
      }
      if (!is.null(plugin$branch_trace) && nrow(plugin$branch_trace)) {
        branch_trace_rows[[length(branch_trace_rows) + 1L]] <- plugin$branch_trace
      }
      if (!is.null(plugin$handoff_trace) && nrow(plugin$handoff_trace)) {
        handoff_trace_rows[[length(handoff_trace_rows) + 1L]] <- plugin$handoff_trace
      }
      if (!is.null(plugin$law_integration_diagnostics) && nrow(plugin$law_integration_diagnostics)) {
        branch_law_integration_rows[[length(branch_law_integration_rows) + 1L]] <-
          plugin$law_integration_diagnostics
      }
      if (!is.null(plugin$root_plugin_diagnostics) && nrow(plugin$root_plugin_diagnostics)) {
        root_plugin_rows[[length(root_plugin_rows) + 1L]] <- plugin$root_plugin_diagnostics
      }
      root_returned_mean <- mean(plugin$subject_q, na.rm = TRUE)
      root_plugin_targeted_mean <- plugin$targeted_mean
      root_reference_mean <- mean(as.numeric(qstar_source_eval$value), na.rm = TRUE)
      root_reference_valid <- is.finite(root_reference_mean) &&
        length(qstar_source_eval$value) == nrow(obs_for_boundary) &&
        !isTRUE(qstar_source_eval$local_qstar_source_prediction_used) &&
        !isTRUE(qstar_source_eval$terminal_full_recursion_used)
      root_vs_branch_rows[[length(root_vs_branch_rows) + 1L]] <- data.frame(
        component = comp,
        root_task_id = task_id,
        returned_mean = root_returned_mean,
        root_plugin_targeted_mean = root_plugin_targeted_mean,
        independent_Qroot_star_mean = root_reference_mean,
        returned_minus_independent_Qroot_star = root_returned_mean - root_reference_mean,
        root_plugin_minus_independent_Qroot_star = root_plugin_targeted_mean - root_reference_mean,
        same_targeted_Q_state_used = TRUE,
        used_subject_matrix_as_reference = FALSE,
        used_root_plugin_diagnostics_as_reference = FALSE,
        used_returned_mean_as_reference = FALSE,
        reference_source = "root_training_source_dp_eval",
        reference_source_eval_mode = as.character(qstar_source_eval$source_eval_mode),
        cached_targeted_continuation_used = isTRUE(qstar_source_eval$cached_targeted_continuation_used),
        terminal_full_recursion_used = isTRUE(qstar_source_eval$terminal_full_recursion_used),
        local_qstar_source_prediction_used = isTRUE(qstar_source_eval$local_qstar_source_prediction_used),
        passed = root_reference_valid,
        failure_class = if (root_reference_valid) {
          "no_failure"
        } else {
          "root_training_source_dp_eval_invalid"
        },
        stringsAsFactors = FALSE
      )
      if (!is.null(plugin$observed_vs_source_row_role_check) &&
          nrow(plugin$observed_vs_source_row_role_check)) {
        observed_vs_source_rows[[length(observed_vs_source_rows) + 1L]] <-
          plugin$observed_vs_source_row_role_check
      }

      rm(plugin)
    }

    if (length(upd$train) != nrow(task_data$H_observed_data) ||
        any(!is.finite(upd$train))) {
      .stop("Invalid targeted continuation storage value for task=", task_id)
    }

    task_graph <- .ltmle_exact_store_targeted_continuation(
      task_graph = task_graph,
      task_id = task_id,
      qstar_observed = upd$train,
      qstar_generated = NULL,
      store_generated = FALSE
    )
    stored_observed_value <- as.numeric(task_graph$observed_values[[task_id]])
    if (length(stored_observed_value) != length(upd$train) ||
        any(!is.finite(stored_observed_value))) {
      .stop("Invalid stored targeted continuation value for task=", task_id)
    }
    stored_minus_upd_train <- stored_observed_value - as.numeric(upd$train)
    stored_minus_qstar_source_eval <- stored_observed_value - as.numeric(qstar_source_eval$value)
    training_outcome_minus_dp_source_eval <- as.numeric(task_data$U_observed) -
      as.numeric(qstar_source_eval$value)
    stored_equals_upd_train <- max(abs(stored_minus_upd_train), na.rm = TRUE) <= 1e-12
    stored_equals_qstar_source_eval <- max(abs(stored_minus_qstar_source_eval), na.rm = TRUE) <= 1e-12
    training_matches_dp_source <- max(abs(training_outcome_minus_dp_source_eval), na.rm = TRUE) <= 1e-12
    targeted_storage_passed <- isTRUE(stored_equals_upd_train) &&
      isTRUE(training_matches_dp_source) &&
      !isTRUE(stored_equals_qstar_source_eval && max(abs(upd$train - qstar_source_eval$value), na.rm = TRUE) > 1e-12)
    targeted_storage_rows[[length(targeted_storage_rows) + 1L]] <- data.frame(
      component = comp,
      task_id = task_id,
      source_task_id = as.character(source_info_for_training$source_task_id %||% NA_character_),
      t = as.integer(task$t %||% NA_integer_),
      node = task$node %||% NA_character_,
      process_type = .ltmle_exact_process_type(task),
      mean_upd_train = mean(upd$train, na.rm = TRUE),
      mean_qstar_source_eval = mean(qstar_source_eval$value, na.rm = TRUE),
      mean_stored_observed_value = mean(stored_observed_value, na.rm = TRUE),
      stored_value_source = "current_task_targeted_continuation",
      stored_minus_upd_train = mean(stored_minus_upd_train, na.rm = TRUE),
      stored_minus_qstar_source_eval = mean(stored_minus_qstar_source_eval, na.rm = TRUE),
      stored_equals_upd_train = isTRUE(stored_equals_upd_train),
      stored_equals_qstar_source_eval = isTRUE(stored_equals_qstar_source_eval),
      training_outcome_source = training_outcome_source,
      training_outcome_mean = mean(task_data$U_observed, na.rm = TRUE),
      training_outcome_hash = .ltmle_exact_numeric_vector_hash(task_data$U_observed),
      dp_source_eval_mean = mean(qstar_source_eval$value, na.rm = TRUE),
      dp_source_eval_hash = .ltmle_exact_numeric_vector_hash(qstar_source_eval$value),
      training_outcome_minus_dp_source_eval = mean(training_outcome_minus_dp_source_eval, na.rm = TRUE),
      virtual_mixed_covariate_set_role = task_covariate_set_role,
      conditioning_covariates = paste(task_data$conditioning_covariates, collapse = "|"),
      design_matrix_columns = paste(design_matrix_columns, collapse = "|"),
      conditioning_covariates_include_current_M1 = "M1" %in% task_data$conditioning_covariates,
      conditioning_covariates_include_current_M2 = "M2" %in% task_data$conditioning_covariates,
      design_matrix_has_current_M1 = any(design_matrix_columns %in% "M1"),
      design_matrix_has_current_M2 = any(design_matrix_columns %in% "M2"),
      legacy_graph_observed_value_used = FALSE,
      virtual_mixed_direct_continuation_used =
        isTRUE(qstar_source_eval$virtual_mixed_direct_continuation_used),
      direct_continuation_allowed =
        isTRUE(qstar_source_eval$direct_continuation_allowed),
      downstream_source_can_represent_virtual_target =
        isTRUE(qstar_source_eval$downstream_source_can_represent_virtual_target),
      dedicated_virtual_Q_target_used =
        isTRUE(qstar_source_eval$dedicated_virtual_Q_target_used),
      passed = targeted_storage_passed,
      failure_class = if (targeted_storage_passed) {
        "no_failure"
      } else if (!isTRUE(stored_equals_upd_train)) {
        "targeted_continuation_storage_mismatch"
      } else if (!isTRUE(training_matches_dp_source)) {
        "training_outcome_not_dp_source_eval"
      } else {
        "targeted_continuation_storage_check_failed"
      },
      stringsAsFactors = FALSE
    )

    if (!is.null(qstar_source_eval$diagnostics$source_eval_training_path_check) &&
        nrow(qstar_source_eval$diagnostics$source_eval_training_path_check)) {
      path_row <- qstar_source_eval$diagnostics$source_eval_training_path_check
      initial_value <- qstar_source_eval$initial_value %||% numeric(0)
      initial_available <- length(initial_value) == length(qstar_source_eval$value) &&
        length(initial_value) > 0L &&
        all(is.finite(initial_value))
      initial_mean <- if (isTRUE(initial_available)) {
        mean(initial_value, na.rm = TRUE)
      } else if ("initial_dp_source_eval_mean" %in% names(path_row)) {
        as.numeric(path_row$initial_dp_source_eval_mean[1L])
      } else {
        NA_real_
      }
      path_row$used_as_training_outcome <- TRUE
      path_row$training_outcome_source <- training_outcome_source
      path_row$training_outcome_mean <- mean(task_data$U_observed, na.rm = TRUE)
      path_row$training_outcome_hash <- .ltmle_exact_numeric_vector_hash(task_data$U_observed)
      path_row$dp_source_eval_mean <- mean(qstar_source_eval$value, na.rm = TRUE)
      path_row$dp_source_eval_hash <- .ltmle_exact_numeric_vector_hash(qstar_source_eval$value)
      path_row$initial_dp_source_eval_available <- isTRUE(initial_available) ||
        isTRUE(path_row$initial_dp_source_eval_available[1L] %||% FALSE)
      path_row$initial_dp_source_eval_mean <- initial_mean
      path_row$initial_dp_source_eval_hash <- if (isTRUE(initial_available)) {
        .ltmle_exact_numeric_vector_hash(initial_value)
      } else {
        NA_character_
      }
      path_row$initial_dp_source_eval_source <- if (isTRUE(initial_available)) {
        "same_target_source_with_targeting_update_epsilon_zero"
      } else {
        path_row$initial_dp_source_eval_source[1L] %||% "not_available"
      }
      path_row$targeted_dp_source_eval_mean <- mean(qstar_source_eval$value, na.rm = TRUE)
      path_row$targeted_dp_source_eval_hash <- .ltmle_exact_numeric_vector_hash(qstar_source_eval$value)
      path_row$targeted_minus_initial_dp_source_eval <- if (is.finite(initial_mean)) {
        mean(qstar_source_eval$value, na.rm = TRUE) - initial_mean
      } else {
        NA_real_
      }
      path_row$legacy_graph_observed_value_used <- FALSE
      path_row$legacy_graph_observed_value_mean <- 0
      path_row$training_outcome_minus_dp_source_eval <- mean(training_outcome_minus_dp_source_eval, na.rm = TRUE)
      path_row$stored_graph_value_quantity <- "current_task_targeted_Qstar_observed"
      path_row$stored_graph_value_mean <- mean(stored_observed_value, na.rm = TRUE)
      path_row$stored_graph_value_minus_current_qstar_observed <- mean(stored_minus_upd_train, na.rm = TRUE)
      path_row$virtual_mixed_covariate_set_role <- task_covariate_set_role
      path_row$conditioning_covariates <- paste(task_data$conditioning_covariates, collapse = "|")
      path_row$design_matrix_columns <- paste(design_matrix_columns, collapse = "|")
      path_row$conditioning_covariates_include_current_M1 <-
        "M1" %in% task_data$conditioning_covariates
      path_row$conditioning_covariates_include_current_M2 <-
        "M2" %in% task_data$conditioning_covariates
      path_row$design_matrix_has_current_M1 <- any(design_matrix_columns %in% "M1")
      path_row$design_matrix_has_current_M2 <- any(design_matrix_columns %in% "M2")
      source_eval_training_rows[[length(source_eval_training_rows) + 1L]] <- path_row
    }

    score_rows[[length(score_rows) + 1L]] <- upd$diagnostics

    term_rows <- .ltmle_exact_make_eif_term_rows(
      task, task_data, upd$D_task, H_obs, task_data$U_observed, upd$train
    )
    if (!isTRUE(is_virtual_mixed_task)) {
      component_D_accum <- component_D_accum + .ltmle_exact_subject_sum(term_rows$D_task, term_rows$id)
    }

    if (identical(diagnostics_level, "full") && !isTRUE(is_virtual_mixed_task)) {
      eif_rows[[length(eif_rows) + 1L]] <- term_rows
    }

    rm(term_rows, task_data, cont_fit, q0_obs, source_eval, source_rows, source_long,
       source_branch_state, source_state_key,
       q0_source, H_source, H_source_finite, source_density_diag, qstar_source_diagnostic,
       H_obs, obs_density_diag, upd, U_task, obs_current, obs_for_boundary,
       source_boundary_state, qstar_source_eval, source_info_for_training,
       stored_observed_value, stored_minus_upd_train, stored_minus_qstar_source_eval,
       training_outcome_minus_dp_source_eval, semantic_cache_key, target_label,
       task_covariate_set_role, design_matrix_columns)
    if (length(score_rows) %% 4L == 0L) invisible(gc(verbose = FALSE))
  }

  invisible(gc(verbose = FALSE))

  terminal_task_id <- task_graph$terminal_task_id
  if (is.null(terminal_subject_q) || length(terminal_subject_q) != nrow(dat_wide) ||
      any(!is.finite(terminal_subject_q))) {
    .stop("deterministic terminal/root plug-in did not produce finite subject_q for component: ", comp)
  }
  if (is.null(terminal_subject_q_initial) || length(terminal_subject_q_initial) != nrow(dat_wide) ||
      any(!is.finite(terminal_subject_q_initial))) {
    .stop("deterministic terminal/root plug-in did not produce finite initial subject_q for component: ", comp)
  }
  if (is.null(terminal_mc_diag)) {
    .stop("deterministic terminal/root plug-in did not produce diagnostics for component: ", comp)
  }

  subject_q <- terminal_subject_q
  mc_diag <- terminal_mc_diag
  task_graph_slim <- .ltmle_exact_slim_task_graph(task_graph)

  baseline_term <- subject_q - mean(subject_q, na.rm = TRUE)
  component_D_accum <- component_D_accum + baseline_term

  baseline_rows <- data.frame(
    id = seq_along(baseline_term),
    component = comp,
    task_id = paste(comp, "baseline", sep = "::"),
    parent_task_id = terminal_task_id,
    t = 0L,
    time = 0L,
    node = "baseline",
    process_type = "baseline_target_parameter",
    role = NA_character_,
    D_task = baseline_term,
    H = 1,
    U = subject_q,
    Q_star = mean(subject_q, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  if (identical(diagnostics_level, "full")) {
    eif_rows[[length(eif_rows) + 1L]] <- baseline_rows
    component_eif_terms <- do.call(rbind, eif_rows)
  } else {
    component_eif_terms <- NULL
  }

  component_D <- component_D_accum

  .ltmle_exact_log(verbose, "[ltmle_exact] finished component=", comp,
                   " mean=", signif(mean(subject_q, na.rm = TRUE), 6))

  list(
    mean = mean(subject_q, na.rm = TRUE),
    initial_mean = mean(terminal_subject_q_initial, na.rm = TRUE),
    subject_q = subject_q,
    subject_q_initial = terminal_subject_q_initial,
    mean_targeted_minus_initial = mean(terminal_subject_q - terminal_subject_q_initial, na.rm = TRUE),
    component_D = component_D,
    score_diagnostics = .ltmle_exact_rbind_fill(score_rows),
    component_eif_terms = component_eif_terms,
    density_ratio_diagnostics = .ltmle_exact_rbind_fill(density_ratio_rows),
    truncation_diagnostics = .ltmle_exact_rbind_fill(truncation_rows),
    clever_covariate_decomposition_diagnostics = .ltmle_exact_rbind_fill(clever_decomposition_rows),
    separate_clever_covariate_identity_check = .ltmle_exact_rbind_fill(
      separate_clever_covariate_identity_rows
    ),
    second_M2_marginal_reference_check = .ltmle_exact_rbind_fill(
      second_M2_marginal_reference_rows
    ),
    mc_integration_diagnostics = mc_diag,
    law_task_integration_diagnostics = if (length(law_task_integration_rows)) {
      .ltmle_exact_rbind_fill(law_task_integration_rows)
    } else {
      data.frame()
    },
    branch_history_state_check = if (length(branch_trace_rows)) {
      .ltmle_exact_rbind_fill(branch_trace_rows)
    } else {
      data.frame()
    },
    branch_state_handoff_trace = if (length(handoff_trace_rows)) {
      .ltmle_exact_rbind_fill(handoff_trace_rows)
    } else {
      data.frame()
    },
    law_integration_diagnostics = if (length(branch_law_integration_rows)) {
      .ltmle_exact_rbind_fill(branch_law_integration_rows)
    } else {
      data.frame()
    },
    source_eval_training_path_check = if (length(source_eval_training_rows)) {
      .ltmle_exact_rbind_fill(source_eval_training_rows)
    } else {
      data.frame()
    },
    cross_regimen_source_boundary_trace = if (length(cross_regimen_source_boundary_rows)) {
      .ltmle_exact_rbind_fill(cross_regimen_source_boundary_rows)
    } else {
      data.frame()
    },
    source_pseudooutcome_independent_recursion_check = if (length(source_pseudooutcome_independent_rows)) {
      .ltmle_exact_rbind_fill(source_pseudooutcome_independent_rows)
    } else {
      data.frame()
    },
    targeted_continuation_storage_check = if (length(targeted_storage_rows)) {
      .ltmle_exact_rbind_fill(targeted_storage_rows)
    } else {
      data.frame()
    },
    dedicated_L_transition_L_fit_training_row_trace = if (length(dedicated_L_training_row_trace_rows)) {
      .ltmle_exact_rbind_fill(dedicated_L_training_row_trace_rows)
    } else {
      data.frame()
    },
    dedicated_L_transition_L_fit_prediction_row_trace = if (length(dedicated_L_prediction_row_trace_rows)) {
      .ltmle_exact_rbind_fill(dedicated_L_prediction_row_trace_rows)
    } else {
      data.frame()
    },
    dedicated_L_transition_L_fit_coefficient_contribution_check =
      if (length(dedicated_L_coefficient_contribution_rows)) {
        .ltmle_exact_rbind_fill(dedicated_L_coefficient_contribution_rows)
      } else {
        data.frame()
      },
    dedicated_L_transition_prediction_weighting_collapse_check =
      if (length(dedicated_L_prediction_weighting_collapse_rows)) {
        .ltmle_exact_rbind_fill(dedicated_L_prediction_weighting_collapse_rows)
      } else {
        data.frame()
      },
    dedicated_L_transition_mediator_path_weight_trace =
      if (length(dedicated_L_mediator_path_weight_trace_rows)) {
        .ltmle_exact_rbind_fill(dedicated_L_mediator_path_weight_trace_rows)
      } else {
        data.frame()
      },
    dedicated_L_transition_training_response_oracle_alignment_check =
      if (length(dedicated_L_training_response_oracle_rows)) {
        .ltmle_exact_rbind_fill(dedicated_L_training_response_oracle_rows)
      } else {
        data.frame()
      },
    observed_vs_source_row_role_check = if (length(observed_vs_source_rows)) {
      .ltmle_exact_rbind_fill(observed_vs_source_rows)
    } else {
      data.frame()
    },
    root_plugin_diagnostics = if (length(root_plugin_rows)) {
      .ltmle_exact_rbind_fill(root_plugin_rows)
    } else {
      data.frame()
    },
    root_vs_targeted_branch_recursion_check = if (length(root_vs_branch_rows)) {
      .ltmle_exact_rbind_fill(root_vs_branch_rows)
    } else {
      data.frame()
    },
    task_graph = task_graph_slim
  )
}

.ltmle_exact_component_summary <- function(eif_matrix, tolerance) {
  n <- nrow(eif_matrix)
  mean_D <- as.numeric(colMeans(eif_matrix, na.rm = TRUE))
  sd_D <- apply(eif_matrix, 2L, stats::sd, na.rm = TRUE)
  se_D <- sd_D / sqrt(n)
  scaled_Z <- abs(mean_D) / se_D
  data.frame(
    component = colnames(eif_matrix),
    mean_D = mean_D,
    sd_D = sd_D,
    se_D = se_D,
    scaled_Z = scaled_Z,
    component_tolerance = tolerance,
    component_equation_solved = is.finite(scaled_Z) & abs(mean_D) <= tolerance,
    stringsAsFactors = FALSE
  )
}

.ltmle_exact_normalize_density_ratio_diagnostics <- function(density_ratio_diagnostics,
                                                             requested_density_ratio_mc_n,
                                                             requested_law_integration_n) {
  if (is.null(density_ratio_diagnostics) || !nrow(density_ratio_diagnostics)) {
    return(data.frame())
  }
  n <- nrow(density_ratio_diagnostics)
  add_if_missing <- function(dat, nm, value) {
    if (!nm %in% names(dat)) {
      dat[[nm]] <- rep(value, n)
    } else {
      dat[[nm]][is.na(dat[[nm]])] <- value
    }
    dat
  }
  density_ratio_diagnostics <- add_if_missing(
    density_ratio_diagnostics,
    "requested_density_ratio_mc_n",
    as.integer(requested_density_ratio_mc_n)
  )
  density_ratio_diagnostics <- add_if_missing(
    density_ratio_diagnostics,
    "effective_density_ratio_mc_n",
    as.integer(requested_density_ratio_mc_n)
  )
  density_ratio_diagnostics <- add_if_missing(
    density_ratio_diagnostics,
    "density_ratio_mc_cap_applied",
    FALSE
  )
  density_ratio_diagnostics <- add_if_missing(
    density_ratio_diagnostics,
    "requested_law_integration_n",
    as.integer(requested_law_integration_n)
  )
  density_ratio_diagnostics <- add_if_missing(
    density_ratio_diagnostics,
    "effective_law_integration_n",
    as.integer(requested_law_integration_n)
  )
  density_ratio_diagnostics <- add_if_missing(
    density_ratio_diagnostics,
    "law_integration_cap_applied",
    FALSE
  )
  density_ratio_diagnostics <- add_if_missing(
    density_ratio_diagnostics,
    "particle_pruning_applied",
    FALSE
  )
  density_ratio_diagnostics <- add_if_missing(
    density_ratio_diagnostics,
    "max_generated_particles",
    0L
  )
  density_ratio_diagnostics <- add_if_missing(
    density_ratio_diagnostics,
    "run_mode",
    NA_character_
  )
  density_ratio_diagnostics <- add_if_missing(
    density_ratio_diagnostics,
    "is_acceptance_gate",
    FALSE
  )
  required <- c(
    "requested_density_ratio_mc_n",
    "effective_density_ratio_mc_n",
    "density_ratio_mc_cap_applied",
    "requested_law_integration_n",
    "effective_law_integration_n",
    "law_integration_cap_applied",
    "particle_pruning_applied",
    "max_generated_particles"
  )
  required_missing <- Reduce(`|`, lapply(required, function(nm) is.na(density_ratio_diagnostics[[nm]])))
  finite_ratios <- is.finite(as.numeric(density_ratio_diagnostics$max_ratio %||% rep(1, n)))
  density_ratio_diagnostics$passed <- !required_missing &
    !as.logical(density_ratio_diagnostics$density_ratio_mc_cap_applied) &
    !as.logical(density_ratio_diagnostics$law_integration_cap_applied) &
    !as.logical(density_ratio_diagnostics$particle_pruning_applied) &
    finite_ratios
  density_ratio_diagnostics$failure_class <- ifelse(
    density_ratio_diagnostics$passed,
    "no_failure",
    "density_ratio_diagnostic_failed"
  )
  density_ratio_diagnostics
}

.ltmle_exact_separate_product_join_diagnostics <- function(branch_state_handoff_trace,
                                                           scenario_id = NA_character_,
                                                           run_mode = NA_character_,
                                                           detail = c("row", "event")) {
  detail <- match.arg(detail)
  empty_trace <- data.frame(
    scenario_id = character(0),
    run_mode = character(0),
    component = character(0),
    task_id = character(0),
    t = integer(0),
    join_type = character(0),
    id0 = integer(0),
    n_first_particles = integer(0),
    n_second_particles = integer(0),
    n_joined_particles = integer(0),
    expected_product_count = integer(0),
    product_count_matches_expected = logical(0),
    row_order_copy_used = logical(0),
    first_weight = numeric(0),
    second_weight = numeric(0),
    joined_weight = numeric(0),
    expected_joined_weight = numeric(0),
    weight_product_correct = logical(0),
    sum_joined_weight_by_id0 = numeric(0),
    sum_weight_close_to_one_by_id0 = logical(0),
    first_auxiliary_mediator_copied_to_outcome = logical(0),
    second_auxiliary_mediator_copied_to_outcome = logical(0),
    target_M1_copied_to_outcome = logical(0),
    target_M2_copied_to_outcome = logical(0),
    event_id = integer(0),
    first_law_id = character(0),
    second_law_id = character(0),
    first_branch_id = integer(0),
    second_branch_id = integer(0),
    actual_product_count = integer(0),
    joined_weight_raw = numeric(0),
    joined_weight_normalized = numeric(0),
    sum_joined_weight_raw_by_event = numeric(0),
    sum_joined_weight_normalized_by_event = numeric(0),
    outer_product_weights_used = logical(0),
    target_m1_source = character(0),
    target_m2_source = character(0),
    auxiliary_mediator_copied_from_wrong_branch = logical(0),
    event_weight_normalization_passed = logical(0),
    product_join_nondegeneracy_gate_applicable = logical(0),
    product_join_nondegenerate_production_rows = integer(0),
    product_join_nondegeneracy_gate_status = character(0),
    synthetic_helper_used_as_production_coverage = logical(0),
    passed = logical(0),
    failure_class = character(0),
    stringsAsFactors = FALSE
  )
  empty_pre <- data.frame(
    scenario_id = character(0),
    run_mode = character(0),
    event_id = integer(0),
    component = character(0),
    task_id = character(0),
    t = integer(0),
    id0 = integer(0),
    n_first_particles_pre_collapse = integer(0),
    n_second_particles_pre_collapse = integer(0),
    n_outcome_particles_pre_collapse = integer(0),
    expected_product_count = integer(0),
    product_count_matches_expected = logical(0),
    row_order_copy_used = logical(0),
    first_auxiliary_M2_copied_to_outcome = logical(0),
    second_auxiliary_M1_copied_to_outcome = logical(0),
    target_M1_copied_to_outcome = logical(0),
    target_M2_copied_to_outcome = logical(0),
    sum_weight_before_collapse = numeric(0),
    sum_weight_after_collapse = numeric(0),
    max_abs_sum_weight_after_collapse_minus_1 = numeric(0),
    product_join_nondegeneracy_gate_applicable = logical(0),
    product_join_nondegenerate_production_rows = integer(0),
    product_join_nondegeneracy_gate_status = character(0),
    tolerance = numeric(0),
    passed = logical(0),
    failure_class = character(0),
    stringsAsFactors = FALSE
  )
  if (is.null(branch_state_handoff_trace) || !nrow(branch_state_handoff_trace) ||
      !"handoff_type" %in% names(branch_state_handoff_trace)) {
    return(list(trace = empty_trace, pre_collapse_trace = empty_pre))
  }

  product_rows <- branch_state_handoff_trace[
    branch_state_handoff_trace$handoff_type == "separate_world_product_join",
    ,
    drop = FALSE
  ]
  product_rows <- product_rows[
    !is.na(product_rows$expected_n_outcome_particles) &
      !is.na(product_rows$n_outcome_particles_after_join),
    ,
    drop = FALSE
  ]
  if (!nrow(product_rows)) {
    return(list(trace = empty_trace, pre_collapse_trace = empty_pre))
  }

  val <- function(dat, nm, default) if (nm %in% names(dat)) dat[[nm]] else rep(default, nrow(dat))
  joined_weight <- as.numeric(val(product_rows, "joined_weight", val(product_rows, "weight_product", NA_real_)))
  expected_joined_weight <- as.numeric(val(
    product_rows,
    "expected_joined_weight",
    as.numeric(val(product_rows, "weight_first", NA_real_)) *
      as.numeric(val(product_rows, "weight_second", NA_real_))
  ))
  product_event_start <- c(
    TRUE,
    product_rows$component[-1L] != product_rows$component[-nrow(product_rows)] |
      product_rows$task_id[-1L] != product_rows$task_id[-nrow(product_rows)] |
      product_rows$id0[-1L] != product_rows$id0[-nrow(product_rows)]
  )
  product_event_id <- cumsum(product_event_start)
  product_group <- interaction(
    product_rows$component,
    product_rows$task_id,
    product_rows$id0,
    product_event_id,
    drop = TRUE
  )
  sum_by_id0 <- ave(joined_weight, product_group, FUN = function(x) sum(x, na.rm = TRUE))
  product_join_nondegenerate_production_rows <-
    sum(product_rows$expected_n_outcome_particles > 1L, na.rm = TRUE)
  product_join_nondegeneracy_gate_applicable <- TRUE
  product_join_nondegeneracy_gate_status <- if (product_join_nondegenerate_production_rows > 0L) {
    "passed"
  } else {
    "failed"
  }

  if (identical(detail, "event")) {
    event_starts <- which(product_event_start)
    event_ends <- c(event_starts[-1L] - 1L, nrow(product_rows))
    event_seq <- seq_along(event_starts)
    event_sum <- function(x) {
      vapply(event_seq, function(ii) sum(x[event_starts[ii]:event_ends[ii]], na.rm = TRUE), numeric(1))
    }
    event_all <- function(x) {
      vapply(event_seq, function(ii) all(x[event_starts[ii]:event_ends[ii]] %in% TRUE), logical(1))
    }
    event_any <- function(x) {
      vapply(event_seq, function(ii) any(x[event_starts[ii]:event_ends[ii]] %in% TRUE), logical(1))
    }
    first_idx <- event_starts
    joined_weight_raw <- as.numeric(val(product_rows, "joined_weight_raw", joined_weight))
    joined_weight_normalized <- as.numeric(val(product_rows, "joined_weight_normalized", joined_weight))
    joined_sum <- event_sum(joined_weight)
    joined_raw_sum <- event_sum(joined_weight_raw)
    joined_normalized_sum <- event_sum(joined_weight_normalized)
    trace <- data.frame(
      scenario_id = scenario_id,
      run_mode = run_mode,
      event_id = product_event_id[first_idx],
      component = product_rows$component[first_idx],
      task_id = product_rows$task_id[first_idx],
      first_law_id = paste0(product_rows$component[first_idx], "::first_law"),
      second_law_id = paste0(product_rows$component[first_idx], "::second_law"),
      first_branch_id = as.integer(val(product_rows, "first_particle_id", NA_integer_))[first_idx],
      second_branch_id = as.integer(val(product_rows, "second_particle_id", NA_integer_))[first_idx],
      t = product_rows$t[first_idx],
      join_type = "separate_world_product_join",
      id0 = product_rows$id0[first_idx],
      n_first_particles = product_rows$n_first_particles[first_idx],
      n_second_particles = product_rows$n_second_particles[first_idx],
      n_joined_particles = product_rows$n_outcome_particles_after_join[first_idx],
      expected_product_count = product_rows$expected_n_outcome_particles[first_idx],
      actual_product_count = product_rows$n_outcome_particles_after_join[first_idx],
      product_count_matches_expected = product_rows$n_outcome_particles_after_join[first_idx] ==
        product_rows$expected_n_outcome_particles[first_idx],
      row_order_copy_used = event_any(as.logical(product_rows$row_order_copy_used)),
      first_weight = as.numeric(product_rows$weight_first)[first_idx],
      second_weight = as.numeric(product_rows$weight_second)[first_idx],
      joined_weight = joined_weight[first_idx],
      joined_weight_raw = joined_weight_raw[first_idx],
      joined_weight_normalized = joined_weight_normalized[first_idx],
      expected_joined_weight = expected_joined_weight[first_idx],
      weight_product_correct = event_all(abs(joined_weight - expected_joined_weight) <= 1e-10),
      sum_joined_weight_by_id0 = joined_sum,
      sum_joined_weight_raw_by_event = joined_raw_sum,
      sum_joined_weight_normalized_by_event = joined_normalized_sum,
      sum_weight_close_to_one_by_id0 = abs(joined_sum - 1) <= 1e-8,
      outer_product_weights_used = TRUE,
      first_auxiliary_mediator_copied_to_outcome =
        event_any(as.logical(product_rows$first_auxiliary_M2_copied_to_outcome)),
      second_auxiliary_mediator_copied_to_outcome =
        event_any(as.logical(product_rows$second_auxiliary_M1_copied_to_outcome)),
      target_M1_copied_to_outcome = event_all(as.logical(val(product_rows, "target_M1_copied_to_outcome", TRUE))),
      target_M2_copied_to_outcome = event_all(as.logical(val(product_rows, "target_M2_copied_to_outcome", TRUE))),
      target_m1_source = "first_law_target_path",
      target_m2_source = "second_law_target_path",
      auxiliary_mediator_copied_from_wrong_branch =
        event_any(as.logical(product_rows$first_auxiliary_M2_copied_to_outcome) |
                    as.logical(product_rows$second_auxiliary_M1_copied_to_outcome)),
      event_weight_normalization_passed = abs(joined_normalized_sum - 1) <= 1e-8,
      product_join_nondegeneracy_gate_applicable = product_join_nondegeneracy_gate_applicable,
      product_join_nondegenerate_production_rows = product_join_nondegenerate_production_rows,
      product_join_nondegeneracy_gate_status = product_join_nondegeneracy_gate_status,
      synthetic_helper_used_as_production_coverage = FALSE,
      stringsAsFactors = FALSE
    )
    trace$passed <- as.logical(trace$product_count_matches_expected) &
      !as.logical(trace$row_order_copy_used) &
      as.logical(trace$weight_product_correct) &
      as.logical(trace$event_weight_normalization_passed) &
      !as.logical(trace$first_auxiliary_mediator_copied_to_outcome) &
      !as.logical(trace$second_auxiliary_mediator_copied_to_outcome) &
      as.logical(trace$target_M1_copied_to_outcome) &
      as.logical(trace$target_M2_copied_to_outcome) &
      trace$product_join_nondegenerate_production_rows > 0L
    trace$failure_class <- ifelse(
      trace$passed,
      "no_failure",
      ifelse(trace$product_join_nondegenerate_production_rows <= 0L,
             "separate_product_join_degenerate_trace_only",
             "separate_product_join_trace_failed")
    )

    pre <- data.frame(
      scenario_id = scenario_id,
      run_mode = run_mode,
      event_id = product_event_id[first_idx],
      component = product_rows$component[first_idx],
      task_id = product_rows$task_id[first_idx],
      t = product_rows$t[first_idx],
      id0 = product_rows$id0[first_idx],
      n_first_particles_pre_collapse = product_rows$n_first_particles[first_idx],
      n_second_particles_pre_collapse = product_rows$n_second_particles[first_idx],
      n_outcome_particles_pre_collapse = product_rows$n_outcome_particles_after_join[first_idx],
      expected_product_count = product_rows$expected_n_outcome_particles[first_idx],
      product_count_matches_expected = product_rows$n_outcome_particles_after_join[first_idx] ==
        product_rows$expected_n_outcome_particles[first_idx],
      row_order_copy_used = event_any(as.logical(product_rows$row_order_copy_used)),
      first_auxiliary_M2_copied_to_outcome =
        event_any(as.logical(product_rows$first_auxiliary_M2_copied_to_outcome)),
      second_auxiliary_M1_copied_to_outcome =
        event_any(as.logical(product_rows$second_auxiliary_M1_copied_to_outcome)),
      target_M1_copied_to_outcome = event_all(as.logical(val(product_rows, "target_M1_copied_to_outcome", TRUE))),
      target_M2_copied_to_outcome = event_all(as.logical(val(product_rows, "target_M2_copied_to_outcome", TRUE))),
      sum_weight_before_collapse = event_sum(expected_joined_weight),
      sum_weight_after_collapse = event_sum(joined_weight),
      product_join_nondegeneracy_gate_applicable = product_join_nondegeneracy_gate_applicable,
      product_join_nondegenerate_production_rows = product_join_nondegenerate_production_rows,
      product_join_nondegeneracy_gate_status = product_join_nondegeneracy_gate_status,
      tolerance = 1e-8,
      stringsAsFactors = FALSE
    )
    pre$max_abs_sum_weight_after_collapse_minus_1 <- abs(as.numeric(pre$sum_weight_after_collapse) - 1)
    pre$passed <- as.logical(pre$product_count_matches_expected) &
      !as.logical(pre$row_order_copy_used) &
      !as.logical(pre$first_auxiliary_M2_copied_to_outcome) &
      !as.logical(pre$second_auxiliary_M1_copied_to_outcome) &
      as.logical(pre$target_M1_copied_to_outcome) &
      as.logical(pre$target_M2_copied_to_outcome) &
      pre$max_abs_sum_weight_after_collapse_minus_1 <= pre$tolerance &
      pre$product_join_nondegenerate_production_rows > 0L
    pre$failure_class <- ifelse(
      pre$passed,
      "no_failure",
      ifelse(pre$product_join_nondegenerate_production_rows <= 0L,
             "separate_product_join_degenerate_trace_only",
             "separate_product_join_pre_collapse_failed")
    )

    return(list(trace = trace, pre_collapse_trace = pre))
  }

  trace <- data.frame(
    scenario_id = scenario_id,
    run_mode = run_mode,
    event_id = product_event_id,
    component = product_rows$component,
    task_id = product_rows$task_id,
    first_law_id = paste0(product_rows$component, "::first_law"),
    second_law_id = paste0(product_rows$component, "::second_law"),
    first_branch_id = as.integer(val(product_rows, "first_particle_id", NA_integer_)),
    second_branch_id = as.integer(val(product_rows, "second_particle_id", NA_integer_)),
    t = product_rows$t,
    join_type = "separate_world_product_join",
    id0 = product_rows$id0,
    n_first_particles = product_rows$n_first_particles,
    n_second_particles = product_rows$n_second_particles,
    n_joined_particles = product_rows$n_outcome_particles_after_join,
    expected_product_count = product_rows$expected_n_outcome_particles,
    actual_product_count = product_rows$n_outcome_particles_after_join,
    product_count_matches_expected = product_rows$n_outcome_particles_after_join ==
      product_rows$expected_n_outcome_particles,
    row_order_copy_used = as.logical(product_rows$row_order_copy_used),
    first_weight = as.numeric(product_rows$weight_first),
    second_weight = as.numeric(product_rows$weight_second),
    joined_weight = joined_weight,
    joined_weight_raw = as.numeric(val(product_rows, "joined_weight_raw", joined_weight)),
    joined_weight_normalized = as.numeric(val(product_rows, "joined_weight_normalized", joined_weight)),
    expected_joined_weight = expected_joined_weight,
    weight_product_correct = abs(joined_weight - expected_joined_weight) <= 1e-10,
    sum_joined_weight_by_id0 = as.numeric(sum_by_id0),
    sum_joined_weight_raw_by_event = ave(
      as.numeric(val(product_rows, "joined_weight_raw", joined_weight)),
      product_group,
      FUN = function(x) sum(x, na.rm = TRUE)
    ),
    sum_joined_weight_normalized_by_event = ave(
      as.numeric(val(product_rows, "joined_weight_normalized", joined_weight)),
      product_group,
      FUN = function(x) sum(x, na.rm = TRUE)
    ),
    sum_weight_close_to_one_by_id0 = abs(as.numeric(sum_by_id0) - 1) <= 1e-8,
    outer_product_weights_used = TRUE,
    first_auxiliary_mediator_copied_to_outcome =
      as.logical(product_rows$first_auxiliary_M2_copied_to_outcome),
    second_auxiliary_mediator_copied_to_outcome =
      as.logical(product_rows$second_auxiliary_M1_copied_to_outcome),
    target_M1_copied_to_outcome = as.logical(val(product_rows, "target_M1_copied_to_outcome", TRUE)),
    target_M2_copied_to_outcome = as.logical(val(product_rows, "target_M2_copied_to_outcome", TRUE)),
    target_m1_source = "first_law_target_path",
    target_m2_source = "second_law_target_path",
    auxiliary_mediator_copied_from_wrong_branch =
      as.logical(product_rows$first_auxiliary_M2_copied_to_outcome) |
        as.logical(product_rows$second_auxiliary_M1_copied_to_outcome),
    event_weight_normalization_passed = abs(as.numeric(sum_by_id0) - 1) <= 1e-8,
    product_join_nondegeneracy_gate_applicable = product_join_nondegeneracy_gate_applicable,
    product_join_nondegenerate_production_rows = product_join_nondegenerate_production_rows,
    product_join_nondegeneracy_gate_status = product_join_nondegeneracy_gate_status,
    synthetic_helper_used_as_production_coverage = FALSE,
    stringsAsFactors = FALSE
  )
  trace$passed <- as.logical(trace$product_count_matches_expected) &
    !as.logical(trace$row_order_copy_used) &
    as.logical(trace$weight_product_correct) &
    as.logical(trace$event_weight_normalization_passed) &
    !as.logical(trace$first_auxiliary_mediator_copied_to_outcome) &
    !as.logical(trace$second_auxiliary_mediator_copied_to_outcome) &
    as.logical(trace$target_M1_copied_to_outcome) &
    as.logical(trace$target_M2_copied_to_outcome) &
    trace$product_join_nondegenerate_production_rows > 0L
  trace$failure_class <- ifelse(
    trace$passed,
    "no_failure",
    ifelse(trace$product_join_nondegenerate_production_rows <= 0L,
           "separate_product_join_degenerate_trace_only",
           "separate_product_join_trace_failed")
  )

  pre <- data.frame(
    scenario_id = scenario_id,
    run_mode = run_mode,
    event_id = product_event_id,
    component = product_rows$component,
    task_id = product_rows$task_id,
    t = product_rows$t,
    id0 = product_rows$id0,
    n_first_particles_pre_collapse = product_rows$n_first_particles,
    n_second_particles_pre_collapse = product_rows$n_second_particles,
    n_outcome_particles_pre_collapse = product_rows$n_outcome_particles_after_join,
    expected_product_count = product_rows$expected_n_outcome_particles,
    product_count_matches_expected = product_rows$n_outcome_particles_after_join ==
      product_rows$expected_n_outcome_particles,
    row_order_copy_used = as.logical(product_rows$row_order_copy_used),
    first_auxiliary_M2_copied_to_outcome =
      as.logical(product_rows$first_auxiliary_M2_copied_to_outcome),
    second_auxiliary_M1_copied_to_outcome =
      as.logical(product_rows$second_auxiliary_M1_copied_to_outcome),
    target_M1_copied_to_outcome = as.logical(val(product_rows, "target_M1_copied_to_outcome", TRUE)),
    target_M2_copied_to_outcome = as.logical(val(product_rows, "target_M2_copied_to_outcome", TRUE)),
    sum_weight_before_collapse = ave(
      as.numeric(val(product_rows, "expected_joined_weight", product_rows$joined_weight)),
      product_group,
      FUN = function(x) sum(x, na.rm = TRUE)
    ),
    sum_weight_after_collapse = ave(
      as.numeric(val(product_rows, "joined_weight", product_rows$weight_product)),
      product_group,
      FUN = function(x) sum(x, na.rm = TRUE)
    ),
    product_join_nondegeneracy_gate_applicable = product_join_nondegeneracy_gate_applicable,
    product_join_nondegenerate_production_rows = product_join_nondegenerate_production_rows,
    product_join_nondegeneracy_gate_status = product_join_nondegeneracy_gate_status,
    tolerance = 1e-8,
    stringsAsFactors = FALSE
  )
  pre$max_abs_sum_weight_after_collapse_minus_1 <- abs(as.numeric(pre$sum_weight_after_collapse) - 1)
  pre$passed <- as.logical(pre$product_count_matches_expected) &
    !as.logical(pre$row_order_copy_used) &
    !as.logical(pre$first_auxiliary_M2_copied_to_outcome) &
    !as.logical(pre$second_auxiliary_M1_copied_to_outcome) &
    as.logical(pre$target_M1_copied_to_outcome) &
    as.logical(pre$target_M2_copied_to_outcome) &
    pre$max_abs_sum_weight_after_collapse_minus_1 <= pre$tolerance &
    pre$product_join_nondegenerate_production_rows > 0L
  pre$failure_class <- ifelse(
    pre$passed,
    "no_failure",
    ifelse(pre$product_join_nondegenerate_production_rows <= 0L,
           "separate_product_join_degenerate_trace_only",
           "separate_product_join_pre_collapse_failed")
  )

  list(trace = trace, pre_collapse_trace = pre)
}

.ltmle_exact_role_key <- function(role, n) {
  if (is.null(role)) role <- rep(NA_character_, n)
  role <- as.character(role)
  role[is.na(role)] <- ""
  role
}

.ltmle_exact_required_checks <- function(score_diagnostics,
                                         component_summary,
                                         factor_tasks,
                                         scaled_z_tolerance) {
  key_cols <- c("component", "t", "node", "process_type", ".role_key")
  factor_tasks$.role_key <- .ltmle_exact_role_key(factor_tasks$role, nrow(factor_tasks))
  score_diagnostics$.role_key <- .ltmle_exact_role_key(
    if ("role" %in% names(score_diagnostics)) score_diagnostics$role else NULL,
    nrow(score_diagnostics)
  )

  task_key <- unique(factor_tasks[, key_cols, drop = FALSE])
  score_key <- unique(score_diagnostics[, key_cols, drop = FALSE])
  task_id <- do.call(paste, task_key[, key_cols, drop = FALSE])
  score_id <- do.call(paste, score_key[, key_cols, drop = FALSE])

  missing_tasks <- task_key[!task_id %in% score_id, , drop = FALSE]
  extra_scores <- score_key[!score_id %in% task_id, , drop = FALSE]

  bad_score <- !isTRUE(all(score_diagnostics$score_equation_solved))
  bad_eif <- !isTRUE(all(component_summary$component_equation_solved))
  bad_z <- !isTRUE(all(is.finite(component_summary$scaled_Z))) ||
    max(component_summary$scaled_Z, na.rm = TRUE) > scaled_z_tolerance
  base_prediction_cols <- intersect(
    c("epsilon", "score_before", "score_after", "H_mean", "H_max"),
    names(score_diagnostics)
  )

  prediction_checks_passed <- length(base_prediction_cols) > 0L &&
    all(vapply(
      score_diagnostics[, base_prediction_cols, drop = FALSE],
      function(x) all(is.finite(x)),
      logical(1)
    ))

  generated_evaluated <- if ("generated_evaluated" %in% names(score_diagnostics)) {
    as.logical(score_diagnostics$generated_evaluated)
  } else {
    rep(TRUE, nrow(score_diagnostics))
  }
  generated_evaluated[is.na(generated_evaluated)] <- FALSE

  if (any(generated_evaluated)) {
    generated_prediction_cols <- intersect(
      c("H_new_mean", "H_new_max"),
      names(score_diagnostics)
    )
    prediction_checks_passed <- prediction_checks_passed &&
      length(generated_prediction_cols) == 2L &&
      all(vapply(
        score_diagnostics[generated_evaluated, generated_prediction_cols, drop = FALSE],
        function(x) all(is.finite(x)),
        logical(1)
      ))
  }

  density_checks_passed <- all(is.finite(score_diagnostics$H_max)) &&
    all(score_diagnostics$H_max >= 0)

  if (any(generated_evaluated)) {
    density_checks_passed <- density_checks_passed &&
      all(is.finite(score_diagnostics$H_new_max[generated_evaluated])) &&
      all(score_diagnostics$H_new_max[generated_evaluated] >= 0)
  }
  list(
    missing_tasks = missing_tasks,
    extra_scores = extra_scores,
    bad_score = bad_score,
    bad_eif = bad_eif,
    bad_z = bad_z,
    all_required_factor_tasks_present = nrow(missing_tasks) == 0L && nrow(extra_scores) == 0L,
    all_predictions_strictly_valid = prediction_checks_passed,
    all_density_ratios_strictly_valid = density_checks_passed,
    ok = nrow(missing_tasks) == 0L && nrow(extra_scores) == 0L &&
      !bad_score && !bad_eif && !bad_z
  )
}

.ltmle_exact_missing_law_specific_equations <- function(score_diagnostics, registry, T,
                                                        node_spec = NULL) {
  L_nodes <- .ltmle_exact_L_nodes(node_spec)
  score_role <- .ltmle_exact_role_key(
    if ("role" %in% names(score_diagnostics)) score_diagnostics$role else NULL,
    nrow(score_diagnostics)
  )
  score_key <- paste(
    score_diagnostics$component,
    as.integer(score_diagnostics$t),
    score_diagnostics$node,
    score_diagnostics$process_type,
    score_role,
    sep = "\r"
  )

  rows <- list()
  add_expected <- function(component, world_type, tt, node, process_type, role = NA_character_) {
    rows[[length(rows) + 1L]] <<- data.frame(
      component = component,
      world_type = world_type,
      t = as.integer(tt),
      node = node,
      process_type = process_type,
      role = role,
      stringsAsFactors = FALSE
    )
  }

  first_M2_role <- "auxiliary second mediator in the first-mediator stochastic intervention law"
  law_Y_role <- "outcome history for stochastic mediator intervention law"
  outcome_type <- .ltmle_exact_outcome_type(node_spec)
  for (ii in seq_len(nrow(registry))) {
    comp <- registry$component[ii]
    wt <- registry$world_type[ii]
    if (identical(wt, "joint")) {
      for (tt in seq_len(max(T - 1L, 0L))) {
        for (L_node in L_nodes) {
          add_expected(comp, wt, tt, L_node, "joint_stochastic_mediator_intervention_law")
        }
        if (identical(outcome_type, "longitudinal")) {
          add_expected(
            comp, wt, tt, "Y", "joint_stochastic_mediator_intervention_law",
            role = law_Y_role
          )
        }
      }
    }
    if (identical(wt, "separate")) {
      for (tt in seq_len(max(T - 1L, 0L))) {
        for (L_node in L_nodes) {
          add_expected(comp, wt, tt, L_node, "first_mediator_stochastic_intervention_law")
          add_expected(comp, wt, tt, L_node, "second_mediator_stochastic_intervention_law")
        }
        if (identical(outcome_type, "longitudinal")) {
          add_expected(
            comp, wt, tt, "Y", "first_mediator_stochastic_intervention_law",
            role = law_Y_role
          )
          add_expected(
            comp, wt, tt, "Y", "second_mediator_stochastic_intervention_law",
            role = law_Y_role
          )
        }
      }
      if (T >= 3L) {
        for (tt in seq.int(2L, T - 1L)) {
          add_expected(
            comp, wt, tt, "M2", "first_mediator_stochastic_intervention_law",
            role = first_M2_role
          )
        }
      }
    }
  }

  if (!length(rows)) {
    return(data.frame(
      component = character(0), world_type = character(0), t = integer(0),
      node = character(0), process_type = character(0), role = character(0),
      stringsAsFactors = FALSE
    ))
  }
  expected <- do.call(rbind, rows)
  expected_role <- .ltmle_exact_role_key(expected$role, nrow(expected))
  expected_key <- paste(
    expected$component,
    expected$t,
    expected$node,
    expected$process_type,
    expected_role,
    sep = "\r"
  )
  expected[!expected_key %in% score_key, , drop = FALSE]
}

.ltmle_exact_validate_stochastic_law_registry <- function(factor_tasks, T, node_spec = NULL) {
  law_process_types <- c(
    "joint_stochastic_mediator_intervention_law",
    "first_mediator_stochastic_intervention_law",
    "second_mediator_stochastic_intervention_law"
  )
  law_tasks <- factor_tasks[factor_tasks$process_type %in% law_process_types, , drop = FALSE]
  if (!nrow(law_tasks)) return(invisible(TRUE))

  baseline_mediator_tasks <- law_tasks[
    law_tasks$t == 1L & law_tasks$node %in% c("M1", "M2"),
    ,
    drop = FALSE
  ]
  if (nrow(baseline_mediator_tasks)) {
    .stop(
      "ltmle_exact must not register baseline mediator stochastic intervention law tasks."
    )
  }

  first_M2_role <- "auxiliary second mediator in the first-mediator stochastic intervention law"
  final_first_law_auxiliary_M2 <- law_tasks[
    law_tasks$t == T &
      law_tasks$node == "M2" &
      law_tasks$process_type == "first_mediator_stochastic_intervention_law" &
      law_tasks$role %in% first_M2_role,
    ,
    drop = FALSE
  ]
  if (nrow(final_first_law_auxiliary_M2)) {
    .stop(
      "ltmle_exact must not register a final-visit auxiliary second-mediator task ",
      "inside the first-mediator stochastic intervention law."
    )
  }

  final_first_law_M1 <- law_tasks[
    law_tasks$t == T &
      law_tasks$node == "M1" &
      law_tasks$process_type == "first_mediator_stochastic_intervention_law",
    ,
    drop = FALSE
  ]
  if (nrow(final_first_law_M1)) {
    bad <- grepl(
      "first_mediator_stochastic_intervention_law",
      final_first_law_M1$observed_pseudooutcome_source_task_id,
      fixed = TRUE
    ) | grepl(
      "auxiliary second mediator",
      final_first_law_M1$observed_pseudooutcome_source_task_id,
      fixed = TRUE
    )
    if (any(bad)) {
      .stop(
        "The final first-mediator stochastic intervention law M1 task must point ",
        "directly to the outcome-process post-mediator covariate continuation, not ",
        "to an auxiliary second-mediator task."
      )
    }
  }

  if (any(law_tasks$observed_pseudooutcome_source_task_id == "observed_terminal_outcome")) {
    .stop("ltmle_exact stochastic mediator intervention law task points to observed_terminal_outcome.")
  }

  final_law_L <- law_tasks[
    law_tasks$t == T & law_tasks$node %in% .ltmle_exact_L_nodes(node_spec),
    ,
    drop = FALSE
  ]
  if (nrow(final_law_L)) {
    .stop("ltmle_exact registered final-visit post-mediator covariate tasks for a stochastic mediator intervention law.")
  }

  if (identical(.ltmle_exact_outcome_type(node_spec), "longitudinal") && T > 1L) {
    rows <- unique(factor_tasks[, c("component", "world_type"), drop = FALSE])
    expected <- list()
    for (ii in seq_len(nrow(rows))) {
      process_types <- if (identical(rows$world_type[ii], "joint")) {
        "joint_stochastic_mediator_intervention_law"
      } else if (identical(rows$world_type[ii], "separate")) {
        c("first_mediator_stochastic_intervention_law",
          "second_mediator_stochastic_intervention_law")
      } else {
        character(0)
      }
      for (process_type in process_types) {
        for (tt in seq_len(T - 1L)) {
          expected[[length(expected) + 1L]] <- paste(rows$component[ii], tt, process_type, sep = "\r")
        }
      }
    }
    if (length(expected)) {
      have <- paste(
        law_tasks$component[law_tasks$t < T & law_tasks$node == "Y"],
        law_tasks$t[law_tasks$t < T & law_tasks$node == "Y"],
        law_tasks$process_type[law_tasks$t < T & law_tasks$node == "Y"],
        sep = "\r"
      )
      missing <- setdiff(unlist(expected, use.names = FALSE), have)
      if (length(missing)) {
        .stop("ltmle_exact is missing law-specific outcome-history tasks for longitudinal mediator histories.")
      }
    }
  }

  invisible(TRUE)
}

.ltmle_exact_fail_closed <- function(score_diagnostics, component_summary, factor_tasks,
                                     scaled_z_tolerance, task_graph_diagnostics = NULL) {
  if (!is.null(task_graph_diagnostics) && nrow(task_graph_diagnostics)) {
    required_graph_cols <- c("source_is_row_order_default", "source_is_explicit")
    missing_graph_cols <- setdiff(required_graph_cols, names(task_graph_diagnostics))
    if (length(missing_graph_cols)) {
      .stop("ltmle_exact task graph diagnostics are missing dependency validation columns: ",
            paste(missing_graph_cols, collapse = ", "))
    }
    if (any(as.logical(task_graph_diagnostics$source_is_row_order_default) %in% TRUE)) {
      .stop("ltmle_exact task graph uses row-order pseudo-outcome dependencies.")
    }
    if (!all(as.logical(task_graph_diagnostics$source_is_explicit) %in% TRUE)) {
      .stop("ltmle_exact task graph contains non-explicit pseudo-outcome dependencies.")
    }
  }
  chk <- .ltmle_exact_required_checks(score_diagnostics, component_summary, factor_tasks,
                                      scaled_z_tolerance)
  fatal_registry_mismatch <- nrow(chk$missing_tasks) > 0L || nrow(chk$extra_scores) > 0L
  if (isTRUE(fatal_registry_mismatch)) {
    bad_cols <- intersect(c("component", "t", "node", "process_type", "role", "score_after"),
                          names(score_diagnostics))
    bad_nodes <- score_diagnostics[!score_diagnostics$score_equation_solved, bad_cols, drop = FALSE]
    msg <- "Full ltmle_exact empirical-equation registry is inconsistent."
    if (nrow(chk$missing_tasks)) {
      msg <- paste0(msg, " Missing required empirical equations: ",
                    paste(apply(chk$missing_tasks[1:min(5, nrow(chk$missing_tasks)), ], 1, paste, collapse = "/"), collapse = "; "), ".")
    }
    if (nrow(chk$extra_scores)) {
      msg <- paste0(msg, " Score diagnostics contain tasks not in registry: ",
                    paste(apply(chk$extra_scores[1:min(5, nrow(chk$extra_scores)), ], 1, paste, collapse = "/"), collapse = "; "), ".")
    }
    if (nrow(bad_nodes)) {
      msg <- paste0(msg, " First unsolved score: ",
                    paste(bad_nodes[1L, ], collapse = " / "))
    }
    .stop(msg)
  }
  chk
}

.ltmle_exact_pick_numeric <- function(dat, cols, default = 0) {
  cols <- intersect(as.character(cols), names(dat))
  if (!length(cols)) return(rep(as.numeric(default)[1L], nrow(dat)))
  x <- dat[[cols[1L]]]
  if (is.numeric(x)) return(as.numeric(x))
  as.numeric(factor(x))
}

.ltmle_exact_normalize_censoring_vars <- function(censoring_vars, T) {
  if (is.null(censoring_vars) || !length(censoring_vars)) return(NULL)
  if (is.list(censoring_vars) && !is.data.frame(censoring_vars)) {
    visit <- as.character(censoring_vars$visit %||% character(0))
    final <- censoring_vars$final %||% NA_character_
    final <- if (length(final) && !is.na(final[1L])) as.character(final[1L]) else NA_character_
    if (length(visit) != T) .stop("visit censoring vars must have length T. Got ", length(visit), ".")
    return(list(visit = visit, final = final))
  }
  censoring_vars <- as.character(censoring_vars)
  visit <- grep("^R_visit", censoring_vars, value = TRUE)
  if (length(visit) >= T) {
    visit_num <- suppressWarnings(as.integer(sub("^R_visit", "", visit)))
    visit <- visit[order(ifelse(is.na(visit_num), seq_along(visit), visit_num))]
    final <- intersect(censoring_vars, "R_final")
    final <- if (length(final)) final[1L] else NA_character_
    return(list(visit = visit[seq_len(T)], final = final))
  }
  if (length(censoring_vars) == T) return(list(visit = censoring_vars, final = NA_character_))
  if (length(censoring_vars) == T + 1L) {
    final <- intersect(censoring_vars, "R_final")
    if (!length(final)) final <- censoring_vars[T + 1L]
    visit <- setdiff(censoring_vars, final[1L])
    return(list(visit = visit[seq_len(T)], final = final[1L]))
  }
  .stop("censoring_vars must have length T or T+1. Got ", length(censoring_vars), " with T=", T, ".")
}

.ltmle_exact_canonicalize_node_spec <- function(dat_wide, node_spec, T) {
  if (is.null(node_spec)) {
    validate_wide_data(dat_wide, T)
    return(dat_wide)
  }

  baseline_vars <- as.character(node_spec$baseline_vars %||% character(0))
  treatment_vars <- as.character(node_spec$treatment_vars %||% character(0))
  mediator1_vars <- as.character(node_spec$mediator1_vars %||% character(0))
  mediator2_vars <- as.character(node_spec$mediator2_vars %||% character(0))
  L_blocks <- node_spec$L_blocks %||% list()
  L_lag_init <- node_spec$L_lag_init %||% list()
  L_lag_init_cols <- unlist(
    L_lag_init[vapply(L_lag_init, is.character, logical(1))],
    use.names = FALSE
  )
  outcome_vars <- as.character(node_spec$outcome_vars %||% character(0))
  censoring_vars <- as.character(node_spec$censoring_vars %||% character(0))

  if (length(treatment_vars) == 1L) treatment_vars <- rep(treatment_vars, T)
  if (length(treatment_vars) != T) .stop("node_spec$treatment_vars must have length 1 or T.")
  if (length(mediator1_vars) != T) .stop("node_spec$mediator1_vars must have length T.")
  if (length(mediator2_vars) != T) .stop("node_spec$mediator2_vars must have length T.")
  outcome_type <- if (length(outcome_vars) == 1L) "terminal_only" else "longitudinal"
  if (!length(outcome_vars) %in% c(1L, T)) .stop("node_spec$outcome_vars must have length 1 or T.")

  L_by_time <- vector("list", T)
  if (length(L_blocks)) {
    if (is.atomic(L_blocks)) L_blocks <- list(L = as.character(L_blocks))
    for (tt in seq_len(T)) {
      cols_tt <- vapply(L_blocks, function(x) as.character(x)[tt], character(1))
      L_by_time[[tt]] <- cols_tt
    }
  } else {
    for (tt in seq_len(T)) L_by_time[[tt]] <- character(0)
  }

  required <- unique(c(
    baseline_vars, treatment_vars, mediator1_vars, mediator2_vars,
    unlist(L_by_time, use.names = FALSE), L_lag_init_cols, outcome_vars, censoring_vars
  ))
  assert_cols(dat_wide, required, "node_spec ltmle_exact input")

  out <- data.frame(.row_id = seq_len(nrow(dat_wide)), stringsAsFactors = FALSE)
  for (base_nm in baseline_vars) {
    out[[base_nm]] <- .ltmle_exact_pick_numeric(dat_wide, base_nm)
  }
  for (init_col in unique(as.character(L_lag_init_cols))) {
    if (nzchar(init_col) && !init_col %in% names(out)) {
      out[[init_col]] <- .ltmle_exact_pick_numeric(dat_wide, init_col)
    }
  }
  out$W1 <- .ltmle_exact_pick_numeric(dat_wide, baseline_vars[1L] %||% character(0))
  out$W2 <- .ltmle_exact_pick_numeric(dat_wide, baseline_vars[2L] %||% character(0))
  out$Y0 <- .ltmle_exact_pick_numeric(dat_wide, baseline_vars[3L] %||% baseline_vars[1L] %||% character(0))
  out$M1_0 <- as.numeric(dat_wide[[mediator1_vars[1L]]])
  out$M2_0 <- as.numeric(dat_wide[[mediator2_vars[1L]]])
  for (tt in seq_len(T)) {
    idx <- tt - 1L
    out[[paste0("A", idx)]] <- as.numeric(dat_wide[[treatment_vars[tt]]])
    out[[paste0("M1_", idx)]] <- as.numeric(dat_wide[[mediator1_vars[tt]]])
    out[[paste0("M2_", idx)]] <- as.numeric(dat_wide[[mediator2_vars[tt]]])
    L_cols_tt <- L_by_time[[tt]]
    for (jj in seq_along(L_cols_tt)) {
      L_name <- names(L_blocks)[jj] %||% paste0("L", jj)
      if (!is.character(L_name) || is.na(L_name) || !nzchar(L_name)) L_name <- paste0("L", jj)
      L_node <- if (startsWith(L_name, "L_") || identical(L_name, "L")) L_name else paste0("L_", L_name)
      out[[paste0(L_node, "_", idx)]] <- as.numeric(dat_wide[[L_cols_tt[jj]]])
    }
    if (length(L_cols_tt) <= 1L) {
      out[[paste0("L", idx)]] <- if (length(L_cols_tt)) as.numeric(dat_wide[[L_cols_tt[1L]]]) else 0
    }
    if (identical(outcome_type, "longitudinal")) {
      out[[paste0("Y_", tt)]] <- as.numeric(dat_wide[[outcome_vars[tt]]])
    } else if (tt == T) {
      out[[paste0("Y_", tt)]] <- as.numeric(dat_wide[[outcome_vars[1L]]])
    }
  }
  for (R_col in censoring_vars) out[[R_col]] <- as.numeric(dat_wide[[R_col]])
  out$.row_id <- NULL
  attr(out, "ltmle_exact_node_spec") <- c(node_spec, list(outcome_type = outcome_type))
  assert_no_na(out, names(out), "canonicalized node_spec ltmle_exact input")
  out
}

.ltmle_exact_node_spec_to_long <- function(dat_wide, T, node_spec = NULL) {
  if (is.null(node_spec)) {
    return(.ltmle_exact_add_terms(wide_to_long(dat_wide, T)))
  }
  L_nodes <- .ltmle_exact_L_nodes(node_spec)
  outcome_type <- .ltmle_exact_outcome_type(node_spec)
  n <- nrow(dat_wide)
  rows <- vector("list", T)
  for (tt in seq_len(T)) {
    idx <- tt - 1L
    row <- data.frame(
      id = seq_len(n),
      t = tt,
      W1 = as.numeric(dat_wide$W1),
      W2 = as.numeric(dat_wide$W2),
      Y0 = as.numeric(dat_wide$Y0),
      M1_0 = as.numeric(dat_wide$M1_0),
      M2_0 = as.numeric(dat_wide$M2_0),
      A = as.numeric(dat_wide[[paste0("A", idx)]]),
      M1 = as.numeric(dat_wide[[paste0("M1_", idx)]]),
      M2 = as.numeric(dat_wide[[paste0("M2_", idx)]]),
      Y = NA_real_,
      M1_lag = if (tt == 1L) as.numeric(dat_wide$M1_0) else as.numeric(dat_wide[[paste0("M1_", tt - 2L)]]),
      M2_lag = if (tt == 1L) as.numeric(dat_wide$M2_0) else as.numeric(dat_wide[[paste0("M2_", tt - 2L)]]),
      stringsAsFactors = FALSE
    )
    for (base_nm in as.character(node_spec$baseline_vars %||% character(0))) {
      if (base_nm %in% names(dat_wide)) row[[base_nm]] <- as.numeric(dat_wide[[base_nm]])
    }
    for (L_node in L_nodes) {
      cur_col <- paste0(L_node, "_", idx)
      if (!cur_col %in% names(dat_wide) && identical(L_node, "L")) cur_col <- paste0("L", idx)
      if (!cur_col %in% names(dat_wide)) .stop("Missing post-mediator covariate column for ", L_node, " t=", tt)
      row[[L_node]] <- as.numeric(dat_wide[[cur_col]])
      lag_nm <- .ltmle_exact_lag_name(L_node)
      if (tt == 1L) {
        init_spec <- node_spec$L_lag_init[[sub("^L_", "", L_node)]] %||% 0
        if (is.character(init_spec) && length(init_spec) == 1L) {
          if (!init_spec %in% names(dat_wide)) .stop("Missing L_lag_init column for ", L_node, ": ", init_spec)
          row[[lag_nm]] <- as.numeric(dat_wide[[init_spec]])
        } else if (is.numeric(init_spec) && length(init_spec) == 1L) {
          row[[lag_nm]] <- rep(as.numeric(init_spec), n)
        } else {
          .stop("Invalid L_lag_init for ", L_node)
        }
      } else {
        lag_col <- paste0(L_node, "_", tt - 2L)
        if (!lag_col %in% names(dat_wide) && identical(L_node, "L")) lag_col <- paste0("L", tt - 2L)
        row[[lag_nm]] <- as.numeric(dat_wide[[lag_col]])
      }
    }
    if (identical(outcome_type, "longitudinal")) {
      row$Y <- as.numeric(dat_wide[[paste0("Y_", tt)]])
      row$Y_lag <- if (tt == 1L) as.numeric(dat_wide$Y0) else as.numeric(dat_wide[[paste0("Y_", tt - 1L)]])
    } else if (tt == T) {
      row$Y <- as.numeric(dat_wide[[paste0("Y_", T)]])
    }
    rows[[tt]] <- row
  }
  long <- do.call(rbind, rows)
  attr(long, "ltmle_exact_node_spec") <- node_spec
  .ltmle_exact_add_terms(long)
}

.ltmle_exact_attach_censoring_to_long <- function(long, dat_wide, T, censoring_vars) {
  censoring_vars <- .ltmle_exact_normalize_censoring_vars(censoring_vars, T)
  if (is.null(censoring_vars)) return(long)
  censor_cols <- c(censoring_vars$visit, censoring_vars$final)
  censor_cols <- censor_cols[!is.na(censor_cols)]
  assert_cols(dat_wide, censor_cols, "ltmle_exact censoring input")
  for (R_col in censoring_vars$visit) long[[R_col]] <- NA_real_
  for (tt in seq_len(T)) {
    idx <- long$t == tt
    long[[censoring_vars$visit[tt]]][idx] <- as.numeric(dat_wide[[censoring_vars$visit[tt]]])[long$id[idx]]
  }
  final_var <- censoring_vars$final %||% NA_character_
  if (is.character(final_var) && length(final_var) == 1L && !is.na(final_var)) {
    long[[final_var]] <- as.numeric(dat_wide[[final_var]])[long$id]
  }
  long
}

.ltmle_exact_make_folds <- function(n, V, seed) {
  V <- as.integer(V)
  if (!is.finite(V) || V < 2L) .stop("V must be at least 2 for ltmle_exact cross-fitting.")
  V <- min(V, n)
  if (!is.null(seed)) set.seed(seed)
  sample(rep(seq_len(V), length.out = n))
}

.ltmle_exact_crossfit_initial_nuisance <- function(dat, node_spec, T, fold_id, learner, sl_library, ...) {
  fold_id <- as.integer(fold_id)
  if (length(fold_id) != nrow(dat)) .stop("fold_id must have one entry per subject.")
  if (any(!is.finite(fold_id)) || length(unique(fold_id)) < 2L) {
    .stop("fold_id must contain at least two folds for ltmle_exact cross-fitting.")
  }
  list(
    used = identical(learner, "sl"),
    fold_id = fold_id,
    V = length(unique(fold_id)),
    sl_library = sl_library
  )
}

.ltmle_exact_stack_out_of_fold_predictions <- function(crossfit_result) {
  crossfit_result
}

fit_ltmle_exact <- function(dat = NULL, dat_wide = NULL, node_spec = NULL, T,
                            reg_a, reg_as,
                            learner = c("glm", "sl"),
                            sl_library = c("SL.glm", "SL.mean"),
                            Q_model = c("correct", "wrong"),
                            seed = 202405L,
                            probability_bounds = c(0.01, 0.99),
                            truncation_enabled = TRUE,
                            truncation_policy = "quantile",
                            truncation_quantile_lower = 0.01,
                            truncation_quantile_upper = 0.99,
                            truncation_target = "clever_covariate_H",
                            ltmle_exact_density_ratio_mc_n = 2000L,
                            ltmle_exact_law_integration_n = 5L,
                            y_bounds_mode = c("train_fold", "fixed"),
                            y_bounds = NULL,
                            treat_mech = c("observational", "baseline_rct"),
                            p_rct = 0.5,
                            censoring_vars = NULL,
                            censoring_mech = c("none", "estimated"),
                            score_tolerance = 1e-4,
                            component_tolerance = 1e-5,
                            scaled_z_tolerance = 3,
                            V = 5L,
                            fold_id = NULL,
                            diagnostics_level = c("summary", "full"),
                            compute_auxiliary_nuisance_gcomp = FALSE,
                            auxiliary_nuisance_gcomp_mc_n = 2000L,
                            force_epsilon_zero = FALSE,
                            component_subset = NULL,
                            verbose = FALSE,
                            ...) {
  .ltmle_exact_assert_no_legacy_truncation_args(list(...))
  learner <- match.arg(learner)
  Q_model <- match.arg(Q_model)
  diagnostics_level <- match.arg(diagnostics_level)
  compute_auxiliary_nuisance_gcomp <- isTRUE(compute_auxiliary_nuisance_gcomp)
  auxiliary_nuisance_gcomp_mc_n <- as.integer(auxiliary_nuisance_gcomp_mc_n)
  if (!is.finite(auxiliary_nuisance_gcomp_mc_n) || auxiliary_nuisance_gcomp_mc_n < 1L) {
    .stop("auxiliary_nuisance_gcomp_mc_n must be an integer >= 1.")
  }
  ltmle_exact_density_ratio_mc_n <- .ltmle_exact_normalize_density_ratio_mc_n(
    ltmle_exact_density_ratio_mc_n
  )
  ltmle_exact_law_integration_n <- as.integer(ltmle_exact_law_integration_n)
  if (!is.finite(ltmle_exact_law_integration_n) || ltmle_exact_law_integration_n < 3L) {
    .stop("ltmle_exact_law_integration_n must be an integer >= 3.")
  }
  y_bounds_mode <- match.arg(y_bounds_mode)
  if (length(treat_mech) == 1L && identical(treat_mech, "known_randomized")) treat_mech <- "baseline_rct"
  treat_mech <- match.arg(treat_mech)
  censoring_mech <- match.arg(censoring_mech)
  truncation <- .ltmle_exact_normalize_truncation(
    truncation_enabled = truncation_enabled,
    truncation_policy = truncation_policy,
    truncation_quantile_lower = truncation_quantile_lower,
    truncation_quantile_upper = truncation_quantile_upper,
    truncation_target = truncation_target
  )
  estimator_variant <- if (isTRUE(truncation$enabled)) {
    "ltmle_exact_quantile_truncated"
  } else {
    "ltmle_exact_untruncated"
  }
  if (is.null(dat_wide)) dat_wide <- dat
  if (is.null(dat_wide)) .stop("Provide dat or dat_wide.")
  dat_wide <- .ltmle_exact_canonicalize_node_spec(dat_wide, node_spec, T)
  if (!is.null(node_spec)) {
    node_spec <- attr(dat_wide, "ltmle_exact_node_spec") %||% node_spec
  }
  if (identical(learner, "sl")) {
    fold_id <- fold_id %||% .ltmle_exact_make_folds(nrow(dat_wide), V, seed)
    dat_wide$.fold_id <- fold_id
    crossfit_info <- .ltmle_exact_crossfit_initial_nuisance(
      dat = dat_wide,
      node_spec = node_spec,
      T = T,
      fold_id = fold_id,
      learner = learner,
      sl_library = sl_library
    )
  } else {
    crossfit_info <- list(used = FALSE, fold_id = NULL)
  }
  censoring_vars <- .ltmle_exact_normalize_censoring_vars(censoring_vars, T)
  if (identical(censoring_mech, "none")) censoring_vars <- NULL
  if (identical(censoring_mech, "estimated") && is.null(censoring_vars)) {
    .stop("censoring_mech='estimated' requires censoring_vars.")
  }
  if (is.null(node_spec)) validate_wide_data(dat_wide, T)
  reg_a <- normalize_regimen(reg_a, T, "reg_a")
  reg_as <- normalize_regimen(reg_as, T, "reg_as")
  validate_baseline_rct_regimens(reg_a, reg_as, treat_mech)
  probability_bounds <- as.numeric(probability_bounds)
  if (length(probability_bounds) != 2L || any(!is.finite(probability_bounds)) ||
      probability_bounds[1L] <= 0 || probability_bounds[2L] >= 1 ||
      probability_bounds[1L] >= probability_bounds[2L]) {
    .stop("probability_bounds must satisfy 0 < lower < upper < 1.")
  }
  long <- .ltmle_exact_node_spec_to_long(dat_wide, T, node_spec)
  if (identical(learner, "sl")) long$.fold_id <- fold_id[long$id]
  for (nm in as.character(node_spec$baseline_vars %||% character(0))) {
    if (nm %in% names(dat_wide) && !nm %in% names(long)) {
      long[[nm]] <- as.numeric(dat_wide[[nm]])[long$id]
    }
  }
  long <- .ltmle_exact_attach_censoring_to_long(long, dat_wide, T, censoring_vars)
  treatment_models <- .ltmle_exact_fit_treatment_models(long, T, learner, sl_library, treat_mech, p_rct,
                                                        node_spec = node_spec,
                                                        fold_id = if (identical(learner, "sl")) long$.fold_id else NULL)
  censoring_models <- .ltmle_exact_fit_censoring_models(long, T, censoring_vars, learner, sl_library,
                                                        node_spec = node_spec,
                                                        fold_id = if (identical(learner, "sl")) long$.fold_id else NULL)
  models <- .ltmle_exact_fit_node_models(long, T, learner, sl_library, Q_model, node_spec = node_spec,
                                         fold_id = if (identical(learner, "sl")) long$.fold_id else NULL)
  mediator_density_models <- .ltmle_exact_mediator_density_summary(models)
  if (identical(learner, "sl") && !isTRUE(mediator_density_models$crossfit_used)) {
    .stop("Super Learner ltmle_exact requires cross-fitted mediator density models.")
  }
  registry <- .ltmle_exact_component_registry(reg_a, reg_as, T)
  world_specs <- .ltmle_exact_world_spec(reg_a, reg_as, T)
  factor_tasks <- .ltmle_exact_factor_tasks(reg_a, reg_as, T, node_spec = node_spec)
  component_subset <- as.character(component_subset %||% character(0))
  component_subset <- component_subset[nzchar(component_subset)]
  if (length(component_subset)) {
    unknown_components <- setdiff(component_subset, registry$component)
    if (length(unknown_components)) {
      .stop("Unknown ltmle_exact component_subset value(s): ", paste(unknown_components, collapse = ", "))
    }
    registry <- registry[match(component_subset, registry$component), , drop = FALSE]
    world_specs <- world_specs[world_specs$component %in% component_subset, , drop = FALSE]
    factor_tasks <- factor_tasks[factor_tasks$component %in% component_subset, , drop = FALSE]
  }
  .ltmle_exact_validate_stochastic_law_registry(factor_tasks, T, node_spec)
  means <- setNames(rep(NA_real_, nrow(registry)), registry$component)
  means_initial <- setNames(rep(NA_real_, nrow(registry)), registry$component)
  subject_matrix <- matrix(NA_real_, nrow = nrow(dat_wide), ncol = nrow(registry))
  colnames(subject_matrix) <- registry$component
  subject_initial_matrix <- matrix(NA_real_, nrow = nrow(dat_wide), ncol = nrow(registry))
  colnames(subject_initial_matrix) <- registry$component
  component_eif_matrix <- matrix(NA_real_, nrow = nrow(dat_wide), ncol = nrow(registry))
  colnames(component_eif_matrix) <- registry$component
  score_rows <- list()
  eif_term_rows <- list()
  density_ratio_rows <- list()
  truncation_rows <- list()
  clever_decomposition_rows <- list()
  mc_rows <- list()
  law_task_integration_rows <- list()
  branch_history_rows <- list()
  branch_handoff_rows <- list()
  law_integration_rows <- list()
  source_eval_training_rows <- list()
  cross_regimen_source_boundary_rows <- list()
  source_pseudooutcome_independent_rows <- list()
  targeted_continuation_storage_rows <- list()
  observed_vs_source_rows <- list()
  dedicated_L_training_row_trace_rows <- list()
  dedicated_L_prediction_row_trace_rows <- list()
  dedicated_L_coefficient_contribution_rows <- list()
  dedicated_L_prediction_weighting_collapse_rows <- list()
  dedicated_L_mediator_path_weight_trace_rows <- list()
  dedicated_L_training_response_oracle_rows <- list()
  root_plugin_rows <- list()
  root_vs_branch_rows <- list()
  separate_clever_covariate_identity_rows <- list()
  second_M2_marginal_reference_rows <- list()
  task_graph_rows <- list()

  for (ii in seq_len(nrow(registry))) {
    comp <- registry$component[ii]
    sp <- world_specs[world_specs$component == comp, , drop = FALSE]
    component_tasks <- factor_tasks[factor_tasks$component == comp, , drop = FALSE]
    one <- .ltmle_exact_component_fit(
      dat_wide = dat_wide,
      long = long,
      T = T,
      spec = sp,
      models = models,
      treatment_models = treatment_models,
      learner = learner,
      sl_library = sl_library,
      Q_model = Q_model,
      seed = seed + ii * 1009L,
      probability_bounds = probability_bounds,
      truncation = truncation,
      estimator_variant = estimator_variant,
      y_bounds_mode = y_bounds_mode,
      y_bounds = y_bounds,
      score_tolerance = score_tolerance,
      treat_mech = treat_mech,
      p_rct = p_rct,
      censoring_models = censoring_models,
      component_tasks = component_tasks,
      node_spec = node_spec,
      ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
      ltmle_exact_law_integration_n = ltmle_exact_law_integration_n,
      auxiliary_nuisance_gcomp_mc_n = auxiliary_nuisance_gcomp_mc_n,
      force_epsilon_zero = force_epsilon_zero,
      diagnostics_level = diagnostics_level,
      verbose = verbose
    )
    means[comp] <- one$mean
    means_initial[comp] <- one$initial_mean
    subject_matrix[, comp] <- one$subject_q
    subject_initial_matrix[, comp] <- one$subject_q_initial
    component_eif_matrix[, comp] <- one$component_D
    score_rows[[length(score_rows) + 1L]] <- one$score_diagnostics
    eif_term_rows[[length(eif_term_rows) + 1L]] <- one$component_eif_terms
    density_ratio_rows[[length(density_ratio_rows) + 1L]] <- one$density_ratio_diagnostics
    truncation_rows[[length(truncation_rows) + 1L]] <- one$truncation_diagnostics
    clever_decomposition_rows[[length(clever_decomposition_rows) + 1L]] <-
      one$clever_covariate_decomposition_diagnostics
    mc_rows[[length(mc_rows) + 1L]] <- one$mc_integration_diagnostics
    law_task_integration_rows[[length(law_task_integration_rows) + 1L]] <- one$law_task_integration_diagnostics
    branch_history_rows[[length(branch_history_rows) + 1L]] <- one$branch_history_state_check
    branch_handoff_rows[[length(branch_handoff_rows) + 1L]] <- one$branch_state_handoff_trace
    law_integration_rows[[length(law_integration_rows) + 1L]] <- one$law_integration_diagnostics
    source_eval_training_rows[[length(source_eval_training_rows) + 1L]] <- one$source_eval_training_path_check
    cross_regimen_source_boundary_rows[[length(cross_regimen_source_boundary_rows) + 1L]] <-
      one$cross_regimen_source_boundary_trace
    source_pseudooutcome_independent_rows[[length(source_pseudooutcome_independent_rows) + 1L]] <-
      one$source_pseudooutcome_independent_recursion_check
    targeted_continuation_storage_rows[[length(targeted_continuation_storage_rows) + 1L]] <-
      one$targeted_continuation_storage_check
    dedicated_L_training_row_trace_rows[[length(dedicated_L_training_row_trace_rows) + 1L]] <-
      one$dedicated_L_transition_L_fit_training_row_trace
    dedicated_L_prediction_row_trace_rows[[length(dedicated_L_prediction_row_trace_rows) + 1L]] <-
      one$dedicated_L_transition_L_fit_prediction_row_trace
    dedicated_L_coefficient_contribution_rows[[length(dedicated_L_coefficient_contribution_rows) + 1L]] <-
      one$dedicated_L_transition_L_fit_coefficient_contribution_check
    dedicated_L_prediction_weighting_collapse_rows[[length(dedicated_L_prediction_weighting_collapse_rows) + 1L]] <-
      one$dedicated_L_transition_prediction_weighting_collapse_check
    dedicated_L_mediator_path_weight_trace_rows[[length(dedicated_L_mediator_path_weight_trace_rows) + 1L]] <-
      one$dedicated_L_transition_mediator_path_weight_trace
    dedicated_L_training_response_oracle_rows[[length(dedicated_L_training_response_oracle_rows) + 1L]] <-
      one$dedicated_L_transition_training_response_oracle_alignment_check
    observed_vs_source_rows[[length(observed_vs_source_rows) + 1L]] <- one$observed_vs_source_row_role_check
    root_plugin_rows[[length(root_plugin_rows) + 1L]] <- one$root_plugin_diagnostics
    root_vs_branch_rows[[length(root_vs_branch_rows) + 1L]] <- one$root_vs_targeted_branch_recursion_check
    separate_clever_covariate_identity_rows[[length(separate_clever_covariate_identity_rows) + 1L]] <-
      one$separate_clever_covariate_identity_check
    second_M2_marginal_reference_rows[[length(second_M2_marginal_reference_rows) + 1L]] <-
      one$second_M2_marginal_reference_check
    graph_order <- one$task_graph$reverse_topological_order
    graph_tasks <- one$task_graph$tasks[graph_order]
    previous_id <- c(NA_character_, graph_order[-length(graph_order)])
    observed_source <- vapply(
      graph_tasks,
      function(x) x$observed_pseudooutcome_source_task_id %||% NA_character_,
      character(1)
    )
    source_is_explicit <- vapply(
      graph_tasks,
      function(x) isTRUE(x$source_is_explicit %||% FALSE) &&
        nzchar(x$observed_pseudooutcome_source_task_id %||% ""),
      logical(1)
    )
    source_matches_previous <- !is.na(previous_id) & observed_source == previous_id
    task_graph_rows[[length(task_graph_rows) + 1L]] <- data.frame(
      component = comp,
      task_id = graph_order,
      parent_task_id = vapply(graph_tasks, function(x) x$parent_task_id %||% NA_character_, character(1)),
      observed_pseudooutcome_source_task_id = observed_source,
      source_is_explicit = source_is_explicit,
      source_is_row_order_default = source_matches_previous & !source_is_explicit,
      node = vapply(graph_tasks, function(x) x$node, character(1)),
      time = vapply(graph_tasks, function(x) as.integer(x$t), integer(1)),
      process_type = vapply(graph_tasks, function(x) x$process_type, character(1)),
      world_type = vapply(graph_tasks, function(x) x$world_type %||% NA_character_, character(1)),
      assigned_treatment_regimen_label = vapply(
        graph_tasks,
        function(x) x$assigned_treatment_regimen_label %||% NA_character_,
        character(1)
      ),
      conditioning_history_type = vapply(
        graph_tasks,
        function(x) x$conditioning_history_type %||% NA_character_,
        character(1)
      ),
      continuation_target_label = vapply(
        graph_tasks,
        .ltmle_exact_continuation_target_label,
        character(1)
      ),
      expected_cache_key = vapply(
        graph_tasks,
        .ltmle_exact_semantic_cache_key,
        character(1)
      ),
      mediator_role = vapply(graph_tasks, function(x) x$mediator_role %||% NA_character_, character(1)),
      source_boundary_type = vapply(graph_tasks, function(x) x$source_boundary_type %||% "pure", character(1)),
      source_boundary_direction = vapply(graph_tasks, function(x) x$source_boundary_direction %||% "not_applicable", character(1)),
      source_boundary_eval_state = vapply(graph_tasks, function(x) x$source_boundary_eval_state %||% "not_applicable", character(1)),
      source_boundary_outcome_history_state = vapply(
        graph_tasks,
        function(x) x$source_boundary_outcome_history_state %||% "not_applicable",
        character(1)
      ),
      source_boundary_m1_history_state = vapply(
        graph_tasks,
        function(x) x$source_boundary_m1_history_state %||% "not_applicable",
        character(1)
      ),
      source_boundary_m2_history_state = vapply(
        graph_tasks,
        function(x) x$source_boundary_m2_history_state %||% "not_applicable",
        character(1)
      ),
      source_boundary_auxiliary_mediator_history_state = vapply(
        graph_tasks,
        function(x) x$source_boundary_auxiliary_mediator_history_state %||% "not_applicable",
        character(1)
      ),
      source_boundary_outer_regimen = vapply(graph_tasks, function(x) x$source_boundary_outer_regimen %||% NA_character_, character(1)),
      source_boundary_m1_regimen = vapply(graph_tasks, function(x) x$source_boundary_m1_regimen %||% NA_character_, character(1)),
      source_boundary_m2_regimen = vapply(graph_tasks, function(x) x$source_boundary_m2_regimen %||% NA_character_, character(1)),
      source_boundary_provenance_rule = vapply(graph_tasks, function(x) x$source_boundary_provenance_rule %||% "not_applicable", character(1)),
      graph_used = isTRUE(one$task_graph$graph_used),
      stringsAsFactors = FALSE
    )
  }

  # Final ltmle_exact estimate is the targeted deterministic-root substitution
  # estimate. Auxiliary fitted-nuisance g-computation must never define it.
  root_plugin_targeted_means <- means
  integrated_worldmean_recursion <- NULL
  if (compute_auxiliary_nuisance_gcomp) {
    integrated_worldmean_recursion <- .ltmle_exact_integrated_worldmeans_from_models(
      dat_wide = dat_wide,
      models = models,
      reg_a = reg_a,
      reg_as = reg_as,
      T = T,
      node_spec = node_spec,
      B_mc = auxiliary_nuisance_gcomp_mc_n,
      seed = seed + 7919L
    )
  }

  score_diagnostics <- .ltmle_exact_rbind_fill(score_rows)
  non_null_eif <- eif_term_rows[!vapply(eif_term_rows, is.null, logical(1))]
  component_eif_terms <- if (length(non_null_eif)) do.call(rbind, non_null_eif) else NULL
  density_ratio_diagnostics <- .ltmle_exact_normalize_density_ratio_diagnostics(
    .ltmle_exact_rbind_fill(density_ratio_rows),
    requested_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
    requested_law_integration_n = ltmle_exact_law_integration_n
  )
  truncation_diagnostics <- .ltmle_exact_rbind_fill(truncation_rows)
  clever_covariate_decomposition_diagnostics <- .ltmle_exact_rbind_fill(clever_decomposition_rows)
  separate_clever_covariate_identity_check <- .ltmle_exact_rbind_fill(
    separate_clever_covariate_identity_rows
  )
  second_M2_marginal_reference_check <- .ltmle_exact_rbind_fill(
    second_M2_marginal_reference_rows
  )
  mc_integration_diagnostics <- .ltmle_exact_rbind_fill(mc_rows)
  non_empty_law_integration <- law_task_integration_rows[
    vapply(law_task_integration_rows, function(x) !is.null(x) && nrow(x) > 0L, logical(1))
  ]
  law_task_integration_diagnostics <- if (length(non_empty_law_integration)) {
    .ltmle_exact_rbind_fill(non_empty_law_integration)
  } else {
    data.frame()
  }
  non_empty_branch_history <- branch_history_rows[
    vapply(branch_history_rows, function(x) !is.null(x) && nrow(x) > 0L, logical(1))
  ]
  branch_history_state_check <- if (length(non_empty_branch_history)) {
    .ltmle_exact_rbind_fill(non_empty_branch_history)
  } else {
    data.frame()
  }
  non_empty_handoff <- branch_handoff_rows[
    vapply(branch_handoff_rows, function(x) !is.null(x) && nrow(x) > 0L, logical(1))
  ]
  branch_state_handoff_trace <- if (length(non_empty_handoff)) {
    .ltmle_exact_rbind_fill(non_empty_handoff)
  } else {
    data.frame()
  }
  separate_product_join_diagnostics <- .ltmle_exact_separate_product_join_diagnostics(
    branch_state_handoff_trace = branch_state_handoff_trace,
    detail = if (identical(diagnostics_level, "full")) "row" else "event"
  )
  non_empty_law_eval <- law_integration_rows[
    vapply(law_integration_rows, function(x) !is.null(x) && nrow(x) > 0L, logical(1))
  ]
  law_integration_diagnostics <- if (length(non_empty_law_eval)) {
    .ltmle_exact_rbind_fill(non_empty_law_eval)
  } else {
    data.frame()
  }
  non_empty_source_eval <- source_eval_training_rows[
    vapply(source_eval_training_rows, function(x) !is.null(x) && nrow(x) > 0L, logical(1))
  ]
  source_eval_training_path_check <- if (length(non_empty_source_eval)) {
    .ltmle_exact_rbind_fill(non_empty_source_eval)
  } else {
    data.frame()
  }
  non_empty_cross_boundary <- cross_regimen_source_boundary_rows[
    vapply(cross_regimen_source_boundary_rows, function(x) !is.null(x) && nrow(x) > 0L, logical(1))
  ]
  cross_regimen_source_boundary_trace <- if (length(non_empty_cross_boundary)) {
    .ltmle_exact_rbind_fill(non_empty_cross_boundary)
  } else {
    data.frame()
  }
  non_empty_source_independent <- source_pseudooutcome_independent_rows[
    vapply(source_pseudooutcome_independent_rows, function(x) !is.null(x) && nrow(x) > 0L, logical(1))
  ]
  source_pseudooutcome_independent_recursion_check <- if (length(non_empty_source_independent)) {
    .ltmle_exact_rbind_fill(non_empty_source_independent)
  } else {
    data.frame()
  }
  non_empty_targeted_storage <- targeted_continuation_storage_rows[
    vapply(targeted_continuation_storage_rows, function(x) !is.null(x) && nrow(x) > 0L, logical(1))
  ]
  targeted_continuation_storage_check <- if (length(non_empty_targeted_storage)) {
    .ltmle_exact_rbind_fill(non_empty_targeted_storage)
  } else {
    data.frame()
  }
  non_empty_dedicated_L_training_row_trace <- dedicated_L_training_row_trace_rows[
    vapply(dedicated_L_training_row_trace_rows, function(x) !is.null(x) && nrow(x) > 0L, logical(1))
  ]
  dedicated_L_transition_L_fit_training_row_trace <- if (length(non_empty_dedicated_L_training_row_trace)) {
    .ltmle_exact_rbind_fill(non_empty_dedicated_L_training_row_trace)
  } else {
    data.frame()
  }
  non_empty_dedicated_L_prediction_row_trace <- dedicated_L_prediction_row_trace_rows[
    vapply(dedicated_L_prediction_row_trace_rows, function(x) !is.null(x) && nrow(x) > 0L, logical(1))
  ]
  dedicated_L_transition_L_fit_prediction_row_trace <- if (length(non_empty_dedicated_L_prediction_row_trace)) {
    .ltmle_exact_rbind_fill(non_empty_dedicated_L_prediction_row_trace)
  } else {
    data.frame()
  }
  non_empty_dedicated_L_coefficient_contribution <- dedicated_L_coefficient_contribution_rows[
    vapply(dedicated_L_coefficient_contribution_rows, function(x) !is.null(x) && nrow(x) > 0L, logical(1))
  ]
  dedicated_L_transition_L_fit_coefficient_contribution_check <-
    if (length(non_empty_dedicated_L_coefficient_contribution)) {
      .ltmle_exact_rbind_fill(non_empty_dedicated_L_coefficient_contribution)
    } else {
      data.frame()
    }
  non_empty_dedicated_L_weighting_collapse <- dedicated_L_prediction_weighting_collapse_rows[
    vapply(dedicated_L_prediction_weighting_collapse_rows, function(x) !is.null(x) && nrow(x) > 0L, logical(1))
  ]
  dedicated_L_transition_prediction_weighting_collapse_check <-
    if (length(non_empty_dedicated_L_weighting_collapse)) {
      .ltmle_exact_rbind_fill(non_empty_dedicated_L_weighting_collapse)
    } else {
      data.frame()
    }
  non_empty_dedicated_L_mediator_path_weight_trace <- dedicated_L_mediator_path_weight_trace_rows[
    vapply(dedicated_L_mediator_path_weight_trace_rows, function(x) !is.null(x) && nrow(x) > 0L, logical(1))
  ]
  dedicated_L_transition_mediator_path_weight_trace <-
    if (length(non_empty_dedicated_L_mediator_path_weight_trace)) {
      .ltmle_exact_rbind_fill(non_empty_dedicated_L_mediator_path_weight_trace)
    } else {
      data.frame()
    }
  non_empty_dedicated_L_training_response_oracle <- dedicated_L_training_response_oracle_rows[
    vapply(dedicated_L_training_response_oracle_rows, function(x) !is.null(x) && nrow(x) > 0L, logical(1))
  ]
  dedicated_L_transition_training_response_oracle_alignment_check <-
    if (length(non_empty_dedicated_L_training_response_oracle)) {
      .ltmle_exact_rbind_fill(non_empty_dedicated_L_training_response_oracle)
    } else {
      data.frame()
    }
  non_empty_observed_vs_source <- observed_vs_source_rows[
    vapply(observed_vs_source_rows, function(x) !is.null(x) && nrow(x) > 0L, logical(1))
  ]
  observed_vs_source_row_role_check <- if (length(non_empty_observed_vs_source)) {
    .ltmle_exact_rbind_fill(non_empty_observed_vs_source)
  } else {
    data.frame()
  }
  non_empty_root_plugin <- root_plugin_rows[
    vapply(root_plugin_rows, function(x) !is.null(x) && nrow(x) > 0L, logical(1))
  ]
  root_plugin_diagnostics <- if (length(non_empty_root_plugin)) {
    .ltmle_exact_rbind_fill(non_empty_root_plugin)
  } else {
    data.frame()
  }
  non_empty_root_vs_branch <- root_vs_branch_rows[
    vapply(root_vs_branch_rows, function(x) !is.null(x) && nrow(x) > 0L, logical(1))
  ]
  root_vs_targeted_branch_recursion_check <- if (length(non_empty_root_vs_branch)) {
    .ltmle_exact_rbind_fill(non_empty_root_vs_branch)
  } else {
    data.frame()
  }
  task_graph_diagnostics <- do.call(rbind, task_graph_rows)
  component_eif_summary <- .ltmle_exact_component_summary(component_eif_matrix, component_tolerance)
  checks <- .ltmle_exact_fail_closed(score_diagnostics, component_eif_summary, factor_tasks,
                                     scaled_z_tolerance, task_graph_diagnostics)
  missing_law_specific <- .ltmle_exact_missing_law_specific_equations(
    score_diagnostics = score_diagnostics,
    registry = registry,
    T = T,
    node_spec = node_spec
  )
  checks$missing_law_specific_equations <- missing_law_specific
  checks$all_law_specific_equations_present <- nrow(missing_law_specific) == 0L
  if (nrow(missing_law_specific)) {
    .stop(
      "ltmle_exact missing law-specific stochastic mediator intervention empirical equations: ",
      paste(
        apply(
          missing_law_specific[seq_len(min(5L, nrow(missing_law_specific))), , drop = FALSE],
          1,
          paste,
          collapse = "/"
        ),
        collapse = "; "
      ),
      "."
    )
  }
  process_counts <- table(score_diagnostics$process_type)
  process_count_string <- paste(paste(names(process_counts), as.integer(process_counts), sep = "="), collapse = ";")
  treatment_diag <- density_ratio_diagnostics[density_ratio_diagnostics$ratio_name == "treatment likelihood factor", , drop = FALSE]
  censoring_diag <- density_ratio_diagnostics[density_ratio_diagnostics$ratio_name == "censoring factor", , drop = FALSE]
  mediator_diag <- density_ratio_diagnostics[grepl("mediator", density_ratio_diagnostics$ratio_name), , drop = FALSE]
  required_process_types <- c(
    "outcome_process",
    "post_mediator_covariate_transition"
  )
  if (any(registry$world_type == "joint")) {
    required_process_types <- c(
      required_process_types,
      "joint_stochastic_mediator_intervention_law"
    )
  }
  if (any(registry$world_type == "separate")) {
    required_process_types <- c(
      required_process_types,
      "first_mediator_stochastic_intervention_law",
      "second_mediator_stochastic_intervention_law"
    )
  }
  missing_process_types <- setdiff(required_process_types, unique(score_diagnostics$process_type))
  if (length(missing_process_types)) {
    .stop("ltmle_exact missing required process types: ", paste(missing_process_types, collapse = ", "))
  }
  # Empirical equation failures are returned as hard-gate diagnostics. They
  # should not make an otherwise completed diagnostic run look runtime-incomplete.
  if (!all(required_process_types %in% unique(score_diagnostics$process_type))) {
    .stop("ltmle_exact is missing empirical equations for required process types.")
  }
  if (any(!is.finite(component_eif_matrix))) {
    .stop("ltmle_exact component efficient influence function matrix contains non-finite values.")
  }
  task_graph_used <- isTRUE(all(task_graph_diagnostics$graph_used))
  if (!task_graph_used) .stop("ltmle_exact task graph was not used for all components.")
  fake_terminal_outcome_history_used <- identical(.ltmle_exact_outcome_type(node_spec), "terminal_only") &&
    ("Y_lag" %in% names(long) || "Y_lag" %in% unlist(lapply(task_graph_diagnostics$task_id, identity)))
  if (fake_terminal_outcome_history_used) .stop("ltmle_exact terminal-only outcome used fake intermediate outcome history.")
  multiple_L_blocks_collapsed <- .ltmle_exact_multi_L(node_spec) && any(c("L", "L_lag") %in% names(long))
  if (multiple_L_blocks_collapsed) .stop("ltmle_exact multiple post-mediator covariates were collapsed into scalar L.")
  visit_censoring_downstream <- any(c("A", "M1", "M2", "L", .ltmle_exact_L_nodes(node_spec)) %in%
                                      .ltmle_exact_visit_censoring_covariates(node_spec))
  if (visit_censoring_downstream) .stop("ltmle_exact visit censoring model uses downstream current nodes.")
  n_clever_covariate_values_truncated <- if ("n_clever_covariate_values_truncated" %in% names(score_diagnostics)) {
    sum(as.integer(score_diagnostics$n_clever_covariate_values_truncated), na.rm = TRUE)
  } else {
    0L
  }
  n_clever_covariate_values_checked <- if (nrow(score_diagnostics) && "n_rows" %in% names(score_diagnostics)) {
    sum(as.integer(score_diagnostics$n_rows), na.rm = TRUE)
  } else {
    0L
  }

  run_summary <- data.frame(
    method = "ltmle_exact",
    estimator_class = .ltmle_exact_class(),
    n_subjects = nrow(dat_wide),
    T = T,
    learner = learner,
    Q_model = Q_model,
    terminal_plugin_type = "deterministic_root",
    terminal_plugin_mc_active = FALSE,
    effective_terminal_plugin_mc_n = 1L,
    ltmle_exact_density_ratio_mc_n = as.integer(ltmle_exact_density_ratio_mc_n),
    ltmle_exact_law_integration_n = as.integer(ltmle_exact_law_integration_n),
    worldmean_plugin_type = "targeted_deterministic_root",
    worldmean_plugin_mc_active = FALSE,
    effective_worldmean_plugin_mc_n = 1L,
    auxiliary_nuisance_gcomp_computed = !is.null(integrated_worldmean_recursion),
    auxiliary_nuisance_gcomp_mc_n = if (!is.null(integrated_worldmean_recursion)) {
      integrated_worldmean_recursion$B_mc
    } else {
      NA_integer_
    },
    diagnostics_level = diagnostics_level,
    treat_mech = treat_mech,
    probability_bounds_lower = probability_bounds[1L],
    probability_bounds_upper = probability_bounds[2L],
    probability_bounding_enabled = TRUE,
    probability_bound_lower = probability_bounds[1L],
    probability_bound_upper = probability_bounds[2L],
    probability_bounding_is_truncation = FALSE,
    truncation_enabled = isTRUE(truncation$enabled),
    truncation_policy = truncation$policy,
    truncation_target = truncation$target,
    truncation_rule = if (isTRUE(truncation$enabled)) "sample_quantile" else "none",
    requested_truncation_quantile_lower = truncation$quantile_lower,
    requested_truncation_quantile_upper = truncation$quantile_upper,
    effective_truncation_quantile_lower = if (isTRUE(truncation$enabled)) truncation$quantile_lower else NA_real_,
    effective_truncation_quantile_upper = if (isTRUE(truncation$enabled)) truncation$quantile_upper else NA_real_,
    fixed_bound_truncation_enabled = FALSE,
    fixed_bound_truncation_used = FALSE,
    clever_covariate_truncation_enabled = isTRUE(truncation$enabled),
    clever_covariate_truncation_policy = truncation$policy,
    clever_covariate_truncation_target = "clever_covariate_H",
    density_ratio_factor_truncation_used = FALSE,
    n_density_ratio_values_truncated = 0L,
    fraction_density_ratio_values_truncated = 0,
    n_clever_covariate_values_truncated = as.integer(n_clever_covariate_values_truncated),
    fraction_clever_covariate_values_truncated = if (n_clever_covariate_values_checked > 0L) {
      n_clever_covariate_values_truncated / n_clever_covariate_values_checked
    } else {
      0
    },
    estimator_variant = estimator_variant,
    mediator_density_engine = mediator_density_models$engine,
    mediator_density_crossfit_used = isTRUE(mediator_density_models$crossfit_used),
    score_tolerance = score_tolerance,
    component_tolerance = component_tolerance,
    scaled_z_tolerance = scaled_z_tolerance,
    n_score_equations_total = nrow(score_diagnostics),
    n_score_equations_by_process_type = process_count_string,
    n_component_eif_columns = ncol(component_eif_matrix),
    max_abs_score_after = max(abs(score_diagnostics$score_after), na.rm = TRUE),
    max_abs_component_eif_mean = max(abs(component_eif_summary$mean_D), na.rm = TRUE),
    max_component_eif_scaled_Z = max(component_eif_summary$scaled_Z, na.rm = TRUE),
    all_score_equations_solved = all(score_diagnostics$score_equation_solved),
    all_component_equations_solved = all(component_eif_summary$component_equation_solved),
    all_required_factor_tasks_present = isTRUE(checks$all_required_factor_tasks_present),
    all_law_specific_equations_present = isTRUE(checks$all_law_specific_equations_present),
    all_required_censoring_factors_present = identical(censoring_models$type, "none") ||
      (length(censoring_models$visit_fits %||% censoring_models$fits) == T &&
         (is.na(censoring_models$vars$final %||% NA_character_) || !is.null(censoring_models$final_fit))),
    all_predictions_strictly_valid = isTRUE(checks$all_predictions_strictly_valid),
    all_density_ratios_strictly_valid = isTRUE(checks$all_density_ratios_strictly_valid) &&
      all(is.finite(density_ratio_diagnostics$max_ratio)),
    all_equations_solved = all(score_diagnostics$score_equation_solved) &&
      all(component_eif_summary$component_equation_solved) &&
      isTRUE(checks$all_required_factor_tasks_present) &&
      isTRUE(checks$all_law_specific_equations_present) &&
      max(component_eif_summary$scaled_Z, na.rm = TRUE) <= scaled_z_tolerance,
    censoring_factor = censoring_models$type,
    cross_fitting_used = isTRUE(crossfit_info$used),
    task_graph_used = task_graph_used,
    task_graph_dependencies_explicit = all(as.logical(task_graph_diagnostics$source_is_explicit) %in% TRUE),
    task_graph_row_order_dependencies_absent =
      !any(as.logical(task_graph_diagnostics$source_is_row_order_default) %in% TRUE),
    fake_terminal_outcome_history_used = fake_terminal_outcome_history_used,
    multiple_L_blocks_collapsed = multiple_L_blocks_collapsed,
    visit_censoring_uses_downstream_nodes = visit_censoring_downstream,
    max_abs_treatment_factor = if (nrow(treatment_diag)) max(abs(treatment_diag$max_ratio), na.rm = TRUE) else NA_real_,
    max_abs_censoring_factor = if (nrow(censoring_diag)) max(abs(censoring_diag$max_ratio), na.rm = TRUE) else NA_real_,
    max_abs_mediator_density_ratio = if (nrow(mediator_diag)) max(abs(mediator_diag$max_ratio), na.rm = TRUE) else NA_real_,
    minimum_effective_sample_size = min(density_ratio_diagnostics$effective_sample_size, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  removed_ratio_arg_name <- .ltmle_exact_removed_density_ratio_arg()
  run_summary[[paste0(removed_ratio_arg_name, "_present")]] <- FALSE
  run_summary[[paste0(removed_ratio_arg_name, "_used")]] <- FALSE
  run_summary[[removed_ratio_arg_name]] <- "none"

  list(
    method = "ltmle_exact",
    estimator_class = .ltmle_exact_class(),
    learner = learner,
    means = as.list(means),
    means_initial = as.list(means_initial),
    means_targeted = as.list(means),
    subject_component_means = as.data.frame(subject_matrix),
    subject_component_initial_means = as.data.frame(subject_initial_matrix),
    diagnostics = list(
      estimator_class = .ltmle_exact_class(),
      component_registry = registry,
      world_specs = world_specs,
      factor_tasks = factor_tasks,
      ltmle_exact_component_law_registry = registry,
      ltmle_exact_factor_tasks = factor_tasks,
      score_diagnostics = score_diagnostics,
      targeting_empirical_equation_check = score_diagnostics,
      component_eif_matrix = as.data.frame(component_eif_matrix),
      component_eif_terms = component_eif_terms,
      component_eif_summary = component_eif_summary,
      ltmle_exact_component_eif_terms = component_eif_terms,
      ltmle_exact_component_eif_summary = component_eif_summary,
      density_ratio_diagnostics = density_ratio_diagnostics,
      truncation_diagnostics = truncation_diagnostics,
      clever_covariate_decomposition_diagnostics = clever_covariate_decomposition_diagnostics,
      separate_clever_covariate_identity_check = separate_clever_covariate_identity_check,
      second_M2_marginal_reference_check = second_M2_marginal_reference_check,
      mc_integration_diagnostics = mc_integration_diagnostics,
      law_task_integration_diagnostics = law_task_integration_diagnostics,
      law_integration_diagnostics = law_integration_diagnostics,
      branch_history_state_check = branch_history_state_check,
      branch_state_handoff_trace = branch_state_handoff_trace,
      separate_product_join_trace = separate_product_join_diagnostics$trace,
      separate_product_join_pre_collapse_trace = separate_product_join_diagnostics$pre_collapse_trace,
      source_eval_training_path_check = source_eval_training_path_check,
      cross_regimen_source_boundary_trace = cross_regimen_source_boundary_trace,
      source_pseudooutcome_independent_recursion_check = source_pseudooutcome_independent_recursion_check,
      targeted_continuation_storage_check = targeted_continuation_storage_check,
      dedicated_L_transition_L_fit_training_row_trace =
        dedicated_L_transition_L_fit_training_row_trace,
      dedicated_L_transition_L_fit_prediction_row_trace =
        dedicated_L_transition_L_fit_prediction_row_trace,
      dedicated_L_transition_L_fit_coefficient_contribution_check =
        dedicated_L_transition_L_fit_coefficient_contribution_check,
      dedicated_L_transition_prediction_weighting_collapse_check =
        dedicated_L_transition_prediction_weighting_collapse_check,
      dedicated_L_transition_mediator_path_weight_trace =
        dedicated_L_transition_mediator_path_weight_trace,
      dedicated_L_transition_training_response_oracle_alignment_check =
        dedicated_L_transition_training_response_oracle_alignment_check,
      allowed_local_source_prediction_map = .ltmle_exact_empty_allowed_local_source_prediction_map(),
      observed_vs_source_row_role_check = observed_vs_source_row_role_check,
      root_plugin_diagnostics = root_plugin_diagnostics,
      root_vs_targeted_branch_recursion_check = root_vs_targeted_branch_recursion_check,
      root_plugin_targeted_means = as.list(root_plugin_targeted_means),
      integrated_worldmean_recursion = integrated_worldmean_recursion,
      auxiliary_nuisance_gcomp_worldmeans = integrated_worldmean_recursion,
      task_graph_diagnostics = task_graph_diagnostics,
      ltmle_exact_required_checks = checks,
      ltmle_exact_run_summary = run_summary
    ),
    metadata = list(
      method = "ltmle_exact",
      estimator_class = .ltmle_exact_class(),
      learner = learner,
      Q_model = Q_model,
      terminal_plugin_type = "deterministic_root",
      terminal_plugin_mc_active = FALSE,
      effective_terminal_plugin_mc_n = 1L,
      ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
      ltmle_exact_law_integration_n = ltmle_exact_law_integration_n,
      worldmean_plugin_type = "targeted_deterministic_root",
      worldmean_plugin_mc_active = FALSE,
      effective_worldmean_plugin_mc_n = 1L,
      auxiliary_nuisance_gcomp_computed = !is.null(integrated_worldmean_recursion),
      auxiliary_nuisance_gcomp_mc_n = if (!is.null(integrated_worldmean_recursion)) {
        integrated_worldmean_recursion$B_mc
      } else {
        NA_integer_
      },
      diagnostics_level = diagnostics_level,
      y_bounds_mode = y_bounds_mode,
      probability_bounds = probability_bounds,
      probability_bounding_enabled = TRUE,
      probability_bounding_is_truncation = FALSE,
      truncation_enabled = isTRUE(truncation$enabled),
      truncation_policy = truncation$policy,
      truncation_target = truncation$target,
      truncation_quantile_lower = truncation$quantile_lower,
      truncation_quantile_upper = truncation$quantile_upper,
      fixed_bound_truncation_used = FALSE,
      density_ratio_factor_truncation_used = FALSE,
      mediator_density_engine = mediator_density_models$engine,
      mediator_density_crossfit_used = isTRUE(mediator_density_models$crossfit_used),
      treat_mech = treat_mech,
      censoring_mech = censoring_mech,
      censoring_vars = censoring_vars,
      crossfit = crossfit_info
      ,
      component_subset = component_subset
    )
  )
}

ltmle_exact_estimate_worldmeans <- function(dat = NULL, dat_wide = NULL, node_spec = NULL, T,
                                            reg_a, reg_as,
                                            learner = c("glm", "sl"),
                                            sl_library = c("SL.glm", "SL.mean"),
                                            Q_model = c("correct", "wrong"),
                                            seed = 202405L,
                                            probability_bounds = c(0.01, 0.99),
                                            truncation_enabled = TRUE,
                                            truncation_policy = "quantile",
                                            truncation_quantile_lower = 0.01,
                                            truncation_quantile_upper = 0.99,
                                            truncation_target = "clever_covariate_H",
                                            ltmle_exact_density_ratio_mc_n = 2000L,
                                            ltmle_exact_law_integration_n = 5L,
                                            y_bounds_mode = c("train_fold", "fixed"),
                                            y_bounds = NULL,
                                            treat_mech = c("observational", "baseline_rct"),
                                            p_rct = 0.5,
                                            censoring_vars = NULL,
                                            censoring_mech = c("none", "estimated"),
                                            score_tolerance = 1e-4,
                                            component_tolerance = 1e-5,
                                            scaled_z_tolerance = 3,
                                            V = 5L,
                                            fold_id = NULL,
                                            diagnostics_level = c("summary", "full"),
                                            compute_auxiliary_nuisance_gcomp = FALSE,
                                            auxiliary_nuisance_gcomp_mc_n = 2000L,
                                            force_epsilon_zero = FALSE,
                                            component_subset = NULL,
                                            verbose = FALSE,
                                            ...) {
  fit_ltmle_exact(
    dat = dat,
    dat_wide = dat_wide,
    node_spec = node_spec,
    T = T,
    reg_a = reg_a,
    reg_as = reg_as,
    learner = learner,
    sl_library = sl_library,
    Q_model = Q_model,
    seed = seed,
    probability_bounds = probability_bounds,
    truncation_enabled = truncation_enabled,
    truncation_policy = truncation_policy,
    truncation_quantile_lower = truncation_quantile_lower,
    truncation_quantile_upper = truncation_quantile_upper,
    truncation_target = truncation_target,
    ltmle_exact_density_ratio_mc_n = ltmle_exact_density_ratio_mc_n,
    ltmle_exact_law_integration_n = ltmle_exact_law_integration_n,
    y_bounds_mode = y_bounds_mode,
    y_bounds = y_bounds,
    treat_mech = treat_mech,
    p_rct = p_rct,
    censoring_vars = censoring_vars,
    censoring_mech = censoring_mech,
    score_tolerance = score_tolerance,
    component_tolerance = component_tolerance,
    scaled_z_tolerance = scaled_z_tolerance,
    V = V,
    fold_id = fold_id,
    diagnostics_level = diagnostics_level,
    compute_auxiliary_nuisance_gcomp = compute_auxiliary_nuisance_gcomp,
    auxiliary_nuisance_gcomp_mc_n = auxiliary_nuisance_gcomp_mc_n,
    force_epsilon_zero = force_epsilon_zero,
    component_subset = component_subset,
    verbose = verbose,
    ...
  )
}
