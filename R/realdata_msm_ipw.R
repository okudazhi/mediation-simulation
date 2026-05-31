################################################################################
# realdata_msm_ipw.R
#
# MSM + IPW estimator for the real-data layer.
#
# This implementation is written to track the revised manuscript more closely
# than the earlier draft package:
#   - treatment/regimen and censoring weights from pooled logistic models,
#   - mediator density ratios for the outcome MSM from ordered Gaussian models,
#   - a weighted marginal structural model for the final outcome,
#   - plug-in evaluation using mediator laws integrated from conditional models.
#
# The exact closed-form algebra in the manuscript appendix is approximated here
# through subject-specific plug-in prediction under the fitted linear working
# models.
################################################################################

.rd_build_formula <- function(response, rhs) {
  rhs <- unique(as.character(rhs))
  rhs <- rhs[nzchar(rhs)]
  if (!length(rhs)) return(stats::as.formula(paste(response, "~ 1")))
  stats::as.formula(paste(response, "~", paste(rhs, collapse = " + ")))
}


.rd_fallback_mean <- function(x, fallback = 0) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(as.numeric(fallback)[1L])
  mean(x)
}

.rd_normalize_weights <- function(weights, data, component = "model") {
  if (is.null(weights)) return(NULL)
  ww <- as.numeric(weights)
  if (length(ww) != nrow(data)) {
    .stop("Weight length mismatch in ", component, ": expected ",
          nrow(data), ", got ", length(ww))
  }
  ww[!is.finite(ww) | ww < 0] <- 0
  if (!length(ww) || all(ww == 0)) {
    .stop("All weights are zero in ", component, ".")
  }
  ww
}

.rd_fit_glm_or_stop <- function(formula,
                                data,
                                family = stats::binomial(),
                                component = "model",
                                weights = NULL,
                                engine = "glm",
                                sl_library = c("SL.glm", "SL.mean")) {
  fit <- tryCatch({
    ww <- .rd_normalize_weights(weights, data = data, component = component)
    if (!identical(engine, "glm")) {
      fit_nuisance_model(formula, data = data, family = family,
                         engine = engine, sl_library = sl_library,
                         component = component, weights = ww)
    } else {
      args <- list(formula = formula, data = data, family = family)
      if (!is.null(ww)) args$weights <- ww
      do.call(stats::glm, args)
    }
  }, error = function(e) e)
  if (inherits(fit, "error")) {
    .stop("Failed to fit ", component, ": ", fit$message)
  }
  fit
}

.rd_fit_lm_or_stop <- function(formula,
                               data,
                               weights = NULL,
                               component = "model",
                               engine = "glm",
                               sl_library = c("SL.glm", "SL.mean")) {
  fit <- tryCatch({
    ww <- .rd_normalize_weights(weights, data = data, component = component)
    if (!identical(engine, "glm")) {
      fit_nuisance_model(formula, data = data, family = stats::gaussian(),
                         engine = engine, sl_library = sl_library,
                         component = component, weights = ww)
    } else {
      args <- list(formula = formula, data = data)
      if (!is.null(ww)) args$weights <- ww
      do.call(stats::lm, args)
    }
  }, error = function(e) e)
  if (inherits(fit, "error")) {
    .stop("Failed to fit ", component, ": ", fit$message)
  }
  fit
}

.rd_at_risk_rows <- function(long) {
  is.finite(long$R_current) & (long$R_current == 1)
}

.rd_weight_diag_row <- function(name, x, raw = NULL, trunc_prob = c(0.01, 0.99)) {
  x <- as.numeric(x)
  raw_x <- if (is.null(raw)) x else as.numeric(raw)

  raw_fin <- raw_x[is.finite(raw_x)]
  if (length(raw_fin) > 0L) {
    bounds <- as.numeric(stats::quantile(raw_fin, probs = trunc_prob,
                                        na.rm = TRUE, names = FALSE, type = 7))
    prop_trunc_raw <- mean(raw_x < bounds[1L] | raw_x > bounds[2L], na.rm = TRUE)
  } else {
    bounds <- c(NA_real_, NA_real_)
    prop_trunc_raw <- NA_real_
  }

  q <- safe_quantile(x, c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1))
  ess <- ess_from_weights(x)

  data.frame(
    weight_type = name,
    n = sum(is.finite(x)),
    trunc_lo = bounds[1L],
    trunc_hi = bounds[2L],
    prop_trunc_raw = prop_trunc_raw,
    mean = mean(x, na.rm = TRUE),
    sd = stats::sd(x, na.rm = TRUE),
    min = q[1L],
    q01 = q[2L],
    q05 = q[3L],
    q25 = q[4L],
    median = q[5L],
    q75 = q[6L],
    q95 = q[7L],
    q99 = q[8L],
    max = q[9L],
    ess = ess[["ess"]],
    ess_frac = ess[["ess_frac"]],
    stringsAsFactors = FALSE
  )
}

.rd_visit_weight_diag <- function(long,
                                  cols = c("wA_step", "wC_step", "wM_y_step",
                                           "wA_y_cum", "wC_y_cum", "wM_y_cum")) {
  cols <- intersect(cols, names(long))
  if (!length(cols)) return(data.frame())

  out <- list()
  k <- 0L
  for (cc in cols) {
    for (tt in sort(unique(long$t))) {
      idx <- long$t == tt
      x <- as.numeric(long[[cc]][idx])
      q <- safe_quantile(x, c(0, 0.01, 0.05, 0.5, 0.95, 0.99, 1))
      ess <- ess_from_weights(x)
      k <- k + 1L
      out[[k]] <- data.frame(
        weight_type = cc,
        visit = tt,
        n = sum(is.finite(x)),
        mean = mean(x, na.rm = TRUE),
        sd = stats::sd(x, na.rm = TRUE),
        min = q[1L],
        q01 = q[2L],
        q05 = q[3L],
        median = q[4L],
        q95 = q[5L],
        q99 = q[6L],
        max = q[7L],
        ess = ess[["ess"]],
        ess_frac = ess[["ess_frac"]],
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, out)
}

.rd_weight_covars <- function(spec) {
  list(
    A_den = c("visit", spec$baseline_vars, "A_lag", "cumA_prev", "M1_lag", "M2_lag",
              paste0("L_", spec$L_names, "_lag")),
    A_num = c("visit", spec$baseline_vars, "A_lag", "cumA_prev"),
    C_den = c("visit", spec$baseline_vars, "A", "cumA_curr", "M1", "M2",
              paste0("L_", spec$L_names)),
    C_num = c("visit", spec$baseline_vars, "A", "cumA_curr"),
    M1_den = c("visit", spec$baseline_vars, "A", "cumA_curr", "M1_lag", "M2_lag",
               paste0("L_", spec$L_names, "_lag")),
    M1_num = c("visit", spec$baseline_vars, "A", "cumA_curr", "M1_lag", "M2_lag"),
    M2_den = c("visit", spec$baseline_vars, "A", "cumA_curr", "M1", "A_M1", "M1_lag", "M2_lag",
               paste0("L_", spec$L_names, "_lag")),
    M2_sep_num = c("visit", spec$baseline_vars, "A", "cumA_curr", "M1_lag", "M2_lag"),
    M2_joint_num = c("visit", spec$baseline_vars, "A", "cumA_curr", "M1", "A_M1", "M1_lag", "M2_lag")
  )
}

compute_realdata_weight_components <- function(prepared,
                                               control = default_esax_msm_cfg()) {
  spec <- prepared$spec
  long <- prepared$long
  trunc_prob <- control$trunc_prob %||% c(0.01, 0.99)
  eps <- control$eps %||% 1e-8
  engine <- control$nuisance_engine %||% "glm"
  sl_library <- control$sl_library %||% c("SL.glm", "SL.mean")
  covars <- .rd_weight_covars(spec)
  at_risk <- .rd_at_risk_rows(long)

  long$wA_step <- 1
  long$wC_step <- 1
  long$wM_y_step <- 1

  fitA_den <- NULL
  fitA_num <- NULL
  if (!(spec$treat_mech %in% c("known_randomized", "known_rct"))) {
    dat_A <- long[at_risk & !is.na(long$A), , drop = FALSE]
    fA_den <- .rd_build_formula("A", covars$A_den)
    fA_num <- .rd_build_formula("A", covars$A_num)
    fitA_den <- .rd_fit_glm_or_stop(fA_den, dat_A, component = "treatment denominator model", engine = engine, sl_library = sl_library)
    fitA_num <- .rd_fit_glm_or_stop(fA_num, dat_A, component = "treatment numerator model", engine = engine, sl_library = sl_library)
    pA_fallback <- mean(dat_A$A, na.rm = TRUE)
    p_den <- safe_glm_prob(fitA_den, dat_A, eps = eps, fallback = pA_fallback)
    p_num <- safe_glm_prob(fitA_num, dat_A, eps = eps, fallback = pA_fallback)
    pr_den <- ifelse(dat_A$A == 1, p_den, 1 - p_den)
    pr_num <- ifelse(dat_A$A == 1, p_num, 1 - p_num)
    long$wA_step[at_risk & !is.na(long$A)] <- pr_num / pmax(pr_den, eps)
  }

  fitC_den <- NULL
  fitC_num <- NULL
  dat_C <- long[at_risk & !is.na(long$R_next), , drop = FALSE]
  if (length(unique(dat_C$R_next[is.finite(dat_C$R_next)])) > 1L) {
    fC_den <- .rd_build_formula("R_next", covars$C_den)
    fC_num <- .rd_build_formula("R_next", covars$C_num)
    fitC_den <- .rd_fit_glm_or_stop(fC_den, dat_C, component = "censoring denominator model", engine = engine, sl_library = sl_library)
    fitC_num <- .rd_fit_glm_or_stop(fC_num, dat_C, component = "censoring numerator model", engine = engine, sl_library = sl_library)
    pC_fallback <- mean(dat_C$R_next, na.rm = TRUE)
    p_den <- safe_glm_prob(fitC_den, dat_C, eps = eps, fallback = pC_fallback)
    p_num <- safe_glm_prob(fitC_num, dat_C, eps = eps, fallback = pC_fallback)
    keep_idx <- dat_C$R_next == 1
    step_c <- rep(1, nrow(dat_C))
    step_c[keep_idx] <- p_num[keep_idx] / pmax(p_den[keep_idx], eps)
    long$wC_step[at_risk & !is.na(long$R_next)] <- step_c
  }

  dat_M <- long[at_risk & long$t >= 1L & !is.na(long$M1) & !is.na(long$M2), , drop = FALSE]
  fM1_den <- .rd_build_formula("M1", covars$M1_den)
  fM1_num <- .rd_build_formula("M1", covars$M1_num)
  fM2_den <- .rd_build_formula("M2", covars$M2_den)
  fM2_joint_num <- .rd_build_formula("M2", covars$M2_joint_num)

  fitM1_den <- .rd_fit_lm_or_stop(fM1_den, dat_M, component = "M1 density denominator model", engine = engine, sl_library = sl_library)
  fitM1_num <- .rd_fit_lm_or_stop(fM1_num, dat_M, component = "M1 density numerator model", engine = engine, sl_library = sl_library)
  fitM2_den <- .rd_fit_lm_or_stop(fM2_den, dat_M, component = "M2 density denominator model", engine = engine, sl_library = sl_library)
  fitM2_joint_num <- .rd_fit_lm_or_stop(fM2_joint_num, dat_M, component = "M2 common mediator-process numerator model", engine = engine, sl_library = sl_library)

  mu_M1_den <- safe_predict_lm(fitM1_den, dat_M, fallback = .rd_fallback_mean(dat_M$M1))
  mu_M1_num <- safe_predict_lm(fitM1_num, dat_M, fallback = .rd_fallback_mean(dat_M$M1))
  mu_M2_den <- safe_predict_lm(fitM2_den, dat_M, fallback = .rd_fallback_mean(dat_M$M2))
  mu_M2_joint_num <- safe_predict_lm(fitM2_joint_num, dat_M, fallback = .rd_fallback_mean(dat_M$M2))

  sig_M1_den <- model_sigma_safe(fitM1_den, fallback = 1e-6, ml = TRUE)
  sig_M1_num <- model_sigma_safe(fitM1_num, fallback = 1e-6, ml = TRUE)
  sig_M2_den <- model_sigma_safe(fitM2_den, fallback = 1e-6, ml = TRUE)
  sig_M2_joint_num <- model_sigma_safe(fitM2_joint_num, fallback = 1e-6, ml = TRUE)

  dens_den <- normal_density(dat_M$M1, mu_M1_den, sig_M1_den) *
    normal_density(dat_M$M2, mu_M2_den, sig_M2_den)
  dens_joint_num <- normal_density(dat_M$M1, mu_M1_num, sig_M1_num) *
    normal_density(dat_M$M2, mu_M2_joint_num, sig_M2_joint_num)

  idx_M <- at_risk & long$t >= 1L & !is.na(long$M1) & !is.na(long$M2)
  long$wM_y_step[idx_M] <- dens_joint_num / pmax(dens_den, eps)

  long$wA_y_cum <- cumprod_by_id(long$wA_step, long$ID)
  long$wC_y_cum <- cumprod_by_id(long$wC_step, long$ID)
  long$wM_y_cum <- cumprod_by_id(long$wM_y_step, long$ID)

  long$wC_m_cum <- lag_by_id(long$wC_y_cum, long$ID, default = 1)
  long$wA_m_cum <- long$wA_y_cum
  long$w_m_visit <- long$wA_m_cum * long$wC_m_cum

  last_t <- max(long$t)
  end <- long[long$t == last_t, c("ID", "wA_y_cum", "wC_y_cum", "wM_y_cum"), drop = FALSE]
  end$w_nat_raw <- end$wA_y_cum * end$wC_y_cum
  end$w_int_raw <- end$w_nat_raw * end$wM_y_cum

  end$w_nat <- truncate_vec(end$w_nat_raw, trunc_prob[1], trunc_prob[2])
  end$w_int <- truncate_vec(end$w_int_raw, trunc_prob[1], trunc_prob[2])
  end$sw_final <- end$w_int

  diag_weights <- rbind(
    .rd_weight_diag_row("w_nat_raw", end$w_nat_raw, raw = end$w_nat_raw, trunc_prob = trunc_prob),
    .rd_weight_diag_row("w_nat_trunc", end$w_nat, raw = end$w_nat_raw, trunc_prob = trunc_prob),
    .rd_weight_diag_row("w_int_raw", end$w_int_raw, raw = end$w_int_raw, trunc_prob = trunc_prob),
    .rd_weight_diag_row("w_int_trunc", end$w_int, raw = end$w_int_raw, trunc_prob = trunc_prob)
  )

  diag_visit_weights <- .rd_visit_weight_diag(long)

  list(
    long = long,
    end_weights = end,
    diagnostics = diag_weights,
    diagnostics_visit = diag_visit_weights,
    nuisance = list(
      fitA_den = fitA_den,
      fitA_num = fitA_num,
      fitC_den = fitC_den,
      fitC_num = fitC_num,
      fitM1_den = fitM1_den,
      fitM1_num = fitM1_num,
      fitM2_den = fitM2_den,
      fitM2_joint_num = fitM2_joint_num,
      sigma = list(
        M1_den = sig_M1_den,
        M1_num = sig_M1_num,
        M2_den = sig_M2_den,
        M2_joint_num = sig_M2_joint_num
      )
    )
  )
}

.rd_add_outcome_design_cols <- function(wide,
                                        T,
                                        include_AM = TRUE,
                                        include_MM = FALSE) {
  out <- wide
  for (tt in 0:(T - 1L)) {
    Acol <- paste0("A", tt)
    M1col <- paste0("M1_", tt)
    M2col <- paste0("M2_", tt)
    if (isTRUE(include_AM)) {
      out[[paste0("AM1_", tt)]] <- out[[Acol]] * out[[M1col]]
      out[[paste0("AM2_", tt)]] <- out[[Acol]] * out[[M2col]]
    }
    if (isTRUE(include_MM)) {
      out[[paste0("MM_", tt)]] <- out[[M1col]] * out[[M2col]]
    }
  }
  out
}

fit_realdata_outcome_msms <- function(prepared,
                                      weight_obj,
                                      control = default_esax_msm_cfg()) {
  spec <- prepared$spec
  T <- spec$T
  include_AM <- isTRUE(control$include_A_M_interactions %||% TRUE)
  include_MM <- isTRUE(control$include_M1_M2_interactions %||% FALSE)
  engine <- control$nuisance_engine %||% "glm"
  sl_library <- control$sl_library %||% c("SL.glm", "SL.mean")

  wide <- .rd_add_outcome_design_cols(prepared$wide, T = T,
                                      include_AM = include_AM,
                                      include_MM = include_MM)
  ew <- weight_obj$end_weights
  idx <- match(wide$ID, ew$ID)
  w_nat <- ew$w_nat[idx]
  w_int <- ew$w_int[idx]

  A_terms <- paste0("A", 0:(T - 1L))
  if (spec$treat_mech %in% c("known_randomized", "known_rct", "baseline_rct")) {
    A_terms <- "A0"
  }
  post_idx <- if (T >= 2L) 1:(T - 1L) else integer(0)

  rhs_nat <- unique(c(spec$baseline_vars, A_terms))
  fN <- .rd_build_formula("Y_final", rhs_nat)

  rhs_int <- c(
    spec$baseline_vars,
    A_terms,
    paste0("M1_", post_idx),
    paste0("M2_", post_idx)
  )
  if (isTRUE(include_AM)) rhs_int <- c(rhs_int, paste0("AM1_", post_idx), paste0("AM2_", post_idx))
  if (isTRUE(include_MM)) rhs_int <- c(rhs_int, paste0("MM_", post_idx))
  fY <- .rd_build_formula("Y_final", rhs_int)

  keep_y <- is.finite(wide$Y_final) & (wide$R_final == 1)

  fit_nat <- .rd_fit_lm_or_stop(
    fN, wide[keep_y, , drop = FALSE],
    weights = w_nat[keep_y],
    component = "natural treatment-only outcome MSM",
    engine = engine,
    sl_library = sl_library
  )
  fit_int <- .rd_fit_lm_or_stop(
    fY, wide[keep_y, , drop = FALSE],
    weights = w_int[keep_y],
    component = "controlled interventional outcome MSM",
    engine = engine,
    sl_library = sl_library
  )

  list(
    formula_nat = deparse(fN),
    formula_int = deparse(fY),
    fit_nat = fit_nat,
    fit_int = fit_int,
    include_AM = include_AM,
    include_MM = include_MM
  )
}

fit_realdata_mediator_msms <- function(prepared,
                                       weight_obj,
                                       control = default_esax_msm_cfg()) {
  spec <- prepared$spec
  long <- prepared$long
  at_risk <- .rd_at_risk_rows(long)
  engine <- control$nuisance_engine %||% "glm"
  sl_library <- control$sl_library %||% c("SL.glm", "SL.mean")
  idx <- at_risk & long$t >= 1L & !is.na(long$M1) & !is.na(long$M2)
  dat <- long[idx, , drop = FALSE]
  w_m <- long$w_m_visit[idx]

  fM1 <- .rd_build_formula("M1", c("visit", spec$baseline_vars, "A", "A_lag", "cumA_prev", "M1_lag"))
  fM2_sep <- .rd_build_formula("M2", c("visit", spec$baseline_vars, "A", "A_lag", "cumA_prev", "M2_lag"))
  fM2_joint <- .rd_build_formula("M2", c("visit", spec$baseline_vars, "A", "A_lag", "cumA_prev", "M2_lag", "M1", "A_M1"))

  fitM1 <- .rd_fit_lm_or_stop(fM1, dat, weights = w_m, component = "mediator MSM for M1", engine = engine, sl_library = sl_library)
  fitM2_sep <- .rd_fit_lm_or_stop(fM2_sep, dat, weights = w_m, component = "mediator MSM for M2 (separate)", engine = engine, sl_library = sl_library)
  fitM2_joint <- .rd_fit_lm_or_stop(fM2_joint, dat, weights = w_m, component = "mediator MSM for M2 (joint)", engine = engine, sl_library = sl_library)

  list(
    fitM1 = fitM1,
    fitM2_sep = fitM2_sep,
    fitM2_joint = fitM2_joint,
    formula_M1 = deparse(fM1),
    formula_M2_sep = deparse(fM2_sep),
    formula_M2_joint = deparse(fM2_joint)
  )
}

.rd_baseline_frame <- function(prepared) {
  req <- c("ID", prepared$spec$baseline_vars, "M1_0", "M2_0")
  miss <- setdiff(req, names(prepared$wide))
  if (length(miss)) {
    .stop("prepared wide data are missing baseline mediator columns required by MSM evaluation: ",
          paste(miss, collapse = ", "))
  }
  prepared$wide[, req, drop = FALSE]
}

.rd_predict_path_m1 <- function(prepared,
                                fitM1,
                                regimen) {
  spec <- prepared$spec
  T <- spec$T
  reg <- normalize_regimen(regimen, T, "regimen")
  base <- .rd_baseline_frame(prepared)
  n <- nrow(base)

  out <- matrix(NA_real_, nrow = n, ncol = T)
  out[, 1L] <- as.numeric(base$M1_0)
  M1_lag <- as.numeric(base$M1_0)

  if (T >= 2L) {
    for (tt in 1:(T - 1L)) {
      nd <- data.frame(
        visit = factor(tt, levels = 0:(T - 1L)),
        A = rep(reg[tt + 1L], n),
        A_lag = rep(reg[tt], n),
        cumA_prev = rep(sum(reg[1:tt]), n),
        M1_lag = M1_lag,
        stringsAsFactors = FALSE
      )
      for (nm in spec$baseline_vars) nd[[nm]] <- base[[nm]]
      M1_t <- safe_predict_lm(fitM1, nd, fallback = .rd_fallback_mean(M1_lag))
      out[, tt + 1L] <- M1_t
      M1_lag <- M1_t
    }
  }

  out
}

.rd_predict_path_m2 <- function(prepared,
                                fitM2,
                                regimen,
                                M1_path = NULL,
                                joint = FALSE) {
  spec <- prepared$spec
  T <- spec$T
  reg <- normalize_regimen(regimen, T, "regimen")
  base <- .rd_baseline_frame(prepared)
  n <- nrow(base)

  out <- matrix(NA_real_, nrow = n, ncol = T)
  out[, 1L] <- as.numeric(base$M2_0)
  M2_lag <- as.numeric(base$M2_0)

  if (T >= 2L) {
    for (tt in 1:(T - 1L)) {
      nd <- data.frame(
        visit = factor(tt, levels = 0:(T - 1L)),
        A = rep(reg[tt + 1L], n),
        A_lag = rep(reg[tt], n),
        cumA_prev = rep(sum(reg[1:tt]), n),
        M2_lag = M2_lag,
        stringsAsFactors = FALSE
      )
      for (nm in spec$baseline_vars) nd[[nm]] <- base[[nm]]
      if (isTRUE(joint)) {
        if (is.null(M1_path)) .stop("M1_path is required when joint=TRUE.")
        nd$M1 <- M1_path[, tt + 1L]
        nd$A_M1 <- nd$A * nd$M1
      }
      M2_t <- safe_predict_lm(fitM2, nd, fallback = .rd_fallback_mean(M2_lag))
      out[, tt + 1L] <- M2_t
      M2_lag <- M2_t
    }
  }

  out
}

.rd_build_outcome_newdata <- function(prepared,
                                      outer_regimen,
                                      M1_path,
                                      M2_path,
                                      outcome_fits) {
  spec <- prepared$spec
  T <- spec$T
  reg <- normalize_regimen(outer_regimen, T, "outer_regimen")
  base <- .rd_baseline_frame(prepared)

  nd <- data.frame(
    ID = base$ID,
    stringsAsFactors = FALSE
  )
  for (nm in spec$baseline_vars) nd[[nm]] <- base[[nm]]
  for (tt in 0:(T - 1L)) {
    nd[[paste0("A", tt)]] <- rep(reg[tt + 1L], nrow(base))
    nd[[paste0("M1_", tt)]] <- M1_path[, tt + 1L]
    nd[[paste0("M2_", tt)]] <- M2_path[, tt + 1L]
  }
  .rd_add_outcome_design_cols(
    nd,
    T = T,
    include_AM = isTRUE(outcome_fits$include_AM),
    include_MM = isTRUE(outcome_fits$include_MM)
  )
}

.rd_fit_conditional_models_for_msm_component_means <- function(prepared,
                                                               control = default_esax_msm_cfg(),
                                                               seed = NULL) {
  if (!is.null(seed)) set.seed(as.integer(seed)[1L])
  input <- .rd_ltmle_exact_input_from_prepared(prepared)
  T <- input$node_spec$T
  learner <- if (identical(control$nuisance_engine %||% "glm", "superlearner")) "sl" else "glm"
  model_set <- control$model_set %||% "main"
  q_model_for_ltmle_engine <- control$ltmle_exact_Q_model %||% "correct"

  dat_wide <- .ltmle_exact_canonicalize_node_spec(input$data, input$node_spec, T)
  node_spec <- attr(dat_wide, "ltmle_exact_node_spec")
  long <- .ltmle_exact_node_spec_to_long(dat_wide, T, node_spec)

  for (nm in as.character(node_spec$baseline_vars %||% character(0))) {
    if (nm %in% names(dat_wide) && !nm %in% names(long)) {
      long[[nm]] <- as.numeric(dat_wide[[nm]])[long$id]
    }
  }

  models <- .ltmle_exact_fit_node_models(
    long = long,
    T = T,
    learner = learner,
    sl_library = control$sl_library %||% c("SL.glm", "SL.mean"),
    Q_model = q_model_for_ltmle_engine,
    node_spec = node_spec
  )

  list(
    dat_wide = dat_wide,
    node_spec = node_spec,
    models = models,
    model_set = model_set,
    component_mean_model_engine = "ltmle_exact_node_models_for_msm_component_mean_integration"
  )
}

.rd_extract_mediator_path_from_ltmle_histories <- function(histories, T, mediator = c("M1", "M2")) {
  mediator <- match.arg(mediator)
  first_row <- histories$outcome_process[[1L]]$row
  n <- nrow(first_row)
  out <- matrix(NA_real_, nrow = n, ncol = T)

  base_col <- paste0(mediator, "_0")
  if (base_col %in% names(first_row)) {
    out[, 1L] <- as.numeric(first_row[[base_col]])
  } else if (mediator %in% names(first_row)) {
    out[, 1L] <- as.numeric(first_row[[mediator]])
  } else {
    .stop("Cannot extract baseline ", mediator, " from generated history.")
  }

  if (T >= 2L) {
    for (tt in 2:T) {
      row_tt <- histories$outcome_process[[tt]]$row
      if (!mediator %in% names(row_tt)) {
        .stop("Generated history is missing ", mediator, " at time index ", tt - 1L, ".")
      }
      out[, tt] <- as.numeric(row_tt[[mediator]])
    }
  }
  out
}

.rd_world_spec_for_msm_target <- function(spec, target, component_name, T) {
  world_specs <- .ltmle_exact_world_spec(spec$reg_a, spec$reg_as, T)
  sp <- world_specs[world_specs$component == component_name, , drop = FALSE]
  if (nrow(sp) == 1L) return(sp)

  outer <- normalize_regimen(target$outer_regimen, T, "target$outer_regimen")
  if (identical(target$type, "natural")) {
    world_type <- "natural"
    m1 <- outer
    m2 <- outer
  } else if (identical(target$type, "joint")) {
    world_type <- "joint"
    m1 <- normalize_regimen(target$med_regimen, T, "target$med_regimen")
    m2 <- m1
  } else if (target$type %in% c("separate", "product")) {
    world_type <- "separate"
    m1 <- normalize_regimen(target$med_regimen1, T, "target$med_regimen1")
    m2 <- normalize_regimen(target$med_regimen2, T, "target$med_regimen2")
  } else {
    .stop("Unsupported target type: ", target$type)
  }

  data.frame(
    component_id = NA_integer_,
    component = component_name,
    world_type = world_type,
    t = seq_len(T),
    outer_A = outer,
    m1_A = m1,
    m2_A = m2,
    stringsAsFactors = FALSE
  )
}

.evaluate_realdata_target_mean_msm_from_conditional_models <- function(prepared,
                                                                      outcome_fits,
                                                                      component_mean_models,
                                                                      target,
                                                                      mc_n,
                                                                      seed = NULL) {
  spec <- prepared$spec
  T <- component_mean_models$node_spec$T
  dat_wide <- component_mean_models$dat_wide
  node_spec <- component_mean_models$node_spec
  models <- component_mean_models$models

  if (!is.null(seed)) set.seed(as.integer(seed)[1L])
  mc_n <- as.integer(mc_n)
  if (!is.finite(mc_n) || mc_n < 1L) .stop("mc_n must be a positive integer.")

  component_name <- target$component %||% target$name %||% NULL
  if (is.null(component_name)) {
    .stop("Each real-data target must provide a component or name field.")
  }
  sp <- .rd_world_spec_for_msm_target(spec, target, component_name, T)

  fit <- if (identical(target$type, "natural")) {
    outcome_fits$fit_nat
  } else if (target$type %in% c("separate", "product")) {
    outcome_fits$fit_int
  } else if (identical(target$type, "joint")) {
    outcome_fits$fit_int
  } else {
    .stop("Unsupported target type: ", target$type)
  }

  pred_sum <- 0
  n_pred <- 0L

  for (b in seq_len(mc_n)) {
    histories <- .ltmle_exact_draw_component_histories(
      dat_wide = dat_wide,
      T = T,
      spec = sp,
      models = models,
      mc_n = 1L,
      seed = if (is.null(seed)) NULL else as.integer(seed)[1L] + b * 1009L,
      node_spec = node_spec
    )

    M1_path <- .rd_extract_mediator_path_from_ltmle_histories(histories, T, mediator = "M1")
    M2_path <- .rd_extract_mediator_path_from_ltmle_histories(histories, T, mediator = "M2")

    nd <- .rd_build_outcome_newdata(
      prepared = prepared,
      outer_regimen = as.numeric(sp$outer_A),
      M1_path = M1_path,
      M2_path = M2_path,
      outcome_fits = outcome_fits
    )

    pred <- safe_predict_lm(fit, nd, fallback = NA_real_)
    pred_sum <- pred_sum + sum(pred, na.rm = TRUE)
    n_pred <- n_pred + sum(is.finite(pred))
  }

  pred_sum / n_pred
}

.compute_realdata_effects <- function(registry_means, spec) {
  fixed_keys <- names(make_fixed_horizon_registry(spec$reg_a, spec$reg_as, spec$T))
  fixed_means <- as.list(registry_means[fixed_keys])

  fixed_main <- compute_fixed_horizon_main_effects_from_means(fixed_means)
  joint <- compute_joint_draw_effects_from_means(fixed_means)

  time_resolved <- data.frame()
  time_registry <- make_time_resolved_sep_registry(spec$reg_a, spec$reg_as, cuts = spec$cuts)
  if (length(time_registry)) {
    tmp <- compute_time_resolved_sep_effects_from_means(as.list(registry_means[names(time_registry)]))
    time_resolved <- tmp$time_resolved
  }

  list(
    fixed = fixed_main,
    joint = joint,
    time_resolved = time_resolved
  )
}

estimate_realdata_msm_ipw <- function(prepared,
                                      control = default_esax_msm_cfg(),
                                      seed = NULL) {
  weight_obj <- compute_realdata_weight_components(prepared, control = control)
  outcome_fits <- fit_realdata_outcome_msms(prepared, weight_obj, control = control)
  component_mean_models <- .rd_fit_conditional_models_for_msm_component_means(
    prepared = prepared,
    control = control,
    seed = seed
  )

  fixed_registry <- make_fixed_horizon_registry(prepared$spec$reg_a, prepared$spec$reg_as, prepared$spec$T)
  time_registry <- make_time_resolved_sep_registry(prepared$spec$reg_a, prepared$spec$reg_as, cuts = prepared$spec$cuts)
  registry <- c(fixed_registry, time_registry)

  mc_n <- as.integer(control$mediator_mc_n %||% control$mc_B_main %||% 2000L)
  means <- vapply(seq_along(registry), function(ii) {
    tgt <- registry[[ii]]
    if (is.null(tgt$component) && is.null(tgt$name)) {
      tgt$component <- names(registry)[ii]
    }
    .evaluate_realdata_target_mean_msm_from_conditional_models(
      prepared = prepared,
      outcome_fits = outcome_fits,
      component_mean_models = component_mean_models,
      target = tgt,
      mc_n = mc_n,
      seed = if (is.null(seed)) NULL else as.integer(seed)[1L] + ii * 1009L
    )
  }, numeric(1))
  names(means) <- names(registry)

  effects <- .compute_realdata_effects(means, prepared$spec)

  endw <- weight_obj$end_weights
  idx <- match(prepared$wide$ID, endw$ID)
  subject_weights <- data.frame(
    ID = prepared$wide$ID,
    sw_final = endw$sw_final[idx],
    w_nat_raw = endw$w_nat_raw[idx],
    w_nat_trunc = endw$w_nat[idx],
    w_int_raw = endw$w_int_raw[idx],
    w_int_trunc = endw$w_int[idx],
    stringsAsFactors = FALSE
  )

  out <- list(
    estimator = "msm_ipw",
    analysis_name = prepared$spec$analysis_name,
    spec = prepared$spec,
    means = as.list(means[names(fixed_registry)]),
    registry_means = as.list(means),
    effects = effects,
    weight_models = list(
      subject_weights = subject_weights,
      long = weight_obj$long
    ),
    diagnostics = list(
      weight_summary = weight_obj$diagnostics,
      visit_weight_summary = weight_obj$diagnostics_visit,
      outcome_formula_nat = outcome_fits$formula_nat,
      outcome_formula_int = outcome_fits$formula_int,
      mediator_intervention_distribution_source =
        "conditional model integration for interventional mediator distributions",
      model_set = component_mean_models$model_set,
      component_mean_model_engine = component_mean_models$component_mean_model_engine
    ),
    model_fits = list(
      outcome = outcome_fits,
      component_mean_conditional_models = component_mean_models$models
    ),
    mc_diagnostics = NULL,
    targeting_diagnostics = NULL,
    control = control,
    registry = list(fixed = fixed_registry, time = time_registry)
  )
  class(out) <- c("realdata_msm_ipw_result", "list")
  out
}

# Backward-compatible name used by the original runner.
realdata_msm_ipw_estimate <- estimate_realdata_msm_ipw
