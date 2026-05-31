################################################################################
# estimator_msm_ipw.R
#
# MSMs + IPW for interventional effects (Yamamuro marginal & Tai joint)
#
# This file implements the *coefficient-based* MSM approach:
#
#   (A) Fit (stabilized) weights for the outcome MSM:
#         SW_Y(t) = W_A^Y(t) * W_M(t)
#
#       - W_A^Y(t): treatment stabilized weight where the numerator excludes
#         time-varying confounders and the denominator includes them.
#
#       - W_M(t): mediator stabilized weight built from Gaussian density ratios.
#
#         The mediator-process numerator is a stabilization choice for fitting
#         one controlled-outcome MSM; target distinctions enter at plug-in.
#
#   (B) For both Q_model = "correct" and Q_model = "wrong", evaluate
#       interventional mediator distributions by numerical integration under the
#       corresponding conditional node models.
#
#   (C) Fit the outcome MSM (linear regression) and use coefficient plug-in
#       to build the 9 world means in the package's canonical order.
#
# IMPORTANT:
#   - This estimator targets the terminal outcome Y_T.
#   - Baseline covariates V are plugged in at their sample means. Because all
#     MSMs are linear in V, this is equivalent to averaging over the empirical
#     distribution of V.
#
# Dependencies:
#   - core_utils.R (normalize_regimen, wide_to_long, truncate_vec, etc.)
#   - core_targets.R (compute_estimands_from_means)
################################################################################

# NOTE: This repository is used both as a small package and as standalone files.
# We keep a lightweight dependency check helper. Use a file-specific name to
# avoid accidental overwrites when multiple scripts are sourced.

.msm_source_if_needed <- function() {
  need_stop <- !exists(".stop", mode = "function")
  need_estimands <- !exists("compute_estimands_from_means", mode = "function")
  if (!need_stop && !need_estimands) return(invisible(TRUE))

  candidates <- c(
    file.path(getwd(), "core_utils.R"),
    file.path(getwd(), "core_targets.R"),
    file.path(getwd(), "R", "core_utils.R"),
    file.path(getwd(), "R", "core_targets.R")
  )
  candidates <- unique(normalizePath(candidates[file.exists(candidates)], winslash = "/", mustWork = FALSE))

  if (length(candidates)) {
    for (fp in candidates) source(fp, chdir = TRUE)
  }

  if (!exists(".stop", mode = "function") || !exists("compute_estimands_from_means", mode = "function")) {
    stop(
      paste0(
        "core_utils.R and core_targets.R must be sourced before estimator_msm_ipw.R. ",
        "Either run run_simulation.R from the package root or source ",
        "R/core_utils.R and R/core_targets.R explicitly."
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

.clamp01 <- function(p, eps = 1e-9) {
  p <- pmax(eps, pmin(1 - eps, p))
  p
}

# ML sigma for Gaussian lm: sqrt(RSS / n)
.ml_sigma <- function(lm_fit) {
  r <- stats::residuals(lm_fit)
  n <- length(r)
  if (n <= 0) return(NA_real_)
  s2 <- sum(r^2) / n
  sqrt(max(s2, 0))
}

.safe_sigma <- function(sig, fallback = 1e-6) {
  if (!is.finite(sig) || sig <= 0) return(fallback)
  sig
}

.dens_norm <- function(x, mu, sd) {
  stats::dnorm(x, mean = mu, sd = sd)
}

.msm_gcomp_source_if_needed <- function() {
  need <- !exists("fit_nuisance_models_gcomp", mode = "function") ||
    !exists("simulate_full_paths_once", mode = "function")
  if (!need) return(invisible(TRUE))

  candidates <- c(
    file.path(getwd(), "estimator_gcomp.R"),
    file.path(getwd(), "R", "estimator_gcomp.R")
  )
  candidates <- unique(normalizePath(candidates[file.exists(candidates)], winslash = "/", mustWork = FALSE))
  for (fp in candidates) source(fp, chdir = TRUE)

  if (!exists("fit_nuisance_models_gcomp", mode = "function") ||
      !exists("simulate_full_paths_once", mode = "function")) {
    .stop("estimator_gcomp.R must be sourced before MSM mediator-law evaluation.")
  }
  invisible(TRUE)
}

.msm_as_positive_integer <- function(x, default = 5000L, name = "B_mc") {
  if (is.null(x)) x <- default
  x <- as.integer(x)[1L]
  if (!is.finite(x) || x < 1L) .stop(name, " must be a positive integer.")
  x
}

.msm_normalize_truncation <- function(trunc = NULL,
                                      truncation_enabled = TRUE,
                                      truncation_policy = "quantile",
                                      truncation_quantile_lower = 0.01,
                                      truncation_quantile_upper = 0.99,
                                      truncation_target = "final_cumulative_weight") {
  if (!is.null(trunc)) {
    trunc <- as.numeric(trunc)
    if (length(trunc) != 2L || any(!is.finite(trunc))) {
      .stop("trunc must be a numeric vector of length 2 when supplied.")
    }
    truncation_quantile_lower <- trunc[1L]
    truncation_quantile_upper <- trunc[2L]
  }
  truncation_policy <- as.character(truncation_policy %||% "quantile")[1L]
  if (!truncation_policy %in% c("quantile", "none")) {
    .stop("MSM/IPW truncation_policy must be either 'quantile' or 'none'.")
  }
  truncation_enabled <- isTRUE(truncation_enabled) && !identical(truncation_policy, "none")
  if (!isTRUE(truncation_enabled)) truncation_policy <- "none"
  truncation_target <- as.character(truncation_target %||% "final_cumulative_weight")[1L]
  if (!identical(truncation_target, "final_cumulative_weight")) {
    .stop("MSM/IPW truncation_target must be 'final_cumulative_weight'.")
  }
  q_lower <- as.numeric(truncation_quantile_lower)[1L]
  q_upper <- as.numeric(truncation_quantile_upper)[1L]
  if (identical(truncation_policy, "quantile")) {
    if (!is.finite(q_lower) || !is.finite(q_upper) ||
        q_lower < 0 || q_upper > 1 || q_lower >= q_upper) {
      .stop("MSM/IPW truncation quantiles must satisfy 0 <= lower < upper <= 1.")
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

.msm_ess <- function(w) {
  w <- as.numeric(w)
  w <- w[is.finite(w)]
  if (!length(w)) return(NA_real_)
  sw <- sum(w)
  sw2 <- sum(w * w)
  if (!is.finite(sw) || !is.finite(sw2) || sw2 <= 0) return(NA_real_)
  (sw * sw) / sw2
}

.msm_quantile_value <- function(x, p) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x) || !is.finite(p)) return(NA_real_)
  as.numeric(stats::quantile(x, probs = p, na.rm = TRUE, names = FALSE, type = 7))
}

.msm_truncate_final_weight <- function(x, truncation) {
  x <- as.numeric(x)
  if (isTRUE(truncation$enabled) && identical(truncation$policy, "quantile")) {
    qs <- stats::quantile(
      x,
      probs = c(truncation$quantile_lower, truncation$quantile_upper),
      na.rm = TRUE,
      names = FALSE,
      type = 7
    )
    return(list(
      value = pmin(pmax(x, qs[1L]), qs[2L]),
      lower = as.numeric(qs[1L]),
      upper = as.numeric(qs[2L])
    ))
  }
  list(value = x, lower = NA_real_, upper = NA_real_)
}

.msm_truncation_summary_row <- function(raw,
                                        truncated,
                                        truncation,
                                        estimator_variant,
                                        component,
                                        task_id,
                                        q_lower_value = NA_real_,
                                        q_upper_value = NA_real_) {
  raw <- as.numeric(raw)
  truncated <- as.numeric(truncated)
  raw_finite <- raw[is.finite(raw)]
  trunc_finite <- truncated[is.finite(truncated)]
  if (!length(raw_finite)) raw_finite <- NA_real_
  if (!length(trunc_finite)) trunc_finite <- NA_real_
  n_values <- sum(is.finite(raw))
  n_lower <- if (isTRUE(truncation$enabled) && is.finite(q_lower_value)) {
    sum(is.finite(raw) & raw < q_lower_value)
  } else {
    0L
  }
  n_upper <- if (isTRUE(truncation$enabled) && is.finite(q_upper_value)) {
    sum(is.finite(raw) & raw > q_upper_value)
  } else {
    0L
  }
  data.frame(
    estimator = "msm_ipw",
    estimator_variant = estimator_variant,
    component = component,
    task_id = task_id,
    t = NA_integer_,
    node = "final_weight",
    process_type = "outcome_msm_weight",
    truncation_enabled = isTRUE(truncation$enabled),
    truncation_policy = truncation$policy,
    truncation_target = truncation$target,
    truncation_rule = if (isTRUE(truncation$enabled)) "sample_quantile" else "none",
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
    p01_raw_value = .msm_quantile_value(raw, 0.01),
    median_raw_value = stats::median(raw_finite, na.rm = TRUE),
    p99_raw_value = .msm_quantile_value(raw, 0.99),
    max_raw_value = max(raw_finite, na.rm = TRUE),
    min_truncated_value = min(trunc_finite, na.rm = TRUE),
    max_truncated_value = max(trunc_finite, na.rm = TRUE),
    mean_raw_value = mean(raw_finite, na.rm = TRUE),
    sd_raw_value = stats::sd(raw_finite, na.rm = TRUE),
    mean_truncated_value = mean(trunc_finite, na.rm = TRUE),
    sd_truncated_value = stats::sd(trunc_finite, na.rm = TRUE),
    ess_raw = .msm_ess(raw),
    ess_truncated = .msm_ess(truncated),
    passed = TRUE,
    failure_class = "no_failure",
    stringsAsFactors = FALSE
  )
}

.msm_predict_lm_numeric <- function(fit, newdata, fallback = 0) {
  n <- nrow(newdata)
  out <- tryCatch(as.numeric(stats::predict(fit, newdata = newdata)),
                  error = function(e) rep(NA_real_, n))
  if (length(out) != n) out <- rep(NA_real_, n)
  bad <- !is.finite(out)
  if (any(bad)) out[bad] <- as.numeric(fallback)[1L]
  out
}

.msm_fit_conditional_models_for_interventional_mediator_distributions <- function(dat_wide, T, Q_model = c("correct", "wrong")) {
  Q_model <- match.arg(Q_model)
  .msm_gcomp_source_if_needed()
  fit_nuisance_models_gcomp(wide_to_long(dat_wide, T), T, Q_model = Q_model)
}

.msm_gcomp_fits_for_mediator_law <- .msm_fit_conditional_models_for_interventional_mediator_distributions

.msm_mediator_density_numerators_from_conditional_models <- function(dat_wide, T, fits, t,
                                                                     B_mc = 5000L,
                                                                     seed = NULL) {
  B_mc <- .msm_as_positive_integer(B_mc, name = "mediator_density_mc_n")
  if (!is.null(seed)) set.seed(as.integer(seed)[1L])

  varmap <- default_varmap(T)
  n <- nrow(dat_wide)
  obs_m1 <- as.numeric(dat_wide[[varmap$M1[t]]])
  obs_m2 <- as.numeric(dat_wide[[varmap$M2[t]]])

  d_m1_sum <- numeric(n)
  d_m2_marg_sum <- numeric(n)
  d_joint_sum <- numeric(n)

  for (b in seq_len(B_mc)) {
    L_lag <- rep(0, n)
    Y_lag <- as.numeric(dat_wide$Y0)
    M1_lag <- as.numeric(dat_wide$M1_0)
    M2_lag <- as.numeric(dat_wide$M2_0)

    if (t > 1L) {
      for (s in seq_len(t - 1L)) {
        A_s <- as.numeric(dat_wide[[varmap$A[s]]])
        M1_s <- as.numeric(dat_wide[[varmap$M1[s]]])
        M2_s <- as.numeric(dat_wide[[varmap$M2[s]]])

        ndL <- data.frame(
          A = A_s, W1 = dat_wide$W1, W2 = dat_wide$W2,
          L_lag = L_lag, Y_lag = Y_lag, M1 = M1_s, M2 = M2_s
        )
        muL <- .msm_predict_lm_numeric(fits$L[[s]], ndL, fallback = mean(ndL$M1 + ndL$M2, na.rm = TRUE))
        L_s <- muL + stats::rnorm(n, mean = 0, sd = fits$sigma_L[s])

        ndY <- data.frame(
          W1 = dat_wide$W1, W2 = dat_wide$W2,
          Y0 = dat_wide$Y0, Y_lag = Y_lag, L_lag = L_lag,
          M1_lag = M1_lag, M2_lag = M2_lag,
          A = A_s, M1 = M1_s, M2 = M2_s, L = L_s
        )
        muY <- .msm_predict_lm_numeric(fits$Y[[s]], ndY, fallback = mean(Y_lag, na.rm = TRUE))
        Y_s <- muY + stats::rnorm(n, mean = 0, sd = fits$sigma_Y[s])

        L_lag <- L_s
        Y_lag <- Y_s
        M1_lag <- M1_s
        M2_lag <- M2_s
      }
    }

    A_t <- as.numeric(dat_wide[[varmap$A[t]]])
    nd1 <- data.frame(
      A = A_t, W1 = dat_wide$W1, W2 = dat_wide$W2,
      L_lag = L_lag, Y_lag = Y_lag, M1_lag = M1_lag, M2_lag = M2_lag
    )
    mu1 <- .msm_predict_lm_numeric(fits$M1[[t]], nd1, fallback = mean(obs_m1, na.rm = TRUE))
    d1 <- .dens_norm(obs_m1, mu1, .safe_sigma(fits$sigma_M1[t]))

    nd2_joint <- data.frame(
      A = A_t, W1 = dat_wide$W1, W2 = dat_wide$W2,
      L_lag = L_lag, Y_lag = Y_lag, M1_lag = M1_lag, M2_lag = M2_lag,
      M1 = obs_m1
    )
    mu2_joint <- .msm_predict_lm_numeric(fits$M2[[t]], nd2_joint, fallback = mean(obs_m2, na.rm = TRUE))
    d2_joint <- .dens_norm(obs_m2, mu2_joint, .safe_sigma(fits$sigma_M2[t]))

    m1_draw <- mu1 + stats::rnorm(n, mean = 0, sd = fits$sigma_M1[t])
    nd2_marg <- nd2_joint
    nd2_marg$M1 <- m1_draw
    mu2_marg <- .msm_predict_lm_numeric(fits$M2[[t]], nd2_marg, fallback = mean(obs_m2, na.rm = TRUE))
    d2_marg <- .dens_norm(obs_m2, mu2_marg, .safe_sigma(fits$sigma_M2[t]))

    d_m1_sum <- d_m1_sum + d1
    d_m2_marg_sum <- d_m2_marg_sum + d2_marg
    d_joint_sum <- d_joint_sum + d1 * d2_joint
  }

  list(
    M1_marginal = d_m1_sum / B_mc,
    M2_marginal = d_m2_marg_sum / B_mc,
    M12_joint = d_joint_sum / B_mc
  )
}

.msm_predict_mediator_moments_from_conditional_models <- function(dat_wide, T, reg_a, reg_as,
                                                                  fits,
                                                                  B_mc = 5000L,
                                                                  seed = NULL) {
  B_mc <- .msm_as_positive_integer(B_mc, name = "mediator_mc_n")
  if (!is.null(seed)) set.seed(as.integer(seed)[1L])

  baseline <- extract_baseline(dat_wide, default_varmap(T))
  V <- .baseline_V_names()
  base_means <- colMeans(dat_wide[, V, drop = FALSE])

  zero <- rep(0, T)
  acc <- list(
    sep = list(
      mu1_a = zero, mu1_as = zero, mu2_a = zero, mu2_as = zero,
      mu12_a = zero, mu12_as = zero, mu12_a_aas = zero
    ),
    joint = list(
      mu1_a = zero, mu1_as = zero, mu2_a = zero, mu2_as = zero,
      mu12_a = zero, mu12_as = zero
    )
  )

  for (b in seq_len(B_mc)) {
    joint_a <- simulate_full_paths_once(baseline, fits, reg_a)
    joint_as <- simulate_full_paths_once(baseline, fits, reg_as)

    sep_m1_a <- simulate_full_paths_once(baseline, fits, reg_a)$M1
    sep_m2_a <- simulate_full_paths_once(baseline, fits, reg_a)$M2
    sep_m1_as <- simulate_full_paths_once(baseline, fits, reg_as)$M1
    sep_m2_as <- simulate_full_paths_once(baseline, fits, reg_as)$M2

    acc$joint$mu1_a <- acc$joint$mu1_a + colMeans(joint_a$M1)
    acc$joint$mu2_a <- acc$joint$mu2_a + colMeans(joint_a$M2)
    acc$joint$mu12_a <- acc$joint$mu12_a + colMeans(joint_a$M1 * joint_a$M2)
    acc$joint$mu1_as <- acc$joint$mu1_as + colMeans(joint_as$M1)
    acc$joint$mu2_as <- acc$joint$mu2_as + colMeans(joint_as$M2)
    acc$joint$mu12_as <- acc$joint$mu12_as + colMeans(joint_as$M1 * joint_as$M2)

    acc$sep$mu1_a <- acc$sep$mu1_a + colMeans(sep_m1_a)
    acc$sep$mu2_a <- acc$sep$mu2_a + colMeans(sep_m2_a)
    acc$sep$mu12_a <- acc$sep$mu12_a + colMeans(sep_m1_a * sep_m2_a)
    acc$sep$mu1_as <- acc$sep$mu1_as + colMeans(sep_m1_as)
    acc$sep$mu2_as <- acc$sep$mu2_as + colMeans(sep_m2_as)
    acc$sep$mu12_as <- acc$sep$mu12_as + colMeans(sep_m1_as * sep_m2_as)
    acc$sep$mu12_a_aas <- acc$sep$mu12_a_aas + colMeans(sep_m1_a * sep_m2_as)
  }

  acc$sep <- lapply(acc$sep, function(x) x / B_mc)
  acc$joint <- lapply(acc$joint, function(x) x / B_mc)
  c(
    list(
      engine = "conditional_model_integration_for_interventional_mediator_distributions",
      base_means = base_means
    ),
    acc
  )
}

.msm_predict_mediator_moments_sequential <- .msm_predict_mediator_moments_from_conditional_models


# Return the value of a coefficient if present, else 0.
.get_coef0 <- function(coefs, name) {
  if (is.null(coefs) || is.null(names(coefs))) return(0)

  val <- NA_real_

  if (name %in% names(coefs)) {
    val <- unname(coefs[[name]])
  } else if (grepl(":", name, fixed = TRUE)) {
    # R may order interaction terms as M1:A instead of A:M1.
    parts <- strsplit(name, ":", fixed = TRUE)[[1]]
    alt <- paste(rev(parts), collapse = ":")
    if (alt %in% names(coefs)) val <- unname(coefs[[alt]])
  } else {
    return(0)
  }

  if (!is.finite(val)) return(0)
  val
}

# Baseline covariates V used in all MSMs
.baseline_V_names <- function() c("W1", "W2", "Y0", "M1_0", "M2_0")

# Convenience: safe subset of a varmap vector by 1:k, returning character(0) if k<1
.take_first <- function(x, k) {
  if (k < 1) return(character(0))
  x[seq_len(min(k, length(x)))]
}

# -----------------------------------------------------------------------------
# Weight estimation
# -----------------------------------------------------------------------------

# Treatment stabilized weights.
#
# - wA_y: used for OUTCOME MSM; numerator excludes (L,Y) history.
# - wA_m: used for MEDIATOR MSM; numerator excludes mediator history and (L,Y).
.compute_treatment_weights <- function(long, T,
                                       treat_mech = c("baseline_rct","sequential_rct","observational"),
                                       p_rct = 0.5,
                                       eps = 1e-9) {
  treat_mech <- match.arg(treat_mech)

  n <- length(unique(long$id))
  ids <- sort(unique(long$id))

  # Default: all ones
  long$wA_y_step <- 1
  long$wA_m_step <- 1

  if (treat_mech != "observational") {
    # Under (sequential|baseline) RCT mechanisms, the stabilized ratios are 1.
    # (baseline_rct additionally implies static A, but we don’t rely on that here.)
    long$wA_y_cum <- ave(long$wA_y_step, long$id, FUN = cumprod)
    long$wA_m_cum <- ave(long$wA_m_step, long$id, FUN = cumprod)
    return(long)
  }

  varmap <- default_varmap(T)
  V <- .baseline_V_names()

  # Build time-specific models for A (binary) and compute stabilized ratios.
  for (t in 1:T) {
    dt <- long[long$t == t, , drop = FALSE]
    if (nrow(dt) != n) {
      # Reorder/align if needed
      dt <- dt[match(ids, dt$id), , drop = FALSE]
    }

    k <- t - 1L  # covariate time

    # Past A: A0..A_{k-1}  (length k)
    A_past <- .take_first(varmap$A, k)

    # Past mediators (exclude baseline): M1_1..M1_{k-1}
    M1_past <- if (k <= 1) character(0) else varmap$M1[2:k]
    M2_past <- if (k <= 1) character(0) else varmap$M2[2:k]

    # Past confounders: L0..L_{k-1} and Y_1..Y_k
    L_past <- .take_first(varmap$L, k)
    Y_past <- .take_first(varmap$Y, k)

    # Denominator (for both outcome and mediator MSM weights)
    rhs_den <- c(V, A_past, M1_past, M2_past, L_past, Y_past)
    rhs_den <- unique(rhs_den)
    f_den <- stats::as.formula(paste("A ~", paste(rhs_den, collapse = " + ")))

    # Numerator for OUTCOME MSM treatment weight: exclude (L,Y) history
    rhs_num_y <- unique(c(V, A_past, M1_past, M2_past))
    f_num_y <- stats::as.formula(paste("A ~", paste(rhs_num_y, collapse = " + ")))

    # Numerator for MEDIATOR MSM treatment weight: exclude mediator history and (L,Y)
    rhs_num_m <- unique(c(V, A_past))
    f_num_m <- stats::as.formula(paste("A ~", paste(rhs_num_m, collapse = " + ")))

    fit_den <- tryCatch(stats::glm(f_den, data = dt, family = stats::binomial()),
                        error = function(e) NULL)
    fit_num_y <- tryCatch(stats::glm(f_num_y, data = dt, family = stats::binomial()),
                          error = function(e) NULL)
    fit_num_m <- tryCatch(stats::glm(f_num_m, data = dt, family = stats::binomial()),
                          error = function(e) NULL)

    if (is.null(fit_den) || is.null(fit_num_y) || is.null(fit_num_m)) {
      .warn_once(
        key = 'msm_treatment_fallback',
        '[MSM+IPW] One or more treatment models failed; using constant 0.5 fallback probabilities for the affected visit(s).',
        call. = FALSE
      )
      p_den <- rep(0.5, nrow(dt))
      p_num_y <- rep(0.5, nrow(dt))
      p_num_m <- rep(0.5, nrow(dt))
    } else {
      p_den   <- .clamp01(stats::predict(fit_den, type = "response"), eps)
      p_num_y <- .clamp01(stats::predict(fit_num_y, type = "response"), eps)
      p_num_m <- .clamp01(stats::predict(fit_num_m, type = "response"), eps)
    }

    A_obs <- dt$A
    # Individual conditional probabilities of the observed treatment.
    pr_den   <- ifelse(A_obs == 1, p_den,   1 - p_den)
    pr_num_y <- ifelse(A_obs == 1, p_num_y, 1 - p_num_y)
    pr_num_m <- ifelse(A_obs == 1, p_num_m, 1 - p_num_m)

    w_step_y <- pr_num_y / pr_den
    w_step_m <- pr_num_m / pr_den

    long$wA_y_step[long$t == t] <- w_step_y
    long$wA_m_step[long$t == t] <- w_step_m
  }

  long$wA_y_cum <- ave(long$wA_y_step, long$id, FUN = cumprod)
  long$wA_m_cum <- ave(long$wA_m_step, long$id, FUN = cumprod)
  long
}


# Stabilized mediator-process weights for fitting the single controlled-outcome MSM.
# The numerator mediator density is a stabilization choice and is not the
# interventional mediator distribution. Separate-draw and joint-draw distinctions
# enter only through the plug-in mediator moments and product moments.
.compute_mediator_weights <- function(long, T,
                                      trunc = c(0.01, 0.99),
                                      eps = 1e-12,
                                      sigma_method = c("ml", "lm"),
                                      Q_model = c("correct", "wrong"),
                                      dat_wide = NULL,
                                      mediator_density_mc_n = 5000L,
                                      seed = NULL,
                                      gcomp_fits = NULL) {
  sigma_method <- match.arg(sigma_method)
  Q_model <- match.arg(Q_model)

  varmap <- default_varmap(T)
  V <- .baseline_V_names()

  if (is.null(dat_wide)) {
    .stop("dat_wide is required for computing mediator density ratios in the marginal structural model estimator.")
  }
  if (is.null(gcomp_fits)) {
    .stop("gcomp_fits is required for computing mediator density numerators from conditional models in the marginal structural model estimator.")
  }
  mediator_density_mc_n <- .msm_as_positive_integer(
    mediator_density_mc_n,
    name = "mediator_density_mc_n"
  )

  n <- length(unique(long$id))
  ids <- sort(unique(long$id))

  long$wM_y_step <- 1

  # Only post-baseline mediators: times t = 2..T correspond to M*_1..M*_{T-1}
  if (T >= 2L) for (t in 2:T) {
    dt <- long[long$t == t, , drop = FALSE]
    if (nrow(dt) != n) {
      dt <- dt[match(ids, dt$id), , drop = FALSE]
    }

    k <- t - 1L  # mediator time index in wide columns

    # Exposure history up to current time: A0..A_k  (length k+1)
    A_hist <- .take_first(varmap$A, k + 1L)

    # Previous mediators (exclude baseline): M*_1..M*_{k-1}
    M1_prev <- if (k <= 1) character(0) else varmap$M1[2:k]
    M2_prev <- if (k <= 1) character(0) else varmap$M2[2:k]

    # Confounder history for mediator models: L0..L_{k-1} and Y_1..Y_k
    L_hist <- .take_first(varmap$L, k)
    Y_hist <- .take_first(varmap$Y, k)

    # -----------------------
    # M1 density models
    # -----------------------
    rhs_m1_den <- unique(c(V, A_hist, M1_prev, M2_prev, L_hist, Y_hist))
    f_m1_den <- stats::as.formula(paste("M1 ~", paste(rhs_m1_den, collapse = " + ")))

    fit_m1_den <- tryCatch(stats::lm(f_m1_den, data = dt), error = function(e) NULL)

    mu_m1_den <- if (is.null(fit_m1_den)) rep(mean(dt$M1), nrow(dt)) else as.numeric(stats::predict(fit_m1_den, newdata = dt))

    if (sigma_method == "ml") {
      sd_m1_den <- if (is.null(fit_m1_den)) stats::sd(dt$M1) else .ml_sigma(fit_m1_den)
    } else {
      sd_m1_den <- if (is.null(fit_m1_den)) stats::sd(dt$M1) else summary(fit_m1_den)$sigma
    }

    sd_m1_den <- .safe_sigma(sd_m1_den)

    d_m1_den <- .dens_norm(dt$M1, mu_m1_den, sd_m1_den)

    # -----------------------
    # M2 density model in the observed ordered mediator process
    # -----------------------
    if (Q_model == "correct") {
      dt$AM1_curr <- dt$A * dt$M1
      rhs_m2_den_joint <- unique(c(V, A_hist, M1_prev, M2_prev, "M1", "AM1_curr", L_hist, Y_hist))
    } else {
      rhs_m2_den_joint <- unique(c(V, A_hist, M1_prev, M2_prev, L_hist, Y_hist))
    }
    f_m2_den_joint <- stats::as.formula(paste("M2 ~", paste(rhs_m2_den_joint, collapse = " + ")))

    fit_m2_den_joint <- tryCatch(stats::lm(f_m2_den_joint, data = dt), error = function(e) NULL)

    failed_models <- is.null(fit_m1_den) || is.null(fit_m2_den_joint)
    if (failed_models) {
      .warn_once(
        key = 'msm_mediator_fallback',
        '[MSM+IPW] One or more mediator density models failed; using simple mean/empirical-sd fallback(s) for the affected visit(s).',
        call. = FALSE
      )
    }

    mu_m2_den_joint <- if (is.null(fit_m2_den_joint)) rep(mean(dt$M2), nrow(dt)) else as.numeric(stats::predict(fit_m2_den_joint, newdata = dt))

    if (sigma_method == "ml") {
      sd_m2_den_joint <- if (is.null(fit_m2_den_joint)) stats::sd(dt$M2) else .ml_sigma(fit_m2_den_joint)
    } else {
      sd_m2_den_joint <- if (is.null(fit_m2_den_joint)) stats::sd(dt$M2) else summary(fit_m2_den_joint)$sigma
    }

    sd_m2_den_joint <- .safe_sigma(sd_m2_den_joint)

    d_m2_den_joint <- .dens_norm(dt$M2, mu_m2_den_joint, sd_m2_den_joint)

    den_step <- pmax(d_m1_den * d_m2_den_joint, eps)
    num <- .msm_mediator_density_numerators_from_conditional_models(
      dat_wide = dat_wide[dt$id, , drop = FALSE],
      T = T,
      fits = gcomp_fits,
      t = t,
      B_mc = mediator_density_mc_n,
      seed = if (is.null(seed)) NULL else as.integer(seed)[1L] + 1009L * t
    )
    w_step_y <- num$M12_joint / den_step

    # NOTE: In the current implementation, truncation is applied to the
    # end-of-follow-up cumulative weights used in the outcome MSM rather than to
    # the step-wise mediator density ratios at each visit.

    long$wM_y_step[long$t == t] <- w_step_y
  }

  long$wM_y_cum <- ave(long$wM_y_step, long$id, FUN = cumprod)
  long
}


# -----------------------------------------------------------------------------
# MSM fitting
# -----------------------------------------------------------------------------

# Fit outcome MSM (terminal Y_T) given weights.
.fit_outcome_msm <- function(dat_wide, T, weights_end,
                             treat_mech = c("baseline_rct","sequential_rct","observational"),
                             Q_model = c("correct","wrong")) {
  treat_mech <- match.arg(treat_mech)
  Q_model <- match.arg(Q_model)

  varmap <- default_varmap(T)
  V <- .baseline_V_names()

  # Outcome
  y_name <- varmap$Y[T]

  # Treatment history A0..A_{T-1}
  A_terms <- varmap$A

  # NOTE:
  # treat_mech == 'baseline_rct' implies A0==A1==...==A_{T-1} in the DGP.
  # Keeping all A-terms makes the outcome MSM design matrix rank-deficient,
  # which yields NA (aliased) coefficients in lm(). This estimator later
  # uses a coefficient plug-in; NA would propagate to the world means.
  #
  # Using A0 only is sufficient/identifiable under baseline_rct.
  if (treat_mech == "baseline_rct") {
    A_terms <- varmap$A[1]
  }

  # Mediator history excluding baseline: M*_1..M*_{T-1}
  M1_terms <- if (T <= 1) character(0) else varmap$M1[2:T]
  M2_terms <- if (T <= 1) character(0) else varmap$M2[2:T]

  # Interactions between A_k and current mediators at each covariate time k.
  int_terms <- character(0)
  if (Q_model == "correct" && T >= 2L) {
    for (k in 1:(T - 1L)) {
      a_nm <- if (treat_mech == "baseline_rct") varmap$A[1L] else varmap$A[k + 1L]
      int_terms <- c(
        int_terms,
        paste0(a_nm, ":", varmap$M1[k + 1L]),
        paste0(a_nm, ":", varmap$M2[k + 1L])
      )
    }
  }

  # Optional: mediator-mediator interaction terms (M1_k * M2_k) in the outcome MSM.
  # This matches the "Q_model correct vs wrong" definition used by the other estimators.
  dat_msm <- dat_wide
  mm_terms <- character(0)
  if (Q_model == "correct" && T >= 2L) {
    for (k in 1:(T - 1L)) {
      m1_nm <- varmap$M1[k + 1L]
      m2_nm <- varmap$M2[k + 1L]
      mm_nm <- paste0("MM", k)
      dat_msm[[mm_nm]] <- dat_msm[[m1_nm]] * dat_msm[[m2_nm]]
      mm_terms <- c(mm_terms, mm_nm)
    }
  }

  rhs <- c(V, A_terms, M1_terms, M2_terms, int_terms, mm_terms)
  rhs <- unique(rhs)
  f <- stats::as.formula(paste(y_name, "~", paste(rhs, collapse = " + ")))

  fit <- stats::lm(f, data = dat_msm, weights = weights_end)
  th <- coef(fit)
  # Make plug-in robust to rank deficiency (aliased coefficients)
  th[!is.finite(th)] <- 0
  th
}

# Fit a treatment-only outcome MSM (terminal Y_T) for the natural (TE) worlds.
#
# Rationale:
# - The natural mean E{Y_T(a_bar)} does not involve mediator interventions.
# - With post-treatment confounding, the natural total effect (TE) generally
#   differs from the interventional total effect (ITE) defined by mediator
#   interventions (Tai et al., 2023). Therefore we estimate mu_nat_* using a
#   separate treatment MSM/IPW fit, rather than forcing mu_nat == mu_joint.
.fit_outcome_msm_natural <- function(dat_wide, T, weights_end,
                                     treat_mech = c("baseline_rct","sequential_rct","observational")) {
  treat_mech <- match.arg(treat_mech)

  varmap <- default_varmap(T)
  V <- .baseline_V_names()

  # Outcome
  y_name <- varmap$Y[T]

  # Treatment history A0..A_{T-1}
  A_terms <- varmap$A
  if (treat_mech == "baseline_rct") {
    # Under baseline randomization, A0==...==A_{T-1} in the DGP; using A0 only
    # avoids rank deficiency in lm().
    A_terms <- varmap$A[1]
  }

  rhs <- unique(c(V, A_terms))
  f <- stats::as.formula(paste(y_name, "~", paste(rhs, collapse = " + ")))

  fit <- stats::lm(f, data = dat_wide, weights = weights_end)
  th <- coef(fit)
  th[!is.finite(th)] <- 0
  th
}

# Plug-in natural world mean from the treatment-only outcome MSM.
.plugin_world_mean_natural <- function(theta, T, outer, base_means) {
  varmap <- default_varmap(T)
  V <- .baseline_V_names()

  out <- .get_coef0(theta, "(Intercept)")

  # Baseline covariates
  for (v in V) {
    out <- out + .get_coef0(theta, v) * base_means[[v]]
  }

  # Treatment history main effects
  for (k in 0:(T - 1L)) {
    a_nm <- varmap$A[k + 1L]
    out <- out + .get_coef0(theta, a_nm) * outer[k + 1L]
  }

  as.numeric(out)
}





# Fit/evaluate mediator trajectory laws and compute E(M_k | regimen).
#
# Returns a list with expected means for each regimen:
#   list(
#     sep  = list(mu1_a, mu1_as, mu2_a, mu2_as),
#     joint= list(mu1_a, mu1_as, mu2_a, mu2_as)
#   )

.fit_and_predict_mediator_moments <- function(long, dat_wide, T, reg_a, reg_as,
                                             weights_long_m = NULL,
                                             sigma_method = c("ml","lm"),
                                             Q_model = c("correct", "wrong"),
                                             mediator_mc_n = 5000L,
                                             seed = NULL,
                                             gcomp_fits = NULL) {
  sigma_method <- match.arg(sigma_method)
  Q_model <- match.arg(Q_model)
  if (!is.null(weights_long_m)) {
    # Conditional-model integration does not use person-visit mediator MSM weights.
  }

  if (is.null(gcomp_fits)) {
    gcomp_fits <- .msm_gcomp_fits_for_mediator_law(
      dat_wide = dat_wide,
      T = T,
      Q_model = Q_model
    )
  }

  return(.msm_predict_mediator_moments_from_conditional_models(
    dat_wide = dat_wide,
    T = T,
    reg_a = reg_a,
    reg_as = reg_as,
    fits = gcomp_fits,
    B_mc = mediator_mc_n,
    seed = seed
  ))
}

.fit_and_predict_mediator_means <- .fit_and_predict_mediator_moments

# Plug-in world mean based on outcome MSM coefficients and E(M).
#
# Arguments:
#   - theta: named coefficient vector from the outcome MSM.
#   - outer: length-T vector of exposure values (outer regimen).
#   - mu1, mu2: length-T vectors giving expected mediator values at times
#               k=0..T-1 (mu[1] is baseline k=0).
#   - base_means: named baseline means for V.
.plugin_world_mean <- function(theta, T, outer, mu1, mu2, base_means,
                               mu12 = NULL,
                               Q_model = c("correct","wrong"),
                               treat_mech = c("baseline_rct", "sequential_rct", "observational")) {

  Q_model <- match.arg(Q_model)
  treat_mech <- match.arg(treat_mech)

  varmap <- default_varmap(T)
  V <- .baseline_V_names()

  # Start with intercept
  out <- .get_coef0(theta, "(Intercept)")

  # Baseline covariates
  for (v in V) {
    out <- out + .get_coef0(theta, v) * base_means[[v]]
  }

  # Treatment history main effects
  for (k in 0:(T - 1L)) {
    a_nm <- varmap$A[k + 1L]
    out <- out + .get_coef0(theta, a_nm) * outer[k + 1L]
  }

  # Post-baseline times k=1..T-1
  if (T >= 2) {
    for (k in 1:(T - 1L)) {
      m1_nm <- varmap$M1[k + 1L]
      m2_nm <- varmap$M2[k + 1L]
      a_nm <- if (treat_mech == "baseline_rct") varmap$A[1L] else varmap$A[k + 1L]

      out <- out + .get_coef0(theta, m1_nm) * mu1[k + 1L]
      out <- out + .get_coef0(theta, m2_nm) * mu2[k + 1L]

      out <- out + .get_coef0(theta, paste0(a_nm, ":", m1_nm)) * (outer[k + 1L] * mu1[k + 1L])
      out <- out + .get_coef0(theta, paste0(a_nm, ":", m2_nm)) * (outer[k + 1L] * mu2[k + 1L])
    }
  }

  # Mediator-mediator interaction terms (M1_k * M2_k) when Q_model == "correct".
  if (Q_model == "correct") {
    if (is.null(mu12) || length(mu12) != T) {
      .stop("mu12 must be a length-T vector when Q_model='correct'.")
    }

    if (T >= 2) {
      for (k in 1:(T - 1L)) {
        mm_nm <- paste0("MM", k)
        out <- out + .get_coef0(theta, mm_nm) * mu12[k + 1L]
      }
    }
  }

  as.numeric(out)
}



# -----------------------------------------------------------------------------
# Public API
# -----------------------------------------------------------------------------

#' MSM + IPW world-mean estimator (coefficient plug-in)
#'
#' Returns the 9 world means in the package's canonical keys.
#'
#' @param dat observed data in wide format
#' @param dat_wide same as dat (legacy alias)
#' @param T number of follow-up outcome times
#' @param reg_a,reg_as exposure regimens (length-T vectors or scalars)
#' @param treat_mech treatment mechanism used for fitting treatment weights
#' @param p_rct randomization probability (used only when treat_mech != observational)
#' @param trunc legacy alias for final cumulative weight quantile bounds
#' @param eps small constant for probability clamping / division safety
#' @return list(means=..., effects=..., fits=...)
msm_ipw_estimate_worldmeans <- function(dat = NULL, dat_wide = NULL,
                                       T, reg_a, reg_as,
                                       treat_mech = c("baseline_rct","sequential_rct","observational"),
                                       p_rct = 0.5,
                                       trunc = NULL,
                                       truncation_enabled = TRUE,
                                       truncation_policy = "quantile",
                                       truncation_quantile_lower = 0.01,
                                       truncation_quantile_upper = 0.99,
                                       truncation_target = "final_cumulative_weight",
                                       eps = 1e-9,
                                       sigma_method = c("ml","lm"),
                                       Q_model = c("correct","wrong"),
                                       B_mc = 5000L,
                                       mediator_density_mc_n = B_mc,
                                       seed = NULL) {
  .msm_source_if_needed()

  if (is.null(dat_wide)) dat_wide <- dat
  if (is.null(dat_wide)) stop("Provide `dat` (wide data).", call. = FALSE)

  treat_mech <- match.arg(treat_mech)
  sigma_method <- match.arg(sigma_method)
  Q_model <- match.arg(Q_model)
  B_mc <- .msm_as_positive_integer(B_mc, name = "B_mc")
  mediator_density_mc_n <- .msm_as_positive_integer(
    mediator_density_mc_n,
    name = "mediator_density_mc_n"
  )
  truncation <- .msm_normalize_truncation(
    trunc = trunc,
    truncation_enabled = truncation_enabled,
    truncation_policy = truncation_policy,
    truncation_quantile_lower = truncation_quantile_lower,
    truncation_quantile_upper = truncation_quantile_upper,
    truncation_target = truncation_target
  )
  estimator_variant <- if (isTRUE(truncation$enabled)) {
    "msm_ipw_quantile_truncated"
  } else {
    "msm_ipw_untruncated"
  }
  effective_trunc <- c(truncation$quantile_lower, truncation$quantile_upper)

  # Normalize regimens
  reg_a  <- normalize_regimen(reg_a,  T, "reg_a")
  reg_as <- normalize_regimen(reg_as, T, "reg_as")
  validate_baseline_rct_regimens(reg_a, reg_as, treat_mech)

  # Basic validation
  validate_wide_data(dat_wide, T)

  # Long data with history columns
  long <- wide_to_long_with_hist(dat_wide, T)

  gcomp_fits <- .msm_gcomp_fits_for_mediator_law(
    dat_wide = dat_wide,
    T = T,
    Q_model = Q_model
  )

  # -------------------------------------------------------------------------
  # 1) Treatment weights
  # -------------------------------------------------------------------------
  long <- .compute_treatment_weights(long, T, treat_mech = treat_mech, p_rct = p_rct, eps = eps)

  # -------------------------------------------------------------------------
  # 2) Common mediator-process weights for the controlled-outcome MSM
  # -------------------------------------------------------------------------
  long <- .compute_mediator_weights(
    long,
    T,
    trunc = effective_trunc,
    eps = eps,
    sigma_method = sigma_method,
    Q_model = Q_model,
    dat_wide = dat_wide,
    mediator_density_mc_n = mediator_density_mc_n,
    seed = seed,
    gcomp_fits = gcomp_fits
  )

  # Extract end-of-follow-up weights per subject (t == T)
  end <- long[long$t == T, c("id","wA_y_cum","wA_m_cum","wM_y_cum"), drop = FALSE]
  end <- end[order(end$id), ]

  w_out_nat_raw <- end$wA_y_cum
  w_out_int_raw <- w_out_nat_raw * end$wM_y_cum

  # Cumulative truncation: truncate the *end-of-follow-up* weights used in the
  # outcome MSM (rather than truncating step-wise ratios).
  w_out_nat_obj <- .msm_truncate_final_weight(w_out_nat_raw, truncation)
  w_out_int_obj <- .msm_truncate_final_weight(w_out_int_raw, truncation)
  w_out_nat <- w_out_nat_obj$value
  w_out_int <- w_out_int_obj$value

  # -------------------------------------------------------------------------
  # 3) Outcome MSMs (theta)
  #    - theta_nat   : treatment-only MSM for natural (TE) worlds
  #    - theta_int   : single controlled-outcome MSM for interventional worlds
  # -------------------------------------------------------------------------
  theta_nat <- .fit_outcome_msm_natural(
    dat_wide, T,
    weights_end = w_out_nat,
    treat_mech = treat_mech
  )
  theta_int <- .fit_outcome_msm(
    dat_wide, T,
    weights_end = w_out_int,
    treat_mech = treat_mech,
    Q_model = Q_model
  )

  # -------------------------------------------------------------------------
  # 4) Mediator moments and product moments
  # -------------------------------------------------------------------------
  med_obj <- .fit_and_predict_mediator_moments(
    long = long,
    dat_wide = dat_wide,
    T = T,
    reg_a = reg_a,
    reg_as = reg_as,
    weights_long_m = long$wA_m_cum,
    sigma_method = sigma_method,
    Q_model = Q_model,
    mediator_mc_n = B_mc,
    seed = seed,
    gcomp_fits = gcomp_fits
  )
  if (!identical(med_obj$engine, "conditional_model_integration_for_interventional_mediator_distributions")) {
    .stop("MSM/IPW mediator moments must be evaluated from the fitted interventional mediator distributions.")
  }

  base_means <- med_obj$base_means

  # Convenience aliases
  mu1_sep_a   <- med_obj$sep$mu1_a
  mu1_sep_as  <- med_obj$sep$mu1_as
  mu2_sep_a   <- med_obj$sep$mu2_a
  mu2_sep_as  <- med_obj$sep$mu2_as

  mu1_joint_a  <- med_obj$joint$mu1_a
  mu1_joint_as <- med_obj$joint$mu1_as
  mu2_joint_a  <- med_obj$joint$mu2_a
  mu2_joint_as <- med_obj$joint$mu2_as

  mu12_joint_a  <- med_obj$joint$mu12_a
  mu12_joint_as <- med_obj$joint$mu12_as

  mu12_sep_a     <- med_obj$sep$mu12_a
  mu12_sep_as    <- med_obj$sep$mu12_as
  mu12_sep_a_aas <- med_obj$sep$mu12_a_aas

  # -------------------------------------------------------------------------
  # 5) Build the 9 world means via coefficient plug-in
  # -------------------------------------------------------------------------

  # Natural worlds (TE): treatment-only MSM plug-in (no mediator interventions)
  mu_nat_a  <- .plugin_world_mean_natural(theta_nat, T, outer = reg_a,
                                          base_means = base_means)
  mu_nat_as <- .plugin_world_mean_natural(theta_nat, T, outer = reg_as,
                                          base_means = base_means)

  # Joint mediator intervention worlds (Tai)
  mu_joint_aa   <- .plugin_world_mean(theta_int, T, outer = reg_a,
                                      mu1 = mu1_joint_a,  mu2 = mu2_joint_a,
                                      base_means = base_means,
                                      mu12 = if (Q_model == "correct") mu12_joint_a else NULL,
                                      Q_model = Q_model,
                                      treat_mech = treat_mech)
  mu_joint_asas <- .plugin_world_mean(theta_int, T, outer = reg_as,
                                      mu1 = mu1_joint_as, mu2 = mu2_joint_as,
                                      base_means = base_means,
                                      mu12 = if (Q_model == "correct") mu12_joint_as else NULL,
                                      Q_model = Q_model,
                                      treat_mech = treat_mech)
  mu_joint_aas  <- .plugin_world_mean(theta_int, T, outer = reg_a,
                                      mu1 = mu1_joint_as, mu2 = mu2_joint_as,
                                      base_means = base_means,
                                      mu12 = if (Q_model == "correct") mu12_joint_as else NULL,
                                      Q_model = Q_model,
                                      treat_mech = treat_mech)

  # Separate mediator intervention worlds (Yamamuro)
  mu_sep_aaa <- .plugin_world_mean(theta_int, T, outer = reg_a,
                                   mu1 = mu1_sep_a, mu2 = mu2_sep_a,
                                   base_means = base_means,
                                   mu12 = if (Q_model == "correct") mu12_sep_a else NULL,
                                   Q_model = Q_model,
                                   treat_mech = treat_mech)
  mu_sep_asas_asas <- .plugin_world_mean(theta_int, T, outer = reg_as,
                                         mu1 = mu1_sep_as, mu2 = mu2_sep_as,
                                         base_means = base_means,
                                         mu12 = if (Q_model == "correct") mu12_sep_as else NULL,
                                         Q_model = Q_model,
                                         treat_mech = treat_mech)
  mu_sep_a_asas <- .plugin_world_mean(theta_int, T, outer = reg_a,
                                      mu1 = mu1_sep_as, mu2 = mu2_sep_as,
                                      base_means = base_means,
                                      mu12 = if (Q_model == "correct") mu12_sep_as else NULL,
                                      Q_model = Q_model,
                                      treat_mech = treat_mech)
  mu_sep_a_aas <- .plugin_world_mean(theta_int, T, outer = reg_a,
                                     mu1 = mu1_sep_a, mu2 = mu2_sep_as,
                                     base_means = base_means,
                                     mu12 = if (Q_model == "correct") mu12_sep_a_aas else NULL,
                                     Q_model = Q_model,
                                     treat_mech = treat_mech)

  means <- list(
    # keep the exact key order used everywhere else
    mu_nat_a = mu_nat_a,
    mu_nat_as = mu_nat_as,

    mu_joint_aa = mu_joint_aa,
    mu_joint_asas = mu_joint_asas,
    mu_joint_aas = mu_joint_aas,

    mu_sep_aaa = mu_sep_aaa,
    mu_sep_asas_asas = mu_sep_asas_asas,
    mu_sep_a_asas = mu_sep_a_asas,
    mu_sep_a_aas = mu_sep_a_aas
  )

  effects <- compute_estimands_from_means(means)

  # ---- Diagnostics (lightweight) --------------------------------------------
  .diag_stats_local <- function(x, probs = c(0,0.01,0.05,0.5,0.95,0.99,1), type = 8) {
    x <- x[is.finite(x)]
    if (!length(x)) {
      out <- rep(NA_real_, 9)
      names(out) <- c("min","p01","p05","median","p95","p99","max","mean","sd")
      return(out)
    }
    qs <- as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = type))
    c(min=qs[1], p01=qs[2], p05=qs[3], median=qs[4], p95=qs[5], p99=qs[6], max=qs[7],
      mean=mean(x), sd=stats::sd(x))
  }

  .diag_ess_local <- function(w) {
    w <- w[is.finite(w)]
    n <- length(w)
    if (n == 0) return(c(ess=NA_real_, ess_frac=NA_real_))
    sw <- sum(w); sw2 <- sum(w*w)
    if (!is.finite(sw) || !is.finite(sw2) || sw2 <= 0) return(c(ess=NA_real_, ess_frac=NA_real_))
    ess <- (sw*sw)/sw2
    c(ess=ess, ess_frac=ess/n)
  }

  weight_cols <- c("wA_y_step","wA_y_cum","wA_m_step","wA_m_cum",
                   "wM_y_step","wM_y_cum")

  wt_rows <- list()
  idx <- 0L
  for (tt in 1:T) {
    dt <- long[long$t == tt, , drop = FALSE]
    for (wc in weight_cols) {
      idx <- idx + 1L
      s <- .diag_stats_local(dt[[wc]])
      e <- .diag_ess_local(dt[[wc]])
      wt_rows[[idx]] <- data.frame(
        t = tt,
        weight = wc,
        n = length(dt[[wc]]),
        trunc_lo = NA_real_,
        trunc_hi = NA_real_,
        prop_truncated = NA_real_,
        t(s),
        ess = e["ess"],
        ess_frac = e["ess_frac"],
        stringsAsFactors = FALSE
      )
    }
  }
  diag_weights_t <- do.call(rbind, wt_rows)

  b_nat <- c(w_out_nat_obj$lower, w_out_nat_obj$upper)
  b_int <- c(w_out_int_obj$lower, w_out_int_obj$upper)

  w_out_nat_trunc <- w_out_nat
  w_out_int_trunc <- w_out_int

  prop_tr <- function(x, lo, hi) mean(x < lo | x > hi, na.rm = TRUE)

  # Truncation cutpoints + truncation rate (computed on raw end-of-follow-up weights)
  prop_trunc_nat <- if (isTRUE(truncation$enabled)) prop_tr(w_out_nat_raw, b_nat[1], b_nat[2]) else 0
  prop_trunc_int <- if (isTRUE(truncation$enabled)) prop_tr(w_out_int_raw, b_int[1], b_int[2]) else 0

  mk_end <- function(name, x, lo, hi, prop_trunc_raw = NA_real_) {
    s <- .diag_stats_local(x)
    e <- .diag_ess_local(x)
    data.frame(
      t = NA_integer_,
      weight = name,
      n = length(x),
      trunc_lo = lo,
      trunc_hi = hi,
      # Proportion of raw end-of-follow-up weights outside [trunc_lo, trunc_hi].
      # This is the key truncation diagnostic for IPW.
      prop_trunc_raw = prop_trunc_raw,
      t(s),
      ess = e["ess"], ess_frac = e["ess_frac"],
      stringsAsFactors = FALSE
    )
  }

  diag_weights_end <- rbind(
    mk_end("w_out_nat_raw",   w_out_nat_raw,   b_nat[1], b_nat[2], prop_trunc_raw = prop_trunc_nat),
    mk_end("w_out_nat_trunc", w_out_nat_trunc, b_nat[1], b_nat[2], prop_trunc_raw = prop_trunc_nat),
    mk_end("w_out_int_raw",   w_out_int_raw,   b_int[1], b_int[2], prop_trunc_raw = prop_trunc_int),
    mk_end("w_out_int_trunc", w_out_int_trunc, b_int[1], b_int[2], prop_trunc_raw = prop_trunc_int)
  )
  diag_truncation <- rbind(
    .msm_truncation_summary_row(
      raw = w_out_nat_raw,
      truncated = w_out_nat_trunc,
      truncation = truncation,
      estimator_variant = estimator_variant,
      component = "natural_outcome_msm",
      task_id = "msm_ipw::w_out_nat",
      q_lower_value = b_nat[1],
      q_upper_value = b_nat[2]
    ),
    .msm_truncation_summary_row(
      raw = w_out_int_raw,
      truncated = w_out_int_trunc,
      truncation = truncation,
      estimator_variant = estimator_variant,
      component = "interventional_outcome_msm",
      task_id = "msm_ipw::w_out_int",
      q_lower_value = b_int[1],
      q_upper_value = b_int[2]
    )
  )

  diag <- list(
    treat_mech = treat_mech,
    p_rct = p_rct,
    trunc = effective_trunc,
    truncation_enabled = isTRUE(truncation$enabled),
    truncation_policy = truncation$policy,
    truncation_target = truncation$target,
    estimator_variant = estimator_variant,
    mediator_law = "conditional_model_integration_for_interventional_mediator_distributions",
    mediator_intervention_distribution_source = "conditional_model_integration_for_interventional_mediator_distributions",
    mediator_mc_n = B_mc,
    mediator_density_mc_n = mediator_density_mc_n,
    weights_t = diag_weights_t,
    weights_end = diag_weights_end,
    truncation_diagnostics = diag_truncation
  )

  list(
    Q_model = Q_model,
    means = means,
    effects = effects,
    fits = list(
      theta_nat = theta_nat,
      theta_int = theta_int,
      base_means = base_means,
      mu12_joint_a = mu12_joint_a,
      mu12_joint_as = mu12_joint_as,
      mu12_sep_a = mu12_sep_a,
      mu12_sep_as = mu12_sep_as,
      mu12_sep_a_aas = mu12_sep_a_aas
    ),
    diagnostics = diag
  )
}
