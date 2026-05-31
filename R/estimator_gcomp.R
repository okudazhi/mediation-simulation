################################################################################
# estimator_gcomp.R
#
# Monte Carlo g-computation (mediational g-formula) for longitudinal mediation
# with two time-varying continuous mediators.
#
# This implementation is adapted from the provided component code
# `R_g-comp_yamamuro-tai.R`, but standardized to:
#   - use the package's common world-mean keys
#   - use the package's single estimand mapping in core_targets.R
#
# IMPORTANT:
#   - Q_model = "correct" uses the manuscript's correct node-mean models.
#   - Q_model = "wrong" implements the manuscript's combined
#     outcome-plus-mediator misspecification:
#       * remove M1 and A:M1 from the M2 model
#       * remove A:M1, A:M2, and M1:M2 from the Y model
################################################################################

# NOTE: core_utils.R and core_targets.R must be sourced before this file.

# ---- Utilities -------------------------------------------------------------
# This estimator relies on common helper functions in core_utils.R:
#   - assert_no_na(), assert_cols()
#   - default_varmap(), validate_wide_data(), extract_baseline()
#   - wide_to_long(), normalize_regimen()

#-------------------------------#
#  Fit nuisance models (by time t)
#-------------------------------#

fit_nuisance_models_gcomp <- function(long, T, Q_model = c("correct", "wrong")) {
  Q_model <- match.arg(Q_model)

  fits <- list(
    M1 = vector("list", T),
    M2 = vector("list", T),
    L  = vector("list", T),
    Y  = vector("list", T),
    sigma_M1 = rep(NA_real_, T),
    sigma_M2 = rep(NA_real_, T),
    sigma_L  = rep(NA_real_, T),
    sigma_Y  = rep(NA_real_, T),
    Q_model  = Q_model
  )

  for (t in 1:T) {
    dt <- long[long$t == t, , drop = FALSE]

    # Baseline mediators (t=1 -> covariate time 0) are treated as pre-treatment
    # history and are not modeled/simulated as post-treatment variables.
    if (t >= 2L) {
      fits$M1[[t]] <- stats::lm(
        M1 ~ A + W1 + W2 + L_lag + Y_lag + M1_lag + M2_lag,
        data = dt
      )
      fits$sigma_M1[t] <- summary(fits$M1[[t]])$sigma

      if (Q_model == "correct") {
        fits$M2[[t]] <- stats::lm(
          M2 ~ A + W1 + W2 + L_lag + Y_lag + M1_lag + M2_lag + M1 + A:M1,
          data = dt
        )
      } else {
        fits$M2[[t]] <- stats::lm(
          M2 ~ A + W1 + W2 + L_lag + Y_lag + M1_lag + M2_lag,
          data = dt
        )
      }
      fits$sigma_M2[t] <- summary(fits$M2[[t]])$sigma
    } else {
      fits$M1[[t]] <- NULL
      fits$M2[[t]] <- NULL
      fits$sigma_M1[t] <- NA_real_
      fits$sigma_M2[t] <- NA_real_
    }

    fits$L[[t]] <- stats::lm(
      L ~ A + W1 + W2 + L_lag + Y_lag + M1 + M2,
      data = dt
    )
    fits$sigma_L[t] <- summary(fits$L[[t]])$sigma

    if (Q_model == "correct") {
      fits$Y[[t]] <- stats::lm(
        Y ~ W1 + W2 + Y0 + Y_lag + L_lag + M1_lag + M2_lag + A + M1 + M2 + L + A:M1 + A:M2 + I(M1 * M2),
        data = dt
      )
    } else {
      fits$Y[[t]] <- stats::lm(
        Y ~ W1 + W2 + Y0 + Y_lag + L_lag + M1_lag + M2_lag + A + M1 + M2 + L,
        data = dt
      )
    }
    fits$sigma_Y[t] <- summary(fits$Y[[t]])$sigma

    if (t >= 2L) {
      if (!is.finite(fits$sigma_M1[t]) || fits$sigma_M1[t] <= 0) .stop("Invalid sigma_M1 at t=", t)
      if (!is.finite(fits$sigma_M2[t]) || fits$sigma_M2[t] <= 0) .stop("Invalid sigma_M2 at t=", t)
    }
    if (!is.finite(fits$sigma_L[t])  || fits$sigma_L[t]  <= 0) .stop("Invalid sigma_L at t=", t)
    if (!is.finite(fits$sigma_Y[t])  || fits$sigma_Y[t]  <= 0) .stop("Invalid sigma_Y at t=", t)
  }

  fits
}

#-------------------------------#
#  Simulation kernels
#-------------------------------#

simulate_full_paths_once <- function(baseline, fits, regimen) {
  n <- nrow(baseline)
  T <- length(regimen)

  L_lag <- rep(0, n)
  Y_lag <- baseline$Y0
  M1_lag <- baseline$M1_0
  M2_lag <- baseline$M2_0

  M1_path <- matrix(NA_real_, n, T)
  M2_path <- matrix(NA_real_, n, T)
  L_path  <- matrix(NA_real_, n, T)
  Y_path  <- matrix(NA_real_, n, T)

  for (t in 1:T) {
    A_t <- rep(regimen[t], n)

    # Baseline mediators (t=1) are fixed at their baseline values.
    if (t == 1L) {
      M1_t <- baseline$M1_0
      M2_t <- baseline$M2_0
    } else {
      nd1 <- data.frame(
        A = A_t, W1 = baseline$W1, W2 = baseline$W2,
        L_lag = L_lag, Y_lag = Y_lag, M1_lag = M1_lag, M2_lag = M2_lag
      )
      mu1 <- stats::predict(fits$M1[[t]], newdata = nd1)
      M1_t <- mu1 + stats::rnorm(n, mean = 0, sd = fits$sigma_M1[t])

      nd2 <- data.frame(
        A = A_t, W1 = baseline$W1, W2 = baseline$W2,
        L_lag = L_lag, Y_lag = Y_lag, M1_lag = M1_lag, M2_lag = M2_lag,
        M1 = M1_t
      )
      mu2 <- stats::predict(fits$M2[[t]], newdata = nd2)
      M2_t <- mu2 + stats::rnorm(n, mean = 0, sd = fits$sigma_M2[t])
    }

    ndL <- data.frame(
      A = A_t, W1 = baseline$W1, W2 = baseline$W2,
      L_lag = L_lag, Y_lag = Y_lag, M1 = M1_t, M2 = M2_t
    )
    muL <- stats::predict(fits$L[[t]], newdata = ndL)
    L_t <- muL + stats::rnorm(n, mean = 0, sd = fits$sigma_L[t])

    ndY <- data.frame(
      W1 = baseline$W1, W2 = baseline$W2,
      Y0 = baseline$Y0,
      Y_lag = Y_lag,
      L_lag = L_lag,
      M1_lag = M1_lag,
      M2_lag = M2_lag,
      A = A_t, M1 = M1_t, M2 = M2_t, L = L_t
    )
    muY <- stats::predict(fits$Y[[t]], newdata = ndY)
    Y_t <- muY + stats::rnorm(n, mean = 0, sd = fits$sigma_Y[t])

    M1_path[, t] <- M1_t
    M2_path[, t] <- M2_t
    L_path[, t]  <- L_t
    Y_path[, t]  <- Y_t

    L_lag <- L_t
    Y_lag <- Y_t
    M1_lag <- M1_t
    M2_lag <- M2_t
  }

  list(M1 = M1_path, M2 = M2_path, L = L_path, Y = Y_path)
}

simulate_LY_given_mediators_once <- function(baseline, fits, outer_regimen, M1_path, M2_path) {
  n <- nrow(baseline)
  T <- length(outer_regimen)
  if (!all(dim(M1_path) == c(n, T))) .stop("M1_path must be n x T.")
  if (!all(dim(M2_path) == c(n, T))) .stop("M2_path must be n x T.")

  L_lag <- rep(0, n)
  Y_lag <- baseline$Y0
  M1_lag <- baseline$M1_0
  M2_lag <- baseline$M2_0

  L_path <- matrix(NA_real_, n, T)
  Y_path <- matrix(NA_real_, n, T)

  for (t in 1:T) {
    A_t  <- rep(outer_regimen[t], n)
    M1_t <- M1_path[, t]
    M2_t <- M2_path[, t]

    ndL <- data.frame(
      A = A_t, W1 = baseline$W1, W2 = baseline$W2,
      L_lag = L_lag, Y_lag = Y_lag, M1 = M1_t, M2 = M2_t
    )
    muL <- stats::predict(fits$L[[t]], newdata = ndL)
    L_t <- muL + stats::rnorm(n, mean = 0, sd = fits$sigma_L[t])

    ndY <- data.frame(
      W1 = baseline$W1, W2 = baseline$W2,
      Y0 = baseline$Y0,
      Y_lag = Y_lag,
      L_lag = L_lag,
      M1_lag = M1_lag,
      M2_lag = M2_lag,
      A = A_t, M1 = M1_t, M2 = M2_t, L = L_t
    )
    muY <- stats::predict(fits$Y[[t]], newdata = ndY)
    Y_t <- muY + stats::rnorm(n, mean = 0, sd = fits$sigma_Y[t])

    L_path[, t] <- L_t
    Y_path[, t] <- Y_t

    L_lag <- L_t
    Y_lag <- Y_t
    M1_lag <- M1_t
    M2_lag <- M2_t
  }

  list(L = L_path, Y = Y_path)
}

#-------------------------------#
#  World mean estimation (MC)
#-------------------------------#

estimate_world_mean_mc <- function(baseline, fits,
                                   outer_regimen,
                                   world = c("natural", "joint", "separate"),
                                   med_regimen1 = NULL,
                                   med_regimen2 = NULL,
                                   B_mc = 5000,
                                   progress = FALSE) {
  world <- match.arg(world)
  T <- length(outer_regimen)
  if (!is.numeric(B_mc) || length(B_mc) != 1 || B_mc <= 0) .stop("B_mc must be positive.")

  mu_b <- numeric(B_mc)
  pb <- NULL
  if (isTRUE(progress)) pb <- utils::txtProgressBar(min = 0, max = B_mc, style = 3)

  for (b in 1:B_mc) {
    if (world == "natural") {
      sim <- simulate_full_paths_once(baseline, fits, outer_regimen)
      y_end <- sim$Y[, T]
    } else if (world == "joint") {
      if (is.null(med_regimen1)) .stop("joint world requires med_regimen1.")
      medsim <- simulate_full_paths_once(baseline, fits, med_regimen1)
      outsim <- simulate_LY_given_mediators_once(baseline, fits, outer_regimen,
                                                M1_path = medsim$M1,
                                                M2_path = medsim$M2)
      y_end <- outsim$Y[, T]
    } else { # separate
      if (is.null(med_regimen1) || is.null(med_regimen2)) .stop("separate world requires med_regimen1 and med_regimen2.")
      sim_M1 <- simulate_full_paths_once(baseline, fits, med_regimen1)
      sim_M2 <- simulate_full_paths_once(baseline, fits, med_regimen2)
      outsim <- simulate_LY_given_mediators_once(baseline, fits, outer_regimen,
                                                M1_path = sim_M1$M1,
                                                M2_path = sim_M2$M2)
      y_end <- outsim$Y[, T]
    }

    mu_b[b] <- mean(y_end)
    if (isTRUE(progress)) utils::setTxtProgressBar(pb, b)
  }

  if (isTRUE(progress)) close(pb)

  list(
    mean  = mean(mu_b),
    mc_se = stats::sd(mu_b) / sqrt(B_mc),
    mu_b  = mu_b
  )
}

# ---- Public API -------------------------------------------------------------

#' Estimate the 9 world means by Monte Carlo g-computation.
#'
#' @param dat wide dataset from simulate_dgp_wide()
#' @param T number of follow-up times
#' @param reg_a active regimen (length 1 or T)
#' @param reg_as reference regimen (length 1 or T)
#' @param B_mc Monte Carlo size for within-dataset forward simulation
#' @param Q_model "correct" or "wrong"
#' @param varmap optional column map
#' @param seed optional seed
#' @param progress show MC progress bars
#' @param treat_mech optional treatment mechanism; baseline_rct requires static regimens
#' @return list(means, means_mcse, effects)
.gcomp_truncation_diagnostics <- function(Q_model) {
  data.frame(
    estimator = "gcomp",
    estimator_variant = paste0("gcomp_", Q_model),
    component = "all_worldmeans",
    task_id = "gcomp::not_applicable",
    t = NA_integer_,
    node = "not_applicable",
    process_type = "not_applicable",
    truncation_enabled = FALSE,
    truncation_policy = "not_applicable",
    truncation_target = "not_applicable",
    truncation_rule = "not_applicable",
    requested_quantile_lower = NA_real_,
    requested_quantile_upper = NA_real_,
    effective_quantile_lower = NA_real_,
    effective_quantile_upper = NA_real_,
    fixed_bound_truncation_enabled = FALSE,
    fixed_bound_truncation_used = FALSE,
    density_ratio_factor_truncation_used = FALSE,
    n_values = NA_integer_,
    n_values_truncated_lower = 0L,
    n_values_truncated_upper = 0L,
    fraction_truncated_lower = 0,
    fraction_truncated_upper = 0,
    fraction_truncated_total = 0,
    min_raw_value = NA_real_,
    p01_raw_value = NA_real_,
    median_raw_value = NA_real_,
    p99_raw_value = NA_real_,
    max_raw_value = NA_real_,
    min_truncated_value = NA_real_,
    max_truncated_value = NA_real_,
    mean_raw_value = NA_real_,
    sd_raw_value = NA_real_,
    mean_truncated_value = NA_real_,
    sd_truncated_value = NA_real_,
    ess_raw = NA_real_,
    ess_truncated = NA_real_,
    passed = TRUE,
    failure_class = "no_failure",
    stringsAsFactors = FALSE
  )
}

estimate_worldmeans_gcomp <- function(dat, T,
                                     reg_a, reg_as,
                                     B_mc = 5000,
                                     Q_model = c("correct", "wrong"),
                                     varmap = default_varmap(T),
                                     seed = NULL,
                                     progress = FALSE,
                                     treat_mech = NULL) {
  Q_model <- match.arg(Q_model)
  if (!is.null(treat_mech)) {
    treat_mech <- match.arg(treat_mech, c("baseline_rct", "sequential_rct", "observational", "fixed"))
  }
  if (!is.null(seed)) set.seed(seed)

  validate_wide_data(dat, T, varmap)

  reg_a  <- normalize_regimen(reg_a,  T, "reg_a")
  reg_as <- normalize_regimen(reg_as, T, "reg_as")
  validate_baseline_rct_regimens(reg_a, reg_as, treat_mech)

  baseline <- extract_baseline(dat, varmap)
  long <- wide_to_long(dat, T, varmap)

  fits <- fit_nuisance_models_gcomp(long, T, Q_model = Q_model)

  # Natural worlds
  nat_a   <- estimate_world_mean_mc(baseline, fits, outer_regimen = reg_a,
                                   world = "natural", B_mc = B_mc, progress = progress)
  nat_as  <- estimate_world_mean_mc(baseline, fits, outer_regimen = reg_as,
                                   world = "natural", B_mc = B_mc, progress = progress)

  # Joint mediator intervention worlds
  joint_aa    <- estimate_world_mean_mc(baseline, fits, outer_regimen = reg_a,
                                       world = "joint", med_regimen1 = reg_a,
                                       B_mc = B_mc, progress = progress)
  joint_asas  <- estimate_world_mean_mc(baseline, fits, outer_regimen = reg_as,
                                       world = "joint", med_regimen1 = reg_as,
                                       B_mc = B_mc, progress = progress)
  joint_aas   <- estimate_world_mean_mc(baseline, fits, outer_regimen = reg_a,
                                       world = "joint", med_regimen1 = reg_as,
                                       B_mc = B_mc, progress = progress)

  # Separate mediator intervention worlds
  sep_aaa <- estimate_world_mean_mc(baseline, fits, outer_regimen = reg_a,
                                   world = "separate", med_regimen1 = reg_a, med_regimen2 = reg_a,
                                   B_mc = B_mc, progress = progress)
  sep_asas_asas <- estimate_world_mean_mc(baseline, fits, outer_regimen = reg_as,
                                         world = "separate", med_regimen1 = reg_as, med_regimen2 = reg_as,
                                         B_mc = B_mc, progress = progress)
  sep_a_asas <- estimate_world_mean_mc(baseline, fits, outer_regimen = reg_a,
                                       world = "separate", med_regimen1 = reg_as, med_regimen2 = reg_as,
                                       B_mc = B_mc, progress = progress)
  sep_a_aas <- estimate_world_mean_mc(baseline, fits, outer_regimen = reg_a,
                                      world = "separate", med_regimen1 = reg_a, med_regimen2 = reg_as,
                                      B_mc = B_mc, progress = progress)

  means <- list(
    mu_nat_a  = nat_a$mean,
    mu_nat_as = nat_as$mean,
    mu_joint_aa    = joint_aa$mean,
    mu_joint_asas  = joint_asas$mean,
    mu_joint_aas   = joint_aas$mean,
    mu_sep_aaa         = sep_aaa$mean,
    mu_sep_asas_asas   = sep_asas_asas$mean,
    mu_sep_a_asas      = sep_a_asas$mean,
    mu_sep_a_aas       = sep_a_aas$mean
  )

  means_mcse <- list(
    mu_nat_a  = nat_a$mc_se,
    mu_nat_as = nat_as$mc_se,
    mu_joint_aa    = joint_aa$mc_se,
    mu_joint_asas  = joint_asas$mc_se,
    mu_joint_aas   = joint_aas$mc_se,
    mu_sep_aaa         = sep_aaa$mc_se,
    mu_sep_asas_asas   = sep_asas_asas$mc_se,
    mu_sep_a_asas      = sep_a_asas$mc_se,
    mu_sep_a_aas       = sep_a_aas$mc_se
  )

  effects <- compute_estimands_from_means(means)

  list(
    estimator = paste0("gcomp_", Q_model),
    T = T,
    reg_a = reg_a,
    reg_as = reg_as,
    B_mc = B_mc,
    means = means,
    means_mcse = means_mcse,
    effects = effects,
    fits = fits,
    diagnostics = list(
      truncation_diagnostics = .gcomp_truncation_diagnostics(Q_model)
    )
  )
}

# -----------------------------------------------------------------------------
# Public alias (kept for consistency across estimators)
# -----------------------------------------------------------------------------

#' Monte Carlo g-computation world-mean estimator (public name).
#'
#' This is an alias to estimate_worldmeans_gcomp(). The wrapper accepts both
#' dat and dat_wide for convenience.
#'
#' @param dat observed data in wide format
#' @param dat_wide same as dat (alternative argument name)
#' @inheritParams estimate_worldmeans_gcomp
#' @return list with elements: means, means_mcse, effects, fits
gcomp_estimate_worldmeans <- function(dat = NULL, dat_wide = NULL,
                                     T, reg_a, reg_as,
                                     B_mc = 5000,
                                     Q_model = c("correct", "wrong"),
                                     seed = NULL,
                                     progress = FALSE,
                                     treat_mech = NULL) {
  if (is.null(dat_wide)) dat_wide <- dat
  if (is.null(dat_wide)) stop("Provide `dat` (wide data).", call. = FALSE)

  estimate_worldmeans_gcomp(
    dat = dat_wide,
    T = T,
    reg_a = reg_a,
    reg_as = reg_as,
    B_mc = B_mc,
    Q_model = Q_model,
    seed = seed,
    progress = progress,
    treat_mech = treat_mech
  )
}
