################################################################################
# truth.R
#
# Ground-truth computation for world means and estimands.
#
# Monte Carlo truth approximation (large synthetic cohort):
#   1) Draw ONE baseline cohort (W1,W2,Y0,M1_0,M2_0) and keep it fixed.
#   2) For each counterfactual world, simulate longitudinal paths under fixed
#      regimens using the updated time indexing:
#
#         covariate time t = 0..T-1:  A_t, M1_t, M2_t, L_t
#         outcome time k = 1..T:     Y_k where Y_{t+1} uses covariate time t.
#
#      The estimand target is mean(Y_T).
#
# IMPORTANT implementation detail (separate-world mediators):
#   The "separate" (marginal) mediator world requires M1 and M2 to be drawn
#   independently from their *marginal* distributions under the specified
#   regimens, given the shared baseline. Therefore, even when regimen1 ==
#   regimen2 (e.g., a,a), we MUST generate M1 and M2 from two independent
#   simulation runs and then pair them.
################################################################################

truth_generate_baseline <- function(B_truth, params, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  W1 <- stats::rnorm(B_truth)
  W2 <- stats::rnorm(B_truth)

  # Y0
  bY0 <- params$baseline$Y0
  Y0 <- bY0$int + bY0$W1 * W1 + bY0$W2 * W2 + stats::rnorm(B_truth, 0, params$sigma$Y0)

  # Baseline mediators (M1_0, M2_0)
  bM10 <- params$baseline$M10
  M1_0 <- bM10$int + bM10$W1 * W1 + bM10$W2 * W2 + bM10$Y0 * Y0 +
    stats::rnorm(B_truth, 0, params$sigma$M10)

  bM20 <- params$baseline$M20
  M2_0 <- bM20$int + bM20$W1 * W1 + bM20$W2 * W2 + bM20$Y0 * Y0 +
    stats::rnorm(B_truth, 0, params$sigma$M20)

  data.frame(W1 = W1, W2 = W2, Y0 = Y0, M1_0 = M1_0, M2_0 = M2_0)
}

# Internal: correlated mediator errors with exposure-modified correlation.
.truth_gen_mediator_errors <- function(A, sigma1, sigma2, rho0, rho1) {
  n <- length(A)
  e1 <- numeric(n)
  e2 <- numeric(n)

  # A = 0 group
  idx0 <- which(A == 0)
  if (length(idx0) > 0) {
    r <- rho0
    z1 <- stats::rnorm(length(idx0))
    z2 <- stats::rnorm(length(idx0))
    e1[idx0] <- sigma1 * z1
    e2[idx0] <- sigma2 * (r * z1 + sqrt(pmax(1 - r^2, 0)) * z2)
  }

  # A = 1 group
  idx1 <- which(A == 1)
  if (length(idx1) > 0) {
    r <- rho1
    z1 <- stats::rnorm(length(idx1))
    z2 <- stats::rnorm(length(idx1))
    e1[idx1] <- sigma1 * z1
    e2[idx1] <- sigma2 * (r * z1 + sqrt(pmax(1 - r^2, 0)) * z2)
  }

  list(e1 = e1, e2 = e2)
}

# Simulate full paths under the DGP from a *fixed* regimen (A is not stochastic here).
truth_simulate_full_paths <- function(baseline, params, regimen) {
  n <- nrow(baseline)
  T <- length(regimen)

  W1 <- baseline$W1
  W2 <- baseline$W2
  Y0 <- baseline$Y0

  # Lags at outcome time k=1 (covariate time t=0):
  #   - L_{-1} is set to 0
  #   - Y_0 is baseline outcome
  #   - baseline mediators M1_0/M2_0 are treated as pre-treatment history
  M1_lag <- baseline$M1_0
  M2_lag <- baseline$M2_0
  L_lag  <- rep(0, n)
  Y_lag  <- Y0

  A_mat  <- matrix(NA_real_, n, T)
  M1_mat <- matrix(NA_real_, n, T)
  M2_mat <- matrix(NA_real_, n, T)
  L_mat  <- matrix(NA_real_, n, T)
  Y_mat  <- matrix(NA_real_, n, T)

  # Baseline mediators occupy the first column (covariate time 0)
  M1_mat[, 1L] <- baseline$M1_0
  M2_mat[, 1L] <- baseline$M2_0

  for (t in 1:T) {
    # A_{t-1}
    A_t <- rep(regimen[t], n)

    # Mediators at covariate time (t-1)
    # Baseline mediators (t=1) are pre-treatment and are not regenerated.
    if (t == 1L) {
      M1_t <- baseline$M1_0
      M2_t <- baseline$M2_0
    } else {
      b1 <- params$M1
      mu1 <- .get_coef_t(b1$int, t) +
        .get_coef_t(b1$A, t) * A_t +
        .get_coef_t(b1$W1, t) * W1 +
        .get_coef_t(b1$W2, t) * W2 +
        .get_coef_t(b1$L_lag, t) * L_lag +
        .get_coef_t(b1$Y_lag, t) * Y_lag +
        .get_coef_t(b1$M1_lag, t) * M1_lag +
        .get_coef_t(b1$M2_lag, t) * M2_lag

      b2 <- params$M2
      mu2 <- .get_coef_t(b2$int, t) +
        .get_coef_t(b2$A, t) * A_t +
        .get_coef_t(b2$W1, t) * W1 +
        .get_coef_t(b2$W2, t) * W2 +
        .get_coef_t(b2$L_lag, t) * L_lag +
        .get_coef_t(b2$Y_lag, t) * Y_lag +
        .get_coef_t(b2$M1_lag, t) * M1_lag +
        .get_coef_t(b2$M2_lag, t) * M2_lag

      errs <- .truth_gen_mediator_errors(
        A = A_t,
        sigma1 = params$sigma$M1,
        sigma2 = params$sigma$M2,
        rho0 = params$rho0,
        rho1 = params$rho1
      )
      M1_t <- mu1 + errs$e1
      M2_t <- mu2 + (params$gamma0 + params$gammaA * A_t) * M1_t + errs$e2
    }

    # L
    bL <- params$L
    muL <- .get_coef_t(bL$int, t) +
      .get_coef_t(bL$A, t) * A_t +
      .get_coef_t(bL$W1, t) * W1 +
      .get_coef_t(bL$W2, t) * W2 +
      .get_coef_t(bL$L_lag, t) * L_lag +
      .get_coef_t(bL$Y_lag, t) * Y_lag +
      .get_coef_t(bL$M1, t) * M1_t +
      .get_coef_t(bL$M2, t) * M2_t

    L_t <- muL + stats::rnorm(n, 0, params$sigma$L)

    # Y
    bY <- params$Y
    muY <- .get_coef_t(bY$int, t) +
      .get_coef_t(bY$W1, t) * W1 +
      .get_coef_t(bY$W2, t) * W2 +
      .get_coef_t(bY$Y_lag, t) * Y_lag +
      .get_coef_t(bY$A, t) * A_t +
      .get_coef_t(bY$M1, t) * M1_t +
      .get_coef_t(bY$M2, t) * M2_t +
      .get_coef_t(bY$L, t) * L_t +
      .get_coef_t(bY$A_M1, t) * (A_t * M1_t) +
      .get_coef_t(bY$A_M2, t) * (A_t * M2_t) +
      params$delta * (M1_t * M2_t)

    Y_t <- muY + stats::rnorm(n, 0, params$sigma$Y)

    # store + update
    A_mat[, t]  <- A_t
    M1_mat[, t] <- M1_t
    M2_mat[, t] <- M2_t
    L_mat[, t]  <- L_t
    Y_mat[, t]  <- Y_t

    L_lag  <- L_t
    Y_lag  <- Y_t
    M1_lag <- M1_t
    M2_lag <- M2_t
  }

  list(A = A_mat, M1 = M1_mat, M2 = M2_mat, L = L_mat, Y = Y_mat)
}

# Simulate L,Y under an outer regimen while overriding mediators with supplied paths.
truth_simulate_LY_given_mediators <- function(baseline, params, outer_regimen, M1_path, M2_path) {
  n <- nrow(baseline)
  T <- length(outer_regimen)
  if (!all(dim(M1_path) == c(n, T))) .stop("M1_path must be n x T")
  if (!all(dim(M2_path) == c(n, T))) .stop("M2_path must be n x T")

  W1 <- baseline$W1
  W2 <- baseline$W2
  Y0 <- baseline$Y0

  L_lag <- rep(0, n)
  Y_lag <- Y0

  L_mat <- matrix(NA_real_, n, T)
  Y_mat <- matrix(NA_real_, n, T)

  for (t in 1:T) {
    A_t  <- rep(outer_regimen[t], n)
    M1_t <- M1_path[, t]
    M2_t <- M2_path[, t]

    bL <- params$L
    muL <- .get_coef_t(bL$int, t) +
      .get_coef_t(bL$A, t) * A_t +
      .get_coef_t(bL$W1, t) * W1 +
      .get_coef_t(bL$W2, t) * W2 +
      .get_coef_t(bL$L_lag, t) * L_lag +
      .get_coef_t(bL$Y_lag, t) * Y_lag +
      .get_coef_t(bL$M1, t) * M1_t +
      .get_coef_t(bL$M2, t) * M2_t

    L_t <- muL + stats::rnorm(n, 0, params$sigma$L)

    bY <- params$Y
    muY <- .get_coef_t(bY$int, t) +
      .get_coef_t(bY$W1, t) * W1 +
      .get_coef_t(bY$W2, t) * W2 +
      .get_coef_t(bY$Y_lag, t) * Y_lag +
      .get_coef_t(bY$A, t) * A_t +
      .get_coef_t(bY$M1, t) * M1_t +
      .get_coef_t(bY$M2, t) * M2_t +
      .get_coef_t(bY$L, t) * L_t +
      .get_coef_t(bY$A_M1, t) * (A_t * M1_t) +
      .get_coef_t(bY$A_M2, t) * (A_t * M2_t) +
      params$delta * (M1_t * M2_t)

    Y_t <- muY + stats::rnorm(n, 0, params$sigma$Y)

    L_mat[, t] <- L_t
    Y_mat[, t] <- Y_t

    L_lag <- L_t
    Y_lag <- Y_t
  }

  list(L = L_mat, Y = Y_mat)
}

# Compute all 9 world means (truth) for one Monte Carlo batch.
.compute_truth_world_means_once <- function(params, T, reg_a, reg_as, B_truth = 200000, seed = 1) {
  reg_a  <- normalize_regimen(reg_a,  T, "reg_a")
  reg_as <- normalize_regimen(reg_as, T, "reg_as")

  # --------------------------------------------------------------------------
  # Shared baseline across all worlds ("common baseline cohort")
  # --------------------------------------------------------------------------
  baseline <- truth_generate_baseline(B_truth, params, seed = seed)

  # --------------------------------------------------------------------------
  # Natural worlds
  # --------------------------------------------------------------------------
  nat_a  <- truth_simulate_full_paths(baseline, params, reg_a)
  nat_as <- truth_simulate_full_paths(baseline, params, reg_as)

  mu_nat_a  <- mean(nat_a$Y[, T])
  mu_nat_as <- mean(nat_as$Y[, T])

  # --------------------------------------------------------------------------
  # Joint mediator intervention worlds
  #   - draw (M1,M2) jointly under regimen and then simulate (L,Y) under outer
  # --------------------------------------------------------------------------
  joint_med_aa   <- truth_simulate_full_paths(baseline, params, reg_a)
  joint_med_asas <- truth_simulate_full_paths(baseline, params, reg_as)

  joint_aa   <- truth_simulate_LY_given_mediators(
    baseline, params, outer_regimen = reg_a,
    M1_path = joint_med_aa$M1, M2_path = joint_med_aa$M2
  )
  joint_asas <- truth_simulate_LY_given_mediators(
    baseline, params, outer_regimen = reg_as,
    M1_path = joint_med_asas$M1, M2_path = joint_med_asas$M2
  )
  joint_aas  <- truth_simulate_LY_given_mediators(
    baseline, params, outer_regimen = reg_a,
    M1_path = joint_med_asas$M1, M2_path = joint_med_asas$M2
  )

  mu_joint_aa   <- mean(joint_aa$Y[, T])
  mu_joint_asas <- mean(joint_asas$Y[, T])
  mu_joint_aas  <- mean(joint_aas$Y[, T])

  # --------------------------------------------------------------------------
  # Separate mediator intervention worlds (marginal/separate)
  #   - draw M1 and M2 independently (even if regimen1 == regimen2)
  #   - then simulate (L,Y) under outer regimen using paired (M1,M2)
  # --------------------------------------------------------------------------
  # M1 marginal draws
  sep_M1_a   <- truth_simulate_full_paths(baseline, params, reg_a)
  sep_M1_as  <- truth_simulate_full_paths(baseline, params, reg_as)

  # M2 marginal draws (independent runs)
  sep_M2_a   <- truth_simulate_full_paths(baseline, params, reg_a)
  sep_M2_as  <- truth_simulate_full_paths(baseline, params, reg_as)

  # sep( a ; a, a )
  sep_aaa <- truth_simulate_LY_given_mediators(
    baseline, params, outer_regimen = reg_a,
    M1_path = sep_M1_a$M1,
    M2_path = sep_M2_a$M2
  )

  # sep( a* ; a*, a* )
  sep_asas_asas <- truth_simulate_LY_given_mediators(
    baseline, params, outer_regimen = reg_as,
    M1_path = sep_M1_as$M1,
    M2_path = sep_M2_as$M2
  )

  # sep( a ; a*, a* )
  sep_a_asas <- truth_simulate_LY_given_mediators(
    baseline, params, outer_regimen = reg_a,
    M1_path = sep_M1_as$M1,
    M2_path = sep_M2_as$M2
  )

  # sep( a ; a, a* )
  sep_a_aas <- truth_simulate_LY_given_mediators(
    baseline, params, outer_regimen = reg_a,
    M1_path = sep_M1_a$M1,
    M2_path = sep_M2_as$M2
  )

  means <- list(
    mu_nat_a = mu_nat_a,
    mu_nat_as = mu_nat_as,

    mu_joint_aa = mu_joint_aa,
    mu_joint_asas = mu_joint_asas,
    mu_joint_aas = mu_joint_aas,

    mu_sep_aaa = mean(sep_aaa$Y[, T]),
    mu_sep_asas_asas = mean(sep_asas_asas$Y[, T]),
    mu_sep_a_asas = mean(sep_a_asas$Y[, T]),
    mu_sep_a_aas = mean(sep_a_aas$Y[, T])
  )

  effects <- compute_estimands_from_means(means)

  list(means = means, effects = effects)
}

.truth_batch_sizes <- function(B_truth, n_batches) {
  B_truth <- as.integer(B_truth)
  n_batches <- as.integer(n_batches)
  if (!is.finite(B_truth) || B_truth < 1L) .stop("B_truth must be a positive integer.")
  if (!is.finite(n_batches) || n_batches < 1L) .stop("n_batches must be a positive integer.")
  n_batches <- min(n_batches, B_truth)
  base <- B_truth %/% n_batches
  rem <- B_truth %% n_batches
  sizes <- rep(base, n_batches)
  if (rem > 0L) sizes[seq_len(rem)] <- sizes[seq_len(rem)] + 1L
  sizes
}

.truth_mcse <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

.truth_nested_mcse <- function(batch_effects_df) {
  out <- list()
  for (est in unique(batch_effects_df$estimand)) {
    ee <- batch_effects_df[batch_effects_df$estimand == est, , drop = FALSE]
    vals <- stats::aggregate(
      estimate ~ effect,
      data = ee,
      FUN = .truth_mcse
    )
    out[[est]] <- as.list(setNames(vals$estimate, vals$effect))
  }
  out
}

# Compute all 9 world means (truth) using large Monte Carlo split into batches.
compute_truth_world_means <- function(params, T, reg_a, reg_as,
                                      B_truth = 200000,
                                      seed = 1,
                                      n_batches = 20L) {
  batch_sizes <- .truth_batch_sizes(B_truth, n_batches)

  batch_results <- vector("list", length(batch_sizes))
  for (bb in seq_along(batch_sizes)) {
    batch_results[[bb]] <- .compute_truth_world_means_once(
      params = params,
      T = T,
      reg_a = reg_a,
      reg_as = reg_as,
      B_truth = batch_sizes[[bb]],
      seed = as.integer(seed)[1L] + 7919L * bb
    )
  }

  mean_keys <- component_mean_keys()
  batch_means <- do.call(rbind, lapply(seq_along(batch_results), function(bb) {
    vals <- unlist(batch_results[[bb]]$means[mean_keys])
    data.frame(
      batch = bb,
      B_batch = batch_sizes[[bb]],
      t(vals),
      check.names = FALSE,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }))

  means_mean <- as.list(setNames(
    vapply(mean_keys, function(nm) {
      stats::weighted.mean(batch_means[[nm]], w = batch_means$B_batch)
    }, numeric(1)),
    mean_keys
  ))
  means_mcse <- as.list(setNames(
    vapply(mean_keys, function(nm) .truth_mcse(batch_means[[nm]]), numeric(1)),
    mean_keys
  ))

  batch_effects <- do.call(rbind, lapply(seq_along(batch_results), function(bb) {
    eff <- flatten_effects(batch_results[[bb]]$effects)
    data.frame(
      batch = bb,
      B_batch = batch_sizes[[bb]],
      estimand = eff$estimand,
      effect = eff$effect,
      estimate = eff$estimate,
      stringsAsFactors = FALSE
    )
  }))

  effects_mean <- compute_estimands_from_means(means_mean)
  effects_mcse <- .truth_nested_mcse(batch_effects)

  list(
    means = means_mean,
    effects = effects_mean,
    means_mcse = means_mcse,
    effects_mcse = effects_mcse,
    batch_means = batch_means,
    batch_effects = batch_effects
  )
}
