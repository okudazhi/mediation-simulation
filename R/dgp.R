################################################################################
# dgp.R
#
# Data-generating process (DGP) for longitudinal mediation with two time-varying
# continuous mediators.
#
# Package-specific time indexing:
#
#   Pre-treatment baseline variables:
#     W1, W2, Y0, M1_0, M2_0
#
#   Baseline mediators M1_0 and M2_0 are generated before A_0 and are treated as
#   pre-treatment history throughout the package.
#
#   At baseline visit (covariate time 0):
#     M1_0, M2_0  ->  A_0  ->  L_0  ->  Y_1
#
#   At follow-up visits t = 1,...,T-1:
#     A_t  ->  M1_t  ->  M2_t  ->  L_t  ->  Y_{t+1}
#
# Therefore the final outcome is Y_T, generated at covariate time t=T-1.
#
# Treatment assignment mechanisms (treat_mech):
#   - "baseline_rct"   : one baseline randomization A_0 carried forward
#   - "sequential_rct" : independent randomization for each A_t
#   - "observational"  : logistic model in observed history (user-specified)
#   - "fixed"          : deterministic regimen provided by the user
#
# NOTE: Coefficient values below are defaults only. All coefficients are editable
# at run time.
################################################################################


# ---- Default parameter list --------------------------------------------------

#' Create a default DGP parameter list.
#'
#' All coefficients are scalar (time-homogeneous) by default, but you may
#' replace any scalar with a length-T vector if you want time-varying effects.
#'
#' The defaults are intentionally moderate so that the simulation runs without
#' numerical instabilities. Please adapt them to your scientific calibration.
#'
#' @param T number of follow-up times
#' @return list of parameters
default_dgp_params <- function(T) {
  list(
    # --- baseline distributions ---
    sigma = list(
      Y0 = 1.0,
      M10 = 1.0,
      M20 = 1.0,
      M1 = 1.0,
      M2 = 1.0,
      L  = 1.0,
      Y  = 1.0
    ),

    baseline = list(
      # Y0
      Y0 = list(int = 0.0, W1 = 0.5, W2 = 0.5),
      # M1_0 and M2_0 depend on baseline W and Y0
      M10 = list(int = 0.0, W1 = 0.4, W2 = 0.4, Y0 = 0.3),
      M20 = list(int = 0.0, W1 = 0.4, W2 = 0.4, Y0 = 0.3)
    ),

    # --- treatment assignment ---
    treatment = list(
      p_rct = 0.5,
      # Observational A_t model: logit P(A_t=1|H_t) = int + ...
      # (Only used when treat_mech == "observational")
      obs = list(
        int = 0.0,
        W1 = 0.2,
        W2 = -0.2,
        Y_lag = 0.1,
        L_lag = 0.1,
        M1_lag = 0.1,
        M2_lag = 0.1
      )
    ),

    # --- mediator dependence switches ---
    # M2_t includes: (gamma0 + gammaA*A_t) * M1_t
    gamma0 = 0.4,
    gammaA = 0.2,

    # Corr(e1,e2 | A_t) = rho0 + (rho1-rho0)*A_t
    rho0 = 0.0,
    rho1 = 0.0,

    # Outcome includes delta * (M1_t*M2_t)
    delta = 0.05,

    # --- time-varying node models (default: time-homogeneous) ---
    M1 = list(
      int = 0.0,
      A = 0.4,
      W1 = 0.2,
      W2 = 0.2,
      L_lag = 0.2,
      Y_lag = 0.2,
      M1_lag = 0.3,
      M2_lag = 0
    ),

    # mu2 excludes the (gamma0 + gammaA*A)*M1 term, which is added separately.
    M2 = list(
      int = 0.0,
      A = 0.3,
      W1 = 0.2,
      W2 = 0.2,
      L_lag = 0.2,
      Y_lag = 0.2,
      M1_lag = 0,
      M2_lag = 0.3
    ),

    L = list(
      int = 0.0,
      A = 0.2,
      W1 = 0.2,
      W2 = 0.2,
      L_lag = 0.3,
      Y_lag = 0.2,
      M1 = 0.2,
      M2 = 0.2
    ),

    Y = list(
      int = 0.0,
      W1 = 0.2,
      W2 = 0.2,
      Y_lag = 0.3,
      A = 0.2,
      M1 = 0.1,
      M2 = 0.1,
      L  = 0.1,
      A_M1 = 0.1,
      A_M2 = 0.1
      # delta is global
    )
  )
}

# ---- Internal helpers --------------------------------------------------------

.get_coef_t <- function(x, t) {
  # Allow scalar or length-T vector.
  if (length(x) == 1L) return(as.numeric(x))
  as.numeric(x[[t]])
}

# Generate correlated mediator errors e1,e2 with correlation depending on A.
# Vectorized by splitting A==0 and A==1.
.gen_mediator_errors <- function(A, sigma1, sigma2, rho0, rho1) {
  n <- length(A)
  e1 <- numeric(n)
  e2 <- numeric(n)

  # For A=0
  idx0 <- which(A == 0)
  if (length(idx0) > 0) {
    rho <- rho0
    z1 <- stats::rnorm(length(idx0))
    z2 <- stats::rnorm(length(idx0))
    e1[idx0] <- sigma1 * z1
    e2[idx0] <- sigma2 * (rho * z1 + sqrt(pmax(1 - rho^2, 0)) * z2)
  }

  # For A=1
  idx1 <- which(A == 1)
  if (length(idx1) > 0) {
    rho <- rho1
    z1 <- stats::rnorm(length(idx1))
    z2 <- stats::rnorm(length(idx1))
    e1[idx1] <- sigma1 * z1
    e2[idx1] <- sigma2 * (rho * z1 + sqrt(pmax(1 - rho^2, 0)) * z2)
  }

  list(e1 = e1, e2 = e2)
}

.logit <- function(x) 1 / (1 + exp(-x))

# Assign A_t under each treatment mechanism.
.assign_A_t <- function(t, n, treat_mech, params, W1, W2,
                        L_lag, Y_lag, M1_lag, M2_lag,
                        A_baseline, fixed_regimen, A_assign_fn) {
  if (!is.null(A_assign_fn)) {
    # User-provided function takes a data.frame of current history.
    hist <- data.frame(
      t = t,
      W1 = W1, W2 = W2,
      L_lag = L_lag, Y_lag = Y_lag,
      M1_lag = M1_lag, M2_lag = M2_lag,
      A_baseline = A_baseline
    )
    A <- A_assign_fn(hist, params)
    return(as.numeric(A))
  }

  if (treat_mech == "baseline_rct") {
    # A_baseline is generated once per subject (typically length n) and carried
    # forward across all follow-up times.
    #
    # IMPORTANT: do NOT use rep(A_baseline, n), which would create a length n*n
    # vector and can trigger recycling warnings/errors when assigning into an
    # n-length vector.
    if (length(A_baseline) == 1L) {
      return(rep(as.numeric(A_baseline), n))
    }
    if (length(A_baseline) != n) {
      .stop("A_baseline must have length 1 or n. Got length=", length(A_baseline), ", n=", n, ".")
    }
    return(as.numeric(A_baseline))
  }

  if (treat_mech == "sequential_rct") {
    p <- params$treatment$p_rct
    return(stats::rbinom(n, 1, p))
  }

  if (treat_mech == "observational") {
    b <- params$treatment$obs
    eta <- b$int + b$W1 * W1 + b$W2 * W2 + b$Y_lag * Y_lag + b$L_lag * L_lag +
      b$M1_lag * M1_lag + b$M2_lag * M2_lag
    p <- .logit(eta)
    return(stats::rbinom(n, 1, p))
  }

  if (treat_mech == "fixed") {
    if (is.null(fixed_regimen)) .stop("fixed_regimen must be supplied when treat_mech='fixed'.")
    return(rep(fixed_regimen[t], n))
  }

  .stop("Unknown treat_mech: ", treat_mech)
}

# ---- Main DGP generator (wide) ----------------------------------------------

#' Simulate a dataset under the DGP in *wide* format.
#'
#' @param n sample size
#' @param T number of follow-up times
#' @param params DGP parameter list (see default_dgp_params)
#' @param treat_mech treatment mechanism (see header)
#' @param fixed_regimen regimen vector length T (only if treat_mech == 'fixed')
#' @param A_assign_fn optional custom function for treatment assignment
#' @param seed optional seed
#' @return wide data.frame with baseline and time-varying columns
simulate_dgp_wide <- function(n, T, params,
                              treat_mech = c("baseline_rct","sequential_rct","observational","fixed"),
                              fixed_regimen = NULL,
                              A_assign_fn = NULL,
                              seed = NULL) {
  treat_mech <- match.arg(treat_mech)
  if (!is.null(seed)) set.seed(seed)

  # --- Baseline ---
  W1 <- stats::rnorm(n)
  W2 <- stats::rnorm(n)

  # Y0
  bY0 <- params$baseline$Y0
  Y0 <- bY0$int + bY0$W1 * W1 + bY0$W2 * W2 + stats::rnorm(n, 0, params$sigma$Y0)

  # M1_0, M2_0 (pre-treatment baseline mediators; generated before A_0)
  bM10 <- params$baseline$M10
  M1_0 <- bM10$int + bM10$W1 * W1 + bM10$W2 * W2 + bM10$Y0 * Y0 + stats::rnorm(n, 0, params$sigma$M10)

  bM20 <- params$baseline$M20
  M2_0 <- bM20$int + bM20$W1 * W1 + bM20$W2 * W2 + bM20$Y0 * Y0 + stats::rnorm(n, 0, params$sigma$M20)

  # --- Initialize storage ---
  # Indexing:
  #   - Matrix column t=1..T corresponds to covariate time (t-1)
  #   - A_mat[,t] stores A_{t-1}
  #   - M*_mat[,t] stores M*_{t-1}
  #   - L_mat[,t] stores L_{t-1}
  #   - Y_mat[,t] stores Y_t (generated from covariate time t-1)
  A_mat  <- matrix(NA_real_, n, T)
  M1_mat <- matrix(NA_real_, n, T)
  M2_mat <- matrix(NA_real_, n, T)
  L_mat  <- matrix(NA_real_, n, T)
  Y_mat  <- matrix(NA_real_, n, T)

  # Baseline mediators at covariate time 0
  M1_mat[, 1L] <- M1_0
  M2_mat[, 1L] <- M2_0

  # Lags used for generating time t-1 nodes:
  #   - At t=1 (covariate time 0): use L_{-1}=0 and Y_0=Y0.
  L_lag  <- rep(0, n)  # L_{-1}
  Y_lag  <- Y0         # Y_0
  M1_lag <- M1_0       # treat as baseline history for A_0
  M2_lag <- M2_0

  # baseline A_0 for baseline_rct (carried forward)
  A0 <- stats::rbinom(n, 1, params$treatment$p_rct)

  # fixed regimen sanity
  if (treat_mech == "fixed") {
    fixed_regimen <- normalize_regimen(fixed_regimen, T, "fixed_regimen")
  }

  for (t in 1:T) {
    # --- Treatment assignment (A_{t-1}) ---
    A_t <- .assign_A_t(
      t, n, treat_mech, params,
      W1, W2,
      L_lag, Y_lag, M1_lag, M2_lag,
      A_baseline = A0,
      fixed_regimen = fixed_regimen,
      A_assign_fn = A_assign_fn
    )
    A_t <- as.numeric(A_t)

    # --- Mediators (M*_{t-1}) ---
    # Baseline mediators (t=1, covariate time 0) are generated *before* A_0
    # and are treated as fixed history here (do not regenerate).
    if (t == 1L) {
      M1_t <- M1_0
      M2_t <- M2_0
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

      errs <- .gen_mediator_errors(
        A = A_t,
        sigma1 = params$sigma$M1,
        sigma2 = params$sigma$M2,
        rho0 = params$rho0,
        rho1 = params$rho1
      )

      M1_t <- mu1 + errs$e1
      M2_t <- mu2 + (params$gamma0 + params$gammaA * A_t) * M1_t + errs$e2
    }

    # --- Post-mediator covariate L ---
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

    # --- Outcome (Y_t = Y_{(t-1)+1}) ---
    bY <- params$Y
    muY <- .get_coef_t(bY$int, t) +
      .get_coef_t(bY$W1, t) * W1 +
      .get_coef_t(bY$W2, t) * W2 +
      .get_coef_t(bY$Y_lag, t) * Y_lag +
      .get_coef_t(bY$A, t) * A_t +
      .get_coef_t(bY$M1, t) * M1_t +
      .get_coef_t(bY$M2, t) * M2_t +
      .get_coef_t(bY$L, t)  * L_t +
      .get_coef_t(bY$A_M1, t) * (A_t * M1_t) +
      .get_coef_t(bY$A_M2, t) * (A_t * M2_t) +
      params$delta * (M1_t * M2_t)

    Y_t <- muY + stats::rnorm(n, 0, params$sigma$Y)

    # store
    A_mat[, t]  <- A_t
    M1_mat[, t] <- M1_t
    M2_mat[, t] <- M2_t
    L_mat[, t]  <- L_t
    Y_mat[, t]  <- Y_t

    # update lags
    L_lag  <- L_t
    Y_lag  <- Y_t
    M1_lag <- M1_t
    M2_lag <- M2_t
  }

  # wide data.frame
  dat <- data.frame(
    id = seq_len(n),
    W1 = W1,
    W2 = W2,
    Y0 = Y0,
    M1_0 = M1_0,
    M2_0 = M2_0
  )

  for (t in 1:T) {
    # Time-varying nodes at covariate time (t-1)
    dat[[paste0("A", t - 1L)]] <- A_mat[, t]
    dat[[paste0("L", t - 1L)]] <- L_mat[, t]

    # Mediators: baseline columns M1_0/M2_0 already exist; follow-ups start at t=2.
    if (t >= 2L) {
      dat[[paste0("M1_", t - 1L)]] <- M1_mat[, t]
      dat[[paste0("M2_", t - 1L)]] <- M2_mat[, t]
    }

    # Outcomes at time t
    dat[[paste0("Y_", t)]] <- Y_mat[, t]
  }

  dat
}

# ---- Scenario adjustment -----------------------------------------------------

#' Apply scenario switches to a parameter list.
#'
#' Scenario components (manuscript-concordant):
#'   - PM: pathway magnitude via coefficient scaling
#'       * c1: scales M1$A
#'       * c2: scales M2$A
#'       * c3: scales Y$M1 and Y$A_M1
#'       * c4: scales Y$M2, Y$A_M2, and params$delta
#'   - MI_mode: mediator interdependence structure (distinct from residual rho)
#'       * "none"        : gamma0=0, gammaA=0, delta=0
#'       * "gamma0_only" : gammaA=0, delta=0 (gamma0 kept at base value)
#'       * "full"        : keep (gamma0, gammaA, delta)
#'   - rho0/rho1: residual dependence Corr(e1,e2 | A_t=0/1)
#'
#' @param params base parameter list
#' @param PM optional list with multipliers c1,c2,c3,c4
#' @param MI_mode character in {"none","gamma0_only","full"}
#' @param rho0,rho1 residual correlations (can be scalars)
#' @return modified params
apply_scenario_to_params <- function(params,
                                     PM = NULL,
                                     MI_mode = c("full", "gamma0_only", "none"),
                                     rho0 = NULL,
                                     rho1 = NULL) {
  p <- params

  if (!is.null(PM)) {
    if (!is.list(PM)) .stop("PM must be a list with keys c1,c2,c3,c4")
    req <- c("c1", "c2", "c3", "c4")
    miss <- setdiff(req, names(PM))
    if (length(miss) > 0) .stop("PM missing keys: ", paste(miss, collapse = ", "))

    p$M1$A   <- p$M1$A   * PM$c1
    p$M2$A   <- p$M2$A   * PM$c2
    p$Y$M1   <- p$Y$M1   * PM$c3
    p$Y$A_M1 <- p$Y$A_M1 * PM$c3
    p$Y$M2   <- p$Y$M2   * PM$c4
    p$Y$A_M2 <- p$Y$A_M2 * PM$c4
    p$delta  <- p$delta  * PM$c4
  }

  MI_mode <- match.arg(MI_mode)
  if (MI_mode == "none") {
    p$gamma0 <- 0
    p$gammaA <- 0
    p$delta  <- 0
  } else if (MI_mode == "gamma0_only") {
    p$gammaA <- 0
    p$delta  <- 0
  }

  if (!is.null(rho0)) p$rho0 <- rho0
  if (!is.null(rho1)) p$rho1 <- rho1

  p
}
