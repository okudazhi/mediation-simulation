################################################################################
# core_utils.R
#
# Shared utilities for the dual-use package:
#   - simulation helpers that are safe to reuse
#   - real-data analysis helpers
#
# The functions here are intentionally dependency-light so that the package can
# be sourced as plain R scripts without installation.
################################################################################

.stop <- function(..., call. = FALSE) stop(paste0(...), call. = call.)

assert_cols <- function(dat, cols, dat_name = "data") {
  miss <- setdiff(cols, names(dat))
  if (length(miss) > 0L) {
    .stop(dat_name, " is missing required columns: ", paste(miss, collapse = ", "))
  }
  invisible(TRUE)
}

assert_no_na <- function(dat, cols, dat_name = "data") {
  bad <- cols[vapply(cols, function(v) any(is.na(dat[[v]])), logical(1))]
  if (length(bad) > 0L) {
    .stop(dat_name, " has NA in required columns: ", paste(bad, collapse = ", "))
  }
  invisible(TRUE)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

.warn_once <- local({
  seen <- new.env(parent = emptyenv())
  function(key, ..., call. = FALSE) {
    if (exists(key, envir = seen, inherits = FALSE)) return(invisible(FALSE))
    assign(key, TRUE, envir = seen)
    warning(..., call. = call.)
    invisible(TRUE)
  }
})

normalize_regimen <- function(reg, T, name = "regimen") {
  if (length(reg) == 1L) return(rep(as.numeric(reg), T))
  if (length(reg) != T) {
    .stop(name, " must have length 1 or length T. Got length=", length(reg), ", T=", T)
  }
  as.numeric(reg)
}

validate_baseline_rct_regimens <- function(reg_a, reg_as, treat_mech) {
  if (!identical(treat_mech, "baseline_rct")) return(invisible(TRUE))
  if (any(!is.finite(reg_a)) || any(reg_a != reg_a[1L])) {
    .stop("baseline_rct requires reg_a to be constant across visits.")
  }
  if (any(!is.finite(reg_as)) || any(reg_as != reg_as[1L])) {
    .stop("baseline_rct requires reg_as to be constant across visits.")
  }
  invisible(TRUE)
}

make_switching_regimen <- function(primary, comparator, t_switch) {
  primary <- as.numeric(primary)
  comparator <- as.numeric(comparator)
  if (length(primary) != length(comparator)) {
    .stop("primary and comparator must have the same length.")
  }
  T <- length(primary)
  if (!is.numeric(t_switch) || length(t_switch) != 1L ||
      t_switch < 0L || t_switch > (T - 1L)) {
    .stop("t_switch must be an integer in {0, ..., T-1}.")
  }
  out <- primary
  if (t_switch < (T - 1L)) {
    out[(t_switch + 2L):T] <- comparator[(t_switch + 2L):T]
  }
  out
}

safe_div <- function(num, den, eps = 1e-12) {
  if (!is.finite(den) || abs(den) < eps) return(NA_real_)
  num / den
}

clamp01 <- function(p, eps = 1e-8) {
  pmin(1 - eps, pmax(eps, as.numeric(p)))
}

truncate_vec <- function(x, p_lo = 0.01, p_hi = 0.99) {
  if (!length(x) || all(is.na(x))) return(x)
  qs <- stats::quantile(x, probs = c(p_lo, p_hi), na.rm = TRUE, type = 7)
  pmin(pmax(x, qs[[1L]]), qs[[2L]])
}

truncate_by_group <- function(df, wcol, gcols, p_lo = 0.01, p_hi = 0.99) {
  assert_cols(df, c(wcol, gcols), "truncate_by_group(df)")
  g <- interaction(df[gcols], drop = TRUE, lex.order = TRUE)
  out <- df[[wcol]]
  levs <- levels(g)
  for (lv in levs) {
    idx <- which(g == lv)
    wi <- out[idx]
    wi <- wi[is.finite(wi)]
    if (!length(wi)) next
    qs <- stats::quantile(wi, probs = c(p_lo, p_hi), na.rm = TRUE, type = 7)
    out[idx] <- ifelse(is.na(out[idx]), NA_real_,
                       pmin(pmax(out[idx], qs[[1L]]), qs[[2L]]))
  }
  df[[wcol]] <- out
  df
}

safe_quantile <- function(x, probs) {
  x <- x[is.finite(x)]
  if (!length(x)) return(rep(NA_real_, length(probs)))
  as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 7))
}

ess_from_weights <- function(w) {
  w <- as.numeric(w)
  w <- w[is.finite(w)]
  n <- length(w)
  if (!n) return(c(ess = NA_real_, ess_frac = NA_real_))
  sw <- sum(w)
  sw2 <- sum(w * w)
  if (!is.finite(sw) || !is.finite(sw2) || sw2 <= 0) {
    return(c(ess = NA_real_, ess_frac = NA_real_))
  }
  ess <- (sw * sw) / sw2
  c(ess = ess, ess_frac = ess / n)
}

cumprod_by_id <- function(x, id) {
  out <- as.numeric(x)
  split_idx <- split(seq_along(out), id)
  for (ii in split_idx) out[ii] <- cumprod(out[ii])
  out
}

lag_by_id <- function(x, id, default = 0) {
  out <- x
  split_idx <- split(seq_along(x), id)
  for (ii in split_idx) {
    if (length(ii) == 1L) {
      out[ii] <- default
    } else {
      out[ii] <- c(default, x[ii][-length(ii)])
    }
  }
  out
}

clean_dot_na <- function(dat, na_strings = c(".", "NA", "")) {
  out <- dat
  for (nm in names(out)) {
    if (is.factor(out[[nm]])) out[[nm]] <- as.character(out[[nm]])
    if (is.character(out[[nm]])) {
      out[[nm]][trimws(out[[nm]]) %in% na_strings] <- NA_character_
    }
  }
  out
}

coerce_numeric_cols <- function(dat, cols) {
  out <- dat
  cols <- intersect(unique(cols), names(out))
  for (nm in cols) {
    if (is.factor(out[[nm]])) out[[nm]] <- as.character(out[[nm]])
    out[[nm]] <- suppressWarnings(as.numeric(out[[nm]]))
  }
  out
}

default_varmap <- function(T) {
  if (!is.numeric(T) || length(T) != 1L || T < 1L) {
    .stop("T must be a positive integer. Got: ", T)
  }
  T <- as.integer(T)
  .seq_if <- function(from, to) if (from <= to) from:to else integer(0)

  list(
    W1 = "W1",
    W2 = "W2",
    Y0 = "Y0",
    M10 = "M1_0",
    M20 = "M2_0",
    A = paste0("A", 0:(T - 1L)),
    M1 = c("M1_0", paste0("M1_", .seq_if(1L, T - 1L))),
    M2 = c("M2_0", paste0("M2_", .seq_if(1L, T - 1L))),
    L = paste0("L", 0:(T - 1L)),
    Y = paste0("Y_", 1:T)
  )
}

validate_wide_data <- function(dat, T, varmap = default_varmap(T)) {
  req <- unique(c(
    varmap$W1, varmap$W2, varmap$Y0, varmap$M10, varmap$M20,
    varmap$A, varmap$M1, varmap$M2, varmap$L, varmap$Y
  ))
  assert_cols(dat, req, "wide data")
  assert_no_na(dat, req, "wide data")
  invisible(TRUE)
}

extract_baseline <- function(dat, varmap = default_varmap(1L), dat_name = "dat (wide)") {
  req <- c(varmap$W1, varmap$W2, varmap$Y0, varmap$M10, varmap$M20)
  assert_cols(dat, req, dat_name)
  data.frame(
    W1 = as.numeric(dat[[varmap$W1]]),
    W2 = as.numeric(dat[[varmap$W2]]),
    Y0 = as.numeric(dat[[varmap$Y0]]),
    M1_0 = as.numeric(dat[[varmap$M10]]),
    M2_0 = as.numeric(dat[[varmap$M20]]),
    stringsAsFactors = FALSE
  )
}

wide_to_long <- function(dat, T, varmap = default_varmap(T), id_name = "id") {
  validate_wide_data(dat, T, varmap)

  n <- nrow(dat)
  id <- seq_len(n)
  out <- vector("list", T)

  for (t in 1:T) {
    A_tm1 <- as.numeric(dat[[varmap$A[t]]])
    M1_tm1 <- as.numeric(dat[[varmap$M1[t]]])
    M2_tm1 <- as.numeric(dat[[varmap$M2[t]]])
    L_tm1 <- as.numeric(dat[[varmap$L[t]]])
    Y_t <- as.numeric(dat[[varmap$Y[t]]])

    L_lag <- if (t == 1L) rep(0, n) else as.numeric(dat[[varmap$L[t - 1L]]])
    Y_lag <- if (t == 1L) as.numeric(dat[[varmap$Y0]]) else as.numeric(dat[[varmap$Y[t - 1L]]])
    M1_lag <- if (t == 1L) as.numeric(dat[[varmap$M10]]) else as.numeric(dat[[varmap$M1[t - 1L]]])
    M2_lag <- if (t == 1L) as.numeric(dat[[varmap$M20]]) else as.numeric(dat[[varmap$M2[t - 1L]]])

    out[[t]] <- data.frame(
      id = id,
      t = t,
      W1 = as.numeric(dat[[varmap$W1]]),
      W2 = as.numeric(dat[[varmap$W2]]),
      Y0 = as.numeric(dat[[varmap$Y0]]),
      M1_0 = as.numeric(dat[[varmap$M10]]),
      M2_0 = as.numeric(dat[[varmap$M20]]),
      A = A_tm1,
      M1 = M1_tm1,
      M2 = M2_tm1,
      L = L_tm1,
      Y = Y_t,
      L_lag = L_lag,
      Y_lag = Y_lag,
      M1_lag = M1_lag,
      M2_lag = M2_lag
    )
  }

  long <- do.call(rbind, out)
  names(long)[names(long) == "id"] <- id_name
  long
}

wide_to_long_with_hist <- function(dat, T, varmap = default_varmap(T), id_name = "id") {
  validate_wide_data(dat, T, varmap)

  long <- wide_to_long(dat, T, varmap, id_name = id_name)
  id <- id_name
  hist_cols <- unique(c(varmap$A, varmap$M1, varmap$M2, varmap$L, varmap$Y))
  hist_cols <- setdiff(hist_cols, names(long))
  if (!length(hist_cols)) return(long)

  wide_hist <- dat[, hist_cols, drop = FALSE]
  wide_hist[[id]] <- seq_len(nrow(dat))
  idx <- match(long[[id]], wide_hist[[id]])
  cbind(long, wide_hist[idx, hist_cols, drop = FALSE])
}

project_timestamp <- function() {
  format(Sys.time(), "%Y%m%d_%H%M%S")
}

safe_predict_lm <- function(fit, newdata, fallback = NA_real_) {
  if (is.list(fit) && !is.null(fit$type) && fit$type %in% c("glm", "superlearner")) {
    return(predict_nuisance_model(fit, newdata, type = "numeric", fallback = fallback))
  }
  n <- nrow(newdata)
  if (!n) return(numeric(0))
  out <- rep(NA_real_, n)

  vars <- tryCatch(all.vars(stats::delete.response(stats::terms(fit))),
                   error = function(e) character(0))
  vars <- intersect(vars, names(newdata))
  keep <- if (length(vars)) stats::complete.cases(newdata[, vars, drop = FALSE]) else rep(TRUE, n)

  if (any(keep)) {
    pred <- tryCatch(
      stats::predict(fit, newdata = newdata[keep, , drop = FALSE]),
      error = function(e) rep(NA_real_, sum(keep))
    )
    pred <- as.numeric(pred)
    if (length(pred) == sum(keep)) out[keep] <- pred
  }

  bad <- !is.finite(out)
  if (any(bad)) {
    if (length(fallback) == 1L) {
      out[bad] <- as.numeric(fallback)[1L]
    } else if (length(fallback) == n) {
      out[bad] <- as.numeric(fallback)[bad]
    }
  }
  out
}

safe_glm_prob <- function(fit, newdata, eps = 1e-8, fallback = 0.5) {
  if (is.list(fit) && !is.null(fit$type) && fit$type %in% c("glm", "superlearner")) {
    return(predict_nuisance_model(fit, newdata, type = "response", eps = eps, fallback = fallback))
  }
  n <- nrow(newdata)
  if (!n) return(numeric(0))
  out <- rep(NA_real_, n)

  vars <- tryCatch(all.vars(stats::delete.response(stats::terms(fit))),
                   error = function(e) character(0))
  vars <- intersect(vars, names(newdata))
  keep <- if (length(vars)) stats::complete.cases(newdata[, vars, drop = FALSE]) else rep(TRUE, n)

  if (any(keep)) {
    pred <- tryCatch(
      stats::predict(fit, newdata = newdata[keep, , drop = FALSE], type = "response"),
      error = function(e) rep(NA_real_, sum(keep))
    )
    pred <- as.numeric(pred)
    if (length(pred) == sum(keep)) out[keep] <- pred
  }

  bad <- !is.finite(out)
  if (any(bad)) {
    if (length(fallback) == 1L) {
      out[bad] <- as.numeric(fallback)[1L]
    } else if (length(fallback) == n) {
      out[bad] <- as.numeric(fallback)[bad]
    } else {
      out[bad] <- as.numeric(fallback)[1L]
    }
  }
  clamp01(out, eps = eps)
}

fit_gaussian_lm_safe <- function(formula, data, weights = NULL) {
  if (!nrow(data)) return(NULL)
  keep <- stats::complete.cases(model.frame(formula, data = data, na.action = stats::na.pass))
  if (!any(keep)) return(NULL)
  dd <- data[keep, , drop = FALSE]
  ww <- NULL
  if (!is.null(weights)) {
    ww <- as.numeric(weights)[keep]
    ww[!is.finite(ww) | ww < 0] <- 0
    if (all(ww == 0)) ww <- NULL
  }
  tryCatch({
    stats::lm(formula, data = dd, weights = ww, model = FALSE, x = FALSE, y = FALSE)
  }, error = function(e) NULL)
}

fit_binomial_glm_safe <- function(formula, data, weights = NULL) {
  if (!nrow(data)) return(NULL)
  keep <- stats::complete.cases(model.frame(formula, data = data, na.action = stats::na.pass))
  if (!any(keep)) return(NULL)
  dd <- data[keep, , drop = FALSE]
  ww <- NULL
  if (!is.null(weights)) {
    ww <- as.numeric(weights)[keep]
    ww[!is.finite(ww) | ww < 0] <- 0
    if (all(ww == 0)) ww <- NULL
  }
  tryCatch({
    stats::glm(formula, data = dd, family = stats::binomial(), weights = ww,
               model = FALSE, x = FALSE, y = FALSE)
  }, error = function(e) NULL)
}

fit_gaussian_sl_or_glm <- function(formula, data, weights = NULL,
                                   learner = c("glm", "sl"),
                                   sl_library = c("SL.glm", "SL.mean")) {
  learner <- match.arg(learner)
  if (identical(learner, "glm")) {
    return(list(type = "glm", fit = fit_gaussian_lm_safe(formula, data, weights = weights)))
  }

  if (!requireNamespace("SuperLearner", quietly = TRUE)) {
    .warn_once("ltmle_exact_sl_missing_superlearner",
               "Package 'SuperLearner' is not installed; falling back to glm for the LTMLE nuisance backend.",
               call. = FALSE)
    return(list(type = "glm", fit = fit_gaussian_lm_safe(formula, data, weights = weights)))
  }

  mf <- tryCatch(model.frame(formula, data = data, na.action = stats::na.pass), error = function(e) NULL)
  if (is.null(mf) || !nrow(mf)) {
    return(list(type = "glm", fit = NULL))
  }
  keep <- stats::complete.cases(mf)
  if (!any(keep)) return(list(type = "glm", fit = NULL))
  mf <- mf[keep, , drop = FALSE]
  y <- mf[[1L]]
  x <- mf[-1L]
  ww <- if (is.null(weights)) rep(1, length(y)) else as.numeric(weights)[keep]
  ww[!is.finite(ww) | ww < 0] <- 0
  if (all(ww == 0)) ww <- rep(1, length(ww))

  sl_fit <- tryCatch(
    SuperLearner::SuperLearner(
      Y = y,
      X = x,
      SL.library = sl_library,
      family = gaussian(),
      obsWeights = ww,
      cvControl = list(V = 2L)
    ),
    error = function(e) NULL
  )
  if (is.null(sl_fit)) {
    .warn_once("ltmle_exact_sl_fit_failed",
               "SuperLearner fitting failed; falling back to glm for the LTMLE nuisance backend.",
               call. = FALSE)
    return(list(type = "glm", fit = fit_gaussian_lm_safe(formula, data, weights = weights)))
  }
  list(type = "sl", fit = sl_fit, formula = formula)
}

predict_gaussian_sl_or_glm <- function(obj, newdata, fallback = NA_real_) {
  if (is.null(obj) || is.null(obj$fit)) return(rep(fallback, nrow(newdata)))
  if (identical(obj$type, "glm")) {
    return(safe_predict_lm(obj$fit, newdata, fallback = fallback))
  }
  out <- tryCatch({
    as.numeric(SuperLearner::predict.SuperLearner(obj$fit, newdata = newdata)$pred)
  }, error = function(e) rep(fallback, nrow(newdata)))
  if (length(out) != nrow(newdata) || any(!is.finite(out))) {
    return(rep(fallback, nrow(newdata)))
  }
  out
}

model_sigma_safe <- function(fit, fallback = 1e-6, ml = FALSE) {
  if (is.null(fit)) return(fallback)
  if (is.list(fit) && !is.null(fit$sigma) && is.finite(fit$sigma) && fit$sigma > 0) return(fit$sigma)
  if (is.list(fit) && !is.null(fit$fit) && !is.null(fit$type) &&
      fit$type %in% c("glm", "superlearner")) fit <- fit$fit
  sig <- tryCatch({
    if (isTRUE(ml)) {
      r <- stats::residuals(fit)
      sqrt(sum(r * r) / length(r))
    } else {
      summary(fit)$sigma
    }
  }, error = function(e) NA_real_)
  if (!is.finite(sig) || sig <= 0) fallback else sig
}

normal_density <- function(x, mean, sd, log = FALSE, sd_floor = 1e-6) {
  sd <- pmax(as.numeric(sd), sd_floor)
  stats::dnorm(as.numeric(x), mean = as.numeric(mean), sd = sd, log = log)
}

subject_bootstrap_sample <- function(dat, id_col = "ID", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  assert_cols(dat, id_col, "subject_bootstrap_sample(dat)")
  ids <- unique(dat[[id_col]])
  draw <- sample(ids, size = length(ids), replace = TRUE)
  pieces <- vector("list", length(draw))
  for (ii in seq_along(draw)) {
    dd <- dat[dat[[id_col]] == draw[[ii]], , drop = FALSE]
    dd[[id_col]] <- paste0(draw[[ii]], "_b", ii)
    pieces[[ii]] <- dd
  }
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

summarize_weights <- function(w) {
  q <- safe_quantile(w, c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1))
  ess <- ess_from_weights(w)
  data.frame(
    n = sum(is.finite(w)),
    mean = mean(w, na.rm = TRUE),
    sd = stats::sd(w, na.rm = TRUE),
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

effect_order_map <- function() {
  list(
    yamamuro = c("TE_Y", "IDE_Y", "IIE_Y1", "IIE_Y2", "R_Y"),
    tai = c("ITE_J", "IDE_J", "PSE1", "PSE2", "MI", "xi")
  )
}

add_effect_order <- function(df,
                             estimand_col = "estimand",
                             effect_col = "effect",
                             out_col = "effect_order") {
  assert_cols(df, c(estimand_col, effect_col), "add_effect_order(df)")

  om <- effect_order_map()
  est <- as.character(df[[estimand_col]])
  eff <- as.character(df[[effect_col]])
  ord <- rep(Inf, nrow(df))

  for (nm in names(om)) {
    idx <- which(est == nm)
    if (!length(idx)) next
    m <- match(eff[idx], om[[nm]])
    m[is.na(m)] <- length(om[[nm]]) + 1L
    ord[idx] <- m
  }

  ord[!is.finite(ord)] <- 9999L
  df[[out_col]] <- as.integer(ord)
  df
}

order_effects_long <- function(df) {
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(df)

  if (all(c("estimand", "effect") %in% names(df))) {
    df <- add_effect_order(df, estimand_col = "estimand", effect_col = "effect",
                           out_col = "effect_order")
    df$estimand_order <- match(as.character(df$estimand), names(effect_order_map()))
    df$estimand_order[is.na(df$estimand_order)] <- 9999L
  }

  ord_cols <- intersect(
    c("scenario_id", "PM", "rho", "MI_mode", "n", "estimator", "rep",
      "estimand_order", "estimand", "effect_order", "effect"),
    names(df)
  )
  if (length(ord_cols)) {
    df <- df[do.call(order, c(df[ord_cols], list(na.last = TRUE))), , drop = FALSE]
  }

  if ("effect_order" %in% names(df)) df$effect_order <- NULL
  if ("estimand_order" %in% names(df)) df$estimand_order <- NULL
  rownames(df) <- NULL
  df
}
