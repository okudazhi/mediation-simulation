################################################################################
# summarize.R
#
# Performance metrics for simulation results.
#
# Metrics (per scenario x n x estimator x effect):
#   - Bias, MAE, RMSE, SD
#   - PercentBias (guarded near zero)
#   - Coverage (if CI columns are provided)
################################################################################

summarize_performance <- function(est_long, truth_long, eps = 1e-6) {
  if (!exists(".stop", mode = "function")) {
    stop("core_utils.R must be sourced before summarize.R", call. = FALSE)
  }
  req_est <- c("scenario_id","n","rep","estimator","estimand","effect","estimate")
  req_tru <- c("scenario_id","estimand","effect","truth")
  assert_cols(est_long, req_est, "est_long")
  assert_cols(truth_long, req_tru, "truth_long")

  truth_key <- unique(truth_long[, c("scenario_id","estimand","effect","truth"), drop = FALSE])
  dup_truth <- duplicated(truth_key[, c("scenario_id","estimand","effect"), drop = FALSE])
  if (any(dup_truth)) {
    stop(
      "truth_long contains multiple truth values for the same scenario_id / estimand / effect key.",
      call. = FALSE
    )
  }

  merged <- merge(est_long, truth_key, by = c("scenario_id","estimand","effect"), all.x = TRUE)
  if (any(is.na(merged$truth))) {
    warning("Some truth values are missing after merge. Check scenario_id / estimand / effect keys.")
  }
  merged$error <- merged$estimate - merged$truth

  by_keys <- c("scenario_id","n","estimand","effect","estimator")

  group_key <- interaction(merged[, by_keys, drop = FALSE], drop = TRUE, lex.order = TRUE)
  groups <- split(seq_len(nrow(merged)), group_key)
  rows <- lapply(groups, function(ii) {
    x <- merged[ii, , drop = FALSE]
    e <- x$error
    e_fin <- e[is.finite(e)]
    R_total <- nrow(x)
    R_finite <- length(e_fin)
    Bias <- if (R_finite) mean(e_fin) else NA_real_
    MAE <- if (R_finite) mean(abs(e_fin)) else NA_real_
    RMSE <- if (R_finite) sqrt(mean(e_fin^2)) else NA_real_
    SD <- if (R_finite >= 2L) stats::sd(e_fin) else NA_real_
    MCSE_Bias <- if (R_finite >= 2L) SD / sqrt(R_finite) else NA_real_
    data.frame(
      x[1L, by_keys, drop = FALSE],
      K = R_finite,
      R_total = R_total,
      R_finite = R_finite,
      finite_prop = if (R_total > 0L) R_finite / R_total else NA_real_,
      Bias = Bias,
      MAE = MAE,
      RMSE = RMSE,
      SD = SD,
      MCSE_Bias = MCSE_Bias,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  # Percent bias: 100 * Bias / truth, with stability guard
  out <- merge(out, truth_key,
               by = c("scenario_id","estimand","effect"), all.x = TRUE)
  out$PercentBias <- ifelse(is.finite(out$truth) & abs(out$truth) >= eps,
                            100 * out$Bias / out$truth,
                            NA_real_)

  z <- stats::qnorm(0.975)
  
  out$Bias_lo <- out$Bias - z * out$MCSE_Bias
  out$Bias_hi <- out$Bias + z * out$MCSE_Bias
  
  out$MCSE_PercentBias <- ifelse(
    is.finite(out$truth) & abs(out$truth) >= eps & is.finite(out$MCSE_Bias),
    100 * out$MCSE_Bias / abs(out$truth),
    NA_real_
  )
  out$PercentBias_lo <- out$PercentBias - z * out$MCSE_PercentBias
  out$PercentBias_hi <- out$PercentBias + z * out$MCSE_PercentBias
  

  # Coverage (optional)
  if (all(c("lcl","ucl") %in% names(est_long))) {
    covg <- stats::aggregate(
      list(covered = merged$truth >= merged$lcl & merged$truth <= merged$ucl),
      by = merged[by_keys],
      FUN = function(z) {
        z <- z[!is.na(z)]
        if (!length(z)) NA_real_ else mean(z)
      }
    )
    out <- merge(out, covg, by = by_keys, all.x = TRUE)
    names(out)[names(out) == "covered"] <- "Coverage"
  }

  # Keep nice column order
  base_cols <- c(by_keys, "truth", "K", "R_total", "R_finite", "finite_prop",
                 "Bias","MAE","RMSE","SD","MCSE_Bias","PercentBias",
                 "Bias_lo", "Bias_hi", "PercentBias_lo", "PercentBias_hi")
  if ("Coverage" %in% names(out)) base_cols <- c(base_cols, "Coverage")
  out <- out[, base_cols, drop = FALSE]

  out
}

summarize_worldmean_bias <- function(worldmeans_long, truth_worldmeans) {
  if (!exists(".stop", mode = "function")) {
    stop("core_utils.R must be sourced before summarize.R", call. = FALSE)
  }
  req_wm <- c("scenario_id","rep","estimator","world","mean_hat")
  assert_cols(worldmeans_long, req_wm, "worldmeans_long")
  assert_cols(truth_worldmeans, "scenario_id", "truth_worldmeans")

  world_cols <- grep("^mu_", names(truth_worldmeans), value = TRUE)
  if (length(world_cols) == 0L) {
    stop("truth_worldmeans must contain mu_* world mean columns.", call. = FALSE)
  }

  truth_long <- stats::reshape(
    truth_worldmeans[, c("scenario_id", world_cols), drop = FALSE],
    varying = list(world_cols),
    v.names = "truth",
    timevar = "world",
    times = world_cols,
    direction = "long"
  )
  truth_long$id <- NULL
  rownames(truth_long) <- NULL
  truth_long <- unique(truth_long[, c("scenario_id","world","truth"), drop = FALSE])

  dup_truth <- duplicated(truth_long[, c("scenario_id","world"), drop = FALSE])
  if (any(dup_truth)) {
    stop("truth_worldmeans contains multiple truth values for the same scenario_id / world key.", call. = FALSE)
  }

  detail <- merge(worldmeans_long, truth_long,
                  by = c("scenario_id","world"), all.x = TRUE, sort = FALSE)
  if (any(is.na(detail$truth))) {
    warning("Some world mean truth values are missing after merge. Check scenario_id / world keys.")
  }
  detail$bias <- detail$mean_hat - detail$truth
  detail$abs_bias <- abs(detail$bias)

  group_cols <- intersect(c("estimator","world","Q_model"), names(detail))
  if (!length(group_cols)) group_cols <- c("estimator","world")

  split_key <- interaction(detail[, group_cols, drop = FALSE], drop = TRUE, lex.order = TRUE)
  idx <- split(seq_len(nrow(detail)), split_key)
  summary_rows <- lapply(idx, function(ii) {
    x <- detail[ii, , drop = FALSE]
    b <- x$bias
    b <- b[is.finite(b)]
    data.frame(
      x[1L, group_cols, drop = FALSE],
      R = length(b),
      mean_bias = if (length(b)) mean(b) else NA_real_,
      mean_abs_bias = if (length(b)) mean(abs(b)) else NA_real_,
      rmse = if (length(b)) sqrt(mean(b^2)) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  summary <- if (length(summary_rows)) do.call(rbind, summary_rows) else data.frame()
  rownames(summary) <- NULL

  preferred <- c("scenario_id","analysis_tier","pathway_setting","rho_setting",
                 "structure_setting","PM","rho","rho0","rho1","MI_mode","n",
                 "fold_count","rep","estimator","learner","Q_model","world",
                 "truth","mean_hat","bias","abs_bias")
  detail <- detail[, c(intersect(preferred, names(detail)),
                       setdiff(names(detail), preferred)), drop = FALSE]

  list(detail = detail, summary = summary)
}
