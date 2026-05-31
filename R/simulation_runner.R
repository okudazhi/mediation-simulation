################################################################################
# simulation_runner.R
#
# High-level orchestration for the simulation study.
#
# This file provides two user-facing functions:
#   (1) create_scenario_grid(...)
#       Builds a scenario data.frame over (n, PM, rho, MI_mode).
#
#   (2) run_simulation(...)
#       For each scenario and replication:
#         - simulate one observed dataset under the DGP
#         - run each estimator (provided as a registry of functions)
#         - store estimated effect components in long format
#       It also computes Monte Carlo ground-truth for each unique scenario
#       (common baseline cohort shared across all counterfactual worlds).
#
# NOTE:
#   - Scenario switches are applied via apply_scenario_to_params() in dgp.R
#   - Truth is computed via compute_truth_world_means() in truth.R
#   - Effect mapping is done via compute_estimands_from_means() in core_targets.R
################################################################################

.runner_source_if_needed <- function() {
  req_funs <- c(
    "simulate_dgp_wide",
    "apply_scenario_to_params",
    "compute_truth_world_means",
    "compute_estimands_from_means",
    "primary_effect_keys",
    "filter_primary_nonratio_effects",
    "flatten_effects",
    "summarize_performance",
    "%||%"
  )
  missing <- req_funs[!vapply(req_funs, exists, logical(1), mode = "function")]
  if (length(missing) > 0) {
    stop(
      paste0(
        "Missing required functions: ", paste(missing, collapse = ", "),
        ". Please source R/core_utils.R, R/core_targets.R, R/dgp.R, R/truth.R, estimator files, and R/summarize.R before R/simulation_runner.R."
      ),
      call. = FALSE
    )
  }
}

# ---- Scenario grid -----------------------------------------------------------

#' Create a manuscript-concordant scenario grid.
#'
#' The protocol grid is driven by the scenario manifest factors:
#'   - analysis_tier      : primary / supp. near-null / supp. sensitivity
#'   - pathway_setting    : high / low / near-null
#'   - rho_setting        : resid-corr-none / resid-corr-present
#'   - structure_setting  : none / gamma0-only / full
#'   - Q_model            : correct / wrong
#'   - n                  : DGP sample size
#'   - fold_count         : 2 (primary) or 5 (supp. sensitivity)
#'
#' Additional columns are included for implementation convenience:
#'   - PM, MI_mode, rho0, rho1, rho
#'   - R, B_truth_init, B_MC_init, truncation_rule
#'
#' @param n_vec integer vector of sample sizes
#' @param R replication count stored in the manifest
#' @param B_truth_init truth MC size stored in the manifest
#' @param B_MC_init within-dataset MC size stored in the manifest
#' @param truncation_rule label stored in the manifest
#' @param include_near_null logical; include supplementary near-null cells
#' @param include_sensitivity logical; include 5-fold sensitivity sentinel cells
#' @param rho_mode residual-correlation mode: "rho1_only" or "rho0_rho1"
#' @param rho_low residual mediator error correlation for A_t=0 in resid-corr-present cells
#' @param rho_high residual mediator error correlation for A_t=1 in resid-corr-present cells
#' @return data.frame representing the protocol scenario manifest
create_scenario_grid <- function(
  n_vec,
  R = 500L,
  B_truth_init = 1e6,
  B_MC_init = 2000L,
  truncation_rule = "T1",
  include_near_null = TRUE,
  include_sensitivity = TRUE,
  rho_mode = c("rho1_only", "rho0_rho1"),
  rho_low = 0.0,
  rho_high = 0.5
) {
  rho_mode <- match.arg(rho_mode)
  if (!is.numeric(rho_low) || length(rho_low) != 1L ||
      !is.finite(rho_low) || abs(rho_low) > 1) {
    stop("rho_low must be a finite scalar correlation in [-1, 1].", call. = FALSE)
  }
  if (!is.numeric(rho_high) || length(rho_high) != 1L ||
      !is.finite(rho_high) || abs(rho_high) > 1) {
    stop("rho_high must be a finite scalar correlation in [-1, 1].", call. = FALSE)
  }
  rho_low <- as.numeric(rho_low)
  rho_high <- as.numeric(rho_high)

  primary <- expand.grid(
    analysis_tier = "primary",
    pathway_setting = c("high", "low"),
    rho_setting = c("resid-corr-none", "resid-corr-present"),
    structure_setting = c("none", "gamma0-only", "full"),
    Q_model = c("correct", "wrong"),
    n = as.integer(n_vec),
    fold_count = 2L,
    stringsAsFactors = FALSE
  )

  near_null <- if (isTRUE(include_near_null)) {
    expand.grid(
      analysis_tier = "supp. near-null",
      pathway_setting = "near-null",
      rho_setting = "resid-corr-none",
      structure_setting = "none",
      Q_model = c("correct", "wrong"),
      n = as.integer(n_vec),
      fold_count = 2L,
      stringsAsFactors = FALSE
    )
  } else NULL

  sensitivity <- if (isTRUE(include_sensitivity)) {
    expand.grid(
      analysis_tier = "supp. sensitivity",
      pathway_setting = "high",
      rho_setting = "resid-corr-present",
      structure_setting = "full",
      Q_model = c("correct", "wrong"),
      n = as.integer(n_vec),
      fold_count = 5L,
      stringsAsFactors = FALSE
    )
  } else NULL

  grid <- do.call(rbind, Filter(Negate(is.null), list(primary, near_null, sensitivity)))

  grid$PM <- grid$pathway_setting
  grid$MI_mode <- ifelse(grid$structure_setting == "gamma0-only", "gamma0_only", grid$structure_setting)
  grid$rho0 <- if (rho_mode == "rho0_rho1") {
    ifelse(grid$rho_setting == "resid-corr-present", rho_low, 0.0)
  } else {
    0.0
  }
  grid$rho1 <- ifelse(grid$rho_setting == "resid-corr-present", rho_high, 0.0)
  grid$rho <- grid$rho1
  grid$R <- as.integer(R)
  grid$B_truth_init <- B_truth_init
  grid$B_MC_init <- as.integer(B_MC_init)
  grid$truncation_rule <- truncation_rule

  grid$scenario_id <- vapply(
    seq_len(nrow(grid)),
    function(i) .scenario_id_protocol(
      analysis_tier = grid$analysis_tier[i],
      pathway_setting = grid$pathway_setting[i],
      rho_setting = grid$rho_setting[i],
      structure_setting = grid$structure_setting[i],
      Q_model = grid$Q_model[i],
      n = grid$n[i],
      fold_count = grid$fold_count[i]
    ),
    character(1)
  )

  grid
}

.scenario_id_protocol <- function(analysis_tier, pathway_setting, rho_setting,
                                  structure_setting, Q_model, n, fold_count) {
  tier_code <- c("primary" = "PRI", "supp. near-null" = "SUP-NN", "supp. sensitivity" = "SUP-SEN")[[analysis_tier]]
  p_code <- c("high" = "H", "low" = "L", "near-null" = "NN")[[pathway_setting]]
  r_code <- c("resid-corr-none" = "R0", "resid-corr-present" = "R1")[[rho_setting]]
  s_code <- c("none" = "S0", "gamma0-only" = "Sg", "full" = "Sf")[[structure_setting]]
  q_code <- c("correct" = "Qc", "wrong" = "Qw")[[Q_model]]
  sprintf("%s-%s-%s-%s-%s-N%d-F%d", tier_code, p_code, r_code, s_code, q_code,
          as.integer(n), as.integer(fold_count))
}

.truth_id <- function(PM, rho0, rho1, MI_mode) {
  sprintf(
    "PM=%s;rho0=%s;rho1=%s;MI=%s",
    PM,
    format(rho0, digits = 3),
    format(rho1, digits = 3),
    MI_mode
  )
}

# ---- Internal helpers: truth computation cache ------------------------------

.compute_truth_for_unique_scenarios <- function(scenarios_unique, params0, run_cfg) {
  T <- run_cfg$T
  reg_a  <- run_cfg$reg_a
  reg_as <- run_cfg$reg_as
  B_truth <- run_cfg$B_truth
  truth_n_batches <- run_cfg$truth_n_batches %||% 20L

  if (is.null(B_truth) || !is.numeric(B_truth) || B_truth <= 0) {
    stop("run_cfg$B_truth must be a positive number.", call. = FALSE)
  }

  sc_map <- unique(scenarios_unique)
  if (!("PM" %in% names(sc_map)) && "pathway_setting" %in% names(sc_map)) sc_map$PM <- sc_map$pathway_setting
  if (!("MI_mode" %in% names(sc_map)) && "structure_setting" %in% names(sc_map)) {
    sc_map$MI_mode <- ifelse(sc_map$structure_setting == "gamma0-only", "gamma0_only", sc_map$structure_setting)
  }
  if (!("rho0" %in% names(sc_map)) && "rho" %in% names(sc_map)) sc_map$rho0 <- sc_map$rho
  if (!("rho1" %in% names(sc_map)) && "rho" %in% names(sc_map)) sc_map$rho1 <- sc_map$rho
  if (!("rho" %in% names(sc_map)) && all(c("rho0", "rho1") %in% names(sc_map))) sc_map$rho <- sc_map$rho1

  sc_map$truth_id <- vapply(
    seq_len(nrow(sc_map)),
    function(i) .truth_id(sc_map$PM[i], sc_map$rho0[i], sc_map$rho1[i], sc_map$MI_mode[i]),
    character(1)
  )

  dgp_unique <- unique(sc_map[, c("truth_id", "PM", "rho0", "rho1", "MI_mode"), drop = FALSE])
  seed_base <- if (!is.null(run_cfg$seed)) as.integer(run_cfg$seed) else 1L

  truth_worldmeans_rows <- list()
  truth_effects_rows <- list()

  for (k in seq_len(nrow(dgp_unique))) {
    sc <- dgp_unique[k, , drop = FALSE]

    pm_mult <- NULL
    if (!is.null(run_cfg$PM_MAP)) {
      pm_mult <- run_cfg$PM_MAP[[as.character(sc$PM)]]
    }

    params_sc <- apply_scenario_to_params(
      params0,
      PM = pm_mult,
      MI_mode = sc$MI_mode,
      rho0 = sc$rho0,
      rho1 = sc$rho1
    )

    truth_seed <- seed_base + 10000L * k
    truth <- compute_truth_world_means(
      params = params_sc,
      T = T,
      reg_a = reg_a,
      reg_as = reg_as,
      B_truth = B_truth,
      seed = truth_seed,
      n_batches = truth_n_batches
    )

    means_vec <- unlist(truth$means)
    means_mcse_vec <- unlist(truth$means_mcse %||% list())
    if (length(means_mcse_vec)) {
      names(means_mcse_vec) <- paste0(names(means_mcse_vec), "_mcse")
    }
    truth_worldmeans_rows[[k]] <- data.frame(
      truth_id = sc$truth_id,
      PM = sc$PM,
      rho0 = sc$rho0,
      rho1 = sc$rho1,
      rho = sc$rho1,
      MI_mode = sc$MI_mode,
      t(means_vec),
      t(means_mcse_vec),
      row.names = NULL,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    eff_df_all <- flatten_effects(truth$effects)
    eff_df <- filter_primary_nonratio_effects(eff_df_all)
    eff_mcse_df <- if (!is.null(truth$effects_mcse)) {
      filter_primary_nonratio_effects(flatten_effects(truth$effects_mcse))
    } else {
      data.frame(estimand = eff_df$estimand, effect = eff_df$effect, estimate = NA_real_)
    }
    names(eff_mcse_df)[names(eff_mcse_df) == "estimate"] <- "truth_mcse"
    eff_df <- merge(eff_df, eff_mcse_df, by = c("estimand", "effect"), all.x = TRUE, sort = FALSE)
    truth_effects_rows[[k]] <- data.frame(
      truth_id = sc$truth_id,
      PM = sc$PM,
      rho0 = sc$rho0,
      rho1 = sc$rho1,
      rho = sc$rho1,
      MI_mode = sc$MI_mode,
      estimand = eff_df$estimand,
      effect = eff_df$effect,
      truth = eff_df$estimate,
      truth_mcse = eff_df$truth_mcse,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }

  truth_worldmeans_base <- do.call(rbind, truth_worldmeans_rows)
  truth_effects_base <- do.call(rbind, truth_effects_rows)

  truth_worldmeans <- merge(
    sc_map,
    truth_worldmeans_base,
    by = c("truth_id", "PM", "rho0", "rho1", "rho", "MI_mode"),
    all.x = TRUE,
    sort = FALSE
  )
  if ("scenario_id" %in% names(sc_map)) {
    truth_worldmeans <- truth_worldmeans[match(sc_map$scenario_id, truth_worldmeans$scenario_id), , drop = FALSE]
  }

  truth_effects <- merge(
    sc_map,
    truth_effects_base,
    by = c("truth_id", "PM", "rho0", "rho1", "rho", "MI_mode"),
    all.x = TRUE,
    sort = FALSE
  )
  if ("scenario_id" %in% names(sc_map)) {
    truth_effects$order_key <- match(truth_effects$scenario_id, sc_map$scenario_id)
    truth_effects <- truth_effects[order(truth_effects$order_key, truth_effects$estimand, truth_effects$effect), , drop = FALSE]
    truth_effects$order_key <- NULL
  }

  means_cache <- setNames(vector("list", nrow(truth_worldmeans)), truth_worldmeans$scenario_id)
  for (i in seq_len(nrow(truth_worldmeans))) {
    row_i <- truth_worldmeans[i, , drop = FALSE]
    keys <- c(
      "mu_nat_a", "mu_nat_as",
      "mu_joint_aa", "mu_joint_asas", "mu_joint_aas",
      "mu_sep_aaa", "mu_sep_asas_asas", "mu_sep_a_asas", "mu_sep_a_aas"
    )
    means_cache[[row_i$scenario_id]] <- as.list(as.numeric(row_i[1, keys]))
    names(means_cache[[row_i$scenario_id]]) <- keys
  }

  list(
    truth_worldmeans = truth_worldmeans,
    truth_effects = truth_effects,
    means_cache = means_cache
  )
}

# ---- Internal helpers: bootstrap CIs ----------------------------------------

.bootstrap_effect_cis <- function(dat_wide, cfg, estimator_fun,
                                  B_boot = 200,
                                  conf_level = 0.95,
                                  seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(dat_wide)
  alpha <- 1 - conf_level

  # One pass to establish effect keys.
  res0 <- estimator_fun(dat_wide, cfg)
  if (is.null(res0) || is.null(res0$means)) return(NULL)
  eff0 <- filter_primary_nonratio_effects(flatten_effects(compute_estimands_from_means(res0$means)))
  k_names <- paste(eff0$estimand, eff0$effect, sep = "|")
  K <- length(k_names)

  boot_mat <- matrix(NA_real_, nrow = B_boot, ncol = K)
  colnames(boot_mat) <- k_names

  for (b in 1:B_boot) {
    idx <- sample.int(n, size = n, replace = TRUE)
    dat_b <- dat_wide[idx, , drop = FALSE]

    # Seed inside bootstrap (optional): keep reproducible but vary across b.
    if (!is.null(seed)) set.seed(seed + b)

    rb <- tryCatch(estimator_fun(dat_b, cfg), error = function(e) NULL)
    if (is.null(rb) || is.null(rb$means)) next

    effb <- filter_primary_nonratio_effects(flatten_effects(compute_estimands_from_means(rb$means)))
    kb <- paste(effb$estimand, effb$effect, sep = "|")

    # Align by key; leave NA if missing.
    m <- match(k_names, kb)
    boot_mat[b, ] <- effb$estimate[m]
  }

  # Percentile CI (ignoring NA rows per component).
  lcl <- apply(boot_mat, 2, function(x) stats::quantile(x, probs = alpha / 2, na.rm = TRUE, names = FALSE))
  ucl <- apply(boot_mat, 2, function(x) stats::quantile(x, probs = 1 - alpha / 2, na.rm = TRUE, names = FALSE))

  out <- data.frame(
    key = k_names,
    lcl = as.numeric(lcl),
    ucl = as.numeric(ucl),
    stringsAsFactors = FALSE
  )

  # Split key back into (estimand, effect).
  parts <- strsplit(out$key, "\\|", fixed = FALSE)
  out$estimand <- vapply(parts, function(z) z[1], character(1))
  out$effect <- vapply(parts, function(z) z[2], character(1))
  out$key <- NULL

  out
}

# ---- Diagnostics helpers -----------------------------------------------------

.diag_default_cfg <- function(run_cfg) {
  d <- run_cfg$diagnostics %||% list()
  d$enabled   <- isTRUE(d$enabled)
  d$dir       <- d$dir %||% file.path(run_cfg$output_dir %||% "output", "diagnostics")
  d$save_full <- d$save_full %||% "none"   # none/outlier/all
  d$save_score_equations <- d$save_score_equations %||% "outlier"
  if (!d$save_score_equations %in% c("none", "outlier", "all")) d$save_score_equations <- "outlier"
  d$thr <- d$thr %||% list()
  d$thr$max_abs_effect <- d$thr$max_abs_effect %||% 1e6
  d$thr$max_H          <- d$thr$max_H %||% 500
  d$thr$min_ess_frac   <- d$thr$min_ess_frac %||% 0.05
  d$thr$max_abs_eps    <- d$thr$max_abs_eps %||% 25
  d$outlier_triggers <- d$outlier_triggers %||% c(
    "nonfinite estimate",
    "estimator failure",
    "absolute effect estimate exceeds diagnostics$thr$max_abs_effect",
    "ltmle_exact clever covariate H exceeds diagnostics$thr$max_H",
    "ltmle_exact ESS fraction below diagnostics$thr$min_ess_frac",
    "absolute targeting epsilon exceeds diagnostics$thr$max_abs_eps"
  )
  d
}

.diag_safe_id <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

.scenario_label <- function(sc) {
  if ("scenario_label" %in% names(sc) && length(sc$scenario_label) &&
      !is.na(sc$scenario_label[1L]) && nzchar(as.character(sc$scenario_label[1L]))) {
    return(as.character(sc$scenario_label[1L]))
  }
  if ("scenario_id" %in% names(sc) && length(sc$scenario_id)) {
    return(as.character(sc$scenario_id[1L]))
  }
  NA_character_
}

.diag_collapse_messages <- function(x, max_messages = 3L, max_chars = 500L) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(NA_character_)
  out <- paste(utils::head(x, max_messages), collapse = " | ")
  if (length(x) > max_messages) {
    out <- paste0(out, " | ... +", length(x) - max_messages, " more")
  }
  if (nchar(out, type = "chars") > max_chars) {
    out <- paste0(substr(out, 1L, max_chars - 3L), "...")
  }
  out
}

.estimator_variant_from_result <- function(est_name, res = NULL) {
  if (!is.null(res)) {
    if (!is.null(res$metadata$estimator_variant) &&
        length(res$metadata$estimator_variant) &&
        !is.na(res$metadata$estimator_variant[1L])) {
      return(as.character(res$metadata$estimator_variant[1L]))
    }
    run_summary <- res$diagnostics$ltmle_exact_run_summary %||% NULL
    if (!is.null(run_summary) && is.data.frame(run_summary) &&
        nrow(run_summary) && "estimator_variant" %in% names(run_summary)) {
      return(as.character(run_summary$estimator_variant[1L]))
    }
  }
  as.character(est_name)
}

.truncation_metadata_from_result <- function(res = NULL) {
  out <- list(
    truncation_policy = NA_character_,
    truncation_target = NA_character_,
    selected_primary_law_integration_n = NA_integer_
  )
  if (is.null(res) || is.null(res$diagnostics)) return(out)
  tr <- res$diagnostics$truncation_diagnostics %||% NULL
  if (!is.null(tr) && is.data.frame(tr) && nrow(tr)) {
    if ("truncation_policy" %in% names(tr)) {
      out$truncation_policy <- as.character(tr$truncation_policy[1L])
    }
    if ("truncation_target" %in% names(tr)) {
      out$truncation_target <- as.character(tr$truncation_target[1L])
    }
  }
  run_summary <- res$diagnostics$ltmle_exact_run_summary %||% NULL
  if (!is.null(run_summary) && is.data.frame(run_summary) && nrow(run_summary)) {
    if ("truncation_policy" %in% names(run_summary) && is.na(out$truncation_policy)) {
      out$truncation_policy <- as.character(run_summary$truncation_policy[1L])
    }
    if ("truncation_target" %in% names(run_summary) && is.na(out$truncation_target)) {
      out$truncation_target <- as.character(run_summary$truncation_target[1L])
    }
    if ("ltmle_exact_law_integration_n" %in% names(run_summary)) {
      out$selected_primary_law_integration_n <- as.integer(run_summary$ltmle_exact_law_integration_n[1L])
    }
  }
  out
}

.is_ltmle_output_name <- function(est_name) {
  est_name <- as.character(est_name)
  identical(est_name, "ltmle_glm") ||
    identical(est_name, "ltmle_sl") ||
    identical(est_name, "ltmle_exact") ||
    grepl("^ltmle_exact_", est_name)
}

.diag_ess <- function(w) {
  w <- w[is.finite(w)]
  n <- length(w)
  if (n == 0) return(c(ess=NA_real_, ess_frac=NA_real_))
  sw <- sum(w); sw2 <- sum(w*w)
  if (!is.finite(sw) || !is.finite(sw2) || sw2 <= 0) return(c(ess=NA_real_, ess_frac=NA_real_))
  ess <- (sw*sw)/sw2
  c(ess=ess, ess_frac=ess/n)
}


.diag_ltmle_exact_score_equations <- function(res, sc, rep, est_name) {
  if (is.null(res$debug) || is.null(res$debug$cc_sep) || is.null(res$debug$cc_joint)) return(NULL)

  cc_all <- rbind(res$debug$cc_sep, res$debug$cc_joint)

  # eps_log is one row per (phi, t) produced by the targeting recursion.
  eps_all <- NULL
  if (!is.null(res$debug$eps_log_sep) && !is.null(res$debug$eps_log_joint)) {
    eps_all <- rbind(res$debug$eps_log_sep, res$debug$eps_log_joint)

    keep_eps <- intersect(c("phi", "t", "eps_hat",
                            "mc_H_max", "mc_H_ess_frac", "mc_prop_trunc_hi",
                            "raw_log_ratio_max"),
                          names(eps_all))
    eps_all <- eps_all[, keep_eps, drop = FALSE]
  }

  groups <- split(cc_all, interaction(cc_all$phi, cc_all$t, drop = TRUE))

  rows <- lapply(groups, function(df) {
    phi <- unique(df$phi)[1]
    t   <- unique(df$t)[1]

    # Cumulative-H truncation bounds (scalar per phi×t).
    lo <- unique(round(df$H_trunc_lo, 12))
    hi <- unique(round(df$H_trunc_hi, 12))
    lo <- if (length(lo)) as.numeric(lo[1]) else NA_real_
    hi <- if (length(hi)) as.numeric(hi[1]) else NA_real_

    rlo <- unique(round(df$H_trunc_rate_lo, 12))
    rhi <- unique(round(df$H_trunc_rate_hi, 12))
    rlo <- if (length(rlo)) as.numeric(rlo[1]) else NA_real_
    rhi <- if (length(rhi)) as.numeric(rhi[1]) else NA_real_

    H_fin <- df$H_targ[is.finite(df$H_targ)]
    H_max <- if (length(H_fin)) max(H_fin) else NA_real_

    ess <- .diag_ess(df$H_targ)

    data.frame(
      phi = phi,
      t = t,
      is_TE = unique(df$is_TE)[1],
      target_A = unique(df$target_A)[1],
      pat_m1 = unique(df$pat_m1_step)[1],
      pat_m2 = unique(df$pat_m2_step)[1],
      n_id = length(unique(df$id)),
      I_match_mean = mean(df$I_match, na.rm = TRUE),

      H_trunc_lo = lo,
      H_trunc_hi = hi,
      H_trunc_rate_lo = rlo,
      H_trunc_rate_hi = rhi,

      H_max = H_max,
      H_ess = ess["ess"],
      H_ess_frac = ess["ess_frac"],

      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)

  if (!is.null(eps_all) && nrow(eps_all) > 0) {
    out <- merge(out, eps_all, by = c("phi","t"), all.x = TRUE, sort = FALSE)
  }

  out$scenario_id <- sc$scenario_id
  out$n <- sc$n
  out$rep <- rep
  out$estimator <- est_name

  keep_first <- c("scenario_id", "n", "rep", "estimator", "phi", "t")
  out <- out[, c(keep_first, setdiff(names(out), keep_first)), drop = FALSE]
  out[order(out$phi, out$t), , drop = FALSE]
}


.diag_ltmle_exact_run_summary <- function(ltmle_exact_score_equations, sc, rep, est_name, res = NULL) {
  if (is.null(ltmle_exact_score_equations) ||
      !is.data.frame(ltmle_exact_score_equations) ||
      nrow(ltmle_exact_score_equations) == 0) return(NULL)

  tm0 <- ltmle_exact_score_equations

  # Keep only the "nontrivial" part where truncation bounds are not degenerate.
  # This avoids flagging deterministic cases as outliers.
  if (all(c("t", "is_TE", "H_trunc_lo", "H_trunc_hi") %in% names(tm0))) {
    ok <- (tm0$t >= 2) & (tm0$is_TE == 0) &
      is.finite(tm0$H_trunc_lo) & is.finite(tm0$H_trunc_hi) &
      abs(tm0$H_trunc_hi - tm0$H_trunc_lo) > 1e-8
    tm0 <- tm0[ok, , drop = FALSE]
  } else {
    if ("t" %in% names(tm0)) tm0 <- tm0[tm0$t >= 2, , drop = FALSE]
    if ("is_TE" %in% names(tm0)) tm0 <- tm0[tm0$is_TE == 0, , drop = FALSE]
  }

  .max_fin <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    max(x)
  }
  .min_fin <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    min(x)
  }

  out <- data.frame(
    scenario_id = sc$scenario_id,
    n = sc$n,
    rep = rep,
    estimator = est_name,

    max_H_max = if ("H_max" %in% names(tm0)) .max_fin(tm0$H_max) else NA_real_,
    min_H_ess_frac = if ("H_ess_frac" %in% names(tm0)) .min_fin(tm0$H_ess_frac) else NA_real_,
    max_abs_eps = if ("eps_hat" %in% names(tm0)) .max_fin(abs(tm0$eps_hat)) else NA_real_,

    max_mc_H_max = if ("mc_H_max" %in% names(tm0)) .max_fin(tm0$mc_H_max) else NA_real_,
    min_mc_H_ess_frac = if ("mc_H_ess_frac" %in% names(tm0)) .min_fin(tm0$mc_H_ess_frac) else NA_real_,
    max_mc_prop_trunc_hi = if ("mc_prop_trunc_hi" %in% names(tm0)) .max_fin(tm0$mc_prop_trunc_hi) else NA_real_,

    max_raw_log_ratio_max = if ("raw_log_ratio_max" %in% names(tm0)) .max_fin(tm0$raw_log_ratio_max) else NA_real_,

    stringsAsFactors = FALSE
  )

  out
}


.diag_is_outlier <- function(diag_cfg, eff_df, ltmle_exact_score_equations = NULL) {
  if (!isTRUE(diag_cfg$enabled)) return(FALSE)

  thr <- diag_cfg$thr

  # 1) effect estimates contain NA/Inf
  if (any(!is.finite(eff_df$estimate))) return(TRUE)

  # 2) effect is extreme (coarse screen)
  if (is.finite(thr$max_abs_effect) &&
      max(abs(eff_df$estimate), na.rm = TRUE) > thr$max_abs_effect) return(TRUE)

  # 3) LTMLE-exact-specific screens: evaluate only on the "nontrivial" part
  #    (t>=2, non-TE, and non-degenerate H-truncation bounds).
  if (!is.null(ltmle_exact_score_equations) &&
      is.data.frame(ltmle_exact_score_equations) &&
      nrow(ltmle_exact_score_equations) > 0) {
    tm0 <- ltmle_exact_score_equations

    if (all(c("t", "is_TE", "H_trunc_lo", "H_trunc_hi") %in% names(tm0))) {
      ok <- (tm0$t >= 2) & (tm0$is_TE == 0) &
        is.finite(tm0$H_trunc_lo) & is.finite(tm0$H_trunc_hi) &
        abs(tm0$H_trunc_hi - tm0$H_trunc_lo) > 1e-8
      tm0 <- tm0[ok, , drop = FALSE]
    } else {
      if ("t" %in% names(tm0)) tm0 <- tm0[tm0$t >= 2, , drop = FALSE]
      if ("is_TE" %in% names(tm0)) tm0 <- tm0[tm0$is_TE == 0, , drop = FALSE]
    }

    if (nrow(tm0) > 0) {
      if (is.finite(thr$max_H) && ("H_max" %in% names(tm0)) &&
          max(tm0$H_max, na.rm = TRUE) > thr$max_H) return(TRUE)
      if (is.finite(thr$max_H) && ("max_abs_H" %in% names(tm0)) &&
          max(tm0$max_abs_H, na.rm = TRUE) > thr$max_H) return(TRUE)

      if (is.finite(thr$min_ess_frac) && ("H_ess_frac" %in% names(tm0)) &&
          min(tm0$H_ess_frac, na.rm = TRUE) < thr$min_ess_frac) return(TRUE)

      if ("eps_hat" %in% names(tm0) && is.finite(thr$max_abs_eps) &&
          max(abs(tm0$eps_hat), na.rm = TRUE) > thr$max_abs_eps) return(TRUE)
      if ("eps" %in% names(tm0) && is.finite(thr$max_abs_eps) &&
          max(abs(tm0$eps), na.rm = TRUE) > thr$max_abs_eps) return(TRUE)
    }
  }

  FALSE
}

.diag_should_save_full <- function(diag_cfg, is_outlier) {
  if (!isTRUE(diag_cfg$enabled)) return(FALSE)

  mode <- diag_cfg$save_full
  if (identical(mode, "none")) return(FALSE)
  if (identical(mode, "all")) return(TRUE)
  if (identical(mode, "outlier")) return(isTRUE(is_outlier))

  FALSE
}

.diag_save_rds <- function(obj, diag_dir_full, sc, rep, est_name) {
  safe_sc <- .diag_safe_id(sc$scenario_id)
  fname <- sprintf("%s__n%d__rep%03d__%s.rds", safe_sc, as.integer(sc$n), as.integer(rep), est_name)
  path <- file.path(diag_dir_full, fname)
  ok <- tryCatch({ saveRDS(obj, path, compress = "xz"); TRUE }, error = function(e) FALSE)
  if (ok) path else NA_character_
}

.runner_rbind_fill <- function(rows) {
  rows <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, rows)
  if (!length(rows)) return(NULL)
  cols <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(x) {
    missing <- setdiff(cols, names(x))
    for (cc in missing) x[[cc]] <- NA
    x[, cols, drop = FALSE]
  })
  do.call(rbind, rows)
}

.attempt_status_row <- function(sc, rep, seed, est_name, estimator_variant,
                                attempted = TRUE, completed = FALSE,
                                error_condition = NULL, warnings = character(),
                                eff_df = NULL, means_vec = NULL,
                                elapsed_seconds = NA_real_,
                                diagnostics_written = FALSE,
                                output_rows_written = FALSE,
                                truncation_policy = NA_character_,
                                truncation_target = NA_character_) {
  n_effect_rows <- if (is.data.frame(eff_df)) nrow(eff_df) else 0L
  n_finite_effect_rows <- if (is.data.frame(eff_df) && "estimate" %in% names(eff_df)) {
    sum(is.finite(suppressWarnings(as.numeric(eff_df$estimate))))
  } else 0L
  n_worldmean_rows <- length(means_vec %||% numeric(0))
  n_finite_worldmean_rows <- if (n_worldmean_rows) {
    sum(is.finite(suppressWarnings(as.numeric(means_vec))))
  } else 0L
  failed <- attempted && !completed
  status <- if (!attempted) "not_applicable" else if (failed) "failed" else "completed"
  nonfinite_effect <- n_effect_rows > 0L && n_finite_effect_rows < n_effect_rows
  nonfinite_worldmean <- n_worldmean_rows > 0L && n_finite_worldmean_rows < n_worldmean_rows
  passed_basic <- isTRUE(attempted) && isTRUE(completed) &&
    n_effect_rows > 0L && n_worldmean_rows > 0L &&
    !isTRUE(nonfinite_effect) && !isTRUE(nonfinite_worldmean)
  error_message <- if (!is.null(error_condition)) conditionMessage(error_condition) else NA_character_
  error_class <- if (!is.null(error_condition)) class(error_condition)[1L] else NA_character_
  failure_class <- if (!attempted) {
    "not_applicable"
  } else if (!completed) {
    error_class %||% "estimator_failed"
  } else if (!passed_basic) {
    "basic_output_check_failed"
  } else {
    "no_failure"
  }
  data.frame(
    scenario_id = as.character(sc$scenario_id),
    scenario_label = .scenario_label(sc),
    n = as.integer(sc$n),
    rep = as.integer(rep),
    seed = as.integer(seed),
    estimator = as.character(est_name),
    estimator_variant = as.character(estimator_variant),
    truncation_policy = as.character(truncation_policy),
    truncation_target = as.character(truncation_target),
    attempted = isTRUE(attempted),
    completed = isTRUE(completed),
    failed = isTRUE(failed),
    status = status,
    error_class = error_class,
    error_message = error_message,
    warning_count = length(warnings),
    warning_messages_summary = .diag_collapse_messages(warnings),
    n_effect_rows = as.integer(n_effect_rows),
    n_finite_effect_rows = as.integer(n_finite_effect_rows),
    n_worldmean_rows = as.integer(n_worldmean_rows),
    n_finite_worldmean_rows = as.integer(n_finite_worldmean_rows),
    nonfinite_effect_estimate = isTRUE(nonfinite_effect),
    nonfinite_worldmean_estimate = isTRUE(nonfinite_worldmean),
    elapsed_seconds = as.numeric(elapsed_seconds),
    memory_peak_mb = NA_real_,
    output_rows_written = isTRUE(output_rows_written),
    diagnostics_written = isTRUE(diagnostics_written),
    passed_basic_output_check = isTRUE(passed_basic),
    failure_class = failure_class,
    failure_detail = if (!is.na(error_message)) error_message else {
      if (identical(failure_class, "no_failure")) NA_character_ else failure_class
    },
    stringsAsFactors = FALSE
  )
}

.summarize_estimator_runtime <- function(attempt_status) {
  if (is.null(attempt_status) || !is.data.frame(attempt_status) || !nrow(attempt_status)) {
    return(NULL)
  }
  key_cols <- c("scenario_id", "scenario_label", "n", "estimator", "estimator_variant")
  for (cc in setdiff(key_cols, names(attempt_status))) attempt_status[[cc]] <- NA
  split_key <- interaction(attempt_status[, key_cols, drop = FALSE], drop = TRUE, lex.order = TRUE)
  idx <- split(seq_len(nrow(attempt_status)), split_key)
  .q <- function(x, p) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    as.numeric(stats::quantile(x, probs = p, na.rm = TRUE, names = FALSE))
  }
  rows <- lapply(idx, function(ii) {
    x <- attempt_status[ii, , drop = FALSE]
    elapsed <- suppressWarnings(as.numeric(x$elapsed_seconds))
    mem <- suppressWarnings(as.numeric(x$memory_peak_mb))
    n_attempts <- sum(as.logical(x$attempted), na.rm = TRUE)
    n_completed <- sum(as.logical(x$completed), na.rm = TRUE)
    n_failed <- sum(as.logical(x$failed), na.rm = TRUE)
    data.frame(
      x[1L, key_cols, drop = FALSE],
      n_attempts = as.integer(n_attempts),
      n_completed = as.integer(n_completed),
      n_failed = as.integer(n_failed),
      mean_elapsed_seconds = if (any(is.finite(elapsed))) mean(elapsed, na.rm = TRUE) else NA_real_,
      sd_elapsed_seconds = if (sum(is.finite(elapsed)) >= 2L) stats::sd(elapsed, na.rm = TRUE) else NA_real_,
      median_elapsed_seconds = .q(elapsed, 0.50),
      p75_elapsed_seconds = .q(elapsed, 0.75),
      p90_elapsed_seconds = .q(elapsed, 0.90),
      p95_elapsed_seconds = .q(elapsed, 0.95),
      max_elapsed_seconds = if (any(is.finite(elapsed))) max(elapsed, na.rm = TRUE) else NA_real_,
      min_elapsed_seconds = if (any(is.finite(elapsed))) min(elapsed, na.rm = TRUE) else NA_real_,
      total_elapsed_seconds = if (any(is.finite(elapsed))) sum(elapsed, na.rm = TRUE) else NA_real_,
      failure_rate = if (n_attempts > 0L) n_failed / n_attempts else NA_real_,
      nonfinite_effect_rate = if (n_attempts > 0L) {
        sum(as.logical(x$nonfinite_effect_estimate), na.rm = TRUE) / n_attempts
      } else NA_real_,
      nonfinite_worldmean_rate = if (n_attempts > 0L) {
        sum(as.logical(x$nonfinite_worldmean_estimate), na.rm = TRUE) / n_attempts
      } else NA_real_,
      memory_peak_mb_mean = if (any(is.finite(mem))) mean(mem, na.rm = TRUE) else NA_real_,
      memory_peak_mb_p95 = .q(mem, 0.95),
      memory_peak_mb_max = if (any(is.finite(mem))) max(mem, na.rm = TRUE) else NA_real_,
      runtime_summary_status = "completed",
      failure_class = "no_failure",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.summarize_truncation_diagnostics <- function(truncation_diagnostics,
                                              attempt_status = NULL) {
  tr <- truncation_diagnostics
  if (!is.null(tr) && is.data.frame(tr) && nrow(tr)) {
    if (!"estimator_name" %in% names(tr)) {
      tr$estimator_name <- tr$estimator %||% NA_character_
    }
    if (!"estimator_variant" %in% names(tr)) tr$estimator_variant <- tr$estimator_name
    if (!"scenario_label" %in% names(tr)) {
      label_map <- if (!is.null(attempt_status) && "scenario_label" %in% names(attempt_status)) {
        unique(attempt_status[, c("scenario_id", "scenario_label"), drop = FALSE])
      } else NULL
      tr$scenario_label <- if (!is.null(label_map)) {
        label_map$scenario_label[match(tr$scenario_id, label_map$scenario_id)]
      } else tr$scenario_id
    }
    if (!"n" %in% names(tr) && !is.null(attempt_status) && "n" %in% names(attempt_status)) {
      n_map <- unique(attempt_status[, c("scenario_id", "n"), drop = FALSE])
      tr$n <- n_map$n[match(tr$scenario_id, n_map$scenario_id)]
    }
  } else {
    tr <- data.frame()
  }

  rows <- list()
  key_cols <- c("scenario_id", "scenario_label", "n", "estimator", "estimator_variant")
  if (nrow(tr)) {
    tr$estimator <- as.character(tr$estimator_name)
    for (cc in setdiff(key_cols, names(tr))) tr[[cc]] <- NA
    split_key <- interaction(tr[, key_cols, drop = FALSE], drop = TRUE, lex.order = TRUE)
    idx <- split(seq_len(nrow(tr)), split_key)
    .num <- function(x) suppressWarnings(as.numeric(x))
    .wmean <- function(value, weight = NULL) {
      value <- .num(value)
      if (is.null(weight)) weight <- rep(1, length(value))
      weight <- .num(weight)
      ok <- is.finite(value) & is.finite(weight) & weight >= 0
      if (!any(ok)) return(NA_real_)
      if (sum(weight[ok]) <= 0) return(mean(value[ok], na.rm = TRUE))
      stats::weighted.mean(value[ok], weight[ok], na.rm = TRUE)
    }
    .first_chr <- function(x, default = NA_character_) {
      x <- as.character(x)
      x <- x[!is.na(x) & nzchar(x)]
      if (length(x)) x[1L] else default
    }
    .min_fin <- function(x) {
      x <- .num(x)
      if (any(is.finite(x))) min(x, na.rm = TRUE) else NA_real_
    }
    .max_fin <- function(x) {
      x <- .num(x)
      if (any(is.finite(x))) max(x, na.rm = TRUE) else NA_real_
    }
    rows <- lapply(idx, function(ii) {
      x <- tr[ii, , drop = FALSE]
      n_values <- if ("n_values" %in% names(x)) .num(x$n_values) else rep(1, nrow(x))
      frac_total <- if ("fraction_truncated_total" %in% names(x)) .num(x$fraction_truncated_total) else rep(0, nrow(x))
      fixed_used <- if ("fixed_bound_truncation_used" %in% names(x)) {
        any(as.logical(x$fixed_bound_truncation_used), na.rm = TRUE)
      } else FALSE
      dr_used <- if ("density_ratio_factor_truncation_used" %in% names(x)) {
        any(as.logical(x$density_ratio_factor_truncation_used), na.rm = TRUE)
      } else FALSE
      prob_is_trunc <- if ("probability_bounding_is_truncation" %in% names(x)) {
        any(as.logical(x$probability_bounding_is_truncation), na.rm = TRUE)
      } else FALSE
      data.frame(
        x[1L, key_cols, drop = FALSE],
        truncation_enabled = if ("truncation_enabled" %in% names(x)) any(as.logical(x$truncation_enabled), na.rm = TRUE) else FALSE,
        truncation_policy = .first_chr(x$truncation_policy %||% NA_character_, "not_applicable"),
        truncation_target = .first_chr(x$truncation_target %||% NA_character_, "not_applicable"),
        truncation_rule = .first_chr(x$truncation_rule %||% NA_character_, "not_applicable"),
        requested_quantile_lower = .wmean(x$requested_quantile_lower %||% NA_real_, n_values),
        requested_quantile_upper = .wmean(x$requested_quantile_upper %||% NA_real_, n_values),
        effective_quantile_lower = .wmean(x$effective_quantile_lower %||% NA_real_, n_values),
        effective_quantile_upper = .wmean(x$effective_quantile_upper %||% NA_real_, n_values),
        n_diagnostic_rows = nrow(x),
        n_values_total = sum(n_values[is.finite(n_values)], na.rm = TRUE),
        mean_fraction_truncated_lower = .wmean(x$fraction_truncated_lower %||% NA_real_, n_values),
        mean_fraction_truncated_upper = .wmean(x$fraction_truncated_upper %||% NA_real_, n_values),
        mean_fraction_truncated_total = .wmean(frac_total, n_values),
        p95_fraction_truncated_total = if (any(is.finite(frac_total))) {
          as.numeric(stats::quantile(frac_total, 0.95, na.rm = TRUE, names = FALSE))
        } else NA_real_,
        max_fraction_truncated_total = if (any(is.finite(frac_total))) max(frac_total, na.rm = TRUE) else NA_real_,
        mean_raw_value = .wmean(x$mean_raw_value %||% NA_real_, n_values),
        sd_raw_value = .wmean(x$sd_raw_value %||% NA_real_, n_values),
        min_raw_value = .min_fin(x$min_raw_value %||% NA_real_),
        p01_raw_value = .wmean(x$p01_raw_value %||% NA_real_, n_values),
        median_raw_value = .wmean(x$median_raw_value %||% NA_real_, n_values),
        p99_raw_value = .wmean(x$p99_raw_value %||% NA_real_, n_values),
        max_raw_value = .max_fin(x$max_raw_value %||% NA_real_),
        mean_truncated_value = .wmean(x$mean_truncated_value %||% NA_real_, n_values),
        sd_truncated_value = .wmean(x$sd_truncated_value %||% NA_real_, n_values),
        min_truncated_value = .min_fin(x$min_truncated_value %||% NA_real_),
        max_truncated_value = .max_fin(x$max_truncated_value %||% NA_real_),
        mean_ess_raw = .wmean(x$ess_raw %||% NA_real_, n_values),
        mean_ess_truncated = .wmean(x$ess_truncated %||% NA_real_, n_values),
        min_ess_truncated = .min_fin(x$ess_truncated %||% NA_real_),
        n_rows_with_any_truncation = sum(is.finite(frac_total) & frac_total > 0, na.rm = TRUE),
        fraction_rows_with_any_truncation = if (nrow(x) > 0L) {
          sum(is.finite(frac_total) & frac_total > 0, na.rm = TRUE) / nrow(x)
        } else NA_real_,
        fixed_bound_truncation_used = fixed_used,
        density_ratio_factor_truncation_used = dr_used,
        probability_bounding_is_truncation = prob_is_trunc,
        summary_status = "completed",
        failure_class = if (fixed_used) {
          "fixed_bound_truncation_used"
        } else if (dr_used) {
          "density_ratio_factor_truncation_used"
        } else "no_failure",
        stringsAsFactors = FALSE
      )
    })
  }
  out <- if (length(rows)) do.call(rbind, rows) else data.frame()

  if (!is.null(attempt_status) && is.data.frame(attempt_status) && nrow(attempt_status)) {
    attempts <- unique(attempt_status[, intersect(
      c("scenario_id", "scenario_label", "n", "estimator", "estimator_variant",
        "truncation_policy", "truncation_target"),
      names(attempt_status)
    ), drop = FALSE])
    for (cc in setdiff(c("scenario_id", "scenario_label", "n", "estimator",
                         "estimator_variant", "truncation_policy", "truncation_target"),
                       names(attempts))) {
      attempts[[cc]] <- NA
    }
    have_key <- if (nrow(out)) paste(out$scenario_id, out$estimator, out$estimator_variant, sep = "\r") else character(0)
    need <- attempts[!paste(attempts$scenario_id, attempts$estimator, attempts$estimator_variant, sep = "\r") %in% have_key, , drop = FALSE]
    if (nrow(need)) {
      filler <- data.frame(
        scenario_id = need$scenario_id,
        scenario_label = need$scenario_label,
        n = need$n,
        estimator = need$estimator,
        estimator_variant = need$estimator_variant,
        truncation_enabled = FALSE,
        truncation_policy = ifelse(is.na(need$truncation_policy) | !nzchar(as.character(need$truncation_policy)),
                                   "not_applicable", as.character(need$truncation_policy)),
        truncation_target = ifelse(is.na(need$truncation_target) | !nzchar(as.character(need$truncation_target)),
                                   "not_applicable", as.character(need$truncation_target)),
        truncation_rule = "not_applicable",
        requested_quantile_lower = NA_real_,
        requested_quantile_upper = NA_real_,
        effective_quantile_lower = NA_real_,
        effective_quantile_upper = NA_real_,
        n_diagnostic_rows = 0L,
        n_values_total = 0,
        mean_fraction_truncated_lower = 0,
        mean_fraction_truncated_upper = 0,
        mean_fraction_truncated_total = 0,
        p95_fraction_truncated_total = 0,
        max_fraction_truncated_total = 0,
        mean_raw_value = NA_real_,
        sd_raw_value = NA_real_,
        min_raw_value = NA_real_,
        p01_raw_value = NA_real_,
        median_raw_value = NA_real_,
        p99_raw_value = NA_real_,
        max_raw_value = NA_real_,
        mean_truncated_value = NA_real_,
        sd_truncated_value = NA_real_,
        min_truncated_value = NA_real_,
        max_truncated_value = NA_real_,
        mean_ess_raw = NA_real_,
        mean_ess_truncated = NA_real_,
        min_ess_truncated = NA_real_,
        n_rows_with_any_truncation = 0L,
        fraction_rows_with_any_truncation = 0,
        fixed_bound_truncation_used = FALSE,
        density_ratio_factor_truncation_used = FALSE,
        probability_bounding_is_truncation = FALSE,
        summary_status = "completed",
        failure_class = "no_failure",
        stringsAsFactors = FALSE
      )
      out <- .runner_rbind_fill(list(out, filler))
    }
  }
  if (is.null(out) || !nrow(out)) return(NULL)
  rownames(out) <- NULL
  out
}

.add_effect_truth_bias <- function(estimates_long, truth_effects) {
  if (is.null(estimates_long) || !is.data.frame(estimates_long) || !nrow(estimates_long) ||
      is.null(truth_effects) || !is.data.frame(truth_effects) || !nrow(truth_effects)) {
    return(estimates_long)
  }
  truth_key <- unique(truth_effects[, intersect(c("scenario_id", "estimand", "effect", "truth"), names(truth_effects)), drop = FALSE])
  if (!all(c("scenario_id", "estimand", "effect", "truth") %in% names(truth_key))) return(estimates_long)
  out <- merge(estimates_long, truth_key, by = c("scenario_id", "estimand", "effect"), all.x = TRUE, sort = FALSE)
  out$bias <- suppressWarnings(as.numeric(out$estimate)) - suppressWarnings(as.numeric(out$truth))
  out
}

.add_worldmean_truth_bias <- function(worldmeans_long, truth_worldmeans) {
  if (is.null(worldmeans_long) || !is.data.frame(worldmeans_long) || !nrow(worldmeans_long) ||
      is.null(truth_worldmeans) || !is.data.frame(truth_worldmeans) || !nrow(truth_worldmeans)) {
    return(worldmeans_long)
  }
  world_col <- if ("world" %in% names(worldmeans_long)) "world" else if ("component" %in% names(worldmeans_long)) "component" else NULL
  value_col <- if ("mean_hat" %in% names(worldmeans_long)) "mean_hat" else if ("worldmean_estimate" %in% names(worldmeans_long)) "worldmean_estimate" else NULL
  if (is.null(world_col) || is.null(value_col)) return(worldmeans_long)
  world_cols <- grep("^mu_", names(truth_worldmeans), value = TRUE)
  if (!length(world_cols)) return(worldmeans_long)
  truth_long <- stats::reshape(
    truth_worldmeans[, c("scenario_id", world_cols), drop = FALSE],
    varying = list(world_cols),
    v.names = "worldmean_truth_or_reference",
    timevar = "component",
    times = world_cols,
    direction = "long"
  )
  truth_long$id <- NULL
  rownames(truth_long) <- NULL
  out <- worldmeans_long
  if (!"component" %in% names(out)) out$component <- out[[world_col]]
  if (!"worldmean_estimate" %in% names(out)) out$worldmean_estimate <- out[[value_col]]
  out <- merge(out, unique(truth_long[, c("scenario_id", "component", "worldmean_truth_or_reference"), drop = FALSE]),
               by = c("scenario_id", "component"), all.x = TRUE, sort = FALSE)
  out$worldmean_bias <- suppressWarnings(as.numeric(out$worldmean_estimate)) -
    suppressWarnings(as.numeric(out$worldmean_truth_or_reference))
  out
}

.augment_performance_summary <- function(perf, attempt_status = NULL,
                                         truncation_summary = NULL) {
  if (is.null(perf) || !is.data.frame(perf) || !nrow(perf)) return(perf)
  out <- perf
  if (!is.null(attempt_status) && is.data.frame(attempt_status) && nrow(attempt_status)) {
    runtime <- .summarize_estimator_runtime(attempt_status)
    if (!is.null(runtime) && nrow(runtime)) {
      runtime_keep <- runtime[, intersect(
        c("scenario_id", "estimator", "estimator_variant", "n_attempts", "n_completed",
          "n_failed", "failure_rate", "mean_elapsed_seconds",
          "nonfinite_effect_rate", "nonfinite_worldmean_rate"),
        names(runtime)
      ), drop = FALSE]
      names(runtime_keep)[names(runtime_keep) == "n_attempts"] <- "n_reps_requested"
      names(runtime_keep)[names(runtime_keep) == "n_completed"] <- "n_success"
      out <- merge(out, runtime_keep, by = c("scenario_id", "estimator"), all.x = TRUE, sort = FALSE)
    }
  }
  if (!is.null(truncation_summary) && is.data.frame(truncation_summary) && nrow(truncation_summary)) {
    trunc_keep <- truncation_summary[, intersect(
      c("scenario_id", "estimator", "truncation_policy", "truncation_target",
        "mean_fraction_truncated_total"),
      names(truncation_summary)
    ), drop = FALSE]
    trunc_keep <- unique(trunc_keep)
    out <- merge(out, trunc_keep, by = c("scenario_id", "estimator"), all.x = TRUE, sort = FALSE)
  }
  if ("K" %in% names(out) && !"n_success" %in% names(out)) out$n_success <- out$K
  if ("R_total" %in% names(out) && !"n_reps_requested" %in% names(out)) out$n_reps_requested <- out$R_total
  if (!"n_failed" %in% names(out) && all(c("n_reps_requested", "n_success") %in% names(out))) {
    out$n_failed <- out$n_reps_requested - out$n_success
  }
  if (!"failure_rate" %in% names(out) && all(c("n_reps_requested", "n_failed") %in% names(out))) {
    out$failure_rate <- ifelse(out$n_reps_requested > 0, out$n_failed / out$n_reps_requested, NA_real_)
  }
  if ("Bias" %in% names(out) && !"bias" %in% names(out)) out$bias <- out$Bias
  if ("SD" %in% names(out) && !"empirical_sd" %in% names(out)) out$empirical_sd <- out$SD
  if ("RMSE" %in% names(out) && !"rmse" %in% names(out)) out$rmse <- out$RMSE
  if ("Coverage" %in% names(out) && !"coverage" %in% names(out)) out$coverage <- out$Coverage
  out
}

make_failed_effect_rows <- function(sc, r, est_name, reason) {
  keys <- primary_effect_keys()
  data.frame(
    scenario_id = sc$scenario_id,
    analysis_tier = sc$analysis_tier %||% NA_character_,
    pathway_setting = sc$pathway_setting %||% sc$PM,
    rho_setting = sc$rho_setting %||% NA_character_,
    structure_setting = sc$structure_setting %||% NA_character_,
    Q_model = sc$Q_model %||% NA_character_,
    fold_count = sc$fold_count %||% NA_integer_,
    PM = sc$PM,
    rho = sc$rho,
    rho0 = sc$rho0,
    rho1 = sc$rho1,
    MI_mode = sc$MI_mode,
    n = sc$n,
    rep = r,
    estimator = est_name,
    estimand = keys$estimand,
    effect = keys$effect,
    estimate = NA_real_,
    finite = FALSE,
    lcl = NA_real_,
    ucl = NA_real_,
    nonfinite_reason = reason,
    stringsAsFactors = FALSE
  )
}

# ---- Main simulation loop ----------------------------------------------------

#' Run the full simulation.
#'
#' @param scenarios data.frame from create_scenario_grid()
#' @param params0 base DGP parameter list (see default_dgp_params())
#' @param run_cfg list with at least:
#'   - T, reg_a, reg_as
#'   - B_truth (truth cohort size)
#'   - B_mc (within-dataset MC size for g-comp; used by estimator wrappers)
#'   - seed (base seed)
#'   - DGP_treat_mech (baseline_rct / sequential_rct / observational)
#'   - p_rct (randomization probability for baseline_rct / sequential_rct)
#'   - B_boot (optional, default 0)
#'   - conf_level (optional, default 0.95)
#' @param estimators named list of functions: function(dat_wide, cfg) -> list(means = ..., ...)
#' @param R_reps number of replications per scenario row
#' @return list with:
#'   - truth_worldmeans (wide)
#'   - truth_effects (long)
#'   - estimates_long (long)
#'   - performance (summary)
run_simulation <- function(scenarios, params0, run_cfg, estimators, R_reps,
                           show_progress = TRUE,
                           progress_every_dgp = NULL,
                           progress_every_analysis = NULL,
                           save_replicate_data = FALSE,
                           save_nuisance_fits = FALSE,
                           save_rng_seed_stream = FALSE) {
  .runner_source_if_needed()

  .fmt_hms <- function(seconds) {
    if (is.null(seconds) || !is.finite(seconds) || seconds < 0) return(NA_character_)
    s <- as.integer(round(seconds))
    h <- s %/% 3600L
    m <- (s %% 3600L) %/% 60L
    ss <- s %% 60L
    sprintf("%02d:%02d:%02d", h, m, ss)
  }

  .progress_line <- function(label, done, total, elapsed_sec) {
    pct <- if (total > 0) 100 * done / total else NA_real_
    eta <- if (done > 0) (elapsed_sec / done) * (total - done) else NA_real_
    sprintf("[%s] %d/%d (%.1f%%) | elapsed=%s | ETA=%s",
            label, done, total, pct, .fmt_hms(elapsed_sec), .fmt_hms(eta))
  }

  if (!is.data.frame(scenarios) || nrow(scenarios) == 0) {
    stop("scenarios must be a non-empty data.frame.", call. = FALSE)
  }
  if (!is.list(estimators) || length(estimators) == 0) {
    stop("estimators must be a non-empty named list of functions.", call. = FALSE)
  }
  if (!is.numeric(R_reps) || length(R_reps) != 1 || R_reps <= 0) {
    stop("R_reps must be a positive integer.", call. = FALSE)
  }
  R_reps <- as.integer(R_reps)
  rep_indices <- run_cfg$rep_indices
  if (is.null(rep_indices)) {
    rep_indices <- seq_len(R_reps)
  } else {
    rep_indices <- as.integer(rep_indices)
    if (length(rep_indices) == 0L || any(is.na(rep_indices)) ||
        any(rep_indices < 1L) || any(rep_indices > R_reps) ||
        any(duplicated(rep_indices))) {
      stop("run_cfg$rep_indices must contain unique integers in 1:R_reps.", call. = FALSE)
    }
    rep_indices <- sort(rep_indices)
  }

  # Normalize scenario columns for manuscript-concordant and legacy manifests.
  if (!("PM" %in% names(scenarios)) && "pathway_setting" %in% names(scenarios)) scenarios$PM <- scenarios$pathway_setting
  if (!("MI_mode" %in% names(scenarios)) && "structure_setting" %in% names(scenarios)) {
    scenarios$MI_mode <- ifelse(scenarios$structure_setting == "gamma0-only", "gamma0_only", scenarios$structure_setting)
  }
  if (!("rho0" %in% names(scenarios)) && "rho" %in% names(scenarios)) scenarios$rho0 <- scenarios$rho
  if (!("rho1" %in% names(scenarios)) && "rho" %in% names(scenarios)) scenarios$rho1 <- scenarios$rho
  if (!("rho" %in% names(scenarios)) && all(c("rho0", "rho1") %in% names(scenarios))) scenarios$rho <- scenarios$rho1
  if (!("Q_model" %in% names(scenarios))) scenarios$Q_model <- "correct"
  if (!("fold_count" %in% names(scenarios))) scenarios$fold_count <- 2L

  # Ensure required columns exist.
  need_cols <- c("scenario_id", "n", "PM", "rho0", "rho1", "MI_mode", "Q_model")
  miss_cols <- setdiff(need_cols, names(scenarios))
  if (length(miss_cols) > 0) {
    stop("scenarios is missing required columns: ", paste(miss_cols, collapse = ", "), call. = FALSE)
  }

  # Defaults
  if (is.null(run_cfg$B_boot)) run_cfg$B_boot <- 0
  if (is.null(run_cfg$conf_level)) run_cfg$conf_level <- 0.95
  skip_truth <- isTRUE(run_cfg$skip_truth)
  skip_performance <- isTRUE(run_cfg$skip_performance)
  if (skip_truth && !skip_performance) {
    stop("run_cfg$skip_truth requires run_cfg$skip_performance = TRUE.", call. = FALSE)
  }

  # Consistency check: p_rct used in DGP vs estimators
  # - DGP uses params0$treatment$p_rct
  # - IPW/LTMLE-exact use run_cfg$p_rct
  if (!is.null(params0$treatment$p_rct) && !is.null(run_cfg$p_rct)) {
    if (!isTRUE(all.equal(params0$treatment$p_rct, run_cfg$p_rct))) {
      stop("p_rct mismatch: params0$treatment$p_rct != run_cfg$p_rct", call. = FALSE)
    }
  }


  # ---- Diagnostics configuration (optional) ---------------------------------
  diag_cfg <- .diag_default_cfg(run_cfg)
  diag_enabled <- isTRUE(diag_cfg$enabled)

  diag_ltmle_exact_score_equation_rows <- list()
  diag_ltmle_exact_run_rows <- list()
  diag_ltmle_exact_fold_rows <- list()
  diag_ltmle_exact_component_registry_rows <- list()
  diag_ltmle_exact_factor_task_rows <- list()
  diag_ltmle_exact_component_eif_rows <- list()
  diag_ltmle_exact_component_eif_term_rows <- list()
  diag_truncation_rows <- list()
  diag_attempt_rows <- list()
  diag_msm_rows <- list()
  diag_gcomp_rows <- list()
  diag_fail_rows <- list()
  diag_full_files <- character(0)

  if (diag_enabled) {
    dir.create(diag_cfg$dir, showWarnings = FALSE, recursive = TRUE)
    if (.diag_should_save_full(diag_cfg, is_outlier = FALSE) ||
        .diag_should_save_full(diag_cfg, is_outlier = TRUE)) {
      dir.create(file.path(diag_cfg$dir, "full"), showWarnings = FALSE, recursive = TRUE)
    }
  }

  seed_base <- if (!is.null(run_cfg$seed)) as.integer(run_cfg$seed) else 1L

  # Progress configuration (auto-scale if not provided)
  n_scen_rows <- nrow(scenarios)
  n_estimators <- length(estimators)
  total_dgp <- n_scen_rows * length(rep_indices)
  total_analysis <- total_dgp * n_estimators

  if (is.null(progress_every_dgp)) {
    progress_every_dgp <- max(1L, as.integer(floor(total_dgp / 10L)))
  }
  if (is.null(progress_every_analysis)) {
    progress_every_analysis <- max(1L, as.integer(floor(total_analysis / 10L)))
  }

  t_start <- proc.time()[["elapsed"]]
  if (isTRUE(show_progress)) {
    cat(sprintf("[Simulation] Started at %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    cat(sprintf("[Simulation] Scenario rows=%d | R_reps=%d | Active reps=%s | Estimators=%d\n",
                n_scen_rows, R_reps, paste(rep_indices, collapse = ","), n_estimators))
    cat(sprintf("[Simulation] Total DGP draws=%d | Total estimator runs=%d\n\n", total_dgp, total_analysis))
  }

  # ---- 1) Truth cache (computed once per unique DGP scenario) ---------------
  scenarios_unique <- unique(scenarios)

  t_truth0 <- proc.time()[["elapsed"]]
  if (skip_truth) {
    truth_obj <- list(truth_worldmeans = NULL, truth_effects = NULL)
    truth_elapsed <- 0.0
    if (isTRUE(show_progress)) {
      cat("[Truth] Skipped for this run_cfg; performance must be computed externally.\n\n")
    }
  } else {
    truth_obj <- .compute_truth_for_unique_scenarios(scenarios_unique, params0, run_cfg)
    truth_elapsed <- proc.time()[["elapsed"]] - t_truth0

    if (isTRUE(show_progress)) {
      n_truth <- if (!is.null(truth_obj$truth_worldmeans) && "truth_id" %in% names(truth_obj$truth_worldmeans)) {
        length(unique(truth_obj$truth_worldmeans$truth_id))
      } else {
        nrow(scenarios_unique)
      }
      cat(sprintf("[Truth] Computed for %d unique DGP scenarios | elapsed=%s\n\n", n_truth, .fmt_hms(truth_elapsed)))
    }
  }

  # ---- 2) Replicate loop ----------------------------------------------------
  est_rows <- list()
  row_idx <- 0L

  replicate_data_dir <- run_cfg$replicate_data_dir %||% "replicate_data"
  save_any_replicate_artifact <- isTRUE(save_replicate_data) ||
    isTRUE(save_nuisance_fits) ||
    isTRUE(save_rng_seed_stream)
  save_rds_quiet <- function(object, path) {
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    saveRDS(object, path, compress = "xz")
    invisible(path)
  }
  nuisance_file_name <- function(est_name) {
    if (grepl("^gcomp($|_)", est_name)) return("nuisance_fits_gcomp.rds")
    if (.is_ltmle_output_name(est_name)) return("nuisance_fits_ltmle_exact.rds")
    paste0("nuisance_fits_", .diag_safe_id(est_name), ".rds")
  }

  means_rows <- list()
  mean_idx <- 0L

  # Progress counters + timers (separate DGP vs analysis)
  dgp_done <- 0L
  analysis_done <- 0L
  dgp_elapsed <- 0.0
  analysis_elapsed <- 0.0

  for (s in seq_len(nrow(scenarios))) {
    sc <- scenarios[s, , drop = FALSE]

    if (isTRUE(show_progress)) {
      cat(sprintf("[Scenario %d/%d] scenario_id=%s | n=%d | PM=%s | rho=(%s,%s) | MI_mode=%s | Q_model=%s\n",
                  s, nrow(scenarios), sc$scenario_id, sc$n,
                  as.character(sc$PM), format(sc$rho0, digits = 3), format(sc$rho1, digits = 3),
                  as.character(sc$MI_mode), as.character(sc$Q_model)))
    }

    # Apply scenario switches to DGP params for observed data generation.
    pm_mult <- NULL
    if (!is.null(run_cfg$PM_MAP)) {
      pm_mult <- run_cfg$PM_MAP[[as.character(sc$PM)]]
    }
    params_sc <- apply_scenario_to_params(
      params0,
      PM = pm_mult,
      MI_mode = sc$MI_mode,
      rho0 = sc$rho0,
      rho1 = sc$rho1
    )

    for (r in rep_indices) {
      # Seed per (scenario row, replication)
      seed_rep <- seed_base + 100000L * s + r
      set.seed(seed_rep)
      rng_seed_before_dgp <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        get(".Random.seed", envir = .GlobalEnv)
      } else {
        NULL
      }

      rep_artifact_dir <- if (save_any_replicate_artifact) {
        file.path(
          replicate_data_dir,
          paste0("scenario_", .diag_safe_id(sc$scenario_id)),
          paste0("rep_", r)
        )
      } else {
        NULL
      }

      # --- DGP ---------------------------------------------------------------
      t_dgp0 <- proc.time()[["elapsed"]]
      dat_wide <- simulate_dgp_wide(
        n = sc$n,
        T = run_cfg$T,
        params = params_sc,
        treat_mech = run_cfg$DGP_treat_mech
      )
      if (save_any_replicate_artifact) {
        if (isTRUE(save_replicate_data)) {
          save_rds_quiet(dat_wide, file.path(rep_artifact_dir, "dat_wide.rds"))
          save_rds_quiet(wide_to_long(dat_wide, run_cfg$T), file.path(rep_artifact_dir, "dat_long.rds"))
          save_rds_quiet(params_sc, file.path(rep_artifact_dir, "dgp_params.rds"))
        }
        if (isTRUE(save_rng_seed_stream)) {
          rng_seed_after_dgp <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            get(".Random.seed", envir = .GlobalEnv)
          } else {
            NULL
          }
          save_rds_quiet(
            list(
              seed_rep = seed_rep,
              rng_seed_before_dgp = rng_seed_before_dgp,
              rng_seed_after_dgp = rng_seed_after_dgp
            ),
            file.path(rep_artifact_dir, "rng_seed_stream.rds")
          )
        }
      }
      dgp_elapsed <- dgp_elapsed + (proc.time()[["elapsed"]] - t_dgp0)
      dgp_done <- dgp_done + 1L

      if (isTRUE(show_progress) && (dgp_done %% progress_every_dgp == 0L || dgp_done == total_dgp)) {
        cat(.progress_line("DGP", dgp_done, total_dgp, dgp_elapsed), "\n")
      }

      # Estimator config passed to wrappers
      cfg_rep <- run_cfg
      cfg_rep$seed <- seed_rep
      cfg_rep$Q_model <- sc$Q_model

      # --- Analysis (estimators) --------------------------------------------
      for (est_name in names(estimators)) {
        est_fun <- estimators[[est_name]]

        t_est0 <- proc.time()[["elapsed"]]
        failure_reason <- NULL
        estimator_warnings <- character()
        estimator_error <- NULL
        res <- tryCatch(
          withCallingHandlers(
            est_fun(dat_wide, cfg_rep),
            warning = function(w) {
              estimator_warnings <<- c(estimator_warnings, conditionMessage(w))
              invokeRestart("muffleWarning")
            }
          ),
          error = function(e) {
          msg <- sprintf(
            "Estimator '%s' failed (scenario_id=%s, n=%s, rep=%s): %s",
            est_name, sc$scenario_id, sc$n, r, conditionMessage(e)
          )

          # Production simulations count failed estimator attempts as nonfinite
          # primary-effect rows. Set SIM_STOP_ON_ESTIMATOR_FAILURE=true only for
          # debugging.
          stop_on_failure <- isTRUE(run_cfg$stop_on_estimator_failure %||% FALSE)
          if (stop_on_failure) {
            stop(msg, call. = FALSE)
          }

          estimator_error <<- e
          failure_reason <<- conditionMessage(e)
          diag_fail_rows <<- c(
            diag_fail_rows,
            list(data.frame(
              scenario_id = sc$scenario_id,
              scenario_label = .scenario_label(sc),
              n = as.integer(sc$n),
              rep = as.integer(r),
              seed = as.integer(seed_rep),
              estimator = est_name,
              estimator_variant = est_name,
              error_class = class(e)[1L],
              error = conditionMessage(e),
              warning_count = length(estimator_warnings),
              warning_messages_summary = .diag_collapse_messages(estimator_warnings),
              stringsAsFactors = FALSE
            ))
          )
          NULL
        })
        est_elapsed <- proc.time()[["elapsed"]] - t_est0
        analysis_elapsed <- analysis_elapsed + est_elapsed
        analysis_done <- analysis_done + 1L

        if (isTRUE(save_nuisance_fits) && !is.null(res)) {
          nuisance_obj <- if (!is.null(res$fits)) {
            res$fits
          } else if (!is.null(res$diagnostics)) {
            res$diagnostics
          } else {
            res
          }
          save_rds_quiet(nuisance_obj, file.path(rep_artifact_dir, nuisance_file_name(est_name)))
        }

        if (isTRUE(show_progress) && (analysis_done %% progress_every_analysis == 0L || analysis_done == total_analysis)) {
          cat(.progress_line("Analysis", analysis_done, total_analysis, analysis_elapsed), "\n")
        }

        if (is.null(res) || is.null(res$means)) {
          reason <- failure_reason %||% "estimator returned NULL or missing means"
          if (is.null(failure_reason)) {
            estimator_error <- simpleError(reason)
            diag_fail_rows[[length(diag_fail_rows) + 1L]] <- data.frame(
              scenario_id = sc$scenario_id,
              scenario_label = .scenario_label(sc),
              n = as.integer(sc$n),
              rep = as.integer(r),
              seed = as.integer(seed_rep),
              estimator = est_name,
              estimator_variant = est_name,
              error_class = class(estimator_error)[1L],
              error = reason,
              warning_count = length(estimator_warnings),
              warning_messages_summary = .diag_collapse_messages(estimator_warnings),
              stringsAsFactors = FALSE
            )
          }
          row_idx <- row_idx + 1L
          failed_eff <- make_failed_effect_rows(sc, r, est_name, reason)
          failed_eff$scenario_label <- .scenario_label(sc)
          failed_eff$seed <- as.integer(seed_rep)
          failed_eff$estimator_variant <- est_name
          failed_eff$truncation_policy <- NA_character_
          failed_eff$truncation_target <- NA_character_
          failed_eff$selected_primary_law_integration_n <- NA_integer_
          failed_eff$status <- "failed"
          failed_eff$failure_class <- class(estimator_error)[1L] %||% "estimator_failed"
          est_rows[[row_idx]] <- failed_eff
          if (diag_enabled) {
            diag_attempt_rows[[length(diag_attempt_rows) + 1L]] <- .attempt_status_row(
              sc = sc,
              rep = r,
              seed = seed_rep,
              est_name = est_name,
              estimator_variant = est_name,
              attempted = TRUE,
              completed = FALSE,
              error_condition = estimator_error,
              warnings = estimator_warnings,
              eff_df = failed_eff,
              means_vec = numeric(0),
              elapsed_seconds = est_elapsed,
              diagnostics_written = TRUE,
              output_rows_written = nrow(failed_eff) > 0L
            )
          }
          next
        }

        estimator_variant <- .estimator_variant_from_result(est_name, res)
        trunc_meta <- .truncation_metadata_from_result(res)

        # Store world-mean estimates (replicate-level; helps diagnose which mu_* broke)
        means_vec <- unlist(res$means)
        if (is.null(names(means_vec))) {
          names(means_vec) <- paste0("mu_", seq_along(means_vec))
        }
        mean_idx <- mean_idx + 1L
        means_rows[[mean_idx]] <- data.frame(
          scenario_id = sc$scenario_id,
          scenario_label = .scenario_label(sc),
          analysis_tier = sc$analysis_tier %||% NA_character_,
          pathway_setting = sc$pathway_setting %||% sc$PM,
          rho_setting = sc$rho_setting %||% NA_character_,
          structure_setting = sc$structure_setting %||% NA_character_,
          PM = sc$PM,
          rho = sc$rho,
          rho0 = sc$rho0,
          rho1 = sc$rho1,
          MI_mode = sc$MI_mode,
          n = sc$n,
          fold_count = sc$fold_count %||% NA_integer_,
          rep = r,
          seed = as.integer(seed_rep),
          estimator = est_name,
          estimator_variant = estimator_variant,
          truncation_policy = trunc_meta$truncation_policy,
          truncation_target = trunc_meta$truncation_target,
          selected_primary_law_integration_n = trunc_meta$selected_primary_law_integration_n,
          status = "completed",
          failure_class = "no_failure",
          learner = res$learner %||% res$metadata$learner %||% NA_character_,
          Q_model = sc$Q_model %||% res$Q_model %||% NA_character_,
          world = names(means_vec),
          mean_hat = as.numeric(means_vec),
          stringsAsFactors = FALSE
        )

        eff_df_all <- flatten_effects(compute_estimands_from_means(res$means))
        eff_df <- filter_primary_nonratio_effects(eff_df_all)

        # Optional bootstrap percentile CI
        if (is.numeric(run_cfg$B_boot) && run_cfg$B_boot > 0) {
          ci_df <- .bootstrap_effect_cis(
            dat_wide = dat_wide,
            cfg = cfg_rep,
            estimator_fun = est_fun,
            B_boot = run_cfg$B_boot,
            conf_level = run_cfg$conf_level,
            seed = seed_rep + 9999L
          )
          if (!is.null(ci_df)) {
            eff_df <- merge(eff_df, ci_df, by = c("estimand", "effect"), all.x = TRUE, sort = FALSE)
          } else {
            eff_df$lcl <- NA_real_
            eff_df$ucl <- NA_real_
          }
        }


                # --- Diagnostics (optional) ---------------------------------------------
        if (diag_enabled) {
          ltmle_exact_score_equations <- NULL
          ltmle_exact_run <- NULL

          truncation_diagnostics <- if (!is.null(res$diagnostics)) {
            res$diagnostics$truncation_diagnostics
          } else {
            NULL
          }
          if (!is.null(truncation_diagnostics) &&
              is.data.frame(truncation_diagnostics) &&
              nrow(truncation_diagnostics) > 0L) {
            truncation_diagnostics$scenario_id <- sc$scenario_id
            truncation_diagnostics$n <- sc$n
            truncation_diagnostics$rep <- r
            truncation_diagnostics$estimator_name <- est_name
            truncation_diagnostics$estimator_variant <- estimator_variant
            diag_truncation_rows[[length(diag_truncation_rows) + 1L]] <- truncation_diagnostics
          }

          # LTMLE-exact diagnostics supplied by the full estimator.
          if (.is_ltmle_output_name(est_name)) {
            if (!is.null(res$diagnostics) &&
                (!is.null(res$diagnostics$score_diagnostics) ||
                 !is.null(res$diagnostics$ltmle_exact_run_summary) ||
                 !is.null(res$diagnostics$ltmle_exact_component_law_registry) ||
                 !is.null(res$diagnostics$ltmle_exact_factor_tasks))) {
              ltmle_exact_score_equations <- res$diagnostics$score_diagnostics
              if (!is.null(ltmle_exact_score_equations) && nrow(ltmle_exact_score_equations) > 0) {
                ltmle_exact_score_equations$scenario_id <- sc$scenario_id
                ltmle_exact_score_equations$n <- sc$n
                ltmle_exact_score_equations$rep <- r
                ltmle_exact_score_equations$estimator <- est_name
              }
              ltmle_exact_run <- if (!is.null(res$diagnostics$ltmle_exact_run_summary)) res$diagnostics$ltmle_exact_run_summary else NULL
              if (!is.null(ltmle_exact_run) && nrow(ltmle_exact_run) > 0) {
                ltmle_exact_run$scenario_id <- sc$scenario_id
                ltmle_exact_run$n <- sc$n
                ltmle_exact_run$rep <- r
                ltmle_exact_run$estimator <- est_name
                diag_ltmle_exact_run_rows[[length(diag_ltmle_exact_run_rows) + 1L]] <- ltmle_exact_run
              }
              ltmle_exact_fold <- if (!is.null(res$diagnostics$fold_summary)) res$diagnostics$fold_summary else NULL
              if (!is.null(ltmle_exact_fold) && nrow(ltmle_exact_fold) > 0) {
                ltmle_exact_fold$scenario_id <- sc$scenario_id
                ltmle_exact_fold$n <- sc$n
                ltmle_exact_fold$rep <- r
                ltmle_exact_fold$estimator <- est_name
                diag_ltmle_exact_fold_rows[[length(diag_ltmle_exact_fold_rows) + 1L]] <- ltmle_exact_fold
              }
              ltmle_exact_component_registry <- if (!is.null(res$diagnostics$ltmle_exact_component_law_registry)) res$diagnostics$ltmle_exact_component_law_registry else NULL
              if (!is.null(ltmle_exact_component_registry) && nrow(ltmle_exact_component_registry) > 0) {
                ltmle_exact_component_registry$scenario_id <- sc$scenario_id
                ltmle_exact_component_registry$n <- sc$n
                ltmle_exact_component_registry$rep <- r
                ltmle_exact_component_registry$estimator <- est_name
                diag_ltmle_exact_component_registry_rows[[length(diag_ltmle_exact_component_registry_rows) + 1L]] <- ltmle_exact_component_registry
              }
              ltmle_exact_factor_tasks <- if (!is.null(res$diagnostics$ltmle_exact_factor_tasks)) res$diagnostics$ltmle_exact_factor_tasks else NULL
              if (!is.null(ltmle_exact_factor_tasks) && nrow(ltmle_exact_factor_tasks) > 0) {
                ltmle_exact_factor_tasks$scenario_id <- sc$scenario_id
                ltmle_exact_factor_tasks$n <- sc$n
                ltmle_exact_factor_tasks$rep <- r
                ltmle_exact_factor_tasks$estimator <- est_name
                diag_ltmle_exact_factor_task_rows[[length(diag_ltmle_exact_factor_task_rows) + 1L]] <- ltmle_exact_factor_tasks
              }
              ltmle_exact_component_eif <- if (!is.null(res$diagnostics$component_eif_summary)) res$diagnostics$component_eif_summary else NULL
              if (!is.null(ltmle_exact_component_eif) && nrow(ltmle_exact_component_eif) > 0) {
                ltmle_exact_component_eif$scenario_id <- sc$scenario_id
                ltmle_exact_component_eif$n <- sc$n
                ltmle_exact_component_eif$rep <- r
                ltmle_exact_component_eif$estimator <- est_name
                diag_ltmle_exact_component_eif_rows[[length(diag_ltmle_exact_component_eif_rows) + 1L]] <- ltmle_exact_component_eif
              }
              ltmle_exact_component_eif_terms <- if (!is.null(res$diagnostics$component_eif_terms)) res$diagnostics$component_eif_terms else NULL
              if (!is.null(ltmle_exact_component_eif_terms) && nrow(ltmle_exact_component_eif_terms) > 0) {
                ltmle_exact_component_eif_terms$scenario_id <- sc$scenario_id
                ltmle_exact_component_eif_terms$n <- sc$n
                ltmle_exact_component_eif_terms$rep <- r
                ltmle_exact_component_eif_terms$estimator <- est_name
                diag_ltmle_exact_component_eif_term_rows[[length(diag_ltmle_exact_component_eif_term_rows) + 1L]] <- ltmle_exact_component_eif_terms
              }
            } else {
              ltmle_exact_score_equations <- .diag_ltmle_exact_score_equations(res, sc, r, est_name)
              if (!is.null(ltmle_exact_score_equations)) {
                ltmle_exact_run <- .diag_ltmle_exact_run_summary(ltmle_exact_score_equations, sc, r, est_name, res = res)
                if (!is.null(ltmle_exact_run)) {
                  diag_ltmle_exact_run_rows[[length(diag_ltmle_exact_run_rows) + 1L]] <- ltmle_exact_run
                }
              }
            }
          }

          # MSM weights (end-of-follow-up raw and truncated stabilized weights).
          if (grepl("^msm_ipw($|_)", est_name) && !is.null(res$diagnostics)) {
            if (!is.null(res$diagnostics$weights_end)) {
              wm_end <- res$diagnostics$weights_end
              keep_msm_weights <- c(
                "w_out_nat_raw", "w_out_nat_trunc",
                "w_out_int_raw", "w_out_int_trunc"
              )
              wm_end <- wm_end[wm_end$weight %in% keep_msm_weights, , drop = FALSE]
              wm_end$scenario_id <- sc$scenario_id
              wm_end$n <- sc$n
              wm_end$rep <- r
              wm_end$estimator <- est_name
              diag_msm_rows[[length(diag_msm_rows) + 1L]] <- wm_end
            }
          }

          # G-computation MCSE (minimal).
          if (grepl("^gcomp($|_)", est_name) && !is.null(res$means_mcse)) {
            mn <- names(res$means_mcse)
            dg <- data.frame(
              scenario_id = sc$scenario_id,
              n = sc$n,
              rep = r,
              estimator = est_name,
              mean_name = mn,
              mcse = as.numeric(unlist(res$means_mcse[mn])),
              stringsAsFactors = FALSE
            )
            diag_gcomp_rows[[length(diag_gcomp_rows) + 1L]] <- dg
          }

          # Outlier screen and optional outputs.
          outlier_flag <- .diag_is_outlier(
            diag_cfg,
            eff_df,
            ltmle_exact_score_equations = ltmle_exact_score_equations
          )

          if (.is_ltmle_output_name(est_name) && !is.null(ltmle_exact_score_equations)) {
            if (identical(diag_cfg$save_score_equations, "all") ||
                (identical(diag_cfg$save_score_equations, "outlier") && outlier_flag)) {
              diag_ltmle_exact_score_equation_rows[[length(diag_ltmle_exact_score_equation_rows) + 1L]] <-
                ltmle_exact_score_equations
            }
          }

          if (.diag_should_save_full(diag_cfg, outlier_flag)) {
            obj_small <- list(
              seed = cfg_rep$seed %||% NA_integer_,
              scenario_id = sc$scenario_id,
              n = sc$n,
              rep = r,
              estimator = est_name,
              effects_long = eff_df,
              ltmle_exact_run_summary = ltmle_exact_run %||% NULL
            )

            fp <- .diag_save_rds(obj_small,
                                 diag_dir_full = file.path(diag_cfg$dir, "full"),
                                 sc = sc, rep = r, est_name = est_name)
            diag_full_files <- c(diag_full_files, fp)
          }
        }

        row_idx <- row_idx + 1L
        estimate_num <- suppressWarnings(as.numeric(eff_df$estimate))
        finite_num <- is.finite(estimate_num)
        est_rows[[row_idx]] <- data.frame(
          scenario_id = sc$scenario_id,
          scenario_label = .scenario_label(sc),
          analysis_tier = sc$analysis_tier %||% NA_character_,
          pathway_setting = sc$pathway_setting %||% sc$PM,
          rho_setting = sc$rho_setting %||% NA_character_,
          structure_setting = sc$structure_setting %||% NA_character_,
          Q_model = sc$Q_model %||% NA_character_,
          fold_count = sc$fold_count %||% NA_integer_,
          PM = sc$PM,
          rho = sc$rho,
          rho0 = sc$rho0,
          rho1 = sc$rho1,
          MI_mode = sc$MI_mode,
          n = sc$n,
          rep = r,
          seed = as.integer(seed_rep),
          estimator = est_name,
          estimator_variant = estimator_variant,
          truncation_policy = trunc_meta$truncation_policy,
          truncation_target = trunc_meta$truncation_target,
          selected_primary_law_integration_n = trunc_meta$selected_primary_law_integration_n,
          status = "completed",
          failure_class = "no_failure",
          estimand = eff_df$estimand,
          effect = eff_df$effect,
          estimate = estimate_num,
          finite = finite_num,
          lcl = if (!is.null(eff_df$lcl)) eff_df$lcl else NA_real_,
          ucl = if (!is.null(eff_df$ucl)) eff_df$ucl else NA_real_,
          nonfinite_reason = ifelse(
            finite_num,
            NA_character_,
            "nonfinite numeric estimate"
          ),
          stringsAsFactors = FALSE
        )
        if (diag_enabled) {
          diag_attempt_rows[[length(diag_attempt_rows) + 1L]] <- .attempt_status_row(
            sc = sc,
            rep = r,
            seed = seed_rep,
            est_name = est_name,
            estimator_variant = estimator_variant,
            attempted = TRUE,
            completed = TRUE,
            error_condition = NULL,
            warnings = estimator_warnings,
            eff_df = eff_df,
            means_vec = means_vec,
            elapsed_seconds = est_elapsed,
            diagnostics_written = !is.null(res$diagnostics),
            output_rows_written = TRUE,
            truncation_policy = trunc_meta$truncation_policy,
            truncation_target = trunc_meta$truncation_target
          )
        }
      }
    }

    if (isTRUE(show_progress)) cat("\n")
  }

  estimates_long <- if (length(est_rows) > 0) do.call(rbind, est_rows) else {
    data.frame(
      scenario_id = character(0),
      scenario_label = character(0),
      analysis_tier = character(0),
      pathway_setting = character(0),
      rho_setting = character(0),
      structure_setting = character(0),
      Q_model = character(0),
      fold_count = integer(0),
      PM = character(0),
      rho = numeric(0),
      rho0 = numeric(0),
      rho1 = numeric(0),
      MI_mode = character(0),
      n = integer(0),
      rep = integer(0),
      seed = integer(0),
      estimator = character(0),
      estimator_variant = character(0),
      truncation_policy = character(0),
      truncation_target = character(0),
      selected_primary_law_integration_n = integer(0),
      status = character(0),
      failure_class = character(0),
      estimand = character(0),
      effect = character(0),
      estimate = numeric(0),
      finite = logical(0),
      lcl = numeric(0),
      ucl = numeric(0),
      nonfinite_reason = character(0),
      stringsAsFactors = FALSE
    )
  }

  worldmeans_estimates_long <- if (length(means_rows) > 0) do.call(rbind, means_rows) else {
    data.frame(
      scenario_id = character(0),
      scenario_label = character(0),
      analysis_tier = character(0),
      pathway_setting = character(0),
      rho_setting = character(0),
      structure_setting = character(0),
      PM = character(0),
      rho = numeric(0),
      rho0 = numeric(0),
      rho1 = numeric(0),
      MI_mode = character(0),
      n = integer(0),
      fold_count = integer(0),
      rep = integer(0),
      seed = integer(0),
      estimator = character(0),
      estimator_variant = character(0),
      truncation_policy = character(0),
      truncation_target = character(0),
      selected_primary_law_integration_n = integer(0),
      status = character(0),
      failure_class = character(0),
      learner = character(0),
      Q_model = character(0),
      world = character(0),
      mean_hat = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  estimator_attempt_status <- if (length(diag_attempt_rows)) {
    .runner_rbind_fill(diag_attempt_rows)
  } else {
    NULL
  }
  truncation_diagnostics_all <- .runner_rbind_fill(diag_truncation_rows)
  estimator_runtime_summary <- .summarize_estimator_runtime(estimator_attempt_status)
  truncation_diagnostics_summary <- .summarize_truncation_diagnostics(
    truncation_diagnostics_all,
    attempt_status = estimator_attempt_status
  )

  # ---- 3) Performance summary ----------------------------------------------
  if (skip_performance) {
    perf <- data.frame()
  } else {
    perf <- summarize_performance(estimates_long, truth_obj$truth_effects)
    perf <- .augment_performance_summary(
      perf,
      attempt_status = estimator_attempt_status,
      truncation_summary = truncation_diagnostics_summary
    )

    # Add scenario descriptors
    scenario_desc_cols <- intersect(
      c("scenario_id", "scenario_label", "analysis_tier", "pathway_setting", "rho_setting", "structure_setting",
        "Q_model", "fold_count", "R", "B_truth_init", "B_MC_init", "truncation_rule",
        "R_reps", "B_truth", "B_mc", "truth_n_batches", "msm_ipw_trunc",
        "msm_ipw_truncation_enabled", "msm_ipw_truncation_policy",
        "msm_ipw_truncation_quantile_lower", "msm_ipw_truncation_quantile_upper",
        "msm_ipw_truncation_target", "mediator_density_mc_n", "ltmle_exact_density_ratio_mc_n",
        "ltmle_exact_probability_bounds", "ltmle_exact_truncation_enabled",
        "ltmle_exact_truncation_policy", "ltmle_exact_truncation_quantile_lower",
        "ltmle_exact_truncation_quantile_upper", "ltmle_exact_truncation_target",
        "ltmle_exact_score_tolerance",
        "ltmle_exact_component_tolerance", "finite_estimate_rule",
        "PM", "rho", "rho0", "rho1", "MI_mode"),
      names(scenarios)
    )
    perf <- merge(
      perf,
      unique(scenarios[, scenario_desc_cols, drop = FALSE]),
      by = "scenario_id",
      all.x = TRUE,
      sort = FALSE
    )

    # Reorder columns (keep core stats first)
    core_cols <- c("scenario_id", "scenario_label", "analysis_tier", "pathway_setting", "rho_setting", "structure_setting",
                   "Q_model", "fold_count", "PM", "rho", "rho0", "rho1", "MI_mode", "n",
                   "estimator", "estimand", "effect")
    other_cols <- setdiff(names(perf), core_cols)
    perf <- perf[, c(core_cols, other_cols), drop = FALSE]
  }

  total_elapsed <- proc.time()[["elapsed"]] - t_start
  if (isTRUE(show_progress)) {
    cat(sprintf("[Simulation] Finished | total elapsed=%s\n", .fmt_hms(total_elapsed)))
    cat(sprintf("[Simulation] Time split: Truth=%s | DGP=%s | Analysis=%s\n",
                .fmt_hms(truth_elapsed), .fmt_hms(dgp_elapsed), .fmt_hms(analysis_elapsed)))
  }

  diagnostics <- NULL
  if (diag_enabled) {
    diagnostics <- list(
      ltmle_exact_score_equations = if (length(diag_ltmle_exact_score_equation_rows)) do.call(rbind, diag_ltmle_exact_score_equation_rows) else NULL,
      ltmle_exact_run = if (length(diag_ltmle_exact_run_rows)) do.call(rbind, diag_ltmle_exact_run_rows) else NULL,
      ltmle_exact_fold = if (length(diag_ltmle_exact_fold_rows)) do.call(rbind, diag_ltmle_exact_fold_rows) else NULL,
      ltmle_exact_component_law_registry = if (length(diag_ltmle_exact_component_registry_rows)) do.call(rbind, diag_ltmle_exact_component_registry_rows) else NULL,
      ltmle_exact_factor_tasks = if (length(diag_ltmle_exact_factor_task_rows)) do.call(rbind, diag_ltmle_exact_factor_task_rows) else NULL,
      ltmle_exact_component_eif_summary = if (length(diag_ltmle_exact_component_eif_rows)) do.call(rbind, diag_ltmle_exact_component_eif_rows) else NULL,
      ltmle_exact_component_eif_terms = if (length(diag_ltmle_exact_component_eif_term_rows)) do.call(rbind, diag_ltmle_exact_component_eif_term_rows) else NULL,
      estimator_attempt_status = estimator_attempt_status,
      estimator_runtime_summary = estimator_runtime_summary,
      truncation_diagnostics = truncation_diagnostics_all,
      truncation_diagnostics_summary = truncation_diagnostics_summary,
      msm = if (length(diag_msm_rows)) do.call(rbind, diag_msm_rows) else NULL,
      gcomp = if (length(diag_gcomp_rows)) do.call(rbind, diag_gcomp_rows) else NULL,
      failures = if (length(diag_fail_rows)) do.call(rbind, diag_fail_rows) else NULL,
      full_files = if (length(diag_full_files)) diag_full_files else NULL
    )
  }

  list(
    truth_worldmeans = truth_obj$truth_worldmeans,
    truth_effects = truth_obj$truth_effects,
    estimates_long = estimates_long,
    worldmeans_estimates_long = worldmeans_estimates_long,
    performance = perf,
    diagnostics = diagnostics
  )
}
