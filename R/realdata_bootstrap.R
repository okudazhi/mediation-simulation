################################################################################
# realdata_bootstrap.R
#
# Subject-level bootstrap for the real-data analysis layer.
################################################################################

.bootstrap_extract_result <- function(result) {
  means_obj <- result$means %||% result$component_means %||% result$registry_means %||% numeric(0)
  means <- as.numeric(unlist(means_obj))
  names(means) <- names(means_obj)

  fixed_obj <- result$effects$fixed %||% result$fixed_time %||% numeric(0)
  fixed <- as.numeric(unlist(fixed_obj))
  names(fixed) <- names(fixed_obj)

  joint_obj <- result$effects$joint %||% result$joint_effects %||% numeric(0)
  joint <- as.numeric(unlist(joint_obj))
  names(joint) <- names(joint_obj)

  time_resolved <- result$effects$time_resolved %||% result$time_resolved
  if (is.null(time_resolved) || !nrow(time_resolved)) {
    time_resolved <- data.frame()
  }

  list(
    means = means,
    fixed = fixed,
    joint = joint,
    time_resolved = time_resolved
  )
}

.bootstrap_ci_from_matrix <- function(mat, conf_level = 0.95, stat_col = "stat") {
  if (is.null(mat) || !nrow(mat)) return(data.frame())
  alpha <- (1 - conf_level) / 2
  stats <- colnames(mat)
  lo <- apply(mat, 2, function(x) safe_quantile(x, alpha))
  hi <- apply(mat, 2, function(x) safe_quantile(x, 1 - alpha))
  data.frame(
    stat = stats,
    lower = as.numeric(lo),
    upper = as.numeric(hi),
    boot_mean = colMeans(mat, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

.bootstrap_ci_time_resolved <- function(df_long, conf_level = 0.95) {
  if (is.null(df_long) || !nrow(df_long)) return(data.frame())
  alpha <- (1 - conf_level) / 2
  keys <- unique(df_long[, c("cut_t", "component"), drop = FALSE])
  rows <- vector("list", nrow(keys))
  for (ii in seq_len(nrow(keys))) {
    kk <- keys[ii, , drop = FALSE]
    x <- df_long$estimate[df_long$cut_t == kk$cut_t & df_long$component == kk$component]
    rows[[ii]] <- data.frame(
      cut_t = kk$cut_t,
      component = kk$component,
      lower = safe_quantile(x, alpha),
      upper = safe_quantile(x, 1 - alpha),
      boot_mean = mean(x, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

bootstrap_realdata_estimator <- function(prepared,
                                         estimator_fn,
                                         control,
                                         B = 1000L,
                                         conf_level = 0.95,
                                         seed_base = 20260408L,
                                         checkpoint_dir = NULL,
                                         prefix = "bootstrap",
                                         save_every = 25L,
                                         resume = TRUE,
                                         progress = TRUE) {
  B <- as.integer(B)[1L]
  save_every <- max(1L, as.integer(save_every)[1L])
  checkpoint_dir <- checkpoint_dir %||% tempdir()
  ensure_dir(checkpoint_dir)
  checkpoint_file <- file.path(checkpoint_dir, paste0(prefix, "_state.rds"))

  state <- NULL
  if (isTRUE(resume) && file.exists(checkpoint_file)) {
    state <- readRDS(checkpoint_file)
  }

  if (is.null(state)) {
    state <- list(
      next_b = 1L,
      means = list(),
      fixed = list(),
      joint = list(),
      time_resolved = list(),
      failures = data.frame(
        b = integer(0),
        message = character(0),
        stringsAsFactors = FALSE
      )
    )
  }

  if (is.null(state$next_b) || length(state$next_b) != 1L || !is.finite(state$next_b)) {
    .stop("Invalid bootstrap checkpoint: state$next_b must be a finite scalar.")
  }
  state$next_b <- as.integer(state$next_b)[1L]
  if (state$next_b < 1L) {
    .stop("Invalid bootstrap checkpoint: state$next_b must be >= 1.")
  }

  control_boot <- control
  if (!is.null(control$mc_B_boot)) {
    control_boot$mc_B_main <- control$mc_B_boot
  }

  # Do not inherit main-analysis ltmle_exact internal verbose output into
  # bootstrap replicates by default. Outer bootstrap progress is controlled by
  # the `progress` argument.
  if (!is.null(control$bootstrap_verbose)) {
    control_boot$verbose <- isTRUE(control$bootstrap_verbose)
  } else if (!is.null(control$verbose)) {
    control_boot$verbose <- FALSE
  }

  if (state$next_b <= B) {
    for (bb in seq.int(state$next_b, B)) {
      if (isTRUE(progress) && (bb == 1L || bb %% 10L == 0L)) {
        message("[bootstrap] replicate ", bb, " / ", B)
      }

      boot_wide <- tryCatch(
        subject_bootstrap_sample(prepared$wide, id_col = "ID", seed = seed_base + bb),
        error = function(e) e
      )
      if (inherits(boot_wide, "error")) {
        state$failures <- rbind(state$failures,
                                data.frame(b = bb, message = boot_wide$message, stringsAsFactors = FALSE))
        state$next_b <- bb + 1L
        if (bb %% save_every == 0L || bb == B) {
          save_rds_safe(state, checkpoint_file)
        }
        next
      }

      boot_prepared <- prepared_from_wide_and_spec(boot_wide, prepared$spec)
      fit <- tryCatch(
        estimator_fn(boot_prepared, control = control_boot, seed = seed_base + 10000L + bb),
        error = function(e) e
      )
      if (inherits(fit, "error")) {
        state$failures <- rbind(state$failures,
                                data.frame(b = bb, message = fit$message, stringsAsFactors = FALSE))
      } else {
        ext <- .bootstrap_extract_result(fit)
        state$means[[as.character(bb)]] <- ext$means
        state$fixed[[as.character(bb)]] <- ext$fixed
        state$joint[[as.character(bb)]] <- ext$joint
        if (nrow(ext$time_resolved)) {
          td <- ext$time_resolved
          td$b <- bb
          state$time_resolved[[as.character(bb)]] <- td
        }
      }

      state$next_b <- bb + 1L
      if (bb %% save_every == 0L || bb == B) {
        save_rds_safe(state, checkpoint_file)
      }
    }
  } else if (isTRUE(progress)) {
    message("[bootstrap] checkpoint already complete: next_b=", state$next_b, ", B=", B)
  }

  means_mat <- if (length(state$means)) do.call(rbind, lapply(state$means, function(x) matrix(x, nrow = 1, dimnames = list(NULL, names(x))))) else matrix(numeric(0), nrow = 0)
  fixed_mat <- if (length(state$fixed)) do.call(rbind, lapply(state$fixed, function(x) matrix(x, nrow = 1, dimnames = list(NULL, names(x))))) else matrix(numeric(0), nrow = 0)
  joint_mat <- if (length(state$joint)) do.call(rbind, lapply(state$joint, function(x) matrix(x, nrow = 1, dimnames = list(NULL, names(x))))) else matrix(numeric(0), nrow = 0)
  time_df <- if (length(state$time_resolved)) do.call(rbind, state$time_resolved) else data.frame()

  list(
    state = state,
    means_ci = .bootstrap_ci_from_matrix(means_mat, conf_level = conf_level),
    fixed_ci = .bootstrap_ci_from_matrix(fixed_mat, conf_level = conf_level),
    joint_ci = .bootstrap_ci_from_matrix(joint_mat, conf_level = conf_level),
    time_resolved_ci = .bootstrap_ci_time_resolved(time_df, conf_level = conf_level),
    failures = state$failures,
    checkpoint_file = checkpoint_file,
    n_success = length(state$means),
    n_failures = nrow(state$failures)
  )
}
