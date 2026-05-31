################################################################################
# figures.R
#
# Panel-grid figures for the simulation performance summary.
#
# Each PDF page corresponds to ONE condition:
#   sample size (n) x PM x rho x MI_mode x metric
#
# Layout (per page):
#   - facet rows: estimand (yamamuro, tai)
#   - facet cols: effect panel (aligned across estimands)
#   - within each panel: methods side-by-side on the x-axis
#       * colour: Q-model spec (correct=blue, wrong=red; shared across methods)
#       * shape:  method (MSM+IPW, g-comp, LTMLE (GLM), LTMLE (SL))
#
# This file is designed to be sourced from run_simulation.R.
#
# Main entry:
#   create_figures(performance_out, out_pdf = "performance_figures.pdf")
################################################################################

# ---- Package attach (lazy) --------------------------------------------------
.templateB_attach_pkgs <- function() {
  pkgs <- c("dplyr", "tidyr", "ggplot2")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing packages required for figures: ", paste(missing, collapse = ", "),
      "\nInstall them and re-run (e.g., install.packages(c(\"dplyr\",\"tidyr\",\"ggplot2\"))).",
      call. = FALSE
    )
  }

  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
  })

  invisible(TRUE)
}

# ---- Defaults ---------------------------------------------------------------

# Internal estimator IDs -> display labels (used for parsing)
ESTIMATOR_LABEL_MAP_TEMPLATEB <- c(
  msm_ipw = "MSM+IPW",
  msm_ipw_Q_wrong = "MSM+IPW (Q wrong)",
  msm_ipw_Q_correct = "MSM+IPW (Q correct)",
  gcomp = "g-comp",
  gcomp_Q_wrong = "g-comp (Q wrong)",
  gcomp_Q_correct = "g-comp (Q correct)",
  ltmle_glm = "LTMLE (GLM)",
  ltmle_glm_Q_wrong = "LTMLE (GLM, Q wrong)",
  ltmle_glm_Q_correct = "LTMLE (GLM, Q correct)",
  ltmle_sl = "LTMLE (SL)",
  ltmle_sl_Q_wrong = "LTMLE (SL, Q wrong)",
  ltmle_sl_Q_correct = "LTMLE (SL, Q correct)",
  ltmle_exact = "LTMLE",
  ltmle_exact_Q_wrong = "LTMLE (Q wrong)",
  ltmle_exact_Q_correct = "LTMLE (Q correct)"
)

# Effects shown by default (core set used in the sample layout)
CORE_EFFECTS_TEMPLATEB <- list(
  yamamuro = c("TE_Y", "IDE_Y", "IIE_Y1", "IIE_Y2", "R_Y"),
  tai = c("ITE_J", "IDE_J", "PSE1", "PSE2", "MI", "xi")
)

# Methods (x-axis)
METHOD_LEVELS_TEMPLATEB <- c("MSM+IPW", "g-comp", "LTMLE (GLM)", "LTMLE (SL)", "Other")

# Shapes for methods (filled): 16=circle, 15=square, 17=triangle, 18=diamond
METHOD_SHAPES_TEMPLATEB <- c(
  "MSM+IPW" = 17,
  "g-comp" = 16,
  "LTMLE (GLM)" = 15,
  "LTMLE (SL)" = 18,
  "Other" = 4
)

# Spec (colour): shared across methods
SPEC_LEVELS_TEMPLATEB <- c("correct", "wrong")
SPEC_COLORS_TEMPLATEB <- c(correct = "blue", wrong = "red")

# Panel columns (aligned across estimands; matches the sample layout)
EFFECT_PANEL_LEVELS_TEMPLATEB <- c(
  "TE_Y/ITE_J",
  "Direct",
  "M1",
  "M2",
  "R_Y/MI",
  "xi"
)

# ---- Helpers ----------------------------------------------------------------

.check_cols <- function(df, cols, name = "data") {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) {
    stop(sprintf("%s is missing columns: %s", name, paste(miss, collapse = ", ")),
         call. = FALSE)
  }
}

.recode_estimator_labels <- function(x, mapping = ESTIMATOR_LABEL_MAP_TEMPLATEB) {
  x_chr <- as.character(x)
  out <- x_chr
  hit <- x_chr %in% names(mapping)
  out[hit] <- unname(mapping[x_chr[hit]])
  out
}

.templateB_method_from_estimator <- function(estimator_label) {
  x <- as.character(estimator_label)

  dplyr::case_when(
    # Accept both recoded labels ("MSM+IPW (Q ...)") and raw IDs ("msm_ipw_...")
    grepl("^MSM\\+IPW", x) | grepl("^msm_ipw", x, ignore.case = TRUE) ~ "MSM+IPW",
    grepl("^g-comp", x) ~ "g-comp",
    grepl("^LTMLE \\(SL", x) | grepl("^ltmle_sl", x, ignore.case = TRUE) ~ "LTMLE (SL)",
    grepl("^LTMLE \\(GLM", x) | grepl("^ltmle_glm", x, ignore.case = TRUE) ~ "LTMLE (GLM)",
    grepl("^LTMLE", x) | grepl("^ltmle", x, ignore.case = TRUE) ~ "LTMLE (GLM)",
    TRUE ~ "Other"
  )
}

.templateB_spec_from_estimator <- function(estimator_label) {
  x <- as.character(estimator_label)

  dplyr::case_when(
    grepl("wrong", x, ignore.case = TRUE) ~ "wrong",
    grepl("correct", x, ignore.case = TRUE) ~ "correct",
    TRUE ~ "correct"  # e.g., MSM+IPW has no Q_model label
  )
}

# Align effect panels across estimands (matches the sample figure).
.templateB_effect_panel <- function(estimand, effect) {
  est <- as.character(estimand)
  eff <- as.character(effect)

  out <- rep(NA_character_, length(eff))

  i_y <- which(est == "yamamuro")
  if (length(i_y) > 0) {
    out[i_y] <- dplyr::case_when(
      eff[i_y] == "TE_Y"   ~ "TE_Y/ITE_J",
      eff[i_y] == "IDE_Y"  ~ "Direct",
      eff[i_y] == "IIE_Y1" ~ "M1",
      eff[i_y] == "IIE_Y2" ~ "M2",
      eff[i_y] == "R_Y"    ~ "R_Y/MI",
      TRUE ~ eff[i_y]
    )
  }

  i_t <- which(est == "tai")
  if (length(i_t) > 0) {
    out[i_t] <- dplyr::case_when(
      eff[i_t] == "ITE_J" ~ "TE_Y/ITE_J",
      eff[i_t] == "IDE_J" ~ "Direct",
      eff[i_t] == "PSE1"  ~ "M1",
      eff[i_t] == "PSE2"  ~ "M2",
      eff[i_t] == "MI"    ~ "R_Y/MI",
      eff[i_t] == "xi"    ~ "xi",
      TRUE ~ eff[i_t]
    )
  }

  out
}

# Convert the wide performance table into a long table with:
#   metric in {Bias, %Bias, RMSE}
#   value, lo, hi (lo/hi used when available for Bias and %Bias)
.templateB_pivot_metrics <- function(df_wide) {
  # Add optional CI columns if missing
  if (!("Bias_lo" %in% names(df_wide))) df_wide$Bias_lo <- NA_real_
  if (!("Bias_hi" %in% names(df_wide))) df_wide$Bias_hi <- NA_real_
  if (!("PercentBias_lo" %in% names(df_wide))) df_wide$PercentBias_lo <- NA_real_
  if (!("PercentBias_hi" %in% names(df_wide))) df_wide$PercentBias_hi <- NA_real_

  base_cols <- setdiff(names(df_wide), c(
    "Bias","Bias_lo","Bias_hi",
    "PercentBias","PercentBias_lo","PercentBias_hi",
    "RMSE"
  ))

  b1 <- df_wide %>%
    dplyr::select(dplyr::all_of(base_cols), Bias, Bias_lo, Bias_hi) %>%
    dplyr::mutate(metric = "Bias", value = Bias, lo = Bias_lo, hi = Bias_hi) %>%
    dplyr::select(-Bias, -Bias_lo, -Bias_hi)

  b2 <- df_wide %>%
    dplyr::select(dplyr::all_of(base_cols), PercentBias, PercentBias_lo, PercentBias_hi) %>%
    dplyr::mutate(metric = "%Bias", value = PercentBias, lo = PercentBias_lo, hi = PercentBias_hi) %>%
    dplyr::select(-PercentBias, -PercentBias_lo, -PercentBias_hi)

  b3 <- df_wide %>%
    dplyr::select(dplyr::all_of(base_cols), RMSE) %>%
    dplyr::mutate(metric = "RMSE", value = RMSE, lo = NA_real_, hi = NA_real_) %>%
    dplyr::select(-RMSE)

  dplyr::bind_rows(b1, b2, b3)
}

# ---- Data prep --------------------------------------------------------------

prep_perf_for_templateB <- function(performance_out,
                                    scenario_df = NULL,
                                    include_all_effects = FALSE,
                                    core_effects = CORE_EFFECTS_TEMPLATEB) {
  .templateB_attach_pkgs()

  req <- c("scenario_id", "n", "estimator", "estimand", "effect", "Bias", "PercentBias", "RMSE")
  .check_cols(performance_out, req, "performance_out")

  df <- performance_out

  # Ensure scenario descriptors exist (MI_mode, rho, PM, optionally Q_model)
  need_scn <- c("MI_mode", "rho", "PM")
  if (!all(need_scn %in% names(df))) {
    if (is.null(scenario_df)) {
      stop("performance_out lacks MI_mode/rho/PM; please provide scenario_df.", call. = FALSE)
    }
    .check_cols(scenario_df, c("scenario_id", need_scn), "scenario_df")
    join_cols <- c("scenario_id", need_scn, intersect(c("Q_model", "rho_setting", "pathway_setting", "structure_setting"), names(scenario_df)))
    df <- df %>%
      dplyr::left_join(
        scenario_df %>% dplyr::select(dplyr::all_of(join_cols)) %>% dplyr::distinct(),
        by = "scenario_id"
      )
  }

  # Estimator parsing (works whether estimator is internal id or already a label)
  df <- df %>%
    dplyr::mutate(
      estimator_label = .recode_estimator_labels(estimator),
      method = .templateB_method_from_estimator(estimator_label),
      spec   = if ("Q_model" %in% names(df)) as.character(Q_model) else .templateB_spec_from_estimator(estimator_label),

      method = factor(method, levels = METHOD_LEVELS_TEMPLATEB),
      spec   = factor(spec, levels = SPEC_LEVELS_TEMPLATEB),

      estimand = factor(estimand, levels = c("yamamuro", "tai")),
      MI_mode = factor(
        MI_mode,
        levels = c("none", "gamma0_only", "full"),
        labels = c("MI: none", "MI: gamma0_only", "MI: full")
      )
    )

  # Effect panel alignment (sample layout)
  df$effect_panel_raw <- .templateB_effect_panel(df$estimand, df$effect)

  if (isTRUE(include_all_effects)) {
    # Keep all effects: map unknown ones to their original effect name.
    df$effect_panel <- ifelse(
      is.na(df$effect_panel_raw),
      as.character(df$effect),
      as.character(df$effect_panel_raw)
    )
    levs <- unique(c(EFFECT_PANEL_LEVELS_TEMPLATEB, sort(unique(df$effect_panel))))
    df$effect_panel <- factor(df$effect_panel, levels = levs)
  } else {
    # Core-effects only: keep the aligned panel order.
    df <- df %>% dplyr::filter(effect %in% unlist(core_effects))
    df$effect_panel <- factor(df$effect_panel_raw, levels = EFFECT_PANEL_LEVELS_TEMPLATEB)
  }

  df$effect_panel_raw <- NULL

  # Metric long table with CI columns (when available)
  df_long <- .templateB_pivot_metrics(df)

  df_long
}

# ---- Plot -------------------------------------------------------------------

make_panel_grid_plot <- function(df_cond,
                                 metric_name = c("Bias", "%Bias", "RMSE"),
                                 show_ci = TRUE,
                                 point_size = 2.0) {
  .templateB_attach_pkgs()
  metric_name <- match.arg(metric_name)

  req <- c("n", "PM", "rho", "MI_mode",
           "estimand", "effect", "effect_panel",
           "method", "spec",
           "metric", "value", "lo", "hi")
  .check_cols(df_cond, req, "df_cond")

  d <- df_cond %>%
    dplyr::filter(metric == metric_name) %>%
    dplyr::filter(!is.na(effect_panel)) %>%
    dplyr::mutate(
      effect_panel = droplevels(effect_panel),
      method = droplevels(method),
      spec = droplevels(spec)
    )

  if (nrow(d) == 0) {
    stop("No rows for metric=", metric_name, " in the selected condition.", call. = FALSE)
  }

  # Condition labels (unique within df_cond by construction)
  n0  <- unique(d$n)[1]
  pm0 <- unique(d$PM)[1]
  rho0 <- unique(d$rho)[1]
  mi0 <- as.character(unique(d$MI_mode)[1])

  rho_chr <- if (is.numeric(rho0)) format(rho0, digits = 3) else as.character(rho0)
  rho_num <- suppressWarnings(as.numeric(rho0))
  corr_chr <- if (is.finite(rho_num) && abs(rho_num) < 1e-12) "corr: no" else "corr: yes"
  main_title <- sprintf("n=%s | PM=%s | rho=%s (%s) | %s | metric=%s",
                        n0, as.character(pm0), rho_chr, corr_chr, mi0, metric_name)

  pos <- ggplot2::position_dodge(width = 0.55)

  p <- ggplot2::ggplot(
    d,
    ggplot2::aes(x = method, y = value, colour = spec, shape = method, group = spec)
  )

  if (isTRUE(show_ci) && metric_name %in% c("Bias", "%Bias")) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lo, ymax = hi),
      position = pos,
      width = 0.15,
      alpha = 0.7,
      linewidth = 0.4,
      na.rm = TRUE
    )
  }

  p <- p +
    ggplot2::geom_point(position = pos, size = point_size, alpha = 0.9, na.rm = TRUE) +
    {if (metric_name %in% c("Bias", "%Bias")) ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35)} +
    ggplot2::facet_grid(rows = ggplot2::vars(estimand), cols = ggplot2::vars(effect_panel), drop = FALSE) +
    ggplot2::scale_colour_manual(values = SPEC_COLORS_TEMPLATEB, drop = FALSE, na.translate = FALSE) +
    ggplot2::scale_shape_manual(values = METHOD_SHAPES_TEMPLATEB, drop = FALSE) +
    ggplot2::labs(
      title = main_title,
      x = NULL,
      y = metric_name,
      colour = "Q model",
      caption = if (isTRUE(show_ci) && metric_name %in% c("Bias", "%Bias")) "Error bars: approx. 95% MCSE interval" else NULL
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 11, face = "bold"),
      plot.caption = ggplot2::element_text(size = 8),
      strip.text = ggplot2::element_text(size = 9),
      axis.text.x = ggplot2::element_text(size = 8, angle = 35, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 8),
      legend.position = "bottom",
      legend.box = "horizontal",
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::guides(shape = "none")

  p
}

# ---- Export -----------------------------------------------------------------

export_templateB_pdf <- function(perf_long,
                                 out_pdf = file.path("output", "performance_figures.pdf"),
                                 width = 20,
                                 height = 7,
                                 show_ci = TRUE) {
  .templateB_attach_pkgs()

  req <- c("n","PM","rho","MI_mode","metric",
           "estimand","effect","effect_panel",
           "method","spec","value","lo","hi")
  .check_cols(perf_long, req, "perf_long")

  # Sorting helpers
  ns <- sort(unique(perf_long$n))

  pm_vals <- unique(perf_long$PM)
  if (all(c("high","low") %in% as.character(pm_vals))) {
    pms <- c("high","low")
  } else {
    pms <- sort(as.character(pm_vals))
  }

  rhos <- sort(unique(perf_long$rho))

  mis <- levels(perf_long$MI_mode)
  if (is.null(mis) || length(mis) == 0) mis <- unique(perf_long$MI_mode)

  metrics <- c("Bias", "%Bias", "RMSE")

  grDevices::pdf(out_pdf, width = width, height = height, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  for (n0 in ns) {
    for (pm0 in pms) {
      for (rho0 in rhos) {
        for (mi0 in mis) {
          d0 <- perf_long %>%
            dplyr::filter(n == n0, PM == pm0, rho == rho0, MI_mode == mi0)

          if (nrow(d0) == 0) next

          for (m0 in metrics) {
            p <- make_panel_grid_plot(d0, metric_name = m0, show_ci = show_ci)
            print(p)
          }
        }
      }
    }
  }

  message("Wrote: ", out_pdf)
  invisible(out_pdf)
}

# ---- Public wrapper ---------------------------------------------------------

create_figures <- function(performance_out,
                           scenario_df = NULL,
                           out_pdf = "performance_figures.pdf",
                           include_all_effects = FALSE,
                           width = 20,
                           height = 7,
                           show_ci = TRUE,
                           ...) {
  .templateB_attach_pkgs()

  # If caller passed only a filename (no directory), save under ./output/
  if (dirname(out_pdf) == ".") {
    out_pdf <- file.path("output", out_pdf)
  }

  out_dir <- dirname(out_pdf)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  perf_long <- prep_perf_for_templateB(
    performance_out = performance_out,
    scenario_df = scenario_df,
    include_all_effects = include_all_effects
  )

  export_templateB_pdf(
    perf_long = perf_long,
    out_pdf = out_pdf,
    width = width,
    height = height,
    show_ci = show_ci
  )
}

# ---- Truth figures ----------------------------------------------------------

#' Create a multi-panel figure of true estimand values across scenarios.
#'
#' This function is intended for the truth table written by run_simulation.R
#' (typically results/simulation_manuscript/truth_effects.csv). It creates
#' one panel per aligned effect target so Yamamuro and Tai truth values can be
#' compared on the same y-axis:
#'   - x-axis: PM and rho
#'   - colour: estimand family
#'   - shape:  MI_mode
#'   - facets: aligned effect targets
#'
#' @param truth_df Optional data.frame with required columns:
#'   PM, rho, MI_mode, estimand, effect, truth.
#' @param truth_csv CSV path used if truth_df is NULL.
#' @param out_pdf Output PDF path.
#' @return Invisibly returns the ggplot object.
create_truth_figures <- function(truth_df = NULL,
                                 truth_csv = file.path("results", "simulation_manuscript", "truth_effects.csv"),
                                 out_pdf = file.path("results", "simulation_manuscript", "truth_figures.pdf")) {

  # Required packages
  pkgs <- c("dplyr", "ggplot2")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
  })

  # Load data
  if (is.null(truth_df)) {
    if (!file.exists(truth_csv)) {
      stop("truth_csv not found: ", truth_csv, call. = FALSE)
    }
    truth_df <- utils::read.csv(truth_csv, stringsAsFactors = FALSE)
  }

  # Output path: if only a filename, write into the manuscript simulation directory.
  if (dirname(out_pdf) %in% c(".", "")) {
    out_pdf <- file.path("results", "simulation_manuscript", out_pdf)
  }
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

  # Required columns
  req <- c("PM", "rho", "MI_mode", "estimand", "effect", "truth")
  miss <- setdiff(req, names(truth_df))
  if (length(miss) > 0) {
    stop("truth_df must contain columns: ", paste(req, collapse = ", "), call. = FALSE)
  }

  # Recode to aligned effect panels.
  df0 <- truth_df %>%
    mutate(
      effect_panel = dplyr::case_when(
        estimand == "yamamuro" & effect == "TE_Y" ~ "Total: TE_Y / ITE_J",
        estimand == "tai" & effect == "ITE_J" ~ "Total: TE_Y / ITE_J",

        estimand == "yamamuro" & effect == "IDE_Y" ~ "Direct: IDE_Y / IDE_J",
        estimand == "tai" & effect == "IDE_J" ~ "Direct: IDE_Y / IDE_J",

        estimand == "yamamuro" & effect == "IIE_Y1" ~ "M1: IIE_Y1 / PSE1",
        estimand == "tai" & effect == "PSE1" ~ "M1: IIE_Y1 / PSE1",

        estimand == "yamamuro" & effect == "IIE_Y2" ~ "M2: IIE_Y2 / PSE2",
        estimand == "tai" & effect == "PSE2" ~ "M2: IIE_Y2 / PSE2",

        estimand == "yamamuro" & effect == "R_Y" ~ "Interdependence: R_Y / MI",
        estimand == "tai" & effect == "MI"  ~ "Interdependence: R_Y / MI",
        estimand == "tai" & effect == "xi" ~ "Tai residual: xi",

        TRUE ~ NA_character_
      ),

      # Fix factor levels
      PM = factor(PM, levels = c("low", "high")),
      rho = as.numeric(rho),
      pm_rho_chr = paste0("PM=", PM, "\n", "rho=", rho),

      MI_mode = factor(MI_mode, levels = c("none", "gamma0_only", "full")),
      estimand = factor(estimand, levels = c("yamamuro", "tai")),
      effect_panel = factor(
        effect_panel,
        levels = c(
          "Total: TE_Y / ITE_J",
          "Direct: IDE_Y / IDE_J",
          "M1: IIE_Y1 / PSE1",
          "M2: IIE_Y2 / PSE2",
          "Interdependence: R_Y / MI",
          "Tai residual: xi"
        )
      )
    )

  axis_levels <- df0 %>%
    distinct(PM, rho, pm_rho_chr) %>%
    arrange(PM, rho) %>%
    pull(pm_rho_chr)
  df0$pm_rho <- factor(df0$pm_rho_chr, levels = axis_levels)

  # Warn and drop unmapped effects (no "Other" panels)
  dropped <- df0 %>% filter(is.na(effect_panel))
  if (nrow(dropped) > 0) {
    drop_tab <- dropped %>%
      dplyr::count(estimand, effect, name = "n") %>%
      dplyr::arrange(estimand, dplyr::desc(n), effect)
    drop_labels <- paste0(drop_tab$estimand, ":", drop_tab$effect, " (n=", drop_tab$n, ")")
    warning(
      sprintf(
        "Dropped %d truth rows with unmapped (estimand, effect): %s",
        nrow(dropped),
        paste(drop_labels, collapse = "; ")
      ),
      call. = FALSE
    )
  }

  df <- df0 %>%
    filter(!is.na(effect_panel)) %>%
    distinct(
      PM, rho, MI_mode, estimand, effect_panel, effect, truth,
      .keep_all = TRUE
    )

  dodge <- position_dodge(width = 0.7)
  p <- ggplot(
    df,
    aes(
      x = pm_rho,
      y = truth,
      color = estimand,
      shape = MI_mode,
      group = interaction(MI_mode, estimand)
    )
  ) +
    geom_hline(yintercept = 0, linewidth = 0.25, color = "grey70") +
    geom_point(position = dodge, size = 2.2, alpha = 0.95) +
    facet_wrap(~ effect_panel, scales = "free_y", nrow = 2) +
    scale_color_manual(
      values = c(yamamuro = "#0072B2", tai = "#D55E00"),
      labels = c(yamamuro = "Yamamuro", tai = "Tai"),
      drop = FALSE
    ) +
    labs(
      x = NULL,
      y = "Truth",
      color = "Estimand",
      shape = "MI_mode",
      title = "Truth by aligned Yamamuro and Tai effects"
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(size = 7),
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 9)
    )

  ggsave(out_pdf, p, width = 14, height = 7, units = "in")
  message("Wrote: ", out_pdf)
  invisible(p)
}
