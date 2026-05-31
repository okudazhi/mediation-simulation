################################################################################
# analysis_tables.R
#
# Table builders for the ESAX-DN real-data analysis layer.
################################################################################

format_est_ci <- function(estimate, lower = NA_real_, upper = NA_real_, digits = 3) {
  if (!is.finite(estimate)) return(NA_character_)
  est_txt <- format(round(estimate, digits), nsmall = digits, trim = TRUE)
  if (!is.finite(lower) || !is.finite(upper)) return(est_txt)
  lo_txt <- format(round(lower, digits), nsmall = digits, trim = TRUE)
  hi_txt <- format(round(upper, digits), nsmall = digits, trim = TRUE)
  paste0(est_txt, " [", lo_txt, ", ", hi_txt, "]")
}

.lookup_ci <- function(ci_df, stat) {
  if (is.null(ci_df) || !nrow(ci_df)) return(c(lower = NA_real_, upper = NA_real_))
  ii <- match(stat, ci_df$stat)
  if (is.na(ii)) return(c(lower = NA_real_, upper = NA_real_))
  c(lower = ci_df$lower[ii], upper = ci_df$upper[ii])
}

.table_from_named_vector <- function(x,
                                     ci_df = NULL,
                                     labels = NULL,
                                     digits = 3) {
  x_num <- as.numeric(x)
  names(x_num) <- names(x)
  x <- x_num
  stats <- names(x)
  if (is.null(labels)) labels <- stats
  lab_map <- setNames(as.character(labels), names(labels))
  rows <- lapply(stats, function(st) {
    ci <- .lookup_ci(ci_df, st)
    data.frame(
      stat = st,
      label = lab_map[[st]] %||% st,
      estimate = x[[st]],
      lower = ci[["lower"]],
      upper = ci[["upper"]],
      estimate_ci = format_est_ci(x[[st]], ci[["lower"]], ci[["upper"]], digits = digits),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

make_week52_fixed_table <- function(result, digits = 3) {
  x <- unlist(result$effects$fixed)
  ord <- c("TE_Y", "IDE_Y", "IIE_Y1", "IIE_Y2", "R_Y")
  x <- x[ord[ord %in% names(x)]]
  labels <- c(
    TE_Y = "Total effect (TE_Y)",
    IDE_Y = "Yamamuro separate-draw direct effect (IDE_Y)",
    IIE_Y1 = "Yamamuro separate-draw indirect effect via SBP (IIE_Y^1)",
    IIE_Y2 = "Yamamuro separate-draw indirect effect via eGFR (IIE_Y^2)",
    R_Y = "Residual decomposition difference (R_Y)"
  )
  ci_df <- if (!is.null(result$bootstrap)) result$bootstrap$fixed_ci else NULL
  .table_from_named_vector(x, ci_df = ci_df, labels = labels, digits = digits)
}

make_joint_draw_supp_table <- function(result, digits = 3) {
  x <- unlist(result$effects$joint)
  ord <- c("ITE_J", "IDE_J", "PSE1", "PSE2", "MI", "xi")
  x <- x[ord[ord %in% names(x)]]
  labels <- c(
    ITE_J = "Joint-draw treatment-history contrast (ITE_J)",
    IDE_J = "Joint-draw direct effect (IDE_J)",
    PSE1 = "Path-specific effect via SBP (PSE^1)",
    PSE2 = "Path-specific effect via eGFR (PSE^2)",
    MI = "Mediator interdependence (MI)",
    xi = "Treatment-mediator dependence gap (xi)"
  )
  ci_df <- if (!is.null(result$bootstrap)) result$bootstrap$joint_ci else NULL
  .table_from_named_vector(x, ci_df = ci_df, labels = labels, digits = digits)
}

make_component_means_table <- function(result, digits = 3) {
  x <- unlist(result$means)
  ci_df <- if (!is.null(result$bootstrap)) result$bootstrap$means_ci else NULL
  .table_from_named_vector(x, ci_df = ci_df, labels = setNames(names(x), names(x)), digits = digits)
}

make_time_resolved_table <- function(result, digits = 3) {
  tab <- result$effects$time_resolved
  if (is.null(tab) || !nrow(tab)) return(data.frame())
  long <- do.call(rbind, lapply(setdiff(names(tab), "cut_t"), function(nm) {
    data.frame(
      cut_t = tab$cut_t,
      component = nm,
      estimate = tab[[nm]],
      stringsAsFactors = FALSE
    )
  }))
  if (!is.null(result$bootstrap) && !is.null(result$bootstrap$time_resolved_ci) &&
      nrow(result$bootstrap$time_resolved_ci)) {
    long <- merge(long, result$bootstrap$time_resolved_ci,
                  by = c("cut_t", "component"), all.x = TRUE, sort = FALSE)
    long$estimate_ci <- mapply(format_est_ci,
                               estimate = long$estimate,
                               lower = long$lower,
                               upper = long$upper,
                               MoreArgs = list(digits = digits))
  } else {
    long$lower <- NA_real_
    long$upper <- NA_real_
    long$estimate_ci <- vapply(long$estimate, format_est_ci, character(1), digits = digits)
  }
  long
}

.write_tex_fragment <- function(tab, path) {
  if (is.null(tab) || !nrow(tab)) return(invisible(path))
  lines <- c("\\begin{tabular}{ll}",
             "\\hline",
             "Quantity & Estimate [95\\% CI] \\\\",
             "\\hline")
  for (ii in seq_len(nrow(tab))) {
    lab <- gsub("_", "\\\\_", tab$label[ii], fixed = TRUE)
    est <- gsub("_", "\\\\_", tab$estimate_ci[ii], fixed = TRUE)
    lines <- c(lines, paste0(lab, " & ", est, " \\\\"))
  }
  lines <- c(lines, "\\hline", "\\end{tabular}")
  write_lines_safe(lines, path)
  invisible(path)
}

write_analysis_tables <- function(result,
                                  outdir,
                                  prefix = NULL,
                                  digits = 3,
                                  write_tex = TRUE) {
  ensure_dir(outdir)
  prefix <- prefix %||% paste(result$analysis_name, result$estimator, sep = "_")

  comp <- make_component_means_table(result, digits = digits)
  fixed <- make_week52_fixed_table(result, digits = digits)
  joint <- make_joint_draw_supp_table(result, digits = digits)
  tr <- make_time_resolved_table(result, digits = digits)

  write_csv_safe(comp, file.path(outdir, paste0(prefix, "_component_means.csv")))
  write_csv_safe(fixed, file.path(outdir, paste0(prefix, "_week52_fixed.csv")))
  write_csv_safe(joint, file.path(outdir, paste0(prefix, "_joint_draw_supplement.csv")))
  if (nrow(tr)) write_csv_safe(tr, file.path(outdir, paste0(prefix, "_time_resolved.csv")))

  if (isTRUE(write_tex)) {
    .write_tex_fragment(fixed, file.path(outdir, paste0(prefix, "_week52_fixed.tex")))
    .write_tex_fragment(joint, file.path(outdir, paste0(prefix, "_joint_draw_supplement.tex")))
  }

  invisible(list(component_means = comp, fixed = fixed, joint = joint, time_resolved = tr))
}
