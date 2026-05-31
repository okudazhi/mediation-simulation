################################################################################
# analysis_figures.R
#
# Plot builders for the ESAX-DN real-data analysis layer.
################################################################################

plot_time_resolved_decomposition <- function(result,
                                             out_file,
                                             main = NULL) {
  tab <- result$effects$time_resolved
  if (is.null(tab) || !nrow(tab)) return(invisible(NULL))

  ensure_dir(dirname(out_file))
  grDevices::pdf(out_file, width = 9, height = 6)
  on.exit(grDevices::dev.off(), add = TRUE)

  x <- tab$cut_t
  ymat <- as.matrix(tab[, c("IDE_Y_le", "IIE_Y1_le", "IIE_Y2_le"), drop = FALSE])
  matplot(x, ymat, type = "b", pch = 1:3, lty = 1:3,
          xlab = "Cut t",
          ylab = "Accumulated effect through t",
          main = main %||% paste(result$analysis_name, "-", result$estimator))
  legend("topleft",
         legend = c("IDE_Y_{<=t}", "IIE_Y^1_{<=t}", "IIE_Y^2_{<=t}"),
         pch = 1:3, lty = 1:3, bty = "n")
  invisible(out_file)
}

plot_weight_diagnostics <- function(result,
                                    out_file,
                                    main = NULL) {
  if (is.null(result$weight_models) || is.null(result$weight_models$subject_weights)) {
    return(invisible(NULL))
  }
  w <- result$weight_models$subject_weights$sw_final
  w <- w[is.finite(w)]
  if (!length(w)) return(invisible(NULL))

  ensure_dir(dirname(out_file))
  grDevices::pdf(out_file, width = 8, height = 6)
  on.exit(grDevices::dev.off(), add = TRUE)

  hist(w,
       breaks = "FD",
       main = main %||% paste(result$analysis_name, "-", result$estimator, "subject weights"),
       xlab = "Final subject weight")
  invisible(out_file)
}

plot_mc_diagnostics <- function(result,
                                out_file,
                                main = NULL) {
  tab <- result$mc_diagnostics
  if (is.null(tab) || !nrow(tab)) return(invisible(NULL))

  ensure_dir(dirname(out_file))
  grDevices::pdf(out_file, width = 10, height = 6)
  on.exit(grDevices::dev.off(), add = TRUE)

  mc_col <- if ("mcse" %in% names(tab)) "mcse" else if ("mc_se" %in% names(tab)) "mc_se" else NULL
  if (is.null(mc_col)) return(invisible(NULL))

  ord <- order(tab[[mc_col]], decreasing = TRUE)
  dotchart(tab[[mc_col]][ord],
           labels = tab$target[ord],
           xlab = "Monte Carlo SE",
           main = main %||% paste(result$analysis_name, "-", result$estimator, "MC diagnostics"))
  invisible(out_file)
}
