################################################################################
# utils.R
#
# Backward-compatible wrapper. Shared helpers live in core_utils.R.
################################################################################

if (!exists("default_varmap", mode = "function") ||
    !exists("wide_to_long", mode = "function") ||
    !exists("order_effects_long", mode = "function")) {
  candidates <- c(
    file.path(getwd(), "core_utils.R"),
    file.path(getwd(), "R", "core_utils.R")
  )
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates)) {
    source(normalizePath(candidates[[1L]], winslash = "/", mustWork = TRUE),
           chdir = TRUE)
  }
}

if (!exists("default_varmap", mode = "function") ||
    !exists("wide_to_long", mode = "function") ||
    !exists("order_effects_long", mode = "function")) {
  stop("core_utils.R must be sourced before utils.R.", call. = FALSE)
}

invisible(TRUE)
