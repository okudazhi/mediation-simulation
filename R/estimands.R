################################################################################
# estimands.R
#
# Compatibility shim.  The canonical component-mean registry and estimand
# contrasts live in core_targets.R.
################################################################################

if (!exists("component_mean_keys", mode = "function") ||
    !exists("compute_estimands_from_means", mode = "function") ||
    !exists("flatten_effects", mode = "function")) {
  stop("Source R/core_targets.R before R/estimands.R.", call. = FALSE)
}
