################################################################################
# model_specs.R
#
# Model/nuisance specification registry.
#
# Simulation:
#   Q_model = "correct" / "wrong" is carried in the scenario manifest.
#   The label controls which covariate terms are included in the estimator's
#   nuisance regressions.  The DGP remains fixed by the scenario.
#
# Real data:
#   model_set = "main" / "sensitivity_restricted_covariates" / "sensitivity_*".
#   The terms correct/wrong are intentionally not used for real data.
################################################################################

simulation_model_scenarios <- function() {
  list(
    all_correct = list(A = "correct", C = "correct", L = "correct",
                       M1 = "correct", M2 = "correct", Y = "correct"),
    wrong_Q = list(A = "correct", C = "correct", L = "correct",
                   M1 = "wrong", M2 = "wrong", Y = "wrong"),
    wrong_M = list(A = "correct", C = "correct", L = "correct",
                   M1 = "wrong", M2 = "wrong", Y = "correct"),
    wrong_Y = list(A = "correct", C = "correct", L = "correct",
                   M1 = "correct", M2 = "correct", Y = "wrong")
  )
}

normalize_simulation_model_scenario <- function(model = c("correct", "wrong", "all_correct", "wrong_Q")) {
  model <- match.arg(model)
  if (identical(model, "correct")) return(simulation_model_scenarios()$all_correct)
  if (identical(model, "wrong")) return(simulation_model_scenarios()$wrong_Q)
  simulation_model_scenarios()[[model]]
}

realdata_model_sets <- function() {
  default_esax_model_bank()
}
