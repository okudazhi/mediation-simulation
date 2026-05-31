################################################################################
# nuisance_learners.R
#
# Lightweight GLM/SuperLearner nuisance backend used by real-data and simulation
# wrappers. SuperLearner must be installed when the SuperLearner backend is
# requested.
################################################################################

.nuisance_family_name <- function(family) {
  if (is.character(family)) return(family[1L])
  if (!is.null(family$family)) return(family$family)
  "gaussian"
}

fit_nuisance_model <- function(formula,
                               data,
                               family = stats::gaussian(),
                               engine = c("glm", "superlearner", "sl"),
                               sl_library = c("SL.glm", "SL.mean"),
                               component = "nuisance",
                               weights = NULL) {
  engine <- match.arg(engine)
  if (identical(engine, "sl")) engine <- "superlearner"
  fam_name <- .nuisance_family_name(family)
  mf <- stats::model.frame(formula, data = data, na.action = stats::na.omit)
  if (!nrow(mf)) .stop("No complete cases available for ", component, ".")
  y <- mf[[1L]]
  X <- mf[, -1L, drop = FALSE]
  if (is.null(weights)) {
    obs_w <- rep(1, length(y))
  } else {
    ww0 <- as.numeric(weights)
    if (length(ww0) == nrow(data)) names(ww0) <- rownames(data)
    rrn <- rownames(mf)
    if (!is.null(names(ww0)) && all(rrn %in% names(ww0))) {
      obs_w <- ww0[rrn]
    } else if (length(ww0) == length(y)) {
      obs_w <- ww0
    } else {
      obs_w <- rep(1, length(y))
    }
  }
  if (any(!is.finite(obs_w)) || length(obs_w) != length(y)) obs_w <- rep(1, length(y))
  obs_w[obs_w < 0] <- 0
  if (all(obs_w == 0)) obs_w <- rep(1, length(y))
  mf$.obs_w <- obs_w
  rhs_names <- names(X)

  if (identical(engine, "superlearner")) {
    if (!requireNamespace("SuperLearner", quietly = TRUE)) {
      .stop("SuperLearner requested but package is unavailable.")
    }
    if (!"package:SuperLearner" %in% search()) {
      suppressPackageStartupMessages(base::library(SuperLearner))
    }
  }

  if (identical(engine, "superlearner")) {
    fit <- tryCatch(
      SuperLearner::SuperLearner(
        Y = y,
        X = X,
        family = family,
        SL.library = sl_library,
        verbose = FALSE,
        obsWeights = obs_w
      ),
      error = function(e) e
    )
    if (!inherits(fit, "error")) {
      pred <- tryCatch(as.numeric(SuperLearner::predict.SuperLearner(fit, newdata = X)$pred),
                       error = function(e) rep(mean(as.numeric(y), na.rm = TRUE), length(y)))
      sigma <- if (identical(fam_name, "gaussian")) {
        sqrt(mean((as.numeric(y) - pred)^2, na.rm = TRUE))
      } else {
        NA_real_
      }
      return(structure(
        list(type = "superlearner", fit = fit, formula = formula, family = family,
             rhs_names = rhs_names, fallback = mean(as.numeric(y), na.rm = TRUE),
             sigma = sigma, component = component, sl_library = sl_library),
        class = c("nuisance_model", "list")
      ))
    }
    .stop("SuperLearner fit failed for ", component, ": ", fit$message)
  }

  if (identical(fam_name, "gaussian")) {
    fit <- tryCatch(stats::lm(formula, data = mf, weights = .obs_w), error = function(e) e)
  } else {
    fit <- tryCatch(stats::glm(formula, data = mf, family = family, weights = .obs_w), error = function(e) e)
  }
  if (inherits(fit, "error")) .stop("Failed to fit ", component, ": ", fit$message)
  sigma <- if (identical(fam_name, "gaussian")) model_sigma_safe(fit, fallback = 1e-6, ml = TRUE) else NA_real_
  structure(
    list(type = "glm", fit = fit, formula = formula, family = family,
         rhs_names = rhs_names, fallback = mean(as.numeric(y), na.rm = TRUE),
         sigma = sigma, component = component),
    class = c("nuisance_model", "list")
  )
}

predict_nuisance_model <- function(fit, newdata, type = c("auto", "response", "numeric"),
                                   eps = 1e-8, fallback = NULL,
                                   allow_fallback = TRUE) {
  type <- match.arg(type)
  n <- nrow(newdata)
  if (!n) return(numeric(0))
  if (is.null(fallback)) fallback <- fit$fallback %||% NA_real_
  out <- rep(NA_real_, n)
  vars <- fit$rhs_names %||% character(0)
  vars <- intersect(vars, names(newdata))
  keep <- if (length(vars)) stats::complete.cases(newdata[, vars, drop = FALSE]) else rep(TRUE, n)

  if (any(keep)) {
    if (identical(fit$type, "superlearner")) {
      pred <- tryCatch(
        as.numeric(SuperLearner::predict.SuperLearner(fit$fit, newdata = newdata[keep, vars, drop = FALSE])$pred),
        error = function(e) rep(NA_real_, sum(keep))
      )
    } else {
      fam <- .nuisance_family_name(fit$family)
      ptype <- if (identical(fam, "binomial") || identical(type, "response")) "response" else "response"
      pred <- tryCatch(
        as.numeric(stats::predict(fit$fit, newdata = newdata[keep, , drop = FALSE], type = ptype)),
        error = function(e) rep(NA_real_, sum(keep))
      )
    }
    if (length(pred) == sum(keep)) out[keep] <- pred
  }

  bad <- !is.finite(out)
  if (any(bad)) {
    if (!isTRUE(allow_fallback)) {
      .stop("Prediction failed and fallback is disabled for ", fit$component %||% "nuisance", ".")
    }
    if (length(fallback) == 1L) {
      out[bad] <- as.numeric(fallback)[1L]
    } else if (length(fallback) == n) {
      out[bad] <- as.numeric(fallback)[bad]
    } else {
      out[bad] <- as.numeric(fallback)[1L]
    }
  }

  fam <- .nuisance_family_name(fit$family)
  if (identical(fam, "binomial") || identical(type, "response")) {
    return(clamp01(out, eps = eps))
  }
  out
}

nuisance_sigma_safe <- function(fit, fallback = 1e-6) {
  if (is.list(fit) && !is.null(fit$sigma) && is.finite(fit$sigma) && fit$sigma > 0) return(fit$sigma)
  model_sigma_safe(fit, fallback = fallback, ml = TRUE)
}
