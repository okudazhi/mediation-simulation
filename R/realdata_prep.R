################################################################################
# realdata_prep.R
#
# Preparation helpers for the real-data layer.
#
# The authoritative applied-data input is a SAS data set placed directly under
# data/.
################################################################################

read_analysis_data <- function(path,
                               format = c("auto", "sas", "csv"),
                               catalog_path = NULL,
                               encoding = NULL,
                               na_strings = c(".", "NA", ""),
                               stringsAsFactors = FALSE) {
  format <- match.arg(format)
  if (identical(format, "auto")) {
    ext <- tolower(tools::file_ext(path))
    if (identical(ext, "sas7bdat")) {
      format <- "sas"
    } else if (identical(ext, "csv")) {
      format <- "csv"
    } else {
      .stop("The real-data layer expects a SAS data set (.sas7bdat) or an explicitly supplied CSV for dry runs. Got: ", path)
    }
  }

  if (!file.exists(path)) .stop("Input data file not found: ", path)

  if (identical(format, "sas")) {
    if (!requireNamespace("haven", quietly = TRUE)) {
      .stop("Package 'haven' is required to read SAS files. Install it via install.packages('haven').")
    }
    cat_file <- NULL
    if (!is.null(catalog_path) && nzchar(catalog_path) && file.exists(catalog_path)) {
      cat_file <- catalog_path
    }
    dat <- haven::read_sas(data_file = path, catalog_file = cat_file, encoding = encoding)
    dat <- as.data.frame(dat, stringsAsFactors = stringsAsFactors, check.names = FALSE)
  } else {
    dat <- utils::read.csv(path, stringsAsFactors = stringsAsFactors, check.names = FALSE,
                           na.strings = na_strings)
  }

  names(dat) <- tolower(names(dat))
  clean_dot_na(dat, na_strings = na_strings)
}

# Backward-compatible alias used by the original wrapper.
read_medidata <- read_analysis_data

.as_character_flag <- function(x) {
  if (is.factor(x)) return(trimws(as.character(x)))
  trimws(as.character(x))
}

.keep_yes_flag <- function(x) {
  x_chr <- toupper(.as_character_flag(x))
  x_chr %in% c("Y", "YES", "1", "TRUE", "T")
}

.normalize_baseline_map <- function(baseline) {
  if (is.null(baseline) || !length(baseline)) {
    .stop("var_spec$baseline must be a non-empty character vector or named character vector.")
  }
  baseline <- as.character(baseline)
  nm <- names(baseline)
  if (is.null(nm) || any(!nzchar(nm))) nm <- baseline
  names(baseline) <- nm
  baseline
}

.normalize_l_block_spec <- function(L) {
  if (is.null(L) || !length(L)) {
    .stop("var_spec$L must be a non-empty named list.")
  }
  if (is.atomic(L)) L <- list(L = as.character(L))
  if (is.null(names(L)) || any(!nzchar(names(L)))) {
    names(L) <- paste0("block", seq_along(L))
  }
  out <- lapply(L, as.character)
  lens <- vapply(out, length, integer(1))
  if (length(unique(lens)) != 1L) {
    .stop("All L blocks must have the same length.")
  }
  out
}

.resolve_A_cols <- function(var_spec, analysis_cfg, T) {
  if (!identical(analysis_cfg$analysis_type %||% "arm", "arm")) {
    .stop("Only randomized arm comparison is supported in the manuscript-concordant real-data layer.")
  }
  x <- rep(as.character(var_spec$Z), T)
  if (length(x) == 1L) x <- rep(x, T)
  if (length(x) != T) {
    .stop("Resolved treatment columns must have length 1 or T. Got length=", length(x), ", T=", T)
  }
  x
}

.resolve_lag_init <- function(L_lag_init, block_names) {
  if (is.null(L_lag_init)) {
    out <- as.list(rep(0, length(block_names)))
    names(out) <- block_names
    return(out)
  }
  if (!is.list(L_lag_init)) .stop("var_spec$L_lag_init must be a list or NULL.")
  out <- L_lag_init
  miss <- setdiff(block_names, names(out))
  for (nm in miss) out[[nm]] <- 0
  out[block_names]
}

.lag_init_column_names <- function(L_lag_init) {
  if (is.null(L_lag_init) || !length(L_lag_init)) return(character(0))

  keep <- vapply(
    L_lag_init,
    function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x),
    logical(1)
  )

  vals <- as.character(unlist(L_lag_init[keep], use.names = FALSE))
  unique(vals[nzchar(vals) & !is.na(vals)])
}

.copy_lag_init_sources_to_wide <- function(wide, dat, spec) {
  lag_init_cols <- .lag_init_column_names(spec$L_lag_init)
  if (!length(lag_init_cols)) return(wide)

  missing <- setdiff(lag_init_cols, names(dat))
  if (length(missing)) {
    .stop("L_lag_init source columns missing from cleaned raw data: ",
          paste(missing, collapse = ", "))
  }

  for (src in lag_init_cols) {
    if (!src %in% names(wide)) {
      wide[[src]] <- as.numeric(dat[[src]])
    }
  }

  wide
}

.derive_retention_matrix <- function(dat, censor_cols, T) {
  if (is.null(censor_cols) || !length(censor_cols)) {
    r_visit <- matrix(1, nrow = nrow(dat), ncol = T)
    colnames(r_visit) <- paste0("R_visit", 0:(T - 1L))
    return(list(R_visit = r_visit, R_final = rep(1, nrow(dat))))
  }
  assert_cols(dat, censor_cols, "raw data (censor columns)")
  ev <- as.matrix(dat[, censor_cols, drop = FALSE])
  storage.mode(ev) <- "numeric"
  ev[is.na(ev)] <- 0
  ev <- ifelse(ev == 1, 1, 0)

  if (ncol(ev) == (T + 1L)) {
    stay <- 1 - ev
    r_full <- t(apply(stay, 1, cumprod))
    r_visit <- r_full[, 1:T, drop = FALSE]
    r_final <- r_full[, T + 1L]
  } else if (ncol(ev) == T) {
    stay <- 1 - ev
    r_visit <- t(apply(stay, 1, cumprod))
    r_final <- r_visit[, T]
  } else {
    .stop("censor column count must equal T or T+1. Got ", ncol(ev), " with T=", T)
  }

  colnames(r_visit) <- paste0("R_visit", 0:(T - 1L))
  list(R_visit = r_visit, R_final = as.numeric(r_final))
}

.get_lag0_value <- function(wide, spec, block_name) {
  src <- spec$L_lag_init[[block_name]] %||% 0
  if (length(src) == 1L && is.character(src) && src %in% names(wide)) {
    return(as.numeric(wide[[src]]))
  }
  if (length(src) == 1L && is.numeric(src)) {
    return(rep(as.numeric(src), nrow(wide)))
  }
  if (length(src) == nrow(wide)) {
    return(as.numeric(src))
  }
  rep(0, nrow(wide))
}

build_realdata_spec <- function(var_spec,
                                analysis_cfg,
                                T = NULL,
                                analysis_name = NULL) {
  baseline_map <- .normalize_baseline_map(var_spec$baseline)
  T <- as.integer(T %||% length(var_spec$M1))[1L]
  if (!is.finite(T) || T < 1L) .stop("T must be a positive integer.")
  L_cols <- .normalize_l_block_spec(var_spec$L)
  lens <- vapply(L_cols, length, integer(1))
  if (any(lens != T)) {
    .stop("Each L block in var_spec$L must have length T.")
  }

  spec <- list(
    T = T,
    analysis_name = analysis_name %||% (analysis_cfg$name %||% "analysis"),
    analysis_type = as.character(analysis_cfg$analysis_type %||% "arm"),
    id_col = as.character(var_spec$id),
    include_flag = as.character(var_spec$include_flag %||% ""),
    randomization_flag = as.character(var_spec$randomization_flag %||% ""),
    baseline_map = baseline_map,
    baseline_vars = names(baseline_map),
    Z_col = as.character(var_spec$Z %||% ""),
    A_cols = .resolve_A_cols(var_spec, analysis_cfg, T),
    M1_cols = as.character(var_spec$M1),
    M2_cols = as.character(var_spec$M2),
    L_cols = L_cols,
    L_names = names(L_cols),
    L_lag_init = .resolve_lag_init(var_spec$L_lag_init %||% NULL, names(L_cols)),
    Y_final_col = as.character(var_spec$Y_final),
    censor_cols = as.character(var_spec$censor %||% character(0)),
    reg_a = normalize_regimen(analysis_cfg$reg_a, T, "analysis_cfg$reg_a"),
    reg_as = normalize_regimen(analysis_cfg$reg_as, T, "analysis_cfg$reg_as"),
    treat_mech = as.character(analysis_cfg$treat_mech %||% "estimated_regimen"),
    p_rct = analysis_cfg$p_rct %||% NULL,
    subset_expr = analysis_cfg$subset_expr %||% NULL,
    trunc_prob = analysis_cfg$trunc_prob %||% c(0.01, 0.99),
    cuts = as.integer(analysis_cfg$cuts %||% (0:(T - 1L))),
    use_joint_supplement = isTRUE(analysis_cfg$use_joint_supplement %||% TRUE),
    description = analysis_cfg$description %||% NULL
  )
  class(spec) <- c("realdata_spec", "list")
  spec
}

.required_columns_from_spec <- function(spec) {
  unique(c(
    spec$id_col,
    unname(spec$baseline_map),
    spec$Z_col,
    spec$A_cols,
    spec$M1_cols,
    spec$M2_cols,
    unlist(spec$L_cols, use.names = FALSE),
    spec$Y_final_col,
    spec$censor_cols,
    spec$include_flag,
    spec$randomization_flag,
    .lag_init_column_names(spec$L_lag_init)
  ))
}

build_realdata_long_from_wide <- function(wide, spec) {
  T <- spec$T
  rows <- vector("list", T)

  for (tt in 0:(T - 1L)) {
    idx <- tt + 1L
    row <- data.frame(
      ID = as.character(wide$ID),
      t = as.integer(tt),
      visit = factor(tt, levels = 0:(T - 1L)),
      A = as.numeric(wide[[paste0("A", tt)]]),
      M1 = as.numeric(wide[[paste0("M1_", tt)]]),
      M2 = as.numeric(wide[[paste0("M2_", tt)]]),
      Y_final = as.numeric(wide$Y_final),
      stringsAsFactors = FALSE
    )
    row$A_M1 <- row$A * row$M1

    row$A_lag <- if (tt == 0L) 0 else as.numeric(wide[[paste0("A", tt - 1L)]])
    row$cumA_prev <- if (tt == 0L) 0 else rowSums(wide[, paste0("A", 0:(tt - 1L)), drop = FALSE], na.rm = TRUE)
    row$cumA_curr <- rowSums(wide[, paste0("A", 0:tt), drop = FALSE], na.rm = TRUE)
    row$M1_lag <- if (tt == 0L) 0 else as.numeric(wide[[paste0("M1_", tt - 1L)]])
    row$M2_lag <- if (tt == 0L) 0 else as.numeric(wide[[paste0("M2_", tt - 1L)]])
    row$R_current <- as.numeric(wide[[paste0("R_visit", tt)]])
    row$R_next <- if (tt < (T - 1L)) as.numeric(wide[[paste0("R_visit", tt + 1L)]]) else as.numeric(wide$R_final)
    row$uncensored_final <- as.numeric(wide$R_final)

    if ("Z" %in% names(wide)) row$Z <- as.numeric(wide$Z)
    for (nm in spec$baseline_vars) row[[nm]] <- wide[[nm]]

    for (lnm in spec$L_names) {
      row[[paste0("L_", lnm)]] <- as.numeric(wide[[paste0("L_", lnm, "_", tt)]])
      if (tt == 0L) {
        row[[paste0("L_", lnm, "_lag")]] <- as.numeric(.get_lag0_value(wide, spec, lnm))
      } else {
        row[[paste0("L_", lnm, "_lag")]] <- as.numeric(wide[[paste0("L_", lnm, "_", tt - 1L)]])
      }
    }

    rows[[idx]] <- row
  }

  long <- do.call(rbind, rows)
  rownames(long) <- NULL
  long
}

prepared_from_wide_and_spec <- function(wide, spec) {
  long <- build_realdata_long_from_wide(wide, spec)
  baseline <- wide[, c("ID", spec$baseline_vars), drop = FALSE]
  if ("Z" %in% names(wide)) baseline$Z <- wide$Z
  rownames(baseline) <- NULL

  out <- list(
    spec = spec,
    analysis_name = spec$analysis_name,
    wide = wide,
    long = long,
    baseline = baseline,
    raw_clean = NULL,
    meta = list(
      T = spec$T,
      analysis_name = spec$analysis_name,
      n_subjects = nrow(wide)
    )
  )
  class(out) <- c("prepared_realdata", "list")
  out
}

prepare_realdata_data <- function(dat_raw,
                                  var_spec = default_esax_var_spec(),
                                  analysis_cfg = default_esax_arm_analysis_cfg(),
                                  T = NULL,
                                  analysis_name = NULL) {
  spec <- build_realdata_spec(
    var_spec = var_spec,
    analysis_cfg = analysis_cfg,
    T = T %||% length(var_spec$M1),
    analysis_name = analysis_name
  )

  dat <- clean_dot_na(as.data.frame(dat_raw, stringsAsFactors = FALSE))
  if (nzchar(spec$include_flag) && spec$include_flag %in% names(dat)) {
    keep <- .keep_yes_flag(dat[[spec$include_flag]])
    keep[is.na(keep)] <- FALSE
    dat <- dat[keep, , drop = FALSE]
  }

  if (nzchar(spec$Z_col) && spec$Z_col %in% names(dat)) {
    dat$treatn <- dat[[spec$Z_col]]
  }

  if (!is.null(spec$subset_expr) && nzchar(spec$subset_expr)) {
    keep <- tryCatch(with(dat, eval(parse(text = spec$subset_expr))), error = function(e) e)
    if (inherits(keep, "error")) {
      .stop("Failed to evaluate subset_expr: ", keep$message)
    }
    if (!is.logical(keep) || length(keep) != nrow(dat)) {
      .stop("subset_expr must evaluate to a logical vector of length nrow(data).")
    }
    keep[is.na(keep)] <- FALSE
    dat <- dat[keep, , drop = FALSE]
  }

  assert_cols(dat, .required_columns_from_spec(spec), "raw data")

  numeric_cols <- unique(c(
    spec$Z_col, spec$A_cols, spec$M1_cols, spec$M2_cols,
    unlist(spec$L_cols, use.names = FALSE), spec$Y_final_col, spec$censor_cols,
    .lag_init_column_names(spec$L_lag_init)
  ))
  numeric_cols <- intersect(numeric_cols, names(dat))
  dat <- coerce_numeric_cols(dat, numeric_cols)

  wide <- data.frame(
    ID = as.character(dat[[spec$id_col]]),
    Y_final = as.numeric(dat[[spec$Y_final_col]]),
    stringsAsFactors = FALSE
  )
  for (nm in spec$baseline_vars) {
    wide[[nm]] <- dat[[spec$baseline_map[[nm]]]]
  }
  if (nzchar(spec$Z_col) && spec$Z_col %in% names(dat)) wide$Z <- as.numeric(dat[[spec$Z_col]])

  # Keep character-valued L_lag_init source columns in prepared$wide.
  # These may be referenced by build_realdata_long_from_wide() even when
  # they are not baseline variables.
  wide <- .copy_lag_init_sources_to_wide(wide = wide, dat = dat, spec = spec)

  for (tt in 0:(spec$T - 1L)) {
    idx <- tt + 1L
    wide[[paste0("A", tt)]] <- as.numeric(dat[[spec$A_cols[idx]]])
    wide[[paste0("M1_", tt)]] <- as.numeric(dat[[spec$M1_cols[idx]]])
    wide[[paste0("M2_", tt)]] <- as.numeric(dat[[spec$M2_cols[idx]]])
    for (lnm in spec$L_names) {
      wide[[paste0("L_", lnm, "_", tt)]] <- as.numeric(dat[[spec$L_cols[[lnm]][idx]]])
    }
  }

  ret <- .derive_retention_matrix(dat, spec$censor_cols, T = spec$T)
  for (tt in 0:(spec$T - 1L)) {
    wide[[paste0("R_visit", tt)]] <- as.numeric(ret$R_visit[, tt + 1L])
  }
  wide$R_final <- as.numeric(ret$R_final)
  wide$analysis_name <- spec$analysis_name
  rownames(wide) <- NULL

  out <- prepared_from_wide_and_spec(wide, spec)
  out$raw_clean <- dat
  out
}

# Backward-compatible alias used by the original runner.
prepare_realdata_analysis <- prepare_realdata_data

bootstrap_prepared_realdata <- function(prepared, seed = NULL) {
  wide_b <- subject_bootstrap_sample(prepared$wide, id_col = "ID", seed = seed)
  prepared_from_wide_and_spec(wide_b, prepared$spec)
}
