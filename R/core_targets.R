################################################################################
# core_targets.R
#
# Shared target registries and estimand mappings.
#
# The fixed-horizon mapping preserves compatibility with the simulation layer.
# The time-resolved registry/effect builders support the real-data manuscript.
################################################################################

component_mean_keys <- function() {
  c(
    "mu_nat_a", "mu_nat_as",
    "mu_joint_aa", "mu_joint_asas", "mu_joint_aas",
    "mu_sep_aaa", "mu_sep_asas_asas", "mu_sep_a_asas", "mu_sep_a_aas"
  )
}

validate_component_means <- function(means) {
  req <- component_mean_keys()
  miss <- setdiff(req, names(means))
  if (length(miss) > 0L) {
    .stop("means missing keys: ", paste(miss, collapse = ", "))
  }
  invisible(TRUE)
}

primary_effect_keys <- function() {
  data.frame(
    estimand = c(rep("yamamuro", 5L), rep("tai", 6L)),
    effect = c(
      "TE_Y", "IDE_Y", "IIE_Y1", "IIE_Y2", "R_Y",
      "ITE_J", "IDE_J", "PSE1", "PSE2", "MI", "xi"
    ),
    stringsAsFactors = FALSE
  )
}

.compute_safe_div <- function(num, den, eps = 1e-12) {
  if (!is.finite(den) || abs(den) < eps) return(NA_real_)
  num / den
}

compute_estimands_from_means <- function(means, include_descriptive_ratios = FALSE) {
  validate_component_means(means)
  mu <- means

  TE_Y   <- mu$mu_nat_a - mu$mu_nat_as
  IDE_Y  <- mu$mu_sep_a_asas - mu$mu_sep_asas_asas
  IIE_Y1 <- mu$mu_sep_a_aas  - mu$mu_sep_a_asas
  IIE_Y2 <- mu$mu_sep_aaa    - mu$mu_sep_a_aas
  R_Y    <- TE_Y - IDE_Y - IIE_Y1 - IIE_Y2

  ITE_J  <- mu$mu_joint_aa - mu$mu_joint_asas
  IDE_J  <- mu$mu_joint_aas - mu$mu_joint_asas
  PSE1   <- IIE_Y1
  PSE2   <- IIE_Y2
  MI     <- ITE_J - IDE_J - PSE1 - PSE2
  xi     <- TE_Y - ITE_J

  out <- list(
    yamamuro = list(
      TE_Y = TE_Y,
      IDE_Y = IDE_Y,
      IIE_Y1 = IIE_Y1,
      IIE_Y2 = IIE_Y2,
      R_Y = R_Y
    ),
    tai = list(
      ITE_J = ITE_J,
      IDE_J = IDE_J,
      PSE1 = PSE1,
      PSE2 = PSE2,
      MI = MI,
      xi = xi
    )
  )

  if (isTRUE(include_descriptive_ratios)) {
    IIE_Y_total <- IIE_Y1 + IIE_Y2
    tai_indirect_total <- PSE1 + PSE2 + MI
    out$descriptive_ratios <- list(
      PE_Y = .compute_safe_div(IIE_Y_total, TE_Y),
      PE_J = .compute_safe_div(tai_indirect_total, ITE_J),
      PE_J_over_TE_Y = .compute_safe_div(tai_indirect_total, TE_Y)
    )
  }

  out
}

filter_primary_nonratio_effects <- function(effects_long) {
  req <- primary_effect_keys()
  merge(effects_long, req, by = c("estimand", "effect"), all = FALSE, sort = FALSE)
}

flatten_effects <- function(effects) {
  out <- list()
  for (est in names(effects)) {
    comp <- effects[[est]]
    out[[est]] <- data.frame(
      estimand = est,
      effect = names(comp),
      estimate = as.numeric(unlist(comp)),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

make_fixed_horizon_registry <- function(reg_a, reg_as, T = length(reg_a)) {
  reg_a  <- normalize_regimen(reg_a, T, "reg_a")
  reg_as <- normalize_regimen(reg_as, T, "reg_as")
  list(
    mu_nat_a = list(type = "natural", outer_regimen = reg_a),
    mu_nat_as = list(type = "natural", outer_regimen = reg_as),
    mu_joint_aa = list(type = "joint", outer_regimen = reg_a, med_regimen = reg_a),
    mu_joint_asas = list(type = "joint", outer_regimen = reg_as, med_regimen = reg_as),
    mu_joint_aas = list(type = "joint", outer_regimen = reg_a, med_regimen = reg_as),
    mu_sep_aaa = list(type = "separate", outer_regimen = reg_a,
                      med_regimen1 = reg_a, med_regimen2 = reg_a),
    mu_sep_asas_asas = list(type = "separate", outer_regimen = reg_as,
                            med_regimen1 = reg_as, med_regimen2 = reg_as),
    mu_sep_a_asas = list(type = "separate", outer_regimen = reg_a,
                         med_regimen1 = reg_as, med_regimen2 = reg_as),
    mu_sep_a_aas = list(type = "separate", outer_regimen = reg_a,
                        med_regimen1 = reg_a, med_regimen2 = reg_as)
  )
}

make_joint_draw_registry <- function(reg_a, reg_as, T = length(reg_a)) {
  # The supplementary Tai-style summaries require the same nine fixed-horizon
  # means because PSE^1/PSE^2 are constructed from the separate-draw bridge
  # means while IDE_J/ITE_J/MI use the joint-draw means.
  make_fixed_horizon_registry(reg_a = reg_a, reg_as = reg_as, T = T)
}

make_time_resolved_sep_registry <- function(reg_a, reg_as, cuts = NULL) {
  T <- length(reg_a)
  reg_a  <- normalize_regimen(reg_a, T, "reg_a")
  reg_as <- normalize_regimen(reg_as, T, "reg_as")
  if (is.null(cuts)) cuts <- 0:(T - 1L)
  cuts <- sort(unique(as.integer(cuts)))
  if (any(cuts < 0L | cuts > (T - 1L))) {
    .stop("cuts must lie in {0, ..., T-1}.")
  }

  out <- list(
    phi_as_as_as = list(type = "product", outer_regimen = reg_as,
                        med_regimen1 = reg_as, med_regimen2 = reg_as),
    phi_a_as_as = list(type = "product", outer_regimen = reg_a,
                       med_regimen1 = reg_as, med_regimen2 = reg_as),
    phi_a_a_as = list(type = "product", outer_regimen = reg_a,
                      med_regimen1 = reg_a, med_regimen2 = reg_as),
    phi_a_a_a = list(type = "product", outer_regimen = reg_a,
                     med_regimen1 = reg_a, med_regimen2 = reg_a)
  )

  for (tt in cuts) {
    reg_sw <- make_switching_regimen(reg_a, reg_as, tt)
    out[[sprintf("phi_switch_outer_t%d", tt)]] <- list(
      type = "product",
      outer_regimen = reg_sw,
      med_regimen1 = reg_as,
      med_regimen2 = reg_as
    )
    out[[sprintf("phi_switch_m1_t%d", tt)]] <- list(
      type = "product",
      outer_regimen = reg_a,
      med_regimen1 = reg_sw,
      med_regimen2 = reg_as
    )
    out[[sprintf("phi_switch_m2_t%d", tt)]] <- list(
      type = "product",
      outer_regimen = reg_a,
      med_regimen1 = reg_a,
      med_regimen2 = reg_sw
    )
  }
  out
}

compute_fixed_horizon_main_effects_from_means <- function(means) {
  validate_component_means(means)
  te <- means$mu_nat_a - means$mu_nat_as
  ide <- means$mu_sep_a_asas - means$mu_sep_asas_asas
  iie1 <- means$mu_sep_a_aas - means$mu_sep_a_asas
  iie2 <- means$mu_sep_aaa - means$mu_sep_a_aas
  r_y <- te - ide - iie1 - iie2

  list(
    TE_Y = te,
    IDE_Y = ide,
    IIE_Y1 = iie1,
    IIE_Y2 = iie2,
    R_Y = r_y
  )
}

compute_joint_draw_effects_from_means <- function(means) {
  validate_component_means(means)
  ite_j <- means$mu_joint_aa - means$mu_joint_asas
  ide_j <- means$mu_joint_aas - means$mu_joint_asas
  pse1 <- means$mu_sep_a_aas - means$mu_sep_a_asas
  pse2 <- means$mu_sep_aaa - means$mu_sep_a_aas
  mi <- ite_j - ide_j - pse1 - pse2
  te <- means$mu_nat_a - means$mu_nat_as
  xi <- te - ite_j
  list(
    ITE_J = ite_j,
    IDE_J = ide_j,
    PSE1 = pse1,
    PSE2 = pse2,
    MI = mi,
    xi = xi
  )
}

compute_time_resolved_sep_effects_from_means <- function(means) {
  req_fixed <- c("phi_as_as_as", "phi_a_as_as", "phi_a_a_as", "phi_a_a_a")
  miss <- setdiff(req_fixed, names(means))
  if (length(miss) > 0L) {
    .stop("means missing required product-law keys: ", paste(miss, collapse = ", "))
  }

  cuts_outer <- grep("^phi_switch_outer_t[0-9]+$", names(means), value = TRUE)
  cuts_m1 <- grep("^phi_switch_m1_t[0-9]+$", names(means), value = TRUE)
  cuts_m2 <- grep("^phi_switch_m2_t[0-9]+$", names(means), value = TRUE)
  cuts <- sort(unique(as.integer(gsub("^.*_t", "", c(cuts_outer, cuts_m1, cuts_m2)))))

  fixed <- list(
    IDE_Y = means$phi_a_as_as - means$phi_as_as_as,
    IIE_Y1 = means$phi_a_a_as - means$phi_a_as_as,
    IIE_Y2 = means$phi_a_a_a - means$phi_a_a_as,
    IIE_Y_total = (means$phi_a_a_as - means$phi_a_as_as) +
                  (means$phi_a_a_a - means$phi_a_a_as)
  )

  rows <- vector("list", length(cuts))
  for (ii in seq_along(cuts)) {
    tt <- cuts[[ii]]
    key_out <- sprintf("phi_switch_outer_t%d", tt)
    key_m1  <- sprintf("phi_switch_m1_t%d", tt)
    key_m2  <- sprintf("phi_switch_m2_t%d", tt)
    if (!all(c(key_out, key_m1, key_m2) %in% names(means))) next
    rows[[ii]] <- data.frame(
      cut_t = tt,
      IDE_Y_le = means[[key_out]] - means$phi_as_as_as,
      IDE_Y_gt = means$phi_a_as_as - means[[key_out]],
      IIE_Y1_le = means[[key_m1]] - means$phi_a_as_as,
      IIE_Y1_gt = means$phi_a_a_as - means[[key_m1]],
      IIE_Y2_le = means[[key_m2]] - means$phi_a_a_as,
      IIE_Y2_gt = means$phi_a_a_a - means[[key_m2]],
      stringsAsFactors = FALSE
    )
  }

  time_resolved <- if (length(rows)) do.call(rbind, Filter(Negate(is.null), rows)) else data.frame()
  list(fixed = fixed, time_resolved = time_resolved)
}
