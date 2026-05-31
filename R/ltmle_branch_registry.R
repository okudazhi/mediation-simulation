################################################################################
# ltmle_branch_registry.R
#
# Diagnostic registries for the full exact-EIF ltmle_exact implementation.
################################################################################

ltmle_exact_component_registry <- function(reg_a, reg_as, T) {
  reg_a <- normalize_regimen(reg_a, T, "reg_a")
  reg_as <- normalize_regimen(reg_as, T, "reg_as")
  data.frame(
    component_id = seq_along(component_mean_keys()),
    component = component_mean_keys(),
    world_type = c("natural", "natural", "joint", "joint", "joint",
                   "separate", "separate", "separate", "separate"),
    outcome_regimen = c("a", "as", "a", "as", "a", "a", "as", "a", "a"),
    first_mediator_regimen = c("outcome", "outcome", "a", "as", "as",
                               "a", "as", "as", "a"),
    second_mediator_regimen = c("outcome", "outcome", "a", "as", "as",
                                "a", "as", "as", "as"),
    stringsAsFactors = FALSE
  )
}

.make_ltmle_exact_task_id <- function(component, t, node, process_type, role = NA_character_) {
  role2 <- if (is.na(role) || !nzchar(role)) "main" else gsub("[^A-Za-z0-9]+", "_", role)
  paste(component, t, node, process_type, role2, sep = "::")
}

.ltmle_exact_virtual_mixed_task_id <- function(component, previous_y_t, next_L_t) {
  paste0(
    component,
    "::",
    as.integer(previous_y_t),
    "::virtual_mixed_continuation_after_Y::to_outer_L_",
    as.integer(next_L_t)
  )
}

ltmle_exact_factor_task_registry <- function(reg_a, reg_as, T, node_spec = NULL) {
  reg <- ltmle_exact_component_registry(reg_a, reg_as, T)

  L_nodes <- as.character(node_spec$L_order %||% character(0))
  if (!length(L_nodes)) {
    L_blocks <- node_spec$L_blocks %||% NULL
    if (!is.null(L_blocks)) L_nodes <- names(L_blocks)
  }
  if (!length(L_nodes)) L_nodes <- "L"

  outcome_type <- node_spec$outcome_type %||% "longitudinal"
  rows <- list()
  ii <- 0L

  y_task_exists <- function(tt) {
    identical(outcome_type, "longitudinal") || tt == T
  }
  task_id_for <- function(rr, tt, node, process_type, role = NA_character_) {
    .make_ltmle_exact_task_id(reg$component[rr], tt, node, process_type, role)
  }
  y_id <- function(rr, tt) {
    task_id_for(rr, tt, "Y", "outcome_process")
  }
  L_id <- function(rr, tt, L_node, process_type = "post_mediator_covariate_transition") {
    task_id_for(rr, tt, L_node, process_type)
  }
  last_L_downstream_id <- function(rr, tt) {
    if (y_task_exists(tt)) return(y_id(rr, tt))
    if (tt < T) return(time_terminal_id(rr, tt + 1L))
    "observed_terminal_outcome"
  }
  first_L_terminal_id <- function(rr, tt) {
    L_id(rr, tt, L_nodes[1L])
  }
  required_mediator_history_state <- function(rr, node) {
    wt <- as.character(reg$world_type[rr])
    if (identical(wt, "joint")) return("joint_law")
    if (identical(wt, "separate")) {
      if (identical(node, "M1")) return("first_law")
      if (identical(node, "M2")) return("second_law")
    }
    "outcome"
  }
  normalized_mediator_regimen <- function(rr, mediator_regimen) {
    lab <- as.character(mediator_regimen)
    if (identical(lab, "outcome")) as.character(reg$outcome_regimen[rr]) else lab
  }
  required_mixed_continuation_spec <- function(rr) {
    outer_regimen <- as.character(reg$outcome_regimen[rr])
    m1_regimen <- normalized_mediator_regimen(rr, reg$first_mediator_regimen[rr])
    m2_regimen <- normalized_mediator_regimen(rr, reg$second_mediator_regimen[rr])
    m1_state <- required_mediator_history_state(rr, "M1")
    m2_state <- required_mediator_history_state(rr, "M2")
    cross_regimen <- !identical(outer_regimen, m1_regimen) ||
      !identical(outer_regimen, m2_regimen)
    split_mediator_history <- !identical(m1_state, m2_state)
    mixed_target_required <- isTRUE(cross_regimen) || isTRUE(split_mediator_history)
    list(
      outcome_history_state = "outcome",
      m1_history_state = m1_state,
      m2_history_state = m2_state,
      auxiliary_mediator_history_state = "not_applicable",
      outer_regimen = outer_regimen,
      m1_regimen = m1_regimen,
      m2_regimen = m2_regimen,
      continuation_target_label = if (isTRUE(mixed_target_required)) {
        "mixed_outer_LY_history_with_mediator_law_path"
      } else {
        "existing_single_state_continuation"
      },
      existing_task_can_represent_target = !isTRUE(mixed_target_required)
    )
  }
  requires_virtual_mixed_continuation <- function(rr, tt) {
    if (tt >= T) return(FALSE)
    spec <- required_mixed_continuation_spec(rr)
    identical(spec$continuation_target_label, "mixed_outer_LY_history_with_mediator_law_path") &&
      !isTRUE(spec$existing_task_can_represent_target)
  }
  separate_mediator_regimens_differ <- function(rr) {
    identical(as.character(reg$world_type[rr]), "separate") &&
      !identical(
        as.character(reg$first_mediator_regimen[rr]),
        as.character(reg$second_mediator_regimen[rr])
      )
  }
  virtual_mixed_continuation_id <- function(rr, previous_y_t) {
    .ltmle_exact_virtual_mixed_task_id(reg$component[rr], previous_y_t, previous_y_t + 1L)
  }
  law_L_terminal_id <- function(rr, tt, process_type) {
    L_id(rr, tt, L_nodes[1L], process_type)
  }
  law_history_entry_id <- function(rr, tt, process_type) {
    if (tt < 1L || tt >= T) {
      .stop("Invalid law-history entry time for ltmle_exact: ", tt)
    }
    law_L_terminal_id(rr, tt, process_type)
  }
  joint_law_history_entry_id <- function(rr, tt) {
    law_history_entry_id(rr, tt, "joint_stochastic_mediator_intervention_law")
  }
  first_law_history_entry_id <- function(rr, tt) {
    law_history_entry_id(rr, tt, "first_mediator_stochastic_intervention_law")
  }
  second_law_history_entry_id <- function(rr, tt) {
    law_history_entry_id(rr, tt, "second_mediator_stochastic_intervention_law")
  }
  law_Y_history_entry_id <- function(rr, tt, process_type) {
    if (tt < 1L || tt >= T) {
      .stop("Invalid law-Y history entry time for ltmle_exact: ", tt)
    }
    task_id_for(rr, tt, "Y", process_type, role = law_Y_role)
  }
  virtual_mixed_downstream_source_id <- function(rr, previous_y_t) {
    wt <- as.character(reg$world_type[rr])
    if (identical(wt, "joint")) {
      return(law_Y_history_entry_id(
        rr,
        previous_y_t,
        "joint_stochastic_mediator_intervention_law"
      ))
    }
    if (identical(wt, "separate")) {
      if (!isTRUE(separate_mediator_regimens_differ(rr))) {
        return(law_Y_history_entry_id(
          rr,
          previous_y_t,
          "second_mediator_stochastic_intervention_law"
        ))
      }
      return(law_Y_history_entry_id(
        rr,
        previous_y_t,
        "first_mediator_stochastic_intervention_law"
      ))
    }
    .stop("Virtual mixed continuation requested for unsupported world_type: ", wt)
  }
  virtual_mixed_m1_task_id <- function(rr, tt) {
    wt <- as.character(reg$world_type[rr])
    if (identical(wt, "joint")) {
      return(task_id_for(rr, tt, "M1", "joint_stochastic_mediator_intervention_law"))
    }
    if (identical(wt, "separate")) {
      return(task_id_for(rr, tt, "M1", "first_mediator_stochastic_intervention_law"))
    }
    NA_character_
  }
  virtual_mixed_m2_task_id <- function(rr, tt) {
    wt <- as.character(reg$world_type[rr])
    if (identical(wt, "joint")) {
      return(task_id_for(rr, tt, "M2", "joint_stochastic_mediator_intervention_law"))
    }
    if (identical(wt, "separate")) {
      return(task_id_for(rr, tt, "M2", "second_mediator_stochastic_intervention_law"))
    }
    NA_character_
  }
  outcome_Y_source_after_t <- function(rr, tt) {
    if (tt == T) return("observed_terminal_outcome")

    if (requires_virtual_mixed_continuation(rr, tt)) {
      return(virtual_mixed_continuation_id(rr, tt))
    }

    wt <- as.character(reg$world_type[rr])
    if (identical(wt, "joint")) {
      return(joint_law_history_entry_id(rr, tt))
    }
    if (identical(wt, "separate")) {
      return(first_law_history_entry_id(rr, tt))
    }

    task_id_for(rr, tt + 1L, "M1", "observed_mediator_process")
  }
  outcome_L_after_mediator_law_id <- function(rr, tt) {
    first_L_terminal_id(rr, tt)
  }
  second_law_entry_after_first_law_mediator <- function(rr, tt) {
    if (tt <= 1L) .stop("first-law mediator cannot hand off to second-law history at t <= 1.")
    second_law_history_entry_id(rr, tt - 1L)
  }
  second_M1_role <- "preceding first mediator in the second-mediator stochastic intervention law"
  first_M2_role <- "auxiliary second mediator in the first-mediator stochastic intervention law"
  law_Y_role <- "outcome history for stochastic mediator intervention law"
  mediator_role_for <- function(node, process_type, role = NA_character_) {
    if (identical(process_type, "joint_stochastic_mediator_intervention_law") &&
        identical(node, "M1")) {
      return("joint_target_M1")
    }
    if (identical(process_type, "joint_stochastic_mediator_intervention_law") &&
        identical(node, "M2")) {
      return("joint_target_M2")
    }
    if (identical(process_type, "first_mediator_stochastic_intervention_law") &&
        identical(node, "M1")) {
      return("first_target_M1")
    }
    if (identical(process_type, "first_mediator_stochastic_intervention_law") &&
        identical(node, "M2")) {
      return("first_auxiliary_M2")
    }
    if (identical(process_type, "second_mediator_stochastic_intervention_law") &&
        identical(node, "M1")) {
      return("second_auxiliary_M1")
    }
    if (identical(process_type, "second_mediator_stochastic_intervention_law") &&
        identical(node, "M2")) {
      return("second_target_M2")
    }
    NA_character_
  }
  time_terminal_id <- function(rr, tt) {
    if (tt < 1L || tt > T) .stop("Invalid ltmle_exact time-terminal lookup: ", tt)
    first_L_terminal_id(rr, tt)
  }
  law_next_mediator_id <- function(rr, tt, process_type) {
    if (identical(process_type, "second_mediator_stochastic_intervention_law")) {
      return(task_id_for(rr, tt + 1L, "M1", process_type, role = second_M1_role))
    }
    task_id_for(rr, tt + 1L, "M1", process_type)
  }
  law_downstream_after_current_mediators <- function(rr, tt, process_type) {
    outcome_L_after_mediator_law_id(rr, tt)
  }

  add_row <- function(rr, tt, node, process_type, assigned_label,
                      conditioning_history_type,
                      uses_current_mediator_density_ratio,
                      uses_current_censoring_factor,
                      role = NA_character_,
                      observed_source,
                      task_id_override = NA_character_,
                      virtual_mixed_task = FALSE,
                      virtual_mixed_previous_y_t = NA_integer_,
                      virtual_mixed_previous_y_task_id = NA_character_,
                      virtual_mixed_m1_task_id = NA_character_,
                      virtual_mixed_m2_task_id = NA_character_,
                      virtual_mixed_continuation_target_label = NA_character_) {
    if (missing(observed_source)) {
      .stop("Every ltmle_exact task must explicitly define an observed pseudo-outcome source.")
    }
    if (!is.character(observed_source) || length(observed_source) != 1L || !nzchar(observed_source)) {
      .stop("Invalid ltmle_exact pseudo-outcome source for component ", reg$component[rr],
            " node=", node, " t=", tt, ".")
    }
    id <- if (!is.na(task_id_override) && nzchar(as.character(task_id_override))) {
      as.character(task_id_override)
    } else {
      .make_ltmle_exact_task_id(reg$component[rr], tt, node, process_type, role)
    }
    ii <<- ii + 1L
    rows[[ii]] <<- data.frame(
      component_id = reg$component_id[rr],
      component = reg$component[rr],
      world_type = reg$world_type[rr],
      t = tt,
      node = node,
      process_type = process_type,
      task_id = id,
      observed_pseudooutcome_source_task_id = observed_source,
      source_is_explicit = TRUE,
      assigned_treatment_regimen_label = assigned_label,
      conditioning_history_type = conditioning_history_type,
      uses_current_mediator_density_ratio = isTRUE(uses_current_mediator_density_ratio),
      uses_current_censoring_factor = isTRUE(uses_current_censoring_factor),
      role = role,
      mediator_role = mediator_role_for(node, process_type, role),
      virtual_mixed_task = isTRUE(virtual_mixed_task),
      virtual_mixed_previous_y_t = as.integer(virtual_mixed_previous_y_t),
      virtual_mixed_previous_y_task_id = as.character(virtual_mixed_previous_y_task_id %||% NA_character_),
      virtual_mixed_m1_task_id = as.character(virtual_mixed_m1_task_id %||% NA_character_),
      virtual_mixed_m2_task_id = as.character(virtual_mixed_m2_task_id %||% NA_character_),
      virtual_mixed_continuation_target_label =
        as.character(virtual_mixed_continuation_target_label %||% NA_character_),
      stringsAsFactors = FALSE
    )
  }

  add_L_chain <- function(rr,
                          tt,
                          process_type,
                          assigned_label,
                          conditioning_history_type,
                          downstream_source) {
    for (ll in seq_along(L_nodes)) {
      L_node <- L_nodes[ll]
      observed_source <- if (ll < length(L_nodes)) {
        task_id_for(rr, tt, L_nodes[ll + 1L], process_type)
      } else {
        downstream_source
      }
      add_row(
        rr = rr,
        tt = tt,
        node = L_node,
        process_type = process_type,
        assigned_label = assigned_label,
        conditioning_history_type = conditioning_history_type,
        uses_current_mediator_density_ratio = TRUE,
        uses_current_censoring_factor = TRUE,
        observed_source = observed_source
      )
    }
  }
  add_law_history_chain <- function(rr,
                                    tt,
                                    process_type,
                                    assigned_label,
                                    conditioning_history_type) {
    if (tt >= T) return(invisible(NULL))
    next_mediator <- law_next_mediator_id(rr, tt, process_type)
    downstream_source <- if (identical(outcome_type, "longitudinal")) {
      task_id_for(rr, tt, "Y", process_type, role = law_Y_role)
    } else {
      next_mediator
    }
    add_L_chain(
      rr = rr,
      tt = tt,
      process_type = process_type,
      assigned_label = assigned_label,
      conditioning_history_type = conditioning_history_type,
      downstream_source = downstream_source
    )
    if (identical(outcome_type, "longitudinal")) {
      add_row(
        rr = rr,
        tt = tt,
        node = "Y",
        process_type = process_type,
        assigned_label = assigned_label,
        conditioning_history_type = conditioning_history_type,
        uses_current_mediator_density_ratio = TRUE,
        uses_current_censoring_factor = TRUE,
        role = law_Y_role,
        observed_source = next_mediator
      )
    }
    invisible(NULL)
  }

  for (rr in seq_len(nrow(reg))) {
    for (tt in seq_len(T)) {
      if (identical(outcome_type, "longitudinal") || tt == T) {
        observed_source <- outcome_Y_source_after_t(rr, tt)
        add_row(
          rr, tt, "Y", "outcome_process", "outcome_process_regimen",
          "outcome_process_history",
          uses_current_mediator_density_ratio = !identical(reg$world_type[rr], "natural"),
          uses_current_censoring_factor = TRUE,
          observed_source = observed_source
        )
      }

      if (requires_virtual_mixed_continuation(rr, tt)) {
        virtual_next_t <- tt + 1L
        virtual_spec <- required_mixed_continuation_spec(rr)
        add_row(
          rr = rr,
          tt = virtual_next_t,
          node = L_nodes[1L],
          process_type = "virtual_mixed_continuation_task",
          assigned_label = "outcome_process_regimen",
          conditioning_history_type = "virtual_mixed_continuation_history",
          uses_current_mediator_density_ratio = TRUE,
          uses_current_censoring_factor = TRUE,
          role = paste0("after_Y", tt, "_to_outer_L", virtual_next_t),
          observed_source = virtual_mixed_downstream_source_id(rr, tt),
          task_id_override = virtual_mixed_continuation_id(rr, tt),
          virtual_mixed_task = TRUE,
          virtual_mixed_previous_y_t = tt,
          virtual_mixed_previous_y_task_id = y_id(rr, tt),
          virtual_mixed_m1_task_id = virtual_mixed_m1_task_id(rr, virtual_next_t),
          virtual_mixed_m2_task_id = virtual_mixed_m2_task_id(rr, virtual_next_t),
          virtual_mixed_continuation_target_label =
            virtual_spec$continuation_target_label
        )
      }

      add_L_chain(
        rr = rr,
        tt = tt,
        process_type = "post_mediator_covariate_transition",
        assigned_label = "outcome_process_regimen",
        conditioning_history_type = "outcome_process_history",
        downstream_source = last_L_downstream_id(rr, tt)
      )

      if (identical(reg$world_type[rr], "joint")) {
        add_law_history_chain(
          rr = rr,
          tt = tt,
          process_type = "joint_stochastic_mediator_intervention_law",
          assigned_label = "joint_mediator_law_regimen",
          conditioning_history_type = "joint_stochastic_mediator_intervention_law_history"
        )
      }

      if (identical(reg$world_type[rr], "separate")) {
        add_law_history_chain(
          rr = rr,
          tt = tt,
          process_type = "first_mediator_stochastic_intervention_law",
          assigned_label = "first_mediator_law_regimen",
          conditioning_history_type = "first_mediator_stochastic_intervention_law_history"
        )
        add_law_history_chain(
          rr = rr,
          tt = tt,
          process_type = "second_mediator_stochastic_intervention_law",
          assigned_label = "second_mediator_law_regimen",
          conditioning_history_type = "second_mediator_stochastic_intervention_law_history"
        )
      }

      if (tt == 1L) next

      outcome_post_L_source <- first_L_terminal_id(rr, tt)

      if (identical(reg$world_type[rr], "joint")) {
        joint_M2_id <- task_id_for(rr, tt, "M2", "joint_stochastic_mediator_intervention_law")
        add_row(
          rr, tt, "M2", "joint_stochastic_mediator_intervention_law",
          "joint_mediator_law_regimen", "joint_stochastic_mediator_intervention_law_history",
          uses_current_mediator_density_ratio = TRUE,
          uses_current_censoring_factor = TRUE,
          observed_source = law_downstream_after_current_mediators(
            rr, tt, "joint_stochastic_mediator_intervention_law"
          )
        )
        add_row(
          rr, tt, "M1", "joint_stochastic_mediator_intervention_law",
          "joint_mediator_law_regimen", "joint_stochastic_mediator_intervention_law_history",
          uses_current_mediator_density_ratio = FALSE,
          uses_current_censoring_factor = TRUE,
          observed_source = joint_M2_id
        )
      } else if (identical(reg$world_type[rr], "separate")) {
        second_M2_id <- task_id_for(rr, tt, "M2", "second_mediator_stochastic_intervention_law")
        add_row(
          rr, tt, "M2", "second_mediator_stochastic_intervention_law",
          "second_mediator_law_regimen", "second_mediator_stochastic_intervention_law_history",
          uses_current_mediator_density_ratio = TRUE,
          uses_current_censoring_factor = TRUE,
          observed_source = law_downstream_after_current_mediators(
            rr, tt, "second_mediator_stochastic_intervention_law"
          )
        )
        add_row(
          rr, tt, "M1", "second_mediator_stochastic_intervention_law",
          "second_mediator_law_regimen", "second_mediator_stochastic_intervention_law_history",
          uses_current_mediator_density_ratio = FALSE,
          uses_current_censoring_factor = TRUE,
          role = second_M1_role,
          observed_source = second_M2_id
        )
        if (tt < T) {
          first_M2_id <- task_id_for(
            rr, tt, "M2", "first_mediator_stochastic_intervention_law",
            role = first_M2_role
          )
          add_row(
            rr, tt, "M2", "first_mediator_stochastic_intervention_law",
            "first_mediator_law_regimen", "first_mediator_stochastic_intervention_law_history",
            uses_current_mediator_density_ratio = TRUE,
            uses_current_censoring_factor = TRUE,
            role = first_M2_role,
            observed_source = second_law_entry_after_first_law_mediator(rr, tt)
          )
          first_M1_source <- first_M2_id
        } else {
          first_M1_source <- second_M2_id
        }
        add_row(
          rr, tt, "M1", "first_mediator_stochastic_intervention_law",
          "first_mediator_law_regimen", "first_mediator_stochastic_intervention_law_history",
          uses_current_mediator_density_ratio = FALSE,
          uses_current_censoring_factor = TRUE,
          observed_source = first_M1_source
        )
      } else {
        natural_M2_id <- task_id_for(rr, tt, "M2", "observed_mediator_process")
        add_row(
          rr, tt, "M2", "observed_mediator_process",
          "outcome_process_regimen", "observed_history",
          uses_current_mediator_density_ratio = FALSE,
          uses_current_censoring_factor = TRUE,
          observed_source = outcome_post_L_source
        )
        add_row(
          rr, tt, "M1", "observed_mediator_process",
          "outcome_process_regimen", "observed_history",
          uses_current_mediator_density_ratio = FALSE,
          uses_current_censoring_factor = TRUE,
          observed_source = natural_M2_id
        )
      }
    }
  }

  tasks <- do.call(rbind, rows)

  branch_key_for_process <- function(process_type) {
    if (process_type %in% c(
      "outcome_process",
      "post_mediator_covariate_transition",
      "observed_mediator_process"
    )) {
      return("outcome")
    }
    if (identical(process_type, "joint_stochastic_mediator_intervention_law")) {
      return("joint_law")
    }
    if (identical(process_type, "first_mediator_stochastic_intervention_law")) {
      return("first_law")
    }
    if (identical(process_type, "second_mediator_stochastic_intervention_law")) {
      return("second_law")
    }
    if (identical(process_type, "virtual_mixed_continuation_task")) {
      return("outcome")
    }
    "not_applicable"
  }
  target_m1_state_for_world <- function(world_type) {
    if (identical(world_type, "joint")) return("joint_law")
    if (identical(world_type, "separate")) return("first_law")
    "outcome"
  }
  target_m2_state_for_world <- function(world_type) {
    if (identical(world_type, "joint")) return("joint_law")
    if (identical(world_type, "separate")) return("second_law")
    "outcome"
  }
  task_history_states <- function(world_type, process_type, eval_state) {
    if (identical(world_type, "natural")) {
      aux <- if (eval_state %in% c("joint_law", "first_law", "second_law")) eval_state else "not_applicable"
      return(list(
        outcome = eval_state,
        m1 = eval_state,
        m2 = eval_state,
        auxiliary = aux
      ))
    }
    if (identical(world_type, "joint")) {
      if (identical(process_type, "joint_stochastic_mediator_intervention_law")) {
        return(list(outcome = "outcome", m1 = "joint_law", m2 = "joint_law", auxiliary = "joint_law"))
      }
      return(list(outcome = "outcome", m1 = "joint_law", m2 = "joint_law", auxiliary = "not_applicable"))
    }
    if (identical(world_type, "separate")) {
      if (identical(process_type, "first_mediator_stochastic_intervention_law")) {
        return(list(outcome = "outcome", m1 = "first_law", m2 = "first_law", auxiliary = "first_law"))
      }
      if (identical(process_type, "second_mediator_stochastic_intervention_law")) {
        return(list(outcome = "outcome", m1 = "second_law", m2 = "second_law", auxiliary = "second_law"))
      }
      return(list(outcome = "outcome", m1 = "first_law", m2 = "second_law", auxiliary = "not_applicable"))
    }
    list(outcome = eval_state, m1 = eval_state, m2 = eval_state, auxiliary = "not_applicable")
  }
  aux_state_for_source <- function(source_key) {
    if (source_key %in% c("joint_law", "first_law", "second_law")) source_key else "not_applicable"
  }

  task_branch_key <- vapply(tasks$process_type, branch_key_for_process, character(1))
  source_idx <- match(tasks$observed_pseudooutcome_source_task_id, tasks$task_id)
  source_branch_key <- ifelse(
    is.na(source_idx),
    ifelse(tasks$observed_pseudooutcome_source_task_id == "observed_terminal_outcome",
           "observed_terminal_outcome", "not_applicable"),
    task_branch_key[source_idx]
  )
  source_process_type <- ifelse(
    is.na(source_idx),
    ifelse(tasks$observed_pseudooutcome_source_task_id == "observed_terminal_outcome",
           "observed_terminal_outcome", "not_applicable"),
    tasks$process_type[source_idx]
  )
  source_node <- ifelse(
    is.na(source_idx),
    ifelse(tasks$observed_pseudooutcome_source_task_id == "observed_terminal_outcome",
           "observed_terminal_outcome", "not_applicable"),
    tasks$node[source_idx]
  )
  source_t <- ifelse(is.na(source_idx), NA_integer_, tasks$t[source_idx])
  reg_idx <- match(tasks$component, reg$component)
  task_world_type <- as.character(reg$world_type[reg_idx])
  outer_regimen <- as.character(reg$outcome_regimen[reg_idx])
  m1_regimen <- as.character(ifelse(
    reg$first_mediator_regimen[reg_idx] == "outcome",
    reg$outcome_regimen[reg_idx],
    reg$first_mediator_regimen[reg_idx]
  ))
  m2_regimen <- as.character(ifelse(
    reg$second_mediator_regimen[reg_idx] == "outcome",
    reg$outcome_regimen[reg_idx],
    reg$second_mediator_regimen[reg_idx]
  ))
  outcome_history_state <- character(nrow(tasks))
  m1_history_state <- character(nrow(tasks))
  m2_history_state <- character(nrow(tasks))
  aux_history_state <- character(nrow(tasks))
  boundary_direction <- character(nrow(tasks))
  boundary_type <- character(nrow(tasks))
  provenance_rule <- character(nrow(tasks))
  for (jj in seq_len(nrow(tasks))) {
    src_key <- source_branch_key[jj]
    consume_key <- task_branch_key[jj]
    wt <- task_world_type[jj]
    if (identical(src_key, "observed_terminal_outcome") ||
        identical(src_key, "not_applicable")) {
      outcome_history_state[jj] <- "not_applicable"
      m1_history_state[jj] <- "not_applicable"
      m2_history_state[jj] <- "not_applicable"
      aux_history_state[jj] <- "not_applicable"
      boundary_direction[jj] <- "terminal_or_missing_source"
      boundary_type[jj] <- "pure"
      provenance_rule[jj] <- "terminal_or_missing_source"
      next
    }

    if (identical(source_process_type[jj], "virtual_mixed_continuation_task")) {
      outcome_history_state[jj] <- "outcome"
      m1_history_state[jj] <- target_m1_state_for_world(wt)
      m2_history_state[jj] <- target_m2_state_for_world(wt)
      aux_history_state[jj] <- "not_applicable"
      boundary_direction[jj] <- "outcome_consumes_virtual_mixed_continuation"
      boundary_type[jj] <- "virtual_mixed_continuation"
      provenance_rule[jj] <- paste(
        "virtual_mixed_history_from",
        paste(
          c(
            paste0("eval=", src_key),
            paste0("outcome=", outcome_history_state[jj]),
            paste0("m1=", m1_history_state[jj]),
            paste0("m2=", m2_history_state[jj]),
            paste0("aux=", aux_history_state[jj])
          ),
          collapse = ";"
        ),
        sep = ":"
      )
      next
    }

    process_jj <- as.character(tasks$process_type[jj])
    states <- task_history_states(wt, process_jj, consume_key)
    if (identical(wt, "separate") &&
        identical(process_jj, "first_mediator_stochastic_intervention_law") &&
        identical(src_key, "second_law")) {
      states <- list(
        outcome = "outcome",
        m1 = "second_law",
        m2 = "second_law",
        auxiliary = "second_law"
      )
    }
    if (identical(src_key, "outcome")) {
      outcome_history_state[jj] <- if (identical(wt, "natural")) "outcome" else states$outcome
      m1_history_state[jj] <- states$m1
      m2_history_state[jj] <- states$m2
      aux_history_state[jj] <- states$auxiliary
    } else {
      outcome_history_state[jj] <- if (identical(wt, "natural")) src_key else states$outcome
      m1_history_state[jj] <- states$m1
      m2_history_state[jj] <- states$m2
      aux_history_state[jj] <- states$auxiliary %||% aux_state_for_source(src_key)
    }

    boundary_direction[jj] <- if (consume_key == src_key) {
      "within_state"
    } else if (consume_key == "outcome" && src_key != "outcome") {
      "outcome_consumes_law"
    } else if (consume_key != "outcome" && src_key == "outcome") {
      "law_consumes_outcome"
    } else {
      "between_law_states"
    }
    provenance_states <- unique(c(
      src_key,
      outcome_history_state[jj],
      m1_history_state[jj],
      m2_history_state[jj],
      aux_history_state[jj]
    ))
    provenance_states <- provenance_states[
      !is.na(provenance_states) &
        provenance_states != "not_applicable" &
        provenance_states != "observed_terminal_outcome"
    ]
    boundary_type[jj] <- if (length(provenance_states) > 1L) {
      "cross_regimen_mixed"
    } else {
      "pure"
    }
    provenance_rule[jj] <- if (identical(boundary_type[jj], "cross_regimen_mixed")) {
      paste(
        "mixed_history_from",
        paste(
          c(
            paste0("eval=", src_key),
            paste0("outcome=", outcome_history_state[jj]),
            paste0("m1=", m1_history_state[jj]),
            paste0("m2=", m2_history_state[jj]),
            paste0("aux=", aux_history_state[jj])
          ),
          collapse = ";"
        ),
        sep = ":"
      )
    } else {
      paste0("pure_history_from:", src_key)
    }
  }

  tasks$source_boundary_type <- boundary_type
  tasks$source_boundary_eval_state <- source_branch_key
  tasks$source_boundary_outcome_history_state <- outcome_history_state
  tasks$source_boundary_m1_history_state <- m1_history_state
  tasks$source_boundary_m2_history_state <- m2_history_state
  tasks$source_boundary_auxiliary_mediator_history_state <- aux_history_state
  tasks$source_boundary_outer_regimen <- outer_regimen
  tasks$source_boundary_m1_regimen <- m1_regimen
  tasks$source_boundary_m2_regimen <- m2_regimen
  tasks$source_boundary_direction <- boundary_direction
  tasks$source_boundary_provenance_rule <- provenance_rule
  tasks$source_boundary_consuming_state <- task_branch_key
  tasks$source_boundary_source_state <- source_branch_key
  tasks$source_boundary_source_process_type <- source_process_type
  tasks$source_boundary_source_node <- source_node
  tasks$source_boundary_source_t <- source_t

  tasks
}
