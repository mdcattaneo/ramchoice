################################################################################
# Targeted and column-generation computation for H-LAO
################################################################################

#' Structured H-LAO Preference Event
#'
#' @description
#' `hlaoEvent` describes the event that one alternative is strictly preferred
#' to every alternative in a supplied comparison set. Unlike an indicator over
#' enumerated rankings, this representation can be priced directly by the
#' H-LAO column-generation algorithm.
#'
#' @param alternative Integer identifying the focal alternative.
#' @param preferred_to Distinct integers identifying alternatives that the
#'   focal alternative must be preferred to.
#' @param name Optional event label.
#'
#' @return An object of class `ramchoiceHLAOEvent`.
#'
#' @references
#' M. D. Cattaneo, P. H. Y. Cheung, X. Ma, and Y. Masatlioglu (2026).
#' \href{https://arxiv.org/abs/2110.10650}{Attention Overload}. Working paper.
#'
#' @examples
#' hlaoEvent(4, 1, name = "4 above 1")
#' hlaoEvent(1, c(2, 3, 4), name = "1 top ranked")
#'
#' @export
hlaoEvent <- function(alternative, preferred_to, name = NULL) {
  if (!is.numeric(alternative) || length(alternative) != 1L ||
      is.na(alternative) || alternative != as.integer(alternative) ||
      alternative < 1L) {
    stop("'alternative' must be one positive integer.", call. = FALSE)
  }
  if (!is.numeric(preferred_to) || !length(preferred_to) ||
      anyNA(preferred_to) || any(preferred_to != as.integer(preferred_to)) ||
      any(preferred_to < 1L) || anyDuplicated(preferred_to)) {
    stop("'preferred_to' must contain distinct positive integers.", call. = FALSE)
  }
  alternative <- as.integer(alternative)
  preferred_to <- as.integer(preferred_to)
  if (alternative %in% preferred_to) {
    stop("The focal alternative cannot appear in 'preferred_to'.", call. = FALSE)
  }
  if (is.null(name)) {
    name <- paste0(
      alternative, " above ", paste(preferred_to, collapse = " and ")
    )
  }
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("'name' must be one nonempty character string.", call. = FALSE)
  }
  structure(
    list(
      alternative = alternative,
      preferred_to = preferred_to,
      name = name
    ),
    class = "ramchoiceHLAOEvent"
  )
}

.hlao_structured_events <- function(events, universe_size) {
  if (is.null(events)) {
    return(list())
  }
  if (inherits(events, "ramchoiceHLAOEvent")) {
    events <- list(events)
  }
  if (!is.list(events) || !length(events) ||
      !all(vapply(events, inherits, logical(1L), "ramchoiceHLAOEvent"))) {
    return(NULL)
  }
  list_names <- names(events)
  for (index in seq_along(events)) {
    event <- events[[index]]
    if (event$alternative > universe_size ||
        any(event$preferred_to > universe_size)) {
      stop("A structured event refers to an unavailable alternative.", call. = FALSE)
    }
    if (!is.null(list_names) && nzchar(list_names[index])) {
      events[[index]]$name <- list_names[index]
    }
  }
  events
}

.hlao_event_value <- function(ranking, event) {
  positions <- match(c(event$alternative, event$preferred_to), ranking)
  as.numeric(positions[1L] < min(positions[-1L]))
}

.hlao_event_matrix <- function(events, rankings) {
  if (!length(events)) {
    return(matrix(numeric(0L), nrow = nrow(rankings), ncol = 0L))
  }
  result <- vapply(events, function(event) {
    apply(rankings, 1L, .hlao_event_value, event = event)
  }, numeric(nrow(rankings)))
  if (is.null(dim(result))) {
    result <- matrix(result, ncol = 1L)
  }
  colnames(result) <- vapply(events, `[[`, character(1L), "name")
  storage.mode(result) <- "double"
  result
}

.hlao_agreement_indices <- function(agreement, n_menus) {
  if (is.null(agreement) || identical(agreement, FALSE)) {
    return(integer(0L))
  }
  if (identical(agreement, TRUE)) {
    return(seq_len(n_menus))
  }
  if (!is.numeric(agreement) || !length(agreement) || anyNA(agreement) ||
      any(agreement != as.integer(agreement)) ||
      any(agreement < 1L | agreement > n_menus)) {
    stop(
      "'agreement' must be FALSE, TRUE, or valid observed-menu indices.",
      call. = FALSE
    )
  }
  unique(as.integer(agreement))
}

.hlao_independent_spec <- function(menu, prob, domain, attention) {
  n_rows <- 1L + sum(rowSums(menu))
  rhs <- numeric(n_rows)
  rhs[1L] <- 1
  labels <- character(n_rows)
  labels[1L] <- "adding-up"
  row_lookup <- vector("list", nrow(menu))
  row_index <- 1L
  for (menu_index in seq_len(nrow(menu))) {
    items <- domain$items[[menu_index]]
    lookup <- stats::setNames(integer(length(items)), as.character(items))
    for (alternative in items) {
      row_index <- row_index + 1L
      rhs[row_index] <- prob[menu_index, alternative]
      labels[row_index] <- paste0(domain$keys[menu_index], ":", alternative)
      lookup[as.character(alternative)] <- row_index
    }
    row_lookup[[menu_index]] <- lookup
  }
  list(
    menu = menu,
    prob = prob,
    domain = domain,
    attention = attention,
    b = rhs,
    labels = labels,
    row_lookup = row_lookup,
    n_rows = n_rows,
    universe_size = ncol(menu)
  )
}

.hlao_independent_column <- function(ranking, spec) {
  column <- numeric(spec$n_rows)
  column[1L] <- 1
  ranking_positions <- integer(spec$universe_size)
  ranking_positions[ranking] <- seq_along(ranking)
  for (menu_index in seq_len(nrow(spec$menu))) {
    items <- spec$domain$items[[menu_index]]
    masses <- spec$attention$masses[[menu_index]]
    for (prefix_size in seq_along(items)) {
      prefix <- items[seq_len(prefix_size)]
      winner <- prefix[which.min(ranking_positions[prefix])]
      row <- spec$row_lookup[[menu_index]][as.character(winner)]
      column[row] <- column[row] + masses[prefix_size + 1L]
    }
  }
  column
}

.hlao_agreement_value <- function(ranking, spec, menu_index) {
  items <- spec$domain$items[[menu_index]]
  winner <- ranking[match(TRUE, ranking %in% items)]
  list_position <- match(winner, items)
  spec$attention$reach[[menu_index]][list_position]
}

.hlao_target_value <- function(target, ranking, spec) {
  if (target$kind == "event") {
    return(.hlao_event_value(ranking, target$event))
  }
  .hlao_agreement_value(ranking, spec, target$menu_index)
}

.hlao_column_targets <- function(events, agreement_indices, domain) {
  targets <- lapply(events, function(event) {
    list(kind = "event", name = event$name, event = event)
  })
  for (menu_index in agreement_indices) {
    targets[[length(targets) + 1L]] <- list(
      kind = "agreement",
      name = paste0("agreement:", menu_index),
      menu_index = menu_index,
      menu = .hlao_menu_label(domain$items[[menu_index]])
    )
  }
  targets
}

.hlao_pricing_base <- function(spec) {
  n <- spec$universe_size
  z_index <- matrix(NA_integer_, n, n)
  n_variables <- 0L
  for (a in seq_len(n)) {
    for (b in seq_len(n)) {
      if (a != b) {
        n_variables <- n_variables + 1L
        z_index[a, b] <- n_variables
      }
    }
  }
  y_layout <- vector("list", nrow(spec$menu))
  for (menu_index in seq_len(nrow(spec$menu))) {
    items <- spec$domain$items[[menu_index]]
    y_layout[[menu_index]] <- vector("list", length(items))
    for (prefix_size in seq_along(items)) {
      prefix <- items[seq_len(prefix_size)]
      indices <- n_variables + seq_along(prefix)
      n_variables <- max(indices)
      names(indices) <- as.character(prefix)
      y_layout[[menu_index]][[prefix_size]] <- indices
    }
  }

  rows <- list()
  directions <- character(0L)
  rhs <- numeric(0L)
  add_constraint <- function(row, direction, value) {
    rows[[length(rows) + 1L]] <<- row
    directions <<- c(directions, direction)
    rhs <<- c(rhs, value)
  }
  for (a in seq_len(n - 1L)) {
    for (b in (a + 1L):n) {
      row <- numeric(n_variables)
      row[z_index[a, b]] <- 1
      row[z_index[b, a]] <- 1
      add_constraint(row, "=", 1)
    }
  }
  if (n >= 3L) {
    for (a in seq_len(n)) {
      for (b in setdiff(seq_len(n), a)) {
        for (c in setdiff(seq_len(n), c(a, b))) {
          row <- numeric(n_variables)
          row[z_index[a, b]] <- 1
          row[z_index[b, c]] <- 1
          row[z_index[a, c]] <- -1
          add_constraint(row, "<=", 1)
        }
      }
    }
  }
  for (menu_index in seq_len(nrow(spec$menu))) {
    items <- spec$domain$items[[menu_index]]
    for (prefix_size in seq_along(items)) {
      prefix <- items[seq_len(prefix_size)]
      y <- y_layout[[menu_index]][[prefix_size]]
      row <- numeric(n_variables)
      row[y] <- 1
      add_constraint(row, "=", 1)
      if (length(prefix) >= 2L) {
        for (a in prefix) {
          for (b in setdiff(prefix, a)) {
            row <- numeric(n_variables)
            row[y[as.character(a)]] <- 1
            row[z_index[a, b]] <- -1
            add_constraint(row, "<=", 0)
          }
        }
      }
    }
  }
  list(
    A = do.call(rbind, rows),
    directions = directions,
    b = rhs,
    n_variables = n_variables,
    z_index = z_index,
    y_layout = y_layout
  )
}

.hlao_price_ranking <- function(base, spec, row_weights,
                                target = NULL, target_weight = 0,
                                direction = c("min", "max"),
                                exclude = NULL) {
  direction <- match.arg(direction)
  A <- base$A
  directions <- base$directions
  rhs <- base$b
  n_variables <- base$n_variables
  objective <- numeric(n_variables)

  for (menu_index in seq_len(nrow(spec$menu))) {
    items <- spec$domain$items[[menu_index]]
    masses <- spec$attention$masses[[menu_index]]
    for (prefix_size in seq_along(items)) {
      y <- base$y_layout[[menu_index]][[prefix_size]]
      for (alternative in items[seq_len(prefix_size)]) {
        row <- spec$row_lookup[[menu_index]][as.character(alternative)]
        objective[y[as.character(alternative)]] <-
          objective[y[as.character(alternative)]] +
          row_weights[row] * masses[prefix_size + 1L]
      }
    }
  }

  event_index <- NA_integer_
  if (!is.null(target) && target$kind == "event" && target_weight != 0) {
    event_index <- n_variables + 1L
    n_variables <- event_index
    A <- cbind(A, 0)
    objective <- c(objective, target_weight)
    event <- target$event
    comparison_indices <- base$z_index[
      event$alternative,
      event$preferred_to
    ]
    for (index in comparison_indices) {
      row <- numeric(n_variables)
      row[event_index] <- 1
      row[index] <- -1
      A <- rbind(A, row)
      directions <- c(directions, "<=")
      rhs <- c(rhs, 0)
    }
    row <- numeric(n_variables)
    row[event_index] <- 1
    row[comparison_indices] <- -1
    A <- rbind(A, row)
    directions <- c(directions, ">=")
    rhs <- c(rhs, -(length(comparison_indices) - 1L))
  } else if (!is.null(target) && target$kind == "agreement" &&
             target_weight != 0) {
    menu_index <- target$menu_index
    items <- spec$domain$items[[menu_index]]
    y <- base$y_layout[[menu_index]][[length(items)]]
    reach <- spec$attention$reach[[menu_index]]
    for (position in seq_along(items)) {
      objective[y[as.character(items[position])]] <-
        objective[y[as.character(items[position])]] +
        target_weight * reach[position]
    }
  }

  if (!is.null(exclude) && nrow(exclude)) {
    A <- cbind(A, matrix(0, nrow(A), n_variables - ncol(A)))
    for (ranking_index in seq_len(nrow(exclude))) {
      ranking <- exclude[ranking_index, ]
      positions <- integer(spec$universe_size)
      positions[ranking] <- seq_along(ranking)
      row <- numeric(n_variables)
      n_ones <- 0L
      for (a in seq_len(spec$universe_size)) {
        for (b in setdiff(seq_len(spec$universe_size), a)) {
          if (positions[a] < positions[b]) {
            row[base$z_index[a, b]] <- 1
            n_ones <- n_ones + 1L
          } else {
            row[base$z_index[a, b]] <- -1
          }
        }
      }
      A <- rbind(A, row)
      directions <- c(directions, "<=")
      rhs <- c(rhs, n_ones - 1L)
    }
  }

  fit <- lpSolve::lp(
    direction,
    objective,
    A,
    directions,
    rhs,
    all.bin = TRUE
  )
  if (fit$status != 0L) {
    return(list(status = fit$status, ranking = NULL, score = NA_real_))
  }
  z_solution <- fit$solution[base$z_index]
  dim(z_solution) <- dim(base$z_index)
  diag(z_solution) <- NA_real_
  wins <- rowSums(z_solution, na.rm = TRUE)
  ranking <- order(wins, decreasing = TRUE)
  column <- .hlao_independent_column(ranking, spec)
  target_value <- if (is.null(target)) 0 else {
    .hlao_target_value(target, ranking, spec)
  }
  list(
    status = 0L,
    ranking = ranking,
    column = column,
    target_value = target_value,
    score = sum(row_weights * column) + target_weight * target_value
  )
}

.hlao_phase_master <- function(columns, rhs) {
  n_rows <- length(rhs)
  matrix <- cbind(columns, diag(n_rows), -diag(n_rows))
  objective <- c(rep(0, ncol(columns)), rep(1, 2L * n_rows))
  lpSolve::lp(
    "min", objective, matrix, rep("=", n_rows), rhs,
    compute.sens = 1
  )
}

.hlao_restricted_master <- function(columns, rhs, objective) {
  lpSolve::lp(
    "min", objective, columns, rep("=", length(rhs)), rhs,
    compute.sens = 1
  )
}

.hlao_solver_status_class <- function(status) {
  if (length(status) != 1L || is.na(status)) {
    "not-run"
  } else if (status == 0L) {
    "success"
  } else if (status == 2L) {
    "infeasible"
  } else {
    "solver-error"
  }
}

.hlao_primal_residual <- function(columns, solution, rhs) {
  if (is.null(solution) || length(solution) != ncol(columns)) {
    return(NA_real_)
  }
  max(abs(drop(columns %*% solution) - rhs))
}

.hlao_column_generation <- function(menu, prob, domain, attention, targets,
                                    tolerance, max_iterations) {
  spec <- .hlao_independent_spec(menu, prob, domain, attention)
  base <- .hlao_pricing_base(spec)
  rankings <- matrix(seq_len(spec$universe_size), nrow = 1L)
  columns <- matrix(
    .hlao_independent_column(rankings[1L, ], spec),
    ncol = 1L
  )
  iterations <- 0L
  add_column <- function(priced) {
    rankings <<- rbind(rankings, priced$ranking)
    columns <<- cbind(columns, priced$column)
  }

  phase_fit <- NULL
  phase_pricing_status <- NA_integer_
  phase_pricing_score <- NA_real_
  phase_primal_residual <- NA_real_
  phase_termination <- "not-run"
  iteration_limit_reached <- FALSE
  repeat {
    if (iterations >= max_iterations) {
      iteration_limit_reached <- TRUE
      phase_termination <- "iteration-limit"
      break
    }
    iterations <- iterations + 1L
    phase_fit <- .hlao_phase_master(columns, spec$b)
    if (phase_fit$status != 0L) {
      phase_termination <- "master-solver-error"
      break
    }
    phase_primal_residual <- .hlao_primal_residual(
      cbind(columns, diag(spec$n_rows), -diag(spec$n_rows)),
      phase_fit$solution,
      spec$b
    )
    dual <- phase_fit$duals[seq_len(spec$n_rows)]
    priced <- .hlao_price_ranking(
      base, spec, row_weights = dual,
      direction = "max", exclude = rankings
    )
    phase_pricing_status <- priced$status
    if (priced$status == 2L) {
      phase_termination <- "all-rankings-generated"
      break
    }
    if (priced$status != 0L) {
      phase_termination <- "pricing-solver-error"
      break
    }
    phase_pricing_score <- priced$score
    if (priced$score <= tolerance) {
      phase_termination <- "pricing-optimal"
      break
    }
    add_column(priced)
  }

  phase_pricing_certified <- phase_pricing_status %in% c(0L, 2L) &&
    (phase_pricing_status == 2L || phase_pricing_score <= tolerance)
  phase_certified <- !iteration_limit_reached && !is.null(phase_fit) &&
    phase_fit$status == 0L && phase_pricing_certified &&
    is.finite(phase_primal_residual) && phase_primal_residual <= tolerance
  if (!phase_certified) {
    feasible <- NA
    feasibility <- phase_fit
    status <- if (!is.null(phase_fit) && phase_fit$status != 0L) {
      phase_fit$status
    } else if (!is.na(phase_pricing_status) &&
               !phase_pricing_status %in% c(0L, 2L)) {
      phase_pricing_status
    } else {
      NA_integer_
    }
  } else if (phase_fit$objval > tolerance) {
    feasible <- FALSE
    feasibility <- phase_fit
    status <- 2L
  } else {
    feasibility <- .hlao_restricted_master(
      columns, spec$b, numeric(ncol(columns))
    )
    feasibility_class <- .hlao_solver_status_class(feasibility$status)
    feasible <- if (feasibility_class == "success") {
      TRUE
    } else if (feasibility_class == "infeasible") {
      FALSE
    } else {
      NA
    }
    status <- feasibility$status
  }

  bounds <- list()
  if (isTRUE(feasible) && length(targets)) {
    for (target_index in seq_along(targets)) {
      target <- targets[[target_index]]
      endpoints <- rep(NA_real_, 2L)
      statuses <- rep(NA_integer_, 2L)
      pricing_statuses <- rep(NA_integer_, 2L)
      reduced_costs <- rep(NA_real_, 2L)
      dual_residuals <- rep(NA_real_, 2L)
      optimality_gap_bounds <- rep(NA_real_, 2L)
      primal_residuals <- rep(NA_real_, 2L)
      certified <- rep(FALSE, 2L)
      for (endpoint in seq_len(2L)) {
        sign <- if (endpoint == 1L) 1 else -1
        repeat {
          if (iterations >= max_iterations) {
            iteration_limit_reached <- TRUE
            break
          }
          iterations <- iterations + 1L
          values <- apply(rankings, 1L, .hlao_target_value,
                          target = target, spec = spec)
          master <- .hlao_restricted_master(columns, spec$b, sign * values)
          if (master$status != 0L) {
            statuses[endpoint] <- master$status
            endpoints[endpoint] <- NA_real_
            break
          }
          primal_residuals[endpoint] <- .hlao_primal_residual(
            columns, master$solution, spec$b
          )
          dual <- master$duals[seq_len(spec$n_rows)]
          priced <- .hlao_price_ranking(
            base, spec,
            row_weights = -dual,
            target = target,
            target_weight = sign,
            direction = "min",
            exclude = rankings
          )
          pricing_statuses[endpoint] <- priced$status
          if (priced$status == 2L) {
            reduced_costs[endpoint] <- Inf
            dual_residuals[endpoint] <- 0
            optimality_gap_bounds[endpoint] <- 0
            certified[endpoint] <- is.finite(primal_residuals[endpoint]) &&
              primal_residuals[endpoint] <= tolerance
            statuses[endpoint] <- if (certified[endpoint]) 0L else NA_integer_
            endpoints[endpoint] <- if (certified[endpoint]) {
              sign * master$objval
            } else {
              NA_real_
            }
            break
          }
          if (priced$status != 0L) {
            statuses[endpoint] <- priced$status
            endpoints[endpoint] <- NA_real_
            break
          }
          reduced_cost <- priced$score
          reduced_costs[endpoint] <- reduced_cost
          dual_residuals[endpoint] <- max(0, -reduced_cost)
          optimality_gap_bounds[endpoint] <- dual_residuals[endpoint]
          if (reduced_cost >= -tolerance) {
            certified[endpoint] <- is.finite(primal_residuals[endpoint]) &&
              primal_residuals[endpoint] <= tolerance
            statuses[endpoint] <- if (certified[endpoint]) 0L else NA_integer_
            endpoints[endpoint] <- if (certified[endpoint]) {
              sign * master$objval
            } else {
              NA_real_
            }
            break
          }
          add_column(priced)
        }
      }
      bounds[[target_index]] <- data.frame(
        target = target$name,
        kind = target$kind,
        menu_id = if (target$kind == "agreement") target$menu_index else NA_integer_,
        menu = if (target$kind == "agreement") target$menu else NA_character_,
        lower = endpoints[1L],
        upper = endpoints[2L],
        lower_status = statuses[1L],
        upper_status = statuses[2L],
        lower_pricing_status = pricing_statuses[1L],
        upper_pricing_status = pricing_statuses[2L],
        lower_reduced_cost = reduced_costs[1L],
        upper_reduced_cost = reduced_costs[2L],
        lower_dual_residual = dual_residuals[1L],
        upper_dual_residual = dual_residuals[2L],
        lower_optimality_gap_bound = optimality_gap_bounds[1L],
        upper_optimality_gap_bound = optimality_gap_bounds[2L],
        lower_primal_residual = primal_residuals[1L],
        upper_primal_residual = primal_residuals[2L],
        lower_certified = certified[1L],
        upper_certified = certified[2L],
        stringsAsFactors = FALSE
      )
    }
  }

  bound_table <- if (length(bounds)) do.call(rbind, bounds) else {
    data.frame(
      target = character(0L), kind = character(0L),
      menu_id = integer(0L), menu = character(0L),
      lower = numeric(0L), upper = numeric(0L),
      lower_status = integer(0L), upper_status = integer(0L),
      lower_pricing_status = integer(0L), upper_pricing_status = integer(0L),
      lower_reduced_cost = numeric(0L), upper_reduced_cost = numeric(0L),
      lower_dual_residual = numeric(0L), upper_dual_residual = numeric(0L),
      lower_optimality_gap_bound = numeric(0L),
      upper_optimality_gap_bound = numeric(0L),
      lower_primal_residual = numeric(0L), upper_primal_residual = numeric(0L),
      lower_certified = logical(0L), upper_certified = logical(0L),
      stringsAsFactors = FALSE
    )
  }
  endpoint_reduced_costs <- c(
    bound_table$lower_reduced_cost, bound_table$upper_reduced_cost
  )
  endpoint_residuals <- c(
    bound_table$lower_primal_residual, bound_table$upper_primal_residual
  )
  endpoint_dual_residuals <- c(
    bound_table$lower_dual_residual, bound_table$upper_dual_residual
  )
  endpoint_gap_bounds <- c(
    bound_table$lower_optimality_gap_bound,
    bound_table$upper_optimality_gap_bound
  )
  all_endpoints_certified <- if (nrow(bound_table)) {
    all(bound_table$lower_certified & bound_table$upper_certified)
  } else {
    isTRUE(feasible)
  }

  list(
    feasible = feasible,
    status = status,
    rankings = rankings,
    columns = columns,
    bounds = bound_table,
    diagnostics = data.frame(
      algorithm = "column-generation",
      n_columns = nrow(rankings),
      n_possible_rankings = gamma(spec$universe_size + 1L),
      iterations = iterations,
      pricing_binary_variables = base$n_variables,
      pricing_constraints = nrow(base$A),
      tolerance = tolerance,
      phase_master_status = if (is.null(phase_fit)) NA_integer_ else phase_fit$status,
      phase_pricing_status = phase_pricing_status,
      phase_objective = if (is.null(phase_fit) || phase_fit$status != 0L) NA_real_ else phase_fit$objval,
      phase_pricing_score = phase_pricing_score,
      phase_primal_residual = phase_primal_residual,
      phase_termination = phase_termination,
      phase_certified = phase_certified,
      iteration_limit_reached = iteration_limit_reached,
      feasibility_status = status,
      feasibility_status_class = .hlao_solver_status_class(status),
      endpoint_count = 2L * nrow(bound_table),
      minimum_final_reduced_cost = if (length(endpoint_reduced_costs) &&
                                         any(!is.na(endpoint_reduced_costs))) {
        min(endpoint_reduced_costs, na.rm = TRUE)
      } else {
        NA_real_
      },
      maximum_primal_residual = if (length(endpoint_residuals) &&
                                       any(!is.na(endpoint_residuals))) {
        max(endpoint_residuals, na.rm = TRUE)
      } else {
        NA_real_
      },
      maximum_dual_residual = if (length(endpoint_dual_residuals) &&
                                     any(!is.na(endpoint_dual_residuals))) {
        max(endpoint_dual_residuals, na.rm = TRUE)
      } else {
        NA_real_
      },
      maximum_optimality_gap_bound = if (length(endpoint_gap_bounds) &&
                                            any(!is.na(endpoint_gap_bounds))) {
        max(endpoint_gap_bounds, na.rm = TRUE)
      } else {
        NA_real_
      },
      all_endpoints_certified = all_endpoints_certified,
      certified = phase_certified && isTRUE(feasible) && all_endpoints_certified,
      stringsAsFactors = FALSE
    ),
    spec = spec
  )
}

.hlao_agreement_objective <- function(system, rankings, domain, attention,
                                      menu_index, mode) {
  n_rankings <- nrow(rankings)
  objective <- numeric(system$n_variables)
  positions <- .hlao_ranking_positions(rankings, ncol(rankings))
  items <- domain$items[[menu_index]]
  winners <- items[max.col(-positions[, items, drop = FALSE], ties.method = "first")]
  list_positions <- match(winners, items)
  if (mode == "independent") {
    objective[seq_len(n_rankings)] <- attention$reach[[menu_index]][list_positions]
    return(objective)
  }
  q <- system$layout[[menu_index]]$q
  for (ranking_index in seq_len(n_rankings)) {
    first_agreeing_stop <- list_positions[ranking_index]
    objective[q[ranking_index, (first_agreeing_stop + 1L):ncol(q)]] <- 1
  }
  objective
}

.hlao_solve_objective <- function(system, objective) {
  directions <- if (is.null(system$directions)) {
    rep("=", length(system$b))
  } else {
    system$directions
  }
  lower <- lpSolve::lp("min", objective, system$A, directions, system$b)
  upper <- lpSolve::lp("max", objective, system$A, directions, system$b)
  c(
    lower = if (lower$status == 0L) lower$objval else NA_real_,
    upper = if (upper$status == 0L) upper$objval else NA_real_,
    lower_status = lower$status,
    upper_status = upper$status
  )
}

.hlao_direct_agreement <- function(full_attention, attention, domain,
                                   agreement_indices) {
  if (is.null(full_attention) || !length(agreement_indices)) {
    return(NULL)
  }
  do.call(rbind, lapply(agreement_indices, function(menu_index) {
    items <- domain$items[[menu_index]]
    point <- sum(
      attention$reach[[menu_index]] *
        full_attention$matrix[menu_index, items]
    )
    data.frame(
      menu_id = menu_index,
      menu = .hlao_menu_label(items),
      mode = "independent",
      lower = point,
      upper = point,
      point_identified = TRUE,
      method = "full-attention-formula",
      lower_status = 0L,
      upper_status = 0L,
      stringsAsFactors = FALSE
    )
  }))
}

.hlao_enumerated_agreement <- function(fits, modes, rankings, domain, attention,
                                       agreement_indices,
                                       full_attention = NULL) {
  if (!length(agreement_indices)) {
    return(data.frame(
      menu_id = integer(0L), menu = character(0L), mode = character(0L),
      lower = numeric(0L), upper = numeric(0L),
      point_identified = logical(0L), method = character(0L),
      lower_status = integer(0L), upper_status = integer(0L),
      stringsAsFactors = FALSE
    ))
  }
  rows <- list()
  for (mode in modes) {
    system <- fits[[mode]]$system
    if (is.null(system) || !fits[[mode]]$fit$feasible) {
      next
    }
    for (menu_index in agreement_indices) {
      if (mode == "independent" && !is.null(full_attention)) {
        rows[[length(rows) + 1L]] <- .hlao_direct_agreement(
          full_attention, attention, domain, menu_index
        )
        next
      }
      objective <- .hlao_agreement_objective(
        system, rankings, domain, attention, menu_index, mode
      )
      endpoint <- .hlao_solve_objective(system, objective)
      rows[[length(rows) + 1L]] <- data.frame(
        menu_id = menu_index,
        menu = .hlao_menu_label(domain$items[[menu_index]]),
        mode = mode,
        lower = endpoint["lower"],
        upper = endpoint["upper"],
        point_identified = isTRUE(
          abs(endpoint["upper"] - endpoint["lower"]) <= 1e-8
        ),
        method = "enumerated-linear-program",
        lower_status = as.integer(endpoint["lower_status"]),
        upper_status = as.integer(endpoint["upper_status"]),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows)) do.call(rbind, rows) else NULL
}
