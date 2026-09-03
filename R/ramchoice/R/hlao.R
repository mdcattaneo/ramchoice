################################################################################
# Heterogeneous List-Based Attention Overload Model
################################################################################

.hlao_menu_key <- function(menu_row) {
  paste(as.integer(menu_row), collapse = "")
}

.hlao_menu_label <- function(items) {
  paste(items, collapse = ",")
}

.hlao_validate_tolerance <- function(tolerance) {
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance < 0) {
    stop("'tolerance' must be one finite nonnegative number.", call. = FALSE)
  }
  tolerance
}

.hlao_validate_list_order <- function(list_order, universe_size) {
  if (is.null(list_order)) {
    return(seq_len(universe_size))
  }
  if (!is.numeric(list_order) || length(list_order) != universe_size ||
      anyNA(list_order) ||
      !identical(sort(as.integer(list_order)), seq_len(universe_size))) {
    stop(
      "'list_order' must be a permutation of the alternative indices.",
      call. = FALSE
    )
  }
  as.integer(list_order)
}

.hlao_validate_menu <- function(menu, allow_duplicates = FALSE) {
  if (!is.matrix(menu) || !is.numeric(menu) && !is.logical(menu) ||
      nrow(menu) < 1L || ncol(menu) < 1L || anyNA(menu) ||
      any(!(menu %in% c(0, 1)))) {
    stop("'menu' must be a nonempty numeric matrix of zeros and ones.", call. = FALSE)
  }
  menu <- matrix(
    as.integer(menu),
    nrow = nrow(menu),
    ncol = ncol(menu),
    dimnames = dimnames(menu)
  )
  if (any(rowSums(menu) == 0L)) {
    stop("Every observed menu must contain at least one alternative.", call. = FALSE)
  }
  keys <- apply(menu, 1L, .hlao_menu_key)
  if (!allow_duplicates && anyDuplicated(keys)) {
    stop("Population inputs must contain one row per distinct menu.", call. = FALSE)
  }
  menu
}

.hlao_prepare_population <- function(menu, prob, outside_prob, tolerance) {
  menu <- .hlao_validate_menu(menu)
  if (!is.matrix(prob) || !is.numeric(prob) ||
      !identical(dim(prob), dim(menu)) || anyNA(prob) ||
      any(!is.finite(prob))) {
    stop("'prob' must be a finite numeric matrix with the same dimensions as 'menu'.", call. = FALSE)
  }
  if (any(abs(prob[menu == 0L]) > tolerance)) {
    stop("Choice probabilities must be zero for unavailable alternatives.", call. = FALSE)
  }
  if (any(prob < -tolerance) || any(prob > 1 + tolerance)) {
    stop("Choice probabilities must lie between zero and one.", call. = FALSE)
  }
  if (is.null(outside_prob)) {
    outside_prob <- 1 - rowSums(prob)
  }
  if (!is.numeric(outside_prob) || length(outside_prob) != nrow(menu) ||
      anyNA(outside_prob) || any(!is.finite(outside_prob)) ||
      any(outside_prob < -tolerance) || any(outside_prob > 1 + tolerance)) {
    stop("'outside_prob' must contain one probability between zero and one per menu.", call. = FALSE)
  }
  if (any(abs(rowSums(prob) + outside_prob - 1) > tolerance)) {
    stop("Inside and outside choice probabilities must add to one within each menu.", call. = FALSE)
  }

  prob <- matrix(
    pmin(1, pmax(0, as.numeric(prob))),
    nrow = nrow(prob),
    ncol = ncol(prob),
    dimnames = dimnames(prob)
  )
  outside_prob <- pmin(1, pmax(0, as.numeric(outside_prob)))
  list(menu = menu, prob = prob, outside_prob = outside_prob)
}

.hlao_domain <- function(menu, list_order) {
  keys <- apply(menu, 1L, .hlao_menu_key)
  lookup <- stats::setNames(seq_len(nrow(menu)), keys)
  items <- lapply(seq_len(nrow(menu)), function(index) {
    list_order[menu[index, list_order] == 1L]
  })
  suffix_index <- vector("list", nrow(menu))
  missing <- character(0L)

  for (index in seq_len(nrow(menu))) {
    current <- items[[index]]
    suffix_index[[index]] <- vapply(seq_along(current), function(position) {
      suffix <- integer(ncol(menu))
      suffix[current[position:length(current)]] <- 1L
      suffix_key <- .hlao_menu_key(suffix)
      matched <- unname(lookup[suffix_key])
      if (is.na(matched)) {
        missing <<- c(missing, suffix_key)
      }
      as.integer(matched)
    }, integer(1L))
  }

  if (length(missing)) {
    stop(
      "The observed menu domain must be suffix closed; missing suffix keys: ",
      paste(unique(missing), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  list(keys = keys, lookup = lookup, items = items, suffix_index = suffix_index)
}

.hlao_domain_unrestricted <- function(menu, list_order) {
  keys <- apply(menu, 1L, .hlao_menu_key)
  list(
    keys = keys,
    lookup = stats::setNames(seq_len(nrow(menu)), keys),
    items = lapply(seq_len(nrow(menu)), function(index) {
      list_order[menu[index, list_order] == 1L]
    }),
    suffix_index = NULL
  )
}

.hlao_recover_attention <- function(menu, outside_prob, domain, tolerance) {
  n_menu <- nrow(menu)
  reach <- masses <- vector("list", n_menu)
  attention_rows <- vector("list", n_menu)

  for (index in seq_len(n_menu)) {
    suffix <- domain$suffix_index[[index]]
    continuation <- 1 - outside_prob[suffix]
    reach[[index]] <- cumprod(continuation)
    current_reach <- reach[[index]]
    masses[[index]] <- c(
      1 - current_reach[1L],
      if (length(current_reach) > 1L) {
        current_reach[-length(current_reach)] - current_reach[-1L]
      } else {
        numeric(0L)
      },
      current_reach[length(current_reach)]
    )
    attention_rows[[index]] <- data.frame(
      menu_id = index,
      menu = .hlao_menu_label(domain$items[[index]]),
      alternative = domain$items[[index]],
      position = seq_along(domain$items[[index]]),
      reach = current_reach,
      stop_mass = masses[[index]][-1L],
      stringsAsFactors = FALSE
    )
  }

  overload <- list()
  overload_index <- 0L
  for (small in seq_len(n_menu)) {
    for (large in seq_len(n_menu)) {
      if (small == large || sum(menu[small, ]) >= sum(menu[large, ]) ||
          !all(menu[small, ] <= menu[large, ])) {
        next
      }
      for (alternative in domain$items[[small]]) {
        small_position <- match(alternative, domain$items[[small]])
        large_position <- match(alternative, domain$items[[large]])
        overload_index <- overload_index + 1L
        overload[[overload_index]] <- data.frame(
          smaller_menu_id = small,
          larger_menu_id = large,
          alternative = alternative,
          smaller_reach = reach[[small]][small_position],
          larger_reach = reach[[large]][large_position],
          violation = reach[[large]][large_position] - reach[[small]][small_position]
        )
      }
    }
  }
  overload <- if (length(overload)) {
    do.call(rbind, overload)
  } else {
    data.frame(
      smaller_menu_id = integer(0L),
      larger_menu_id = integer(0L),
      alternative = integer(0L),
      smaller_reach = numeric(0L),
      larger_reach = numeric(0L),
      violation = numeric(0L)
    )
  }

  mass_values <- unlist(masses, use.names = FALSE)
  max_overload <- if (nrow(overload)) max(overload$violation) else -Inf
  valid <- min(mass_values) >= -tolerance && max_overload <= tolerance
  diagnostics <- data.frame(
    valid = valid,
    min_prefix_mass = min(mass_values),
    n_negative_prefix_masses = sum(mass_values < -tolerance),
    max_attention_overload_violation = if (is.finite(max_overload)) {
      max(c(max_overload, 0))
    } else {
      0
    },
    n_attention_overload_violations = sum(overload$violation > tolerance)
  )

  list(
    reach = reach,
    masses = masses,
    attention = do.call(rbind, attention_rows),
    overload = overload,
    diagnostics = diagnostics
  )
}

#' Enumerate Strict Preference Rankings
#'
#' @description
#' `hlaoRankings` returns all strict rankings of a supplied alternative set in
#' the deterministic ordering used by the H-LAO linear-programming routines.
#'
#' @param alternatives Vector of distinct alternative labels.
#'
#' @return A matrix with one strict ranking per row, from best to worst.
#'
#' @references
#' M. D. Cattaneo, P. H. Y. Cheung, X. Ma, and Y. Masatlioglu (2026).
#' \href{https://arxiv.org/abs/2110.10650}{Attention Overload}. Working paper.
#'
#' @examples
#' hlaoRankings(1:3)
#'
#' @export
hlaoRankings <- function(alternatives) {
  if (length(alternatives) < 1L || anyNA(alternatives) ||
      anyDuplicated(alternatives)) {
    stop("'alternatives' must contain distinct nonmissing labels.", call. = FALSE)
  }
  recurse <- function(values) {
    if (length(values) == 1L) {
      return(matrix(values, nrow = 1L))
    }
    do.call(rbind, lapply(seq_along(values), function(index) {
      cbind(values[index], recurse(values[-index]))
    }))
  }
  result <- recurse(alternatives)
  colnames(result) <- paste0("rank", seq_len(ncol(result)))
  result
}

.hlao_normalize_events <- function(events, rankings) {
  if (is.null(events)) {
    return(matrix(numeric(0L), nrow = nrow(rankings), ncol = 0L))
  }
  structured <- .hlao_structured_events(events, ncol(rankings))
  if (!is.null(structured)) {
    return(.hlao_event_matrix(structured, rankings))
  }
  if (is.list(events) && !is.data.frame(events)) {
    if (!length(events)) {
      return(matrix(numeric(0L), nrow = nrow(rankings), ncol = 0L))
    }
    event_names <- names(events)
    if (is.null(event_names) || any(!nzchar(event_names))) {
      event_names <- paste0("event", seq_along(events))
    }
    events <- do.call(cbind, events)
    colnames(events) <- event_names
  } else if (is.vector(events) && !is.list(events)) {
    event_name <- names(events)[1L]
    events <- matrix(events, ncol = 1L)
    colnames(events) <- if (!is.null(event_name) && nzchar(event_name)) {
      event_name
    } else {
      "event1"
    }
  } else {
    events <- as.matrix(events)
  }
  if (nrow(events) != nrow(rankings) || !is.numeric(events) && !is.logical(events) ||
      anyNA(events) || any(!(events %in% c(0, 1)))) {
    stop(
      "'events' must provide a zero-one indicator for every enumerated ranking.",
      call. = FALSE
    )
  }
  storage.mode(events) <- "double"
  if (is.null(colnames(events))) {
    colnames(events) <- paste0("event", seq_len(ncol(events)))
  }
  events
}

.hlao_ranking_positions <- function(rankings, universe_size) {
  positions <- matrix(0L, nrow = nrow(rankings), ncol = universe_size)
  for (index in seq_len(nrow(rankings))) {
    positions[index, rankings[index, ]] <- seq_len(universe_size)
  }
  positions
}

.hlao_independent_system <- function(menu, prob, domain, attention, rankings) {
  positions <- .hlao_ranking_positions(rankings, ncol(menu))
  rows <- vector("list", 1L + sum(rowSums(menu)))
  rhs <- numeric(length(rows))
  labels <- character(length(rows))
  rows[[1L]] <- rep(1, nrow(rankings))
  rhs[1L] <- 1
  labels[1L] <- "adding-up"
  row_index <- 1L

  for (menu_index in seq_len(nrow(menu))) {
    items <- domain$items[[menu_index]]
    coefficients <- matrix(0, nrow = length(items), ncol = nrow(rankings))
    for (prefix_size in seq_along(items)) {
      prefix <- items[seq_len(prefix_size)]
      winner_position <- max.col(-positions[, prefix, drop = FALSE], ties.method = "first")
      winners <- prefix[winner_position]
      mass <- attention$masses[[menu_index]][prefix_size + 1L]
      for (alternative_index in seq_along(items)) {
        coefficients[alternative_index, ] <-
          coefficients[alternative_index, ] +
          mass * as.numeric(winners == items[alternative_index])
      }
    }
    for (alternative_index in seq_along(items)) {
      row_index <- row_index + 1L
      alternative <- items[alternative_index]
      rows[[row_index]] <- coefficients[alternative_index, ]
      rhs[row_index] <- prob[menu_index, alternative]
      labels[row_index] <- paste0(domain$keys[menu_index], ":", alternative)
    }
  }
  list(
    A = do.call(rbind, rows),
    b = rhs,
    labels = labels,
    n_rankings = nrow(rankings),
    n_variables = nrow(rankings)
  )
}

.hlao_population_robust_system <- function(menu, prob, domain, attention, rankings) {
  n_rankings <- nrow(rankings)
  positions <- .hlao_ranking_positions(rankings, ncol(menu))
  widths <- vapply(domain$items, function(items) {
    n_rankings * (length(items) + 1L)
  }, integer(1L))
  starts <- n_rankings + c(0L, utils::head(cumsum(widths), -1L)) + 1L
  n_variables <- n_rankings + sum(widths)
  layout <- vector("list", nrow(menu))
  rows <- list()
  rhs <- numeric(0L)

  add_constraint <- function(row, value) {
    rows[[length(rows) + 1L]] <<- row
    rhs <<- c(rhs, value)
  }
  row <- numeric(n_variables)
  row[seq_len(n_rankings)] <- 1
  add_constraint(row, 1)

  for (menu_index in seq_len(nrow(menu))) {
    items <- domain$items[[menu_index]]
    n_stops <- length(items) + 1L
    q_index <- function(ranking, stop) {
      starts[menu_index] + (ranking - 1L) * n_stops + stop
    }
    layout[[menu_index]] <- list(
      q = matrix(
        vapply(seq_len(n_rankings), function(ranking) {
          vapply(0:(n_stops - 1L), function(stop) {
            q_index(ranking, stop)
          }, integer(1L))
        }, integer(n_stops)),
        nrow = n_rankings,
        byrow = TRUE
      )
    )

    for (stop in 0:(n_stops - 1L)) {
      row <- numeric(n_variables)
      row[vapply(seq_len(n_rankings), q_index, integer(1L), stop = stop)] <- 1
      add_constraint(row, attention$masses[[menu_index]][stop + 1L])
    }
    for (ranking in seq_len(n_rankings)) {
      row <- numeric(n_variables)
      row[ranking] <- -1
      row[vapply(0:(n_stops - 1L), function(stop) {
        q_index(ranking, stop)
      }, integer(1L))] <- 1
      add_constraint(row, 0)
    }
    for (alternative in items) {
      row <- numeric(n_variables)
      for (prefix_size in seq_along(items)) {
        prefix <- items[seq_len(prefix_size)]
        winner_position <- max.col(
          -positions[, prefix, drop = FALSE],
          ties.method = "first"
        )
        winners <- prefix[winner_position]
        selected <- which(winners == alternative)
        if (length(selected)) {
          row[vapply(selected, q_index, integer(1L), stop = prefix_size)] <- 1
        }
      }
      add_constraint(row, prob[menu_index, alternative])
    }
  }
  list(
    A = do.call(rbind, rows),
    b = rhs,
    n_rankings = n_rankings,
    n_variables = n_variables,
    layout = layout
  )
}

.hlao_nopi_system <- function(menu, domain, rankings, prob = NULL,
                              outside_prob = NULL, bands = NULL) {
  sample_problem <- !is.null(bands)
  if (sample_problem == (!is.null(prob) || !is.null(outside_prob))) {
    stop(
      "Supply either population probabilities or sample probability bands.",
      call. = FALSE
    )
  }
  n_rankings <- nrow(rankings)
  positions <- .hlao_ranking_positions(rankings, ncol(menu))
  layout <- vector("list", nrow(menu))
  offset <- n_rankings
  for (menu_index in seq_len(nrow(menu))) {
    n_stops <- length(domain$items[[menu_index]]) + 1L
    mass <- offset + seq_len(n_stops)
    offset <- max(mass)
    q <- matrix(
      offset + seq_len(n_rankings * n_stops),
      nrow = n_rankings,
      ncol = n_stops,
      byrow = TRUE
    )
    offset <- max(q)
    layout[[menu_index]] <- list(mass = mass, q = q)
  }
  n_variables <- offset
  rows <- list()
  directions <- character(0L)
  rhs <- numeric(0L)
  add_constraint <- function(row, direction, value) {
    rows[[length(rows) + 1L]] <<- row
    directions <<- c(directions, direction)
    rhs <<- c(rhs, value)
  }

  row <- numeric(n_variables)
  row[seq_len(n_rankings)] <- 1
  add_constraint(row, "=", 1)

  for (menu_index in seq_len(nrow(menu))) {
    current <- layout[[menu_index]]
    items <- domain$items[[menu_index]]
    n_items <- length(items)

    row <- numeric(n_variables)
    row[current$mass] <- 1
    add_constraint(row, "=", 1)

    row <- numeric(n_variables)
    row[current$mass[1L]] <- 1
    if (sample_problem) {
      add_constraint(row, ">=", bands$outside_lower[menu_index])
      add_constraint(row, "<=", bands$outside_upper[menu_index])
    } else {
      add_constraint(row, "=", outside_prob[menu_index])
    }

    for (stop in 0:n_items) {
      row <- numeric(n_variables)
      row[current$q[, stop + 1L]] <- 1
      row[current$mass[stop + 1L]] <- -1
      add_constraint(row, "=", 0)
    }
    for (ranking in seq_len(n_rankings)) {
      row <- numeric(n_variables)
      row[current$q[ranking, ]] <- 1
      row[ranking] <- -1
      add_constraint(row, "=", 0)
    }
    for (alternative in items) {
      row <- numeric(n_variables)
      for (prefix_size in seq_along(items)) {
        prefix <- items[seq_len(prefix_size)]
        winner_position <- max.col(
          -positions[, prefix, drop = FALSE],
          ties.method = "first"
        )
        winners <- prefix[winner_position]
        selected <- which(winners == alternative)
        if (length(selected)) {
          row[current$q[selected, prefix_size + 1L]] <- 1
        }
      }
      if (sample_problem) {
        add_constraint(row, ">=", bands$lower[menu_index, alternative])
        add_constraint(row, "<=", bands$upper[menu_index, alternative])
      } else {
        add_constraint(row, "=", prob[menu_index, alternative])
      }
    }
  }

  for (small in seq_len(nrow(menu))) {
    for (large in seq_len(nrow(menu))) {
      if (small == large || sum(menu[small, ]) >= sum(menu[large, ]) ||
          !all(menu[small, ] <= menu[large, ])) {
        next
      }
      for (alternative in domain$items[[small]]) {
        small_position <- match(alternative, domain$items[[small]])
        large_position <- match(alternative, domain$items[[large]])
        row <- numeric(n_variables)
        row[layout[[small]]$mass[(small_position + 1L):
                                  length(layout[[small]]$mass)]] <- 1
        row[layout[[large]]$mass[(large_position + 1L):
                                  length(layout[[large]]$mass)]] <- -1
        add_constraint(row, ">=", 0)
      }
    }
  }

  list(
    A = do.call(rbind, rows),
    directions = directions,
    b = rhs,
    n_rankings = n_rankings,
    n_variables = n_variables,
    layout = layout
  )
}

.hlao_solve_general_bounds <- function(system, events, mode) {
  feasibility <- lpSolve::lp(
    "min", rep(0, system$n_variables), system$A,
    system$directions, system$b
  )
  feasible <- feasibility$status == 0L
  empty <- data.frame(
    event = character(0L), mode = character(0L),
    lower = numeric(0L), upper = numeric(0L),
    lower_status = integer(0L), upper_status = integer(0L),
    stringsAsFactors = FALSE
  )
  if (!feasible || ncol(events) == 0L) {
    return(list(feasible = feasible, status = feasibility$status, bounds = empty))
  }
  rows <- vector("list", ncol(events))
  for (event_index in seq_len(ncol(events))) {
    objective <- numeric(system$n_variables)
    objective[seq_len(nrow(events))] <- events[, event_index]
    lower <- lpSolve::lp(
      "min", objective, system$A, system$directions, system$b
    )
    upper <- lpSolve::lp(
      "max", objective, system$A, system$directions, system$b
    )
    rows[[event_index]] <- data.frame(
      event = colnames(events)[event_index],
      mode = mode,
      lower = if (lower$status == 0L) lower$objval else NA_real_,
      upper = if (upper$status == 0L) upper$objval else NA_real_,
      lower_status = lower$status,
      upper_status = upper$status,
      stringsAsFactors = FALSE
    )
  }
  list(
    feasible = feasible,
    status = feasibility$status,
    bounds = do.call(rbind, rows)
  )
}

.hlao_solve_bounds <- function(system, events, mode) {
  n_variables <- ncol(system$A)
  directions <- rep("=", length(system$b))
  feasibility <- lpSolve::lp(
    "min",
    rep(0, n_variables),
    system$A,
    directions,
    system$b
  )
  feasible <- feasibility$status == 0L
  bounds <- data.frame(
    event = character(0L),
    mode = character(0L),
    lower = numeric(0L),
    upper = numeric(0L),
    lower_status = integer(0L),
    upper_status = integer(0L),
    stringsAsFactors = FALSE
  )
  if (!feasible || ncol(events) == 0L) {
    return(list(feasible = feasible, status = feasibility$status, bounds = bounds))
  }

  result <- vector("list", ncol(events))
  for (event_index in seq_len(ncol(events))) {
    objective <- numeric(n_variables)
    objective[seq_len(nrow(events))] <- events[, event_index]
    lower <- lpSolve::lp(
      "min", objective, system$A, directions, system$b
    )
    upper <- lpSolve::lp(
      "max", objective, system$A, directions, system$b
    )
    result[[event_index]] <- data.frame(
      event = colnames(events)[event_index],
      mode = mode,
      lower = if (lower$status == 0L) lower$objval else NA_real_,
      upper = if (upper$status == 0L) upper$objval else NA_real_,
      lower_status = lower$status,
      upper_status = upper$status,
      stringsAsFactors = FALSE
    )
  }
  list(
    feasible = feasible,
    status = feasibility$status,
    bounds = do.call(rbind, result)
  )
}

.hlao_is_full_domain <- function(menu) {
  universe_size <- ncol(menu)
  if (universe_size > 20L || nrow(menu) != 2^universe_size - 1L) {
    return(FALSE)
  }
  observed <- sort(apply(menu, 1L, .hlao_menu_key))
  expected <- character(0L)
  for (size in seq_len(universe_size)) {
    combinations <- utils::combn(universe_size, size)
    for (index in seq_len(ncol(combinations))) {
      row <- integer(universe_size)
      row[combinations[, index]] <- 1L
      expected <- c(expected, .hlao_menu_key(row))
    }
  }
  identical(observed, sort(expected))
}

.hlao_full_attention <- function(menu, prob, domain, attention, tolerance) {
  if (!.hlao_is_full_domain(menu)) {
    return(NULL)
  }
  terminal <- vapply(attention$reach, function(value) value[length(value)], numeric(1L))
  if (any(terminal <= tolerance)) {
    return(NULL)
  }

  full <- matrix(0, nrow = nrow(menu), ncol = ncol(menu))
  menu_order <- order(rowSums(menu))
  for (menu_index in menu_order) {
    items <- domain$items[[menu_index]]
    if (length(items) == 1L) {
      full[menu_index, items] <- 1
      next
    }
    reach <- attention$reach[[menu_index]]
    for (position in 2:length(items)) {
      alternative <- items[position]
      adjustment <- 0
      if (position <= length(items) - 1L) {
        for (prefix_size in position:(length(items) - 1L)) {
          prefix_menu <- integer(ncol(menu))
          prefix_menu[items[seq_len(prefix_size)]] <- 1L
          prefix_index <- unname(domain$lookup[.hlao_menu_key(prefix_menu)])
          adjustment <- adjustment +
            (reach[prefix_size] - reach[prefix_size + 1L]) *
            full[prefix_index, alternative]
        }
      }
      full[menu_index, alternative] <-
        (prob[menu_index, alternative] - adjustment) /
        reach[length(items)]
    }
    full[menu_index, items[1L]] <- 1 - sum(full[menu_index, items[-1L]])
  }

  rows <- list()
  row_index <- 0L
  for (menu_index in seq_len(nrow(menu))) {
    for (alternative in domain$items[[menu_index]]) {
      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        menu_id = menu_index,
        menu = .hlao_menu_label(domain$items[[menu_index]]),
        alternative = alternative,
        probability = full[menu_index, alternative]
      )
    }
  }
  list(matrix = full, results = do.call(rbind, rows))
}

.hlao_bm_diagnostics <- function(menu, domain, full_attention, tolerance) {
  if (is.null(full_attention)) {
    return(NULL)
  }
  rows <- list()
  row_index <- 0L
  for (target in seq_len(nrow(menu))) {
    target_items <- domain$items[[target]]
    supersets <- which(apply(menu, 1L, function(row) all(row >= menu[target, ])))
    for (alternative in target_items) {
      value <- sum(vapply(supersets, function(superset) {
        (-1)^(sum(menu[superset, ]) - sum(menu[target, ])) *
          full_attention$matrix[superset, alternative]
      }, numeric(1L)))
      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        menu_id = target,
        menu = .hlao_menu_label(target_items),
        alternative = alternative,
        bm_polynomial = value,
        violated = value < -tolerance
      )
    }
  }
  results <- do.call(rbind, rows)
  list(
    results = results,
    valid = !any(results$violated),
    minimum = min(results$bm_polynomial),
    n_violated = sum(results$violated)
  )
}

.hlao_population_pairwise <- function(menu, prob, domain, attention, tolerance) {
  binary <- which(rowSums(menu) == 2L)
  if (!length(binary)) {
    return(data.frame(
      earlier = integer(0L), later = integer(0L), reach = numeric(0L),
      share_later_preferred = numeric(0L), identified = logical(0L)
    ))
  }
  do.call(rbind, lapply(binary, function(menu_index) {
    items <- domain$items[[menu_index]]
    terminal_reach <- attention$reach[[menu_index]][2L]
    data.frame(
      menu_id = menu_index,
      earlier = items[1L],
      later = items[2L],
      reach = terminal_reach,
      share_later_preferred = if (terminal_reach > tolerance) {
        prob[menu_index, items[2L]] / terminal_reach
      } else {
        NA_real_
      },
      identified = terminal_reach > tolerance
    )
  }))
}

#' Population Analysis for Heterogeneous List-Based Attention Overload
#'
#' @description
#' `hlaoModel` recovers list-based reach probabilities and prefix masses on a
#' suffix-closed menu domain. It evaluates recovered-attention restrictions,
#' constructs sharp independent, dependence-robust, or path-independence-robust
#' preference polytopes, and computes sharp bounds for supplied preference
#' events. The first two modes require a suffix-closed domain because they use
#' Sequential Path Independence to recover attention. The `"noPI"` mode treats
#' prefix masses as latent and is available on any observed-menu domain. With
#' full menu data and positive terminal reach, the SPI modes also recover the
#' full-attention choice rule and report Block--Marschak diagnostics. Optional
#' agreement targets measure whether observed and full-attention choices agree.
#' Under benchmark independence, structured events support status-checked
#' column generation with mixed-integer pricing over linear orders. Returned
#' diagnostics report solver statuses, tolerance, reduced costs, primal and
#' dual residuals, an optimality-gap bound, and whether the numerical
#' certificate checks succeeded.
#'
#' @param menu Numeric matrix of zeros and ones with one row per distinct menu.
#' @param prob Numeric matrix of inside choice probabilities with the same
#'   dimensions as `menu`.
#' @param outside_prob Optional vector of outside-option probabilities. When
#'   omitted, it is computed as one minus the row sum of `prob`.
#' @param list_order Permutation giving the observed presentation order. The
#'   default is the column order of `menu`.
#' @param events Optional zero-one event indicators over the rows returned by
#'   [hlaoRankings()], or one or more structured [hlaoEvent()] objects.
#' @param dependence Which population polytope to construct: `"independent"`,
#'   `"robust"`, `"noPI"`, `"both"`, or `"all"`. For backward compatibility,
#'   `"both"` continues to request the independent and dependence-robust SPI
#'   polytopes; `"all"` adds the no-SPI polytope.
#' @param tolerance Nonnegative numerical tolerance for model diagnostics.
#' @param agreement `FALSE`, `TRUE`, or observed-menu indices. `TRUE` computes
#'   full-attention agreement bounds for every observed menu.
#' @param algorithm Computational method: `"auto"`, `"enumerate"`, or
#'   `"column_generation"`. Column generation currently applies to the
#'   benchmark independent model and structured events.
#' @param max_rankings Maximum number of ranking columns to enumerate.
#' @param max_iterations Maximum number of master and pricing iterations under
#'   column generation.
#'
#' @return An object of class `ramchoiceHLAOModel` containing recovered
#'   `attention`, attention `diagnostics`, population `pairwise` shares,
#'   compatibility by dependence mode, event `bounds`, full-attention
#'   `agreement`, ranking columns used by the selected algorithm, computation
#'   diagnostics, and, when available, `full_attention` and `block_marschak`
#'   results.
#'
#' @references
#' M. D. Cattaneo, P. H. Y. Cheung, X. Ma, and Y. Masatlioglu (2026).
#' \href{https://arxiv.org/abs/2110.10650}{Attention Overload}. Working paper.
#'
#' @examples
#' menu <- rbind(c(1, 0), c(0, 1), c(1, 1))
#' prob <- rbind(c(.8, 0), c(0, .75), c(.56, .24))
#' rankings <- hlaoRankings(1:2)
#' event <- rankings[, 1] == 2
#' hlaoModel(menu, prob, events = list(`2 above 1` = event))
#' hlaoModel(
#'   menu, prob,
#'   events = hlaoEvent(2, 1, name = "2 above 1"),
#'   agreement = TRUE
#' )
#'
#' @export
hlaoModel <- function(menu, prob, outside_prob = NULL, list_order = NULL,
                      events = NULL,
                      dependence = c(
                        "independent", "robust", "noPI", "both", "all"
                      ),
                      tolerance = sqrt(.Machine$double.eps),
                      agreement = FALSE,
                      algorithm = c("auto", "enumerate", "column_generation"),
                      max_rankings = 5000L,
                      max_iterations = 1000L) {
  started <- proc.time()[["elapsed"]]
  tolerance <- .hlao_validate_tolerance(tolerance)
  dependence <- match.arg(dependence)
  algorithm <- match.arg(algorithm)
  if (!is.numeric(max_rankings) || length(max_rankings) != 1L ||
      is.na(max_rankings) || max_rankings < 1) {
    stop("'max_rankings' must be one positive number.", call. = FALSE)
  }
  if (!is.numeric(max_iterations) || length(max_iterations) != 1L ||
      is.na(max_iterations) || max_iterations != as.integer(max_iterations) ||
      max_iterations < 1L) {
    stop("'max_iterations' must be one positive integer.", call. = FALSE)
  }
  max_iterations <- as.integer(max_iterations)
  prepared <- .hlao_prepare_population(menu, prob, outside_prob, tolerance)
  menu <- prepared$menu
  prob <- prepared$prob
  outside_prob <- prepared$outside_prob
  list_order <- .hlao_validate_list_order(list_order, ncol(menu))
  modes <- switch(
    dependence,
    both = c("independent", "robust"),
    all = c("independent", "robust", "noPI"),
    dependence
  )
  needs_spi <- any(modes != "noPI")
  domain <- if (needs_spi) {
    .hlao_domain(menu, list_order)
  } else {
    .hlao_domain_unrestricted(menu, list_order)
  }
  attention <- if (needs_spi) {
    .hlao_recover_attention(menu, outside_prob, domain, tolerance)
  } else {
    NULL
  }
  full_attention <- if (needs_spi) {
    .hlao_full_attention(menu, prob, domain, attention, tolerance)
  } else {
    NULL
  }
  block_marschak <- if (needs_spi) {
    .hlao_bm_diagnostics(menu, domain, full_attention, tolerance)
  } else {
    NULL
  }

  n_rankings <- gamma(ncol(menu) + 1L)
  structured_events <- .hlao_structured_events(events, ncol(menu))
  agreement_indices <- .hlao_agreement_indices(agreement, nrow(menu))
  column_agreement_indices <- if (is.null(full_attention)) {
    agreement_indices
  } else {
    integer(0L)
  }
  use_column_generation <- algorithm == "column_generation" ||
    algorithm == "auto" && (!is.finite(n_rankings) || n_rankings > max_rankings)
  if (use_column_generation && !identical(modes, "independent")) {
    stop(
      "Column generation currently supports dependence = 'independent'.",
      call. = FALSE
    )
  }
  if (use_column_generation && is.null(structured_events)) {
    stop(
      "Column generation requires events created by 'hlaoEvent()'.",
      call. = FALSE
    )
  }
  if (!use_column_generation &&
      (!is.finite(n_rankings) || n_rankings > max_rankings)) {
    stop(
      "The ranking polytope requires ", n_rankings,
      " columns, exceeding 'max_rankings'. Use structured events with column generation.",
      call. = FALSE
    )
  }

  fits <- list()
  compatibility <- vector("list", length(modes))
  bounds <- list()
  agreement_results <- NULL
  computation <- NULL

  if (use_column_generation) {
    targets <- .hlao_column_targets(
      structured_events, column_agreement_indices, domain
    )
    if (!attention$diagnostics$valid) {
      column_fit <- NULL
      rankings <- matrix(integer(0L), nrow = 0L, ncol = ncol(menu))
      event_matrix <- matrix(
        numeric(0L), nrow = 0L, ncol = length(structured_events)
      )
      fit <- list(
        feasible = FALSE, status = NA_integer_,
        bounds = data.frame(
          event = character(0L), mode = character(0L),
          lower = numeric(0L), upper = numeric(0L),
          lower_status = integer(0L), upper_status = integer(0L)
        )
      )
      system <- NULL
    } else {
      column_fit <- .hlao_column_generation(
        menu, prob, domain, attention, targets,
        tolerance = max(tolerance, 1e-9),
        max_iterations = max_iterations
      )
      rankings <- column_fit$rankings
      event_matrix <- .hlao_event_matrix(structured_events, rankings)
      target_bounds <- column_fit$bounds
      event_bounds <- target_bounds[target_bounds$kind == "event", , drop = FALSE]
      fit <- list(
        feasible = column_fit$feasible,
        status = column_fit$status,
        bounds = if (nrow(event_bounds)) {
          data.frame(
            event = event_bounds$target,
            mode = "independent",
            lower = event_bounds$lower,
            upper = event_bounds$upper,
            lower_status = event_bounds$lower_status,
            upper_status = event_bounds$upper_status,
            lower_pricing_status = event_bounds$lower_pricing_status,
            upper_pricing_status = event_bounds$upper_pricing_status,
            lower_reduced_cost = event_bounds$lower_reduced_cost,
            upper_reduced_cost = event_bounds$upper_reduced_cost,
            lower_dual_residual = event_bounds$lower_dual_residual,
            upper_dual_residual = event_bounds$upper_dual_residual,
            lower_optimality_gap_bound = event_bounds$lower_optimality_gap_bound,
            upper_optimality_gap_bound = event_bounds$upper_optimality_gap_bound,
            lower_primal_residual = event_bounds$lower_primal_residual,
            upper_primal_residual = event_bounds$upper_primal_residual,
            lower_certified = event_bounds$lower_certified,
            upper_certified = event_bounds$upper_certified,
            stringsAsFactors = FALSE
          )
        } else {
          data.frame(
            event = character(0L), mode = character(0L),
            lower = numeric(0L), upper = numeric(0L),
            lower_status = integer(0L), upper_status = integer(0L)
          )
        }
      )
      agreement_bounds <- target_bounds[
        target_bounds$kind == "agreement", , drop = FALSE
      ]
      if (nrow(agreement_bounds)) {
        agreement_results <- data.frame(
          menu_id = agreement_bounds$menu_id,
          menu = agreement_bounds$menu,
          mode = "independent",
          lower = agreement_bounds$lower,
          upper = agreement_bounds$upper,
          point_identified = abs(
            agreement_bounds$upper - agreement_bounds$lower
          ) <= 1e-8,
          method = "column-generation",
          lower_status = agreement_bounds$lower_status,
          upper_status = agreement_bounds$upper_status,
          lower_pricing_status = agreement_bounds$lower_pricing_status,
          upper_pricing_status = agreement_bounds$upper_pricing_status,
          lower_reduced_cost = agreement_bounds$lower_reduced_cost,
          upper_reduced_cost = agreement_bounds$upper_reduced_cost,
          lower_dual_residual = agreement_bounds$lower_dual_residual,
          upper_dual_residual = agreement_bounds$upper_dual_residual,
          lower_optimality_gap_bound = agreement_bounds$lower_optimality_gap_bound,
          upper_optimality_gap_bound = agreement_bounds$upper_optimality_gap_bound,
          lower_primal_residual = agreement_bounds$lower_primal_residual,
          upper_primal_residual = agreement_bounds$upper_primal_residual,
          lower_certified = agreement_bounds$lower_certified,
          upper_certified = agreement_bounds$upper_certified,
          stringsAsFactors = FALSE
        )
      }
      system <- list(
        A = column_fit$columns,
        b = column_fit$spec$b,
        n_variables = ncol(column_fit$columns)
      )
      computation <- column_fit$diagnostics
    }
    fits$independent <- list(fit = fit, system = system)
    compatibility[[1L]] <- data.frame(
      mode = "independent",
      attention_required = TRUE,
      attention_valid = attention$diagnostics$valid,
      preference_polytope_nonempty = fit$feasible,
      compatible = attention$diagnostics$valid && fit$feasible,
      lp_status = fit$status,
      n_rankings = nrow(rankings),
      n_variables = if (is.null(system)) NA_integer_ else ncol(system$A),
      n_constraints = if (is.null(system)) NA_integer_ else nrow(system$A)
    )
    if (nrow(fit$bounds)) {
      bounds[[1L]] <- fit$bounds
    }
  } else {
    rankings <- hlaoRankings(seq_len(ncol(menu)))
    event_matrix <- .hlao_normalize_events(events, rankings)
    for (mode_index in seq_along(modes)) {
      mode <- modes[mode_index]
      if (mode != "noPI" && !attention$diagnostics$valid) {
        fit <- list(
          feasible = FALSE,
          status = NA_integer_,
          bounds = data.frame(
            event = character(0L), mode = character(0L),
            lower = numeric(0L), upper = numeric(0L),
            lower_status = integer(0L), upper_status = integer(0L)
          )
        )
        system <- NULL
      } else if (mode == "independent") {
        system <- .hlao_independent_system(menu, prob, domain, attention, rankings)
        fit <- .hlao_solve_bounds(system, event_matrix, mode)
      } else if (mode == "robust") {
        system <- .hlao_population_robust_system(
          menu, prob, domain, attention, rankings
        )
        fit <- .hlao_solve_bounds(system, event_matrix, mode)
      } else {
        system <- .hlao_nopi_system(
          menu, domain, rankings,
          prob = prob,
          outside_prob = outside_prob
        )
        fit <- .hlao_solve_general_bounds(system, event_matrix, mode)
      }
      fits[[mode]] <- list(fit = fit, system = system)
      compatibility[[mode_index]] <- data.frame(
        mode = mode,
        attention_required = mode != "noPI",
        attention_valid = if (mode == "noPI") NA else attention$diagnostics$valid,
        preference_polytope_nonempty = fit$feasible,
        compatible = (mode == "noPI" || attention$diagnostics$valid) && fit$feasible,
        lp_status = fit$status,
        n_rankings = nrow(rankings),
        n_variables = if (is.null(system)) NA_integer_ else ncol(system$A),
        n_constraints = if (is.null(system)) NA_integer_ else nrow(system$A)
      )
      if (nrow(fit$bounds)) {
        bounds[[length(bounds) + 1L]] <- fit$bounds
      }
    }
    computation <- data.frame(
      algorithm = "enumeration",
      n_columns = nrow(rankings),
      n_possible_rankings = n_rankings,
      iterations = NA_integer_,
      pricing_binary_variables = NA_integer_,
      pricing_constraints = NA_integer_,
      stringsAsFactors = FALSE
    )
  }

  if (use_column_generation && !is.null(full_attention)) {
    agreement_results <- .hlao_direct_agreement(
      full_attention, attention, domain, agreement_indices
    )
  }
  if (!use_column_generation) {
    agreement_results <- .hlao_enumerated_agreement(
      fits, modes, rankings, domain, attention, agreement_indices,
      full_attention = full_attention
    )
  }
  empty_attention <- data.frame(
    menu_id = integer(0L), menu = character(0L),
    alternative = integer(0L), position = integer(0L),
    reach = numeric(0L), stop_mass = numeric(0L)
  )
  empty_attention_diagnostics <- data.frame(
    valid = NA, min_prefix_mass = NA_real_,
    n_negative_prefix_masses = NA_integer_,
    max_attention_overload_violation = NA_real_,
    n_attention_overload_violations = NA_integer_
  )
  result <- list(
    compatibility = do.call(rbind, compatibility),
    attention = if (needs_spi) attention$attention else empty_attention,
    attention_diagnostics = if (needs_spi) {
      attention$diagnostics
    } else {
      empty_attention_diagnostics
    },
    attention_overload = if (needs_spi) attention$overload else NULL,
    prefix_masses = if (needs_spi) attention$masses else NULL,
    pairwise = if (needs_spi) {
      .hlao_population_pairwise(menu, prob, domain, attention, tolerance)
    } else {
      NULL
    },
    bounds = if (length(bounds)) do.call(rbind, bounds) else {
      data.frame(
        event = character(0L), mode = character(0L),
        lower = numeric(0L), upper = numeric(0L),
        lower_status = integer(0L), upper_status = integer(0L)
      )
    },
    agreement = agreement_results,
    rankings = rankings,
    events = event_matrix,
    full_attention = full_attention,
    block_marschak = block_marschak,
    menu = menu,
    prob = prob,
    outside_prob = outside_prob,
    list_order = list_order,
    fits = fits,
    computation = computation,
    tolerance = tolerance,
    elapsed = unname(proc.time()[["elapsed"]] - started)
  )
  class(result) <- "ramchoiceHLAOModel"
  result
}

.hlao_sample_summary <- function(menu, choice, outside, cluster = NULL) {
  menu <- .hlao_validate_menu(menu, allow_duplicates = TRUE)
  if (!is.matrix(choice) || !is.numeric(choice) && !is.logical(choice) ||
      !identical(dim(choice), dim(menu)) || anyNA(choice) ||
      any(!(choice %in% c(0, 1)))) {
    stop("'choice' must be a zero-one matrix with the same dimensions as 'menu'.", call. = FALSE)
  }
  choice <- matrix(as.integer(choice), nrow = nrow(choice), ncol = ncol(choice))
  if (any(choice[menu == 0L] != 0L) || any(rowSums(choice) > 1L)) {
    stop("Each observation must choose at most one available inside alternative.", call. = FALSE)
  }
  if (is.null(outside)) {
    outside <- as.integer(rowSums(choice) == 0L)
  }
  if (!is.numeric(outside) && !is.logical(outside) || length(outside) != nrow(menu) ||
      anyNA(outside) || any(!(outside %in% c(0, 1))) ||
      any(rowSums(choice) + outside != 1L)) {
    stop("'outside' must identify the outside choice for every observation.", call. = FALSE)
  }

  keys <- apply(menu, 1L, .hlao_menu_key)
  unique_keys <- unique(keys)
  groups <- match(keys, unique_keys)
  first <- match(unique_keys, keys)
  sum_menu <- menu[first, , drop = FALSE]
  sample_size <- tabulate(groups, nbins = length(unique_keys))
  sum_choice <- matrix(0, nrow = length(unique_keys), ncol = ncol(menu))
  sum_outside <- numeric(length(unique_keys))
  for (group in seq_along(unique_keys)) {
    selected <- groups == group
    sum_choice[group, ] <- colSums(choice[selected, , drop = FALSE])
    sum_outside[group] <- sum(outside[selected])
  }
  cluster_info <- .ram_validate_cluster(cluster, nrow(menu))
  n_cell <- sum(rowSums(sum_menu) + 1L)
  covariance <- matrix(0, nrow = n_cell, ncol = n_cell)
  outside_cell <- integer(nrow(sum_menu))
  inside_cell <- matrix(
    0L,
    nrow = nrow(sum_menu),
    ncol = ncol(sum_menu)
  )
  cluster_scores <- cluster_cell_count <- cluster_menu_count <- NULL
  if (!is.null(cluster_info)) {
    cluster_scores <- matrix(0, nrow = cluster_info$n, ncol = n_cell)
    cluster_cell_count <- matrix(0, nrow = cluster_info$n, ncol = n_cell)
    cluster_menu_count <- matrix(
      0,
      nrow = cluster_info$n,
      ncol = nrow(sum_menu)
    )
  }

  next_cell <- 0L
  for (menu_index in seq_len(nrow(sum_menu))) {
    items <- which(sum_menu[menu_index, ] == 1L)
    current <- next_cell + seq_len(length(items) + 1L)
    next_cell <- max(current)
    outside_cell[menu_index] <- current[1L]
    inside_cell[menu_index, items] <- current[-1L]
    probabilities <- c(
      sum_outside[menu_index] / sample_size[menu_index],
      sum_choice[menu_index, items] / sample_size[menu_index]
    )
    if (is.null(cluster_info)) {
      covariance[current, current] <- (
        diag(probabilities) - tcrossprod(probabilities)
      ) / sample_size[menu_index]
    } else {
      selected_menu <- groups == menu_index
      outcomes <- cbind(outside, choice[, items, drop = FALSE])
      for (cluster_index in seq_len(cluster_info$n)) {
        selected <- selected_menu & cluster_info$id == cluster_index
        cluster_menu_size <- sum(selected)
        cluster_menu_count[cluster_index, menu_index] <- cluster_menu_size
        if (cluster_menu_size) {
          counts <- colSums(outcomes[selected, , drop = FALSE])
          cluster_cell_count[cluster_index, current] <- counts
          cluster_scores[cluster_index, current] <- (
            counts - cluster_menu_size * probabilities
          ) / sample_size[menu_index]
        }
      }
    }
  }
  if (!is.null(cluster_info)) {
    covariance <- .ram_cluster_covariance(cluster_scores)
  }

  list(
    menu = sum_menu,
    prob = sum_choice / sample_size,
    outside_prob = sum_outside / sample_size,
    sample_size = sample_size,
    covariance = covariance,
    outside_cell = outside_cell,
    inside_cell = inside_cell,
    sampling = if (is.null(cluster_info)) "iid" else "cluster",
    n_cluster = if (is.null(cluster_info)) NA_integer_ else cluster_info$n,
    cluster_labels = if (is.null(cluster_info)) NULL else cluster_info$labels,
    cluster_scores = cluster_scores,
    cluster_cell_count = cluster_cell_count,
    cluster_menu_count = cluster_menu_count
  )
}

.hlao_ensure_inference_summary <- function(summary) {
  if (!is.null(summary$covariance) && !is.null(summary$sampling)) {
    return(summary)
  }
  n_cell <- sum(rowSums(summary$menu) + 1L)
  summary$covariance <- matrix(0, nrow = n_cell, ncol = n_cell)
  summary$outside_cell <- integer(nrow(summary$menu))
  summary$inside_cell <- matrix(
    0L,
    nrow = nrow(summary$menu),
    ncol = ncol(summary$menu)
  )
  next_cell <- 0L
  for (menu_index in seq_len(nrow(summary$menu))) {
    items <- which(summary$menu[menu_index, ] == 1L)
    current <- next_cell + seq_len(length(items) + 1L)
    next_cell <- max(current)
    summary$outside_cell[menu_index] <- current[1L]
    summary$inside_cell[menu_index, items] <- current[-1L]
    probabilities <- c(
      summary$outside_prob[menu_index],
      summary$prob[menu_index, items]
    )
    summary$covariance[current, current] <- (
      diag(probabilities) - tcrossprod(probabilities)
    ) / summary$sample_size[menu_index]
  }
  summary$sampling <- "iid"
  summary$n_cluster <- NA_integer_
  summary$cluster_labels <- NULL
  summary$cluster_scores <- NULL
  summary$cluster_cell_count <- NULL
  summary$cluster_menu_count <- NULL
  summary
}

.hlao_probability_bands <- function(summary, alpha,
                                    band_method = c("hoeffding", "gaussian"),
                                    n_band_draws = 2000L,
                                    boundary_count = 5L) {
  summary <- .hlao_ensure_inference_summary(summary)
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be one number strictly between zero and one.", call. = FALSE)
  }
  band_method <- match.arg(band_method)
  if (!is.numeric(n_band_draws) || length(n_band_draws) != 1L ||
      is.na(n_band_draws) || n_band_draws < 100L ||
      n_band_draws != as.integer(n_band_draws)) {
    stop("'n_band_draws' must be an integer of at least 100.", call. = FALSE)
  }
  if (!is.numeric(boundary_count) || length(boundary_count) != 1L ||
      is.na(boundary_count) || boundary_count < 0L ||
      boundary_count != as.integer(boundary_count)) {
    stop("'boundary_count' must be a nonnegative integer.", call. = FALSE)
  }
  n_band_draws <- as.integer(n_band_draws)
  boundary_count <- as.integer(boundary_count)
  n_cells <- sum(rowSums(summary$menu) + 1L)
  hoeffding_width <- function(error_probability) {
    if (summary$sampling == "iid") {
      return(sqrt(
        log(2 * n_cells / error_probability) /
          (2 * summary$sample_size)
      ))
    }
    cluster_weights <- sweep(
      summary$cluster_menu_count,
      2L,
      summary$sample_size,
      "/"
    )
    sqrt(
      log(2 * n_cells / error_probability) / 2 *
        colSums(cluster_weights^2)
    )
  }
  epsilon <- hoeffding_width(alpha)
  cell_standard_error <- sqrt(pmax(0, diag(summary$covariance)))
  standard_error <- matrix(
    0,
    nrow = nrow(summary$prob),
    ncol = ncol(summary$prob)
  )
  outside_standard_error <- cell_standard_error[summary$outside_cell]
  for (menu_index in seq_len(nrow(summary$menu))) {
    items <- which(summary$menu[menu_index, ] == 1L)
    standard_error[menu_index, items] <- cell_standard_error[
      summary$inside_cell[menu_index, items]
    ]
  }
  half_width <- matrix(
    rep(epsilon, ncol(summary$menu)),
    nrow = nrow(summary$menu),
    ncol = ncol(summary$menu)
  )
  half_width[summary$menu == 0L] <- 0
  outside_half_width <- epsilon
  gaussian_active <- matrix(
    FALSE,
    nrow = nrow(summary$menu),
    ncol = ncol(summary$menu)
  )
  outside_gaussian_active <- rep(FALSE, nrow(summary$menu))
  critical_value <- NA_real_
  gaussian_alpha <- NA_real_
  fallback_alpha <- alpha
  n_gaussian_cells <- 0L
  n_fallback_cells <- n_cells

  if (band_method == "gaussian") {
    active_by_menu <- vector("list", nrow(summary$menu))
    active_cells <- rep(FALSE, n_cells)
    for (menu_index in seq_len(nrow(summary$menu))) {
      items <- which(summary$menu[menu_index, ] == 1L)
      cells <- c(
        summary$outside_cell[menu_index],
        summary$inside_cell[menu_index, items]
      )
      if (summary$sampling == "iid") {
        probabilities <- c(
          summary$outside_prob[menu_index],
          summary$prob[menu_index, items]
        )
        counts <- round(probabilities * summary$sample_size[menu_index])
        active <- counts >= boundary_count &
          summary$sample_size[menu_index] - counts >= boundary_count
      } else {
        cell_counts <- summary$cluster_cell_count[, cells, drop = FALSE]
        menu_counts <- summary$cluster_menu_count[, menu_index]
        successes <- colSums(cell_counts > 0)
        failures <- colSums(sweep(cell_counts, 1L, menu_counts, "-") < 0)
        active <- successes >= boundary_count & failures >= boundary_count
      }
      active <- active &
        cell_standard_error[cells] > sqrt(.Machine$double.eps)
      active_by_menu[[menu_index]] <- active
      active_cells[cells] <- active
    }
    n_gaussian_cells <- sum(active_cells)
    n_fallback_cells <- n_cells - n_gaussian_cells
    if (n_gaussian_cells > 0L && n_fallback_cells > 0L) {
      gaussian_alpha <- alpha / 2
      fallback_alpha <- alpha / 2
      epsilon <- hoeffding_width(fallback_alpha)
      half_width <- matrix(
        rep(epsilon, ncol(summary$menu)),
        nrow = nrow(summary$menu),
        ncol = ncol(summary$menu)
      )
      half_width[summary$menu == 0L] <- 0
      outside_half_width <- epsilon
    } else if (n_gaussian_cells > 0L) {
      gaussian_alpha <- alpha
      fallback_alpha <- NA_real_
    }
    if (n_gaussian_cells > 0L) {
      draws <- if (summary$sampling == "cluster") {
        .ram_cluster_multiplier_draws(
          summary$cluster_scores,
          n_band_draws
        )
      } else {
        MASS::mvrnorm(
        n = n_band_draws,
          mu = rep(0, n_cells),
          Sigma = summary$covariance,
        tol = 1e-10
      )
      }
      draws <- matrix(
        draws,
        nrow = n_band_draws,
        ncol = n_cells
      )
      standardized <- abs(sweep(
        draws[, active_cells, drop = FALSE],
        2L,
        cell_standard_error[active_cells],
        "/"
      ))
      critical_value <- as.numeric(stats::quantile(
        apply(standardized, 1L, max),
        probs = 1 - gaussian_alpha,
        names = FALSE,
        type = 8
      ))
      for (menu_index in seq_len(nrow(summary$menu))) {
        items <- which(summary$menu[menu_index, ] == 1L)
        active <- active_by_menu[[menu_index]]
        outside_gaussian_active[menu_index] <- active[1L]
        gaussian_active[menu_index, items] <- active[-1L]
        if (active[1L]) {
          outside_half_width[menu_index] <-
            critical_value * outside_standard_error[menu_index]
        }
        selected_items <- items[active[-1L]]
        if (length(selected_items)) {
          half_width[menu_index, selected_items] <-
            critical_value * standard_error[menu_index, selected_items]
        }
      }
    }
  }

  lower <- matrix(
    pmax(0, as.numeric(summary$prob - half_width)),
    nrow = nrow(summary$prob),
    ncol = ncol(summary$prob)
  )
  upper <- matrix(
    pmin(1, as.numeric(summary$prob + half_width)),
    nrow = nrow(summary$prob),
    ncol = ncol(summary$prob)
  )
  lower[summary$menu == 0L] <- 0
  upper[summary$menu == 0L] <- 0
  outside_lower <- pmax(0, summary$outside_prob - outside_half_width)
  outside_upper <- pmin(1, summary$outside_prob + outside_half_width)
  if (band_method == "gaussian" && n_fallback_cells > 0L &&
      summary$sampling == "iid") {
    cell_alpha <- fallback_alpha / n_cells
    exact_interval <- function(count, sample_size) {
      c(
        lower = if (count == 0L) {
          0
        } else {
          stats::qbeta(cell_alpha / 2, count, sample_size - count + 1L)
        },
        upper = if (count == sample_size) {
          1
        } else {
          stats::qbeta(
            1 - cell_alpha / 2,
            count + 1L,
            sample_size - count
          )
        }
      )
    }
    for (menu_index in seq_len(nrow(summary$menu))) {
      sample_size <- summary$sample_size[menu_index]
      if (!outside_gaussian_active[menu_index]) {
        interval <- exact_interval(
          round(summary$outside_prob[menu_index] * sample_size),
          sample_size
        )
        outside_lower[menu_index] <- interval["lower"]
        outside_upper[menu_index] <- interval["upper"]
        outside_half_width[menu_index] <- max(
          summary$outside_prob[menu_index] - interval["lower"],
          interval["upper"] - summary$outside_prob[menu_index]
        )
      }
      for (alternative in which(
        summary$menu[menu_index, ] == 1L &
          !gaussian_active[menu_index, ]
      )) {
        interval <- exact_interval(
          round(summary$prob[menu_index, alternative] * sample_size),
          sample_size
        )
        lower[menu_index, alternative] <- interval["lower"]
        upper[menu_index, alternative] <- interval["upper"]
        half_width[menu_index, alternative] <- max(
          summary$prob[menu_index, alternative] - interval["lower"],
          interval["upper"] - summary$prob[menu_index, alternative]
        )
      }
    }
  }
  method_label <- if (band_method == "hoeffding") {
    if (summary$sampling == "cluster") "cluster-hoeffding" else "hoeffding"
  } else {
    if (summary$sampling == "cluster") {
      "cluster-multiplier-hoeffding-hybrid"
    } else {
      "correlated-gaussian-exact-hybrid"
    }
  }

  rows <- list()
  row_index <- 0L
  for (menu_index in seq_len(nrow(summary$menu))) {
    row_index <- row_index + 1L
    rows[[row_index]] <- data.frame(
      menu_id = menu_index,
      outcome = 0L,
      estimate = summary$outside_prob[menu_index],
      lower = outside_lower[menu_index],
      upper = outside_upper[menu_index],
      standard_error = outside_standard_error[menu_index],
      half_width = outside_half_width[menu_index],
      epsilon = outside_half_width[menu_index],
      gaussian_active = outside_gaussian_active[menu_index],
      sample_size = summary$sample_size[menu_index],
      band_method = method_label
    )
    for (alternative in which(summary$menu[menu_index, ] == 1L)) {
      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        menu_id = menu_index,
        outcome = alternative,
        estimate = summary$prob[menu_index, alternative],
        lower = lower[menu_index, alternative],
        upper = upper[menu_index, alternative],
        standard_error = standard_error[menu_index, alternative],
        half_width = half_width[menu_index, alternative],
        epsilon = half_width[menu_index, alternative],
        gaussian_active = gaussian_active[menu_index, alternative],
        sample_size = summary$sample_size[menu_index],
        band_method = method_label
      )
    }
  }
  list(
    results = do.call(rbind, rows),
    lower = lower,
    upper = upper,
    outside_lower = outside_lower,
    outside_upper = outside_upper,
    epsilon = epsilon,
    n_cells = n_cells,
    method = method_label,
    critical_value = critical_value,
    gaussian_alpha = gaussian_alpha,
    fallback_alpha = fallback_alpha,
    n_gaussian_cells = n_gaussian_cells,
    n_fallback_cells = n_fallback_cells
  )
}

.hlao_reach_bounds <- function(domain, bands) {
  lower <- upper <- vector("list", length(domain$items))
  for (menu_index in seq_along(domain$items)) {
    suffix <- domain$suffix_index[[menu_index]]
    lower[[menu_index]] <- cumprod(1 - bands$outside_upper[suffix])
    upper[[menu_index]] <- cumprod(1 - bands$outside_lower[suffix])
  }
  list(lower = lower, upper = upper)
}

.hlao_specification_row <- function(restriction, restriction_id, estimate,
                                    lower, upper, available = TRUE,
                                    tolerance = 0) {
  data.frame(
    restriction = restriction,
    restriction_id = restriction_id,
    estimate = estimate,
    lower = lower,
    upper = upper,
    available = available,
    reject = available && is.finite(lower) && lower > tolerance,
    stringsAsFactors = FALSE
  )
}

.hlao_attention_specification <- function(domain, attention, reach_bounds,
                                          tolerance) {
  if (!nrow(attention$overload)) {
    return(data.frame(
      restriction = character(0L), restriction_id = character(0L),
      estimate = numeric(0L), lower = numeric(0L), upper = numeric(0L),
      available = logical(0L), reject = logical(0L)
    ))
  }
  rows <- vector("list", nrow(attention$overload))
  for (index in seq_len(nrow(attention$overload))) {
    current <- attention$overload[index, ]
    small_position <- match(
      current$alternative,
      domain$items[[current$smaller_menu_id]]
    )
    large_position <- match(
      current$alternative,
      domain$items[[current$larger_menu_id]]
    )
    lower <-
      reach_bounds$lower[[current$larger_menu_id]][large_position] -
      reach_bounds$upper[[current$smaller_menu_id]][small_position]
    upper <-
      reach_bounds$upper[[current$larger_menu_id]][large_position] -
      reach_bounds$lower[[current$smaller_menu_id]][small_position]
    rows[[index]] <- .hlao_specification_row(
      restriction = "attention-overload",
      restriction_id = paste0(
        "AO-M", current$smaller_menu_id, "-M", current$larger_menu_id,
        "-A", current$alternative
      ),
      estimate = current$violation,
      lower = lower,
      upper = upper,
      tolerance = tolerance
    )
  }
  do.call(rbind, rows)
}

.hlao_full_attention_bands <- function(summary, domain, bands, reach_bounds,
                                       full_attention, tolerance) {
  if (!.hlao_is_full_domain(summary$menu)) {
    return(list(
      available = FALSE,
      reason = "Full-attention diagnostics require the complete menu domain.",
      lower = NULL,
      upper = NULL,
      probability_diagnostics = data.frame(
        restriction = character(0L), restriction_id = character(0L),
        estimate = numeric(0L), lower = numeric(0L), upper = numeric(0L),
        available = logical(0L), reject = logical(0L)
      )
    ))
  }

  lower <- upper <- matrix(
    0,
    nrow = nrow(summary$menu),
    ncol = ncol(summary$menu)
  )
  probability_rows <- list()
  probability_index <- 0L
  menu_order <- order(rowSums(summary$menu))

  add_probability_diagnostic <- function(menu_index, alternative,
                                         raw_lower, raw_upper) {
    point <- if (is.null(full_attention)) {
      NA_real_
    } else {
      full_attention$matrix[menu_index, alternative]
    }
    probability_index <<- probability_index + 1L
    probability_rows[[probability_index]] <<- .hlao_specification_row(
      restriction = "full-attention-probability",
      restriction_id = paste0("F-M", menu_index, "-A", alternative),
      estimate = if (is.finite(point)) max(-point, point - 1) else NA_real_,
      lower = max(-raw_upper, raw_lower - 1),
      upper = max(-raw_lower, raw_upper - 1),
      tolerance = tolerance
    )
  }

  for (menu_index in menu_order) {
    items <- domain$items[[menu_index]]
    n_items <- length(items)
    if (n_items == 1L) {
      lower[menu_index, items] <- 1
      upper[menu_index, items] <- 1
      add_probability_diagnostic(menu_index, items, 1, 1)
      next
    }

    terminal_lower <- reach_bounds$lower[[menu_index]][n_items]
    terminal_upper <- reach_bounds$upper[[menu_index]][n_items]
    for (position in 2:n_items) {
      alternative <- items[position]
      if (terminal_lower <= tolerance) {
        raw_lower <- 0
        raw_upper <- 1
      } else {
        adjustment_lower <- adjustment_upper <- 0
        if (position <= n_items - 1L) {
          for (prefix_size in position:(n_items - 1L)) {
            prefix <- integer(ncol(summary$menu))
            prefix[items[seq_len(prefix_size)]] <- 1L
            prefix_index <- unname(
              domain$lookup[.hlao_menu_key(prefix)]
            )
            mass_lower <- max(
              0,
              reach_bounds$lower[[menu_index]][prefix_size] -
                reach_bounds$upper[[menu_index]][prefix_size + 1L]
            )
            mass_upper <- min(
              1,
              max(
                0,
                reach_bounds$upper[[menu_index]][prefix_size] -
                  reach_bounds$lower[[menu_index]][prefix_size + 1L]
              )
            )
            adjustment_lower <- adjustment_lower +
              mass_lower * lower[prefix_index, alternative]
            adjustment_upper <- adjustment_upper +
              mass_upper * upper[prefix_index, alternative]
          }
        }
        numerator_lower <-
          bands$lower[menu_index, alternative] - adjustment_upper
        numerator_upper <-
          bands$upper[menu_index, alternative] - adjustment_lower
        ratios <- c(
          numerator_lower / terminal_lower,
          numerator_lower / terminal_upper,
          numerator_upper / terminal_lower,
          numerator_upper / terminal_upper
        )
        raw_lower <- min(ratios)
        raw_upper <- max(ratios)
      }
      add_probability_diagnostic(
        menu_index,
        alternative,
        raw_lower,
        raw_upper
      )
      if (raw_lower <= 1 && raw_upper >= 0) {
        lower[menu_index, alternative] <- max(0, raw_lower)
        upper[menu_index, alternative] <- min(1, raw_upper)
      } else {
        lower[menu_index, alternative] <- 0
        upper[menu_index, alternative] <- 1
      }
    }

    first <- items[1L]
    raw_lower <- 1 - sum(upper[menu_index, items[-1L]])
    raw_upper <- 1 - sum(lower[menu_index, items[-1L]])
    add_probability_diagnostic(
      menu_index,
      first,
      raw_lower,
      raw_upper
    )
    if (raw_lower <= 1 && raw_upper >= 0) {
      lower[menu_index, first] <- max(0, raw_lower)
      upper[menu_index, first] <- min(1, raw_upper)
    } else {
      lower[menu_index, first] <- 0
      upper[menu_index, first] <- 1
    }
  }

  list(
    available = TRUE,
    reason = NA_character_,
    lower = lower,
    upper = upper,
    probability_diagnostics = do.call(rbind, probability_rows)
  )
}

.hlao_bm_specification <- function(summary, domain, full_attention,
                                   full_attention_bands, tolerance) {
  if (!full_attention_bands$available) {
    return(data.frame(
      restriction = character(0L), restriction_id = character(0L),
      estimate = numeric(0L), lower = numeric(0L), upper = numeric(0L),
      available = logical(0L), reject = logical(0L)
    ))
  }
  rows <- list()
  row_index <- 0L
  for (target in seq_len(nrow(summary$menu))) {
    target_items <- domain$items[[target]]
    supersets <- which(apply(summary$menu, 1L, function(row) {
      all(row >= summary$menu[target, ])
    }))
    signs <- (-1)^(rowSums(summary$menu[supersets, , drop = FALSE]) -
                    sum(summary$menu[target, ]))
    for (alternative in target_items) {
      bm_lower <- sum(ifelse(
        signs > 0,
        full_attention_bands$lower[supersets, alternative],
        -full_attention_bands$upper[supersets, alternative]
      ))
      bm_upper <- sum(ifelse(
        signs > 0,
        full_attention_bands$upper[supersets, alternative],
        -full_attention_bands$lower[supersets, alternative]
      ))
      estimate <- if (is.null(full_attention)) {
        NA_real_
      } else {
        sum(signs * full_attention$matrix[supersets, alternative])
      }
      row_index <- row_index + 1L
      rows[[row_index]] <- .hlao_specification_row(
        restriction = "block-marschak",
        restriction_id = paste0("BM-M", target, "-A", alternative),
        estimate = if (is.finite(estimate)) -estimate else NA_real_,
        lower = -bm_upper,
        upper = -bm_lower,
        tolerance = tolerance
      )
    }
  }
  do.call(rbind, rows)
}

.hlao_delta_moments <- function(summary, domain, attention, full_attention,
                                tolerance) {
  if (!.hlao_is_full_domain(summary$menu) || is.null(full_attention)) {
    stop(
      "Direct delta diagnostics require a complete menu domain and positive terminal reach.",
      call. = FALSE
    )
  }

  n_menu <- nrow(summary$menu)
  n_alternative <- ncol(summary$menu)
  n_cell <- sum(rowSums(summary$menu) + 1L)
  covariance <- summary$covariance
  outside_gradient <- matrix(0, n_menu, n_cell)
  probability_gradient <- array(
    0,
    dim = c(n_menu, n_alternative, n_cell)
  )
  next_cell <- 0L
  for (menu_index in seq_len(n_menu)) {
    items <- domain$items[[menu_index]]
    current <- next_cell + seq_len(length(items) + 1L)
    next_cell <- max(current)
    outside_gradient[menu_index, current[1L]] <- 1
    for (position in seq_along(items)) {
      probability_gradient[
        menu_index,
        items[position],
        current[position + 1L]
      ] <- 1
    }
  }

  reach_gradient <- vector("list", n_menu)
  for (menu_index in seq_len(n_menu)) {
    suffix <- domain$suffix_index[[menu_index]]
    reach_gradient[[menu_index]] <- matrix(0, length(suffix), n_cell)
    previous <- 1
    previous_gradient <- numeric(n_cell)
    for (position in seq_along(suffix)) {
      outside_index <- suffix[position]
      continuation <- 1 - summary$outside_prob[outside_index]
      current_gradient <- previous_gradient * continuation -
        previous * outside_gradient[outside_index, ]
      reach_gradient[[menu_index]][position, ] <- current_gradient
      previous <- attention$reach[[menu_index]][position]
      previous_gradient <- current_gradient
    }
  }

  full <- full_attention$matrix
  full_gradient <- array(0, dim = c(n_menu, n_alternative, n_cell))
  for (menu_index in order(rowSums(summary$menu))) {
    items <- domain$items[[menu_index]]
    n_items <- length(items)
    if (n_items == 1L) {
      next
    }
    for (position in 2:n_items) {
      alternative <- items[position]
      adjustment <- 0
      adjustment_gradient <- numeric(n_cell)
      if (position <= n_items - 1L) {
        for (prefix_size in position:(n_items - 1L)) {
          prefix <- integer(n_alternative)
          prefix[items[seq_len(prefix_size)]] <- 1L
          prefix_index <- unname(domain$lookup[.hlao_menu_key(prefix)])
          mass <- attention$reach[[menu_index]][prefix_size] -
            attention$reach[[menu_index]][prefix_size + 1L]
          mass_gradient <- reach_gradient[[menu_index]][prefix_size, ] -
            reach_gradient[[menu_index]][prefix_size + 1L, ]
          adjustment <- adjustment + mass * full[prefix_index, alternative]
          adjustment_gradient <- adjustment_gradient +
            mass_gradient * full[prefix_index, alternative] +
            mass * full_gradient[prefix_index, alternative, ]
        }
      }
      numerator <- summary$prob[menu_index, alternative] - adjustment
      numerator_gradient <- probability_gradient[menu_index, alternative, ] -
        adjustment_gradient
      denominator <- attention$reach[[menu_index]][n_items]
      denominator_gradient <- reach_gradient[[menu_index]][n_items, ]
      full_gradient[menu_index, alternative, ] <-
        numerator_gradient / denominator -
        numerator * denominator_gradient / denominator^2
    }
    first <- items[1L]
    full_gradient[menu_index, first, ] <- -apply(
      full_gradient[menu_index, items[-1L], , drop = FALSE],
      3L,
      sum
    )
  }

  estimates <- numeric(0L)
  gradients <- list()
  restrictions <- restriction_ids <- character(0L)
  add_restriction <- function(restriction, restriction_id, estimate,
                              gradient) {
    restrictions <<- c(restrictions, restriction)
    restriction_ids <<- c(restriction_ids, restriction_id)
    estimates <<- c(estimates, estimate)
    gradients[[length(gradients) + 1L]] <<- gradient
  }

  for (row_index in seq_len(nrow(attention$overload))) {
    current <- attention$overload[row_index, ]
    small_position <- match(
      current$alternative,
      domain$items[[current$smaller_menu_id]]
    )
    large_position <- match(
      current$alternative,
      domain$items[[current$larger_menu_id]]
    )
    add_restriction(
      "attention-overload",
      paste0(
        "AO-M", current$smaller_menu_id, "-M", current$larger_menu_id,
        "-A", current$alternative
      ),
      current$violation,
      reach_gradient[[current$larger_menu_id]][large_position, ] -
        reach_gradient[[current$smaller_menu_id]][small_position, ]
    )
  }

  for (menu_index in seq_len(n_menu)) {
    for (alternative in domain$items[[menu_index]]) {
      add_restriction(
        "full-attention-probability",
        paste0("F-L-M", menu_index, "-A", alternative),
        -full[menu_index, alternative],
        -full_gradient[menu_index, alternative, ]
      )
      add_restriction(
        "full-attention-probability",
        paste0("F-U-M", menu_index, "-A", alternative),
        full[menu_index, alternative] - 1,
        full_gradient[menu_index, alternative, ]
      )
    }
  }

  for (target in seq_len(n_menu)) {
    supersets <- which(apply(summary$menu, 1L, function(row) {
      all(row >= summary$menu[target, ])
    }))
    signs <- (-1)^(rowSums(summary$menu[supersets, , drop = FALSE]) -
                    sum(summary$menu[target, ]))
    for (alternative in domain$items[[target]]) {
      add_restriction(
        "block-marschak",
        paste0("BM-M", target, "-A", alternative),
        -sum(signs * full[supersets, alternative]),
        -apply(
          signs * full_gradient[supersets, alternative, , drop = FALSE],
          3L,
          sum
        )
      )
    }
  }

  list(
    restriction = restrictions,
    restriction_id = restriction_ids,
    estimate = estimates,
    gradient = do.call(rbind, gradients),
    covariance = covariance,
    sampling = summary$sampling,
    cluster_scores = summary$cluster_scores,
    tolerance = tolerance
  )
}

.hlao_delta_specification <- function(summary, domain, attention,
                                      full_attention, alpha, n_draws,
                                      tolerance) {
  moments <- .hlao_delta_moments(
    summary,
    domain,
    attention,
    full_attention,
    tolerance
  )
  standard_error <- sqrt(pmax(
    0,
    rowSums((moments$gradient %*% moments$covariance) * moments$gradient)
  ))
  primitive_draws <- if (moments$sampling == "cluster") {
    .ram_cluster_multiplier_draws(moments$cluster_scores, n_draws)
  } else {
    MASS::mvrnorm(
      n = n_draws,
      mu = numeric(ncol(moments$gradient)),
      Sigma = moments$covariance,
      tol = 1e-10
    )
  }
  functional_draws <- primitive_draws %*% t(moments$gradient)
  active <- standard_error > tolerance
  standardized <- matrix(0, n_draws, length(moments$estimate))
  standardized[, active] <- sweep(
    functional_draws[, active, drop = FALSE],
    2L,
    standard_error[active],
    "/"
  )
  maximum_draws <- apply(standardized, 1L, max)
  critical_value <- as.numeric(stats::quantile(
    maximum_draws,
    probs = 1 - alpha,
    names = FALSE,
    type = 8
  ))
  lower <- moments$estimate - critical_value * standard_error
  observed_standardized <- rep(-Inf, length(moments$estimate))
  observed_standardized[active] <-
    moments$estimate[active] / standard_error[active]
  observed_standardized[!active & moments$estimate > tolerance] <- Inf
  simultaneous_p_value <- function(index) {
    if (!length(index)) return(NA_real_)
    observed <- max(observed_standardized[index])
    if (is.infinite(observed) && observed < 0) return(1)
    (1 + sum(maximum_draws >= observed)) / (length(maximum_draws) + 1)
  }
  details <- data.frame(
    restriction = moments$restriction,
    restriction_id = moments$restriction_id,
    estimate = moments$estimate,
    lower = lower,
    upper = Inf,
    standard_error = standard_error,
    available = TRUE,
    reject = lower > tolerance,
    p_value = vapply(
      seq_along(moments$estimate),
      simultaneous_p_value,
      numeric(1L)
    ),
    alpha = alpha,
    method = if (moments$sampling == "cluster") {
      "direct-delta-cluster-multiplier"
    } else {
      "direct-delta-gaussian"
    },
    stringsAsFactors = FALSE
  )
  finite_max <- function(value) {
    value <- value[is.finite(value)]
    if (length(value)) max(value) else NA_real_
  }
  restrictions <- c(
    "attention-overload",
    "full-attention-probability",
    "block-marschak"
  )
  summaries <- lapply(restrictions, function(restriction) {
    selected_index <- which(details$restriction == restriction)
    selected <- details[selected_index, , drop = FALSE]
    data.frame(
      restriction = restriction,
      available = nrow(selected) > 0L,
      n_restrictions = nrow(selected),
      max_violation_estimate = finite_max(selected$estimate),
      max_violation_lower = finite_max(selected$lower),
      n_rejected = sum(selected$reject),
      reject = any(selected$reject),
      p_value = simultaneous_p_value(selected_index),
      alpha = alpha,
      method = if (moments$sampling == "cluster") {
        "direct-delta-cluster-multiplier"
      } else {
        "direct-delta-gaussian"
      },
      stringsAsFactors = FALSE
    )
  })
  summaries[[length(summaries) + 1L]] <- data.frame(
    restriction = "omnibus",
    available = TRUE,
    n_restrictions = nrow(details),
    max_violation_estimate = finite_max(details$estimate),
    max_violation_lower = finite_max(details$lower),
    n_rejected = sum(details$reject),
    reject = any(details$reject),
    p_value = simultaneous_p_value(seq_len(nrow(details))),
    alpha = alpha,
    method = if (moments$sampling == "cluster") {
      "direct-delta-cluster-multiplier"
    } else {
      "direct-delta-gaussian"
    },
    stringsAsFactors = FALSE
  )
  list(
    summary = do.call(rbind, summaries),
    details = details,
    critical_value = critical_value,
    n_active = sum(active)
  )
}

.hlao_specification_diagnostics <- function(summary, domain, bands, attention,
                                            full_attention, alpha,
                                            tolerance) {
  finite_max <- function(value) {
    value <- value[is.finite(value)]
    if (length(value)) max(value) else NA_real_
  }
  reach_bounds <- .hlao_reach_bounds(domain, bands)
  attention_results <- .hlao_attention_specification(
    domain,
    attention,
    reach_bounds,
    tolerance
  )
  full_attention_bands <- .hlao_full_attention_bands(
    summary,
    domain,
    bands,
    reach_bounds,
    full_attention,
    tolerance
  )
  bm_results <- .hlao_bm_specification(
    summary,
    domain,
    full_attention,
    full_attention_bands,
    tolerance
  )
  details <- rbind(
    attention_results,
    full_attention_bands$probability_diagnostics,
    bm_results
  )
  details$alpha <- alpha
  details$method <- bands$method

  restrictions <- c(
    "attention-overload",
    "full-attention-probability",
    "block-marschak"
  )
  summaries <- lapply(restrictions, function(restriction) {
    selected <- details[details$restriction == restriction, , drop = FALSE]
    available <- nrow(selected) > 0L && all(selected$available)
    data.frame(
      restriction = restriction,
      available = available,
      n_restrictions = nrow(selected),
      max_violation_estimate = if (available) finite_max(selected$estimate) else NA_real_,
      max_violation_lower = if (available) finite_max(selected$lower) else NA_real_,
      n_rejected = if (available) sum(selected$reject) else NA_integer_,
      reject = if (available) any(selected$reject) else NA,
      alpha = alpha,
      method = bands$method,
      stringsAsFactors = FALSE
    )
  })
  available_details <- details[details$available, , drop = FALSE]
  summaries[[length(summaries) + 1L]] <- data.frame(
    restriction = "omnibus",
    available = nrow(available_details) > 0L,
    n_restrictions = nrow(available_details),
    max_violation_estimate = if (nrow(available_details)) {
      finite_max(available_details$estimate)
    } else {
      NA_real_
    },
    max_violation_lower = if (nrow(available_details)) {
      finite_max(available_details$lower)
    } else {
      NA_real_
    },
    n_rejected = if (nrow(available_details)) {
      sum(available_details$reject)
    } else {
      NA_integer_
    },
    reject = if (nrow(available_details)) any(available_details$reject) else NA,
    alpha = alpha,
    method = bands$method,
    stringsAsFactors = FALSE
  )
  list(
    summary = do.call(rbind, summaries),
    details = details,
    reach_bounds = reach_bounds,
    full_attention_bands = full_attention_bands
  )
}

.hlao_pairwise_intervals <- function(summary, domain, bands, tolerance) {
  binary <- which(rowSums(summary$menu) == 2L)
  if (!length(binary)) {
    return(data.frame(
      earlier = integer(0L), later = integer(0L), estimate = numeric(0L),
      reach_estimate = numeric(0L), reach_lower = numeric(0L),
      reach_upper = numeric(0L), lower = numeric(0L), upper = numeric(0L),
      empty = logical(0L)
    ))
  }
  do.call(rbind, lapply(binary, function(menu_index) {
    items <- domain$items[[menu_index]]
    singleton_index <- domain$suffix_index[[menu_index]][2L]
    y <- summary$prob[menu_index, items[2L]]
    y_lower <- bands$lower[menu_index, items[2L]]
    y_upper <- bands$upper[menu_index, items[2L]]
    reach_estimate <-
      (1 - summary$outside_prob[menu_index]) *
      (1 - summary$outside_prob[singleton_index])
    reach_lower <-
      (1 - bands$outside_upper[menu_index]) *
      (1 - bands$outside_upper[singleton_index])
    reach_upper <-
      (1 - bands$outside_lower[menu_index]) *
      (1 - bands$outside_lower[singleton_index])

    if (reach_upper <= tolerance) {
      empty <- y_lower > tolerance
      lower <- if (empty) NA_real_ else 0
      upper <- if (empty) NA_real_ else 1
    } else {
      lower <- max(0, y_lower / reach_upper)
      upper <- if (reach_lower > tolerance) min(1, y_upper / reach_lower) else 1
      empty <- lower > upper + tolerance
      if (empty) {
        lower <- upper <- NA_real_
      }
    }
    data.frame(
      menu_id = menu_index,
      earlier = items[1L],
      later = items[2L],
      estimate = if (reach_estimate > tolerance) y / reach_estimate else NA_real_,
      reach_estimate = reach_estimate,
      reach_lower = reach_lower,
      reach_upper = reach_upper,
      lower = lower,
      upper = upper,
      empty = empty
    )
  }))
}

.hlao_quadratic_set <- function(coefficients, tolerance) {
  coefficients[abs(coefficients) <= tolerance] <- 0
  nonzero <- which(coefficients != 0)
  last_coefficient <- if (length(nonzero)) max(nonzero) else 1L
  if (last_coefficient == 1L) {
    if (coefficients[1L] <= tolerance) {
      return(matrix(c(0, 1), nrow = 1L, dimnames = list(NULL, c("lower", "upper"))))
    }
    return(matrix(numeric(0L), nrow = 0L, ncol = 2L,
                  dimnames = list(NULL, c("lower", "upper"))))
  }
  roots <- polyroot(coefficients[seq_len(last_coefficient)])
  roots <- Re(roots[abs(Im(roots)) <= sqrt(tolerance)])
  roots <- pmin(1, pmax(0, roots[roots >= -tolerance & roots <= 1 + tolerance]))
  knots <- sort(unique(c(0, roots, 1)))
  evaluate <- function(theta) {
    coefficients[1L] + coefficients[2L] * theta +
      coefficients[3L] * theta^2
  }
  pieces <- list()
  if (length(knots) > 1L) {
    for (index in seq_len(length(knots) - 1L)) {
      lower <- knots[index]
      upper <- knots[index + 1L]
      if (evaluate((lower + upper) / 2) <= tolerance) {
        pieces[[length(pieces) + 1L]] <- c(lower, upper)
      }
    }
  }
  for (root in knots[evaluate(knots) <= tolerance]) {
    covered <- length(pieces) && any(vapply(pieces, function(piece) {
      root >= piece[1L] - tolerance && root <= piece[2L] + tolerance
    }, logical(1L)))
    if (!covered) pieces[[length(pieces) + 1L]] <- c(root, root)
  }
  if (!length(pieces)) {
    return(matrix(numeric(0L), nrow = 0L, ncol = 2L,
                  dimnames = list(NULL, c("lower", "upper"))))
  }
  pieces <- pieces[order(vapply(pieces, `[[`, numeric(1L), 1L))]
  merged <- list(pieces[[1L]])
  if (length(pieces) > 1L) {
    for (piece in pieces[-1L]) {
      last <- merged[[length(merged)]]
      if (piece[1L] <= last[2L] + tolerance) {
        merged[[length(merged)]][2L] <- max(last[2L], piece[2L])
      } else {
        merged[[length(merged) + 1L]] <- piece
      }
    }
  }
  result <- do.call(rbind, merged)
  colnames(result) <- c("lower", "upper")
  result
}

.hlao_studentized_pairwise <- function(summary, domain, alpha, tolerance) {
  binary <- which(rowSums(summary$menu) == 2L)
  empty_summary <- data.frame(
    menu_id = integer(0L), earlier = integer(0L), later = integer(0L),
    estimate = numeric(0L), reach_estimate = numeric(0L),
    lower = numeric(0L), upper = numeric(0L), width = numeric(0L),
    n_components = integer(0L), empty = logical(0L),
    studentizable = logical(0L), critical_value = numeric(0L),
    alpha = numeric(0L), method = character(0L)
  )
  empty_components <- data.frame(
    menu_id = integer(0L), earlier = integer(0L), later = integer(0L),
    component = integer(0L), lower = numeric(0L), upper = numeric(0L)
  )
  if (!length(binary)) {
    return(list(summary = empty_summary, components = empty_components))
  }
  critical_value <- stats::qnorm(1 - alpha / (2 * length(binary)))
  summaries <- components <- list()
  component_index <- 0L
  for (pair_index in seq_along(binary)) {
    menu_index <- binary[pair_index]
    items <- domain$items[[menu_index]]
    singleton_index <- domain$suffix_index[[menu_index]][2L]
    y <- summary$prob[menu_index, items[2L]]
    e_ab <- summary$outside_prob[menu_index]
    e_b <- summary$outside_prob[singleton_index]
    reach <- (1 - e_ab) * (1 - e_b)
    cells <- c(
      summary$inside_cell[menu_index, items[2L]],
      summary$outside_cell[menu_index],
      summary$outside_cell[singleton_index]
    )
    covariance <- summary$covariance[cells, cells, drop = FALSE]
    variance_coefficients <- c(
      covariance[1L, 1L],
      2 * (
        (1 - e_b) * covariance[1L, 2L] +
          (1 - e_ab) * covariance[1L, 3L]
      ),
      (1 - e_b)^2 * covariance[2L, 2L] +
        (1 - e_ab)^2 * covariance[3L, 3L] +
        2 * (1 - e_b) * (1 - e_ab) * covariance[2L, 3L]
    )
    quadratic <- c(
      y^2 - critical_value^2 * variance_coefficients[1L],
      -2 * y * reach - critical_value^2 * variance_coefficients[2L],
      reach^2 - critical_value^2 * variance_coefficients[3L]
    )
    accepted <- .hlao_quadratic_set(quadratic, tolerance)
    candidate_variance <- vapply(c(0, 0.5, 1), function(theta) {
      variance_coefficients[1L] + variance_coefficients[2L] * theta +
        variance_coefficients[3L] * theta^2
    }, numeric(1L))
    studentizable <- max(candidate_variance) > tolerance
    if (!studentizable) {
      accepted <- matrix(
        c(0, 1), nrow = 1L,
        dimnames = list(NULL, c("lower", "upper"))
      )
    }
    if (nrow(accepted)) {
      for (component in seq_len(nrow(accepted))) {
        component_index <- component_index + 1L
        components[[component_index]] <- data.frame(
          menu_id = menu_index,
          earlier = items[1L],
          later = items[2L],
          component = component,
          lower = accepted[component, "lower"],
          upper = accepted[component, "upper"]
        )
      }
    }
    summaries[[pair_index]] <- data.frame(
      menu_id = menu_index,
      earlier = items[1L],
      later = items[2L],
      estimate = if (reach > tolerance) y / reach else NA_real_,
      reach_estimate = reach,
      lower = if (nrow(accepted)) min(accepted[, "lower"]) else NA_real_,
      upper = if (nrow(accepted)) max(accepted[, "upper"]) else NA_real_,
      width = if (nrow(accepted)) {
        sum(accepted[, "upper"] - accepted[, "lower"])
      } else {
        0
      },
      n_components = nrow(accepted),
      empty = nrow(accepted) == 0L,
      studentizable = studentizable,
      critical_value = critical_value,
      alpha = alpha,
      method = if (summary$sampling == "cluster") {
        "cluster-studentized-undivided-bonferroni"
      } else {
        "studentized-undivided-bonferroni"
      }
    )
  }
  list(
    summary = do.call(rbind, summaries),
    components = if (length(components)) {
      do.call(rbind, components)
    } else {
      empty_components
    }
  )
}

.hlao_projection_system <- function(summary, domain, bands, rankings) {
  n_rankings <- nrow(rankings)
  positions <- .hlao_ranking_positions(rankings, ncol(summary$menu))
  layout <- vector("list", nrow(summary$menu))
  offset <- n_rankings
  for (menu_index in seq_len(nrow(summary$menu))) {
    n_items <- length(domain$items[[menu_index]])
    v <- offset + seq_len(n_items)
    offset <- max(v)
    mass <- offset + seq_len(n_items + 1L)
    offset <- max(mass)
    q <- matrix(
      offset + seq_len(n_rankings * (n_items + 1L)),
      nrow = n_rankings,
      ncol = n_items + 1L,
      byrow = TRUE
    )
    offset <- max(q)
    layout[[menu_index]] <- list(v = v, mass = mass, q = q)
  }
  n_variables <- offset
  rows <- list()
  directions <- character(0L)
  rhs <- numeric(0L)
  add_constraint <- function(row, direction, value) {
    rows[[length(rows) + 1L]] <<- row
    directions <<- c(directions, direction)
    rhs <<- c(rhs, value)
  }

  row <- numeric(n_variables)
  row[seq_len(n_rankings)] <- 1
  add_constraint(row, "=", 1)

  reach_lower <- reach_upper <- vector("list", nrow(summary$menu))
  for (menu_index in seq_len(nrow(summary$menu))) {
    suffix <- domain$suffix_index[[menu_index]]
    reach_lower[[menu_index]] <- cumprod(1 - bands$outside_upper[suffix])
    reach_upper[[menu_index]] <- cumprod(1 - bands$outside_lower[suffix])
    current <- layout[[menu_index]]
    items <- domain$items[[menu_index]]
    n_items <- length(items)

    for (position in seq_len(n_items)) {
      row <- numeric(n_variables)
      row[current$v[position]] <- 1
      add_constraint(row, ">=", reach_lower[[menu_index]][position])
      add_constraint(row, "<=", reach_upper[[menu_index]][position])
    }
    row <- numeric(n_variables)
    row[current$mass[1L]] <- 1
    row[current$v[1L]] <- 1
    add_constraint(row, "=", 1)
    for (position in seq_len(n_items)) {
      row <- numeric(n_variables)
      row[current$mass[position + 1L]] <- 1
      row[current$v[position]] <- -1
      if (position < n_items) {
        row[current$v[position + 1L]] <- 1
      }
      add_constraint(row, "=", 0)
    }
    for (stop in 0:n_items) {
      row <- numeric(n_variables)
      row[current$q[, stop + 1L]] <- 1
      row[current$mass[stop + 1L]] <- -1
      add_constraint(row, "=", 0)
    }
    for (ranking in seq_len(n_rankings)) {
      row <- numeric(n_variables)
      row[current$q[ranking, ]] <- 1
      row[ranking] <- -1
      add_constraint(row, "=", 0)
    }
    for (alternative in items) {
      row <- numeric(n_variables)
      for (prefix_size in seq_along(items)) {
        prefix <- items[seq_len(prefix_size)]
        winner_position <- max.col(
          -positions[, prefix, drop = FALSE],
          ties.method = "first"
        )
        winners <- prefix[winner_position]
        selected <- which(winners == alternative)
        if (length(selected)) {
          row[current$q[selected, prefix_size + 1L]] <- 1
        }
      }
      add_constraint(row, ">=", bands$lower[menu_index, alternative])
      add_constraint(row, "<=", bands$upper[menu_index, alternative])
    }
  }

  for (small in seq_len(nrow(summary$menu))) {
    for (large in seq_len(nrow(summary$menu))) {
      if (small == large || sum(summary$menu[small, ]) >= sum(summary$menu[large, ]) ||
          !all(summary$menu[small, ] <= summary$menu[large, ])) {
        next
      }
      for (alternative in domain$items[[small]]) {
        row <- numeric(n_variables)
        row[layout[[small]]$v[match(alternative, domain$items[[small]])]] <- 1
        row[layout[[large]]$v[match(alternative, domain$items[[large]])]] <- -1
        add_constraint(row, ">=", 0)
      }
    }
  }

  list(
    A = do.call(rbind, rows),
    directions = directions,
    b = rhs,
    n_rankings = n_rankings,
    n_variables = n_variables,
    reach_lower = reach_lower,
    reach_upper = reach_upper
  )
}

.hlao_solve_projection <- function(system, events) {
  feasibility <- lpSolve::lp(
    "min", rep(0, system$n_variables), system$A,
    system$directions, system$b
  )
  if (feasibility$status != 0L) {
    return(list(
      feasible = FALSE,
      status = feasibility$status,
      intervals = data.frame(
        event = colnames(events), lower = NA_real_, upper = NA_real_,
        lower_status = NA_integer_, upper_status = NA_integer_
      )
    ))
  }
  rows <- vector("list", ncol(events))
  for (event_index in seq_len(ncol(events))) {
    objective <- numeric(system$n_variables)
    objective[seq_len(nrow(events))] <- events[, event_index]
    lower <- lpSolve::lp(
      "min", objective, system$A, system$directions, system$b
    )
    upper <- lpSolve::lp(
      "max", objective, system$A, system$directions, system$b
    )
    rows[[event_index]] <- data.frame(
      event = colnames(events)[event_index],
      lower = if (lower$status == 0L) lower$objval else NA_real_,
      upper = if (upper$status == 0L) upper$objval else NA_real_,
      lower_status = lower$status,
      upper_status = upper$status
    )
  }
  list(
    feasible = TRUE,
    status = feasibility$status,
    intervals = do.call(rbind, rows)
  )
}

#' Inference for Heterogeneous List-Based Attention Overload
#'
#' @description
#' `hlaoTest` forms simultaneous bands for all observed menu--outcome
#' probabilities. The default Hoeffding method is finite-sample valid. The
#' correlated-Gaussian method uses the estimated block-multinomial covariance
#' and retains exact binomial bands for sparse or degenerate cells. When
#' `cluster` is supplied, the Gaussian component instead uses cluster-level
#' influence vectors and multiplier draws, and the fallback is a
#' cluster-Hoeffding band. When both types of cells are present, each component
#' receives half of the common error budget. The function
#' inverts the undivided binary-menu moments to obtain simultaneous pairwise
#' preference-share intervals that remain valid at zero reach. It also reports
#' a Bonferroni-calibrated studentized inversion of the same moments. The
#' studentized set is obtained by exact quadratic inversion and may therefore
#' contain more than one component; exactly degenerate moments return `[0,1]`.
#' For supplied general preference events, the function also computes the
#' dependence-robust outer projection intervals described in the
#' Supplemental Appendix.
#'
#' @param menu Zero-one matrix of menus, with one row per observation.
#' @param choice Zero-one matrix of inside choices. An all-zero row denotes the
#'   outside option unless `outside` is supplied.
#' @param outside Optional zero-one indicator for outside choices.
#' @param cluster Optional vector identifying independent sampling clusters.
#'   When supplied, covariance estimation and Gaussian calibration use
#'   cluster-level influence vectors and multiplier draws.
#' @param list_order Permutation giving the observed presentation order.
#' @param events Optional zero-one event indicators over [hlaoRankings()].
#' @param alpha Nominal error probability for the common simultaneous region.
#' @param band_method Probability-band method, either `"hoeffding"` or
#'   `"gaussian"`.
#' @param diagnostic_method Specification-diagnostic method. `"outer"` uses
#'   the simultaneous probability region and is the finite-sample default.
#'   `"delta"` uses a direct delta-Gaussian approximation on a complete menu
#'   domain with positive terminal reach.
#' @param n_band_draws Number of Gaussian or cluster-multiplier draws used for
#'   the simultaneous band and direct-diagnostic critical values.
#' @param boundary_count Minimum number of successes and failures required for
#'   a cell to use the Gaussian band. Under clustered sampling this counts
#'   clusters with successes and failures. Other cells retain the applicable
#'   simultaneous fallback band.
#' @param tolerance Nonnegative numerical tolerance for zero-reach conventions.
#' @param max_rankings Maximum ranking count used for general-event projection.
#'
#' @return An object of class `ramchoiceHLAOTest` containing aggregated choice
#'   `summary`, plug-in `attention` and `full_attention` estimates, simultaneous
#'   probability `bands`, weak-reach `pairwise` intervals, studentized
#'   `pairwise_studentized` sets and their `pairwise_studentized_components`,
#'   optional general-event `event_intervals`, simultaneous specification
#'   diagnostics, projection dimensions, options, and elapsed time.
#'
#' @references
#' M. D. Cattaneo, P. H. Y. Cheung, X. Ma, and Y. Masatlioglu (2026).
#' \href{https://arxiv.org/abs/2110.10650}{Attention Overload}. Working paper.
#'
#' @examples
#' menu <- rbind(
#'   matrix(rep(c(1, 0), 10), ncol = 2, byrow = TRUE),
#'   matrix(rep(c(0, 1), 10), ncol = 2, byrow = TRUE),
#'   matrix(rep(c(1, 1), 10), ncol = 2, byrow = TRUE)
#' )
#' choice <- matrix(0, nrow = nrow(menu), ncol = 2)
#' choice[1:8, 1] <- 1
#' choice[11:17, 2] <- 1
#' choice[21:26, 1] <- 1
#' choice[27:28, 2] <- 1
#' hlaoTest(menu, choice)
#'
#' @export
hlaoTest <- function(menu, choice, outside = NULL, list_order = NULL,
                     events = NULL, alpha = 0.05,
                     band_method = c("hoeffding", "gaussian"),
                     diagnostic_method = c("outer", "delta"),
                     n_band_draws = 2000L, boundary_count = 5L,
                     tolerance = sqrt(.Machine$double.eps),
                     max_rankings = 5000L, cluster = NULL) {
  started <- proc.time()[["elapsed"]]
  tolerance <- .hlao_validate_tolerance(tolerance)
  diagnostic_method <- match.arg(diagnostic_method)
  summary <- .hlao_sample_summary(menu, choice, outside, cluster)
  list_order <- .hlao_validate_list_order(list_order, ncol(summary$menu))
  domain <- .hlao_domain(summary$menu, list_order)
  bands <- .hlao_probability_bands(
    summary,
    alpha,
    band_method = band_method,
    n_band_draws = n_band_draws,
    boundary_count = boundary_count
  )
  attention <- .hlao_recover_attention(
    summary$menu,
    summary$outside_prob,
    domain,
    tolerance
  )
  full_attention <- .hlao_full_attention(
    summary$menu,
    summary$prob,
    domain,
    attention,
    tolerance
  )
  block_marschak <- .hlao_bm_diagnostics(
    summary$menu,
    domain,
    full_attention,
    tolerance
  )
  outer_specification <- .hlao_specification_diagnostics(
    summary,
    domain,
    bands,
    attention,
    full_attention,
    alpha,
    tolerance
  )
  specification <- if (diagnostic_method == "delta") {
    .hlao_delta_specification(
      summary,
      domain,
      attention,
      full_attention,
      alpha,
      n_band_draws,
      tolerance
    )
  } else {
    outer_specification
  }
  pairwise <- .hlao_pairwise_intervals(summary, domain, bands, tolerance)
  pairwise_studentized <- .hlao_studentized_pairwise(
    summary, domain, alpha, tolerance
  )

  rankings <- NULL
  event_matrix <- NULL
  projection <- NULL
  event_intervals <- data.frame(
    event = character(0L), lower = numeric(0L), upper = numeric(0L),
    lower_status = integer(0L), upper_status = integer(0L)
  )
  if (!is.null(events)) {
    n_rankings <- gamma(ncol(summary$menu) + 1L)
    if (!is.finite(n_rankings) || n_rankings > max_rankings) {
      stop(
        "General-event projection requires ", n_rankings,
        " ranking columns, exceeding 'max_rankings'.",
        call. = FALSE
      )
    }
    rankings <- hlaoRankings(seq_len(ncol(summary$menu)))
    event_matrix <- .hlao_normalize_events(events, rankings)
    projection <- .hlao_projection_system(summary, domain, bands, rankings)
    fit <- .hlao_solve_projection(projection, event_matrix)
    event_intervals <- fit$intervals
    event_intervals$alpha <- alpha
    event_intervals$method <- paste0(
      "dependence-robust-", bands$method, "-projection"
    )
    projection$feasible <- fit$feasible
    projection$status <- fit$status
  }

  pairwise$alpha <- alpha
  pairwise$method <- paste0("weak-reach-", bands$method, "-projection")
  result <- list(
    pairwise = pairwise,
    pairwise_studentized = pairwise_studentized$summary,
    pairwise_studentized_components = pairwise_studentized$components,
    event_intervals = event_intervals,
    attention = attention$attention,
    attention_diagnostics = attention$diagnostics,
    full_attention = full_attention,
    block_marschak = block_marschak,
    specification = specification$summary,
    specification_details = specification$details,
    outer_specification = outer_specification$summary,
    outer_specification_details = outer_specification$details,
    full_attention_bands = outer_specification$full_attention_bands,
    summary = summary,
    bands = bands$results,
    rankings = rankings,
    events = event_matrix,
    projection = projection,
    options = list(
      alpha = alpha,
      band_method = bands$method,
      diagnostic_method = diagnostic_method,
      diagnostic_critical_value = if (diagnostic_method == "delta") {
        specification$critical_value
      } else {
        NA_real_
      },
      n_active_diagnostics = if (diagnostic_method == "delta") {
        specification$n_active
      } else {
        NA_integer_
      },
      n_band_draws = n_band_draws,
      boundary_count = boundary_count,
      band_critical_value = bands$critical_value,
      gaussian_alpha = bands$gaussian_alpha,
      fallback_alpha = bands$fallback_alpha,
      n_gaussian_cells = bands$n_gaussian_cells,
      n_fallback_cells = bands$n_fallback_cells,
      list_order = list_order,
      tolerance = tolerance,
      max_rankings = max_rankings,
      n_cells = bands$n_cells,
      sampling = summary$sampling,
      n_cluster = summary$n_cluster
    ),
    elapsed = unname(proc.time()[["elapsed"]] - started)
  )
  class(result) <- "ramchoiceHLAOTest"
  result
}

#' Path-Independence-Robust Inference for H-LAO
#'
#' @description
#' `hlaoNoPITest` projects a simultaneous confidence region for primitive
#' menu-choice probabilities through the sharp H-LAO model that retains prefix
#' consideration, attention overload, and a stable marginal preference
#' distribution but does not impose Sequential Path Independence. Prefix masses
#' and menu-specific preference--stopping couplings are latent variables. Every
#' reported event endpoint is obtained by linear programming, and the observed
#' menu domain need not be suffix closed.
#'
#' @param menu Zero-one matrix of menus, with one row per observation.
#' @param choice Zero-one matrix of inside choices. An all-zero row denotes the
#'   outside option unless `outside` is supplied.
#' @param outside Optional zero-one indicator for outside choices.
#' @param cluster Optional vector identifying independent sampling clusters.
#'   When supplied, covariance estimation and Gaussian calibration use
#'   cluster-level influence vectors and multiplier draws.
#' @param list_order Permutation giving the observed presentation order.
#' @param events Optional zero-one event indicators over [hlaoRankings()].
#' @param alpha Nominal error probability for the common simultaneous region.
#' @param band_method Probability-band method, either `"hoeffding"` or
#'   `"gaussian"`.
#' @param n_band_draws Number of Gaussian or cluster-multiplier draws used by
#'   the covariance-aware probability band.
#' @param boundary_count Minimum number of successes and failures required for
#'   a cell to use the Gaussian band. Under clustered sampling this counts
#'   clusters with successes and failures.
#' @param tolerance Nonnegative numerical tolerance.
#' @param max_rankings Maximum ranking count used for event projection.
#'
#' @return An object of class `ramchoiceHLAONoPITest` containing event
#'   `intervals`, simultaneous probability `bands`, the LP `projection`,
#'   options, and elapsed time.
#'
#' @references
#' M. D. Cattaneo, P. H. Y. Cheung, X. Ma, and Y. Masatlioglu (2026).
#' \href{https://arxiv.org/abs/2110.10650}{Attention Overload}. Working paper.
#'
#' @examples
#' menu <- rbind(
#'   matrix(rep(c(1, 0), 20), ncol = 2, byrow = TRUE),
#'   matrix(rep(c(1, 1), 20), ncol = 2, byrow = TRUE)
#' )
#' choice <- matrix(0, nrow = nrow(menu), ncol = 2)
#' choice[1:15, 1] <- 1
#' choice[21:30, 1] <- 1
#' choice[31:36, 2] <- 1
#' rankings <- hlaoRankings(1:2)
#' hlaoNoPITest(
#'   menu, choice,
#'   events = list(`2 above 1` = rankings[, 1] == 2)
#' )
#'
#' @export
hlaoNoPITest <- function(menu, choice, outside = NULL, list_order = NULL,
                          events = NULL, alpha = 0.05,
                          band_method = c("hoeffding", "gaussian"),
                          n_band_draws = 2000L, boundary_count = 5L,
                          tolerance = sqrt(.Machine$double.eps),
                          max_rankings = 5000L, cluster = NULL) {
  started <- proc.time()[["elapsed"]]
  tolerance <- .hlao_validate_tolerance(tolerance)
  summary <- .hlao_sample_summary(menu, choice, outside, cluster)
  list_order <- .hlao_validate_list_order(list_order, ncol(summary$menu))
  domain <- .hlao_domain_unrestricted(summary$menu, list_order)
  bands <- .hlao_probability_bands(
    summary,
    alpha,
    band_method = band_method,
    n_band_draws = n_band_draws,
    boundary_count = boundary_count
  )
  n_rankings <- gamma(ncol(summary$menu) + 1L)
  if (!is.finite(n_rankings) || n_rankings > max_rankings) {
    stop(
      "No-SPI event projection requires ", n_rankings,
      " ranking columns, exceeding 'max_rankings'.",
      call. = FALSE
    )
  }
  rankings <- hlaoRankings(seq_len(ncol(summary$menu)))
  event_matrix <- .hlao_normalize_events(events, rankings)
  projection <- .hlao_nopi_system(
    summary$menu,
    domain,
    rankings,
    bands = bands
  )
  fit <- .hlao_solve_general_bounds(projection, event_matrix, "noPI")
  intervals <- fit$bounds
  if (nrow(intervals)) {
    intervals$alpha <- alpha
    intervals$method <- paste0(
      "noPI-", bands$method, "-exact-projection"
    )
  }
  projection$feasible <- fit$feasible
  projection$status <- fit$status
  result <- list(
    intervals = intervals,
    summary = summary,
    bands = bands$results,
    rankings = rankings,
    events = event_matrix,
    projection = projection,
    options = list(
      alpha = alpha,
      band_method = bands$method,
      n_band_draws = n_band_draws,
      boundary_count = boundary_count,
      band_critical_value = bands$critical_value,
      gaussian_alpha = bands$gaussian_alpha,
      fallback_alpha = bands$fallback_alpha,
      n_gaussian_cells = bands$n_gaussian_cells,
      n_fallback_cells = bands$n_fallback_cells,
      list_order = list_order,
      tolerance = tolerance,
      max_rankings = max_rankings,
      n_cells = bands$n_cells,
      sampling = summary$sampling,
      n_cluster = summary$n_cluster
    ),
    elapsed = unname(proc.time()[["elapsed"]] - started)
  )
  class(result) <- "ramchoiceHLAONoPITest"
  result
}

#' @export
summary.ramchoiceHLAONoPITest <- function(object, ...) {
  list(
    feasible = object$projection$feasible,
    intervals = object$intervals
  )
}

#' @export
print.ramchoiceHLAONoPITest <- function(x, ...) {
  cat("\nPath-independence-robust inference for H-LAO\n\n")
  cat("Observations:", sum(x$summary$sample_size), "\n")
  cat("Observed menus:", nrow(x$summary$menu), "\n")
  cat("LP feasible:", x$projection$feasible, "\n")
  cat("Elapsed seconds:", format(round(x$elapsed, 3), nsmall = 3), "\n")
  if (nrow(x$intervals)) {
    cat("\nPreference-event intervals\n\n")
    print(x$intervals, row.names = FALSE)
  }
  invisible(x)
}

#' @export
summary.ramchoiceHLAOModel <- function(object, ...) {
  list(
    compatibility = object$compatibility,
    attention = object$attention_diagnostics,
    bounds = object$bounds,
    pairwise = object$pairwise,
    agreement = object$agreement,
    computation = object$computation
  )
}

#' @export
print.ramchoiceHLAOModel <- function(x, ...) {
  cat("\nPopulation analysis for H-LAO\n\n")
  print(x$compatibility, row.names = FALSE)
  if (nrow(x$bounds)) {
    cat("\nPreference-event bounds\n\n")
    print(x$bounds, row.names = FALSE)
  }
  if (!is.null(x$agreement) && nrow(x$agreement)) {
    cat("\nFull-attention agreement\n\n")
    print(x$agreement, row.names = FALSE)
  }
  if (!is.null(x$computation)) {
    cat("\nComputation\n\n")
    print(x$computation, row.names = FALSE)
  }
  invisible(x)
}

#' @export
summary.ramchoiceHLAOTest <- function(object, ...) {
  list(
    pairwise = object$pairwise,
    pairwise_studentized = object$pairwise_studentized,
    event_intervals = object$event_intervals,
    specification = object$specification
  )
}

#' @export
print.ramchoiceHLAOTest <- function(x, ...) {
  cat("\nFinite-sample inference for H-LAO\n\n")
  cat("Observations:", sum(x$summary$sample_size), "\n")
  cat("Observed menus:", nrow(x$summary$menu), "\n")
  cat("Elapsed seconds:", format(round(x$elapsed, 3), nsmall = 3), "\n\n")
  print(x$pairwise, row.names = FALSE)
  if (nrow(x$pairwise_studentized)) {
    cat("\nStudentized pairwise confidence sets\n\n")
    print(x$pairwise_studentized, row.names = FALSE)
  }
  cat("\nSpecification diagnostics\n\n")
  print(x$specification, row.names = FALSE)
  if (nrow(x$event_intervals)) {
    cat("\nPreference-event intervals\n\n")
    print(x$event_intervals, row.names = FALSE)
  }
  invisible(x)
}
