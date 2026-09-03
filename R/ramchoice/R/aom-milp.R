################################################################################
# Mixed-integer identification for the Homogeneous Attention Overload Model
################################################################################

.aom_status_class <- function(status) {
  if (identical(as.integer(status), 0L)) {
    "success"
  } else if (identical(as.integer(status), 2L)) {
    "infeasible"
  } else {
    "solver-error"
  }
}

.aom_validate_population <- function(menu, prob, tolerance) {
  if (!is.matrix(menu) || !is.numeric(menu) && !is.logical(menu) ||
      min(dim(menu)) == 0L || anyNA(menu) ||
      any(!(menu %in% c(0, 1)))) {
    stop("'menu' must be a nonempty zero-one matrix.", call. = FALSE)
  }
  menu <- matrix(as.integer(menu), nrow = nrow(menu), ncol = ncol(menu))
  if (any(rowSums(menu) == 0L)) {
    stop("Every row of 'menu' must contain an available alternative.", call. = FALSE)
  }
  if (!is.matrix(prob) || !is.numeric(prob) ||
      !identical(dim(prob), dim(menu)) || anyNA(prob) ||
      any(!is.finite(prob)) || any(prob < 0) || any(prob > 1)) {
    stop(
      "'prob' must be a finite probability matrix with the same dimensions as 'menu'.",
      call. = FALSE
    )
  }
  if (any(abs(prob[menu == 0L]) > tolerance)) {
    stop("Choice probabilities must be zero outside each menu.", call. = FALSE)
  }
  if (any(abs(rowSums(prob) - 1) > tolerance)) {
    stop("Choice probabilities must sum to one on every menu.", call. = FALSE)
  }

  keys <- apply(menu, 1L, paste0, collapse = "")
  for (key in unique(keys[duplicated(keys)])) {
    rows <- which(keys == key)
    differences <- abs(
      prob[rows, , drop = FALSE] -
        matrix(prob[rows[1L], ], nrow = length(rows), ncol = ncol(prob), byrow = TRUE)
    )
    if (max(differences) > tolerance) {
      stop("Duplicate menu rows contain conflicting choice probabilities.", call. = FALSE)
    }
  }
  keep <- !duplicated(keys)
  menu <- menu[keep, , drop = FALSE]
  prob <- prob[keep, , drop = FALSE]
  list(menu = menu, prob = prob)
}

.aom_milp_system <- function(menu, prob, tolerance) {
  n_alternatives <- ncol(menu)
  z_index <- matrix(NA_integer_, n_alternatives, n_alternatives)
  current <- 0L
  for (preferred in seq_len(n_alternatives)) {
    for (inferior in seq_len(n_alternatives)) {
      if (preferred != inferior) {
        current <- current + 1L
        z_index[preferred, inferior] <- current
      }
    }
  }

  rows <- list()
  directions <- character(0L)
  rhs <- numeric(0L)
  labels <- character(0L)
  add_constraint <- function(row, direction, value, label) {
    rows[[length(rows) + 1L]] <<- row
    directions <<- c(directions, direction)
    rhs <<- c(rhs, value)
    labels <<- c(labels, label)
  }

  for (a in seq_len(n_alternatives - 1L)) {
    for (b in (a + 1L):n_alternatives) {
      row <- numeric(current)
      row[z_index[a, b]] <- 1
      row[z_index[b, a]] <- 1
      add_constraint(row, "=", 1, paste0("totality:", a, ":", b))
    }
  }

  if (n_alternatives >= 3L) {
    for (a in seq_len(n_alternatives)) {
      for (b in setdiff(seq_len(n_alternatives), a)) {
        for (c in setdiff(seq_len(n_alternatives), c(a, b))) {
          row <- numeric(current)
          row[z_index[a, b]] <- 1
          row[z_index[b, c]] <- 1
          row[z_index[a, c]] <- -1
          add_constraint(
            row, "<=", 1,
            paste0("transitivity:", a, ":", b, ":", c)
          )
        }
      }
    }
  }

  for (small in seq_len(nrow(menu))) {
    for (large in seq_len(nrow(menu))) {
      if (small == large || rowSums(menu)[small] >= rowSums(menu)[large] ||
          !all(menu[small, ] <= menu[large, ])) {
        next
      }
      items <- which(menu[small, ] == 1L)
      for (a in items) {
        row <- numeric(current)
        for (b in setdiff(items, a)) {
          row[z_index[b, a]] <- -prob[small, b]
        }
        add_constraint(
          row,
          "<=",
          prob[small, a] - prob[large, a] + tolerance,
          paste0("regularity:", a, ":", small, ":", large)
        )
      }
    }
  }

  list(
    A = do.call(rbind, rows),
    directions = directions,
    b = rhs,
    labels = labels,
    z_index = z_index,
    n_variables = current
  )
}

.aom_solve_milp <- function(system, fixed_index = NULL, fixed_value = NULL) {
  A <- system$A
  directions <- system$directions
  rhs <- system$b
  if (!is.null(fixed_index)) {
    row <- numeric(system$n_variables)
    row[fixed_index] <- 1
    A <- rbind(A, row)
    directions <- c(directions, "=")
    rhs <- c(rhs, fixed_value)
  }
  lpSolve::lp(
    direction = "min",
    objective.in = numeric(system$n_variables),
    const.mat = A,
    const.dir = directions,
    const.rhs = rhs,
    all.bin = TRUE
  )
}

.aom_solution_ranking <- function(solution, z_index) {
  n_alternatives <- nrow(z_index)
  wins <- numeric(n_alternatives)
  for (alternative in seq_len(n_alternatives)) {
    indices <- z_index[alternative, ]
    wins[alternative] <- sum(solution[indices[!is.na(indices)]])
  }
  order(wins, decreasing = TRUE)
}

#' Population Identification for Homogeneous AOM
#'
#' @description
#' `aomIdentify` implements the mixed-integer characterization of homogeneous
#' Attention Overload. Binary variables encode pairwise comparisons, while
#' totality, transitivity, and the observed `succ`-Regularity inequalities
#' characterize the sharp set of compatible strict preferences. The routine
#' tests model feasibility and all pairwise revealed-preference conclusions
#' without enumerating the factorial collection of rankings. A status of 2 is
#' treated as a solver-certified infeasibility result; any other nonzero status
#' is reported as a solver error rather than as model incompatibility.
#'
#' @param menu Numeric zero-one matrix with one row per observed menu.
#' @param prob Numeric matrix of population choice probabilities with the same
#'   dimensions as `menu`.
#' @param tolerance Nonnegative numerical tolerance added to the population
#'   inequalities.
#' @param pairwise Logical; if `TRUE`, determine whether each direction of
#'   every pairwise comparison occurs in a compatible preference.
#'
#' @return An object of class `ramchoiceAOMIdentification`. It contains model
#'   `compatible`, one feasible `preference` when the model is nonempty,
#'   pairwise possibility and revelation results, solver diagnostics, and the
#'   mixed-integer system used in the calculation.
#'
#' @references
#' M. D. Cattaneo, P. H. Y. Cheung, X. Ma, and Y. Masatlioglu (2026).
#' \href{https://arxiv.org/abs/2110.10650}{Attention Overload}. Working paper.
#'
#' @examples
#' menu <- prob <- matrix(c(
#'   1, 1, 1,
#'   1, 1, 0,
#'   1, 0, 1,
#'   0, 1, 1
#' ), ncol = 3, byrow = TRUE)
#' for (i in seq_len(nrow(prob))) {
#'   prob[i, menu[i, ] == 1] <- logitAtte(sum(menu[i, ]), 2)$choiceProb
#' }
#' aomIdentify(menu, prob)
#'
#' @export
aomIdentify <- function(menu, prob,
                        tolerance = sqrt(.Machine$double.eps),
                        pairwise = TRUE) {
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance < 0) {
    stop("'tolerance' must be one finite nonnegative number.", call. = FALSE)
  }
  if (!is.logical(pairwise) || length(pairwise) != 1L || is.na(pairwise)) {
    stop("'pairwise' must be TRUE or FALSE.", call. = FALSE)
  }
  prepared <- .aom_validate_population(menu, prob, tolerance)
  menu <- prepared$menu
  prob <- prepared$prob
  system <- .aom_milp_system(menu, prob, tolerance)

  started <- proc.time()[["elapsed"]]
  baseline <- .aom_solve_milp(system)
  baseline_class <- .aom_status_class(baseline$status)
  compatible <- if (baseline_class == "success") {
    TRUE
  } else if (baseline_class == "infeasible") {
    FALSE
  } else {
    NA
  }
  preference <- if (isTRUE(compatible)) {
    .aom_solution_ranking(baseline$solution, system$z_index)
  } else {
    integer(0L)
  }

  pairwise_results <- data.frame(
    alternative_a = integer(0L),
    alternative_b = integer(0L),
    a_preferred_possible = logical(0L),
    b_preferred_possible = logical(0L),
    opposite_status = integer(0L),
    opposite_status_class = character(0L),
    revealed_preference = character(0L),
    stringsAsFactors = FALSE
  )
  n_solves <- 1L
  solver_statuses <- baseline$status
  if (isTRUE(compatible) && pairwise && ncol(menu) >= 2L) {
    rows <- vector("list", choose(ncol(menu), 2L))
    row_index <- 0L
    for (a in seq_len(ncol(menu) - 1L)) {
      for (b in (a + 1L):ncol(menu)) {
        row_index <- row_index + 1L
        baseline_a <- baseline$solution[system$z_index[a, b]] > 0.5
        opposite_index <- if (baseline_a) {
          system$z_index[b, a]
        } else {
          system$z_index[a, b]
        }
        opposite <- .aom_solve_milp(system, opposite_index, 1)
        n_solves <- n_solves + 1L
        solver_statuses <- c(solver_statuses, opposite$status)
        opposite_class <- .aom_status_class(opposite$status)
        opposite_possible <- if (opposite_class == "success") {
          TRUE
        } else if (opposite_class == "infeasible") {
          FALSE
        } else {
          NA
        }
        a_possible <- if (baseline_a) TRUE else opposite_possible
        b_possible <- if (baseline_a) opposite_possible else TRUE
        rows[[row_index]] <- data.frame(
          alternative_a = a,
          alternative_b = b,
          a_preferred_possible = a_possible,
          b_preferred_possible = b_possible,
          opposite_status = opposite$status,
          opposite_status_class = opposite_class,
          revealed_preference = if (isTRUE(a_possible) && identical(b_possible, FALSE)) {
            paste0(a, " > ", b)
          } else if (identical(a_possible, FALSE) && isTRUE(b_possible)) {
            paste0(b, " > ", a)
          } else {
            NA_character_
          },
          stringsAsFactors = FALSE
        )
      }
    }
    pairwise_results <- do.call(rbind, rows)
  }

  result <- list(
    compatible = compatible,
    preference = preference,
    pairwise = pairwise_results,
    status = baseline$status,
    diagnostics = data.frame(
      algorithm = "mixed-integer",
      n_alternatives = ncol(menu),
      n_binary_variables = system$n_variables,
      n_constraints = nrow(system$A),
      n_milp_solves = n_solves,
      baseline_status = baseline$status,
      baseline_status_class = baseline_class,
      all_solver_statuses_resolved = all(solver_statuses %in% c(0L, 2L)),
      elapsed = unname(proc.time()[["elapsed"]] - started),
      stringsAsFactors = FALSE
    ),
    menu = menu,
    prob = prob,
    tolerance = tolerance,
    system = system
  )
  class(result) <- "ramchoiceAOMIdentification"
  result
}

#' @export
summary.ramchoiceAOMIdentification <- function(object, ...) {
  list(
    compatible = object$compatible,
    preference = object$preference,
    pairwise = object$pairwise,
    diagnostics = object$diagnostics
  )
}

#' @export
print.ramchoiceAOMIdentification <- function(x, ...) {
  cat("\nPopulation identification for homogeneous AOM\n\n")
  cat("Compatible:", x$compatible, "\n")
  cat("Solver status:", x$diagnostics$baseline_status_class,
      "(", x$diagnostics$baseline_status, ")\n")
  if (isTRUE(x$compatible)) {
    cat("One compatible preference:", paste(x$preference, collapse = " > "), "\n")
  }
  cat("Binary variables:", x$diagnostics$n_binary_variables, "\n")
  cat("Constraints:", x$diagnostics$n_constraints, "\n")
  cat("MILP solves:", x$diagnostics$n_milp_solves, "\n")
  cat("Elapsed seconds:", format(round(x$diagnostics$elapsed, 3), nsmall = 3), "\n")
  if (nrow(x$pairwise)) {
    cat("\nPairwise identified set\n\n")
    print(x$pairwise, row.names = FALSE)
  }
  invisible(x)
}
