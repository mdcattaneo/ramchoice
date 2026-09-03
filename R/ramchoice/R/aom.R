################################################################################
# Homogeneous Attention Overload Model
################################################################################

.aom_preference_labels <- function(preferences) {
  apply(preferences, 1L, paste, collapse = " > ")
}

.aom_constraint_counts <- function(result, n_preferences) {
  counts <- result$constraints$ConstN
  if (is.null(counts)) {
    counts <- rep.int(0L, n_preferences)
  }
  as.integer(counts)
}

.aom_split_inequalities <- function(values, counts) {
  endpoints <- cumsum(c(0L, counts))
  lapply(seq_along(counts), function(index) {
    if (counts[index] == 0L) {
      return(numeric(0L))
    }
    values[(endpoints[index] + 1L):endpoints[index + 1L]]
  })
}

.aom_validate_alpha <- function(alpha) {
  supported <- c(0.10, 0.05, 0.01)
  if (!is.numeric(alpha) || length(alpha) == 0L || any(!is.finite(alpha))) {
    stop("'alpha' must contain one or more finite numeric values.", call. = FALSE)
  }
  matched <- vapply(alpha, function(value) {
    distance <- abs(supported - value)
    if (min(distance) > 1e-12) {
      return(NA_real_)
    }
    supported[which.min(distance)]
  }, numeric(1L))
  if (anyNA(matched)) {
    stop("'alpha' must use 0.10, 0.05, or 0.01.", call. = FALSE)
  }
  unique(matched)
}

.aom_validate_method <- function(method) {
  if (!is.character(method) || length(method) != 1L || is.na(method)) {
    stop("'method' must be one character string.", call. = FALSE)
  }
  method <- toupper(method)
  supported <- c("GMS", "PI", "LF", "2MS", "2UB", "ALL")
  if (!(method %in% supported)) {
    stop(
      "'method' must be GMS, PI, LF, 2MS, 2UB, or ALL.",
      call. = FALSE
    )
  }
  method
}

#' Population Analysis for the Homogeneous Attention Overload Model
#'
#' @description
#' `aomModel` evaluates the population choice-probability inequalities implied
#' by a collection of candidate preference orderings under the homogeneous
#' Attention Overload Model (AOM). It provides a model-specific interface to
#' the AOM restrictions implemented by [revealPrefModel()].
#'
#' @param menu Numeric matrix of zeros and ones. Each row identifies an
#'   observed menu.
#' @param prob Numeric matrix of choice probabilities with the same dimensions
#'   as `menu`.
#' @param pref_list Numeric matrix whose rows are candidate strict preference
#'   orderings. The default is `1, 2, ...`.
#' @param tolerance Nonnegative numerical tolerance used when classifying a
#'   population inequality as violated.
#' @param attBinary Numeric value between one half and one. Values below one
#'   impose the attentive-at-binaries restriction used by the legacy API.
#'
#' @return An object of class `ramchoiceAOMModel`. Its `results` component has
#'   one row per candidate preference, including compatibility, inequality
#'   counts, and violation magnitudes. The object also contains `preferences`,
#'   candidate-specific `inequalities`, the classification `tolerance`, and the
#'   complete legacy [revealPrefModel()] result.
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
#' aomModel(menu, prob, pref_list = rbind(1:3, 3:1))
#'
#' @export
aomModel <- function(menu, prob, pref_list = NULL,
                     tolerance = sqrt(.Machine$double.eps),
                     attBinary = 1) {
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance < 0) {
    stop("'tolerance' must be one finite nonnegative number.", call. = FALSE)
  }

  legacy <- revealPrefModel(
    menu = menu,
    prob = prob,
    pref_list = pref_list,
    RAM = FALSE,
    AOM = TRUE,
    limDataCorr = TRUE,
    attBinary = attBinary
  )

  preferences <- legacy$pref
  n_preferences <- nrow(preferences)
  counts <- .aom_constraint_counts(legacy, n_preferences)
  inequalities <- .aom_split_inequalities(
    as.numeric(legacy$inequalities$R),
    counts
  )
  labels <- .aom_preference_labels(preferences)

  results <- do.call(rbind, lapply(seq_len(n_preferences), function(index) {
    values <- inequalities[[index]]
    n_violated <- sum(values > tolerance)
    data.frame(
      preference_id = index,
      preference = labels[index],
      compatible = n_violated == 0L,
      n_inequalities = counts[index],
      n_violated = n_violated,
      max_inequality = if (length(values)) max(values) else NA_real_,
      max_violation = if (length(values)) max(c(values, 0)) else 0,
      stringsAsFactors = FALSE
    )
  }))

  result <- list(
    results = results,
    preferences = preferences,
    inequalities = inequalities,
    tolerance = tolerance,
    legacy = legacy
  )
  class(result) <- "ramchoiceAOMModel"
  result
}

.aom_cluster_summary <- function(menu, choice, cluster, summary) {
  cluster <- .ram_validate_cluster(cluster, nrow(menu))
  menu_key <- apply(menu, 1L, paste0, collapse = "")
  summary_key <- apply(summary$sumMenu, 1L, paste0, collapse = "")
  menu_group <- match(menu_key, summary_key)
  if (anyNA(menu_group)) {
    stop("Could not align observations with the aggregated menus.", call. = FALSE)
  }

  n_cell <- sum(summary$sumMsize)
  scores <- matrix(0, nrow = cluster$n, ncol = n_cell)
  next_cell <- 0L
  for (menu_index in seq_len(nrow(summary$sumMenu))) {
    items <- which(summary$sumMenu[menu_index, ] == 1L)
    current <- next_cell + seq_along(items)
    next_cell <- max(current)
    selected_menu <- menu_group == menu_index
    menu_size <- sum(selected_menu)
    probabilities <- summary$sumProb[menu_index, items]
    for (cluster_index in seq_len(cluster$n)) {
      selected <- selected_menu & cluster$id == cluster_index
      cluster_menu_size <- sum(selected)
      if (cluster_menu_size) {
        scores[cluster_index, current] <- (
          colSums(choice[selected, items, drop = FALSE]) -
            cluster_menu_size * probabilities
        ) / menu_size
      }
    }
  }

  covariance <- .ram_cluster_covariance(scores)
  summary$Sigma <- covariance * cluster$n
  summary$covariance <- covariance
  summary$cluster_scores <- scores
  summary$n_cluster <- cluster$n
  summary$cluster_labels <- cluster$labels
  summary$sampling <- "cluster"
  summary
}

.aom_cluster_results <- function(summary, constraints, preferences, method,
                                 alpha, n_draws, baratio_ms, baratio_ub,
                                 mnratio_gms) {
  n_preferences <- nrow(preferences)
  counts <- as.integer(constraints$ConstN)
  methods <- if (method == "ALL") {
    c("GMS", "PI", "LF", "2MS", "2UB")
  } else {
    method
  }
  if (is.null(mnratio_gms)) {
    mnratio_gms <- 1 / log(max(summary$n_cluster, 3L))
  }
  primitive_draws <- .ram_cluster_multiplier_draws(
    summary$cluster_scores,
    n_draws
  )
  labels <- .aom_preference_labels(preferences)
  rows <- list()
  row_index <- 0L
  constraint_start <- 0L

  for (preference_index in seq_len(n_preferences)) {
    n_constraint <- counts[preference_index]
    if (!n_constraint) {
      moment <- standard_error <- numeric(0L)
      statistic <- 0
      standardized_draws <- matrix(0, nrow = n_draws, ncol = 0L)
      t_statistic <- numeric(0L)
    } else {
      indices <- constraint_start + seq_len(n_constraint)
      restriction <- constraints$R[indices, , drop = FALSE]
      moment <- as.numeric(restriction %*% summary$sumProbVec)
      moment_covariance <- restriction %*% summary$covariance %*%
        t(restriction)
      standard_error <- sqrt(pmax(0, diag(moment_covariance)))
      active <- standard_error > sqrt(.Machine$double.eps)
      t_statistic <- numeric(n_constraint)
      t_statistic[active] <- moment[active] / standard_error[active]
      t_statistic[!active & moment > 0] <- Inf
      statistic <- max(c(0, t_statistic))
      standardized_draws <- matrix(0, nrow = n_draws, ncol = n_constraint)
      if (any(active)) {
        moment_draws <- primitive_draws %*% t(restriction[active, , drop = FALSE])
        standardized_draws[, active] <- sweep(
          moment_draws,
          2L,
          standard_error[active],
          "/"
        )
      }
    }
    constraint_start <- constraint_start + n_constraint

    least_favorable <- pmax(0, .ram_row_maximum(standardized_draws))
    distributions <- list(
      LF = least_favorable,
      PI = pmax(0, .ram_row_maximum(sweep(
        standardized_draws,
        2L,
        pmin(t_statistic, 0),
        "+"
      ))),
      GMS = pmax(0, .ram_row_maximum(sweep(
        standardized_draws,
        2L,
        pmin(t_statistic * sqrt(mnratio_gms), 0),
        "+"
      )))
    )

    for (method_index in methods) {
      for (level in alpha) {
        if (method_index %in% c("LF", "PI", "GMS")) {
          distribution <- distributions[[method_index]]
          critical_value <- as.numeric(stats::quantile(
            distribution,
            1 - level,
            names = FALSE
          ))
          p_value <- mean(distribution > statistic)
        } else if (method_index == "2MS") {
          beta_critical <- as.numeric(stats::quantile(
            least_favorable,
            1 - level * baratio_ms,
            names = FALSE
          ))
          selected <- t_statistic > -2 * beta_critical
          distribution <- if (any(selected)) {
            pmax(
              0,
              .ram_row_maximum(standardized_draws[, selected, drop = FALSE])
            )
          } else {
            numeric(n_draws)
          }
          critical_value <- as.numeric(stats::quantile(
            distribution,
            1 - level + 2 * level * baratio_ms,
            names = FALSE
          ))
          p_value <- NA_real_
        } else {
          beta_critical <- as.numeric(stats::quantile(
            least_favorable,
            1 - level * baratio_ub,
            names = FALSE
          ))
          center <- pmin(t_statistic + beta_critical, 0)
          distribution <- pmax(0, .ram_row_maximum(sweep(
            standardized_draws,
            2L,
            center,
            "+"
          )))
          critical_value <- as.numeric(stats::quantile(
            distribution,
            1 - level + level * baratio_ub,
            names = FALSE
          ))
          p_value <- NA_real_
        }
        row_index <- row_index + 1L
        rows[[row_index]] <- data.frame(
          preference_id = preference_index,
          preference = labels[preference_index],
          method = method_index,
          alpha = level,
          statistic = statistic,
          critical_value = critical_value,
          p_value = p_value,
          reject = statistic > critical_value,
          n_inequalities = n_constraint,
          n_positive_sample_inequalities = sum(moment > 0),
          max_sample_inequality = if (n_constraint) max(moment) else NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

#' Sample Inference for the Homogeneous Attention Overload Model
#'
#' @description
#' `aomTest` tests candidate preference orderings under the homogeneous AOM and
#' returns a tidy inference table. Row-i.i.d. calculations delegate to
#' [revealPref()] with AOM restrictions only. When `cluster` is supplied, the
#' function uses cluster-level influence vectors and multiplier critical
#' values while preserving the legacy result for backward-compatible auditing.
#'
#' @param menu Numeric matrix of zeros and ones containing observed menus.
#' @param choice Numeric matrix of zeros and ones containing observed choices.
#' @param pref_list Numeric matrix whose rows are candidate strict preference
#'   orderings. The default is `1, 2, ...`.
#' @param method Critical-value method: `"GMS"`, `"PI"`, `"LF"`, `"2MS"`,
#'   `"2UB"`, or `"ALL"`.
#' @param alpha One or more nominal test levels chosen from `0.10`, `0.05`, and
#'   `0.01`.
#' @param nCritSimu Number of Gaussian or cluster-multiplier simulations used
#'   for critical values.
#' @param BARatio2MS Beta-to-alpha ratio for two-step moment selection.
#' @param BARatio2UB Beta-to-alpha ratio for the two-step upper-bound method.
#' @param MNRatioGMS Generalized moment-selection tuning parameter. `NULL`
#'   uses `1/log(N)`, where `N` is the total sample size under row-i.i.d.
#'   sampling and the number of clusters under clustered sampling.
#' @param attBinary Numeric value between one half and one. Values below one
#'   impose the attentive-at-binaries restriction used by the legacy API.
#' @param cluster Optional vector identifying independent sampling clusters.
#'   When supplied, covariance estimation and Gaussian critical values use
#'   cluster-level influence vectors and multiplier draws.
#'
#' @return An object of class `ramchoiceAOMTest`. Its `results` component has
#'   one row per preference, method, and nominal level. The object also contains
#'   `preferences`, candidate-specific `inequalities`, menu-level `summary`
#'   estimates, `constraints`, inference `options`, elapsed computation time,
#'   and the complete legacy [revealPref()] result.
#'
#' @references
#' M. D. Cattaneo, P. H. Y. Cheung, X. Ma, and Y. Masatlioglu (2026).
#' \href{https://arxiv.org/abs/2110.10650}{Attention Overload}. Working paper.
#'
#' @examples
#' set.seed(42)
#' simulated <- lapply(4:2, function(size) {
#'   logitSimu(n = 10, uSize = 4, mSize = size, a = 2)
#' })
#' menu <- do.call(rbind, lapply(simulated, `[[`, "menu"))
#' choice <- do.call(rbind, lapply(simulated, `[[`, "choice"))
#' aomTest(
#'   menu,
#'   choice,
#'   pref_list = rbind(1:4, 4:1),
#'   nCritSimu = 100
#' )
#'
#' @export
aomTest <- function(menu, choice, pref_list = NULL, method = "GMS",
                    alpha = 0.05, nCritSimu = 2000,
                    BARatio2MS = 0.1, BARatio2UB = 0.1,
                    MNRatioGMS = NULL, attBinary = 1, cluster = NULL) {
  method <- .aom_validate_method(method)
  alpha <- .aom_validate_alpha(alpha)

  started <- proc.time()[["elapsed"]]
  legacy <- revealPref(
    menu = menu,
    choice = choice,
    pref_list = pref_list,
    method = method,
    nCritSimu = nCritSimu,
    BARatio2MS = BARatio2MS,
    BARatio2UB = BARatio2UB,
    MNRatioGMS = MNRatioGMS,
    RAM = FALSE,
    AOM = TRUE,
    limDataCorr = TRUE,
    attBinary = attBinary
  )
  elapsed <- unname(proc.time()[["elapsed"]] - started)

  if (!is.null(cluster)) {
    cluster_summary <- .aom_cluster_summary(
      menu,
      choice,
      cluster,
      legacy$sumStats
    )
    cluster_results <- .aom_cluster_results(
      cluster_summary,
      legacy$constraints,
      legacy$pref,
      method,
      alpha,
      nCritSimu,
      BARatio2MS,
      BARatio2UB,
      MNRatioGMS
    )
    counts <- .aom_constraint_counts(legacy, nrow(legacy$pref))
    inequality_values <- as.numeric(
      legacy$constraints$R %*% legacy$sumStats$sumProbVec
    )
    result <- list(
      results = cluster_results,
      preferences = legacy$pref,
      inequalities = .aom_split_inequalities(inequality_values, counts),
      summary = cluster_summary,
      constraints = legacy$constraints,
      options = list(
        method = method,
        alpha = alpha,
        nCritSimu = nCritSimu,
        BARatio2MS = BARatio2MS,
        BARatio2UB = BARatio2UB,
        MNRatioGMS = MNRatioGMS,
        attBinary = attBinary,
        sampling = "cluster",
        n_cluster = cluster_summary$n_cluster
      ),
      elapsed = unname(proc.time()[["elapsed"]] - started),
      legacy = legacy
    )
    class(result) <- "ramchoiceAOMTest"
    return(result)
  }

  preferences <- legacy$pref
  n_preferences <- nrow(preferences)
  counts <- .aom_constraint_counts(legacy, n_preferences)
  inequality_values <- as.numeric(
    legacy$constraints$R %*% legacy$sumStats$sumProbVec
  )
  inequalities <- .aom_split_inequalities(inequality_values, counts)
  n_positive <- vapply(inequalities, function(values) sum(values > 0), integer(1L))
  max_inequality <- vapply(inequalities, function(values) {
    if (length(values)) max(values) else NA_real_
  }, numeric(1L))

  methods <- if (method == "ALL") {
    c("GMS", "PI", "LF", "2MS", "2UB")
  } else {
    method
  }
  critical_fields <- c(GMS = "GMS", PI = "PI", LF = "LF", `2MS` = "MS", `2UB` = "UB")
  p_value_fields <- c(GMS = "GMS", PI = "PI", LF = "LF")
  labels <- .aom_preference_labels(preferences)

  rows <- list()
  row_index <- 0L
  for (method_index in methods) {
    critical <- legacy$critVal[[critical_fields[[method_index]]]]
    available_alpha <- as.numeric(colnames(critical))
    p_values <- if (method_index %in% names(p_value_fields)) {
      as.numeric(legacy$pVal[[p_value_fields[[method_index]]]])
    } else {
      rep.int(NA_real_, n_preferences)
    }

    for (level in alpha) {
      column <- which.min(abs(available_alpha - level))
      critical_value <- as.numeric(critical[, column])
      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        preference_id = seq_len(n_preferences),
        preference = labels,
        method = method_index,
        alpha = level,
        statistic = as.numeric(legacy$Tstat),
        critical_value = critical_value,
        p_value = p_values,
        reject = as.numeric(legacy$Tstat) > critical_value,
        n_inequalities = counts,
        n_positive_sample_inequalities = n_positive,
        max_sample_inequality = max_inequality,
        stringsAsFactors = FALSE
      )
    }
  }

  result <- list(
    results = do.call(rbind, rows),
    preferences = preferences,
    inequalities = inequalities,
    summary = legacy$sumStats,
    constraints = legacy$constraints,
    options = list(
      method = method,
      alpha = alpha,
      nCritSimu = nCritSimu,
      BARatio2MS = BARatio2MS,
      BARatio2UB = BARatio2UB,
      MNRatioGMS = MNRatioGMS,
      attBinary = attBinary,
      sampling = "iid",
      n_cluster = NA_integer_
    ),
    elapsed = elapsed,
    legacy = legacy
  )
  class(result) <- "ramchoiceAOMTest"
  result
}

#' Sample Inference for the Random Attention Model
#'
#' @description
#' `ramTest` provides tidy candidate-ranking inference for the Random Attention
#' Model of Cattaneo, Ma, Masatlioglu, and Suleymanov (2020). Row-i.i.d.
#' calculations use [revealPref()]. When `cluster` is supplied, the function
#' retains the same RAM inequalities but estimates their joint covariance from
#' cluster influence vectors and uses multiplier critical values.
#'
#' @inheritParams aomTest
#' @param limDataCorr Logical indicating whether to use the limited-menu-domain
#'   correction from the legacy RAM implementation.
#'
#' @return An object of class `ramchoiceRAMTest` with the same tidy components
#'   as [aomTest()] and a complete legacy [revealPref()] result.
#'
#' @references
#' M. D. Cattaneo, X. Ma, Y. Masatlioglu, and E. Suleymanov (2020).
#' A Random Attention Model. Journal of Political Economy 128(7): 2796--2836.
#' \doi{10.1086/706861}
#'
#' @export
ramTest <- function(menu, choice, pref_list = NULL, method = "GMS",
                    alpha = 0.05, nCritSimu = 2000,
                    BARatio2MS = 0.1, BARatio2UB = 0.1,
                    MNRatioGMS = NULL, attBinary = 1,
                    limDataCorr = TRUE, cluster = NULL) {
  method <- .aom_validate_method(method)
  alpha <- .aom_validate_alpha(alpha)
  if (!is.logical(limDataCorr) || length(limDataCorr) != 1L ||
      is.na(limDataCorr)) {
    stop("'limDataCorr' must be one nonmissing logical value.", call. = FALSE)
  }

  started <- proc.time()[["elapsed"]]
  legacy <- revealPref(
    menu = menu,
    choice = choice,
    pref_list = pref_list,
    method = method,
    nCritSimu = nCritSimu,
    BARatio2MS = BARatio2MS,
    BARatio2UB = BARatio2UB,
    MNRatioGMS = MNRatioGMS,
    RAM = TRUE,
    AOM = FALSE,
    limDataCorr = limDataCorr,
    attBinary = attBinary
  )
  preferences <- legacy$pref
  n_preferences <- nrow(preferences)
  counts <- .aom_constraint_counts(legacy, n_preferences)
  inequality_values <- as.numeric(
    legacy$constraints$R %*% legacy$sumStats$sumProbVec
  )
  inequalities <- .aom_split_inequalities(inequality_values, counts)

  if (!is.null(cluster)) {
    summary <- .aom_cluster_summary(menu, choice, cluster, legacy$sumStats)
    results <- .aom_cluster_results(
      summary,
      legacy$constraints,
      preferences,
      method,
      alpha,
      nCritSimu,
      BARatio2MS,
      BARatio2UB,
      MNRatioGMS
    )
  } else {
    summary <- legacy$sumStats
    n_positive <- vapply(
      inequalities,
      function(values) sum(values > 0),
      integer(1L)
    )
    max_inequality <- vapply(inequalities, function(values) {
      if (length(values)) max(values) else NA_real_
    }, numeric(1L))
    methods <- if (method == "ALL") {
      c("GMS", "PI", "LF", "2MS", "2UB")
    } else {
      method
    }
    critical_fields <- c(
      GMS = "GMS", PI = "PI", LF = "LF", `2MS` = "MS", `2UB` = "UB"
    )
    p_value_fields <- c(GMS = "GMS", PI = "PI", LF = "LF")
    labels <- .aom_preference_labels(preferences)
    rows <- list()
    row_index <- 0L
    for (method_index in methods) {
      critical <- legacy$critVal[[critical_fields[[method_index]]]]
      available_alpha <- as.numeric(colnames(critical))
      p_values <- if (method_index %in% names(p_value_fields)) {
        as.numeric(legacy$pVal[[p_value_fields[[method_index]]]])
      } else {
        rep.int(NA_real_, n_preferences)
      }
      for (level in alpha) {
        column <- which.min(abs(available_alpha - level))
        critical_value <- as.numeric(critical[, column])
        row_index <- row_index + 1L
        rows[[row_index]] <- data.frame(
          preference_id = seq_len(n_preferences),
          preference = labels,
          method = method_index,
          alpha = level,
          statistic = as.numeric(legacy$Tstat),
          critical_value = critical_value,
          p_value = p_values,
          reject = as.numeric(legacy$Tstat) > critical_value,
          n_inequalities = counts,
          n_positive_sample_inequalities = n_positive,
          max_sample_inequality = max_inequality,
          stringsAsFactors = FALSE
        )
      }
    }
    results <- do.call(rbind, rows)
  }

  result <- list(
    results = results,
    preferences = preferences,
    inequalities = inequalities,
    summary = summary,
    constraints = legacy$constraints,
    options = list(
      method = method,
      alpha = alpha,
      nCritSimu = nCritSimu,
      BARatio2MS = BARatio2MS,
      BARatio2UB = BARatio2UB,
      MNRatioGMS = MNRatioGMS,
      attBinary = attBinary,
      limDataCorr = limDataCorr,
      sampling = if (is.null(cluster)) "iid" else "cluster",
      n_cluster = if (is.null(cluster)) NA_integer_ else summary$n_cluster
    ),
    elapsed = unname(proc.time()[["elapsed"]] - started),
    legacy = legacy
  )
  class(result) <- "ramchoiceRAMTest"
  result
}

#' @export
summary.ramchoiceAOMModel <- function(object, ...) {
  object$results
}

#' @export
print.ramchoiceAOMModel <- function(x, ...) {
  cat("\nPopulation analysis for homogeneous AOM\n\n")
  print(x$results, row.names = FALSE)
  invisible(x)
}

#' @export
summary.ramchoiceAOMTest <- function(object, ...) {
  object$results
}

#' @export
print.ramchoiceAOMTest <- function(x, ...) {
  cat("\nSample inference for homogeneous AOM\n\n")
  cat("Observations:", sum(x$summary$sumN), "\n")
  cat("Observed menus:", nrow(x$summary$sumMenu), "\n")
  cat("Elapsed seconds:", format(round(x$elapsed, 3), nsmall = 3), "\n\n")
  print(x$results, row.names = FALSE)
  invisible(x)
}

#' @export
summary.ramchoiceRAMTest <- function(object, ...) {
  object$results
}

#' @export
print.ramchoiceRAMTest <- function(x, ...) {
  cat("\nSample inference for RAM\n\n")
  cat("Observations:", sum(x$summary$sumN), "\n")
  cat("Observed menus:", nrow(x$summary$sumMenu), "\n")
  cat("Elapsed seconds:", format(round(x$elapsed, 3), nsmall = 3), "\n\n")
  print(x$results, row.names = FALSE)
  invisible(x)
}
