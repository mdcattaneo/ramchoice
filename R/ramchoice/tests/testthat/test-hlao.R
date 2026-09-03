hlao_test_rankings <- function(values) {
  if (length(values) == 1L) {
    return(matrix(values, nrow = 1L))
  }
  do.call(rbind, lapply(seq_along(values), function(index) {
    cbind(values[index], hlao_test_rankings(values[-index]))
  }))
}

hlao_test_key <- function(items) paste(items, collapse = "")

hlao_test_suffix_close <- function(menus) {
  result <- menus
  for (items in menus) {
    for (position in seq_along(items)) {
      suffix <- items[position:length(items)]
      result[[hlao_test_key(suffix)]] <- suffix
    }
  }
  codes <- vapply(result, function(items) sum(2^(items - 1L)), numeric(1L))
  result[order(codes)]
}

hlao_test_population <- function() {
  rankings <- hlao_test_rankings(1:4)
  seed_menus <- list(
    c(1), c(2), c(3), c(4), c(1, 2), c(1, 3),
    c(2, 4), c(1, 2, 4), c(1, 3, 4)
  )
  names(seed_menus) <- vapply(seed_menus, hlao_test_key, character(1L))
  menus <- hlao_test_suffix_close(seed_menus)
  continuation <- c(0.86, 0.73, 0, 0.64)
  raw_tau <- c(
    31, 17, 23, 41, 29, 13, 37, 19, 43, 11, 47, 7,
    5, 53, 3, 59, 2, 61, 67, 71, 73, 79, 83, 89
  )
  tau <- raw_tau / sum(raw_tau)
  menu <- prob <- matrix(0, nrow = length(menus), ncol = 4)
  outside <- numeric(length(menus))

  for (menu_index in seq_along(menus)) {
    items <- menus[[menu_index]]
    menu[menu_index, items] <- 1
    reach <- cumprod(continuation[items])
    masses <- c(
      1 - reach[1L],
      if (length(reach) > 1L) reach[-length(reach)] - reach[-1L] else numeric(0L),
      reach[length(reach)]
    )
    outside[menu_index] <- masses[1L]
    for (ranking_index in seq_len(nrow(rankings))) {
      ranking <- rankings[ranking_index, ]
      for (prefix_size in seq_along(items)) {
        prefix <- items[seq_len(prefix_size)]
        winner <- ranking[match(TRUE, ranking %in% prefix)]
        prob[menu_index, winner] <- prob[menu_index, winner] +
          tau[ranking_index] * masses[prefix_size + 1L]
      }
    }
  }
  event <- apply(rankings, 1L, function(ranking) {
    match(4L, ranking) < match(1L, ranking)
  })
  list(
    menu = menu,
    prob = prob,
    outside = outside,
    rankings = rankings,
    tau = tau,
    event = event
  )
}

hlao_test_count_sample <- function() {
  menu_rows <- choice_rows <- list()
  add_counts <- function(menu_row, outside, first, second = 0L) {
    n <- outside + first + second
    menu_piece <- matrix(rep(menu_row, n), nrow = n, byrow = TRUE)
    choice_piece <- matrix(0, nrow = n, ncol = 2)
    if (first > 0L) {
      choice_piece[outside + seq_len(first), which(menu_row == 1L)[1L]] <- 1
    }
    if (second > 0L) {
      choice_piece[outside + first + seq_len(second), which(menu_row == 1L)[2L]] <- 1
    }
    menu_rows[[length(menu_rows) + 1L]] <<- menu_piece
    choice_rows[[length(choice_rows) + 1L]] <<- choice_piece
  }
  add_counts(c(1, 0), 20L, 80L)
  add_counts(c(0, 1), 25L, 75L)
  add_counts(c(1, 1), 20L, 56L, 24L)
  list(menu = do.call(rbind, menu_rows), choice = do.call(rbind, choice_rows))
}

hlao_test_nopi_population <- function() {
  menu <- rbind(c(1, 0), c(1, 1))
  prob <- rbind(c(0.8, 0), c(0.5, 0.2))
  outside <- c(0.2, 0.3)
  rankings <- hlaoRankings(1:2)
  list(
    menu = menu,
    prob = prob,
    outside = outside,
    event = rankings[, 1L] == 2L,
    truth = 0.4
  )
}

hlao_test_population_sample <- function(population, n = 100L) {
  menu_rows <- choice_rows <- list()
  for (menu_index in seq_len(nrow(population$menu))) {
    probabilities <- c(
      population$outside[menu_index],
      population$prob[menu_index, ]
    )
    counts <- as.integer(round(n * probabilities))
    counts[1L] <- counts[1L] + n - sum(counts)
    menu_piece <- matrix(
      rep(population$menu[menu_index, ], n),
      nrow = n,
      byrow = TRUE
    )
    choice_piece <- matrix(0, nrow = n, ncol = ncol(population$menu))
    offset <- counts[1L]
    for (alternative in seq_len(ncol(population$menu))) {
      if (counts[alternative + 1L] > 0L) {
        rows <- offset + seq_len(counts[alternative + 1L])
        choice_piece[rows, alternative] <- 1
        offset <- offset + counts[alternative + 1L]
      }
    }
    menu_rows[[menu_index]] <- menu_piece
    choice_rows[[menu_index]] <- choice_piece
  }
  list(menu = do.call(rbind, menu_rows), choice = do.call(rbind, choice_rows))
}

hlao_test_exact_specification <- function(menu, prob, outside) {
  summary <- list(
    menu = menu,
    prob = prob,
    outside_prob = outside,
    sample_size = rep(Inf, nrow(menu))
  )
  domain <- ramchoice:::.hlao_domain(menu, seq_len(ncol(menu)))
  bands <- ramchoice:::.hlao_probability_bands(summary, alpha = 0.05)
  attention <- ramchoice:::.hlao_recover_attention(
    menu,
    outside,
    domain,
    tolerance = 1e-12
  )
  full_attention <- ramchoice:::.hlao_full_attention(
    menu,
    prob,
    domain,
    attention,
    tolerance = 1e-12
  )
  ramchoice:::.hlao_specification_diagnostics(
    summary,
    domain,
    bands,
    attention,
    full_attention,
    alpha = 0.05,
    tolerance = 1e-12
  )
}

hlao_test_bm_alternative <- function(delta = 0.2) {
  menus <- lapply(1:7, function(code) {
    which(as.logical(intToBits(code)[seq_len(3L)]))
  })
  keys <- vapply(menus, paste, character(1L), collapse = "-")
  menu <- full <- matrix(0, nrow = length(menus), ncol = 3L)
  masses <- vector("list", length(menus))
  outside <- rep(0.25, length(menus))
  for (menu_index in seq_along(menus)) {
    items <- menus[[menu_index]]
    menu[menu_index, items] <- 1
    reach <- cumprod(rep(0.75, length(items)))
    masses[[menu_index]] <- c(
      1 - reach[1L],
      if (length(reach) > 1L) {
        reach[-length(reach)] - reach[-1L]
      } else {
        numeric(0L)
      },
      reach[length(reach)]
    )
    full[menu_index, min(items)] <- 1
  }
  target <- match("1-2-3", keys)
  full[target, 1L] <- 1 - delta
  full[target, 2L] <- delta

  prob <- matrix(0, nrow = nrow(menu), ncol = ncol(menu))
  for (menu_index in seq_along(menus)) {
    items <- menus[[menu_index]]
    for (prefix_size in seq_along(items)) {
      prefix_index <- match(
        paste(items[seq_len(prefix_size)], collapse = "-"),
        keys
      )
      prob[menu_index, ] <- prob[menu_index, ] +
        masses[[menu_index]][prefix_size + 1L] * full[prefix_index, ]
    }
  }
  list(menu = menu, prob = prob, outside = outside)
}

test_that("hlaoRankings uses a complete deterministic ordering", {
  rankings <- hlaoRankings(1:4)
  expect_identical(dim(rankings), c(24L, 4L))
  expect_identical(rankings[1L, ], c(rank1 = 1L, rank2 = 2L, rank3 = 3L, rank4 = 4L))
  expect_identical(rankings[24L, ], c(rank1 = 4L, rank2 = 3L, rank3 = 2L, rank4 = 1L))
  expect_identical(sort(apply(rankings, 1L, paste, collapse = "")),
                   sort(unique(apply(rankings, 1L, paste, collapse = ""))))
})

test_that("hlaoModel reproduces the supplement incomplete-menu bounds", {
  population <- hlao_test_population()
  fit <- hlaoModel(
    population$menu,
    population$prob,
    outside_prob = population$outside,
    events = list(`4 above 1` = population$event),
    dependence = "both"
  )

  expect_s3_class(fit, "ramchoiceHLAOModel")
  expect_true(fit$attention_diagnostics$valid)
  expect_true(all(fit$compatibility$compatible))
  expect_equal(
    fit$bounds[fit$bounds$mode == "independent", c("lower", "upper")],
    data.frame(lower = 0.545171, upper = 0.779855),
    tolerance = 1e-5,
    ignore_attr = TRUE
  )
  expect_equal(
    fit$bounds[fit$bounds$mode == "robust", c("lower", "upper")],
    data.frame(lower = 0.219045, upper = 1),
    tolerance = 1e-5,
    ignore_attr = TRUE
  )
  expect_null(fit$full_attention)
  expect_true(any(!fit$pairwise$identified))
})

test_that("structured events reproduce enumerated H-LAO bounds", {
  population <- hlao_test_population()
  enumerated <- hlaoModel(
    population$menu,
    population$prob,
    outside_prob = population$outside,
    events = hlaoEvent(4, 1, name = "4 above 1"),
    dependence = "independent",
    algorithm = "enumerate"
  )

  expect_equal(enumerated$bounds$lower, 0.545171, tolerance = 1e-5)
  expect_equal(enumerated$bounds$upper, 0.779855, tolerance = 1e-5)
  expect_identical(colnames(enumerated$events), "4 above 1")
})

test_that("column generation matches exhaustive H-LAO optimization", {
  population <- hlao_test_population()
  event <- hlaoEvent(4, 1, name = "4 above 1")
  enumerated <- hlaoModel(
    population$menu,
    population$prob,
    outside_prob = population$outside,
    events = event,
    dependence = "independent",
    agreement = c(7, 10),
    algorithm = "enumerate"
  )
  generated <- hlaoModel(
    population$menu,
    population$prob,
    outside_prob = population$outside,
    events = event,
    dependence = "independent",
    agreement = c(7, 10),
    algorithm = "column_generation",
    max_rankings = 1L
  )

  expect_equal(generated$bounds$lower, enumerated$bounds$lower, tolerance = 1e-7)
  expect_equal(generated$bounds$upper, enumerated$bounds$upper, tolerance = 1e-7)
  expect_equal(
    generated$agreement[, c("lower", "upper")],
    enumerated$agreement[, c("lower", "upper")],
    tolerance = 1e-7,
    ignore_attr = TRUE
  )
  expect_identical(generated$computation$algorithm, "column-generation")
  expect_lt(generated$computation$n_columns, 24)
  expect_true(generated$computation$phase_certified)
  expect_true(generated$computation$all_endpoints_certified)
  expect_true(generated$computation$certified)
  expect_true(all(generated$bounds$lower_certified))
  expect_true(all(generated$bounds$upper_certified))
  expect_lte(
    max(
      generated$bounds$lower_primal_residual,
      generated$bounds$upper_primal_residual
    ),
    generated$computation$tolerance
  )
  expect_gte(
    min(
      generated$bounds$lower_reduced_cost,
      generated$bounds$upper_reduced_cost
    ),
    -generated$computation$tolerance
  )
  expect_lte(
    max(
      generated$bounds$lower_dual_residual,
      generated$bounds$upper_dual_residual
    ),
    generated$computation$tolerance
  )
  expect_lte(
    max(
      generated$bounds$lower_optimality_gap_bound,
      generated$bounds$upper_optimality_gap_bound
    ),
    generated$computation$tolerance
  )
})

test_that("automatic column generation operates beyond the enumeration cap", {
  menu <- prob <- matrix(0, nrow = 6L, ncol = 6L)
  outside <- rep(0.2, 6L)
  for (menu_index in seq_len(6L)) {
    items <- menu_index:6L
    menu[menu_index, items] <- 1
    prob[menu_index, items[1L]] <- 0.8
  }
  fit <- hlaoModel(
    menu,
    prob,
    outside_prob = outside,
    events = hlaoEvent(1, 2:6, name = "1 top ranked"),
    dependence = "independent",
    algorithm = "auto",
    max_rankings = 100L
  )

  expect_true(fit$compatibility$compatible)
  expect_identical(fit$computation$algorithm, "column-generation")
  expect_equal(fit$computation$n_possible_rankings, 720)
  expect_lt(fit$computation$n_columns, 720)
  expect_equal(fit$bounds$lower, 1, tolerance = 1e-8)
  expect_equal(fit$bounds$upper, 1, tolerance = 1e-8)
  expect_true(fit$computation$certified)
})

test_that("H-LAO solver statuses are classified explicitly", {
  expect_identical(ramchoice:::.hlao_solver_status_class(0L), "success")
  expect_identical(ramchoice:::.hlao_solver_status_class(2L), "infeasible")
  expect_identical(ramchoice:::.hlao_solver_status_class(5L), "solver-error")
  expect_identical(ramchoice:::.hlao_solver_status_class(NA_integer_), "not-run")
})

test_that("hlaoModel computes agreement under all robustness envelopes", {
  population <- hlao_test_population()
  fit <- hlaoModel(
    population$menu,
    population$prob,
    outside_prob = population$outside,
    dependence = "all",
    agreement = 10L
  )

  expect_identical(fit$agreement$mode, c("independent", "robust", "noPI"))
  expect_true(all(fit$agreement$lower <= fit$agreement$upper + 1e-10))
  expect_true(all(fit$agreement$lower >= -1e-10))
  expect_true(all(fit$agreement$upper <= 1 + 1e-10))
})

test_that("hlaoModel adds a sharp no-SPI sensitivity envelope", {
  population <- hlao_test_population()
  fit <- hlaoModel(
    population$menu,
    population$prob,
    outside_prob = population$outside,
    events = list(`4 above 1` = population$event),
    dependence = "all"
  )

  expect_identical(fit$compatibility$mode, c("independent", "robust", "noPI"))
  expect_true(all(fit$compatibility$compatible))
  bounds <- fit$bounds
  independent <- bounds[bounds$mode == "independent", ]
  robust <- bounds[bounds$mode == "robust", ]
  nopi <- bounds[bounds$mode == "noPI", ]
  expect_lte(nopi$lower, robust$lower + 1e-8)
  expect_gte(nopi$upper, robust$upper - 1e-8)
  expect_lte(robust$lower, independent$lower + 1e-8)
  expect_gte(robust$upper, independent$upper - 1e-8)
})

test_that("no-SPI population analysis permits arbitrary menu domains", {
  population <- hlao_test_nopi_population()
  fit <- hlaoModel(
    population$menu,
    population$prob,
    outside_prob = population$outside,
    events = list(`2 above 1` = population$event),
    dependence = "noPI"
  )

  expect_true(fit$compatibility$compatible)
  expect_lte(fit$bounds$lower, population$truth)
  expect_gte(fit$bounds$upper, population$truth)
  expect_equal(nrow(fit$attention), 0L)
  expect_null(fit$prefix_masses)
})

test_that("hlaoModel recovers pairwise and full-attention objects", {
  menu <- rbind(c(1, 0), c(0, 1), c(1, 1))
  prob <- rbind(c(0.8, 0), c(0, 0.75), c(0.56, 0.24))
  rankings <- hlaoRankings(1:2)
  fit <- hlaoModel(
    menu,
    prob,
    events = list(`2 above 1` = rankings[, 1L] == 2L),
    dependence = "both",
    agreement = TRUE
  )

  expect_true(all(fit$compatibility$compatible))
  expect_equal(fit$pairwise$reach, 0.6)
  expect_equal(fit$pairwise$share_later_preferred, 0.4)
  expect_equal(fit$bounds$lower, c(0.4, 0.24), tolerance = 1e-7)
  expect_equal(fit$bounds$upper, c(0.4, 0.64), tolerance = 1e-7)
  expect_false(is.null(fit$full_attention))
  expect_true(fit$block_marschak$valid)
  expect_equal(
    fit$full_attention$matrix[which(rowSums(menu) == 2L), ],
    c(0.6, 0.4),
    tolerance = 1e-8
  )
  independent_agreement <- fit$agreement[
    fit$agreement$mode == "independent", , drop = FALSE
  ]
  expect_equal(independent_agreement$lower, c(0.8, 0.75, 0.72))
  expect_equal(independent_agreement$upper, c(0.8, 0.75, 0.72))
  expect_true(all(independent_agreement$method == "full-attention-formula"))
})

test_that("hlaoModel detects invalid recovered attention", {
  menu <- rbind(c(1, 0), c(0, 1), c(1, 1))
  prob <- rbind(c(0.5, 0), c(0, 0.75), c(0.56, 0.24))
  fit <- hlaoModel(menu, prob)

  expect_false(fit$attention_diagnostics$valid)
  expect_gt(fit$attention_diagnostics$max_attention_overload_violation, 0)
  expect_false(fit$compatibility$compatible)
})

test_that("hlaoTest returns weak-reach and event projection intervals", {
  sample <- hlao_test_count_sample()
  rankings <- hlaoRankings(1:2)
  fit <- hlaoTest(
    sample$menu,
    sample$choice,
    events = list(`2 above 1` = rankings[, 1L] == 2L),
    alpha = 0.05
  )

  expect_s3_class(fit, "ramchoiceHLAOTest")
  expect_equal(fit$pairwise$estimate, 0.4)
  expect_equal(
    fit$attention$reach[fit$attention$menu_id == 3L &
                          fit$attention$alternative == 2L],
    0.6
  )
  expect_false(is.null(fit$full_attention))
  expect_true(fit$block_marschak$valid)
  expect_lte(fit$pairwise$lower, 0.4)
  expect_gte(fit$pairwise$upper, 0.4)
  expect_false(fit$pairwise$empty)
  expect_true(fit$projection$feasible)
  expect_lte(fit$event_intervals$lower, 0.4)
  expect_gte(fit$event_intervals$upper, 0.4)
  expect_true(all(c(
    "attention-overload", "full-attention-probability",
    "block-marschak", "omnibus"
  ) %in% fit$specification$restriction))
  expect_false(any(fit$specification$reject, na.rm = TRUE))
  studentized <- fit$pairwise_studentized_components
  expect_true(any(
    studentized$lower <= 0.4 & studentized$upper >= 0.4
  ))
  expect_true(fit$pairwise_studentized$studentizable)
})

test_that("studentized pairwise inversion is uninformative at degeneracy", {
  menu <- rbind(
    matrix(rep(c(1, 0), 100), ncol = 2, byrow = TRUE),
    matrix(rep(c(0, 1), 100), ncol = 2, byrow = TRUE),
    matrix(rep(c(1, 1), 100), ncol = 2, byrow = TRUE)
  )
  choice <- matrix(0, nrow = 300, ncol = 2)
  choice[21:100, 1L] <- 1
  choice[221:300, 1L] <- 1
  fit <- hlaoTest(menu, choice)

  expect_false(fit$pairwise_studentized$studentizable)
  expect_equal(fit$pairwise_studentized$lower, 0)
  expect_equal(fit$pairwise_studentized$upper, 1)
  expect_equal(fit$pairwise_studentized$width, 1)
})

test_that("no-SPI sample projection covers events without suffix closure", {
  population <- hlao_test_nopi_population()
  sample <- hlao_test_population_sample(population, n = 100L)
  fit <- hlaoNoPITest(
    sample$menu,
    sample$choice,
    events = list(`2 above 1` = population$event)
  )

  expect_s3_class(fit, "ramchoiceHLAONoPITest")
  expect_true(fit$projection$feasible)
  expect_lte(fit$intervals$lower, population$truth)
  expect_gte(fit$intervals$upper, population$truth)
  expect_identical(fit$intervals$mode, "noPI")
})

test_that("H-LAO diagnostics separate attention and Block-Marschak failures", {
  attention_menu <- rbind(c(1, 0), c(0, 1), c(1, 1))
  attention_prob <- rbind(c(0.7, 0), c(0, 0.75), c(0.6, 0.3))
  attention <- hlao_test_exact_specification(
    attention_menu,
    attention_prob,
    outside = c(0.3, 0.25, 0.1)
  )
  attention_summary <- attention$summary[
    attention$summary$restriction == "attention-overload",
  ]
  expect_true(attention_summary$reject)
  expect_equal(attention_summary$max_violation_lower, 0.2, tolerance = 1e-12)

  bm_population <- hlao_test_bm_alternative()
  bm <- hlao_test_exact_specification(
    bm_population$menu,
    bm_population$prob,
    bm_population$outside
  )
  expect_false(bm$summary$reject[
    bm$summary$restriction == "attention-overload"
  ])
  expect_true(bm$summary$reject[
    bm$summary$restriction == "block-marschak"
  ])
  expect_equal(
    bm$summary$max_violation_lower[
      bm$summary$restriction == "block-marschak"
    ],
    0.2,
    tolerance = 1e-12
  )
})

test_that("population event bounds remain sharp at universal zero reach", {
  menu <- rbind(c(1, 0), c(0, 1), c(1, 1))
  prob <- matrix(0, nrow = 3L, ncol = 2L)
  rankings <- hlaoRankings(1:2)
  fit <- hlaoModel(
    menu,
    prob,
    outside_prob = rep(1, 3L),
    events = list(`2 above 1` = rankings[, 1L] == 2L),
    dependence = "all"
  )

  expect_true(all(fit$compatibility$compatible))
  expect_equal(fit$bounds$lower, rep(0, 3L))
  expect_equal(fit$bounds$upper, rep(1, 3L))
})

test_that("hlaoTest remains defined when terminal reach is zero", {
  menu <- rbind(
    matrix(rep(c(1, 0), 100), ncol = 2, byrow = TRUE),
    matrix(rep(c(0, 1), 100), ncol = 2, byrow = TRUE),
    matrix(rep(c(1, 1), 100), ncol = 2, byrow = TRUE)
  )
  choice <- matrix(0, nrow = 300, ncol = 2)
  choice[21:100, 1L] <- 1
  choice[221:300, 1L] <- 1
  fit <- hlaoTest(menu, choice)

  expect_equal(fit$pairwise$estimate, NA_real_)
  expect_equal(fit$pairwise$lower, 0)
  expect_equal(fit$pairwise$upper, 1)
  expect_false(fit$pairwise$empty)
})

test_that("correlated Gaussian bands use covariance and tighten regular cells", {
  sample <- hlao_test_count_sample()
  rankings <- hlaoRankings(1:2)
  hoeffding <- hlaoTest(
    sample$menu,
    sample$choice,
    events = list(`2 above 1` = rankings[, 1L] == 2L)
  )
  set.seed(20260713)
  gaussian <- hlaoTest(
    sample$menu,
    sample$choice,
    events = list(`2 above 1` = rankings[, 1L] == 2L),
    band_method = "gaussian",
    n_band_draws = 500L
  )

  expect_true(is.finite(gaussian$options$band_critical_value))
  expect_gt(gaussian$options$n_gaussian_cells, 0L)
  expect_identical(
    gaussian$pairwise$method,
    "weak-reach-correlated-gaussian-exact-hybrid-projection"
  )
  expect_lte(
    gaussian$pairwise$upper - gaussian$pairwise$lower,
    hoeffding$pairwise$upper - hoeffding$pairwise$lower
  )
  expect_lte(
    gaussian$event_intervals$upper - gaussian$event_intervals$lower,
    hoeffding$event_intervals$upper - hoeffding$event_intervals$lower
  )
  expect_lte(gaussian$pairwise$lower, 0.4)
  expect_gte(gaussian$pairwise$upper, 0.4)
})

test_that("Gaussian inference keeps Hoeffding protection at zero reach", {
  menu <- rbind(
    matrix(rep(c(1, 0), 100), ncol = 2, byrow = TRUE),
    matrix(rep(c(0, 1), 100), ncol = 2, byrow = TRUE),
    matrix(rep(c(1, 1), 100), ncol = 2, byrow = TRUE)
  )
  choice <- matrix(0, nrow = 300, ncol = 2)
  choice[21:100, 1L] <- 1
  choice[221:300, 1L] <- 1
  set.seed(20260714)
  fit <- hlaoTest(
    menu,
    choice,
    band_method = "gaussian",
    n_band_draws = 500L
  )
  hoeffding <- hlaoTest(menu, choice)

  expect_gt(fit$options$n_fallback_cells, 0L)
  expect_identical(
    fit$options$band_method,
    "correlated-gaussian-exact-hybrid"
  )
  expect_equal(fit$options$gaussian_alpha, 0.025)
  expect_equal(fit$options$fallback_alpha, 0.025)
  expect_equal(fit$pairwise$lower, 0)
  expect_equal(fit$pairwise$upper, 1)
  fallback <- merge(
    fit$bands[!fit$bands$gaussian_active, c("menu_id", "outcome", "lower", "upper")],
    hoeffding$bands[, c("menu_id", "outcome", "lower", "upper")],
    by = c("menu_id", "outcome"),
    suffixes = c("_exact", "_hoeffding")
  )
  expect_true(any(
    fallback$upper_exact - fallback$lower_exact <
      fallback$upper_hoeffding - fallback$lower_hoeffding
  ))
})

test_that("direct delta diagnostics use one joint Gaussian critical value", {
  sample <- hlao_test_count_sample()
  set.seed(20260715)
  fit <- hlaoTest(
    sample$menu,
    sample$choice,
    diagnostic_method = "delta",
    n_band_draws = 500L
  )

  expect_identical(fit$options$diagnostic_method, "delta")
  expect_true(is.finite(fit$options$diagnostic_critical_value))
  expect_gt(fit$options$n_active_diagnostics, 0L)
  expect_true(all(fit$specification$method == "direct-delta-gaussian"))
  expect_true(all(fit$outer_specification$method == "hoeffding"))
  expect_false(any(fit$specification$reject))
  expect_true(all(is.finite(fit$specification$p_value)))
  expect_true(all(
    fit$specification$p_value >= 0 & fit$specification$p_value <= 1
  ))
  expect_equal(
    fit$specification$p_value[fit$specification$restriction == "omnibus"],
    min(fit$specification_details$p_value)
  )
})

test_that("direct diagnostic derivatives match finite differences", {
  sample <- hlao_test_count_sample()
  summary <- ramchoice:::.hlao_sample_summary(
    sample$menu,
    sample$choice,
    NULL
  )
  domain <- ramchoice:::.hlao_domain(summary$menu, 1:2)
  recover <- function(current) {
    attention <- ramchoice:::.hlao_recover_attention(
      current$menu,
      current$outside_prob,
      domain,
      tolerance = 1e-12
    )
    full_attention <- ramchoice:::.hlao_full_attention(
      current$menu,
      current$prob,
      domain,
      attention,
      tolerance = 1e-12
    )
    ramchoice:::.hlao_delta_moments(
      current,
      domain,
      attention,
      full_attention,
      tolerance = 1e-12
    )
  }
  moments <- recover(summary)
  epsilon <- 1e-6
  plus <- minus <- summary
  plus$outside_prob[3L] <- plus$outside_prob[3L] - epsilon
  plus$prob[3L, 2L] <- plus$prob[3L, 2L] + epsilon
  minus$outside_prob[3L] <- minus$outside_prob[3L] + epsilon
  minus$prob[3L, 2L] <- minus$prob[3L, 2L] - epsilon
  numerical <- (recover(plus)$estimate - recover(minus)$estimate) /
    (2 * epsilon)
  analytic <- moments$gradient[, 7L] - moments$gradient[, 5L]

  expect_equal(analytic, numerical, tolerance = 1e-6)
})

test_that("direct diagnostics require identified full-attention recovery", {
  menu <- rbind(
    matrix(rep(c(1, 0), 100), ncol = 2, byrow = TRUE),
    matrix(rep(c(0, 1), 100), ncol = 2, byrow = TRUE),
    matrix(rep(c(1, 1), 100), ncol = 2, byrow = TRUE)
  )
  choice <- matrix(0, nrow = 300, ncol = 2)
  choice[21:100, 1L] <- 1
  choice[221:300, 1L] <- 1

  expect_error(
    hlaoTest(menu, choice, diagnostic_method = "delta"),
    "complete menu domain and positive terminal reach"
  )
})

test_that("H-LAO inference supports subject-level clustering", {
  sample <- hlao_test_count_sample()
  set.seed(720)
  cluster <- sample(rep(seq_len(20L), length.out = nrow(sample$menu)))

  set.seed(721)
  iid <- hlaoTest(
    sample$menu,
    sample$choice,
    list_order = 1:2,
    band_method = "gaussian",
    diagnostic_method = "delta",
    n_band_draws = 200L
  )
  set.seed(721)
  clustered <- hlaoTest(
    sample$menu,
    sample$choice,
    list_order = 1:2,
    band_method = "gaussian",
    diagnostic_method = "delta",
    n_band_draws = 200L,
    cluster = cluster
  )

  expect_identical(clustered$options$sampling, "cluster")
  expect_identical(clustered$options$n_cluster, 20L)
  expect_match(clustered$options$band_method, "cluster-multiplier")
  expect_true(all(grepl("cluster", clustered$specification$method)))
  expect_equal(clustered$pairwise$estimate, iid$pairwise$estimate)
  expect_equal(
    clustered$summary$covariance,
    t(clustered$summary$covariance)
  )
  expect_gt(max(abs(
    clustered$summary$covariance[
      clustered$summary$outside_cell[1L],
      clustered$summary$outside_cell[2L]
    ]
  )), 0)
})

test_that("cluster identifiers are validated", {
  sample <- hlao_test_count_sample()
  expect_error(
    hlaoTest(sample$menu, sample$choice, cluster = rep(1L, nrow(sample$menu))),
    "at least two clusters"
  )
  expect_error(
    hlaoTest(sample$menu, sample$choice, cluster = seq_len(2L)),
    "one nonmissing identifier"
  )
})

test_that("H-LAO interfaces enforce suffix closure", {
  expect_error(
    hlaoModel(
      rbind(c(1, 0), c(1, 1)),
      rbind(c(0.8, 0), c(0.56, 0.24))
    ),
    "suffix closed"
  )
})

test_that("H-LAO Gaussian options are validated", {
  sample <- hlao_test_count_sample()
  expect_error(
    hlaoTest(
      sample$menu,
      sample$choice,
      band_method = "gaussian",
      n_band_draws = 99L
    ),
    "at least 100"
  )
  expect_error(
    hlaoTest(
      sample$menu,
      sample$choice,
      band_method = "gaussian",
      boundary_count = -1L
    ),
    "nonnegative"
  )
})
