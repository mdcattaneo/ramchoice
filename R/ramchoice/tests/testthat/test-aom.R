make_aom_sample_fixture <- function() {
  set.seed(101)
  simulated <- lapply(4:2, function(size) {
    logitSimu(n = 20, uSize = 4, mSize = size, a = 2)
  })
  list(
    menu = do.call(rbind, lapply(simulated, `[[`, "menu")),
    choice = do.call(rbind, lapply(simulated, `[[`, "choice"))
  )
}

make_aom_population_fixture <- function() {
  pieces <- lapply(4:2, function(size) {
    menus <- t(utils::combn(4, size))
    probabilities <- matrix(0, nrow = nrow(menus), ncol = 4)
    choice_probability <- logitAtte(size, 2)$choiceProb
    for (index in seq_len(nrow(menus))) {
      probabilities[index, menus[index, ]] <- choice_probability
    }
    list(menu = {
      result <- matrix(0, nrow = nrow(menus), ncol = 4)
      for (index in seq_len(nrow(menus))) {
        result[index, menus[index, ]] <- 1
      }
      result
    }, prob = probabilities)
  })
  list(
    menu = do.call(rbind, lapply(pieces, `[[`, "menu")),
    prob = do.call(rbind, lapply(pieces, `[[`, "prob"))
  )
}

test_that("aomTest preserves the legacy AOM calculation", {
  data <- make_aom_sample_fixture()
  preferences <- rbind(1:4, 4:1)

  set.seed(102)
  result <- aomTest(
    data$menu,
    data$choice,
    pref_list = preferences,
    method = "GMS",
    alpha = 0.05,
    nCritSimu = 200
  )
  set.seed(102)
  legacy <- revealPref(
    data$menu,
    data$choice,
    pref_list = preferences,
    method = "GMS",
    nCritSimu = 200,
    RAM = FALSE,
    AOM = TRUE,
    limDataCorr = TRUE
  )

  expect_s3_class(result, "ramchoiceAOMTest")
  expect_identical(result$legacy, legacy)
  expect_identical(
    names(result$results),
    c(
      "preference_id", "preference", "method", "alpha", "statistic",
      "critical_value", "p_value", "reject", "n_inequalities",
      "n_positive_sample_inequalities", "max_sample_inequality"
    )
  )
  expect_identical(result$results$preference_id, 1:2)
  expect_identical(result$results$method, rep("GMS", 2))
  expect_equal(result$results$alpha, rep(0.05, 2))
  expect_equal(result$results$statistic, legacy$Tstat)
  expect_equal(result$results$critical_value, unname(legacy$critVal$GMS[, 2]))
  expect_equal(result$results$p_value, as.numeric(legacy$pVal$GMS))
  expect_identical(result$results$reject, c(FALSE, FALSE))
  expect_identical(result$results$n_inequalities, c(26L, 26L))
  expect_equal(summary(result), result$results)
})

test_that("aomTest supports cluster-robust multiplier inference", {
  data <- make_aom_sample_fixture()
  cluster <- rep(seq_len(20L), length.out = nrow(data$menu))
  preferences <- rbind(1:4, 4:1)

  set.seed(1021)
  iid <- aomTest(
    data$menu,
    data$choice,
    pref_list = preferences,
    method = "LF",
    nCritSimu = 200L
  )
  set.seed(1021)
  clustered <- aomTest(
    data$menu,
    data$choice,
    pref_list = preferences,
    method = "LF",
    nCritSimu = 200L,
    cluster = cluster
  )

  expect_identical(clustered$options$sampling, "cluster")
  expect_identical(clustered$options$n_cluster, 20L)
  expect_equal(clustered$inequalities, iid$inequalities)
  expect_equal(clustered$summary$covariance, t(clustered$summary$covariance))
  expect_true(all(is.finite(clustered$results$critical_value)))
  expect_true(all(clustered$results$p_value >= 0 &
                    clustered$results$p_value <= 1))
})

test_that("ramTest preserves RAM restrictions under iid and clustered sampling", {
  data <- make_aom_sample_fixture()
  preferences <- rbind(1:4, 4:1)
  cluster <- rep(seq_len(20L), length.out = nrow(data$menu))

  set.seed(1022)
  iid <- ramTest(
    data$menu,
    data$choice,
    pref_list = preferences,
    method = "GMS",
    nCritSimu = 200L
  )
  set.seed(1022)
  legacy <- revealPref(
    data$menu,
    data$choice,
    pref_list = preferences,
    method = "GMS",
    nCritSimu = 200L,
    RAM = TRUE,
    AOM = FALSE,
    limDataCorr = TRUE
  )
  set.seed(1023)
  clustered <- ramTest(
    data$menu,
    data$choice,
    pref_list = preferences,
    method = "LF",
    nCritSimu = 200L,
    cluster = cluster
  )

  expect_s3_class(iid, "ramchoiceRAMTest")
  expect_equal(iid$results$statistic, legacy$Tstat)
  expect_equal(iid$constraints, legacy$constraints)
  expect_identical(clustered$options$sampling, "cluster")
  expect_true(all(is.finite(clustered$results$critical_value)))
})

test_that("aomTest validates new interface options", {
  data <- make_aom_sample_fixture()
  expect_error(aomTest(data$menu, data$choice, alpha = 0.025), "0.10")
  expect_error(aomTest(data$menu, data$choice, method = 1), "character")
  expect_error(aomTest(data$menu, data$choice, method = "bootstrap"), "must be")
})

test_that("aomTest reports the all-inequality baseline alongside GMS", {
  data <- make_aom_sample_fixture()
  preferences <- rbind(1:4, 4:1)

  set.seed(103)
  result <- aomTest(
    data$menu,
    data$choice,
    pref_list = preferences,
    method = "ALL",
    alpha = 0.05,
    nCritSimu = 50
  )

  expect_identical(
    unique(result$results$method),
    c("GMS", "PI", "LF", "2MS", "2UB")
  )
  expect_equal(nrow(result$results), 5L * nrow(preferences))
  expect_equal(
    result$results$critical_value[result$results$method == "LF"],
    unname(result$legacy$critVal$LF[, 2])
  )

  set.seed(104)
  implicit <- aomTest(
    data$menu,
    data$choice,
    pref_list = preferences,
    method = "GMS",
    alpha = 0.05,
    nCritSimu = 50
  )
  set.seed(104)
  explicit <- aomTest(
    data$menu,
    data$choice,
    pref_list = preferences,
    method = "GMS",
    alpha = 0.05,
    nCritSimu = 50,
    MNRatioGMS = 1 / log(nrow(data$menu))
  )
  expect_equal(implicit$results, explicit$results)
})

test_that("two-step moment selection handles one effective inequality", {
  set.seed(28)
  simulated <- lapply(c(6L, 2L), function(size) {
    logitSimu(n = 200L, uSize = 6L, mSize = size, a = 2)
  })
  menu <- do.call(rbind, lapply(simulated, `[[`, "menu"))
  choice <- do.call(rbind, lapply(simulated, `[[`, "choice"))
  preferences <- matrix(
    c(
      1, 2, 3, 4, 5, 6,
      2, 3, 4, 5, 6, 1,
      1, 2, 6, 5, 4, 3,
      1, 6, 5, 4, 3, 2
    ),
    ncol = 6,
    byrow = TRUE
  )

  set.seed(10028)
  result <- revealPref(
    menu,
    choice,
    pref_list = preferences,
    method = "2MS",
    nCritSimu = 100L,
    RAM = FALSE,
    AOM = TRUE,
    limDataCorr = TRUE
  )

  expect_s3_class(result, "ramchoiceRevealPref")
  expect_true(all(is.finite(result$critVal$MS)))
})

test_that("aomModel reports population compatibility and violations", {
  population <- make_aom_population_fixture()
  preferences <- rbind(1:4, 4:1)
  result <- aomModel(
    population$menu,
    population$prob,
    pref_list = preferences
  )
  legacy <- revealPrefModel(
    population$menu,
    population$prob,
    pref_list = preferences,
    RAM = FALSE,
    AOM = TRUE,
    limDataCorr = TRUE
  )

  expect_s3_class(result, "ramchoiceAOMModel")
  expect_identical(result$legacy, legacy)
  expect_identical(
    names(result$results),
    c(
      "preference_id", "preference", "compatible", "n_inequalities",
      "n_violated", "max_inequality", "max_violation"
    )
  )
  expect_identical(result$results$n_inequalities, c(26L, 26L))
  expect_true(result$results$compatible[1])
  expect_false(result$results$compatible[2])
  expect_identical(result$results$n_violated[1], 0L)
  expect_gt(result$results$n_violated[2], 0L)
  expect_lte(result$results$max_violation[1], result$tolerance)
  expect_gt(result$results$max_violation[2], result$tolerance)
  expect_equal(summary(result), result$results)
})

test_that("aomIdentify matches exhaustive ranking enumeration", {
  population <- make_aom_population_fixture()
  rankings <- hlaoRankings(1:4)
  exhaustive <- aomModel(
    population$menu,
    population$prob,
    pref_list = rankings
  )
  identified <- aomIdentify(population$menu, population$prob)
  compatible <- rankings[exhaustive$results$compatible, , drop = FALSE]

  expect_s3_class(identified, "ramchoiceAOMIdentification")
  expect_true(identified$compatible)
  expect_true(any(apply(compatible, 1L, function(ranking) {
    all(ranking == identified$preference)
  })))
  expect_identical(identified$diagnostics$n_binary_variables, 12L)
  expect_identical(identified$diagnostics$n_milp_solves, 7L)

  for (row_index in seq_len(nrow(identified$pairwise))) {
    row <- identified$pairwise[row_index, ]
    a_above_b <- apply(compatible, 1L, function(ranking) {
      match(row$alternative_a, ranking) < match(row$alternative_b, ranking)
    })
    expect_identical(row$a_preferred_possible, any(a_above_b))
    expect_identical(row$b_preferred_possible, any(!a_above_b))
  }
})

test_that("aomIdentify works on incomplete menu domains", {
  population <- make_aom_population_fixture()
  keep <- rowSums(population$menu) >= 3L
  rankings <- hlaoRankings(1:4)
  exhaustive <- aomModel(
    population$menu[keep, , drop = FALSE],
    population$prob[keep, , drop = FALSE],
    pref_list = rankings
  )
  identified <- aomIdentify(
    population$menu[keep, , drop = FALSE],
    population$prob[keep, , drop = FALSE]
  )

  expect_true(identified$compatible)
  for (row_index in seq_len(nrow(identified$pairwise))) {
    row <- identified$pairwise[row_index, ]
    compatible <- rankings[exhaustive$results$compatible, , drop = FALSE]
    a_above_b <- apply(compatible, 1L, function(ranking) {
      match(row$alternative_a, ranking) < match(row$alternative_b, ranking)
    })
    expect_identical(row$a_preferred_possible, any(a_above_b))
    expect_identical(row$b_preferred_possible, any(!a_above_b))
  }
})

test_that("aomIdentify validates adding up and duplicate menus", {
  population <- make_aom_population_fixture()
  invalid <- population$prob
  invalid[1L, which(population$menu[1L, ] == 1L)[1L]] <-
    invalid[1L, which(population$menu[1L, ] == 1L)[1L]] - 0.05
  expect_error(
    aomIdentify(population$menu, invalid),
    "sum to one"
  )

  duplicate_menu <- rbind(population$menu, population$menu[1L, ])
  duplicate_prob <- rbind(population$prob, population$prob[1L, ])
  collapsed <- aomIdentify(duplicate_menu, duplicate_prob, pairwise = FALSE)
  baseline <- aomIdentify(population$menu, population$prob, pairwise = FALSE)
  expect_equal(collapsed$menu, baseline$menu)
  expect_identical(collapsed$compatible, baseline$compatible)

  duplicate_prob[nrow(duplicate_prob), which(duplicate_menu[nrow(duplicate_menu), ] == 1L)[1:2]] <-
    duplicate_prob[nrow(duplicate_prob), which(duplicate_menu[nrow(duplicate_menu), ] == 1L)[1:2]] + c(0.01, -0.01)
  expect_error(
    aomIdentify(duplicate_menu, duplicate_prob),
    "conflicting"
  )
})

test_that("aomIdentify distinguishes infeasibility from solver errors", {
  expect_identical(ramchoice:::.aom_status_class(0L), "success")
  expect_identical(ramchoice:::.aom_status_class(2L), "infeasible")
  expect_identical(ramchoice:::.aom_status_class(5L), "solver-error")
})
