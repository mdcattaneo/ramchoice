make_v2_fixture_data <- function() {
  set.seed(101)
  simulated <- lapply(4:2, function(size) {
    logitSimu(n = 20, uSize = 4, mSize = size, a = 2)
  })

  list(
    menu = do.call(rbind, lapply(simulated, `[[`, "menu")),
    choice = do.call(rbind, lapply(simulated, `[[`, "choice"))
  )
}

test_that("version 2.2 exports and signatures remain available", {
  version_2_exports <- c(
    "genMat", "logitAtte", "logitSimu", "rAtte", "revealAtte",
    "revealPref", "revealPrefModel", "sumData"
  )
  expect_true(
    all(version_2_exports %in% getNamespaceExports("ramchoice")),
    info = "All version 2.2 exports must remain available"
  )

  expect_identical(names(formals(sumData)), c("menu", "choice"))
  expect_identical(
    names(formals(genMat)),
    c(
      "sumMenu", "sumMsize", "pref_list", "RAM", "AOM",
      "limDataCorr", "attBinary"
    )
  )
  expect_identical(names(formals(logitAtte)), c("mSize", "a"))
  expect_identical(
    names(formals(logitSimu)),
    c("n", "uSize", "mSize", "a")
  )
  expect_identical(formals(rAtte), formals(revealPref))
  expect_identical(rAtte, revealPref)
  expect_identical(
    names(formals(revealAtte)),
    c(
      "menu", "choice", "alternative", "S", "lower", "upper", "pref",
      "nCritSimu", "level"
    )
  )
  expect_identical(
    names(formals(revealPrefModel)),
    c(
      "menu", "prob", "pref_list", "RAM", "AOM", "limDataCorr",
      "attBinary"
    )
  )
})

test_that("version 2.2 S3 methods remain registered", {
  classes <- c(
    "ramchoiceRevealAtte",
    "ramchoiceRevealPref",
    "ramchoiceRevealPrefModel"
  )

  for (class in classes) {
    expect_true(is.function(getS3method("print", class)))
    expect_true(is.function(getS3method("summary", class)))
  }
})

test_that("logit attention and simulation retain seeded outputs", {
  expect_equal(
    logitAtte(3, 2),
    list(
      choiceProb = c(0.75, 0.208333333333333, 0.0416666666666667),
      atteFreq = 0.75
    ),
    tolerance = 1e-14
  )

  set.seed(42)
  simulated <- logitSimu(n = 2, uSize = 4, mSize = 3, a = 2)
  expected_menu <- matrix(
    c(
      1, 1, 1, 1, 1, 1, 0, 0,
      1, 1, 1, 1, 0, 0, 1, 1,
      1, 1, 0, 0, 1, 1, 1, 1,
      0, 0, 1, 1, 1, 1, 1, 1
    ),
    nrow = 8,
    ncol = 4
  )
  expected_choice <- matrix(
    c(
      0, 1, 0, 1, 1, 1, 0, 0,
      0, 0, 1, 0, 0, 0, 1, 1,
      1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0
    ),
    nrow = 8,
    ncol = 4
  )

  expect_identical(simulated$menu, expected_menu)
  expect_identical(simulated$choice, expected_choice)

  summary <- sumData(simulated$menu, simulated$choice)
  expect_identical(
    names(summary),
    c("sumMenu", "sumProb", "sumN", "sumMsize", "sumProbVec", "Sigma")
  )
  expect_identical(summary$sumN, rep(2L, 4))
  expect_equal(summary$sumMsize, rep(3, 4))
  expect_identical(dim(summary$sumProbVec), c(12L, 1L))
  expect_identical(dim(summary$Sigma), c(12L, 12L))
})

test_that("preference inference retains version 2.2 contracts and values", {
  data <- make_v2_fixture_data()
  preferences <- rbind(1:4, 4:1)
  summary <- sumData(data$menu, data$choice)
  constraints <- genMat(
    summary$sumMenu,
    summary$sumMsize,
    pref_list = preferences,
    RAM = TRUE,
    AOM = FALSE,
    limDataCorr = FALSE
  )

  expect_identical(names(constraints), c("R", "ConstN"))
  expect_identical(dim(constraints$R), c(36L, 28L))
  expect_identical(constraints$ConstN, c(18L, 18L))

  set.seed(102)
  result <- revealPref(
    data$menu,
    data$choice,
    pref_list = preferences,
    method = "GMS",
    nCritSimu = 200,
    RAM = TRUE,
    AOM = FALSE,
    limDataCorr = FALSE
  )

  expect_s3_class(result, "ramchoiceRevealPref")
  expect_identical(
    names(result),
    c("sumStats", "constraints", "Tstat", "critVal", "pVal", "pref", "method", "opt")
  )
  expect_identical(names(result$critVal), c("GMS", "PI", "LF", "MS", "UB"))
  expect_identical(names(result$pVal), c("GMS", "PI", "LF"))
  expect_equal(
    result$Tstat,
    c(1.0690449676497, 1.87867287325545),
    tolerance = 1e-12
  )
  gms_critical_values <- unname(result$critVal$GMS)
  gms_p_values <- unname(result$pVal$GMS)
  expect_identical(dim(gms_critical_values), c(2L, 3L))
  expect_true(all(is.finite(gms_critical_values)))
  expect_true(all(gms_critical_values >= 0))
  expect_true(all(apply(gms_critical_values, 1, diff) > 0))
  expect_identical(dim(gms_p_values), c(2L, 1L))
  expect_true(all(is.finite(gms_p_values)))
  expect_true(all(gms_p_values >= 0 & gms_p_values <= 1))

  set.seed(102)
  alias_result <- rAtte(
    data$menu,
    data$choice,
    pref_list = preferences,
    method = "GMS",
    nCritSimu = 200,
    RAM = TRUE,
    AOM = FALSE,
    limDataCorr = FALSE
  )
  expect_identical(alias_result, result)
})

test_that("population and attention results retain version 2.2 contracts", {
  data <- make_v2_fixture_data()
  preferences <- rbind(1:4, 4:1)
  summary <- sumData(data$menu, data$choice)

  population <- revealPrefModel(
    summary$sumMenu,
    summary$sumProb,
    pref_list = preferences,
    RAM = TRUE,
    AOM = FALSE,
    limDataCorr = FALSE
  )
  expect_s3_class(population, "ramchoiceRevealPrefModel")
  expect_identical(
    names(population),
    c("constraints", "probVec", "inequalities", "pref", "opt", "sumStats")
  )
  expect_identical(names(population$inequalities), c("R", "ConstN"))
  expect_identical(population$inequalities$ConstN, c(18L, 18L))

  set.seed(103)
  attention <- revealAtte(
    data$menu,
    data$choice,
    alternative = 1,
    S = matrix(c(1, 1, 0, 0), nrow = 1),
    lower = TRUE,
    upper = TRUE,
    pref = matrix(1:4, nrow = 1),
    nCritSimu = 200,
    level = 0.95
  )
  expect_s3_class(attention, "ramchoiceRevealAtte")
  expect_identical(
    names(attention),
    c("sumStats", "lowerBound", "upperBound", "critVal", "opt")
  )
  expect_true(all(is.finite(c(
    attention$lowerBound,
    attention$upperBound,
    attention$critVal
  ))))
  expect_gt(unname(attention$critVal), 0)
  expect_lt(unname(attention$lowerBound), unname(attention$upperBound))
  expect_equal(
    (attention$lowerBound + attention$upperBound) / 2,
    matrix(0.8),
    tolerance = 1e-12
  )

  set.seed(103)
  repeated_attention <- revealAtte(
    data$menu,
    data$choice,
    alternative = 1,
    S = matrix(c(1, 1, 0, 0), nrow = 1),
    lower = TRUE,
    upper = TRUE,
    pref = matrix(1:4, nrow = 1),
    nCritSimu = 200,
    level = 0.95
  )
  expect_identical(repeated_attention, attention)
})
