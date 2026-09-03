make_genmat_cache_fixture <- function() {
  sum_menu <- rbind(
    c(1, 1, 1),
    c(1, 1, 0),
    c(1, 0, 1),
    c(0, 1, 1)
  )
  list(sumMenu = sum_menu, sumMsize = rowSums(sum_menu))
}

clear_genmat_cache <- function() {
  cache <- getFromNamespace(".genMatCache", "ramchoice")
  rm(list = ls(envir = cache, all.names = TRUE), envir = cache)
  cache
}

test_that("genMat reuses only an identical most-recent construction", {
  cache <- clear_genmat_cache()
  fixture <- make_genmat_cache_fixture()

  first <- genMat(
    fixture$sumMenu,
    fixture$sumMsize,
    pref_list = matrix(1:3, nrow = 1),
    RAM = FALSE,
    AOM = TRUE
  )
  first_key <- cache$key
  second <- genMat(
    fixture$sumMenu,
    fixture$sumMsize,
    pref_list = matrix(1:3, nrow = 1),
    RAM = FALSE,
    AOM = TRUE
  )

  expect_identical(second, first)
  expect_identical(cache$key, first_key)

  changed <- genMat(
    fixture$sumMenu,
    fixture$sumMsize,
    pref_list = matrix(3:1, nrow = 1),
    RAM = FALSE,
    AOM = TRUE
  )

  expect_false(identical(cache$key, first_key))
  expect_identical(cache$value, changed)
  expect_false(identical(changed$R, first$R))
})

test_that("genMat cache preserves copy-on-modify behavior", {
  clear_genmat_cache()
  fixture <- make_genmat_cache_fixture()
  arguments <- list(
    sumMenu = fixture$sumMenu,
    sumMsize = fixture$sumMsize,
    pref_list = matrix(1:3, nrow = 1),
    RAM = FALSE,
    AOM = TRUE
  )

  original <- do.call(genMat, arguments)
  modified <- do.call(genMat, arguments)
  modified$R[1, 1] <- modified$R[1, 1] + 1

  expect_identical(do.call(genMat, arguments), original)
})

test_that("genMat retains its legacy one-menu result without caching", {
  cache <- clear_genmat_cache()

  result <- genMat(matrix(c(1, 1), nrow = 1), 2)

  expect_identical(result, list(R = matrix(NA, nrow = 0, ncol = 2), constN = 0))
  expect_false(exists("key", envir = cache, inherits = FALSE))
})
