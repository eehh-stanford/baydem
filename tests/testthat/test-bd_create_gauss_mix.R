test_that("non-truncated works", {
  testthat::expect_s4_class(
    bd_create_gauss_mix(x = tibble::tibble(
      pi = c(0.2, 0.8),
      mu = c(775, 1000),
      sig = c(35, 45)
    )),
    "AbscontDistribution"
  )
})

test_that("truncated works", {
  testthat::expect_s4_class(
    bd_create_gauss_mix(
      x = tibble::tibble(
        pi = c(0.2, 0.8),
        mu = c(775, 1000),
        sig = c(35, 45)
      ),
      taumin = 600,
      taumax = 1280
    ),
    "AbscontDistribution"
  )
})

test_that("big mixtures work", {
  n <- 100
  testthat::expect_s4_class(
    bd_create_gauss_mix(x = tibble::tibble(
      pi = runif(n) / n,
      mu = runif(n, 1, 1000),
      sig = runif(n, 20, 120)
    )),
    "AbscontDistribution"
  )
})
