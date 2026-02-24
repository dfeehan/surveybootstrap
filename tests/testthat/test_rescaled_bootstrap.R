
library(surveybootstrap)

set.seed(12345)
survey <- MU284.complex.surveys[[1]]

# ── rescaled.bootstrap.sample ─────────────────────────────────────────────────

context("rescaled.bootstrap.sample - structure")

test_that("returns a list of length num.reps", {
  res <- rescaled.bootstrap.sample(survey.data   = survey,
                                   survey.design = ~ CL,
                                   num.reps      = 5)
  expect_length(res, 5)
})

test_that("each replicate has 'index' and 'weight.scale' columns", {
  res <- rescaled.bootstrap.sample(survey.data   = survey,
                                   survey.design = ~ CL,
                                   num.reps      = 3)
  for (rep in res) {
    expect_true("index"        %in% names(rep))
    expect_true("weight.scale" %in% names(rep))
  }
})

test_that("index column is a plain numeric/integer vector (not a nested tibble)", {
  res <- rescaled.bootstrap.sample(survey.data   = survey,
                                   survey.design = ~ CL,
                                   num.reps      = 3)
  for (rep in res) {
    expect_true(is.numeric(rep$index) || is.integer(rep$index))
    expect_false(is.data.frame(rep$index))
  }
})

test_that("weight.scale values are non-negative", {
  res <- rescaled.bootstrap.sample(survey.data   = survey,
                                   survey.design = ~ CL,
                                   num.reps      = 10)
  for (rep in res) {
    expect_true(all(rep$weight.scale >= 0))
  }
})

test_that("each replicate has the same number of rows as the survey", {
  res <- rescaled.bootstrap.sample(survey.data   = survey,
                                   survey.design = ~ CL,
                                   num.reps      = 3)
  for (rep in res) {
    expect_equal(nrow(rep), nrow(survey))
  }
})

# ── Regression: weight.scale column accessible (no tibble $ warning) ──────────

context("rescaled.bootstrap.sample - regression: weight.scale column")

test_that("accessing $weight.scale on each replicate emits no warnings", {
  res <- rescaled.bootstrap.sample(survey.data   = survey,
                                   survey.design = ~ CL,
                                   num.reps      = 5)
  for (rep in res) {
    expect_no_warning(rep$weight.scale)
  }
})

# ── bootstrap.estimates with rescaled.bootstrap.sample ───────────────────────

context("bootstrap.estimates + rescaled.bootstrap.sample")

test_that("produces no warnings about weight.scale", {
  estimator <- function(survey.data, weights) {
    data.frame(TS82 = sum(survey.data$S82 * weights))
  }
  expect_no_warning(
    bootstrap.estimates(
      survey.data   = survey,
      survey.design = ~ CL,
      num.reps      = 20,
      estimator.fn  = estimator,
      weights       = "sample_weight",
      bootstrap.fn  = "rescaled.bootstrap.sample"
    )
  )
})

test_that("returns finite estimates (regression against NaN/Inf from naming bug)", {
  estimator <- function(survey.data, weights) {
    data.frame(TS82 = sum(survey.data$S82 * weights))
  }
  res <- bootstrap.estimates(
    survey.data   = survey,
    survey.design = ~ CL,
    num.reps      = 20,
    estimator.fn  = estimator,
    weights       = "sample_weight",
    bootstrap.fn  = "rescaled.bootstrap.sample"
  )
  estimates <- sapply(res, `[[`, "TS82")
  expect_true(all(is.finite(estimates)),
              info = "All bootstrap estimates should be finite (not NaN/Inf)")
})
