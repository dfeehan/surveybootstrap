


## tests for get.var and get.weights
context("helpers - get.var")

gv.tmp <- data.frame('A'=1:5, 'B'=5:1, 'C'=c(1,2,NA,3,NA), weights=10:14)

expect_that(get.var(gv.tmp, 'A'), equals(1:5))
expect_that(get.var(gv.tmp, 1), equals(1:5))

expect_that(get.var(gv.tmp, 'B'), equals(5:1))
expect_that(get.var(gv.tmp, 2), equals(5:1))

expect_that(get.var(gv.tmp, 'C'), equals(c(1,2,NA,3,NA)))
expect_that(get.var(gv.tmp, 3), equals(c(1,2,NA,3,NA)))

expect_that(get.var(gv.tmp, NULL), equals(rep(NA,5)))
expect_that(get.var(gv.tmp, NULL, default=-1), equals(rep(-1,5)))

expect_that(get.var(gv.tmp, 10), throws_error())
expect_that(get.var(gv.tmp, c(1,2)), throws_error())
expect_that(get.var(gv.tmp, NA), throws_error())

expect_that(get.weights(gv.tmp, "weights"), equals(10:14))
expect_that(get.weights(gv.tmp, NULL), equals(rep(1,5)))

## TODO -- write test for df.to.kpvec


## TODO -- write test for parse_design
d1 <- parse_design( ~ a + b + strata(c + d))
expect_that(d1$psu.formula, is_identical_to( ~ a + b))
expect_that(d1$strata.formula, is_identical_to( ~ c + d))

d2 <- parse_design( ~ a + b)
expect_that(d2$psu.formula, is_identical_to( ~ a + b))
expect_that(d2$strata.formula, is_identical_to(NULL))

d3 <- parse_design( ~ a + strata(b))
expect_that(d3$psu.formula, is_identical_to(~ a))
expect_that(d3$strata.formula, is_identical_to(~ b))

## TODO -- write test for topcodev.var
## TODO -- write test for topcode.data

context("helpers - weighted.mean")

test_that("equal weights give the same result as mean()", {
  x <- c(1, 2, 3, 4, 5)
  w <- rep(1, 5)
  expect_equal(surveybootstrap:::weighted.mean(x, w), mean(x))
})

test_that("a zero-weight observation does not affect the result", {
  x <- c(1, 2, 3, 999)
  w <- c(1, 1, 1, 0)
  expect_equal(surveybootstrap:::weighted.mean(x, w), mean(c(1, 2, 3)))
})

test_that("na.rm=FALSE with NA in x returns NA", {
  x <- c(1, 2, NA, 4)
  w <- c(1, 1, 1, 1)
  expect_true(is.na(surveybootstrap:::weighted.mean(x, w, na.rm = FALSE)))
})

test_that("na.rm=TRUE with NA in x ignores missing values", {
  x <- c(1, 2, NA, 4)
  w <- c(1, 1, 1, 1)
  expect_equal(surveybootstrap:::weighted.mean(x, w, na.rm = TRUE), mean(c(1, 2, 4)))
})

## TODO -- do we need weighted mean fn?

## TODO -- write test for parse.total.popn.size

## TODO -- write test for estimate.error

## TODO -- write test for attributes.to.long (or, move to dhstools?)

