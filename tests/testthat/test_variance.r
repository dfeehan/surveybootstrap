## LEFT OFF HERE:
## - I think mu284.RData (which we're not loading anymore) has the reference
##   values for the unit test. figure out how to fix this!
##   -> these are MU284.boot.res.summ


## TODO -- test that there is some variance across a set of bootstrap
##         samples

## TODO -- for ratios, test also that there is some variance in
##         numerator and denominator

## TODO -- test paired ego / alter datasets

## TODO -- test that calling bootstrap.estimates works when
##         total.popn.size is an argument and not an attribute of
##         the data frame (had to use parent.frame(2)) to fix
##         a bug about this

## TODO -- test cases where estimates should never be negative

## TODO -- look at
## http://stackoverflow.com/questions/8898469/is-it-possible-to-use-r-package-data-in-testthat-tests-or-run-examples
## to try and figure out the real way to include package data in
## unit tests...

set.seed(12345)

#########################################
## setup
## NB: see this for possible alternatives
## http://stackoverflow.com/questions/8898469/is-it-possible-to-use-r-package-data-in-testthat-tests-or-run-examples
#load("mu284.RData")

#########################################
## rescaled (Rao / Wu) bootstrap

## TODO -- think about how many bootstrap reps
## we should use in the unit tests...
##M <- 1000
M <- 2000



context(paste0("variance estimators - rescaled bootstrap - correctness (M=", M, ")"))

rbsfn <- functional::Curry(bootstrap.estimates,
                           survey.design= ~ CL,
                           num.reps=M,
                           estimator.fn=MU284.estimator.fn,
                           weights="sample_weight",
                           bootstrap.fn="rescaled.bootstrap.sample")

test.boot <- plyr::llply(MU284.complex.surveys,
                   function(svy) { do.call("rbind", rbsfn(survey.data=svy)) })

test.boot.summ <- plyr::ldply(test.boot,
                        plyr::summarize,
                        mean.TS82.hat=mean(TS82.hat),
                        mean.R.RMT85.P85.hat=mean(R.RMT85.P85.hat),
                        sd.TS82.hat=sd(TS82.hat),
                        sd.R.RMT85.P85.hat=sd(R.RMT85.P85.hat))

## TODO -- how to figure out what tolerance to use?
## for now, using .1 for everything, but this is pretty big; also,
## we'd expect different tolerances for different qois, i think.
qoi <- colnames(test.boot.summ)

plyr::l_ply(qoi,
      function(this.qoi) {
          plyr::l_ply(1:nrow(test.boot.summ),
                function(idx) {
                    expect_that(test.boot.summ[idx,this.qoi],
                                equals(MU284.boot.res.summ[idx,this.qoi],
                                       tolerance=.05),
                                label=paste0("qty: ", this.qoi,
                                             "; MU284 survey #", idx))
                })
      })

#########################################
## rescaled.bootstrap.sample output structure

context("rescaled.bootstrap.sample - output structure")

rbs_svy <- MU284.complex.surveys[[1]]
rbs_result <- rescaled.bootstrap.sample(survey.data = rbs_svy,
                                        survey.design = ~ CL,
                                        num.reps = 5)

test_that("returns a list of length num.reps", {
  expect_equal(length(rbs_result), 5)
})

test_that("each rep has columns 'index' and 'weight.scale'", {
  for (rep in rbs_result) {
    expect_true("index" %in% names(rep))
    expect_true("weight.scale" %in% names(rep))
  }
})

test_that("index values are valid row indices", {
  for (rep in rbs_result) {
    expect_true(all(rep$index >= 1 & rep$index <= nrow(rbs_svy)))
  }
})

test_that("weight.scale values are non-negative", {
  for (rep in rbs_result) {
    expect_true(all(rep$weight.scale >= 0))
  }
})

rbs_cc_result <- rescaled.bootstrap.sample(survey.data = rbs_svy,
                                           survey.design = ~ CL,
                                           num.reps = 5,
                                           include_cc = TRUE)

test_that("include_cc=TRUE returns list with weight_factors and cluster_counts", {
  expect_true("weight_factors" %in% names(rbs_cc_result))
  expect_true("cluster_counts" %in% names(rbs_cc_result))
})

test_that("weight_factors has same structure as default output", {
  expect_equal(length(rbs_cc_result$weight_factors), 5)
  for (rep in rbs_cc_result$weight_factors) {
    expect_true("index" %in% names(rep))
    expect_true("weight.scale" %in% names(rep))
  }
})

test_that("cluster_counts entries have expected columns", {
  for (cc in rbs_cc_result$cluster_counts) {
    expect_true(".cluster_id" %in% names(cc))
    expect_true("CL" %in% names(cc))
    expect_true("cluster_count" %in% names(cc))
  }
})

#########################################
## srs.bootstrap.sample output structure

context("srs.bootstrap.sample - output structure")

srs_svy <- MU284.surveys[[1]]
srs_result <- srs.bootstrap.sample(survey.data = srs_svy,
                                   num.reps = 5)

test_that("returns a list of length num.reps", {
  expect_equal(length(srs_result), 5)
})

test_that("each rep has columns 'index' and 'weight.scale'", {
  for (rep in srs_result) {
    expect_true("index" %in% names(rep))
    expect_true("weight.scale" %in% names(rep))
  }
})

test_that("each rep has nrow(survey.data) rows", {
  for (rep in srs_result) {
    expect_equal(nrow(rep), nrow(srs_svy))
  }
})

test_that("weight.scale is always 1 (SRS does not rescale)", {
  for (rep in srs_result) {
    expect_true(all(rep$weight.scale == 1))
  }
})

test_that("index values are valid row indices", {
  for (rep in srs_result) {
    expect_true(all(rep$index >= 1 & rep$index <= nrow(srs_svy)))
  }
})

#########################################
## get.rescaled.bootstrap.weights output structure

context("get.rescaled.bootstrap.weights - output structure")

grbs_svy <- MU284.complex.surveys[[1]]
grbs_result <- get.rescaled.bootstrap.weights(survey.data = grbs_svy,
                                              survey.design = ~ CL,
                                              idvar = 'LABEL',
                                              weights = 'sample_weight',
                                              num.reps = 5)

test_that("result has orig_weights and boot_weights", {
  expect_true("orig_weights" %in% names(grbs_result))
  expect_true("boot_weights" %in% names(grbs_result))
})

test_that("boot_weights has idvar column and one column per rep", {
  expected_cols <- c("LABEL", paste0("boot_weight_", 1:5))
  expect_true(all(expected_cols %in% names(grbs_result$boot_weights)))
})

test_that("boot_weights has same number of rows as survey", {
  expect_equal(nrow(grbs_result$boot_weights), nrow(grbs_svy))
})

test_that("all boot weights are non-negative", {
  for (col in paste0("boot_weight_", 1:5)) {
    expect_true(all(grbs_result$boot_weights[[col]] >= 0))
  }
})

grbs_sf_result <- get.rescaled.bootstrap.weights(survey.data = grbs_svy,
                                                 survey.design = ~ CL,
                                                 idvar = 'LABEL',
                                                 weights = 'sample_weight',
                                                 num.reps = 5,
                                                 include_scaling_factors = TRUE)

test_that("include_scaling_factors=TRUE adds weight_scaling_factor with boot_rep_* cols", {
  expect_true("weight_scaling_factor" %in% names(grbs_sf_result))
  expected_cols <- c("LABEL", paste0("boot_rep_", 1:5))
  expect_true(all(expected_cols %in% names(grbs_sf_result$weight_scaling_factor)))
})

grbs_cc_result <- get.rescaled.bootstrap.weights(survey.data = grbs_svy,
                                                 survey.design = ~ CL,
                                                 idvar = 'LABEL',
                                                 weights = 'sample_weight',
                                                 num.reps = 5,
                                                 include_cc = TRUE)

test_that("include_cc=TRUE adds cluster_counts with PSU col and boot_rep_* cols", {
  expect_true("cluster_counts" %in% names(grbs_cc_result))
  expected_cols <- c("CL", paste0("boot_rep_", 1:5))
  expect_true(all(expected_cols %in% names(grbs_cc_result$cluster_counts)))
})

#########################################
## bootstrap.estimates + rescaled.bootstrap.sample integration (regression)

context("bootstrap.estimates - rescaled bootstrap integration (regression)")

## Regression test: verifies that bootstrap.estimates correctly reads the
## weight.scale column from rescaled.bootstrap.sample output.
## The bug (weight_scale vs weight.scale naming mismatch) would cause
## all estimates to be NA/NaN/Inf due to numeric(0) weights.

reg_svy <- MU284.complex.surveys[[1]]
reg_result <- bootstrap.estimates(survey.data = reg_svy,
                                  survey.design = ~ CL,
                                  num.reps = 10,
                                  estimator.fn = MU284.estimator.fn,
                                  weights = 'sample_weight',
                                  bootstrap.fn = 'rescaled.bootstrap.sample')

test_that("returns a list of length num.reps", {
  expect_equal(length(reg_result), 10)
})

test_that("each result has expected estimator columns", {
  for (est in reg_result) {
    expect_true("TS82.hat" %in% names(est))
    expect_true("R.RMT85.P85.hat" %in% names(est))
  }
})

test_that("all estimates are finite (regression: weight.scale column name)", {
  all_ts82  <- sapply(reg_result, function(x) x$TS82.hat)
  all_ratio <- sapply(reg_result, function(x) x$R.RMT85.P85.hat)
  expect_true(all(is.finite(all_ts82)))
  expect_true(all(is.finite(all_ratio)))
})

#########################################
## rescaled.bootstrap.weights output structure

context("rescaled.bootstrap.weights - output structure")

rbw_svy <- MU284.complex.surveys[[1]]
rbw_result <- rescaled.bootstrap.weights(survey.data = rbw_svy,
                                         survey.design = ~ CL,
                                         idvar = 'LABEL',
                                         weights = 'sample_weight',
                                         num.reps = 5)

test_that("returns a data frame", {
  expect_true(is.data.frame(rbw_result))
})

test_that("has idvar column and one boot weight column per rep", {
  expected_cols <- c("LABEL", paste0("boot_weight_", 1:5))
  expect_true(all(expected_cols %in% names(rbw_result)))
})

test_that("has same number of rows as survey", {
  expect_equal(nrow(rbw_result), nrow(rbw_svy))
})

test_that("all boot weights are non-negative", {
  for (col in paste0("boot_weight_", 1:5)) {
    expect_true(all(rbw_result[[col]] >= 0))
  }
})

#########################################
## get.rescaled.bootstrap.weights semantic content

context("get.rescaled.bootstrap.weights - semantic content")

sem_svy <- MU284.complex.surveys[[1]]

## cluster_counts semantics
sem_cc <- get.rescaled.bootstrap.weights(survey.data = sem_svy,
                                         survey.design = ~ CL,
                                         idvar = 'LABEL',
                                         weights = 'sample_weight',
                                         num.reps = 5,
                                         include_cc = TRUE)

test_that("cluster_counts has one row per unique PSU", {
  n_psus <- dplyr::n_distinct(sem_svy$CL)
  expect_equal(nrow(sem_cc$cluster_counts), n_psus)
})

test_that("cluster_counts boot_rep_* values are non-negative", {
  for (col in paste0("boot_rep_", 1:5)) {
    expect_true(all(sem_cc$cluster_counts[[col]] >= 0))
  }
})

test_that("orig_weights has a weight column with positive values", {
  expect_true("weight" %in% names(sem_cc$orig_weights))
  expect_true(all(sem_cc$orig_weights$weight > 0))
})

test_that("orig_weights has same number of rows as survey", {
  expect_equal(nrow(sem_cc$orig_weights), nrow(sem_svy))
})

## scaling factor relationship: boot_weight == orig_weight * scaling_factor
sem_sf <- get.rescaled.bootstrap.weights(survey.data = sem_svy,
                                         survey.design = ~ CL,
                                         idvar = 'LABEL',
                                         weights = 'sample_weight',
                                         num.reps = 5,
                                         include_scaling_factors = TRUE)

test_that("boot_weight equals orig_weight times scaling_factor for each rep", {
  orig_w <- sem_sf$orig_weights[, c("LABEL", "weight")]
  for (i in 1:5) {
    bw_col <- paste0("boot_weight_", i)
    sf_col <- paste0("boot_rep_", i)
    joined <- merge(orig_w,
                    sem_sf$boot_weights[, c("LABEL", bw_col)],
                    by = "LABEL")
    joined <- merge(joined,
                    sem_sf$weight_scaling_factor[, c("LABEL", sf_col)],
                    by = "LABEL")
    expect_equal(joined[[bw_col]], joined[["weight"]] * joined[[sf_col]])
  }
})

