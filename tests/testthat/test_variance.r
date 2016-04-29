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
load("mu284.RData")


## size of the entire population
tot.pop.size <- 10718378

## column names for connections to groups of known size
kp.q <- paste(example.knownpop.dat$known.popn)


boot.example <- example.survey
attr(boot.example, 'total.popn.size') <- tot.pop.size

boot.example <- add.kp(boot.example,
                       df.to.kpvec(example.knownpop.dat,
                                   kp.var='known.popn',
                                   kp.value='size'))

knownpop.tot <- sum(example.knownpop.dat$size)

boot.example$d.hat <- kp.individual.estimator_(boot.example, 
                               known.populations=kp.q,
                               total.kp.size=knownpop.tot,
                               alter.popn.size=tot.pop.size)$dbar.Fcell.F

M1 <- 10
M2 <- 50

test.bootfn <- function(bfn) {

  boot1 <- bootstrap.estimates(survey.data=boot.example,
                               survey.design= ~ cluster + strata(region),
                               num.reps=M1,
                               estimator.fn="nsum.estimator",
                               ## these args below all go to
                               ## the estimator
                               kp.method=TRUE,
                               return.plot=FALSE,
                               weights="indweight",
                               missing="complete.obs",
                               y.vals="clients",
                               verbose=FALSE,
                               bootstrap.fn=bfn,
                               d.hat.vals="d.hat")

  boot2 <- bootstrap.estimates(survey.data=boot.example,
                               survey.design= ~ cluster + strata(region),
                               num.reps=M2,
                               estimator.fn="nsum.estimator",
                               ## these args below all go to
                               ## the estimator
                               kp.method=TRUE,
                               return.plot=FALSE,
                               weights="indweight",
                               missing="complete.obs",
                               y.vals="clients",
                               verbose=FALSE,
                               bootstrap.fn=bfn,
                               d.hat.vals="d.hat")

  ests1 <- laply(boot1, function(x) { x$estimate })
  nums1 <- laply(boot1, function(x) { x$tot.connections })
  denoms1 <- laply(boot1, function(x) { x$sum.d.hat })

  ests2 <- laply(boot2, function(x) { x$estimate })
  nums2 <- laply(boot2, function(x) { x$tot.connections })
  denoms2 <- laply(boot2, function(x) { x$sum.d.hat })

  ## be sure that there is variation in the numerator
  ## the denominator and the estimates
  expect_that(var(ests1) == 0, is_false())
  expect_that(var(nums1) == 0, is_false())
  expect_that(var(denoms1) == 0, is_false())

  expect_that(var(ests2) == 0, is_false())
  expect_that(var(nums2) == 0, is_false())
  expect_that(var(denoms2) == 0, is_false())

  ## be sure that the estimates, the numerator,
  ## and the denominator are always nonnegative
  expect_that(all(ests1 >= 0), is_true())
  expect_that(all(nums1 >= 0), is_true())
  expect_that(all(denoms1 >= 0), is_true())

  expect_that(all(ests2 >= 0), is_true())
  expect_that(all(nums2 >= 0), is_true())
  expect_that(all(denoms2 >= 0), is_true())

  ## eventually, we might want to do more with the resamples,
  ## so return them
  return(list(boot1, boot2))
}

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

test.boot <- llply(MU284.surveys,
                   function(svy) { do.call("rbind", rbsfn(survey.data=svy)) })

test.boot.summ <- ldply(test.boot,
                        summarize,
                        mean.TS82.hat=mean(TS82.hat),
                        mean.R.RMT85.P85.hat=mean(R.RMT85.P85.hat),
                        sd.TS82.hat=sd(TS82.hat),
                        sd.R.RMT85.P85.hat=sd(R.RMT85.P85.hat))

## TODO -- how to figure out what tolerance to use?
## for now, using .1 for everything, but this is pretty big; also,
## we'd expect different tolerances for different qois, i think.
qoi <- colnames(test.boot.summ)

l_ply(qoi,
      function(this.qoi) {
          l_ply(1:nrow(test.boot.summ),
                function(idx) {
                    expect_that(test.boot.summ[idx,this.qoi],
                                equals(MU284.boot.res.summ[idx,this.qoi],
                                       tolerance=.05),
                                label=paste0("qty: ", this.qoi,
                                             "; MU284 survey #", idx))
                })
      })

#########################################
## simple random sample (SRS) bootstrap
context("variance estimators - srs bootstrap - sanity checks")

tmp <- test.bootfn("srs.bootstrap.sample")


#########################################
## rescaled (Rao / Wu) bootstrap
context("variance estimators - rescaled bootstrap - sanity checks")

tmp <- test.bootfn("rescaled.bootstrap.sample")


## TODO -- LEFT OFF HERE...
##  * consider increasing M
##  * test that mean is close to raw value
##  * draft the rest of the tests...


## raw estimates
## rw.raw <- nsum.internal.validation(rw.data,
##                                    kp.method=TRUE,
##                                    na.rm=TRUE,
##                                    weights="indweight")

## conduct the bootstrap resamples
## rw.rbs <- rw.befn(bootstrap.fn="rescaled.bootstrap.sample")
## rw.srs <- rw.befn(bootstrap.fn="srs.bootstrap.sample")
