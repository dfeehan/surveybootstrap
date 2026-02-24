# bootstrap.estimates

Use a given bootstrap method to estimate sampling uncertainty from a
given estimator.

## Usage

``` r
bootstrap.estimates(
  survey.data,
  survey.design,
  bootstrap.fn,
  estimator.fn,
  num.reps,
  weights = NULL,
  ...,
  summary.fn = NULL,
  verbose = TRUE,
  parallel = FALSE,
  paropts = NULL
)
```

## Arguments

- survey.data:

  The dataset to use

- survey.design:

  A formula describing the design of the survey (see Details below)

- bootstrap.fn:

  Name of the method to be used to take bootstrap resamples

- estimator.fn:

  The name of a function which, given a dataset like `survey.data` and
  arguments in `...`, will produce an estimate of interest

- num.reps:

  The number of bootstrap replication samples to draw

- weights:

  Weights to use in estimation (or `NULL`, if none)

- ...:

  additional arguments which will be passed on to `estimator.fn`

- summary.fn:

  (Optional) Name of a function which, given the set of estimates
  produced by `estimator.fn`, summarizes them. If not specified, all of
  the estimates are returned in a list.

- verbose:

  If `TRUE`, produce lots of feedback about what is going on

- parallel:

  If `TRUE`, use the `plyr` library's `.parallel` argument to produce
  bootstrap resamples and estimates in parallel

- paropts:

  If not `NULL`, additional arguments to pass along to the
  parallelization routine

## Value

If `summary.fn` is not specified, then return the list of estimates
produced by `estimator.fn`; if `summary.fn` is specified, then return
its output

## Details

The formula describing the survey design should have the form
`~ psu_v1 + psu_v2 + ... + strata(strata_v1 + strata_v2 + ...)`, where
`psu_v1, ...` are the variables identifying primary sampling units
(PSUs) and `strata_v1, ...` identifies the strata

## Examples

``` r
# example using a simple random sample

survey <- MU284.surveys[[1]]

estimator <- function(survey.data, weights) {
  plyr::summarise(survey.data,
                  T82.hat = sum(S82 * weights))
}

ex.mu284 <- bootstrap.estimates(
   survey.design = ~1,
   num.reps = 10,
   estimator.fn = estimator,
   weights='sample_weight',
   bootstrap.fn = 'srs.bootstrap.sample',
   survey.data=survey)

if (FALSE) { # \dontrun{
idu.est <- bootstrap.estimates(
  ## this describes the sampling design of the
  ## survey; here, the PSUs are given by the
  ## variable cluster, and the strata are given
  ## by the variable region
  survey.design = ~ cluster + strata(region),
  ## the number of bootstrap resamples to obtain
  num.reps=1000,
  ## this is the name of the function
  ## we want to use to produce an estimate
  ## from each bootstrapped dataset
  estimator.fn="our.estimator",
  ## these are the sampling weights
  weights="indweight",
  ## this is the name of the type of bootstrap
  ## we wish to use
  bootstrap.fn="rescaled.bootstrap.sample",
  ## our dataset
  survey.data=example.survey,
  ## other parameters we need to pass
  ## to the estimator function
  d.hat.vals=d.hat,
  total.popn.size=tot.pop.size,
  y.vals="clients",
  missing="complete.obs")

} # }
```
