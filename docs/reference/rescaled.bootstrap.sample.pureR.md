<div id="main" class="col-md-9" role="main">

rescaled.bootstrap.sample.pureR
===============================

<div class="ref-description section level2">

(this is the pure R version; it has been supplanted by
`rescaled.bootstrap.sample`, which is partially written in C++)

</div>

<div class="section level2">

Usage
-----

<div class="sourceCode">

``` r
rescaled.bootstrap.sample.pureR(
  survey.data,
  survey.design,
  parallel = FALSE,
  paropts = NULL,
  num.reps = 1
)
```

</div>

</div>

<div class="section level2">

Arguments
---------

-   survey.data:

    the dataset to use

-   survey.design:

    a formula describing the design of the survey (see below - TODO)

-   parallel:

    if TRUE, use parallelization (via `plyr`)

-   paropts:

    an optional list of arguments passed on to `plyr` to control details
    of parallelization

-   num.reps:

    the number of bootstrap replication samples to draw

</div>

<div class="section level2">

Value
-----

a list with `num.reps` entries. each entry is a dataset which has at
least the variables `index` (the row index of the original dataset that
was resampled) and `weight.scale` (the factor by which to multiply the
sampling weights in the original dataset).

</div>

<div class="section level2">

Details
-------

given a survey dataset and a description of the survey design (ie, which
combination of vars determines primary sampling units, and which
combination of vars determines strata), take a bunch of bootstrap
samples for the rescaled bootstrap estimator (see, eg, Rust and Rao
1996).

Note that we assume that the formula uniquely specifies PSUs. This will
always be true if the PSUs were selected without replacement. If they
were selected with replacement, then it will be necessary to make each
realization of a given PSU in the sample a unique id. Bottom line: the
code below assumes that all observations within each PSU (as identified
by the design formula) are from the same draw of the PSU.

The rescaled bootstrap technique works by adjusting the estimation
weights based on the number of times each row is included in the
resamples. If a row is never selected, it is still included in the
returned results, but its weight will be set to 0. It is therefore
important to use estimators that make use of the estimation weights on
the resampled datasets.

We always take m\_i = n\_i - 1, according to the advice presented in Rao
and Wu (1988) and Rust and Rao (1996).

`survey.design` is a formula of the form  
weight \~ psu\_vars + strata(strata\_vars), where weight is the variable
with the survey weights and psu is the variable denoting the primary
sampling unit

</div>

</div>
