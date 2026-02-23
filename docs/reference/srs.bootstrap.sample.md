<div id="main" class="col-md-9" role="main">

srs.bootstrap.sample
====================

<div class="ref-description section level2">

Given a survey dataset and a description of the survey design (ie, which
combination of vars determines primary sampling units, and which
combination of vars determines strata), take a bunch of bootstrap
samples under a simple random sampling (with repetition) scheme

</div>

<div class="section level2">

Usage
-----

<div class="sourceCode">

``` r
srs.bootstrap.sample(
  survey.data,
  num.reps = 1,
  parallel = FALSE,
  paropts = NULL,
  ...
)
```

</div>

</div>

<div class="section level2">

Arguments
---------

-   survey.data:

    The dataset to use

-   num.reps:

    The number of bootstrap replication samples to draw

-   parallel:

    If `TRUE`, use parallelization (via `plyr`)

-   paropts:

    An optional list of arguments passed on to `plyr` to control details
    of parallelization

-   ...:

    Ignored, but useful because it allows params like `survey.design`
    which are used in other bootstrap designs, to be passed in without
    error

</div>

<div class="section level2">

Value
-----

A list with `num.reps` entries. Each entry is a dataset which has at
least the variables `index` (the row index of the original dataset that
was resampled) and `weight.scale` (the factor by which to multiply the
sampling weights in the original dataset).

</div>

<div class="section level2">

Examples
--------

<div class="sourceCode">

``` r
survey <- MU284.surveys[[1]]
boot_surveys <- srs.bootstrap.sample(survey, num.reps = 2)
```

</div>

</div>

</div>
