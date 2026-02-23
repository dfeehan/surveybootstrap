<div id="main" class="col-md-9" role="main">

rescaled.bootstrap.weights
==========================

<div class="ref-description section level2">

This function creates a dataset with rescaled bootstrap weights; it can
be a helpful alternative to `bootstrap.estimates` in some situations

</div>

<div class="section level2">

Usage
-----

<div class="sourceCode">

``` r
rescaled.bootstrap.weights(
  survey.data,
  survey.design,
  num.reps,
  weights = NULL,
  idvar,
  verbose = TRUE,
  include_cc = FALSE,
  parallel = FALSE,
  paropts = NULL
)
```

</div>

</div>

<div class="section level2">

Arguments
---------

-   survey.data:

    The dataset to use

-   survey.design:

    A formula describing the design of the survey (see Details in
    `bootstrap.estimates()` help page)

-   num.reps:

    the number of bootstrap replication samples to draw

-   weights:

    weights to use in estimation (or NULL, if none)

-   idvar:

    the name of the column in `survey.data` that has the respondent id

-   verbose:

    if TRUE, produce lots of feedback about what is going on

-   include\_cc:

    if TRUE, include cluster counts (see `rescaled_bootstrap`)

-   parallel:

    if TRUE, use the plyr library's .parallel argument to produce
    bootstrap resamples and estimates in parallel

-   paropts:

    if not NULL, additional arguments to pass along to the
    parallelization routine

</div>

<div class="section level2">

Value
-----

if no summary.fn is specified, then return the list of estimates
produced by estimator.fn; if summary.fn is specified, then return its
output

</div>

<div class="section level2">

Details
-------

The formula describing the survey design should have the form
`~ psu_v1 + psu_v2 + ... + strata(strata_v1 + strata_v2 + ...)`, where
`psu_v1, ...` are the variables identifying primary sampling units
(PSUs) and `strata_v1, ...` identify the strata

</div>

<div class="section level2">

Examples
--------

<div class="sourceCode">

``` r
survey <- MU284.complex.surveys[[1]]
rescaled.bootstrap.weights(survey.data = survey,
                          survey.design = ~ CL,
                          weights='sample_weight',
                          idvar='LABEL',
                          num.reps = 2)
#>    LABEL boot_weight_rep.1 boot_weight_rep.2
#> 1     96                 0          33.33333
#> 2     99                 0          33.33333
#> 3    101                 0          33.33333
#> 4     16                 0          20.83333
#> 5     18                 0          20.83333
#> 6     19                 0          20.83333
#> 7    221                 0           0.00000
#> 8    223                 0           0.00000
#> 9    224                 0           0.00000
#> 10   188                50           0.00000
#> 11   190                50           0.00000
#> 12   191                50           0.00000
#> 13   199                50          50.00000
#> 14   200                50          50.00000
#> 15   203                50          50.00000


if (FALSE) { # \dontrun{
bootweights <- rescaled.bootstrap.weights(
                                         # formula describing survey design:
                                         # psu and strata
                                         survey.design = ~ psu +
                                                           stratum(stratum_analysis),
                                         num.reps=10000,
                                         # column with respondent ids
                                         idvar='caseid',
                                         # column with sampling weight
                                         weights='wwgt',
                                         # survey dataset
                                         survey.data=mw.ego)

} # }
```

</div>

</div>

</div>
