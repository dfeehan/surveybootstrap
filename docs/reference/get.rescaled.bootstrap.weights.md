<div id="main" class="col-md-9" role="main">

rescaled.bootstrap.weights
==========================

<div class="ref-description section level2">

Get a dataset of rescaled bootstrap weights.

</div>

<div class="section level2">

Usage
-----

<div class="sourceCode">

``` r
get.rescaled.bootstrap.weights(
  survey.data,
  survey.design,
  idvar,
  weights = NULL,
  parallel = FALSE,
  paropts = NULL,
  num.reps = 1,
  include_scaling_factors = FALSE,
  include_cc = FALSE
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

    A formula describing the design of the survey (see Details)

-   idvar:

    the name of the column in `survey.data` that has the respondent id

-   weights:

    Survey weights to be rescaled (or `NULL`, if none)

-   parallel:

    If `TRUE`, use parallelization (via `plyr`)

-   paropts:

    An optional list of arguments passed on to `plyr` to control details
    of parallelization

-   num.reps:

    The number of bootstrap replication samples to draw

-   include\_scaling\_factors:

    If `TRUE`, include `weight_scaling_factor` in the output

-   include\_cc:

    If `TRUE`, include cluster counts in the output. These cluster
    counts are needed for jackknife after bootstrap calculations

</div>

<div class="section level2">

Value
-----

<div class="sourceCode">

    A list with entries:

    * `orig_weights` - A dataframe with rows corresponding to the original
       `survey.data` and a column with the original (un rescaled) weight for the obs.
       Can be useful for debugging.

    * `boot_weights` - A dataframe with rows corresponding to the original
       `survey.data` and one column for each bootstrap rep. The entries have
       the rescaled bootstrap weights for each row and bootstrap rep.

    * `weight_scaling_factor` - (if `include_scaling_factors` is `TRUE`). A dataframe with rows corresponding to the original
       `survey.data` and one column for each bootstrap rep. The entries have
       the values by which the original weights were scaled to produce `boot_weights`.
       You typically won't need to use this dataframe, but it can be helpful for debugging.

    * `cluster_counts` - (if `include_cc` is `TRUE`). A list with `num.reps` entries.
      Each entry is a dataset which has the variable '.cluster_id',
      the values of columns that specify the PSU (from the `survey.design` formula),
      and `cluster_count` (the number of times the given PSU was resampled))

</div>

</div>

<div class="section level2">

Details
-------

This is a new function being added as part of a refactor in v.0.2. This
function focuses on being able to return a dataframe with bootstrap
weights, and optionally also a dataframe with cluster/PSU inclusion
counts.

`survey.design` is a formula of the form

`weight ~ psu_vars + strata(strata_vars)`

where:

-   `weight` is the variable with the survey weights

-   `psu_vars` has the form `psu_v1 + psu_v2 + ...`, where primary
    sampling units (PSUs) are determined by `psu_v1`, etc

-   `strata_vars` has the form `strata_v1 + strata_v2 + ...`, which
    determine strata

Note that we assume that the formula uniquely specifies PSUs. This will
always be true if the PSUs were selected without replacement. If they
were selected with replacement, then it will be necessary to make each
realization of a given PSU in the sample a unique id. The code below
assumes that all observations within each PSU (as identified by the
design formula) are from the same draw of the PSU.

The rescaled bootstrap technique works by adjusting the estimation
weights based on the number of times each row is included in the
resamples. If a row is never selected, it is still included in the
returned results, but its weight will be set to 0. It is therefore
important to use estimators that make use of the estimation weights on
the resampled datasets.

We always take m\_i = n\_i - 1, according to the advice presented in Rao
and Wu (1988) and Rust and Rao (1996).

(This is a C++ version; a previous version, written in pure R, is called
`rescaled.bootstrap.sample.pureR()` )

References:

-   Rust, Keith F., and J. N. K. Rao. "Variance estimation for complex
    surveys using replication techniques." *Statistical methods in
    medical research* 5.3 (1996): 283-310.

-   Rao, Jon NK, and C. F. J. Wu. "Resampling inference with complex
    survey data." *Journal of the American Statistical Association*
    83.401 (1988): 231-241.

</div>

<div class="section level2">

Examples
--------

<div class="sourceCode">

``` r
survey <- MU284.complex.surveys[[1]]
boot_surveys <- get.rescaled.bootstrap.weights(survey.data = survey,
                                               survey.design = ~ CL,
                                               num.reps = 2)
#> Error in select(., .internal_id, .cluster_id, any_of(psu.var.names), one_of(idvar)): ℹ In argument: `one_of(idvar)`.
#> Caused by error:
#> ! argument "idvar" is missing, with no default
```

</div>

</div>

</div>
