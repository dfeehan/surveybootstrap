<div id="main" class="col-md-9" role="main">

Draw RDS bootstrap resamples
============================

<div class="ref-description section level2">

Draw bootstrap resamples for an RDS dataset, using the algorithm
described in the supporting online material of Weir et al 2012 (see
`rds.boot.draw.chain()` ).

</div>

<div class="section level2">

Usage
-----

<div class="sourceCode">

``` r
rds.chain.boot.draws(chains, mm, dd, num.reps, keep.vars = NULL)
```

</div>

</div>

<div class="section level2">

Arguments
---------

-   chains:

    A list whose entries are the chains we want to resample

-   mm:

    The mixing model

-   dd:

    The degree distributions

-   num.reps:

    The number of bootstrap resamples we want

-   keep.vars:

    If not `NULL`, then the names of variables from the original dataset
    we want appended to each bootstrap resampled dataset (default is
    `NULL`)

</div>

<div class="section level2">

Value
-----

A list of length `num.reps`; each entry in the list has one
bootstrap-resampled dataset

</div>

</div>
