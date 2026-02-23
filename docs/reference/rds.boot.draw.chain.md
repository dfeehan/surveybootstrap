<div id="main" class="col-md-9" role="main">

Draw RDS bootstrap resamples for one chain
==========================================

<div class="ref-description section level2">

This function uses the algorithm described in the supporting online
material for Weir et al 2012 (see Details) to take bootstrap resamples
of one chain from an RDS dataset.

</div>

<div class="section level2">

Usage
-----

<div class="sourceCode">

``` r
rds.boot.draw.chain(chain, mm, dd, parent.trait, idvar = "uid")
```

</div>

</div>

<div class="section level2">

Arguments
---------

-   chain:

    The chain to draw resamples for

-   mm:

    The mixing model to use

-   dd:

    The degree distns to use

-   parent.trait:

    A vector whose length is the number of bootstrap reps we want

-   idvar:

    The name of the variable used to label the columns of the output
    (presumably some id identifying the row in the original dataset they
    come from)

</div>

<div class="section level2">

Value
-----

A list of dataframes with one entry for each respondent in the chain.
each dataframe has one row for each bootstrap replicate. so if we take
10 bootstrap resamples of a chain of length 50, there will be 50 entries
in the list that is returned. each entry will be a dataframe with 10
rows.

</div>

<div class="section level2">

Details
-------

See

-   Weir, Sharon S., et al. "A comparison of respondent-driven and
    venue-based sampling of female sex workers in Liuzhou, China."
    *Sexually transmitted infections* 88.Suppl 2 (2012): i95-i101.

</div>

</div>
