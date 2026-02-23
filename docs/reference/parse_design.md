<div id="main" class="col-md-9" role="main">

Parse a formula that describes the design of a survey
=====================================================

<div class="ref-description section level2">

Parse a formula of the form
`~ psu_v1 + psu_v2 + ... + strata(strata_v1 + strata_v2 + ...)` into a
PSU formula and a strata formula. Note that we only document
`strata(...)` but this will also work with `stratum(...)` because it is
annoying to have to remember which one to use

</div>

<div class="section level2">

Usage
-----

<div class="sourceCode">

``` r
parse_design(formula)
```

</div>

</div>

<div class="section level2">

Arguments
---------

-   formula:

    a formula describing the sample design (see Description of
    `bootstrap.estimates()`)

</div>

<div class="section level2">

Value
-----

a list with entries `psu.formula` and `strata.formula`

</div>

</div>
