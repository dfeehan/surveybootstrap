<div id="main" class="col-md-9" role="main">

Estimate degree distributions by trait
======================================

<div class="ref-description section level2">

Break down RDS degree distributions by trait, and return an object which
has the degrees for each trait as well as functions to draw degrees from
each trait.

</div>

<div class="section level2">

Usage
-----

<div class="sourceCode">

``` r
estimate.degree.distns(survey.data, d.hat.vals, traits, keep.vars = NULL)
```

</div>

</div>

<div class="section level2">

Arguments
---------

-   survey.data:

    The respondent info

-   d.hat.vals:

    The variable that contains the degrees for each respondent

-   traits:

    A vector of the names of the columns of `survey.data` which refer to
    the traits

-   keep.vars:

    Additional vars to return along with degrees

</div>

<div class="section level2">

Value
-----

An object with

-   `distns` a list with one entry per trait value; each

-   `draw.degrees.fn` a function which gets called with one

-   `keep.vars` the name of the other vars that are kept (if any)

</div>

<div class="section level2">

Details
-------

One of the items returned as a result is a function, `draw.degrees.fn`,
which takes one argument, `traits`. This is a vector of traits and, for
each entry in this vector, `draw.degress.fn` returns a draw from the
empirical distribution of degrees among respondents with that trait. So,
`draw.degrees.fn(c("0.0", "0.1", "0.1")` would return a degree drawn
uniformly at random from among the observed degrees of respondents with
trait "0.0" and then two degrees from respondents with trait "0.1"

</div>

</div>
