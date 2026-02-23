<div id="main" class="col-md-9" role="main">

Take a set of traits and turn into a string
===========================================

<div class="ref-description section level2">

This is a helper function that is useful when we wish to make several
traits into one variable

</div>

<div class="section level2">

Usage
-----

<div class="sourceCode">

``` r
traits.to.string(data, traits, na.action = "drop", sep = ".")
```

</div>

</div>

<div class="section level2">

Arguments
---------

-   data:

    The respondent info

-   traits:

    The names of the traits to build the model on

-   na.action:

    Defaults to 'drop' (meaning all rows of data with any missingness on
    the traits are dropped). Anything else means `NA`s are treated like
    any other value.

-   sep:

    The separator character used to combine values

</div>

<div class="section level2">

Value
-----

A list whose entries are

-   `used.idx`, which indicates which rows from the original dataset
    were used (may not be all of them if there is missingness); and

-   `traits`, which has the string version of the traits

</div>

</div>
