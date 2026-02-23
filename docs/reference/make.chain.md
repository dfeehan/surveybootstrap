<div id="main" class="col-md-9" role="main">

Build an RDS seed's chain from the dataset
==========================================

<div class="ref-description section level2">

Note that this assumes that the chain is a tree (no loops)

</div>

<div class="section level2">

Usage
-----

<div class="sourceCode">

``` r
make.chain(seed.id, survey.data, is.child.fn = is.child.ct)
```

</div>

</div>

<div class="section level2">

Arguments
---------

-   seed.id:

    The id of the seed whose chain we wish to build from the dataset

-   survey.data:

    The dataset

-   is.child.fn:

    A function which takes two ids as arguments; it is expected to
    return `TRUE` if the second argument is the parent of the first, and
    `FALSE` otherwise. it defaults to `is.child.ct()`

</div>

<div class="section level2">

Value
-----

info

</div>

</div>
