<div id="main" class="col-md-9" role="main">

Get all of the values of the given variable found among members of a chain
==========================================================================

<div class="ref-description section level2">

Get all of the values of the given variable found among members of a
chain

</div>

<div class="section level2">

Usage
-----

<div class="sourceCode">

``` r
chain.vals(chain, qoi.var = "uid")
```

</div>

</div>

<div class="section level2">

Arguments
---------

-   chain:

    The chain to get values from

-   qoi.var:

    The name of the variable to get from each member of the chain

</div>

<div class="section level2">

Value
-----

A vector with all of the values of `qoi.var` found in this chain.
(Currently, the order of the values in the vector is not guaranteed.)

</div>

</div>
