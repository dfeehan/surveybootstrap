<div id="main" class="col-md-9" role="main">

Run a markov model
==================

<div class="ref-description section level2">

Run a given markov model for n time steps, starting at a specified
state.

</div>

<div class="section level2">

Usage
-----

<div class="sourceCode">

``` r
mc.sim(mm, start, n)
```

</div>

</div>

<div class="section level2">

Arguments
---------

-   mm:

    The markov model object returned by `estimate.mixing()`

-   start:

    The name of the state to start in

-   n:

    The number of time-steps to run through

</div>

<div class="section level2">

Value
-----

A vector with the state visited at each time step. the first entry has
the starting state

</div>

<div class="section level2">

Details
-------

This uses the markov model produced by `estimate.mixing()`

</div>

</div>
