<div id="main" class="col-md-9" role="main">

Grab a function based on its name
=================================

<div class="ref-description section level2">

Helper to grab a function that is passed in as an argument

</div>

<div class="section level2">

Usage
-----

<div class="sourceCode">

``` r
get.fn(fn, env = parent.frame())
```

</div>

</div>

<div class="section level2">

Arguments
---------

-   fn:

    The function to search for

-   env:

    The environment to start searching in

</div>

<div class="section level2">

Value
-----

`fn`, if `fn` is already a function; otherwise, the first function found
in env or one of its parents whose name is `fn`

</div>

<div class="section level2">

Details
-------

This is based on Hadley Wickham's response to an SO post:
<https://stackoverflow.com/questions/14183766/match-fun-provide-error-with-functions-defined-inside-functions>
with some minor modifications

</div>

</div>
