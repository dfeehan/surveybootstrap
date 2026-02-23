<div id="main" class="col-md-9" role="main">

Determine whether or not one id is a parent of another
======================================================

<div class="ref-description section level2">

This function allows us to determine which ids are directly descended
from which other ones. It is the only part of the code that relies on
the ID format used by the Curitiba study (see Details); by modifying
this function, it should be possible to adapt this code to another
study.

</div>

<div class="section level2">

Usage
-----

<div class="sourceCode">

``` r
is.child.ct(id, seed.id)
```

</div>

</div>

<div class="section level2">

Arguments
---------

-   id:

    the id of the potential child

-   seed.id:

    the id of the potential parent

</div>

<div class="section level2">

Value
-----

TRUE if `id` is the direct descendant of `seed.id` and FALSE otherwise

</div>

<div class="section level2">

Details
-------

See:

-   Salganik, M. J., Fazito, D., Bertoni, N., Abdo, A. H., Mello, M. B.,
    & Bastos, F. I. (2011). Assessing network scale-up estimates for
    groups most at risk of HIV/AIDS: evidence from a multiple-method
    study of heavy drug users in Curitiba, Brazil. *American journal of
    epidemiology*, 174(10), 1190-1196.

</div>

</div>
