<div id="main" class="col-md-9" role="main">

Changelog
=========

<div class="section level2">

surveybootstrap 0.2.0
---------------------

-   Added `get.rescaled.bootstrap.weights()`: new function that returns
    rescaled bootstrap weights (and optionally cluster counts and
    scaling factors) in a wide data frame, making it straightforward to
    use the rescaled bootstrap for jackknife-after-bootstrap (JAB)
    variance estimation
-   `rescaled.bootstrap.sample()` gains an `include_cc` argument; when
    `TRUE` it returns per-rep cluster inclusion counts alongside the
    weight-scaling factors, which are needed for JAB calculations
-   Fixed a bug where `bootstrap.estimates()` failed silently when used
    with `rescaled.bootstrap.sample()`: the weight-scaling column was
    named `weight_scale` (underscore) in the sampler but read as
    `weight.scale` (dot) by the estimator, causing weights to collapse
    to `numeric(0)` and all estimates to be `NaN` / `Inf`
-   Expanded test suite: added structural and semantic tests for
    `rescaled.bootstrap.sample()`, `srs.bootstrap.sample()`,
    `get.rescaled.bootstrap.weights()`, and
    `rescaled.bootstrap.weights()`; added a regression test that
    directly catches the `weight.scale` naming bug; added tests for the
    `weighted.mean()` helper

</div>

<div class="section level2">

surveybootstrap 0.1.0.9000
--------------------------

</div>

<div class="section level2">

surveybootstrap 0.0.3
---------------------

CRAN release: 2023-01-09

-   fixed some dependency errors that suddenly arose (explicitly use
    plyr:: in a few places)
-   dropped now-deprecated group\_indices\_() in favor of
    cur\_group\_index()
-   changed docs to use Markdown roxygen
-   added MU284 dataset for examples and tests
-   added examples to exported function docs

</div>

<div class="section level2">

surveybootstrap 0.0.2
---------------------

-   moved unit tests for goc and helper functions in from
    networkreporting
-   added rescaled.bootstrap.weights helper function for case when a
    tibble with all of the rescaled bootstrap weights is useful
-   added dependencies to tidyverse, magrittr, and purrr

</div>

<div class="section level2">

surveybootstrap 0.0.1
---------------------

CRAN release: 2016-05-02

-   first CRAN version, May 2016

</div>

</div>
