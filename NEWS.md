# surveybootstrap 0.2.1

* Fixed a warning ("Unknown or uninitialised column: `weight.scale`") that
  appeared when using `bootstrap.estimates()` with `rescaled.bootstrap.sample()`:
  the index column was being extracted as a one-column tibble instead of a plain
  vector, causing tibble's `$` accessor to misreport the column as missing
* Fixed a deprecation warning from rlang ("Unquoting language objects with `!!!`
  is deprecated"): replaced `!!!psu.vars` with `across(all_of(...))` / `all_of(...)`
  throughout `rescaled.bootstrap.sample()` and `get.rescaled.bootstrap.weights()`
* Fixed `rescaled.bootstrap.weights()` producing column names `boot_weight_rep.1`,
  `boot_weight_rep.2`, … instead of the documented `boot_weight_1`, `boot_weight_2`, …
* Added vignette "The rescaled bootstrap" illustrating the main workflow with the
  MU284 example data
* Added pkgdown site
* Added regression tests covering the above bug fixes

# surveybootstrap 0.2.0

* Added `get.rescaled.bootstrap.weights()`: new function that returns rescaled
  bootstrap weights (and optionally cluster counts and scaling factors) in a
  wide data frame, making it straightforward to use the rescaled bootstrap for
  jackknife-after-bootstrap (JAB) variance estimation
* `rescaled.bootstrap.sample()` gains an `include_cc` argument; when `TRUE` it
  returns per-rep cluster inclusion counts alongside the weight-scaling factors,
  which are needed for JAB calculations
* Fixed a bug where `bootstrap.estimates()` failed silently when used with
  `rescaled.bootstrap.sample()`: the weight-scaling column was named
  `weight_scale` (underscore) in the sampler but read as `weight.scale` (dot)
  by the estimator, causing weights to collapse to `numeric(0)` and all
  estimates to be `NaN` / `Inf`
* Expanded test suite: added structural and semantic tests for
  `rescaled.bootstrap.sample()`, `srs.bootstrap.sample()`,
  `get.rescaled.bootstrap.weights()`, and `rescaled.bootstrap.weights()`;
  added a regression test that directly catches the `weight.scale` naming bug;
  added tests for the `weighted.mean()` helper

# surveybootstrap 0.1.0.9000


# surveybootstrap 0.0.3

* fixed some dependency errors that suddenly arose (explicitly use plyr:: in a few places)
* dropped now-deprecated group_indices_() in favor of cur_group_index()
* changed docs to use Markdown roxygen
* added MU284 dataset for examples and tests
* added examples to exported function docs

# surveybootstrap 0.0.2

* moved unit tests for goc and helper functions in from networkreporting
* added rescaled.bootstrap.weights helper function for case when a tibble with all of the rescaled bootstrap weights is useful
* added dependencies to tidyverse, magrittr, and purrr

# surveybootstrap 0.0.1

* first CRAN version, May 2016
