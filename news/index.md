# Changelog

## surveybootstrap 0.2.1

- Fixed a warning (“Unknown or uninitialised column: `weight.scale`”)
  that appeared when using
  [`bootstrap.estimates()`](http://dennisfeehan.org/surveybootstrap/reference/bootstrap.estimates.md)
  with
  [`rescaled.bootstrap.sample()`](http://dennisfeehan.org/surveybootstrap/reference/rescaled.bootstrap.sample.md):
  the index column was being extracted as a one-column tibble instead of
  a plain vector, causing tibble’s `$` accessor to misreport the column
  as missing
- Fixed a deprecation warning from rlang (“Unquoting language objects
  with `!!!` is deprecated”): replaced `!!!psu.vars` with
  `across(all_of(...))` / `all_of(...)` throughout
  [`rescaled.bootstrap.sample()`](http://dennisfeehan.org/surveybootstrap/reference/rescaled.bootstrap.sample.md)
  and
  [`get.rescaled.bootstrap.weights()`](http://dennisfeehan.org/surveybootstrap/reference/get.rescaled.bootstrap.weights.md)
- Fixed
  [`rescaled.bootstrap.weights()`](http://dennisfeehan.org/surveybootstrap/reference/rescaled.bootstrap.weights.md)
  producing column names `boot_weight_rep.1`, `boot_weight_rep.2`, …
  instead of the documented `boot_weight_1`, `boot_weight_2`, …
- Added vignette “The rescaled bootstrap” illustrating the main workflow
  with the MU284 example data
- Added pkgdown site
- [`rescaled.bootstrap.sample()`](http://dennisfeehan.org/surveybootstrap/reference/rescaled.bootstrap.sample.md)
  and
  [`get.rescaled.bootstrap.weights()`](http://dennisfeehan.org/surveybootstrap/reference/get.rescaled.bootstrap.weights.md)
  now issue a warning when any stratum contains fewer than 2 PSUs
  (variance cannot be estimated in that case; bootstrap weights for the
  stratum will be 0)
- Added regression tests covering the above bug fixes and the single-PSU
  warning

## surveybootstrap 0.2.0

- Added
  [`get.rescaled.bootstrap.weights()`](http://dennisfeehan.org/surveybootstrap/reference/get.rescaled.bootstrap.weights.md):
  new function that returns rescaled bootstrap weights (and optionally
  cluster counts and scaling factors) in a wide data frame, making it
  straightforward to use the rescaled bootstrap for
  jackknife-after-bootstrap (JAB) variance estimation
- [`rescaled.bootstrap.sample()`](http://dennisfeehan.org/surveybootstrap/reference/rescaled.bootstrap.sample.md)
  gains an `include_cc` argument; when `TRUE` it returns per-rep cluster
  inclusion counts alongside the weight-scaling factors, which are
  needed for JAB calculations
- Fixed a bug where
  [`bootstrap.estimates()`](http://dennisfeehan.org/surveybootstrap/reference/bootstrap.estimates.md)
  failed silently when used with
  [`rescaled.bootstrap.sample()`](http://dennisfeehan.org/surveybootstrap/reference/rescaled.bootstrap.sample.md):
  the weight-scaling column was named `weight_scale` (underscore) in the
  sampler but read as `weight.scale` (dot) by the estimator, causing
  weights to collapse to `numeric(0)` and all estimates to be `NaN` /
  `Inf`
- Expanded test suite: added structural and semantic tests for
  [`rescaled.bootstrap.sample()`](http://dennisfeehan.org/surveybootstrap/reference/rescaled.bootstrap.sample.md),
  [`srs.bootstrap.sample()`](http://dennisfeehan.org/surveybootstrap/reference/srs.bootstrap.sample.md),
  [`get.rescaled.bootstrap.weights()`](http://dennisfeehan.org/surveybootstrap/reference/get.rescaled.bootstrap.weights.md),
  and
  [`rescaled.bootstrap.weights()`](http://dennisfeehan.org/surveybootstrap/reference/rescaled.bootstrap.weights.md);
  added a regression test that directly catches the `weight.scale`
  naming bug; added tests for the
  [`weighted.mean()`](http://dennisfeehan.org/surveybootstrap/reference/weighted.mean.md)
  helper

## surveybootstrap 0.1.0.9000

## surveybootstrap 0.0.3

CRAN release: 2023-01-09

- fixed some dependency errors that suddenly arose (explicitly use
  plyr:: in a few places)
- dropped now-deprecated group_indices\_() in favor of cur_group_index()
- changed docs to use Markdown roxygen
- added MU284 dataset for examples and tests
- added examples to exported function docs

## surveybootstrap 0.0.2

- moved unit tests for goc and helper functions in from networkreporting
- added rescaled.bootstrap.weights helper function for case when a
  tibble with all of the rescaled bootstrap weights is useful
- added dependencies to tidyverse, magrittr, and purrr

## surveybootstrap 0.0.1

CRAN release: 2016-05-02

- first CRAN version, May 2016
