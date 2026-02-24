# Package index

## Main bootstrap interface

High-level function for obtaining bootstrap estimates of sampling
uncertainty from any estimator.

- [`bootstrap.estimates()`](https://dfeehan.github.io/surveybootstrap/reference/bootstrap.estimates.md)
  : bootstrap.estimates

## Bootstrap sampling methods

Functions that draw bootstrap resamples from a survey dataset. Pass the
name of one of these to the `bootstrap.fn` argument of
[`bootstrap.estimates()`](https://dfeehan.github.io/surveybootstrap/reference/bootstrap.estimates.md).

- [`rescaled.bootstrap.sample()`](https://dfeehan.github.io/surveybootstrap/reference/rescaled.bootstrap.sample.md)
  : rescaled.bootstrap.sample
- [`srs.bootstrap.sample()`](https://dfeehan.github.io/surveybootstrap/reference/srs.bootstrap.sample.md)
  : srs.bootstrap.sample

## Bootstrap weight utilities

Functions that return rescaled bootstrap weights as data frames, useful
for downstream calculations such as jackknife-after-bootstrap (JAB)
variance estimation.

- [`get.rescaled.bootstrap.weights()`](https://dfeehan.github.io/surveybootstrap/reference/get.rescaled.bootstrap.weights.md)
  : rescaled.bootstrap.weights
- [`rescaled.bootstrap.weights()`](https://dfeehan.github.io/surveybootstrap/reference/rescaled.bootstrap.weights.md)
  : rescaled.bootstrap.weights

## RDS bootstrap methods

Bootstrap functions for respondent-driven sampling (RDS) data.

- [`rds.boot.draw.chain()`](https://dfeehan.github.io/surveybootstrap/reference/rds.boot.draw.chain.md)
  : Draw RDS bootstrap resamples for one chain
- [`rds.chain.boot.draws()`](https://dfeehan.github.io/surveybootstrap/reference/rds.chain.boot.draws.md)
  : Draw RDS bootstrap resamples
- [`rds.mc.boot.draws()`](https://dfeehan.github.io/surveybootstrap/reference/rds.mc.boot.draws.md)
  : Draw RDS bootstrap resamples using the algorithm in Salganik 2006
  (see Details below)
- [`estimate.degree.distns()`](https://dfeehan.github.io/surveybootstrap/reference/estimate.degree.distns.md)
  : Estimate degree distributions by trait
- [`estimate.mixing()`](https://dfeehan.github.io/surveybootstrap/reference/estimate.mixing.md)
  : Construct a mixing model from GoC/RDS data

## RDS chain utilities

Internal helper functions for working with RDS recruitment chains.

- [`chain.data()`](https://dfeehan.github.io/surveybootstrap/reference/chain.data.md)
  : Get a dataset from a chain
- [`chain.size()`](https://dfeehan.github.io/surveybootstrap/reference/chain.size.md)
  : Get the size of a chain
- [`chain.vals()`](https://dfeehan.github.io/surveybootstrap/reference/chain.vals.md)
  : Get all of the values of the given variable found among members of a
  chain
- [`is.child.ct()`](https://dfeehan.github.io/surveybootstrap/reference/is.child.ct.md)
  : Determine whether or not one id is a parent of another
- [`make.chain()`](https://dfeehan.github.io/surveybootstrap/reference/make.chain.md)
  : Build an RDS seed's chain from the dataset
- [`max(`*`<depth>`*`)`](https://dfeehan.github.io/surveybootstrap/reference/max.depth.md)
  : Get the height (maximum depth) of a chain
- [`mc.sim()`](https://dfeehan.github.io/surveybootstrap/reference/mc.sim.md)
  : Run a markov model

## Example data

The MU284 dataset (284 Swedish municipalities) and simulated samples
drawn from it, used in package examples and tests.

- [`MU284`](https://dfeehan.github.io/surveybootstrap/reference/MU284.md)
  : The MU284 Population dataset
- [`MU284.surveys`](https://dfeehan.github.io/surveybootstrap/reference/MU284.surveys.md)
  : Simulated sample surveys drawn from the MU284 Population
- [`MU284.complex.surveys`](https://dfeehan.github.io/surveybootstrap/reference/MU284.complex.surveys.md)
  : Simulated sample surveys drawn from the MU284 Population using a
  complex design
- [`MU284.boot.res.summ`](https://dfeehan.github.io/surveybootstrap/reference/MU284.boot.res.summ.md)
  : Benchmarks for unit tests
- [`MU284.estimator.fn()`](https://dfeehan.github.io/surveybootstrap/reference/MU284.estimator.fn.md)
  : MU284.estimator.fn
- [`MU284.estimator.summary.fn()`](https://dfeehan.github.io/surveybootstrap/reference/MU284.estimator.summary.fn.md)
  : MU284.estimator.summary.fn

## Legacy / internal

Older implementations retained for reference.

- [`rescaled.bootstrap.sample.pureR()`](https://dfeehan.github.io/surveybootstrap/reference/rescaled.bootstrap.sample.pureR.md)
  : rescaled.bootstrap.sample.pureR
- [`surveybootstrap-package`](https://dfeehan.github.io/surveybootstrap/reference/surveybootstrap.md)
  [`surveybootstrap`](https://dfeehan.github.io/surveybootstrap/reference/surveybootstrap.md)
  [`package-surveybootstrap`](https://dfeehan.github.io/surveybootstrap/reference/surveybootstrap.md)
  : Survey bootstrap variance estimators
