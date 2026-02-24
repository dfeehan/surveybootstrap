# Package index

## Main bootstrap interface

High-level function for obtaining bootstrap estimates of sampling
uncertainty from any estimator.

- [`bootstrap.estimates()`](http://dennisfeehan.org/surveybootstrap/reference/bootstrap.estimates.md)
  : bootstrap.estimates

## Bootstrap sampling methods

Functions that draw bootstrap resamples from a survey dataset. Pass the
name of one of these to the `bootstrap.fn` argument of
[`bootstrap.estimates()`](http://dennisfeehan.org/surveybootstrap/reference/bootstrap.estimates.md).

- [`rescaled.bootstrap.sample()`](http://dennisfeehan.org/surveybootstrap/reference/rescaled.bootstrap.sample.md)
  : rescaled.bootstrap.sample
- [`srs.bootstrap.sample()`](http://dennisfeehan.org/surveybootstrap/reference/srs.bootstrap.sample.md)
  : srs.bootstrap.sample

## Bootstrap weight utilities

Functions that return rescaled bootstrap weights as data frames, useful
for downstream calculations such as jackknife-after-bootstrap (JAB)
variance estimation.

- [`get.rescaled.bootstrap.weights()`](http://dennisfeehan.org/surveybootstrap/reference/get.rescaled.bootstrap.weights.md)
  : rescaled.bootstrap.weights
- [`rescaled.bootstrap.weights()`](http://dennisfeehan.org/surveybootstrap/reference/rescaled.bootstrap.weights.md)
  : rescaled.bootstrap.weights

## RDS bootstrap methods

Bootstrap functions for respondent-driven sampling (RDS) data.

- [`rds.boot.draw.chain()`](http://dennisfeehan.org/surveybootstrap/reference/rds.boot.draw.chain.md)
  : Draw RDS bootstrap resamples for one chain
- [`rds.chain.boot.draws()`](http://dennisfeehan.org/surveybootstrap/reference/rds.chain.boot.draws.md)
  : Draw RDS bootstrap resamples
- [`rds.mc.boot.draws()`](http://dennisfeehan.org/surveybootstrap/reference/rds.mc.boot.draws.md)
  : Draw RDS bootstrap resamples using the algorithm in Salganik 2006
  (see Details below)
- [`estimate.degree.distns()`](http://dennisfeehan.org/surveybootstrap/reference/estimate.degree.distns.md)
  : Estimate degree distributions by trait
- [`estimate.mixing()`](http://dennisfeehan.org/surveybootstrap/reference/estimate.mixing.md)
  : Construct a mixing model from GoC/RDS data

## RDS chain utilities

Internal helper functions for working with RDS recruitment chains.

- [`chain.data()`](http://dennisfeehan.org/surveybootstrap/reference/chain.data.md)
  : Get a dataset from a chain
- [`chain.size()`](http://dennisfeehan.org/surveybootstrap/reference/chain.size.md)
  : Get the size of a chain
- [`chain.vals()`](http://dennisfeehan.org/surveybootstrap/reference/chain.vals.md)
  : Get all of the values of the given variable found among members of a
  chain
- [`is.child.ct()`](http://dennisfeehan.org/surveybootstrap/reference/is.child.ct.md)
  : Determine whether or not one id is a parent of another
- [`make.chain()`](http://dennisfeehan.org/surveybootstrap/reference/make.chain.md)
  : Build an RDS seed's chain from the dataset
- [`max(`*`<depth>`*`)`](http://dennisfeehan.org/surveybootstrap/reference/max.depth.md)
  : Get the height (maximum depth) of a chain
- [`mc.sim()`](http://dennisfeehan.org/surveybootstrap/reference/mc.sim.md)
  : Run a markov model

## Example data

The MU284 dataset (284 Swedish municipalities) and simulated samples
drawn from it, used in package examples and tests.

- [`MU284`](http://dennisfeehan.org/surveybootstrap/reference/MU284.md)
  : The MU284 Population dataset
- [`MU284.surveys`](http://dennisfeehan.org/surveybootstrap/reference/MU284.surveys.md)
  : Simulated sample surveys drawn from the MU284 Population
- [`MU284.complex.surveys`](http://dennisfeehan.org/surveybootstrap/reference/MU284.complex.surveys.md)
  : Simulated sample surveys drawn from the MU284 Population using a
  complex design
- [`MU284.boot.res.summ`](http://dennisfeehan.org/surveybootstrap/reference/MU284.boot.res.summ.md)
  : Benchmarks for unit tests
- [`MU284.estimator.fn()`](http://dennisfeehan.org/surveybootstrap/reference/MU284.estimator.fn.md)
  : MU284.estimator.fn
- [`MU284.estimator.summary.fn()`](http://dennisfeehan.org/surveybootstrap/reference/MU284.estimator.summary.fn.md)
  : MU284.estimator.summary.fn

## Legacy / internal

Older implementations retained for reference.

- [`rescaled.bootstrap.sample.pureR()`](http://dennisfeehan.org/surveybootstrap/reference/rescaled.bootstrap.sample.pureR.md)
  : rescaled.bootstrap.sample.pureR
- [`surveybootstrap-package`](http://dennisfeehan.org/surveybootstrap/reference/surveybootstrap.md)
  [`surveybootstrap`](http://dennisfeehan.org/surveybootstrap/reference/surveybootstrap.md)
  [`package-surveybootstrap`](http://dennisfeehan.org/surveybootstrap/reference/surveybootstrap.md)
  : Survey bootstrap variance estimators
