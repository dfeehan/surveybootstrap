<div id="main" class="col-md-9" role="main">

Package index
=============

<div class="section level2">

Main bootstrap interface
------------------------

<div class="section-desc">

High-level function for obtaining bootstrap estimates of sampling
uncertainty from any estimator.

</div>

</div>

<div class="section level2">

-   `bootstrap.estimates()` : bootstrap.estimates

</div>

<div class="section level2">

Bootstrap sampling methods
--------------------------

<div class="section-desc">

Functions that draw bootstrap resamples from a survey dataset. Pass the
name of one of these to the `bootstrap.fn` argument of
`bootstrap.estimates()`.

</div>

</div>

<div class="section level2">

-   `rescaled.bootstrap.sample()` : rescaled.bootstrap.sample
-   `srs.bootstrap.sample()` : srs.bootstrap.sample

</div>

<div class="section level2">

Bootstrap weight utilities
--------------------------

<div class="section-desc">

Functions that return rescaled bootstrap weights as data frames, useful
for downstream calculations such as jackknife-after-bootstrap (JAB)
variance estimation.

</div>

</div>

<div class="section level2">

-   `get.rescaled.bootstrap.weights()` : rescaled.bootstrap.weights
-   `rescaled.bootstrap.weights()` : rescaled.bootstrap.weights

</div>

<div class="section level2">

RDS bootstrap methods
---------------------

<div class="section-desc">

Bootstrap functions for respondent-driven sampling (RDS) data.

</div>

</div>

<div class="section level2">

-   `rds.boot.draw.chain()` : Draw RDS bootstrap resamples for one chain
-   `rds.chain.boot.draws()` : Draw RDS bootstrap resamples
-   `rds.mc.boot.draws()` : Draw RDS bootstrap resamples using the
    algorithm in Salganik 2006 (see Details below)
-   `estimate.degree.distns()` : Estimate degree distributions by trait
-   `estimate.mixing()` : Construct a mixing model from GoC/RDS data

</div>

<div class="section level2">

RDS chain utilities
-------------------

<div class="section-desc">

Internal helper functions for working with RDS recruitment chains.

</div>

</div>

<div class="section level2">

-   `chain.data()` : Get a dataset from a chain
-   `chain.size()` : Get the size of a chain
-   `chain.vals()` : Get all of the values of the given variable found
    among members of a chain
-   `is.child.ct()` : Determine whether or not one id is a parent of
    another
-   `make.chain()` : Build an RDS seed's chain from the dataset
-   `max(<depth>)` : Get the height (maximum depth) of a chain
-   `mc.sim()` : Run a markov model

</div>

<div class="section level2">

Example data
------------

<div class="section-desc">

The MU284 dataset (284 Swedish municipalities) and simulated samples
drawn from it, used in package examples and tests.

</div>

</div>

<div class="section level2">

-   `MU284` : The MU284 Population dataset
-   `MU284.surveys` : Simulated sample surveys drawn from the MU284
    Population
-   `MU284.complex.surveys` : Simulated sample surveys drawn from the
    MU284 Population using a complex design
-   `MU284.boot.res.summ` : Benchmarks for unit tests
-   `MU284.estimator.fn()` : MU284.estimator.fn
-   `MU284.estimator.summary.fn()` : MU284.estimator.summary.fn

</div>

<div class="section level2">

Legacy / internal
-----------------

<div class="section-desc">

Older implementations retained for reference.

</div>

</div>

<div class="section level2">

-   `rescaled.bootstrap.sample.pureR()` :
    rescaled.bootstrap.sample.pureR
-   `surveybootstrap-package` `surveybootstrap`
    `package-surveybootstrap` : Survey bootstrap variance estimators

</div>

</div>
