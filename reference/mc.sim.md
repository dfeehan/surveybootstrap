# Run a markov model

Run a given markov model for n time steps, starting at a specified
state.

## Usage

``` r
mc.sim(mm, start, n)
```

## Arguments

- mm:

  The markov model object returned by
  [`estimate.mixing()`](http://dennisfeehan.org/surveybootstrap/reference/estimate.mixing.md)

- start:

  The name of the state to start in

- n:

  The number of time-steps to run through

## Value

A vector with the state visited at each time step. the first entry has
the starting state

## Details

This uses the markov model produced by
[`estimate.mixing()`](http://dennisfeehan.org/surveybootstrap/reference/estimate.mixing.md)
