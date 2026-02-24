# Draw RDS bootstrap resamples

Draw bootstrap resamples for an RDS dataset, using the algorithm
described in the supporting online material of Weir et al 2012 (see
[`rds.boot.draw.chain()`](http://dennisfeehan.org/surveybootstrap/reference/rds.boot.draw.chain.md)
).

## Usage

``` r
rds.chain.boot.draws(chains, mm, dd, num.reps, keep.vars = NULL)
```

## Arguments

- chains:

  A list whose entries are the chains we want to resample

- mm:

  The mixing model

- dd:

  The degree distributions

- num.reps:

  The number of bootstrap resamples we want

- keep.vars:

  If not `NULL`, then the names of variables from the original dataset
  we want appended to each bootstrap resampled dataset (default is
  `NULL`)

## Value

A list of length `num.reps`; each entry in the list has one
bootstrap-resampled dataset
