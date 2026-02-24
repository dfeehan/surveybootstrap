# Draw RDS bootstrap resamples using the algorithm in Salganik 2006 (see Details below)

This algorithm picks a respondent from the survey to be a seed uniformly
at random. it then generates a bootstrap draw by simulating the markov
process forward for n steps, where n is the size of the draw required.

If you wish the bootstrap dataset to end up with variables from the
original dataset other than the traits and degree, then you must specify
this when you construct `dd` using the '`estimate.degree.distns`
function.

## Usage

``` r
rds.mc.boot.draws(chains, mm, dd, num.reps)
```

## Arguments

- chains:

  A list with the chains constructed from the survey using `make.chain`

- mm:

  The mixing model

- dd:

  The degree distributions

- num.reps:

  The number of bootstrap resamples we want

## Value

A list of length `num.reps`; each entry in the list has one
bootstrap-resampled dataset

## Details

See:

- Salganik, Matthew J. "Variance estimation, design effects, and sample
  size calculations for respondent-driven sampling." *Journal of Urban
  Health* 83.1 (2006): 98-112.
