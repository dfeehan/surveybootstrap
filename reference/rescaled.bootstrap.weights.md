# rescaled.bootstrap.weights

This function creates a dataset with rescaled bootstrap weights; it can
be a helpful alternative to `bootstrap.estimates` in some situations

## Usage

``` r
rescaled.bootstrap.weights(
  survey.data,
  survey.design,
  num.reps,
  weights = NULL,
  idvar,
  verbose = TRUE,
  parallel = FALSE,
  paropts = NULL
)
```

## Arguments

- survey.data:

  The dataset to use

- survey.design:

  A formula describing the design of the survey (see Details in
  [`bootstrap.estimates()`](http://dennisfeehan.org/surveybootstrap/reference/bootstrap.estimates.md)
  help page)

- num.reps:

  the number of bootstrap replication samples to draw

- weights:

  weights to use in estimation (or NULL, if none)

- idvar:

  the name of the column in `survey.data` that has the respondent id

- verbose:

  if TRUE, produce lots of feedback about what is going on

- parallel:

  if TRUE, use the plyr library's .parallel argument to produce
  bootstrap resamples and estimates in parallel

- paropts:

  if not NULL, additional arguments to pass along to the parallelization
  routine

## Value

if no summary.fn is specified, then return the list of estimates
produced by estimator.fn; if summary.fn is specified, then return its
output

## Details

The formula describing the survey design should have the form
`~ psu_v1 + psu_v2 + ... + strata(strata_v1 + strata_v2 + ...)`, where
`psu_v1, ...` are the variables identifying primary sampling units
(PSUs) and `strata_v1, ...` identify the strata

## Examples

``` r
survey <- MU284.complex.surveys[[1]]
rescaled.bootstrap.weights(survey.data = survey,
                          survey.design = ~ CL,
                          weights='sample_weight',
                          idvar='LABEL',
                          num.reps = 2)
#>    LABEL boot_weight_1 boot_weight_2
#> 1     96       0.00000      33.33333
#> 2     99       0.00000      33.33333
#> 3    101       0.00000      33.33333
#> 4     16      41.66667      20.83333
#> 5     18      41.66667      20.83333
#> 6     19      41.66667      20.83333
#> 7    221      41.66667       0.00000
#> 8    223      41.66667       0.00000
#> 9    224      41.66667       0.00000
#> 10   188       0.00000      25.00000
#> 11   190       0.00000      25.00000
#> 12   191       0.00000      25.00000
#> 13   199       0.00000      25.00000
#> 14   200       0.00000      25.00000
#> 15   203       0.00000      25.00000


if (FALSE) { # \dontrun{
bootweights <- rescaled.bootstrap.weights(
                                         # formula describing survey design:
                                         # psu and strata
                                         survey.design = ~ psu +
                                                           stratum(stratum_analysis),
                                         num.reps=10000,
                                         # column with respondent ids
                                         idvar='caseid',
                                         # column with sampling weight
                                         weights='wwgt',
                                         # survey dataset
                                         survey.data=mw.ego)

} # }
```
