# Parse a formula that describes the design of a survey

Parse a formula of the form
`~ psu_v1 + psu_v2 + ... + strata(strata_v1 + strata_v2 + ...)` into a
PSU formula and a strata formula. Note that we only document
`strata(...)` but this will also work with `stratum(...)` because it is
annoying to have to remember which one to use

## Usage

``` r
parse_design(formula)
```

## Arguments

- formula:

  a formula describing the sample design (see Description of
  [`bootstrap.estimates()`](https://dfeehan.github.io/surveybootstrap/reference/bootstrap.estimates.md))

## Value

a list with entries `psu.formula` and `strata.formula`
