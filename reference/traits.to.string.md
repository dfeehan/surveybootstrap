# Take a set of traits and turn into a string

This is a helper function that is useful when we wish to make several
traits into one variable

## Usage

``` r
traits.to.string(data, traits, na.action = "drop", sep = ".")
```

## Arguments

- data:

  The respondent info

- traits:

  The names of the traits to build the model on

- na.action:

  Defaults to 'drop' (meaning all rows of data with any missingness on
  the traits are dropped). Anything else means `NA`s are treated like
  any other value.

- sep:

  The separator character used to combine values

## Value

A list whose entries are

- `used.idx`, which indicates which rows from the original dataset were
  used (may not be all of them if there is missingness); and

- `traits`, which has the string version of the traits
