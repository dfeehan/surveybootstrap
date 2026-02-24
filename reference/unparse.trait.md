# Unparse a collapsed trait string

For a few of the RDS-related functions, it is useful to combine several
traits into one variable as a string; for example, "male" and "young"
might become "male.young". this function takes a string with combined
traits and explodes it back into several variables.

## Usage

``` r
unparse.trait(trait.string, names, sep = "\\.")
```

## Arguments

- trait.string:

  A vector whose values are collapsed traits

- names:

  A vector with the names of each trait (in order)

- sep:

  The character used to separate the traits in their collapsed string
  representation

## Value

A dataframe whose rows correspond to the entries in `trait.string`, with
one column per trait
