# Get a variable from a dataframe or vector

This function was written because a few of the estimator functions need
to use weights, and there are several cases to handle: the user could
pass in a column name, a vector of weights, or nothing (in which case,
the weights should default to 1 for each row in the dataset). For the
special case of getting weights, look at the curried function
[`get.weights()`](http://dennisfeehan.org/surveybootstrap/reference/get.weights.md)

## Usage

``` r
get.var(survey.data, var, default = NA)
```

## Arguments

- survey.data:

  The survey dataset

- var:

  Either `NULL`, a column name, or a vector of values

- default:

  The default value to fill in if the variable is not found

## Value

A vector of values whose length is the same as the number of rows in
`survey.data`; if `var` is `NULL`, this has the default values
