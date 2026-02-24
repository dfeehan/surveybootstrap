# Get all of the values of the given variable found among members of a chain

Get all of the values of the given variable found among members of a
chain

## Usage

``` r
chain.vals(chain, qoi.var = "uid")
```

## Arguments

- chain:

  The chain to get values from

- qoi.var:

  The name of the variable to get from each member of the chain

## Value

A vector with all of the values of `qoi.var` found in this chain.
(Currently, the order of the values in the vector is not guaranteed.)
