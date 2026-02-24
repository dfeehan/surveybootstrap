# Build an RDS seed's chain from the dataset

Note that this assumes that the chain is a tree (no loops)

## Usage

``` r
make.chain(seed.id, survey.data, is.child.fn = is.child.ct)
```

## Arguments

- seed.id:

  The id of the seed whose chain we wish to build from the dataset

- survey.data:

  The dataset

- is.child.fn:

  A function which takes two ids as arguments; it is expected to return
  `TRUE` if the second argument is the parent of the first, and `FALSE`
  otherwise. it defaults to
  [`is.child.ct()`](http://dennisfeehan.org/surveybootstrap/reference/is.child.ct.md)

## Value

info
