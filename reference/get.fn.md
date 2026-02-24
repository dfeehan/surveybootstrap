# Grab a function based on its name

Helper to grab a function that is passed in as an argument

## Usage

``` r
get.fn(fn, env = parent.frame())
```

## Arguments

- fn:

  The function to search for

- env:

  The environment to start searching in

## Value

`fn`, if `fn` is already a function; otherwise, the first function found
in env or one of its parents whose name is `fn`

## Details

This is based on Hadley Wickham's response to an SO post:
<https://stackoverflow.com/questions/14183766/match-fun-provide-error-with-functions-defined-inside-functions>
with some minor modifications
