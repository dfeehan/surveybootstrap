<div id="main" class="col-md-9" role="main">

Construct a mixing model from GoC/RDS data
==========================================

<div class="ref-description section level2">

Given a dataset with the respondents and a dataset on the parents (in
many cases the same individuals), and a set of relevant traits, estimate
mixing parameters and return a markov model.

</div>

<div class="section level2">

Usage
-----

<div class="sourceCode">

``` r
estimate.mixing(survey.data, parent.data, traits)
```

</div>

</div>

<div class="section level2">

Arguments
---------

-   survey.data:

    The respondent info

-   parent.data:

    The parent info

-   traits:

    The names of the traits to build the model on

</div>

<div class="section level2">

Value
-----

A list with entries:

-   `mixing.df` the data used to estimate the mixing function

-   `choose.next.state.fn` a function which can be passed a vector of
    states and will return a draw of a subsequent state for each entry
    in the vector

-   `mixing.df` a dataframe (long-form) representation of the transition
    counts used to estimate the transition probabilities

-   `states` a list with an entry for each state. within each state's
    entry are

    -   `trans.probs` a vector of estimated transition probabilities

    -   `trans.fn` a function which, when called, randomly chooses a
        next state with probabilities given by the transition probs.

</div>

</div>
