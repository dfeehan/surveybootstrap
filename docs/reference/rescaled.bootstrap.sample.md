# rescaled.bootstrap.sample

Given a survey dataset and a description of the survey design (ie, which
combination of variables determines primary sampling units, and which
combination of variables determines strata), take a bunch of bootstrap
samples for the rescaled bootstrap estimator (see Details).

## Usage

``` r
rescaled.bootstrap.sample(
  survey.data,
  survey.design,
  parallel = FALSE,
  paropts = NULL,
  num.reps = 1,
  include_cc = FALSE
)
```

## Arguments

- survey.data:

  The dataset to use

- survey.design:

  A formula describing the design of the survey (see Details)

- parallel:

  If `TRUE`, use parallelization (via `plyr`)

- paropts:

  An optional list of arguments passed on to `plyr` to control details
  of parallelization

- num.reps:

  The number of bootstrap replication samples to draw

- include_cc:

  If `TRUE`, include cluster counts in the output. These cluster counts
  are needed for jackknife after bootstrap calculations

## Value

    The return value will depend on whether or not you ask for the cluster counts (using `include_cc`).
    If `include_cc` is `FALSE` (the default), return
      A list with `num.reps` entries. Each entry is a dataset which
      has at least the variables `index` (the row index of the original
      dataset that was resampled) and `weight.scale`
      (the factor by which to multiply the sampling weights
      in the original dataset).

     If `include_cc` is `TRUE`, then return a list with two entries:

    * `weight_scaling_factor` - A list with `num.reps` entries.
       Each entry is a dataset which
      has at least the variables `index` (the row index of the original
      dataset that was resampled) and `weight.scale`
      (the factor by which to multiply the sampling weights
      in the original dataset).
    * `cluster_counts` - A list with `num.reps` entries.
      Each entry is a dataset which has the variable '.cluster_id',
      the values of columns that specify the PSU (from the `survey.design` formula),
      and `cluster_count` (the number of times the given PSU was resampled))

## Details

`survey.design` is a formula of the form

`weight ~ psu_vars + strata(strata_vars)`

where:

- `weight` is the variable with the survey weights

- `psu_vars` has the form `psu_v1 + psu_v2 + ...`, where primary
  sampling units (PSUs) are determined by `psu_v1`, etc

- `strata_vars` has the form `strata_v1 + strata_v2 + ...`, which
  determine strata

Note that we assume that the formula uniquely specifies PSUs. This will
always be true if the PSUs were selected without replacement. If they
were selected with replacement, then it will be necessary to make each
realization of a given PSU in the sample a unique id. The code below
assumes that all observations within each PSU (as identified by the
design formula) are from the same draw of the PSU.

The rescaled bootstrap technique works by adjusting the estimation
weights based on the number of times each row is included in the
resamples. If a row is never selected, it is still included in the
returned results, but its weight will be set to 0. It is therefore
important to use estimators that make use of the estimation weights on
the resampled datasets.

We always take m_i = n_i - 1, according to the advice presented in Rao
and Wu (1988) and Rust and Rao (1996).

(This is a C++ version; a previous version, written in pure R, is called
[`rescaled.bootstrap.sample.pureR()`](https://dfeehan.github.io/surveybootstrap/reference/rescaled.bootstrap.sample.pureR.md)
)

References:

- Rust, Keith F., and J. N. K. Rao. "Variance estimation for complex
  surveys using replication techniques." *Statistical methods in medical
  research* 5.3 (1996): 283-310.

- Rao, Jon NK, and C. F. J. Wu. "Resampling inference with complex
  survey data." *Journal of the American Statistical Association* 83.401
  (1988): 231-241.

## Examples

``` r
survey <- MU284.complex.surveys[[1]]
boot_surveys <- rescaled.bootstrap.sample(survey.data = survey,
                                          survey.design = ~ CL,
                                          num.reps = 2)
```
