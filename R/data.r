#' The MU284 Population dataset
#'
#' A dataset containing information about Sweden's 284 municipalities. This
#' dataset comes from Model-Assisted Survey Sampling by
#' Sarndal, Swensson, and Wretman (2003, ISBN:0387406204). The columns are:
#'
#' @name MU284
#'
#' @format
#' A data frame with 284 rows and 11 columns:
#' \describe{
#'   \item{LABEL}{ID}
#'   \item{P85}{Population in 1985}
#'   \item{P75}{Population in 1975}
#'   \item{RMT85}{Municipal tax revenue in 1985}
#'   \item{CS82}{Number of Conservative seats in municipal council}
#'   \item{SS82}{Number of Social-Democratic seats in municipal council}
#'   \item{S82}{Total number of seats in municipal council}
#'   \item{ME84}{Number of municipal employees}
#'   \item{REV84}{Real estate values according to 1984 assessment}
#'   \item{REG}{Geographic location indicator}
#'   \item{CL}{Cluster indicator (neighboring municipalities are clustered together)}
#' }
#' @source 'Model Assisted Survey Sampling' by Sarndal, Swensson, and Wretman (2003, ISBN:0387406204)
NULL

#' Simulated sample surveys drawn from the MU284 Population
#'
#' A list with 10 sample surveys with sample size 15 drawn from the [MU284]
#' dataset using simple random sampling with replacment.
#'
#' @name MU284.surveys
#'
#' @format
#' A list with 10 data frames, each with 15 rows and 11 columns:
#' \describe{
#'   \item{LABEL, ..., CL}{Same as MU284 dataset}
#'   \item{sample_weight}{The sampling weight for the row}
#' }
NULL

#' Simulated sample surveys drawn from the MU284 Population using a complex design
#'
#' A list with 10 sample surveys with sample size 15 drawn from the [MU284]
#' dataset using a complex sampling design.
#'
#' @details
#' The sampling design comes from Ex. 4.3.2 (pg 142-3) of
#' 'Model Assisted Survey Sampling' by Sarndal, Swensson, and Wretman
#' (2003, ISBN:0387406204).
#'
#' The design is a two-stage sample:
#'  * stage I: the primary sampling units (PSUs) are the standard clusters from
#'    [MU284]; we take a simple random sample without replacement of `n_I = 5`
#'    out of `N_I = 50` of these
#'  * stage II: within each sampled PSU, we take a simple random sample without
#'    replacement of `n_i = 3` out of `N_i` municipalities
#'
#' @name MU284.complex.surveys
#'
#' @format
#' A list with 10 data frames, each with 15 rows and 11 columns:
#' \describe{
#'   \item{LABEL, ..., CL}{Same as MU284 dataset}
#'   \item{sample_weight}{The sampling weight for the row}
#' }
NULL

#' Benchmarks for unit tests
#'
#' Benchmark results to use in unit tests; these are based on
#' [MU284.complex.surveys].
#'
#' @name MU284.boot.res.summ
#'
#' @format
#' A list with 10 data frames, each with 15 rows and 11 columns:
#' \describe{
#'   \item{mean.TS82.hat, ..., sd.R.RMT85.P85.hat}{summaries for each estimand}
#' }
NULL

#' MU284.estimator.fn
#'
#' Produce estimates from a simulated sample survey of the [MU284] population.
#' Used in package tests and examples.
#'
#' @param survey.data the survey dataset
#' @param weights a vector with the survey weights
#' @return a data.frame with one row and two columns:
#'   * `TS82.hat` - the estimated total of `S82`
#'   * `R.RMT85.P85.hat` - the estimated ratio of `RMT85` / `P85`
MU284.estimator.fn <- function(survey.data, weights) {
  survey.data$weight <- weights

  res <- plyr::summarise(survey.data,
                         TS82.hat=sum(survey.data$S82*survey.data$weight),
                         R.RMT85.P85.hat=sum(survey.data$RMT85*survey.data$weight)/
                                         sum(survey.data$P85*survey.data$weight))
  return(res)
}

#' MU284.estimator.summary.fn
#'
#' Summarize results from [MU284.estimator.fn()] applied to many surveys.
#' (This is a dummy function, used for tests)
#'
#' @param res a dataframe whose rows are the results of calling [MU284.estimator.fn()]
#' @return the same dataframe
MU284.estimator.summary.fn <- function(res) {
  plyr::ldply(res, identity)
}
