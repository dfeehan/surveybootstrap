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
#'   \item{ME84}{Number of municipal employees}
#'   \item{REV84}{Real estate values according to 1984 assessment}
#'   \item{REG}{Geographic location indicator}
#'   \item{CL}{Cluster indicator (neighboring municipalities are clustered together)}
#' }
#' @source 'Model Assisted Survey Sampling' by Sarndal, Swensson, and Wretman (2003, ISBN:0387406204)
NULL

#' Simulated sample surveys drawn from the MU284 Population
#'
#' A list with 15 sample surveys drawn from the [MU284] dataset.
#'
#' @name MU284.surveys
#'
#' @format
#' A data frame with 284 rows and 11 columns:
#' \describe{
#'   \item{LABEL, ..., CL}{Same as MU284 dataset}
#'   \item{sample_weight}{The sampling weight for the row}
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
                         TS82.hat=sum(S82*weight),
                         R.RMT85.P85.hat=sum(RMT85*weight)/sum(P85*weight))
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
