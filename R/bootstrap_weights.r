#####################################################
##' rescaled.bootstrap.weights
##'
##' this helper function creates a dataset with rescaled bootstrap weights;
##' it can be a helpful alternative to \code{bootstrap.estimates} in some situations\cr
##'
##' @param survey.data the dataset to use
##' @param survey.design a formula describing the design of the survey
##'                      (see below - TODO)
##' @param num.reps the number of bootstrap replication samples to draw
##' @param weights weights to use in estimation (or NULL, if none)
##' @param idvar the name of the column in \code{survey.data} that has the respondent id
##' @param parallel if TRUE, use the plyr library's .parallel argument to
##'                 produce bootstrap resamples and estimates in parallel
##' @param paropts if not NULL, additional arguments to pass along to the
##'                parallelization routine
##' @param verbose if TRUE, produce lots of feedback about what is going on
##' @return if no summary.fn is specified, then return the list of estimates
##'         produced by estimator.fn; if summary.fn is specified, then return
##'         its output
##' @export
##' @examples
##' \donttest{
##' # code goes here
##'
##' }
rescaled.bootstrap.weights <- function(survey.data,
                                       survey.design,
                                       num.reps,
                                       weights=NULL,
                                       idvar,
                                       verbose=TRUE,
                                       parallel=FALSE,
                                       paropts=NULL)
{

  #
  boot_weights <- function(survey.data, weights) {

    res <- survey.data %>% select(one_of(idvar))
    res$.bootstrapweights <- weights

    return(list(boot.data=res))
  }

  raw.res <- bootstrap.estimates(survey.design=survey.design,
                                 num.reps=num.reps,
                                 weights=weights,
                                 estimator.fn=boot_weights,
                                 bootstrap.fn='rescaled.bootstrap.sample',
                                 survey.data=survey.data)

  # for the tibble generated from the ith boostrap replicate, rename the weight column
  # 'boot_weight_i' (so that the columns each have distinct names)
  res <- purrr::imap(raw.res,
                     ~.x$boot.data %>% rename(!!paste0('boot_weight_', .y) := .bootstrapweights))

  # take the list of tibbles and turn it into one big tibble
  # NB: I think this could be faster if we bind_cols instead of using left_join.
  #     But that approach relies on the row order being stable; I think it is now, but it seems like
  #     the code will be more fragile if we make that assumption here
  #dfres <- res %>% purrr::reduce(left_join, by='caseid')
  dfres <- res %>% purrr::reduce(left_join, by=idvar)

  return(dfres)
}
