
#####################################################
##' srs.bootstrap.sample
##'
##' Given a survey dataset and a description of the survey
##' design (ie, which combination of vars determines primary sampling
##' units, and which combination of vars determines strata), take
##' a bunch of bootstrap samples under a simple random sampling
##' (with repetition) scheme
##'
##' @param survey.data The dataset to use
##' @param num.reps The number of bootstrap replication samples to draw
##' @param parallel If `TRUE`, use parallelization (via `plyr`)
##' @param paropts An optional list of arguments passed on to `plyr` to control
##'   details of parallelization
##' @param ... Ignored, but useful because it allows params like `survey.design`
##'   which are used in other bootstrap designs, to be passed in without error
##' @return A list with `num.reps` entries. Each entry is a dataset which has
##'   at least the variables `index` (the row index of the original dataset that
##'   was resampled) and `weight.scale` (the factor by which to multiply the
##'   sampling weights in the original dataset).
##'
##' @export
##' @examples
##'
##' survey <- MU284.surveys[[1]]
##' boot_surveys <- srs.bootstrap.sample(survey, num.reps = 2)
##'
srs.bootstrap.sample <- function(survey.data,
                                 num.reps=1,
                                 parallel=FALSE,
                                 paropts=NULL,
                                 ...)
{

  survey.data$.internal_id <- 1:nrow(survey.data)

  res <- plyr::llply(1:num.reps,
               function(rep.idx) {

                 these.samples <- sample(1:nrow(survey.data),
                                         nrow(survey.data),
                                         replace=TRUE)

                 this.rep <- data.frame(index=these.samples,
                                        weight.scale=1)

                 return(this.rep)
               },
               .parallel=parallel,
               .paropts=paropts)

  return(res)
}
