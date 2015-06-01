#####################################################
## variance_estimators.R
##
## this file has the functions that produce
## variance estimates

#####################################################
##' killworth.se
##'
##' compute standard errors for scale-up estimates
##' based on the Killworth estimator
##'
##' note that this is provided for comparison, but
##' that we do not generally recommend using this
##' strategy for estimating variance
##'
##' @param estimates TODO
##' @param d.hat TODO
##' @param total.popn.size TODO
##' @param total TODO
##' @param missing TODO
##' @return the estimated standard error
##' @keywords internal
killworth.se <- function(estimates,
                         d.hat,
                         total.popn.size=NULL,
                         total=TRUE,
                         missing="ignore") {

  stop("killworth.se is not yet implemented.")

  ## TODO -- code below this point not yet altered...

  if (total & is.null(total.popn.size)) {
    stop("you must pass in total.popn.size to get the Killworth variance estimates of the total (rather than proportions).")
  }

  na.rm <- ifelse(missing == "ignore", TRUE, FALSE)

  sum.d.hat <- sum(d.hat, na.rm=na.rm)

  est.props <- ifelse(rep(total, length(estimates)),
                      estimates / total.popn.size,
                      estimates)

  res <- aaply(est.props,
               1,
               function(est) {
                 return(sqrt((est*(1-est))/sum.d.hat))
               })

  names(res) <- names(estimates)

  if (total) {
    res <- res*total.popn.size
  }

  return(res)
}


#####################################################
##' bootstrap.estimates
##'
##' this function contains the core of the rescaled bootstrap
##' method for estimating uncertainty in our estimates
##' it should be designed so that it can be passed in to
##' estimation functions as an argument\cr
##' OR\cr
##' \cr
##' @section TODO:
##' \itemize{
##'   \item{estimator.fn/bootstrap.fn and summary.fn are treated differently
##'         (one expects characters, one expects an actual fn. fix!)}
##'   \item{write description block, including estimator.fn, bootstrap.fn,
##'         summary.fn, more?}
##'}
##'
##' @param survey.data the dataset to use
##' @param survey.design a formula describing the design of the survey
##'                      (see below - TODO)
##' @param estimator.fn name of a function which, given a dataset like
##'                     survey data and arguments in \code{...},
##'                     will produce an estimate of interest
##' @param bootstrap.fn name of the method to be used to take
##'                     bootstrap resamples; see below
##' @param num.reps the number of bootstrap replication samples to draw
##' @param weights weights to use in estimation (or NULL, if none)
##' @param summary.fn (optional) name of a function which, given the set of estimates
##'                   produced by estimator.fn, summarizes them. if not specified, all of
##'                   the estimates are returned in a list
##' @param parallel if TRUE, use the plyr library's .parallel argument to
##'                 produce bootstrap resamples and estimates in parallel
##' @param paropts if not NULL, additional arguments to pass along to the
##'                parallelization routine
##' @param verbose if TRUE, produce lots of feedback about what is going on
##' @param ... additional arguments which will be passed on to the estimator fn
##' @return if no summary.fn is specified, then return the list of estimates
##'         produced by estimator.fn; if summary.fn is specified, then return
##'         its output
##' @export
##' @examples
##' \donttest{
##' # code goes here
##' ...
##' }
bootstrap.estimates <- function(survey.data,
                                survey.design,
                                bootstrap.fn,
                                estimator.fn,
                                num.reps,
                                weights=NULL,
                                ...,
                                summary.fn=NULL,
                                verbose=TRUE,
                                parallel=FALSE,
                                paropts=NULL)
{

  ## get the weights
  weights <- get.weights(survey.data, weights)

  ## build up a single call to obtain an actual bootstrap
  ## replicate; we'll call this once for each one...
  boot.call <- match.call()

  boot.arg.idx <- match(c("survey.data",
                          "survey.design",
                          "num.reps",
                          "parallel",
                          "paropts"),
                        names(boot.call),
                        0L)
  boot.call <- boot.call[c(1,boot.arg.idx)]
  
  boot.call[[1]] <- get.fn(bootstrap.fn)

  ## also build up a call to obtain an estimate from the data
  est.call <- match.call(expand.dots=TRUE)

  ## these are the args we *won't* use when we call the estimator
  ## (ie, we use them here or in the bootstrap fn instead)
  est.arg.idx <- match(c("survey.design",
                         "estimator.fn",
                         "num.reps",
                         ##"weights",
                         "summary.fn",
                         "bootstrap.fn",
                         "parallel",
                         "paropts",
                         ## don't use survey.data b/c we're going to pass
                         ## in a bootstrap resample each time instead...
                         "survey.data"),
                       names(est.call),
                       0L)
  est.call <- est.call[-est.arg.idx]
  est.call[[1]] <- get.fn(estimator.fn)

  ## get the bootstrap samples
  boot.idx <- eval(boot.call, parent.frame())

  ## produce our estimate for each one
  res <- llply(boot.idx,

               function(this.rep) {

                 ## use the resampled indices to construct
                 ## a full resampled dataset
                 tmpdat <- survey.data[this.rep$index,]
                 tmpweights <- weights[this.rep$index]

                 ## apply the weight.scale to the estimation weights
                 tmpweights <- tmpweights * this.rep$weight.scale

                 ## add the information about which rows in the individual dataset
                 ## these resamples come from as an attribute
                 attr(tmpdat, "resampled.rows.orig.idx") <- this.rep$index

                 est.call[["survey.data"]] <- tmpdat
                 est.call[["weights"]] <- tmpweights

                 ## call estimator.fn to produce an estimate from
                 ## the bootstrap-resampled dataset
                 this.est <- eval(est.call, parent.frame(2))

                 return(this.est)
               },
               .parallel=parallel,
               .paropts=paropts)

  ## if the user specified a summary function, use it
  if (! is.null(summary.fn)) {

    this.sf <- get.fn(summary.fn)
    res <- do.call(this.sf, list(res))

  }

  return(res)

}

#####################################################
##' rescaled.bootstrap.sample.cpp
##'
##' C++ version: given a survey dataset and a description of the survey
##' design (ie, which combination of vars determines primary sampling
##' units, and which combination of vars determines strata), take
##' a bunch of bootstrap samples for the rescaled bootstrap estimator
##' (see, eg, Rust and Rao 1996).
##'
##' Note that we assume that the formula uniquely specifies PSUs.
##' This will always be true if the PSUs were selected without replacement.
##' If they were selected with replacement, then it will be necessary
##' to make each realization of a given PSU in the sample a unique id.
##' Bottom line: the code below assumes that all observations within
##' each PSU (as identified by the design formula) are from the same draw
##' of the PSU.
##'
##' The rescaled bootstrap technique works by adjusting the
##' estimation weights based on the number of times each
##' row is included in the resamples. If a row is never selected,
##' it is still included in the returned results, but its weight
##' will be set to 0. It is therefore important to use estimators
##' that make use of the estimation weights on the resampled
##' datasets.
##'
##' We always take m_i = n_i - 1, according to the advice presented
##' in Rao and Wu (1988) and Rust and Rao (1996).
##'
##' @param survey.data the dataset to use
##' @param survey.design a formula describing the design of the survey (see below - TODO)
##' @param num.reps the number of bootstrap replication samples to draw
##' @param parallel if TRUE, use parallelization (via \code{plyr})
##' @param paropts an optional list of arguments passed on to \code{plyr} to control
##'        details of parallelization
##' @return a list with \code{num.reps} entries. each entry is a dataset which
##' has at least the variables \code{index} (the row index of the original
##' dataset that was resampled) and \code{weight.scale}
##' (the factor by which to multiply the sampling weights
##' in the original dataset).
##' @details \code{survey.design} is a formula of the form\cr
##'    weight ~ psu_vars + strata(strata_vars),
##' where weight is the variable with the survey weights and psu
##' is the variable denoting the primary sampling unit
##' @export
rescaled.bootstrap.sample <- function(survey.data,
                                      survey.design,
                                      parallel=FALSE,
                                      paropts=NULL,
                                      num.reps=1)
{

  survey.data$.internal_id <- 1:nrow(survey.data)

  design <- parse_design(survey.design)

  ## drop the "~" at the start of the formula
  psu.vars <- design$psu.formula[c(-1)][[1]]

  ## in the special case where there are no PSU vars, treat each row as
  ## its own PSU
  if (length(psu.vars)==1 & psu.vars=="1") {
    psu.vars <- as.name(".internal_id")
  }

  ## create a single variable with an id number for each PSU
  ## (we need this to use the C++ code, below)
  ## TODO - check this more; I believe it works
  survey.data$.cluster_id <- group_indices_(survey.data, .dots=all.vars(psu.vars))

  ## if no strata are specified, enclose the entire survey all in
  ## one stratum
  if (is.null(design$strata.formula)) {
    strata <- list(survey.data)
  } else {
    strata <- dlply(survey.data, design$strata.formula, identity)
  }

  ## get num.reps bootstrap resamples within each stratum,
  ## according to the rescaled bootstrap scheme
  ## (see, eg, Rust and Rao 1996)

  ## this llply call returns a list, with one entry for each stratum
  ## each stratum's entry contains a list with the bootstrap resamples
  ## (see the note for the inner llply call below)
  bs <- llply(strata,
              function(stratum.data) {

                ## (this part is written in c++)
                res <- resample_stratum(stratum.data$.cluster_id,
                                        num.reps)

                colnames(res) <- paste0("rep.", 1:ncol(res))
                res <- cbind("index"=stratum.data$.internal_id,
                             res)
                
                return(res)
              })

  ## bs: list, one entry for each stratum
  ## each stratum's entry is a matrix.
  ## first column of the matrix is called 'index',
  ## which is the row number for the observation in the
  ## original dataset; there is one remaining column for each
  ## bootstrap resample. the entries of each column are the factors
  ## by which the original weights should be scaled

  bs.all <- do.call("rbind", bs)

  res <- alply(bs.all[,-1],
               2,
               function(this_col) {
                   return(data.frame(index=bs.all[,1],
                                     weight.scale=this_col))
               })

  return(res)

}

#####################################################
##' rescaled.bootstrap.sample.pureR
##'
##' (this is the pure R version; it has been supplanted by
##'  \code{rescaled.bootstrap.sample}, which is partially written in C++)
##' 
##' given a survey dataset and a description of the survey
##' design (ie, which combination of vars determines primary sampling
##' units, and which combination of vars determines strata), take
##' a bunch of bootstrap samples for the rescaled bootstrap estimator
##' (see, eg, Rust and Rao 1996).
##'
##' Note that we assume that the formula uniquely specifies PSUs.
##' This will always be true if the PSUs were selected without replacement.
##' If they were selected with replacement, then it will be necessary
##' to make each realization of a given PSU in the sample a unique id.
##' Bottom line: the code below assumes that all observations within
##' each PSU (as identified by the design formula) are from the same draw
##' of the PSU.
##'
##' The rescaled bootstrap technique works by adjusting the
##' estimation weights based on the number of times each
##' row is included in the resamples. If a row is never selected,
##' it is still included in the returned results, but its weight
##' will be set to 0. It is therefore important to use estimators
##' that make use of the estimation weights on the resampled
##' datasets.
##'
##' We always take m_i = n_i - 1, according to the advice presented
##' in Rao and Wu (1988) and Rust and Rao (1996).
##'
##' @param survey.data the dataset to use
##' @param survey.design a formula describing the design of the survey (see below - TODO)
##' @param num.reps the number of bootstrap replication samples to draw
##' @param parallel if TRUE, use parallelization (via \code{plyr})
##' @param paropts an optional list of arguments passed on to \code{plyr} to control
##'        details of parallelization
##' @return a list with \code{num.reps} entries. each entry is a dataset which
##' has at least the variables \code{index} (the row index of the original
##' dataset that was resampled) and \code{weight.scale}
##' (the factor by which to multiply the sampling weights
##' in the original dataset).
##' @details \code{survey.design} is a formula of the form\cr
##'    weight ~ psu_vars + strata(strata_vars),
##' where weight is the variable with the survey weights and psu
##' is the variable denoting the primary sampling unit
##' @export
rescaled.bootstrap.sample.pureR <- function(survey.data,
                                      survey.design,
                                      parallel=FALSE,
                                      paropts=NULL,
                                      num.reps=1)
{

  survey.data$.internal_id <- 1:nrow(survey.data)

  design <- parse_design(survey.design)

  ## drop the "~" at the start of the formula
  psu.vars <- design$psu.formula[c(-1)][[1]]

  ## in the special case where there are no PSU vars, treat each row as
  ## its own PSU
  if (length(psu.vars)==1 & psu.vars=="1") {
    psu.vars <- as.name(".internal_id")
  }

  ## if no strata are specified, enclose the entire survey all in
  ## one stratum
  if (is.null(design$strata.formula)) {
    strata <- list(survey.data)
  } else {
    strata <- dlply(survey.data, design$strata.formula, identity)
  }

  ## get num.reps bootstrap resamples within each stratum,
  ## according to the rescaled bootstrap scheme
  ## (see, eg, Rust and Rao 1996)

  ## this llply call returns a list, with one entry for each stratum
  ## each stratum's entry contains a list with the bootstrap resamples
  ## (see the note for the inner llply call below)
  bs <- llply(strata,
              function(stratum.data) {

                ## figure out how many PSUs we have in our sample
                psu.count <- count(stratum.data,
                                   psu.vars)

                n.h <- nrow(psu.count)

                ## take m_h = n_h - 1, which is what the literature
                ## most commonly recommends
                m.h <- n.h - 1

                ## this llply call returns a list, with one entry for each bootstrap rep.
                ## each list entry has a data frame with the same number of rows as
                ## stratum.data, 
                ## and with colums for 
                ## the survey design variables, the .internal_id,
                ## r.hi, and weight.scale
                resamples <- llply(1:num.reps,
                                   ## for each bootstrap rep, this function returns
                                   function(rep) {

                                     ## sample m.h PSUs with replacement
                                     these.psu.samples <- sample(1:nrow(psu.count),
                                                                 m.h,
                                                                 replace=TRUE)

                                     ## r.hi is the number of times PSU i in stratum
                                     ## h was chosen in our resample
                                     r.hi <- count(data.frame(psu.row=these.psu.samples))

                                     r.hi$weight.scale <- r.hi$freq * (n.h / m.h)

                                     psu.count$freq <- NULL

                                     psu.count$r.hi <- 0
                                     psu.count$weight.scale <- 0

                                     psu.count[ r.hi$psu.row, "r.hi" ] <- r.hi$freq

                                     ## this is the factor by which we
                                     ## need to multiply
                                     ## the sampling weights
                                     ## (again, see, eg, Rust and Rao 1996, pg 292)
                                     psu.count[ r.hi$psu.row,
                                                "weight.scale" ] <- r.hi$weight.scale

                                     this.resample <- merge(stratum.data[,c(all.vars(psu.vars),
                                                                            ".internal_id")],
                                                            psu.count,
                                                            by=all.vars(psu.vars))
                                     this.resample$.internal_id.1 <- NULL

                                     return(this.resample)
                                   },
                                   .parallel=parallel,
                                   .paropts=paropts)

                return(resamples)
              })

  # now reassemble each stratum...
  res <- llply(1:num.reps,
               function(rep.idx) {
                 this.rep <- ldply(bs,
                                   function(this.stratum.samples) {
                                     return(this.stratum.samples[[rep.idx]])
                                   })

                 this.rep <- rename(this.rep,
                                    c(".internal_id"="index"))

                 return(this.rep)
               })

  return(res)
}

#####################################################
##' srs.bootstrap.sample
##'
##' given a survey dataset and a description of the survey
##' design (ie, which combination of vars determines primary sampling
##' units, and which combination of vars determines strata), take
##' a bunch of bootstrap samples under a simple random sampling
##' (with repetition) scheme
##'
##' @param survey.data the dataset to use
##' @param num.reps the number of bootstrap replication samples to draw
##' @param parallel if TRUE, use parallelization (via \code{plyr})
##' @param paropts an optional list of arguments passed on to \code{plyr} to control
##'        details of parallelization
##' @param ... ignored, but useful because it allows params like
##\code{survey.design},
##' which are used in other bootstrap designs, to be passed in without error
##' @return a list with \code{num.reps} entries. each entry is a dataset which has
##' at least the variables \code{index} (the row index of the original dataset that
##' was resampled) and \code{weight.scale} (the factor by which to multiply the
##' sampling weights in the original dataset).
##'
##' @details \code{survey.design} is not needed; it's included as an argument
##' to make it easier to drop \code{srs.bootstrap.sample} into the place of
##' other bootstrap functions, which do require information about the survey design
##' @export
srs.bootstrap.sample <- function(survey.data,
                                 num.reps=1,
                                 parallel=FALSE,
                                 paropts=NULL,
                                 ...)
{

  survey.data$.internal_id <- 1:nrow(survey.data)

  res <- llply(1:num.reps,
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

#####################################################
##' draw RDS bootstrap resamples for one chain
##'
##' this function uses the algorithm described in the
##' supporting online material for Weir et al 2012
##' (TODO PROPER CITE) to take bootstrap resamples
##' of one chain from an RDS dataset
##'
##' @param chain the chain to draw resamples for
##' @param mm the mixing model to use
##' @param dd the degree distns to use
##' @param parent.trait a vector whose length is the number
##' of bootstrap reps we want
##' @param idvar the name of the variable used to label the
##' columns of the output (presumably some id identifying the
##' row in the original dataset they come from -- see below)
##' @return a list of dataframes with one entry for each respondent in the chain.
##' each dataframe has one row for each bootstrap replicate. so if we take 10
##' bootstrap resamples of a chain of length 50, there will be 50 entries in
##' the list that is returned. each entry will be a dataframe with 10 rows.
rds.boot.draw.chain <- function(chain, mm, dd, parent.trait, idvar="uid") {

    thisid <- paste0(idvar, ".", chain$data[,idvar])

    ## choose the traits of the referrals from this set of parents
    trait.draws <- mm$choose.next.state.fn(parent.trait)

    ## choose the degree (and other vars) of each referral based on
    ## the drawn set of traits
    deg.draws <- dd$draw.degrees.fn(trait.draws)

    ## if no children, this is the last wave
    ## return a list whose only entry is dataframe of draws for this wave
    if (is.null(chain$children)) {
        return(setNames(list(deg.draws), thisid))
    }

    ## otherwise, if there are children, recursively get bootstrap
    ## draws for them
    child.draws <- llply(chain$children,
                         rds.boot.draw.chain,
                         mm=mm,
                         dd=dd,
                         parent.trait=trait.draws)
    child.draws <- unlist(child.draws, recursive=FALSE)

    
    return(setNames(c(list(deg.draws), child.draws), 
                    c(thisid, names(child.draws))))

}

#####################################################
##' draw RDS bootstrap resamples
##'
##' draw boostrap resamples for an RDS dataset, using
##' the algorithm described in the supporting online
##' material of Weir et al 2012 (TODO PROPER CITE)
##'
##' TODO -- consider constructing chains, mm from other args
##'
##' TODO be sure to comment the broken-out trait variables
##'      (ie these could all be different from the originals)
##'
##' @param chains a list whose entries are the chains
##' we want to resample
##' @param mm the mixing model
##' @param dd the degree distributions
##' @param num.reps the number of bootstrap resamples we want
##' @param keep.vars if not NULL, then the names of variables
##' from the original dataset we want appended to each bootstrap
##' resampled dataset (default is NULL)
##' @return a list of length \code{num.reps}; each entry in
##' the list has one bootstrap-resampled dataset
rds.chain.boot.draws <- function(chains,
                                 mm,
                                 dd,
                                 num.reps,
                                 keep.vars=NULL) {

    traits <- mm$traits

    ## get the bootstrap resamples for each chain
    res <- llply(chains,
                 function(this.chain) {

                     ## TODO -- should handle the case where there is
                     ## missingness in a seed trait?

                     seed.trait <- traits.to.string(this.chain$data[,traits,drop=FALSE],
                                                    traits)$traits

                     return(rds.boot.draw.chain(this.chain,
                                                mm,
                                                dd,
                                                rep(seed.trait,num.reps)))
                 })

    ## within each chain, we now have a list whose entries are dataframes,
    ## one for each repsondent, with one row for each bootstrap resample
    ## convert this into a list whose entries are dataframes, one for each
    ## bootstrap resample, and whose rows are respondents
    res.byboot <- llply(res,
                        function(this.chain.res) {
                            by.rep <- llply(1:num.reps,
                                            function(this.rep.id) {
                                                ldply(this.chain.res,
                                                      function(x) x[this.rep.id,])
                                            })
                            return(by.rep)
                        })

    ## assemble the bootstrap resamples from each chain together
    ## to end up with a list whose entries are datasets, one for each bootstrap
    ## resample (across all chains)
    res.dat <- llply(1:num.reps,
                     function(this.rep.id) {
                        ldply(res.byboot,
                              function(x) { 
                                  this.dat <- x[[this.rep.id]] 
                                  ## remove the .id column
                                  this.dat$.id <- NULL
                                  ## add the trait columns back in...
                                  return(cbind(this.dat,
                                               unparse.trait(this.dat$trait,
                                                             traits)))
                              })
                     })

    return(res.dat)

}

#####################################################
##' draw RDS bootstrap resamples using the
##' algorithm in Salganik 2006 (TODO PROPER CITE)
##'
##' this algorithm picks a respondent from the survey
##' to be a seed uniformly at random. it then generates
##' a bootstrap draw by simulating the markov process
##' forward for n steps, where n is the size of the draw
##' required.
##'
##' if you wish the bootstrap dataset to end up with
##' variables from the original dataset other than the
##' traits and degree, then you must specify this when
##' you construct \code{dd} using the 
##' '\code{\link{estimate.degree.distns}} function.
##'
##'
##' TODO be sure to comment the broken-out trait variables
##'      (ie these could all be different from the originals)
##'
##' @param chains a list with the chains constructed from the survey
##' using \code{\link{make.chain}}
##' @param mm the mixing model
##' @param dd the degree distributions
##' @param num.reps the number of bootstrap resamples we want
##' @return a list of length \code{num.reps}; each entry in
##' the list has one bootstrap-resampled dataset
rds.mc.boot.draws <- function(chains,
                              mm,
                              dd,
                              num.reps) {

    traits <- mm$traits

    all.data <- ldply(chains,
                      chain.data)

    ## NB: only using traits that are nonmissing
    all.traits <- traits.to.string(all.data[,traits,drop=FALSE], traits)
    all.traits.str <- all.traits$traits

    chain.sizes <- laply(chains,
                         chain.size)

    num.chains <- length(chain.sizes)

    res <- llply(1:num.reps,
                 function(this.rep) {

                     ## NB: an alternative would be do to the bootstrap this way,
                     ## but this is not what has been done in the literature so far
                     ## so it's commented out
                     ##
                     ## draw chains with the same length as the
                     ## empirically observed ones
                     ## draw seeds (with replacement, right?)
                     ## these.seeds <- sample(all.traits.str, size=num.chains, replace=TRUE)
                     ##
                     ## these.traits <- unlist(llply(1:num.chains,
                     ##                              function(x) {
                     ##                                  res <- mc.sim(mm,
                     ##                                                these.seeds[x],
                     ##                                                chain.sizes[x])
                     ##                                  return(res)
                     ##                              }))

                     ## use only one chain
                     this.seed <- sample(all.traits.str, size=1)

                     ## draw the chain of respondent traits, starting from our seed,
                     ## using the markov model we were given
                     these.traits <- mc.sim(mm,
                                            this.seed,
                                            sum(chain.sizes))

                     ## draw degrees for each of the respondents in the chain we
                     ## just drew
                     these.degs <- dd$draw.degrees.fn(these.traits)

                     return(cbind(these.degs,
                                  unparse.trait(these.traits, traits)))
                 })

    return(res)

}
