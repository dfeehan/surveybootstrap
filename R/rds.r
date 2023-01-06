#####################################################
## rds.r
##
## limited implementation of RDS estimator

#####################################################
##' Take a set of traits and turn into a string
##'
##' This is a helper function that is useful when we wish
##' to make several traits into one variable
##'
##' @param data The respondent info
##' @param traits The names of the traits to build the model on
##' @param na.action Defaults to 'drop' (meaning all rows of data
##'   with any missingness on the traits are dropped). Anything else
##'   means `NA`s are treated like any other value.
##' @param sep The separator character used to combine values
##' @return A list whose entries are
##'   * `used.idx`, which indicates which rows from the original dataset were used
##'     (may not be all of them if there is missingness); and
##'   * `traits`, which has the string version of the traits
##' @keywords internal
traits.to.string <- function(data, traits, na.action="drop", sep=".") {

  if (na.action == "drop") {
    touse.idx <- plyr::aaply(data[, traits],
                       1,
                       function(x) { return(! any(is.na(x))) },
                       .expand=FALSE)
    touse.idx <- which(touse.idx)
  } else {
    touse.idx <- 1:nrow(data)
  }

  traits.str <- plyr::aaply(data[touse.idx, traits],
                      1,
                      paste0,
                      collapse=sep,
                      .expand=FALSE)

  return(list(used.idx=touse.idx,
              traits=traits.str,
              sep=sep,
              names=traits))

}

#####################################################
##' Unparse a collapsed trait string
##'
##' For a few of the RDS-related functions, it is useful
##' to combine several traits into one variable as a string;
##' for example, "male" and "young" might become
##' "male.young". this function takes a string with
##' combined traits and explodes it back into
##' several variables.
##'
##' @param trait.string A vector whose values are collapsed
##' traits
##' @param names A vector with the names of each trait (in order)
##' @param sep The character used to separate the traits in their
##' collapsed string representation
##' @return A dataframe whose rows correspond to the entries in
##' `trait.string`, with one column per trait
##' @keywords internal
unparse.trait <- function(trait.string, names, sep="\\.") {

    if (sep == ".") {
        sep <- "\\."
    }

    vals <- stringr::str_split(trait.string, sep)

    vals <- plyr::ldply(vals, as.numeric)

    colnames(vals) <- names

    return(vals)

}

#####################################################
##' Estimate degree distributions by trait
##'
##' Break down RDS degree distributions by trait,
##' and return an object which has the degrees
##' for each trait as well as functions to draw
##' degrees from each trait.
##'
##' @details One of the items returned as a result is a function,
##' `draw.degrees.fn`, which takes one argument,
##' `traits`. This is a vector of traits and,
##' for each entry in this vector, `draw.degress.fn`
##' returns a draw from the empirical distribution of
##' degrees among respondents with that trait. So,
##' `draw.degrees.fn(c("0.0", "0.1", "0.1")` would
##' return a degree drawn uniformly at random from among
##' the observed degrees of respondents with trait "0.0"
##' and then two degrees from respondents with trait "0.1"
##'
##' @param survey.data The respondent info
##' @param d.hat.vals The variable that contains
##'    the degrees for each respondent
##' @param traits A vector of the names of the columns
##'    of `survey.data` which refer to the traits
##' @param keep.vars Additional vars to return along with degrees
##' @return An object with
##'   * `distns` a list with one entry per trait value; each
##entry has a dataframe with all of the degrees from respondents with
##the given trait
##'   * `draw.degrees.fn` a function which gets called with one
##argument, \code{traits}. See description above.
##'   * `keep.vars` the name of the other vars that are kept (if any)
##'
estimate.degree.distns <- function(survey.data,
                                   d.hat.vals,
                                   traits,
                                   keep.vars=NULL) {

  st <- traits.to.string(survey.data,
                         traits)

  degs <- get.var(survey.data, d.hat.vals)

  ## NOTE: we need to guarantee that the order of deg.dat's columns is
  ##      1: trait
  ##      2: degree
  ## [3...]: keep.vars
  if (! is.null(keep.vars)) {
      ## TODO -- should eventually make grabbing these others vars more robust
      other.vars <- survey.data[st$used.idx, keep.vars]
      deg.dat <- data.frame(trait=st$traits,
                            degree=degs[st$used.idx],
                            other.vars)
  } else {
      deg.dat <- data.frame(trait=st$traits, degree=degs[st$used.idx])
  }

  ## TODO -- for now, we can just represent the degrees with duplicates
  ## and use SI sampling to pick one when we need to. if there are
  ## huge datasets, this might need to be changed later

  ## to placate R CMD CHECK
  trait <- NULL

  deg.distns <- plyr::dlply(deg.dat,
                      plyr::.(trait),
                      identity)

  deg.fns <- unlist(plyr::llply(deg.distns,
                          function(this.trait.deg) {
                              ## if we don't force evaluation here,
                              ## R's lazy evaluation implies that only the
                              ## last version of this.trait.deg will
                              ## get used
                              ## (see, eg,
                              ##  http://adv-r.had.co.nz/Functions.html)
                              force(this.trait.deg)
                              return(function(n=1) {
                                  idx <- sample(1:nrow(this.trait.deg), size=n, replace=TRUE)
                                  return(this.trait.deg[idx,])
                              })
                          }))

  draw.degrees.fn <- function(traits) {
      tocall <- deg.fns[traits]
      degs <- plyr::llply(tocall, do.call, args=list())
      degs <- do.call("rbind", degs)
      rownames(degs) <- NULL
      return(degs)
  }

  return(list(distns=deg.distns, draw.degrees.fn=draw.degrees.fn, keep.vars=keep.vars))

}

#####################################################
##' Construct a mixing model from GoC/RDS data
##'
##' Given a dataset with the respondents and a dataset
##' on the parents (in many cases the same individuals),
##' and a set of relevant traits,
##' estimate mixing parameters and return a markov model.
##'
##' @param survey.data The respondent info
##' @param parent.data The parent info
##' @param traits The names of the traits to build the model on
##' @return A list with entries:
##'   * `mixing.df` the data used to estimate the mixing function
##'   * `choose.next.state.fn` a function which can be passed
##'      a vector of states and will return a draw of a subsequent state
##'      for each entry in the vector
##'   * `mixing.df` a dataframe (long-form) representation of
##'      the transition counts used to estimate the transition probabilities
##'    * `states` a list with an entry for each state. within
##'       each state's entry are
##'        - `trans.probs` a vector of estimated transition probabilities
##'        - `trans.fn` a function which, when called, randomly chooses a
##'           next state with probabilities given by the transition probs.
##'
estimate.mixing <- function(survey.data, parent.data, traits) {

  ## reduce to dataframe with [ child trait, parent trait ]
  ## then basically do a cross tab

  pkey <- attr(parent.data, "key")
  ckey <- attr(survey.data, "key")

  if (is.null(pkey) || is.null(ckey)) {
    stop("parent and survey datasets need to have an attribute which indicates what their keys are")
  }

  st <- traits.to.string(survey.data,
                         traits)
  pt <- traits.to.string(parent.data,
                         traits)

  parent.tomix <- data.frame(key=parent.data[pt$used.idx,pkey],
                             parent.trait=pt$traits)
  survey.tomix <- data.frame(key=survey.data[st$used.idx,ckey],
                             child.trait=st$traits)

  mix.data <- merge(survey.tomix,
                    parent.tomix,
                    by="key",
                    all.x=TRUE)

  res <- list()

  ## we'll reutrn a dataframe which has the parameters we use to estimate the mixing pattern
  res$mixing.df <- as.data.frame(xtabs(~ child.trait + parent.trait, data=mix.data))

  # to placate R CMD CHECK
  parent.trait <- NULL

  ## return a fn which will give us the next step in the chain from each state
  ## based on these transition probabilities
  res$states <- plyr::dlply(res$mixing.df,
                    plyr::.(parent.trait),
                    function(this.trait) {

                      if (all(this.trait$Freq == 0)) {
                        return(NULL)
                      }

                      probs <- this.trait$Freq / sum(this.trait$Freq)
                      names(probs) <- this.trait$child.trait

                      return(list(trans.probs=probs,
                                  trans.fn=function() {
                                    draw <- stats::rmultinom(1, 1, probs)
                                    return(names(probs)[as.logical(draw)])
                                  }))
                    })

  ## given a list of preceding states, this function obtains a
  ## succeeding state for each one
  res$choose.next.state.fn <- function(prev.states) {
      tfns <- plyr::llply(res$states, function(x) { x$trans.fn })

      tocall <- tfns[prev.states]

      next.states <- unlist(plyr::llply(tocall, do.call, args=list()))
      return(next.states)

  }

  res$traits <- traits

  return(res)
}


#####################################################
##' Run a markov model
##'
##' Run a given markov model for n time steps, starting
##' at a specified state.
##'
##' This uses the markov model produced by [estimate.mixing()]
##'
##' @param mm The markov model object returned by [estimate.mixing()]
##' @param start The name of the state to start in
##' @param n The number of time-steps to run through
##' @return A vector with the state visited at each time step. the first entry
##'         has the starting state
mc.sim <- function(mm, start, n) {

  if (n <= 1) {
      return(start)
  }

  if (is.null(mm) || is.null(mm$state) || is.null(mm$state[[ start ]])) {
    stop(paste("there was a problem starting at state", start))
  }

  path <- rep(NA, n)

  path[1] <- start

  for(i in 2:n) {
    path[i] <- mm$state[[ path[i-1] ]]$trans.fn()
  }

  return(path)

}

###########################################################
##' Determine whether or not one id is a parent of another
##'
##' This function allows us to determine which ids are
##' directly descended from which other ones. It is the only part
##' of the code that relies on the ID format used by the
##' Curitiba study (see Details); by modifying this function,
##' it should be possible to adapt this code to another study.
##'
##' @details See:
##'   * Salganik, M. J., Fazito, D., Bertoni, N., Abdo, A. H., Mello, M. B., &
##'     Bastos, F. I. (2011). Assessing network scale-up estimates for groups
##'     most at risk of HIV/AIDS: evidence from a multiple-method study of heavy
##'     drug users in Curitiba, Brazil. *American journal of epidemiology*,
##'     174(10), 1190-1196.
##'
##' @param id the id of the potential child
##' @param seed.id the id of the potential parent
##' @return TRUE if `id` is the direct descendant of `seed.id`
##' and FALSE otherwise
is.child.ct <- function(id, seed.id) {

    res <- stringr::str_locate(paste(id), paste(seed.id))

    if (! any(is.na(res)) &&
        res[1] == 1 &&
        res[2] == (nchar(id)-1)) {
        return(TRUE)
    }

    return(FALSE)
}

#####################################################
##' Build an RDS seed's chain from the dataset
##'
##' Note that this assumes that the chain is a tree (no loops)
##'
##' @param seed.id The id of the seed whose chain we
##' wish to build from the dataset
##' @param survey.data The dataset
##' @param is.child.fn A function which takes two ids as arguments;
##' it is expected to return `TRUE` if the second argument is the parent of the
##' first, and `FALSE` otherwise. it defaults to [is.child.ct()]
##' @return info
make.chain <- function(seed.id, survey.data, is.child.fn=is.child.ct) {

    key.var <- attr(survey.data, "key")
    keys <- survey.data[,key.var]

    is.child <- get.fn(is.child.fn)

    these.data <- survey.data[keys==seed.id,]

    child.ids <- keys[which(plyr::laply(keys, is.child, seed.id=seed.id))]

    ## if no children of this seed, finish
    if (length(child.ids)==0) {
        return (list(data=these.data,
                     children=NULL))
    }

    return(list(data=these.data,
                children=plyr::llply(child.ids,
                               make.chain,
                               survey.data=survey.data)))

}

#####################################################
##' Get the height (maximum depth) of a chain
##'
##' Get the height (maximum depth) of a chain
##'
##' @param chain The chain object
##' @return The maximum depth of the chain
max.depth <- function(chain) {
    if (is.null(chain$children)) {
        return(1)
    }

    return(1 + max(plyr::laply(chain$children,
                         max.depth)))
}

#####################################################
##' Get the size of a chain
##'
##' Count the total number of respondents in the chain
##' and return it
##'
##' @param chain The chain object
##' @return The number of respondents involved in
##' the chain
chain.size <- function(chain) {

    if (is.null(chain$children)) {
        return(1)
    }

    return(1 + sum(unlist(lapply(chain$children,
                                 chain.size))))
}

#####################################################
##' Get all of the values of the given variable
##' found among members of a chain
##'
##' @param chain The chain to get values from
##' @param qoi.var The name of the variable to
##' get from each member of the chain
##' @return A vector with all of the values of `qoi.var`
##' found in this chain. (Currently, the order of the values
##' in the vector is not guaranteed.)
chain.vals <- function(chain, qoi.var="uid") {

    if (is.null(chain$children)) {
        return(chain$data[,qoi.var])
    }

    return(c(chain$data[,qoi.var],
             unlist(lapply(chain$children,
                           chain.vals,
                           qoi.var=qoi.var))))
}

#####################################################
##' Get a dataset from a chain
##'
##' Take the data for each member of the given chain
##' and assemble it together in a dataset.
##'
##' @param chain The chain to build a dataset from
##' @return A dataset comprised of all of the chain's
##' members' data put together. The order of the rows
##' in the dataset is not specified.
chain.data <- function(chain) {
    if (is.null(chain$children)) {
        return(chain$data)
    }

    return(rbind(chain$data,
                 plyr::ldply(chain$children,
                       chain.data)))
}
