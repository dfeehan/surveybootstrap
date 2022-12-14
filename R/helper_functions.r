#######################################################################
## helper_functions.R
##
## various helper functions for use in the networkreporting
## package
##
## these are all internal
##


##########################################################################
##' Grab a function based on its name
##'
##' Helper to grab a function that is passed in as an argument
##'
##' This is based on Hadley Wickham's response to an SO
##' post:
##' <https://stackoverflow.com/questions/14183766/match-fun-provide-error-with-functions-defined-inside-functions>
##' with some minor modifications
##'
##' @param fn The function to search for
##' @param env The environment to start searching in
##' @return `fn`, if `fn` is already a function; otherwise, the first function found
##'         in env or one of its parents whose name is `fn`
##' @keywords internal
get.fn <- function(fn, env = parent.frame()) {

    ## base case: fn is already a function
    if (is.function(fn)) {
      return(fn)
    }

    ## base case: nothing left to search through
    if (identical(env, emptyenv())) {
        stop("Could not find function ", fn, "!")
    }

    ## base case: found function in env
    if (exists(fn, env, inherits=FALSE) &&
        is.function(env[[fn]])) {
        return(env[[fn]])

    ## recursive case: look through the environment
    ## above env
    } else {
        return(get.fn(fn, parent.env(env)))
    }

}

##########################################################################
##' Get a variable from a dataframe or vector
##'
##' This function was written because a few of the estimator functions
##' need to use weights, and there are several cases to handle:
##' the user could pass in a column name, a vector of weights, or
##' nothing (in which case, the weights should default to 1 for each
##' row in the dataset). For the special case of getting weights, look
##' at the curried function [get.weights()]
##'
##' @param survey.data The survey dataset
##' @param var Either `NULL`, a column name, or a vector of values
##' @param default The default value to fill in if the variable
##'        is not found
##' @return A vector of values whose length is the same as the
##'         number of rows in `survey.data`; if `var` is `NULL`, this has
##'         the default values
##' @keywords internal
get.var <- function(survey.data, var, default=NA) {

  ## weights will default to 1 for everyone, unless the user specified
  ## a weights variable
  if (is.null(var)) {

    return(rep(default, nrow(survey.data)))

  } else if (length(var) == 1) {

    ## ... otherwise, see if the weights variable is referring
    ## to a column of the dataframe; try to
    ## grab sampling weights from survey dataframe
    var.vals <- try(survey.data[,var,drop=FALSE],
                    silent=TRUE)

    ## if( inherits(var.vals, "try-error") ||
    ##    ncol(var.vals) != 1 ||
    ##    ! is.numeric(var.vals[,1]) ) {
    if( inherits(var.vals, "try-error") ||
       ncol(var.vals) != 1) {

      stop(paste(var,
                 " does not identify a valid column in the data.\n"))
    }

    var <- var.vals[,1]

    if (inherits(var.vals, "tbl_df")) {
        var <- collect(select(var.vals, 1))[[1]]
    }

    return(var)

  } else if (length(var) == nrow(survey.data)) {

    ## if var a vector with one entry per row, then these
    ## are our values
    return(var)
  } else {
    stop("can't determine what the values should be for ", var, ".")
  }


}

##########################################################################
##' Get the weights column from a dataframe
##'
##' This is the same as [get.var()] with the default value set to 1
##' instead of `NA`
##' @param ... (this is a function curried from [get.var()])
##' @keywords internal
get.weights <- functional::Curry(get.var, default=1)

##########################################################################
##' Only prints things out in verbose mode
##'
##' @param verbose If `TRUE`, print things out; otherwise, do nothing
##' @param ... Arguments to pass to cat if verbose is `TRUE`
##' @keywords internal
vcat <- function(verbose=TRUE, ...) {

  if(verbose) {
    message(...)
  }

  invisible()
}

##########################################################################
##' Parse a formula that describes the design of a survey
##'
##' Parse a formula of the form
##' `~ psu_v1 + psu_v2 + ... + strata(strata_v1 + strata_v2 + ...)`
##' into a PSU formula and a strata formula.
##'
##' @param formula a formula describing the sample design (see Description of [bootstrap.estimates()])
##' @return a list with entries `psu.formula` and `strata.formula`
##' @keywords internal
parse_design <- function(formula) {

  ## see https://stackoverflow.com/questions/10224805/how-to-select-a-part-of-formula-in-formula-in-r
  ## for some helpful info

  ## TODO / wishlist:
  ##  - check to be sure no response is included (or warn)
  ##  - check formulas for strata more carefully

  psu.formula <- formula
  strata.formula <- NULL

  these.labels <- attr(terms(formula), "term.labels")

  strata.idx <- grep("strata\\(", these.labels)

  if (length(strata.idx) == 1) {

    # grab the expression in the strata(...) part of the formula
    strata.text <- stringr::str_match(these.labels[strata.idx],
                                      "strata\\((.+)\\)")[2]

    ## updating instead of creating a new formula b/c this preserves
    ## the environment that the original formula was created in...
    strata.formula <- update(formula,
                             paste("~ ", strata.text))

    psu.formula <- update.formula(formula,
                                  paste("~ . - strata(",strata.text,")"))

  } else if (length(strata.idx > 1)) {

    stop("Cannot have more than one strata() specification in the design formula.")
  }

  return(list(psu.formula=psu.formula,
              strata.formula=strata.formula))

}


##########################################################################
##' Compute the weighted mean
##'
##' Given a vector of values and a vector of weights, compute the
##' weighted mean
##'
##' @param x The vector of values
##' @param w The vector of weights
##' @param na.rm if `TRUE`, only consider elements of `x` that are not missing
##'              (and their corresponding entries in `w`). Defaults to `FALSE`.
##' @return The weighted mean
##' @keywords internal
weighted.mean <- function(x, w, na.rm=FALSE) {

  if (na.rm) {
    idx <- (1:length(x))[!is.na(x)]
  } else {
    idx <- 1:length(x)
  }

  return(sum(x[idx]*w[idx])/sum(w[idx]))
}

