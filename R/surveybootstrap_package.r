##' Survey bootstrap variance estimators
##'
##' `surveybootstrap` has methods for analyzing data that were collected
##' using network reporting techniques. It includes estimators appropriate for
##' the simple bootstrap and the rescaled bootstrap.
##'
##' @docType package
##' @name surveybootstrap
##' @aliases surveybootstrap package-surveybootstrap
##' @import dplyr
##' @importFrom functional Curry
##' @importFrom stats rmultinom setNames terms update update.formula xtabs
NULL

##' @useDynLib surveybootstrap
##' @importFrom Rcpp sourceCpp
NULL

##' MU284 population
##'
##' Data used in unit tests for variance estimation.
##' See TODO-Sarndal TODO-sampling package
##' TODO-doc describing unit tests
##'
##' @name MU284
##' @docType data
NULL

