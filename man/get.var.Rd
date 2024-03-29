% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.r
\name{get.var}
\alias{get.var}
\title{Get a variable from a dataframe or vector}
\usage{
get.var(survey.data, var, default = NA)
}
\arguments{
\item{survey.data}{The survey dataset}

\item{var}{Either \code{NULL}, a column name, or a vector of values}

\item{default}{The default value to fill in if the variable
is not found}
}
\value{
A vector of values whose length is the same as the
number of rows in \code{survey.data}; if \code{var} is \code{NULL}, this has
the default values
}
\description{
This function was written because a few of the estimator functions
need to use weights, and there are several cases to handle:
the user could pass in a column name, a vector of weights, or
nothing (in which case, the weights should default to 1 for each
row in the dataset). For the special case of getting weights, look
at the curried function \code{\link[=get.weights]{get.weights()}}
}
\keyword{internal}
