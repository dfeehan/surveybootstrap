% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rds.r
\name{unparse.trait}
\alias{unparse.trait}
\title{unparse a collapsed trait string}
\usage{
unparse.trait(trait.string, names, sep = "\\\\.")
}
\arguments{
\item{trait.string}{a vector whose values are collapsed
traits}

\item{names}{a vector with the names of each trait (in order)}

\item{sep}{the character used to separate the traits in their
collpased string representation}
}
\value{
a dataframe whose rows correspond to the entries in
\code{trait.string}, with one column per trait
}
\description{
for a few of the RDS-related functions, it is useful
to combine several traits into one variable as a string;
for example, "male" and "young" might become
"male.young". this function takes a string with
combined traits and explodes it back into
several variables
}
\keyword{internal}
