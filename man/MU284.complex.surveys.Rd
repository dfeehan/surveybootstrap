% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.r
\name{MU284.complex.surveys}
\alias{MU284.complex.surveys}
\title{Simulated sample surveys drawn from the MU284 Population using a complex design}
\format{
A list with 10 data frames, each with 15 rows and 11 columns:
\describe{
\item{LABEL, ..., CL}{Same as MU284 dataset}
\item{sample_weight}{The sampling weight for the row}
}
}
\description{
A list with 10 sample surveys with sample size 15 drawn from the \link{MU284}
dataset using a complex sampling design.
}
\details{
The sampling design comes from Ex. 4.3.2 (pg 142-3) of
'Model Assisted Survey Sampling' by Sarndal, Swensson, and Wretman
(2003, ISBN:0387406204).

The design is a two-stage sample:
\itemize{
\item stage I: the primary sampling units (PSUs) are the standard clusters from
\link{MU284}; we take a simple random sample without replacement of \code{n_I = 5}
out of \code{N_I = 50} of these
\item stage II: within each sampled PSU, we take a simple random sample without
replacement of \code{n_i = 3} out of \code{N_i} municipalities
}
}
