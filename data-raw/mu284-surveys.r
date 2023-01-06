#####################################################################
##
## mu284-surveys.r
##
## take Sarndal et al's MU284 dataset and generate a handful
## of surveys from a known sampling design. then compute a
## bunch of RBS bootstrap samples from each one, based on
## the assumption that the current version of the RBS code
## (20130719) is correct. this conclusion comes from the
## results of our RJ simulation study.
## these results will form the basis of unit tests in
## the networkreporting package
##

library(networkreporting)
library(stringr)
library(plyr)
library(reshape2)
library(ggplot2)
library(functional)
library(sampling)

root.dir <- "~/Dropbox/dennis_matt"

code.dir <- paste(root.dir, "/code/gnsum-variance", sep="")
data.dir <- paste(root.dir, "/data/processed", sep="")
out.dir <- paste(root.dir, "/out/gnsum-variance", sep="")

## set the seed so that we can reproduce these results
set.seed(12345)

data(MU284)

## step 1: take a sample of K surveys using the same design
##         save these for use in the package
## design: described in Ex 4.3.2 (pg 142-3) of SSW
##   2-stage sample:
##     PSUs are the standard clusters for MU284;
##          we take an SI sample of n.I = 5 out of N.I = 50 of these
##     w/in each PSU, we take an SI sample of n.i = 3 out of N.i
##          municipalities

K <- 10
n.stage1 <- 5
n.stage2 <- 3
psus <- unique(MU284$CL)
ssus <- dlply(MU284, .(CL), function(x) { x$LABEL })

draw.survey <- function(all.psus=psus, all.ssus=ssus,
                        n1=n.stage1, n2=n.stage2) {
    these.psus <- sample(all.psus, size=n1, replace=FALSE)
    these.mun <- ldply(these.psus,
                       function(this.psu) {
                           all.mun <- all.ssus[[this.psu]]
                           selected.mun <- sample(all.mun, size=n2, replace=FALSE)
                           selected.data <- MU284[MU284$LABEL %in% selected.mun,]

                           ## add sampling weights
                           selected.data$sample_weight <- 1 / ((n1/length(all.psus))*
                                                               (n2/length(all.mun)))
                           return(selected.data)
                       })
    return(these.mun)
}

surveys <- llply(1:K,
                 function(x) { draw.survey() })

## step 2: for each survey, take a lot of RBS bootstrap samples
##         for estimators of
##         - a total (TODO)
##         - a ratio (TODO)

## TODO -- write / sketch a vignette describing what this fn needs to
## take as args
svy.estimators <- function(survey.data, weights) {
    survey.data$weight <- weights

    res <- summarise(survey.data,
                     TS82.hat=sum(S82*weight),
                     R.RMT85.P85.hat=sum(RMT85*weight)/sum(P85*weight))
    return(res)
}

num.bootstrap.samples <- 5000

rbsfn <- Curry(bootstrap.estimates,
               survey.design= ~ CL,
               num.reps=num.bootstrap.samples,
               estimator.fn="svy.estimators",
               weights="sample_weight",
               bootstrap.fn="rescaled.bootstrap.sample")

boot.res <- llply(surveys,
                  function(x) {
                      br <- rbsfn(survey.data=x)
                      return(do.call("rbind", br))
                  })

## step 3: compute summary values from the bootstrap resamples
##         for each survey, and save them to be used as unit
##         tests in the networkreporting package
boot.res.summ <- llply(boot.res,
                       function(svy.boot.res) {
                           this.summ <- summarise(svy.boot.res,
                                                  mean.TS82.hat=mean(TS82.hat),
                                                  mean.R.RMT85.P85.hat=mean(R.RMT85.P85.hat),
                                                  sd.TS82.hat=sd(TS82.hat),
                                                  sd.R.RMT85.P85.hat=sd(R.RMT85.P85.hat))
                           return(this.summ)
                       })

MU284.surveys <- surveys
MU284.boot.res.summ <- do.call("rbind", boot.res.summ)
MU284.estimator.fn <- svy.estimators

save(MU284.boot.res.summ,
     MU284.surveys,
     MU284,
     MU284.estimator.fn,
     file=file.path(out.dir, "mu284.RData"))

save.image(file.path(out.dir, "mu284-bootstrapped-surveys.RData"))
##load(file.path(out.dir, "mu284-bootstrapped-surveys.RData"))

## TODO NEXT
##  - write actual unit tests
##  - sketch vignette of how to add bootstrap fn...
