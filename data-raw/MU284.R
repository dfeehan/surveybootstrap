# NOTE: the file mu284-surveys.r in this directory has the code
# that was originally used to generate these data

library(plyr)

load('data-raw/MU284.RData')


## set the seed so that we can reproduce these results
set.seed(12345)

# number of surveys we will take
K <- 10

# params for complex design
n.stage1 <- 5
n.stage2 <- 3
psus <- unique(MU284$CL)
ssus <- dlply(MU284, .(CL), function(x) { x$LABEL })

# params for simple design
n.simple <- n.stage1 * n.stage2

## SIMPLE DESIGN:
##
## draw K surveys from MU284 using a simple random sampling
## without replacement design
draw.simple.survey <- function(n=n.simple) {
  idx <- sample(1:nrow(MU284),
                n,
                replace=FALSE)

  res <- MU284[idx,]
  res$sample_weight <- nrow(MU284) / n

  return(res)
}

## COMPLEX DESIGN:
##
## draw K surveys from MU284 using the same design
## described in Ex 4.3.2 (pg 142-3) of SSW
##   2-stage sample:
##     PSUs are the standard clusters for MU284;
##          we take an SI sample of n.I = 5 out of N.I = 50 of these
##     w/in each PSU, we take an SI sample of n.i = 3 out of N.i
##          municipalities
draw.complex.survey <- function(all.psus=psus, all.ssus=ssus,
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

MU284.surveys <- llply(1:K,
                       function(x) { draw.simple.survey() })

MU284.complex.surveys <- llply(1:K,
                               function(x) { draw.complex.survey() })


#### Calculate results from bootstrap for use in unit tests
num.bootstrap.samples <- 5000

rbsfn <- Curry(bootstrap.estimates,
               survey.design= ~ CL,
               num.reps=num.bootstrap.samples,
               estimator.fn="MU284.estimator.fn",
               weights="sample_weight",
               bootstrap.fn="rescaled.bootstrap.sample")

boot.res <- llply(MU284.complex.surveys,
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
MU284.boot.res.summ <- do.call("rbind", boot.res.summ)


usethis::use_data(MU284, overwrite = TRUE)
usethis::use_data(MU284.surveys, overwrite = TRUE)
usethis::use_data(MU284.complex.surveys, overwrite = TRUE)
usethis::use_data(MU284.boot.res.summ, overwrite = TRUE)
