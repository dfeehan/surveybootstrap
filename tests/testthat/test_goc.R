## TODO -- test that there is some variance across a set of bootstrap
##         samples

## TODO -- for ratios, test also that there is some variance in
##         numerator and denominator

## TODO -- test paired ego / alter datasets

## TODO -- test that calling bootstrap.estimates works when
##         total.popn.size is an argument and not an attribute of
##         the data frame (had to use parent.frame(2)) to fix
##         a bug about this

## TODO -- test cases where estimates should never be negative

## TODO -- look at
## http://stackoverflow.com/questions/8898469/is-it-possible-to-use-r-package-data-in-testthat-tests-or-run-examples
## to try and figure out the real way to include package data in
## unit tests...

load("goc.RData")

set.seed(12345)

#########################################
## RDS - markov chain bootstrap (matt's algorithm)
context("goc / rds - build markov model")

mm <- estimate.mixing(survey.data, parent.data, c("use.crack"))

## idea: if we take a bunch of draws starting from each
## state, we should end up in neighboring states in proportion
## to the transition probabilities
this.state <- mm$states[["0"]]
tmp <- plyr::laply(1:10000, function(x) { this.state$trans.fn() })
tmptab <- as.numeric(table(tmp)/sum(table(tmp)))
expect_that(tmptab, equals(as.vector(this.state$trans.probs),
                           tolerance=.01, scale=1),
             label="trans.fn from state '0'")

this.state <- mm$states[["1"]]
tmp <- plyr::laply(1:10000, function(x) { this.state$trans.fn() })
tmptab <- as.numeric(table(tmp)/sum(table(tmp)))
expect_that(tmptab, equals(as.vector(this.state$trans.probs),
                           tolerance=.01, scale=1),
             label="trans.fn from state '1'")


## similar test to above, but using choose.next.state.fn
this.state <- mm$states[["0"]]
parents <- rep("0", 10000)
tmp <- mm$choose.next.state.fn(parents)
tmptab <- as.numeric(table(tmp)/sum(table(tmp)))
expect_that(tmptab, equals(as.vector(this.state$trans.probs),
                           tolerance=.01, scale=1),
             label="choose.next.state.fn from state '0'")

this.state <- mm$states[["1"]]
parents <- rep("1", 10000)
tmp <- mm$choose.next.state.fn(parents)
tmptab <- as.numeric(table(tmp)/sum(table(tmp)))
expect_that(tmptab, equals(as.vector(this.state$trans.probs),
                           tolerance=.01, scale=1),
             label="choose.next.state.fn from state '1'")

## mixture of both states...
parents <- c(rep("1", 5000), rep("0", 5000))
tmptab <- as.numeric(table(tmp)/sum(table(tmp)))
tp <- mean(c(mm$states[["0"]]$trans.probs[1], mm$states[["1"]]$trans.probs[1]))
tp <- c(tp, mean(c(mm$states[["0"]]$trans.probs[2], mm$states[["1"]]$trans.probs[2])))
expect_that(tmptab, equals(as.vector(this.state$trans.probs),
                           tolerance=.01, scale=1),
             label="choose.next.state.fn from mixture of states")

#########################################
## RDS - markov chain bootstrap (matt's algorithm)
context("goc / rds - degree estimation")

these.traits <- c("use.crack", "female")

dd <- estimate.degree.distns(survey.data,
                             d.hat.vals="netsize.5.bss",
                             traits=these.traits)

tt <- traits.to.string(survey.data, these.traits)
survey.data$tt[tt$used.idx] <- tt$traits

dmeans <- plyr::ddply(survey.data,
                plyr::.(tt),
                summarise,
                mean.degree=mean(netsize.5.bss))

for(cur.trait in c("0.0", "0.1", "1.0", "1.1")) {
    res <- dd$draw.degrees.fn(rep(cur.trait, 10000))
    expect_that(mean(res[,'degree']),
                equals(dmeans[paste(dmeans$tt)==cur.trait, "mean.degree"],
                       tol=.1),
                label=paste0("draw from degree distn for trait ", cur.trait))
}


#########################################
## RDS - static chain bootstrap (Weir et al algorithm)
context("variance estimators - rds static chain bootstrap - sanity checks")

## TODO
