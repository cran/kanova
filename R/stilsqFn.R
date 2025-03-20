stilsqFn <- function(s2j,wtm) {
# wtm is an a x b matrix of weights
m <- do.call(cbind,s2j) # m is (rsteps+1) x b.
x <- apply(wtm,2,sum)
wtm <- t(wtm)/apply(wtm,2,sum)
m%*%wtm
}
