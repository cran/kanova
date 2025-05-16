oEngine <- function(sumFns,divByVar) {
#
# Calculate the test statistic in the "one way" case.`
#
# Form the relevant means.
#
A     <- attr(sumFns,"A")
a     <- length(levels(A))
xi    <- split(sumFns,f=A)
barxi <- lapply(xi,wtdMean)
barx  <- wtdMean(sumFns)

# Form the test statistic numerators.
M <- (do.call(cbind,barxi) - barx)^2

# Form s^2. 
s2 <- estSigsq(sumFns) # vector of length = length(attr(sumFns,"r"))

# Form the multiplier.
if(divByVar) {
    wtv   <- sapply(split(getWts(sumFns),f=A),sum) # weight vector
    wtsum <- sum(wtv)
    vmlt  <- 1/wtv - 1/wtsum
    V     <- outer(s2,vmlt,"*")
    MoV   <- M/V
} else {
    MoV <- M
}

# Integrate and sum up.
r    <- attr(sumFns,"r")
ens  <- table(A)
sum(ens*apply(MoV,2,trapint,r=r))
}
