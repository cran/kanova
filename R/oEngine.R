oEngine <- function(sumFns) {
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
wtv   <- sapply(split(getWts(sumFns),f=A),sum) # weight vector
wtsum <- sum(wtv)
vmlt  <- 1/wtv - 1/wtsum
V     <- outer(s2,vmlt,"*")

# Take the ratio.
MoV <- M/V

# Integrate and sum up.
r   <- attr(sumFns,"r")
ens <- table(A)
rslt <- sum(ens*apply(MoV,2,trapint,r=r))
rslt
}
