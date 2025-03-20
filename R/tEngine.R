tEngine <- function(sumFns) {
#
# Calculate the test statistic in the "twoway" case;
# testing for an A effect, allowing for a B effect.
#
# Form the relevant means.
#
A     <- attr(sumFns,"A")
a     <- length(levels(A))
xi    <- split(sumFns,f=A)
barxi <- lapply(xi,wtdMean)
barx  <- wtdMean(sumFns)
# Form the test statistic numerators in the form of a matrix M which
# is of dimension ((length(r)) x a), where r is attr(sumFns,"r").
M1 <- do.call(cbind,barxi) - barx
# The numerator of the test statistic is M1^2 (square the entries of M1)

# Form the variance estimates.  Each column of s2j is a vector 
# whose l-th entry is an (unbiased) estimate of sigma^2_j(r[l]) =
# var(K_{ijk}(r[l]), where r[l] is the l-th entry of the vector r
# that constitutes the argument of the diagnostic/summary function
# which is being used.
s2j <- do.call(cbind,estSigsq(sumFns))

# Form stilsq, the matrix whose (i,l)-th entry is the weighted
# average of {s^2_1(r[l]), s^2_2(r[l]), ... s^2_b(r[l])} with
# weights c(w[l,1],...,w[l,b]) normalised to sum to 1, i.e.
# (w_{l1}, ..., w_{lb})/w_{l.} .
wtm    <- matrix(sapply(split(getWts(sumFns),f=attr(sumFns,"AB")),sum),nrow=a)
stilsq <- s2j%*%t(wtm/apply(wtm,1,sum))
nu <- apply(wtm,1,sum)/sum(wtm)
zeta <- matrix(nu,nrow=length(nu),ncol=length(nu))
diag(zeta) <- (nu-1)^2/nu
V <- stilsq%*%zeta/sum(wtm) # t(zeta) ??? (seems OK; zeta seems to be symmetric).
# V is a matrix of dimension (length(r) x a).
MoV  <- M1^2/V
r    <- attr(sumFns,"r")
ens  <- table(A)
rslt <- sum(ens*apply(MoV,2,trapint,r=r))
rslt
}
