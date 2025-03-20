iEngine <- function(sumFns) {
#
# Calculate the test statistic in the "interaction" case.`
#
A      <- attr(sumFns,"A")
B      <- attr(sumFns,"B")
AB     <- attr(sumFns,"AB")
Khat   <- wtdMean(sumFns)
Khati  <- lapply(split(sumFns,f=A),wtdMean)
Khatj  <- lapply(split(sumFns,f=B),wtdMean)
Khatij <- lapply(split(sumFns,f=AB),wtdMean)
s2ij   <- estSigsq(sumFns)
nr     <- length(attr(sumFns,"r"))
a      <- length(levels(A))
b      <- length(levels(B))
s2ij   <- array(unlist(s2ij),dim=c(a,b,nr))
M1     <- buildM1(Khat,Khati,Khatj,Khatij)
wts    <- getWts(sumFns)
wtm    <- matrix(sapply(split(wts,f=AB),sum),nrow=a,ncol=b)
V      <- buildV(s2ij,wtm)
MoV    <- M1^2/V
r      <- attr(sumFns,"r")
ens    <- table(A,B)
rslt   <- sum(ens*apply(MoV,c(1,2),trapint,r=r))
rslt
}
