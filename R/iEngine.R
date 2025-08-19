iEngine <- function(sumFns,divByVar) {
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
M      <- buildM1(Khat,Khati,Khatj,Khatij)^2
if(divByVar) {
    s2ij   <- estSigsq(sumFns)
    nr     <- length(attr(sumFns,"r"))
    a      <- length(levels(A))
    b      <- length(levels(B))
    s2ij   <- aperm(array(unlist(s2ij),dim=c(nr,a,b)),c(2,3,1))
    wts    <- getWts(sumFns)
    wtm    <- matrix(sapply(split(wts,f=AB),sum),nrow=a,ncol=b)
    V      <- buildV(s2ij,wtm)
    MoV    <- M/V
} else {
    Mov <- M
}

r      <- attr(sumFns,"r")
ens    <- table(A,B)
sum(ens*apply(MoV,c(1,2),trapint,r=r))
}
