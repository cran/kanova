resAndFitCmpnts <- function(sumFns) {
type <- attr(sumFns,"type")
switch(EXPR=type,
    oneway = {
        Khat  <- wtdMean(sumFns)
        rslt  <- list(Khat=Khat)
    },
    twoway = {
        B     <- attr(sumFns,"B")
        xxx   <- split(sumFns,f=B)
        Khatj <- lapply(xxx,wtdMean)
        Khat  <- wtdMean(sumFns)
        rslt  <- list(Khat=Khat,Khatj=Khatj)
   },
   interac = {
        A      <- attr(sumFns,"A")
        B      <- attr(sumFns,"B")
        AB     <- attr(sumFns,"AB")
        xxx    <- split(sumFns,f=AB)
        Khatij <- lapply(xxx,wtdMean)
        xxx    <- split(sumFns,f=A)
        Khati  <- lapply(xxx,wtdMean)
        xxx    <- split(sumFns,f=B)
        Khatj  <- lapply(xxx,wtdMean)
        Khat   <- wtdMean(sumFns)
        rslt   <- list(Khat=Khat,Khati=Khati,Khatj=Khatj,Khatij=Khatij)
   })
   rslt
}
