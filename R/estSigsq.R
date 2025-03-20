estSigsq <- function(sumFns,satMod=FALSE) {
#
# If satMod is FALSE, the estimated variances ("sigma-squared")  are
# be used used for the purpose of "standardising" the test statistic,
# which is analogous to Hahn's "studentisation" procedure.  If
# satMod is TRUE they are used for standardising the residuals.
# In this case they must be the variances appropriate to the
# saturated model.  Note that in the single classification setting
# the full model and the saturated model are the same, so there is
# only one sort of variance.  Only in the two classification setting
# is there a difference, whereby we need to choose between the
# "twoway" model and the saturated model.
#
type    <- attr(sumFns,"type")
if(type=="oneway") {
    if(satMod) stop("Do not set satMod=TRUE when \"type\" is \"oneway\".\n")
} else if(type=="twoway") {
    if(satMod) type <- "interac"
}
rslt <- switch(EXPR=type,
    oneway={
        xxx  <- split(sumFns,f=attr(sumFns,"A"))
        SS   <- lapply(xxx,wtdSS)
        ndot <- length(sumFns)
        a    <- length(SS)
        Reduce("+",SS)/(ndot - a)
    },
    twoway={
        A   <- attr(sumFns,"A")
        B   <- attr(sumFns,"B")
        xxx <- split(sumFns,f=B)
        asv <- split(A,f=B)[[1]]
        yyy <- lapply(xxx,function(x,f){split(x,f=f)},f=asv) 
        zzz <- lapply(yyy,function(y){lapply(y,wtdSS)})
        lapply(zzz,function(z){
                      urk <- sapply(z,function(v){attr(v,"n")})
                      ndot <- sum(urk)
                      Reduce("+",z)/(ndot - length(z))
                   })
        },
    interac={
        AB  <- attr(sumFns,"AB")
        xxx <- split(sumFns,f=AB)
        zzz <- lapply(xxx,wtdSS)
        sss <- lapply(zzz,function(w){w/(attr(w,"n")-1)})
    })
rslt
}
