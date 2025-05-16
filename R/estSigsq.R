estSigsq <- function(sumFns,satMod=FALSE) {
#
# We may wish to use the estimated variances for the for the
# purpose of "normalising" or "homongenising" the summands of the
# test statistic.  Doing so is analogous to Hahn's "studentisation"
# procedure.  We may also wish to use them for the purpose of
# standardising the residuals.  In the latter case they must be
# the variances appropriate to the saturated model.  In the single
# classification setting, and (of course) in the interaction
# setting, the full model and the saturated model are the same.
# Only in the two-classification ("twoway") setting is there a
# difference, and only then does the satMod argument have an impact.
# The impact is to change "type" to "interac" (which is the saturated
# model in the two-classification setting).

type <- attr(sumFns,"type")
if(type=="twoway" & satMod) type <- "interac"

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
        knt <- table(asv)
        if(any(knt <= 1)) {
             Anm <- attr(sumFns,"Anm")
             iii <- which(knt <= 1)
             badlev <- paste(levels(asv)[iii],collapse=" ")
             whinge <- paste0("Levels ",badlev," of factor ",Anm," have too few\n",
                              "  observations for variances to be calculated.\n")
             stop(whinge)
        }
        yyy <- lapply(xxx,function(x,f){split(x,f=f)},f=asv) 
        zzz <- lapply(yyy,function(y){lapply(y,wtdSS)})
        lapply(zzz,function(z){
                      ens  <- sapply(z,function(v){attr(v,"n")})
                      ndot <- sum(ens)
                      den  <- ndot - length(z)
                      sss  <- Reduce("+",z)/den
                   })
        },
    interac={
        AB  <- attr(sumFns,"AB")
        xxx <- split(sumFns,f=AB)
        zzz <- lapply(xxx,wtdSS)
        sss <- lapply(zzz,function(z){z/(attr(z,"n")-1)})
    })
rslt
}
