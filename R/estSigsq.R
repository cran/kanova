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
        knt <- table(A,B)
        if(any(knt <= 1)) {
             iii <- which(knt <= 1,arr.ind=TRUE)
             badcells <- apply(iii,1,paste,collapse=",")
             badcells <- paste0("(",badcells)
             badcells <- paste0(badcells,")")
             badcells <- paste(badcells,collapse=", ")
             whinge <- paste0("Cells ",badcells," of the model have too few\n",
                              "  observations for variances to be calculated.\n")
             stop(whinge)
        }
        xxx <- split(sumFns,f=B)
        asv <- split(A,f=B)
        yyy <- lapply(1:length(xxx),function(k,x,f){split(x[[k]],f=f[[k]])},
                      x=xxx,f=asv) 
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
        lapply(zzz,function(z){z/(attr(z,"n")-1)})
    })
rslt
}
