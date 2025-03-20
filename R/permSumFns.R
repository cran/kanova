permSumFns <- function(sumFns,rAndF,permtype) {
    N <- length(sumFns)
    type <- attr(sumFns,"type")
    if(permtype=="data") {
        if(type=="oneway") {
            ip <- sample(1:N,N)
            psF <- sumFns[ip]
        } else if(type=="twoway") {
            ip      <- permWithin(attr(sumFns,"B"))
            psF <- sumFns[ip]
        } else if(type=="interac") {
            whinge <- paste0("Cannot use \"permtype = \"data\" when",
                             " there is interaction in the model.\n")
            stop(whinge)
        }
    } else if(permtype == "stdres") {
        fV  <- rAndF$fitVals
        rez <- rAndF$resids
        wts <- getWts(sumFns)
        if(type=="oneway") {
# Here the full model and the saturated model are the same, and we
# *don't* want to use the variances from the two-way saturated model.
# (So leave satMod=FALSE in the call to estSigsq().)
            sd   <- sqrt(estSigsq(sumFns))
            ip   <- sample(1:N)
            prez <- lapply(1:N,function(k,x,sd,w){x[,k]/(sd*sqrt(w[k]))},
                                        x=rez,sd=sd,w=wts)[ip]
            psF  <- lapply(1:N,function(k,fv,sd,w,prez){
                              fvstar <- fv[[k]] + sd*sqrt(w[k])*prez[[k]]
                              attr(fvstar,"weight") <- w[k]
                              fvstar},fv=fV,sd=sd,w=wts,prez=prez)
        } else if(type=="twoway") {
            sd   <- lapply(estSigsq(sumFns,satMod=TRUE),sqrt)
            ip   <- sample(1:N)
            prez <- lapply(1:N,function(k,x,sd,w){
                               sdloc <- sd[[colnames(x)[k]]]
                               x[,k]/sdloc*sqrt(w[k])
                               },x=rez,sd=sd,w=wts)[ip]
            psF  <- lapply(1:N,function(k,fv,sd,w,prez){
                               sdloc <- sd[[names(fv)[k]]]
                               fvstar <- fv[[k]] + sdloc*sqrt(w[k])*prez[[k]]
                               attr(fvstar,"weight") <- w[k]
                               fvstar},fv=fV,sd=sd,w=wts,prez=prez)
       } else if(type=="interac") {
            sd   <- lapply(estSigsq(sumFns,satMod=TRUE),sqrt)
            ip   <- sample(1:N)
            prez <- lapply(1:N,function(k,x,sd,w){
                               sdloc <- sd[[colnames(x)[k]]]
                               x[,k]/sdloc*sqrt(w[k])},x=rez,sd=sd,w=wts)
            prez <- prez[ip]
            psF  <- lapply(1:N,function(k,fv,sd,w,prez){
                               sdloc <- sd[[names(fv)[k]]]
                               fvstar <- fv[[k]] + sdloc*sqrt(w[k])*prez[[k]]
                               attr(fvstar,"weight") <- w[k]
                               fvstar
                               },fv=fV,sd=sd,w=wts,prez=prez)
       }
    }
    attributes(psF) <- attributes(sumFns)
    psF
}
