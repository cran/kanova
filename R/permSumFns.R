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
        eps <- .Machine$double.eps
        if(type=="oneway") {
# We *do* want the variance of the residuals from the saturated
# model, but here this is the same as the full model, so we needn't
# set satMod=TRUE in the call to estSigsq().

            sdpb <- sqrt(estSigsq(sumFns)) # pb = possibly bad
            nok  <- sdpb <= .Machine$double.eps
            sdok <- ifelse(nok,1,sdpb)
            ip   <- sample(1:N)
            prez <- lapply(1:N,function(k,x,sd){x[,k]/sd},
                                        x=rez,sd=sdok)[ip]
            psF  <- lapply(1:N,function(k,fv,sd,w,prez){
                              fvstar <- fv[[k]] + sd*prez[[k]]
                              attr(fvstar,"weight") <- w[k]
                              fvstar},fv=fV,sd=sdok,w=wts,prez=prez)
        } else if(type=="twoway") {
# Here we need to set satMod=TRUE, so as to get variance of the residuals
# from the saturated (interaction) model.
            sdpb <- lapply(estSigsq(sumFns,satMod=TRUE),sqrt) # pb = possibly bad
            sdok <- lapply(sdpb,function(x){
                              nok <- sapply(x,function(xx){xx <= eps})
                              ifelse(nok,1,x)
                              })
            ip   <- sample(1:N)
            prez <- lapply(1:N,function(k,x,sd){
                               sdloc <- sd[[colnames(x)[k]]]
                               x[,k]/sdloc
                               },x=rez,sd=sdok)[ip]
            psF  <- lapply(1:N,function(k,fv,sd,w,prez){
                               sdloc <- sd[[names(fv)[k]]]
                               fvstar <- fv[[k]] + sdloc*prez[[k]]
                               attr(fvstar,"weight") <- w[k]
                               fvstar},fv=fV,sd=sdok,w=wts,prez=prez)
       } else if(type=="interac") {
# As in the oneway setting, we *do* want the variance of the
# residuals from the saturated model, but again this is the same
# as the full model, so we needn't set satMod=TRUE in the call
# to estSigsq().
            sdpb <- lapply(estSigsq(sumFns),sqrt) # pb = possibly bad
            sdok <- lapply(sdpb,function(x){nok <- x <= eps
                             ifelse(nok,1,x)})
            ip   <- sample(1:N)
            prez <- lapply(1:N,function(k,x,sd){
                               sdloc <- sd[[colnames(x)[k]]]
                               x[,k]/sdloc},x=rez,sd=sdok)
            prez <- prez[ip]
            psF  <- lapply(1:N,function(k,fv,sd,w,prez){
                               sdloc <- sd[[names(fv)[k]]]
                               fvstar <- fv[[k]] + sdloc*prez[[k]]
                               attr(fvstar,"weight") <- w[k]
                               fvstar
                               },fv=fV,sd=sdok,w=wts,prez=prez)
            junk  <- sapply(1:N,function(k,fv,sd,w,prez){
                               sdloc <- sd[[names(fv)[k]]]
                               mean(fv[[k]] + sdloc*prez[[k]])
                               },fv=fV,sd=sdok,w=wts,prez=prez)
            crap <- sapply(fV,mean)
       }
    }
    attributes(psF) <- attributes(sumFns)
    psF
}
