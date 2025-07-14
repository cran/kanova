kanova <- function(fmla,data,expo=0,rsteps=128,r=NULL,sumFnNm=NULL,
                   warnSFN=TRUE,test=TRUE,bylevel=FALSE,
                   permtype=c("stdres","data"),nperm=99,brief=TRUE,
                   verb=TRUE,keepdata=FALSE,divByVar=TRUE) {
#
# Function to conduct one or two-way analysis of variance of
# summary functions (Kest, Fest, Gest, or Jest)  of replicated point
# patterns classified by a grouping factor A or two grouping factors
# A and B.
#

# Dig out the response name.
if(length(fmla) == 2) {
    rspNm <- names(data)[1]
    cform <- paste(rspNm,"~",as.character(fmla[2]))
    fmla  <- stats::as.formula(cform)
} else {
    rspNm <- as.character(fmla[2])
    if(!rspNm %in% names(data)) {
        whinge <- paste0("The response name is ",rspNm,", which is not the name\n",
                         "  of any of the columns of the hyperframe \"data\".\n")
        stop(whinge)
    }
}

# Set the summary function name.
if(inherits(data[,rspNm,drop=TRUE],"ppplist")) {
    if(is.null(sumFnNm)) sumFnNm <- "Kest"
    if(!(sumFnNm %in% c("Kest","Fest","Gest","Jest")) & warnSFN) {
        whinge <- paste0("Argument \"sumFnNm\" is \"",sumFnNm,"\", which is not\n",
                         "  one of the standard four.  The results may be fragile.\n")
        warning(whinge)
    }
} else {
    sumFnNm <- NA
}

# Augment object "data" from the parent frame if necessary.
preds  <- attr(terms(fmla),"term.labels")
npreds <- length(preds)
if(npreds == 0) stop("No predictors. (???)\n")
for(i in 1:npreds) {
    if(grepl(":",preds[i])) next
    if(!preds[i] %in% names(data)) {
       xp <- try(get(preds[i],envir=parent.frame()),silent=TRUE)
       if(inherits(xp,"try-error"))
           stop(paste0("Cannot find predictor \"",preds[i],"\".\n"))
       newcol <- data.frame(xp)
       names(newcol) <- preds[i]
       data <- cbind(data,newcol)
    }
}
if(bylevel) {
    if(npreds!=2) {
       stop("The number of predictors must be 2 when \"bylevel\" is TRUE.\n")
    }
    cform <- paste(rspNm,"~",preds[1])
    fmla  <- stats::as.formula(cform)
    B     <- factor(data[[preds[2]]])
    sdata <- split(data,f=B)
    b     <- length(sdata)
    rslt  <- vector("list",b)
    for(i in 1:b) {
        rslt[[i]] <- kanova(fmla,data=sdata[[i]],expo=expo,rsteps=rsteps,r=r,
                            sumFnNm=sumFnNm,warnSFN=FALSE,test=test,bylevel=FALSE,
                            permtype=permtype,nperm=nperm,brief=brief,
                            verb=verb,keepdata=FALSE)
    }
    class(rslt) <- "multi.kanova"
    if(test) {
# Form the "overall p-value".
        pvmin <- min(sapply(rslt,function(x){x[["pvalue"]]}))
        oapv  <- 1 - (1 - pvmin)^b
        attr(rslt,"oapv") <- oapv
    }
    if(keepdata) attr(rslt,"data") <- data
    return(rslt)
} else {
    permtype <- match.arg(permtype)
    if(npreds > 3) {
        whinge <- paste0("The length of the vector of predictor names, (",
                         paste(preds,collapse=", "), "), is ",npreds,".\n",
                         "  It must be at most 3.\n")
        stop(whinge)
    }
    switch(EXPR=npreds,
        {Anm  <- preds[1]
         Bnm  <- NULL
         type <- "oneway"
         effNm <- Anm
        },
        {Anm <- preds[1]
         Bnm <- preds[2]
         type <- "twoway"
         effNm <- paste0(Anm," allowing for ",Bnm)
        },
        {Anm <- preds[1]
         Bnm <- preds[2]
         if(preds[3] != paste0(Anm,":",Bnm)) {
             stop("Argument \"fmla\" is of an incorrect form.\n")
         }
         if(permtype != "stdres") {
             stop("Must use permtype = \"stdres\", in a model with interaction.\n")
         }
         type  <- "interac"
         effNm <- paste0("interaction of ",Anm," with ",Bnm)
        }
    )
    # Initial (real) data:
    sumFns <- initPrep(data,rspNm=rspNm,Anm=Anm,Bnm=Bnm,sumFnNm=sumFnNm,
                       type=type,expo=expo,rsteps=rsteps,r=r)
    
    # Calculate the statistic.
    Tobs <- testStat(sumFns,divByVar=divByVar)
    if(!test) {
       if(brief) {
           rslt <- list(EffectName=effNm,stat=Tobs)
       } else {
           rslt <- list(EffectName=effNm,fmla=fmla,sumFnNm=sumFnNm,stat=Tobs)
           if(keepdata) {
               rslt <- c(rslt,list(data=data))
           }
       }
       class(rslt) <- "kanova"
       return(rslt)
    }
    
    # Testing;  carry out the Monte Carlo test.
    # If permtype is "stdres", create the fitted values and residuals.
    if(permtype=="stdres") {
       rAndF <- resAndFit(sumFns) # List with components "resids" and "fitVals".
    } else {
       rAndF <- NULL
    }
    
    Tstar <- numeric(nperm)
    for(i in 1:nperm){
        pSumFns  <- permSumFns(sumFns,rAndF,permtype)
        Tstar[i] <- testStat(pSumFns,divByVar=divByVar)
        if(verb) cat(i,"")
        if(verb & i%%10 == 0) cat("\n")
    }
    if(verb & i%%10 != 0) cat("\n")
    m    <- sum(Tstar >= Tobs)
    pv   <- (m+1)/(nperm+1)
    bres <- list(EffectName=effNm,stat=Tobs,pvalue=pv) # brief result
    if(brief) {
        rslt <- bres
    } else {
        rslt <- c(bres,list(nperm=nperm,permtype=permtype,Tstar=Tstar,
                            fmla=fmla,sumFnNm=sumFnNm))
       if(keepdata) {
           rslt <- c(rslt,list(data=data))
       }
    }
    class(rslt) <- "kanova"
    rslt
    }
}
