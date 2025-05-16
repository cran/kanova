resAndFit <- function(sumFns) {
#
# Note:  The residuals are those from the saturated model.  The
# fitted values are those from under the appropriate null hypothesis.
#

N <- length(sumFns)
type <- attr(sumFns,"type")
switch(EXPR=type,
    oneway = { # Fitted values from one factor model; A only.  Null
               # hypothesis is that there is no A effect, so K_{ij}-hat = Khat.
        Khat  <- resAndFitCmpnts(sumFns)$Khat
        fitz  <- rep(list(Khat),N)
        A     <- attr(sumFns,"A")
        sfitz <- lapply(split(sumFns,f=A),wtdMean) # saturated "fitz"
        sfitz <- sfitz[match(A,names(sfitz))]
    },
    twoway = { # Fitted values from a two factor data structure. Null
               # hypothesis is that there is no A effect, so
               # K_{ijk}-hat = K_{.j.}-bar
        ufitz <- resAndFitCmpnts(sumFns)$Khatj
        B     <- attr(sumFns,"B")
        AB    <- attr(sumFns,"AB")
        fitz  <- reenlist(ufitz,B)
        names(fitz) <- AB
        sfitz <- lapply(split(sumFns,f=AB),wtdMean) # saturated "fitz"
        sfitz <- reenlist(sfitz,AB)
        names(sfitz) <- AB
    },
    interac = { # Fitted values from model with interaction, A * B.
                # Null hypothesisis that there is no interaction,
                # i.e. that the model is additive, whence the fitted values
                # are K_{ijk}-hat = K_{i..}-bar + K_{.j.}-bar - Khat
        xxx    <- resAndFitCmpnts(sumFns)
        ufitzA <- xxx$Khati
        ufitzB <- xxx$Khatj
        Khat   <- xxx$Khat
        A      <- attr(sumFns,"A")
        B      <- attr(sumFns,"B")
        fitzA  <- ufitzA[match(A,names(ufitzA))]
        fitzB  <- ufitzB[match(B,names(ufitzB))]
        fitz   <- lapply(1:N,function(i,a,b,c){a[[i]] + b[[i]] - c},
                             a=fitzA,b=fitzB,c=Khat)
        nms    <- sapply(1:N,function(i,A,B){paste(A[i],B[i],sep=".")},A=A,B=B)
        names(fitz) <- nms
        AB    <- attr(sumFns,"AB")
        sfitz <- lapply(split(sumFns,f=AB),wtdMean) # saturated "fitz"
        sfitz <- reenlist(sfitz,f=AB)
        names(sfitz) <- AB
    }
)
rez <- do.call(cbind,sumFns) - do.call(cbind,as.list(sfitz))
dimnames(rez) <- list(seq(along=attr(sumFns,"r")),names(sfitz))
list(resids=rez,fitVals=fitz)
}
