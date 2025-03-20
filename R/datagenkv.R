datagenkv <- function(av,bv=NULL,interac=NULL,nrep=5,rlen=129,sigma=1,seed=NULL,
                      pseudonoise=NULL){
#
# Function to generate examples of artificial data in the form of numeric
# vectors that are designed to "look like" K-functions.
#
if(is.null(seed)) seed <- sample(1:1e5,1)
set.seed(seed)
if(!requireNamespace("spatstat.geom"))
    stop("Required package \"spatstat.geom\" is not available.\n")
if(!(is.null(pseudonoise) | length(pseudonoise) == rlen))
    stop("Argument \"pseudonoise\" must have length \"rlen\".\n")

if(is.null(interac)) {
    if(is.null(av)) {
        stop("If the interaction is NULL, then \"av\" must be supplied.\n")
    }
    if(is.null(bv)) {
# One-way
        if(length(sigma) != 1) {
        whinge <-paste0("The length of \"sigma\" is ",length(sigma),"; it should\n",
                         "  be 1.\n")
        } 
        a <- length(av)
        b <- 1
        A  <- factor(rep(letters[1:a],nrep))
        B  <- factor(rep(1,a*nrep))
        mn <- rep(av,nrep)
        names(mn) <- A
    } else {
# Two-way
        a <- length(av)
        b <- length(bv)
        if(!(length(sigma) %in% c(1,b))) {
            whinge <-(paste0("The length of \"sigma\" is ",length(sigma),"; it should\n",
                         "  be either 1 or ",b,".\n"))
            stop(whinge)
        }
        A   <- factor(rep(letters[1:a],b*nrep))
        B   <- factor(rep(1:b,each=a,nrep))
        rav <- rep(av,b*nrep)
        rbv <- rep(bv,each=a,nrep)
        mn  <- rav + rbv
        names(mn) <- interaction(A,B)
        }
    } else {
# Model with interaction
    if(!inherits(interac,"matrix")) {
        stop("Argument \"interac\", if specified, must be a matrix.\n")
    }
    a <- nrow(interac)
    b <- ncol(interac)
    A <- factor(rep(letters[1:a],b*nrep))
    B <- factor(rep(1:b,each=a,nrep))
    if(!(length(sigma) %in% c(1,a*b))) {
        whinge <-(paste0("The length of \"sigma\" is ",length(sigma),"; it should\n",
                         "  be either 1 or ",a*b,".\n"))
        stop(whinge)
    }
    mn <- rep(as.vector(interac),nrep)
    names(mn) <- interaction(A,B)
}
sigma           <- matrix(sigma,nrow=a,ncol=b)
dimnames(sigma) <- list(levels(a),levels(b))

if(rlen > 1) {
    rrr   <- seq(0,0.25,length=rlen)
    kbase <- rrr^2
} else {
    rrr   <- 1
    kbase <- 0
}
rslt <- as.hyperframe(matrix(nrow=0,ncol=4))
names(rslt) <- c("y","rep","A","B")
lmv  <- length(mn) # length of mean vector
krow <- 0
for(k in 1:nrep) {
    for(i in 1:a) {
        for(j in 1:b) {
            krow <- krow + 1
            ko   <- if(krow%%lmv == 0) lmv else krow%%lmv
            if(is.null(pseudonoise)) {
                sbase <- sigma[i,j]*(if(rlen>1) kbase else 1)
                noise <- stats::rnorm(rlen,0,sbase)
            } else {
                noise <- pseudonoise
            }
            xxx  <- noise+mn[ko]+kbase
            yyy  <- spatstat.geom::hyperframe(y=list(pmax(xxx,0)),
                                              rep=k,A=A[krow],B=B[krow])
            rslt <- rbind(rslt,yyy)
        }
    }
}
row.names(rslt) <- 1:nrow(rslt)
if(rlen > 1) attr(rslt,"r") <- rrr
attr(rslt,"seed") <- seed
asig <- array(NA,dim=c(a,b,rlen))
for(k in 1:rlen) asig[,,k] <- sigma
attr(rslt,"sigma") <- asig
rslt
}
