datagenpp <- function(kapmin=30,kapmax=45,nkap=6,scale=0.1,
                      lambda=300,nrep=10,mdlEff=TRUE,seed=NULL) {
#
# Function datagenpp to generate a hyperframe with respons column
# equal to a list of point patterns from Thomas, Matern and Cauchy
# cluster models with varying values of the "kappa" parameter.
# (The other parameters are held fixed or determined from the kappa
# parameter and the fixed parameters.)  The predictor columns
# correspond to the kappa values and to the models.
#
# To generate data with no kappa effect, choose kapmin and kapmax
# to be equal.
# To generate data with no model effect, set mdlEff=FALSE.

if(kapmax < kapmin)
    stop("Argument \"kapmax\" must be at least as large as \"kapmin\".\n")
if(!requireNamespace("spatstat.random"))
    stop("Required package \"spatstat.random\" is not available.\n")
if(!requireNamespace("spatstat.geom"))
    stop("Required package \"spatstat.geom\" is not available.\n")
if(is.null(seed)) seed <- sample(1:1e5)
y    <- vector("list",nkap)
kvec <- seq(kapmin,kapmax,length=nkap)
set.seed(seed)
if(mdlEff) {
   rfun1 <- spatstat.random::rThomas
   rfun2 <- spatstat.random::rMatClust
   rfun3 <- spatstat.random::rCauchy
} else {
   rfun1 <- rfun2 <- rfun3 <- spatstat.random::rThomas
}
for(k in 1:nkap) {
    mu <- lambda/kvec[k]
    y1 <- rfun1(kappa=kvec[k],scale=scale,mu=mu,nsim=nrep,drop=FALSE)
    y2 <- rfun2(kappa=kvec[k],scale=scale,mu=mu,nsim=nrep,drop=FALSE)
    y3 <- rfun3(kappa=kvec[k],scale=scale,mu=mu,nsim=nrep,drop=FALSE)
    y[[k]]  <- c(y1,y2,y3)
}
yy <- do.call(c,y)
X  <- spatstat.geom::hyperframe(y=yy,kappa=factor(rep(letters[1:nkap],each=3*nrep)),
                 model=factor(rep(c("T","M","C"),each=nrep,nkap)))
attr(X,"seed") <- seed
X
}
