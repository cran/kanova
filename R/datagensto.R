datagensto <- function(nrep=10,sigma=1,interac=0,seed=NULL,
                       meansonly=FALSE,perturbLayer=0) {
#
# Function datagensto(), to generate/simulate "artificial" point
# pattern data, based on the stomata data.  The results take the
# form of pseudo K-functions, structured so that an additive model
# makes sense.
#

stomata <- kanova::stomata
X <- stomata[,"patterns",drop=TRUE]
Y <- lapply(X,function(x){spatstat.geom::affine(x,mat=diag(1/c(1200,900)))})
r <- seq(0,0.25,length=129)
Kays  <- lapply(Y,function(y,r){kkk <- Kest(y,r=r)
                                wt <- npoints(y)
                                rslt <- kkk[[attr(kkk,"valu")]]
                                attr(rslt,"weight") <- wt
                                rslt
                                },r=r)
K.Pos <- split(Kays,f=stomata$Pos)
K.Lay <- split(Kays,f=stomata$Layer)
base.Pos <- lapply(K.Pos,wtdMean)
base.Lay <- lapply(K.Lay,wtdMean)
if(meansonly) {
   return(list(Pos=base.Pos,Layer=base.Lay,r=r))
}
if(!requireNamespace("spatstat.geom"))
    stop("Required package \"spatstat.geom\" is not available.\n")
if(is.null(seed)) seed <- sample(1:1e5,1)
set.seed(seed)
N        <- nrep*3*6
simKays  <- vector("list",N)
Pos      <- numeric(N)
Layer    <- numeric(N)
sigma    <- 1
kount    <- 0
for(i in 1:3) {
    for(j in 1:6) {
        xxx <- base.Pos[[i]] + base.Lay[[j]]+j*perturbLayer*r^2
        for(k in 1:nrep) {
            kount <- kount + 1
            Pos[[kount]] <- i
            Layer[[kount]] <- j
            yyy <- xxx + rnorm(129,0,sigma*r^2)
            simKays[[kount]] <- pmax(0,yyy)
            if(i==3 & j==6 & interac>0) {
                simKays[[kount]] <- simKays[[kount]] + interac*r
            }
        }
    }
}
Pos   <- factor(Pos)
Layer <- factor(Layer)
Xsim  <- hyperframe(y=simKays,Pos=Pos,Layer=Layer)
attr(Xsim,"r") <- r
attr(Xsim,"seed") <- seed
Xsim
}
