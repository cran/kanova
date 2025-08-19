buildV <- function(s2ij,wtm){
    a  <- nrow(wtm)
    b  <- ncol(wtm)
    nr <- dim(s2ij)[3]
    V  <- array(NA,dim=dim(s2ij))
    Wi <- wtm/apply(wtm,1,sum)
    Wj <- t(wtm)/apply(wtm,2,sum)
    Z  <- lapply(1:nr,function(k,a,b){a[,,k]*b},a=s2ij,b=wtm/sum(wtm))
    Z  <- array(unlist(Z),dim=dim(s2ij))
    Z  <- apply(Z,3,sum)
    for(i in 1:a) {
        si. <- Wi %*% s2ij[i,,]
        for(j in 1:b) {
            s.j <- Wj %*% s2ij[,j,]
            mult1 <- 1/wtm[i,j] - 2/sum(wtm[i,]) - 2/sum(wtm[,j]) +
                     2*wtm[i,j]/(sum(wtm[i,])*sum(wtm[,j])) + 2/sum(wtm)
            part1 <- mult1*s2ij[i,j,]
            mult2 <- 1/sum(wtm[i,]) - 2/sum(wtm)
            part2 <- mult2*si.[i,]
            mult3 <- 1/sum(wtm[,j]) - 2/sum(wtm)
            part3 <- mult3*s.j[j,]
            mult4 <- 1/sum(wtm)
            part4 <- mult4*Z
            V[i,j,] <- part1 + part2 + part3 + part4
       }
    }
    V
}
