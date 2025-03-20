buildM1 <- function(Khat,Khati,Khatj,Khatij) {
    a <- length(Khati)
    b <- length(Khatj)
    Khatij <- matrix(Khatij,nrow=a,ncol=b)
    M1 <- array(dim=c(a,b,length(Khat)))
    for(i in 1:a) {
        for(j in 1:b) {
            part1 <- Khatij[[i,j]]
            part2 <- Khati[[i]] #apply(do.call(cbind,Khati[[i]]),1,sum)
            part3 <- Khatj[[j]] #apply(do.call(cbind,Khatj[[j]]),1,sum)
            M1[i,j,] <- part1 - part2 - part3 + Khat
        }
    }
    M1
}
