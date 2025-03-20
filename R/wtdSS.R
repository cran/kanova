wtdSS <- function(x) {
    mn   <- wtdMean(x)
    n    <- length(x)
    sqdiff <- vector("list",n)
    for(i in 1:n) {
        sqdiff[[i]] <- (x[[i]] - mn)^2*attr(x[[i]],"weight")
    }
    rslt <- Reduce("+",sqdiff)
    attr(rslt,"n") <- n
    rslt
}
