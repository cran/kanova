wtdVar <- function(x) {
    mn   <- wtdMean(x)
    n    <- length(x)
    sqdiff <- vector("list",n)
    for(i in 1:n) {
        sqdiff[[i]] <- (x[[i]] - mn)^2
        attr(sqdiff[[i]],"weight") <- attr(x[[i]],"weight")
    }
    wtdMean(sqdiff)
}
