wtdMean <- function(x) {
    wts <- getWts(x)
    wm  <- Reduce("+",lapply(x,function(u){u*attr(u,"weight")}))/sum(wts)
    attr(wm,"weight") <- NULL
    wm
}
