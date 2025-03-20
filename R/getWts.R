getWts   <- function(x){unlist(sapply(x,function(y){attr(y,"weight")}))}
