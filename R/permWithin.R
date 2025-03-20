permWithin <- function(G) {
n <- length(G)
splind <- split(1:n,f=G)
ipl  <- lapply(splind,sample)
unsplit(ipl,G)
}
