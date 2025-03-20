reenlist <- function(x,f) {
if(!is.vector(x)) stop("Argument \"x\" must be a vector (possibly a list).\n")
f <- as.factor(f)
if(is.null(names(x)))
    stop("Argument \"x\" must be a named vector (list).\n")
if(!all(names(x) %in% levels(f)))
    stop("All names of \"x\" must appear amongst the levels of \"f\".\n")
x[match(f,names(x))]
}
