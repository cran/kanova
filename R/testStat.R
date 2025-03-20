testStat <- function(sumFns) {
# Test statistic.
type <- attr(sumFns,"type")
#
if(type == "oneway") {
    return(oEngine(sumFns))
} else if(type == "twoway") {
    return(tEngine(sumFns))
} else if(type=="interac") {
    return(iEngine(sumFns))
} else {
    stop(paste0("Value of \"type\" ",type," not recognised.\n"))
}
}
