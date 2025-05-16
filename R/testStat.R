testStat <- function(sumFns,divByVar=divByVar) {
# Test statistic.
type <- attr(sumFns,"type")
#
if(type == "oneway") {
    return(oEngine(sumFns,divByVar=divByVar))
} else if(type == "twoway") {
    return(tEngine(sumFns,divByVar=divByVar))
} else if(type=="interac") {
    return(iEngine(sumFns,divByVar=divByVar))
} else {
    stop(paste0("Value of \"type\" ",type," not recognised.\n"))
}
}
