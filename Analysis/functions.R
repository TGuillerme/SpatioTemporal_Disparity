#library(devtools)
#install("~/PhD/Projects/SpatioTemporal_Disparity/Functions/STD")
#library(STD)

#sourceDir function (from man source)
sourceDir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
        source(file.path(path, nm), ...)
        if(trace) cat("\n")
    }
}

sourceDir(path="../Functions/")
sourceDir(path="../Functions/test/")