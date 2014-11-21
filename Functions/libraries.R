#Installing packages if necessary
list.of.packages <- c("diversitree", "Claddis", "paleotree", "strap", "ape", "caper", "phytools", "vegan", "grDevices", "testthat")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#Load the packages
suppressMessages({
    library(diversitree)
    library(Claddis)
    library(paleotree)
    library(strap)
    library(ape)
    library(caper)
    library(phytools)
    library(vegan)
    library(grDevices)
    library(testthat)
})