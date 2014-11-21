#Installing packages if necessary
list.of.packages <- c("diversitree", "paleotree", "strap", "ape", "caper", "phytools", "vegan", "grDevices", "plotrix", "testthat", "devtools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#Load the packages
suppressMessages({
    library(diversitree)
    library(devtools)
    library(paleotree)
    library(strap)
    library(ape)
    library(caper)
    library(phytools)
    library(vegan)
    library(grDevices)
    library(plotrix)
    library(testthat)
})

#Install Claddis from github

list.of.packages <- c("Claddis")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install_github("graemetlloyd/Claddis")

suppressMessages({
library(Claddis)
})