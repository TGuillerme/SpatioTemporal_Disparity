#Library
setwd('~/PhD/Projects/SpatioTemporal_Disparity/Analysis')
library(ape)

#Plotting the polygon outlines per type

plot.pco<-function(pco.scores ,clade.1 , clade.2, ...) {
    library(grDevices)
    plot(1,1, col="white", xlab="PC1", ylab="PC2", xlim=c(min(pco.scores[,1]), max(pco.scores[,1])), ylim=c(min(pco.scores[,2]), max(pco.scores[,2])), ...)
    points(pco.scores[,1], pco.scores[,2], col=c("red", "blue")[as.factor(pco.scores[,clade.col])])
    polygon(clade.1[chull(clade.1),], border="red")
    polygon(clade.2[chull(clade.2),], border="blue")
}

source('treeAge.R')

#Data
Slater.table<-read.dna('../Empirical_mammals/Matrices/2013-Slater-MEE-morpho.phylip', as.character=TRUE, as.matrix=TRUE)