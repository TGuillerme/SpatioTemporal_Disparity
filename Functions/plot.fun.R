plot.pco<-function(pco.scores ,clade.1 , clade.2, ...) {
    library(grDevices)
    plot(1,1, col="white", xlab="PC1", ylab="PC2", xlim=c(min(pco.scores[,1]), max(pco.scores[,1])), ylim=c(min(pco.scores[,2]), max(pco.scores[,2])), ...)
    points(pco.scores[,1], pco.scores[,2], col=c("red", "blue")[as.factor(pco.scores[,clade.col])])
    polygon(clade.1[chull(clade.1),], border="red")
    polygon(clade.2[chull(clade.2),], border="blue")
}

pcoSlice<-function(tree, main="slice") {
    pcoSlice<-pco.scores[tree$tip.label,]
    clade.1<-pcoSlice[match(clade1, row.names(pcoSlice)),1:2]
    clade.2<-pcoSlice[match(clade2, row.names(pcoSlice)),1:2]
    plot.pco(pcoSlice, clade.1[-which(is.na(clade.1)),], clade.2[-which(is.na(clade.2)),], main=main)
}
