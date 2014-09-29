
plot.pco<-function(pco.scores, tax.col, legend=FALSE, xlim='default', ylim='default',...) {
    library(grDevices)
    #empty plot
    if(xlim=='default') {
        xlim=c(min(pco.scores[,1]), max(pco.scores[,1]))
    }
    if(ylim=='default') {
        ylim=c(min(pco.scores[,2]), max(pco.scores[,2]))
    }
    plot(1,1, col="white", xlab="PC1", ylab="PC2", xlim=xlim, ylim=ylim, ...)
    groups=length(levels(as.factor(pco.scores[,tax.col])))
    #points
    points(pco.scores[,1], pco.scores[,2], col=palette()[1:groups][as.factor(pco.scores[,tax.col])])
    #polygon
    for (group in 1:groups) {
        clade<-pco.scores[which(pco.scores[,tax.col] == levels(as.factor(pco.scores[,tax.col]))[group]),1:2]
        polygon(clade[chull(clade),], border=palette()[group])
    }
    if(legend == TRUE){
        legend(min(pco.scores[,1]), max(pco.scores[,2]), levels(as.factor(pco.scores[,tax.col])), col=palette()[1:groups], pch=21, cex=0.7)
    }
}

#fix the scale
pco.slice<-function(tree, pco.scores, slices, method, tax.col, legend=FALSE, pars=c(3,3), xlim, ylim, ...) {
    #Setting the age slices (can be one value (time is split equitably) or a vector of values containing the age)
    if(length(slices) == 1) {
        age<-seq(from=0, to=max(tree.age(tree)[,1]), length.out=slices)
    } else {
        age<-slices
        slices<-length(slices)
    }

    #plot window setting
    op<-par(mfrow=pars) 

    #calculating and plotting each slice
    for (slice in 1:slices) {
        subtree<-slice.tree(tree, age=age[slice], method)
        pco.scores[,tax.col]<-as.factor(pco.scores[,tax.col])
        pco_slice<-subset(pco.scores[subtree$tip.label,])
        plot.pco(pco_slice, tax.col, legend, main=paste(round(age[slice]), "Mya"), xlim=xlim, ylim=ylim, ...)
    }
    par(op)
}

#the subplots must keep the same levels 
