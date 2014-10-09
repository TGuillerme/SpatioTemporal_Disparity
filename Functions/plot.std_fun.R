#plotting the pco axis load
plot.load<-function(data, legend, pos.leg, pars=c(1,2), ...) {
    #Window settings
    op<-par(mfrow=pars) 

    #Selecting the axis load and their cumulative value from the pco
    cum<-data$values$Cum_corr_eig
    load<-data$values$Rel_corr_eig

    #plotting the load
    barplot(load, main="Relative variance per axis")
    if(legend==TRUE) {
        text(round(length(load)/2), (max(load)-0.1*max(load)), paste("1st axis = ", round(load[1]*100, digit=2), "% variance", sep=""), cex=0.5)
        text(round(length(load)/2), (max(load)-0.125*max(load)), paste("2nd axis = ", round(load[2]*100, digit=2), "% variance", sep=""), cex=0.5)
        text(round(length(load)/2), (max(load)-0.15*max(load)), paste("3rd axis = ", round(load[3]*100, digit=2), "% variance", sep=""), cex=0.5)
    }

    #plotting the cumulative variance
    barplot(cum, main="Cumulative variance per axis")
    if(legend==TRUE) {
        abline(0.95,0)
        text(round(length(load)/3), 0.95, paste("0.95 cumulative variance (", length(which(cum <= 0.95)), " axis)", sep=""), cex=0.5 , pos=1)
    }

    par(op)
}

#plotting a global pco
plot.pco<-function(pco.scores, tax.col, legend=FALSE, pos.leg, xlim='default', ylim='default',...) {
    
    #empty plot
    suppressWarnings(
        if(xlim=='default') {
            xlim=c(min(pco.scores[,1]), max(pco.scores[,1]))
        }
    )
    suppressWarnings(
        if(ylim=='default') {
            ylim=c(min(pco.scores[,2]), max(pco.scores[,2]))
        }
    )
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
        if(missing(pos.leg)){
            legend(min(pco.scores[,1]), max(pco.scores[,2]), levels(as.factor(pco.scores[,tax.col])), col=palette()[1:groups], pch=21, cex=0.7)
        } else {
            legend(pos.leg[1], pos.leg[2], levels(as.factor(pco.scores[,tax.col])), col=palette()[1:groups], pch=21, cex=0.7)
        }
    }
}

#plotting a pco through time
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
