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
plot.pco<-function(data, legend=FALSE, pos.leg, xlim='default', ylim='default', col, ...) {
    
    #empty plot
    suppressWarnings(
        if(xlim=='default') {
            xlim=c(min(data[[1]][,1]), max(data[[1]][,1]))
        }
    )
    suppressWarnings(
        if(ylim=='default') {
            ylim=c(min(data[[1]][,2]), max(data[[1]][,2]))
        }
    )
    plot(1,1, col="white", xlim=xlim, ylim=ylim, ...)

    groups=length(levels(as.factor(data[[2]][,1])))

    #points
    points(data[[1]][,1], data[[1]][,2], col=col[1:groups][as.factor(data[[2]][,1])])

    #polygon
    for (group in 1:groups) {
        n<-which(data[[2]][,1] == levels(as.factor(data[[2]][,1]))[group])
        clade<-data[[1]][n,]
        #which(row.names(data[[1]]) == levels((as.factor(data[[2]][,1]))[group]),1:2]
        polygon(clade[chull(clade),], border=palette()[group])
    }
    if(legend == TRUE){
        if(pos.leg == 'default'){
            legend(min(data[[1]][,1]), max(data[[1]][,2]), levels(as.factor(data[[2]][,1])), col=col[1:groups], pch=21, cex=0.7)
        } else {
            legend(pos.leg[1], pos.leg[2], levels(as.factor(data[[2]][,1])), col=col[1:groups], pch=21, cex=0.7)
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
