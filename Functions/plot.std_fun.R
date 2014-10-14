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
plot.pco<-function(data, legend, pos.leg, xlim, ylim, col, ...) {
    
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
        #Draw Convec Hull if clade > 1 (2 coordinates)
        if(length(clade) > 2) {
            polygon(clade[chull(clade),], border=palette()[group])
        }
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
pco.slice<-function(data, legend, pos.leg, xlim, ylim, col, pars, ...) {
    
    #main (default/missing)
    if(!exists("main")){
        main=paste("Mya", data$method, sep=" - ")
    }

    #xlim (default)
    if(xlim=='default') {
        X<-NULL
        for(slice in 1:(length(data)-1)){ #-1 for method
            X<-c(X, data[[slice]][[3]][,1])
        }
        xlim=c(min(X), max(X))
    }

    #ylim (default)
    if(ylim=='default') {
        Y<-NULL
        for(slice in 1:(length(data)-1)){ #-1 for method
            Y<-c(Y, data[[slice]][[3]][,2])
        }
        ylim=c(min(Y), max(Y))
    }   

    #pars
    if(missing(pars)) {
        stop("Plotting layout parameters must be defined.")
    } else {
        check.class(pars, 'numeric', ' must be a vector of two numerical values.')
        check.length(pars, 2, ' must be a vector of two numerical values.')
    }

    #plot window setting
    op<-par(mfrow=pars)

    #Plotting the first slice (optional legend)
    #Creating pco.scores object
    data_tmp<-data[[1]][3:4]
    main_tmp<-paste(data[[1]][[1]], main)
    plot.pco(data_tmp, legend, pos.leg, xlim, ylim, col, main=main_tmp, ...)

    if (length(data)-1 > 1) {
        for (slice in 3:length(data)-1) {
            data_tmp<-data[[slice]][3:4]
            main_tmp<-paste(data[[slice]][[1]], main)
            plot.pco(data_tmp, legend=FALSE, pos.leg, xlim, ylim, col, main=main_tmp, ...)
        }
    }

    par(op)
}
