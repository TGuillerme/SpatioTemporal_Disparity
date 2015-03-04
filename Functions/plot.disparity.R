##########################
#Plotting disparity results
##########################
#Plots the disparity results
#v0.1
##########################
#SYNTAX :
#<disparity> disparity data
#<measure> the name of the column containing the disparity measurement. If set to 'default' the measure will be the first measure (second column) of the table.
#<rarefaction> whether to plot the rarefaction results or not
##########################
#----
#guillert(at)tcd.ie 03/03/2015
##########################

plot.disparity<-function(disparity_data, measure="default", rarefaction=FALSE, ...){
    #SANITIZING
    #Disparity
    check.class(disparity_data, 'data.frame', " must be a disparity data.frame")
    if(length(disparity_data) < 4) {
        stop("Disparity data.frame must have in the following order:\na 'rarefaction' column, a 'measurement' column and at least two 'Confidence Interval' columns.\nUse the disparity() function to generate the proper formatted data.frame.")
    }

    #Measure
    check.class(measure, 'character', " must be 'default' or one of the names of the columns in the disparity data.frame.")
    check.length(measure, 1, " must be 'default' or one of the names of the columns in the disparity data.frame.", errorif=FALSE)
    #Get the right column number
    if(measure == 'default') {
        measure_col<-2
    } else {
        measure_col<-grep(measure, colnames(disparity_data))
        if(length(measure_col) != 1) {
            stop("measure column not found in disparity_data.\nUse the disparity() function to generate the proper formatted data.frame.")
        }
    }
    #Get the Confidence intervals columns
    measure_col_tmp<-measure_col
    while(length(grep("%", colnames(disparity_data)[(measure_col_tmp+1)]))==1) {
        measure_col_tmp<-measure_col_tmp+1
    }
    #Get the column position for the different CI_values
    CI_length<-(measure_col_tmp-measure_col)
    CI_min<-measure_col+1
    CI_max<-measure_col+CI_length
    #Get the eventual intermediate CIs
    #number of CI values
    n_CI<-CI_length/2
    #extracting all the CI column number pairs in a table
    CI_pairs<-matrix(nrow=n_CI, ncol=2)
    for(n in 1:n_CI) {
        CI_pairs[n,1]<-CI_min+n-1
        CI_pairs[n,2]<-CI_max-n+1
    }
    #Setting the line types for the CIs
    lty_list<-c(44,33,22,21,12)

    #Rarefaction
    check.class(rarefaction, 'logical', " must be logical.")
    #Check if rarefaction data is available
    options(warn=-1)
    if(rarefaction == TRUE & is.na(disparity_data[,1])) {
        stop("No rarefaction data available.\nUse the disparity() function to generate the proper formatted data.frame with the option 'rarefaction=TRUE'.")
    }
    #If no rarefaction, check if disparity data is > 1
    if(rarefaction == FALSE & nrow(disparity_data) < 2) {
        warning("Only one disparity point is available.")
    }
    options(warn=0)
    #Plotting the disparity results
    if(rarefaction == TRUE) {
        #Plotting the rarefaction curve
        plot(disparity_data[,1], disparity_data[,measure_col], type='l', ylim=c(min(disparity_data[,CI_min]),max(disparity_data[,CI_max])) , ...)
        #Add the CIs
        for (n in 1:(CI_length/2)) {
            #Add both lines
            lines(disparity_data[,1], disparity_data[,CI_pairs[n,1]], type='l', lty=lty_list[n+1])
            lines(disparity_data[,1], disparity_data[,CI_pairs[n,2]], type='l', lty=lty_list[n+1])
        }

    } else {
        #Plotting the disparity curve
        if(nrow(disparity_data) == 1) {
            #If only one data point is available, do box plot style
            plot(1,1, xlab='', ylab='', ylim=c(min(disparity_data[,CI_min]),max(disparity_data[,CI_max])), type='n', xaxt='n')
            points(1,disparity_data[,measure_col], pch=19)
            #line types for this one
            lty_list2<-c(44,1,1,1,1,1)
            for (n in 1:(CI_length/2)) {
                #Add CIs lines
                lines(c(1,1), c(disparity_data[,CI_pairs[n,1]], disparity_data[,CI_pairs[n,2]]), lwd=1+(n-1)*3, lty=lty_list2[n])
            }
        } else {
            #Plotting the curve
            plot(seq(from=1, to=nrow(disparity_data)), disparity_data[,measure_col], type='l', ylim=c(min(disparity_data[,CI_min]),max(disparity_data[,CI_max])) ,col="white", ...)
            #Set the polygon colors
            polygon_colors<-c("lightgrey", "grey")
            #Add the polygons
            for (n in 1:(CI_length/2)) {
                polygon(c(seq(from=1, to=nrow(disparity_data)), seq(from=nrow(disparity_data), to=1)),
                    c(disparity_data[,CI_pairs[n,1]], rev(disparity_data[,CI_pairs[n,2]])),
                    col=polygon_colors[n], density=100-(100/(CI_length/2)/2.5*n))
            }
            #Add the central tendency line
            lines(seq(from=1, to=nrow(disparity_data)), disparity_data[,measure_col], type='l', ylim=c(min(disparity_data[,CI_min]),max(disparity_data[,CI_max])))
        }
    }
}