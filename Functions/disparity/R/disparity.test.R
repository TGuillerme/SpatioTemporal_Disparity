##########################
#Disparity testing function
##########################
#Calculates the differences between PCO intervals based on Anderson & Friedman 2012 test.
#v0.0.1
##########################
#SYNTAX :
#<time_pco> a time_pco matrix
#<method> the method for calculating the disparity. Can be any of the following: "volume", "centroid", "sum.range", "product.range", "sum.variance", "product.variance"
#<bootstraps> the number of bootstrap replicates (default=1000).
#<correction> one of the methods from p.adjust function to correct the p-values. If NULL, no correction will be applied.
#<rarefaction> a rarefaction value. If NULL, is ignored.
#<verbose> whether to be verbose or not
#<rm.last.axis> whether to remove the last axis from the pco time_pco. Can be a threshold value.
##########################
#----
#guillert(at)tcd.ie 09/07/2015
##########################

disparity.test<-function(time_pco, method, bootstraps=1000, correction="bonferroni", rarefaction=NULL, verbose=FALSE, rm.last.axis=FALSE) {
    #-----------------------------
    #SANITIZING
    #-----------------------------
    #distance
    check.class(time_pco, "matrix")

    #Test if applicable (> 2 rows)
    if(nrow(time_pco) < 2) {
        stop("Disparity can not be calculated because less than two taxa are present in the time_pco!")
    } 

    #method
    check.class(method, "character", " can be 'volume', 'centroid', 'sum.range', 'product.range', 'sum.variance' or/and 'product.variance'.")
    methods_list<-c('volume', "centroid", "sum.range", "product.range", "sum.variance", "product.variance")
    check.length(method, 1, ": only one method can be input")
    if(all(is.na(match(method, methods_list)))) {
        stop("method can be 'volume', 'centroid', 'sum.range', 'product.range', 'sum.variance' or/and 'product.variance'.")
    }

    #Bootstrap
    check.class(bootstraps, "numeric", " must be a single (entire) numerical value.")
    check.length(bootstraps, 1, " must be a single (entire) numerical value.")
    #Make sure the bootstrap is a whole number
    bootstraps<-round(abs(bootstraps))
    
    #Central tendency
    central_tendency=mean
    
    #rarefaction
    if(is.null(rarefaction)) {
        rarefaction<-FALSE
    } else {
        check.class(rarefaction, "numeric")
    }

    #verbose
    check.class(verbose, "logical")

    #rm.last.axis
    if(class(rm.last.axis) == "logical") {
        if(rm.last.axis == FALSE) {
            rm.axis<-FALSE
        } else {
            rm.axis<-TRUE
            last.axis<-0.95
        }
    } else {
        check.class(rm.last.axis, "numeric", " must be logical or a probability threshold value.")
        check.length(rm.last.axis, 1, " must be logical or a probability threshold value.", errorif=FALSE)
        if(rm.last.axis < 0) {
            stop("rm.last.axis must be logical or a probability threshold value.")
        } else {
            if(rm.last.axis > 1) {
                stop("rm.last.axis must be logical or a probability threshold value.")
            } else {
                rm.axis<-TRUE
                last.axis<-rm.last.axis
            }
        }
    }

    #centroid - LEAVE OR NOT?
    centroid.type_function<-cen.apply.mea
    
    #-----------------------------
    #CLEANING / BOOTSTRAPING
    #-----------------------------

    #Removing the last pco axis
    if(rm.axis==TRUE) {
        #calculate the cumulative variance per axis
        scree_data<-cumsum(apply(time_pco, 2, var) / sum(apply(time_pco, 2, var)))
        #extract the axis  below the threshold value
        axis_selected<-length(which(scree_data < last.axis))
        #remove the extra axis
        time_pco<-time_pco[,c(1:axis_selected)]
        #warning
        message(paste("The", length(scree_data)-axis_selected, "last axis have been removed from the pco time_pco."))
    }

    #Bootstraping the matrix
    #verbose
    if(verbose==TRUE) {
        message("Bootstraping...", appendLF=FALSE)
    }
    BSresult<-lapply(time_pco, Bootstrap.rarefaction, bootstraps, rarefaction)
    if(verbose==TRUE) {
        message("Done.", appendLF=TRUE)
    }

    #Selecting the method function
    if(method == 'volume') {
        method.fun<-volume.calc
    }
    if(method == 'centroid') {
        method.fun<-centroid.calc
    }
    if(method == 'range') {
        method.fun<-range.calc
    }
    if(method == 'variance') {
        method.fun<-variance.calc
    }

    #sum of the sums?

    BS_tmp<-matrix(NA, nrow=bootstraps, ncol=length(time_pco))
    for (int in 1:length(time_pco)) {
        BS_tmp[,int]<-apply(lapply(BSresult[[int]],centroid.calc)[[1]], 1, sum)
    }
    BS_result<-BS_tmp




    #-----------------------------
    #VOLUME
    #-----------------------------
    #Hyperspace volume calculation
    if(any(method == 'volume')) {
        #Calculate the hyperspace volume
        if(verbose==TRUE) {
            message("Calculating hyperspace volume...", appendLF=FALSE)
        }
        volumes<-lapply(BSresult, volume.calc)
        #Volumes table
        Volume_table<-Disparity.measure.table(type_function=no.apply, volumes, central_tendency, CI, save.all)
        #Results type
        if(save.all == FALSE) {
            #Renaming the column
            colnames(Volume_table)[1]<-"Volume"
        } else {
            #Isolating the time_pco parts
            Volume_values<-Volume_table[[2]]
            Volume_table<-Volume_table[[1]]
            colnames(Volume_table)[1]<-"Volume"
        }
        if(verbose==TRUE) {
            message("Done.", appendLF=TRUE)
        }
    }

    #-----------------------------
    #CENTROID
    #-----------------------------
    if(any(method == 'centroid')) {
        #Calculate the distance from centroid for the rarefaction and the bootstrapped matrices
        if(verbose==TRUE) {
            message("Calculating distance from centroid...", appendLF=FALSE)
        }
        centroids<-lapply(BSresult, centroid.calc)
        #Distance to centroid
        Centroid_dist_table<-Disparity.measure.table(type_function=centroid.type_function, centroids, central_tendency, CI, save.all)
        #Results type
        if(save.all == FALSE) {
            #Renaming the column
            colnames(Centroid_dist_table)[1]<-"Cent.dist"
        } else {
            #Isolating the time_pco parts
            Centroid_values<-Centroid_dist_table[[2]]
            Centroid_dist_table<-Centroid_dist_table[[1]]
            colnames(Centroid_dist_table)[1]<-"Cent.dist"
        }
        if(verbose==TRUE) {
            message("Done.", appendLF=TRUE)
        }
    }

    #-----------------------------
    #RANGES
    #-----------------------------
    #Wills 1994 range calculations
    if(any(grep("range", method))) {
        #Calculate the range for the rarefaction and the bootstrapped matrices
        if(verbose==TRUE) {
            message("Calculating ranges...", appendLF=FALSE)
        }
        ranges<-lapply(BSresult, range.calc)

        #Sum of ranges
        if(any(method == 'sum.range')) {
            #Sum of range
            Sum_range_table<-Disparity.measure.table(type_function=sum.apply, ranges, central_tendency, CI, save.all)

            #Results type
            if(save.all == FALSE) {
                #Renaming the column
                colnames(Sum_range_table)[1]<-"Sum.range"
            } else {
                #Isolating the time_pco parts
                Sum_range_values<-Sum_range_table[[2]]
                Sum_range_table<-Sum_range_table[[1]]
                colnames(Sum_range_table)[1]<-"Sum.range"
            }

        }

        #Product of ranges
        if(any(method == 'product.range')) {
            #Product of range
            Product_range_table<-Disparity.measure.table(type_function=prod.apply, ranges, central_tendency, CI, save.all)

            #Results type
            if(save.all == FALSE) {
                #Renaming the column
                colnames(Product_range_table)[1]<-"Prod.range"
            } else {
                #Isolating the time_pco parts
                Prod_range_values<-Product_range_table[[2]]
                Product_range_table<-Product_range_table[[1]]
                colnames(Product_range_table)[1]<-"Prod.range"
            }

        }
        if(verbose==TRUE) {
            message("Done.", appendLF=TRUE)
        }
    }

    #-----------------------------
    #VARIANCE
    #-----------------------------
    #Wills 1994 variance calculations
    if(any(grep("variance", method))) {
        #Calculate the variance for the rarefaction and the bootstrapped matrices
        if(verbose==TRUE) {
            message("Calculating variance...", appendLF=FALSE)
        }
        variances<-lapply(BSresult, variance.calc)

        #Sum of variance
        if(any(method == 'sum.variance')) {
            #Sum of variance
            Sum_variance_table<-Disparity.measure.table(type_function=sum.apply, variances, central_tendency, CI, save.all)

            #Results type
            if(save.all == FALSE) {
                #Renaming the column
                colnames(Sum_variance_table)[1]<-"Sum.var"
            } else {
                #Isolating the time_pco parts
                Sum_variance_values<-Sum_variance_table[[2]]
                Sum_variance_table<-Sum_variance_table[[1]]
                colnames(Sum_variance_table)[1]<-"Sum.var"
            }
  
        }

        #Product of variance
        if(any(method == 'product.variance')) {
            #Product of variance
            Product_variance_table<-Disparity.measure.table(type_function=prod.apply, variances, central_tendency, CI, save.all)

            #Results type
            if(save.all == FALSE) {
                #Renaming the column
                colnames(Product_variance_table)[1]<-"Prod.var"
            } else {
                #Isolating the time_pco parts
                Prod_variance_values<-Product_variance_table[[2]]
                Product_variance_table<-Product_variance_table[[1]]
                colnames(Product_variance_table)[1]<-"Prod.var"
            }
        
        }
        if(verbose==TRUE) {
            message("Done.", appendLF=TRUE)
        }
    }

    #-----------------------------
    #OUTPUT
    #-----------------------------
    #Empty output table
    if(rarefaction==FALSE) {
        output<-matrix(nrow=1, time_pco=rep(NA, 1))
        colnames(output)[1]<-"rarefaction"
    } else {
        output<-matrix(nrow=(nrow(time_pco)-2), time_pco=seq(from=3, to=nrow(time_pco)))
        colnames(output)[1]<-"rarefaction"
    }
    #Volume
    if(any(method == 'volume')) {
        #Combine the results
        output<-cbind(output, Volume_table)
    }
    #Distance form centroid
    if(any(method == 'centroid')) {
        #Combine the results
        output<-cbind(output, Centroid_dist_table)
    }
    #Sum of ranges
    if(any(method == 'sum.range')) {
        #Combine the results
        output<-cbind(output, Sum_range_table)
    }
    #Product of ranges
    if(any(method == 'product.range')) {
        #Combine the results
        output<-cbind(output, Product_range_table)
    }
    #Sum of variance
    if(any(method == 'sum.variance')) {
        #Combine the results
        output<-cbind(output, Sum_variance_table)  
    }
    #Product of variance
    if(any(method == 'product.variance')) {
        #Combine the results
        output<-cbind(output, Product_variance_table)   
    }

    if(save.all == FALSE) {
        #Quantiles only
        return(output)
    } else {
        #Quantiles and values
        output<-list("table"=output)
        #Add the values of each metric
        #centroid
        if(any(method == 'volume')) {
            output[[length(output)+1]]<-Volume_values
            names(output)[length(output)]<-"volume"
        }
        #centroid
        if(any(method == 'centroid')) {
            output[[length(output)+1]]<-Centroid_values
            names(output)[length(output)]<-"centroid"
        }
        #Sum of sum.range
        if(any(method == 'sum.range')) {
            output[[length(output)+1]]<-Sum_range_values
            names(output)[length(output)]<-"sum.range"
        }
        #Product of ranges
        if(any(method == 'product.range')) {
            output[[length(output)+1]]<-Prod_range_values
            names(output)[length(output)]<-"product.range"
        }
        #Sum of variance
        if(any(method == 'sum.variance')) {
            output[[length(output)+1]]<-Sum_variance_values
            names(output)[length(output)]<-"sum.variance"
        }
        #Product of variance
        if(any(method == 'product.variance')) {
            output[[length(output)+1]]<-Prod_variance_values
            names(output)[length(output)]<-"product.variance"
        }

        return(output)
    }

#End
}