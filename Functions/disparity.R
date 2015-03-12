##########################
#Disparity functions
##########################
#Calculate the disparity as the distance from centroid
#This function is based on DisparityCalc() from Smith et al. 2014 - Evolution (http://dx.doi.org/10.1111/evo.12435) http://datadryad.org/resource/doi:10.5061/dryad.d380g 
#v0.1.1
##########################
#SYNTAX :
#<distance> the distance matrix
#<method> the method for calculating the disparity. Can be any of the following: "centroid", "sum.range", "product.range", "sum.variance", "product.variance"
#<CI> confidence intervals (default=c(50,95))
#<bootstraps> the number of boostrap replicates (default=1000)
#<central_tendency> any function for calculating the central tendency
#<verbose> whether to be verbose or not
##########################
#----
#guillert(at)tcd.ie 04/03/2015
##########################

disparity<-function(data, method=c("centroid", "sum.range", "product.range", "sum.variance", "product.variance"), CI=c(50, 95), bootstraps=1000, central_tendency=median, rarefaction=FALSE, verbose=FALSE) {

    #SANITIZING
    #distance
    check.class(data, "matrix", " must be a distance matrix.")

    #method
    check.class(method, "character", " must be 'centroid', 'sum.range', 'product.range', 'sum.variance' or/and 'product.variance'.")
    methods_list<-c("centroid", "sum.range", "product.range", "sum.variance", "product.variance")
    if(all(is.na(match(method, methods_list)))) {
        stop("method must be 'centroid', 'sum.range', 'product.range', 'sum.variance' or/and 'product.variance'.")
    }

    #CI
    check.class(CI, "numeric", " must be any value between 1 and 100.")
    #remove warnings
    options(warn=-1)
    if(any(CI) < 1) {
        stop("CI must be any value between 1 and 100.")
    }
    if(any(CI) > 100) {
        stop("CI must be any value between 1 and 100.")
    }
    options(warn=0)
    #Bootstrap
    check.class(bootstraps, "numeric", " must be a single (entire) numerical value.")
    check.length(bootstraps, 1, " must be a single (entire) numerical value.")
    #Make sure the bootstrap is a whole number
    bootstraps<-round(abs(bootstraps))

    #Central tendency
    check.class(central_tendency, "function", " must be either a function (e.g. 'mean' or 'median'.")

    #rarefaction
    check.class(rarefaction, "logical", " must be logical.")

    #verbose
    check.class(verbose, "logical", " must be logical.")

    #CALCULATING THE DISPARITY
    #Bootstraping the matrix
    #verbose
    if(verbose==TRUE) {
        message("Bootstraping...", appendLF=FALSE)
    }
    BSresult<-Bootstrap.rarefaction(data, bootstraps, rarefaction)
    if(verbose==TRUE) {
        message("Done.", appendLF=TRUE)
    }

    #CENTROID
    #Distance form centroid
    if(any(method == 'centroid')) {
        #Calculate the distance from centroid for the rarefaction and the bootstrapped matrices
        if(verbose==TRUE) {
            message("Calculating distance from centroid...", appendLF=FALSE)
        }
        centroids<-lapply(BSresult, centroid.calc)
        #Distance to centroid
        Centroid_dist_table<-Disparity.measure.table(type_function=no.apply, centroids, central_tendency, CI)
        #Renaming the column
        colnames(Centroid_dist_table)[1]<-"Cent.dist"
        if(verbose==TRUE) {
            message("Done.", appendLF=TRUE)
        }
    }

    #RANGES
    if(any(grep("range", method))) {
        #Calculate the range for the rarefaction and the bootstrapped matrices
        if(verbose==TRUE) {
            message("Calculating ranges...", appendLF=FALSE)
        }
        ranges<-lapply(BSresult, range.calc)

        #Sum of ranges
        if(any(method == 'sum.range')) {
            #Sum of range
            Sum_range_table<-Disparity.measure.table(type_function=sum.apply, ranges, central_tendency, CI)
            #Renaming the column
            colnames(Sum_range_table)[1]<-"Sum.range"
        }

        #Product of ranges
        if(any(method == 'product.range')) {
            #Product of range
            Product_range_table<-Disparity.measure.table(type_function=prod.apply, ranges, central_tendency, CI)
            #Renaming the column
            colnames(Product_range_table)[1]<-"Prod.range"           
        }
        if(verbose==TRUE) {
            message("Done.", appendLF=TRUE)
        }
    }

    #VARIANCE
    if(any(grep("variance", method))) {
        #Calculate the variance for the rarefaction and the bootstrapped matrices
        if(verbose==TRUE) {
            message("Calculating variance...", appendLF=FALSE)
        }
        variances<-lapply(BSresult, variance.calc)

        #Sum of variance
        if(any(method == 'sum.variance')) {
            #Sum of variance
            Sum_variance_table<-Disparity.measure.table(type_function=sum.apply, variances, central_tendency, CI)
            #Renaming the column
            colnames(Sum_variance_table)[1]<-"Sum.var"   
        }

        #Product of variance
        if(any(method == 'product.variance')) {
            #Product of variance
            Product_variance_table<-Disparity.measure.table(type_function=prod.apply, variances, central_tendency, CI)
            #Renaming the column
            colnames(Product_variance_table)[1]<-"Prod.var"            
        }
        if(verbose==TRUE) {
            message("Done.", appendLF=TRUE)
        }
    }

    #OUTPUT
    #Empty output table
    if(rarefaction==FALSE) {
        output<-matrix(nrow=1, data=rep(NA, 1))
        colnames(output)[1]<-"rarefaction"
    } else {
        output<-matrix(nrow=(ncol(data)-1), data=seq(from=2, to=ncol(data)))
        colnames(output)[1]<-"rarefaction"
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

    return(output)

#End
}