##########################
#Disparity functions
##########################
#Calculate the disparity as the distance from centroid
#This function is based on DisparityCalc() from Smith et al. 2014 - Evolution (http://dx.doi.org/10.1111/evo.12435) http://datadryad.org/resource/doi:10.5061/dryad.d380g 
#v0.1
##########################
#SYNTAX :
#<distance> the distance matrix
#<method> the method for calculating the disparity. Can be any of the following: "centroid", "sum.range", "product.range", "sum.variance", "product.variance"
#<CI> confidence intervals (default=c(50,95))
#<bootstraps> the number of boostrap replicates (default=1000)
#<central_tendency> any function for calculating the central tendency
##########################
#----
#guillert(at)tcd.ie 02/03/2015
##########################

#DEBUG
distance<-dist.data$max.dist.matrix

disparity<-function(data, method=c("centroid", "sum.range", "product.range", "sum.variance", "product.variance"), CI=c(50, 95), bootstraps=1000, central_tendency=median, rarefaction=FALSE) {

    #SANITIZING
    #distance
    check.class(distance, "matrix", " must be a distance matrix.")

    #method
    check.class(method, "character", " must be 'centroid', 'sum.range', 'product.range', 'sum.variance' or/and 'product.variance'.")
    methods_list<-c("centroid", "sum.range", "product.range", "sum.variance", "product.variance")
    if(all(is.na(match(method, methods_list)))) {
        stop("method must be 'centroid', 'sum.range', 'product.range', 'sum.variance' or/and 'product.variance'.")
    }

    #CI
    check.class(CI, "numeric", " must be any value between 1 and 100.")
    if(any(CI) < 1) {
        stop("CI must be any value between 1 and 100.")
    }
    if(any(CI) > 100) {
        stop("CI must be any value between 1 and 100.")
    }

    #Bootstrap
    check.class(bootstraps, "numeric", " must be a single (entire) numerical value.")
    check.length(bootstraps, 1, " must be a single (entire) numerical value.")
    #Make sure the bootstrap is a whole number
    bootstraps<-round(abs(bootstraps))

    #Central tendency
    check.class(central_tendency, "function", " must be either a function (e.g. 'mean' or 'median'.")

    #rarefaction
    check.class(rarefaction, "logical", " must be logical.")

    #FUNCTIONS

    #Performs bootstrap and eventual rarefaction
    Bootstrap.rarefaction<-function(data, iterations, rarefaction) {

        #Set rarefaction (or not)
        if(rarefaction == TRUE) {
            rarefaction_max<-seq(1:ncol(data))
        } else {
            rarefaction_max<-ncol(data)
        }

        #Rarefaction
        for(rare in rarefaction_max){
            #Bootstraps
            for(BS in 1:iterations){ #iterations -> bootstraps
                #Bootstrap
                output<-as.matrix(data[sample(1:nrow(data),rare,TRUE),])
                result[BS] <- list(output)
            }
            #Rarefaction + BS results
            BSresult[rare]<-list(result)
        }

        #Remove first element if rarefaction
        if(rarefaction == TRUE) {
            BSresult[[1]]=NULL
        } else {
            #Removing the n-1 first elements
            BSresult<-BSresult[-c(1:(rarefaction_max-1))]
        }

        return(BSresult)
    }

    #Range Calculations
    range.calc<-function(list_table) {
        #Empty matrix (for output)
        output<-matrix(nrow=length(list_table), ncol=ncol(list_table[[1]]))
        #Looping through columns and rows
        for(row in 1:length(list_table)) { #Rows are bootstraps
            for(column in 1:ncol(list_table[[1]])) { #Columns are axes
                output[row,column]<-max(list_table[[row]][,column])-min(list_table[[row]][,column])
            }
        }
        return(output)
    }

    #Variance calculation
    variance.calc<-function(list_table) {
        #Empty matrix (for output)
        output<-matrix(nrow=length(list_table), ncol=ncol(list_table[[1]]))
        #Looping through columns and rows
        for(row in 1:length(list_table)) { #Rows are bootstraps
            for(column in 1:ncol(list_table[[1]])) { #Columns are axes
                output[row,column]<-var(list_table[[row]][,column])
            }
        }
        return(output)
    }

    #Centroid Apply
    centroid.apply<-function(X) {
        #FUNCTION FROM SIVE, ADD DOI
        #Centroid (mean score of each PC axis)
        centroid<-apply(X, 2, mean)
        #Euclidean distances to the centroid
        cent.dist<-NULL
        for (j in 1:nrow(X)){
            cent.dist[j] <- dist(rbind(X[j,], centroid), method="euclidean")
        }
        return(cent.dist)
    }

    #Set-up for the NthRoot function in order to scale your products correctly.
    NthRoot<-function(x=data, n=ncol(data)){
       x^(1/n)
    }

    #Apply loop for calculating the product
    prod.apply<-function(X) {
        output<-NthRoot(apply(X,1,prod))
        return(output)
    }

    #Apply loop for calculating the sum
    sum.apply<-function(X) {
        output<-apply(X,1,sum)
        return(output)
    }

    #No apply (does nothin)
    no.apply<-function(X) {
        return(X)
    }

    #Apply loop for calculating the centroid
    centroid.apply<-function(X) {
        #FUNCTION FROM SIVE, ADD DOI
        #Centroid (mean score of each PC axis)
        centroid<-apply(X, 2, mean)
        #Euclidean distances to the centroid
        cent.dist<-NULL
        for (j in 1:nrow(X)){
            cent.dist[j] <- dist(rbind(X[j,], centroid), method="euclidean")
        }
        return(cent.dist)
    }

    #Lapply loop for calculating the centroid
    centroid.calc<-function(X) {
        Y<-lapply(X, centroid.apply)
        return(matrix(nrow=10, data=unlist(Y), byrow=TRUE))
    }

    #Converts one or more CI into a quantile probabilities
    CI.converter<-function(CI) {
        sort(c(50-CI/2, 50+CI/2)/100)
    }

    #Calculate product for the central tendency and the CIs for variance or range
    Disparity.measure.table<-function(type_function, distribution_variable, central_tendency, CI) {

        #Products/Sum of distribution_variable (correct by the NthRoot)
        Disparity_measure<-lapply(distribution_variable, type_function)
        #Confidence intervals for the products/sum of distribution_variable
        Disparity_measure_quantile<-lapply(Disparity_measure, quantile, probs=CI.converter(CI))
        #Calculate the central tendency
        Disparity_measure_central<-lapply(Disparity_measure, central_tendency)

        #Transform the results into a data.frame
        #Add just the first column (central tendency)
        Disparity_measure_table<-data.frame("Central"=unlist(Disparity_measure_central))
        #Add the CIs
        Disparity_measure_table<-cbind(Disparity_measure_table, matrix(ncol=length(CI)*2, data=unlist(Disparity_measure_quantile), byrow=TRUE))
        #Add the CIs names
        colnames(Disparity_measure_table)[-1]<-paste(CI.converter(CI)*100, "%", sep="")

        return(Disparity_measure_table)
    }

    #CALCULATING THE DISPARITY
    #Bootstraping the matrix
    BSresults<-Bootstrap.rarefaction(data, iterations, rarefaction)

    #Calculate the variance for the rarefaction and the bootstrapped matrices
    variances<-lapply(BSresult, variance.calc)

    #CENTROID
    #Distance form centroid
    if(any(method == 'centroid')) {
        #Calculate the distance from centroid for the rarefaction and the bootstrapped matrices
        centroids<-lapply(BSresult, centroid.calc)
        #Distance to centroid
        Centroid_dist_table<-Disparity.measure.table(type_function=no.apply, centroids, central_tendency, CI)
        #Renaming the column
        colnames(Centroid_dist_table)[1]<-"Cent.dist"
    }

    #RANGES
    if(any(grep("range", method))) {
        #Calculate the range for the rarefaction and the bootstrapped matrices
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
    }

    #VARIANCE
    if(any(grep("variance", method))) {
        #Calculate the variance for the rarefaction and the bootstrapped matrices
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