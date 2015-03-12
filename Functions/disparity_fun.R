#FUNCTIONS FOR DISPARITY

#Performs bootstrap and eventual rarefaction
Bootstrap.rarefaction<-function(data, bootstraps, rarefaction) {

    #Set rarefaction (or not)
    if(rarefaction == TRUE) {
        rarefaction_max<-seq(1:nrow(data))
    } else {
        rarefaction_max<-nrow(data)
    }
    #Rarefaction
    result<-NULL
    BSresult<-NULL
    for(rare in rarefaction_max){
        #Bootstraps
        for(BS in 1:bootstraps){ #bootstraps -> bootstraps
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

#Set-up for the NthRoot function in order to scale your products correctly.
nth.root<-function(x, n){
    x^(1/n)
}

#Apply loop for calculating the product
prod.apply<-function(X) {
    output<-nth.root(apply(X,1,prod), n=ncol(X))
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