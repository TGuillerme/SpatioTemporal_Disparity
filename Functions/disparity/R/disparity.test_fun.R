#FUNCTIONS FOR DISPARITY.TEST

#Making the univariate table (e.g. for aov)
make.univariate.table<-function(time_disparity) {
    #Extracting the list of values
    list_val<-time_disparity$values

    #Extracting the values
    val_uni<-unlist(list_val)
    #Removing the names
    names(val_uni)<-NULL

    #Extracting the time
    time_counts<-unlist(lapply(list_val, lapply, length))
    #Creating the time vector
    time<-as.factor(rep(names(list_val), time_counts))

    #Returning the results
    output<-cbind(val_uni, as.factor(time))
    output<-as.data.frame(output)
    output[,2]<-as.factor(time)
    names(output)<-c('Values', 'Time')
    return(output)
}
#next: summary(aov(bla$Values~bla$Time))


#Making the multivariate table (e.g. for adonis)
make.multivariate.table<-function(time_pco) {
    #Binding the first two elements
    pco_mat<-rbind(time_pco[[1]], time_pco[[2]])

    #Binding the other elements
    if(length(time_pco) > 2) {
        for (interval in 3:length(time_pco)) {
            pco_mat<-rbind(pco_mat, time_pco[[interval]])
        }
    }
    
    #Extracting the time
    time_counts<-unlist(lapply(time_pco, nrow))
    #Creating the time vector
    time<-as.factor(rep(names(time_pco), time_counts))

    #Returning the table
    output<-list(pco_mat, time)
    names(output)<-c("Data", "Time")
    return(output)
}
#next: adonis(bla~bla, method='euclidean', permutations=1000)

#Testing assumptions of anova
test.anova<-function(data) {
    #anova
    anova<-aov(data[,1]~data[,2])

    #residuals
    normal_res<-shapiro.test(rstandard(anova))

    #variance
    var_test<-bartlett.test(data[,1]~data[,2])

    #Test
    if(normal_res$p.value > 0.05 & var_test$p.value > 0.05) {
        test_results<-TRUE
    } else {
        test_results<-FALSE
    }

    #Saving the test table
    test.save<-data.frame(row.names=c("Shapiro","Bartlett"), statistic=c(normal_res$statistic[[1]], var_test$statistic[[1]]), df=c(NA, var_test$parameter[[1]]), p.value=c(normal_res$p.value[[1]], var_test$p.value[[1]]))

    #output
    output<-list("Pass"=test_results, "Test"=test.save)
    return(output)
}