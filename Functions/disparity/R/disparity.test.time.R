##########################
#disparity.test.time
##########################
#function for testing the effect of time on data
#v0.1
##########################
#SYNTAX :
#<data> can be either a list of matrices (e.g. for multivariate pco through time using permanova) or a time.disparity output (for univariate analysis, e.g. aov).
#<force.test> optional, whether to force the test to be parametric or non parametric (can be "parametric", "non-parametric", or "NULL" (default))
#<method> optional, distance method for the permanova (see adonis). If NULL method="euclidean".
#<permutations> optional, number of permutations for the permanova (see adonis). If NULL permutations=1000.
##########################
#----
#guillert(at)tcd.ie 09/07/2015
##########################

disparity.test.time<-function(data, force.test=NULL, method=NULL, permutations=NULL, ...) {

    #SANITIZING
    message("Warning: no sanitizing in this version yet!\nUse this function at your own risks...")
    #Add data sanitizing input
    if(class(data[[1]]) == "matrix") {
        data_type<-"Multivariate"
    } else {
        data_type<-"Univariate"
    }
    
    #type
    type<-"Time"

    #force.test
    if(!is.null(force.test) == TRUE) {
        check.class(force.test, "character", " can be either 'parametric', 'non-parametric' or 'NULL'.")
        check.length(force.test, 1, " can be either 'parametric', 'non-parametric' or 'NULL'.")
        if(force.test != "parametric") {
            if(force.test != "non-parametric") {
                stop("force.test must be either 'parametric', 'non-parametric' or 'NULL'.")
            }
        }     
    } else {
        force.test<-FALSE
    }

    #setting the test type to null (default)
    test_type=NULL

    #method
    if(is.null(method)) {
        method="euclidean"
    } else {
        check.class(method, "character")
        test_type<-"Multivariate"
        message("'method' options is not null. The test is assumed to be multivariate.")
    }

    #method
    if(is.null(permutations)) {
        permutations=1000
    } else {
        check.class(permutations, "numeric")
        test_type<-"Multivariate"
        message("'permutations' options is not null. The test is assumed to be multivariate.")
    }

    #TEST DIFFERENCES IN DISPARITY

    #TIME: testing the effect of time on data
    if(type == "Time") {

        if(data_type == "Univariate") {
            #Univariate analysis (Kruskall Wallis or ANOVA)
            univ_data<-make.univariate.table(data)
            data.values<-univ_data[,1]
            time<-univ_data[,2]
        
            #Testing assumptions of anova
            test<-test.anova(univ_data)

            #Running the univariate test
            if(test[[1]] == TRUE) {
                if(force.test == "non-parametric") {
                    #KW test
                    univ_test<-kruskal.test(data.values~time)
                    message("Tested input values ~ input time using Kruskal-Wallis rank sum test.")
                } else {
                    #ANOVA
                    univ_test<-summary(aov(data.values~time))
                    message("Tested input values ~ input time using ANOVA.")
                    message("Residuals are normally distributed and variance within groups is homogeneous.")
                    message("Is the data independent? If not use the 'force.test=\"non-parametric\"' options.")
                }
            } else {
                if(force.test == "parametric") {
                    #ANOVA
                    univ_test<-summary(aov(data.values~time))
                    message("Tested input values ~ input time using ANOVA (parametric test was enforced using 'force.test=\"parametric\"' options).")
                    message("Residuals are NOT normally distributed and variance within groups is NOT homogeneous.")                
                } else {
                    #KW test
                    univ_test<-kruskal.test(data.values~time)
                    message("Tested input values ~ input time using Kruskal-Wallis rank sum test.")
                }

            }

            #Do test output
        }

        if(is.null(test_type) & data_type == "Multivariate") {
            test_type<-"Multivariate"
        } else {
            test_type<-"Univariate"
        }

        if(test_type == "Multivariate") {
            #Is the data univariate?
            if(data_type == "Univariate") {
                multi_data<-make.univariate.table(data)
                data.values<-univ_data[,1]
                time<-univ_data[,2]
                message("Multivariate test is performed on univariate data. Remove the multivariate options to enforce univariate test.")
            } else {
                multi_data<-make.multivariate.table(data)
                data.values<-multi_data[[1]]
                time<-multi_data[[2]]
            }

            #Running PERMANOVA
            multi_test<-adonis(data.values~time, method=method, permutations=permutations)
            message(paste("Tested input values ~ input time using PERMANOVA with ", method," distance and ", permutations, " permutations.", sep=""))

        }

        #Output
        if(test_type == "Multivariate") {
            return(multi_test)
        } else {
            if(test[[1]] == TRUE) {
                return(list("Disparity test"=univ_test, "Parametric test"=test[[2]]))
            } else {
                return(univ_test)
            }
        }

    }
}


