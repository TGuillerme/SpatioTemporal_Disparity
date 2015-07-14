#Loading the package
library(devtools)
install_github("TGuillerme/SpatioTemporal_Disparity/Functions/disparity")
library(disparity)


#KT test
#KTdata_pro_slater<-int.pco(pco_slater, tree_slater , c(170,65,0), FAD_LAD=FADLADbeck)#
#KTtest_pro_slater<-disparity.test(KTdata_pro_slater, method="centroid", test="pairwise", bootstraps=1000)

#KTdata_pro_beck  <-int.pco(pco_beck  , tree_beck  , c(170,65,0), FAD_LAD=FADLADbeck)
#KTtest_pro_beck  <-disparity.test(KTdata_pro_beck, method="centroid", test="pairwise", bootstraps=1000)

######################################
# Making the xtable
######################################
library(xtable)

make.table.permanova<-function() {
    #Fixed rows
    #Empty matrix
    permanova_terms<-as.data.frame(matrix(" ", nrow=8, ncol=3))
    #Column names
    colnames(permanova_terms)<-c("Data", "model", "terms")#, c(colnames(permanova_pro_beck[[1]][1,])))
    #Data
    permanova_terms$Data<-c("Eutherian", rep(" ", 3), "Mammaliformes", rep(" ", 3))
    permanova_terms$model<-rep(c(c("gradual", rep(" ", 1)),c("punctuate", rep(" ", 1))),2)
    permanova_terms$terms<-rep(c("time", "Residuals"), 4)

    #Fixed rows
    #Empty matrix
    permanova_results<-as.data.frame(matrix(NA, nrow=8, ncol=6))
    #Column names
    colnames(permanova_results)<-c(colnames(permanova_pro_beck[[1]][1,]))

    #Variable rows
    permanova_results[1, ]<-as.vector(permanova_pro_beck[[1]][1,])
    permanova_results[2, ]<-as.vector(permanova_pro_beck[[1]][2,])
    permanova_results[3, ]<-as.vector(permanova_ran_beck[[1]][1,])
    permanova_results[4, ]<-as.vector(permanova_ran_beck[[1]][2,])
    permanova_results[5, ]<-as.vector(permanova_pro_slater[[1]][1,])
    permanova_results[6, ]<-as.vector(permanova_pro_slater[[1]][2,])
    permanova_results[7, ]<-as.vector(permanova_ran_slater[[1]][1,])
    permanova_results[8, ]<-as.vector(permanova_ran_slater[[1]][2,])

    #Rounding
    for (n in 1:6) {
        permanova_results[,n]<-round(permanova_results[,n], digit=n)
    }

    return(cbind(permanova_terms, permanova_results))
}

#permanova table
xtable(make.table.permanova())

#lag test table (to modify manually in LaTeX)
xtable(cbind(reftest_pro_beck[[1]][,-5], reftest_ran_beck[[1]][,-5]))
xtable(cbind(reftest_pro_slater[[1]][,-5], reftest_ran_slater[[1]][,-5]))
cat("add: & & gradual & & & & punctuated & & \\\\' in the header.")
cat("replace 'rrrrrrrrr' by 'rrrrr|rrrr'.")
