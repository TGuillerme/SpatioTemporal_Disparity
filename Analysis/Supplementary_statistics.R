#Loading the package
#library(devtools)
#install_github("TGuillerme/SpatioTemporal_Disparity/Functions/disparity")
#library(disparity)

###################################
#
# TESTING THE DIFFERENCES THROUGH TIME
#
###################################

######################################
# functions from disparity.test_fun.R for significance tokens
######################################

#P-value significance levels
apply.signif<-function(p.value) {
    if(p.value < 0.001) {
        p<-"***"
    } else {
        if(p.value < 0.01) {
            p<-"**"
        } else {
            if(p.value < 0.05) {
                p<-"*"
            } else {
                if(p.value < 0.1) {
                    p<-"."
                } else {
                    if(p.value > 0.1) {
                        p<-" "
                    }
                }
            }
        }
    }
    return(p)
}

#Getting the significance tokens
signif.token<-function(p.values) {
    output<-vector()
    for(val in 1:length(p.values)) {
        output<-c(output, apply.signif(p.values[val]))
    }
    return(output)
}

#Combining tests function
combine.tests<-function(methods_names, method_numbers, test_file, test_metric, data_name) {
    n_rows<-0
    n_methods<-vector()
    row_names<-vector()
    for(meth in 1:length(methods_names)) {
        n_row_tmp<-nrow(test_file[[method_numbers[meth]]][[1]])
        n_rows<-c(n_rows+n_row_tmp)
        n_methods<-c(n_methods, rep(methods_names[meth], n_row_tmp))
        row_names<-c(row_names, rownames(test_file[[method_numbers[meth]]][[1]]))
    }
    n_methods<-as.factor(n_methods)

    #Generating the csv file for each metric
    csv_out_tmp<-as.data.frame(matrix(NA, ncol=9, nrow=n_rows))
    colnames(csv_out_tmp)<-c("Data", "Method", "metric", "slice", "difference", "Df", "T", "p-value"," ")
    csv_out_tmp$Data<-rep(data_name, n_rows)
    csv_out_tmp$Method<-n_methods
    csv_out_tmp$metric<-rep(test_metric, n_rows)
    csv_out_tmp$slice<-row_names

    for(meth in 1:length(levels(n_methods))) {
        meth_length<-which(n_methods == levels(n_methods)[meth])
        for(test in 1:length(meth_length)) {
            csv_out_tmp[meth_length[test], 5:9]<-test_file[[method_numbers[meth]]][[1]][test,]
        }
    }

    return(csv_out_tmp)
}


######################################
# permanova results
######################################

#Loading the tests
permanovas_slater<-get(load("../Data/Statistics/permanova_slater.Rda"))
permanovas_beck  <-get(load("../Data/Statistics/permanova_beck.Rda"))

#Generating the csv file
csv_out<-as.data.frame(matrix(NA, ncol=10, nrow=24))
colnames(csv_out)<-c("Data", "Method", "terms", "Df", "SumsOfSqs", "MeanSqs", "F.Model", "R2", "Pr(>F)"," ")
csv_out$Data<-c(rep("Mammaliaformes", 12), rep("Eutherians", 12))
csv_out$Method<-c(rep(c(rep("Intervals:tips",2), rep("Intervals:tips+nodes",2),rep("Slices:punctuated",2), rep("Slices:ACCTRAN",2), rep("Slices:DELTRAN",2), rep("Slices:gradual",2)), 2))
csv_out$terms<-c(rep(c("time", "residuals"), 12))

for(test in 1:length(permanovas_slater)) {
    #slater
    csv_out[(test+(test-1)),4:9]<-as.vector(as.data.frame(permanovas_slater[[test]][[1]][1,]))
    csv_out[(test+(test-1)),10]<-signif.token(csv_out[(test+(test-1)),9])
    csv_out[(test+(test-1)+1),4:9]<-as.vector(as.data.frame(permanovas_slater[[test]][[1]][2,]))
    csv_out[(test+(test-1)+1),10]<-" "
    #beck
    csv_out[(test+(test-1)+12),4:9]<-as.vector(as.data.frame(permanovas_beck[[test]][[1]][1,]))
    csv_out[(test+(test-1)+12),10]<-signif.token(csv_out[(test+(test-1)+12),9])
    csv_out[(test+(test-1)+1+12),4:9]<-as.vector(as.data.frame(permanovas_beck[[test]][[1]][2,]))
    csv_out[(test+(test-1)+1+12),10]<-" "
}

write.csv(csv_out, "../Data/Permanova_results.csv")

######################################
# reference results
######################################

#loading the tests
#slater
reftests_centroid_slater<-get(load("../Data/Statistics/reftests_centroid_slater.Rda"))
reftests_sum.range_slater<-get(load("../Data/Statistics/reftests_sum.range_slater.Rda"))
reftests_product.range_slater<-get(load("../Data/Statistics/reftests_product.range_slater.Rda"))
reftests_sum.variance_slater<-get(load("../Data/Statistics/reftests_sum.variance_slater.Rda"))
reftests_product.variance_slater<-get(load("../Data/Statistics/reftests_product.variance_slater.Rda"))
#beck
reftests_centroid_beck<-get(load("../Data/Statistics/reftests_centroid_beck.Rda"))
reftests_sum.range_beck<-get(load("../Data/Statistics/reftests_sum.range_beck.Rda"))
reftests_product.range_beck<-get(load("../Data/Statistics/reftests_product.range_beck.Rda"))
reftests_sum.variance_beck<-get(load("../Data/Statistics/reftests_sum.variance_beck.Rda"))
reftests_product.variance_beck<-get(load("../Data/Statistics/reftests_product.variance_beck.Rda"))

#Generating the csv file
#slater
csv_out<-combine.tests(c("Intervals:tips", "Slices:DELTRAN", "Slices:gradual"), c(1,5,6), reftests_centroid_slater, "centroid", "Mammaliaformes")
csv_out<-rbind(csv_out, combine.tests(c("Intervals:tips", "Slices:DELTRAN", "Slices:gradual"), c(1,5,6), reftests_sum.range_slater, "sum.range", "Mammaliaformes"))
csv_out<-rbind(csv_out, combine.tests(c("Intervals:tips", "Slices:DELTRAN", "Slices:gradual"), c(1,5,6), reftests_product.range_slater, "product.range", "Mammaliaformes"))
csv_out<-rbind(csv_out, combine.tests(c("Intervals:tips", "Slices:DELTRAN", "Slices:gradual"), c(1,5,6), reftests_sum.variance_slater, "sum.variance", "Mammaliaformes"))
csv_out<-rbind(csv_out, combine.tests(c("Intervals:tips", "Slices:DELTRAN", "Slices:gradual"), c(1,5,6), reftests_product.variance_slater, "product.variance", "Mammaliaformes"))
#beck
csv_out<-rbind(csv_out, combine.tests(c("Intervals:tips", "Intervals:tips+nodes", "Slices:punctuated", "Slices:ACCTRAN", "Slices:DELTRAN", "Slices:gradual"), c(1,2,3,4,5,6), reftests_centroid_beck, "centroid", "Eutherians"))
csv_out<-rbind(csv_out, combine.tests(c("Intervals:tips", "Intervals:tips+nodes", "Slices:punctuated", "Slices:ACCTRAN", "Slices:DELTRAN", "Slices:gradual"), c(1,2,3,4,5,6), reftests_sum.range_beck, "sum.range", "Eutherians"))
csv_out<-rbind(csv_out, combine.tests(c("Intervals:tips", "Intervals:tips+nodes", "Slices:punctuated", "Slices:ACCTRAN", "Slices:DELTRAN", "Slices:gradual"), c(1,2,3,4,5,6), reftests_product.range_beck, "product.range", "Eutherians"))
csv_out<-rbind(csv_out, combine.tests(c("Intervals:tips", "Intervals:tips+nodes", "Slices:punctuated", "Slices:ACCTRAN", "Slices:DELTRAN", "Slices:gradual"), c(1,2,3,4,5,6), reftests_sum.variance_beck, "sum.variance", "Eutherians"))
csv_out<-rbind(csv_out, combine.tests(c("Intervals:tips", "Intervals:tips+nodes", "Slices:punctuated", "Slices:ACCTRAN", "Slices:DELTRAN", "Slices:gradual"), c(1,2,3,4,5,6), reftests_product.variance_beck, "product.variance", "Eutherians"))

write.csv(csv_out, "../Data/reftests_results.csv")

######################################
# sequential results
######################################

#loading the tests
#slater
seqtests_centroid_slater<-get(load("../Data/Statistics/seqtests_centroid_slater.Rda"))
seqtests_sum.range_slater<-get(load("../Data/Statistics/seqtests_sum.range_slater.Rda"))
seqtests_product.range_slater<-get(load("../Data/Statistics/seqtests_product.range_slater.Rda"))
seqtests_sum.variance_slater<-get(load("../Data/Statistics/seqtests_sum.variance_slater.Rda"))
seqtests_product.variance_slater<-get(load("../Data/Statistics/seqtests_product.variance_slater.Rda"))
#beck
seqtests_centroid_beck<-get(load("../Data/Statistics/seqtests_centroid_beck.Rda"))
seqtests_sum.range_beck<-get(load("../Data/Statistics/seqtests_sum.range_beck.Rda"))
seqtests_product.range_beck<-get(load("../Data/Statistics/seqtests_product.range_beck.Rda"))
seqtests_sum.variance_beck<-get(load("../Data/Statistics/seqtests_sum.variance_beck.Rda"))
seqtests_product.variance_beck<-get(load("../Data/Statistics/seqtests_product.variance_beck.Rda"))

#Generating the csv file
#slater
csv_out<-combine.tests(c("Intervals:tips", "Slices:DELTRAN", "Slices:gradual"), c(1,5,6), seqtests_centroid_slater, "centroid", "Mammaliaformes")
csv_out<-rbind(csv_out, combine.tests(c("Intervals:tips", "Slices:DELTRAN", "Slices:gradual"), c(1,5,6), seqtests_sum.range_slater, "sum.range", "Mammaliaformes"))
csv_out<-rbind(csv_out, combine.tests(c("Intervals:tips", "Slices:DELTRAN", "Slices:gradual"), c(1,5,6), seqtests_product.range_slater, "product.range", "Mammaliaformes"))
csv_out<-rbind(csv_out, combine.tests(c("Intervals:tips", "Slices:DELTRAN", "Slices:gradual"), c(1,5,6), seqtests_sum.variance_slater, "sum.variance", "Mammaliaformes"))
csv_out<-rbind(csv_out, combine.tests(c("Intervals:tips", "Slices:DELTRAN", "Slices:gradual"), c(1,5,6), seqtests_product.variance_slater, "product.variance", "Mammaliaformes"))
#beck
csv_out<-rbind(csv_out, combine.tests(c("Intervals:tips", "Intervals:tips+nodes", "Slices:punctuated", "Slices:ACCTRAN", "Slices:DELTRAN", "Slices:gradual"), c(1,2,3,4,5,6), seqtests_centroid_beck, "centroid", "Eutherians"))
csv_out<-rbind(csv_out, combine.tests(c("Intervals:tips", "Intervals:tips+nodes", "Slices:punctuated", "Slices:ACCTRAN", "Slices:DELTRAN", "Slices:gradual"), c(1,2,3,4,5,6), seqtests_sum.range_beck, "sum.range", "Eutherians"))
csv_out<-rbind(csv_out, combine.tests(c("Intervals:tips", "Intervals:tips+nodes", "Slices:punctuated", "Slices:ACCTRAN", "Slices:DELTRAN", "Slices:gradual"), c(1,2,3,4,5,6), seqtests_product.range_beck, "product.range", "Eutherians"))
csv_out<-rbind(csv_out, combine.tests(c("Intervals:tips", "Intervals:tips+nodes", "Slices:punctuated", "Slices:ACCTRAN", "Slices:DELTRAN", "Slices:gradual"), c(1,2,3,4,5,6), seqtests_sum.variance_beck, "sum.variance", "Eutherians"))
csv_out<-rbind(csv_out, combine.tests(c("Intervals:tips", "Intervals:tips+nodes", "Slices:punctuated", "Slices:ACCTRAN", "Slices:DELTRAN", "Slices:gradual"), c(1,2,3,4,5,6), seqtests_product.variance_beck, "product.variance", "Eutherians"))

write.csv(csv_out, "../Data/seqtests_results.csv")