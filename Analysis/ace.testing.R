#Script for testing the difference between my anc.state function and Graemes AncStateEstMatrix function

#Setwd
if(length(grep("TGuillerme", getwd()))) {
    setwd('~/PhD/Projects/SpatioTemporal_Disparity/Analysis')
} else {
    warning("You might have to change the directory!")
}
if(length(grep("SpatioTemporal_Disparity/Analysis", getwd()))==0) {
    if(length(grep("SpatioTemporal_Disparity-master/Analysis", getwd()))==0) {
        stop("Wrong directory!\nThe current directory must be:\nSpatioTemporal_Disparity/Analysis/ OR SpatioTemporal_Disparity-master/Analysis/\nYou can clone the whole repository from:\nhttps://github.com/TGuillerme/SpatioTemporal_Disparity")
    }
}

#Load the functions and the packages
source("functions.R")

#Dummy data
set.seed(0)
table<-matrix(data=sample(c(0,1), 100, replace=TRUE), nrow=10)
set.seed(0)
tree<-rtree(10)
rownames(table)<-tree$tip.label

#Creating the nexus data table (as in ReadMorphNexus output)
table.nex<-list()
#header
table.nex$header<-"test"
#matrix
table.nex$matrix<-table
#ordering
table.nex$ordering<-rep("unord", 10) #unordered
#weighting
table.nex$weights<-rep(1, 10) #none
#max values
table.nex$max.vals<-apply(table, 2, max)
#min values
table.nex$min.vals<-apply(table, 2, min)
#step.matrices
table.nex$step.matrices<-NULL
#symbols
table.nex$symbols<-c("0","1")


#anc.state (ML-ape)
set.seed(0)
time_ape<-system.time(test_ape<-anc.state(tree, table.nex, method='ML-ape', verbose=TRUE))

#anc.state (ML-claddis)
set.seed(0)
time_claddis<-system.time(test_claddis<-anc.state(tree, table.nex, method='ML-claddis', verbose=TRUE))

#Time difference
cat(paste("ML-ape is", round((time_claddis/time_ape)[3], digit=2), "times faster than ML-claddis!"))

#Results difference
#Results for tips must be equal!
all(test_ape[[1]][c(1:Ntip(tree)),] == test_claddis[[1]][c(1:Ntip(tree)),])
all(test_ape[[2]][c(1:Ntip(tree)),] == test_claddis[[2]][c(1:Ntip(tree)),])
#Results for nodes?
diff_matrix<-test_ape[[1]][-c(1:Ntip(tree)),] == test_claddis[[1]][-c(1:Ntip(tree)),]
#Number of difference
diff_cells<-length(which(diff_matrix == FALSE))
#Which characters are different
diff_char<-which((apply(diff_matrix, 2, all)) == FALSE)

#Check the probabilities (rounded to 4 digit)
diff_prob<-round(test_ape[[2]][-c(1:Ntip(tree)),], digit=4) == round(test_claddis[[2]][-c(1:Ntip(tree)),], digit=4)

#Results are too ambiguous when the Max.lik is nearly equal to 0.5. Use the anc.unc function to set a threshold of uncertainty
test_ape_unc<-anc.unc(test_ape, 0.95)
test_claddis_unc<-anc.unc(test_claddis, 0.95)

#Now the results are equal !
all(test_ape_unc[[1]][-c(1:Ntip(tree)),] == test_claddis_unc[[1]][-c(1:Ntip(tree)),])
#But we have a lot of missing data... But at least we're being conservative


#Random NAs (25%)
set.seed(0)
data=sample(c(0,1), 100, replace=TRUE)
set.seed(0)
data[sample(c(1:100), 25)]<-NA
table<-matrix(data=data, nrow=10)
rownames(table)<-tree$tip.label
table.nex$matrix<-table

#Testing with NAs
#test_ape<-anc.state(tree, table.nex, method='ML-ape', verbose=TRUE)
#test_claddis<-anc.state(tree, table.nex, method='ML-claddis', verbose=TRUE)


#Running the simulations
try(result_20t_50c_010na<-testing.ace(20, 50, 0.10, 100) , silent=TRUE)
try(save(result_20t_50c_010na, file="../Data/ace_test/result_20t_50c_010na.Rda"), silent=TRUE)

try(result_20t_50c_025na<-testing.ace(20, 50, 0.25, 100) , silent=TRUE)
try(save(result_20t_50c_025na, file="../Data/ace_test/result_20t_50c_025na.Rda"), silent=TRUE)

try(result_20t_50c_050na<-testing.ace(20, 50, 0.50, 100) , silent=TRUE)
try(save(result_20t_50c_050na, file="../Data/ace_test/result_20t_50c_050na.Rda"), silent=TRUE)

try(result_20t_100c_010na<-testing.ace(20, 100, 0.10, 100) , silent=TRUE)
try(save(result_20t_100c_010na, file="../Data/ace_test/result_20t_100c_010na.Rda"), silent=TRUE)

try(result_20t_100c_025na<-testing.ace(20, 100, 0.25, 100) , silent=TRUE)
try(save(result_20t_100c_025na, file="../Data/ace_test/result_20t_100c_025na.Rda"), silent=TRUE)

try(result_20t_100c_050na<-testing.ace(20, 100, 0.50, 100) , silent=TRUE)
try(save(result_20t_100c_050na, file="../Data/ace_test/result_20t_100c_050na.Rda"), silent=TRUE)

try(result_20t_200c_010na<-testing.ace(20, 200, 0.10, 100) , silent=TRUE)
try(save(result_20t_200c_010na, file="../Data/ace_test/result_20t_200c_010na.Rda"), silent=TRUE)

try(result_20t_200c_025na<-testing.ace(20, 200, 0.25, 100) , silent=TRUE)
try(save(result_20t_200c_025na, file="../Data/ace_test/result_20t_200c_025na.Rda"), silent=TRUE)

try(result_20t_200c_050na<-testing.ace(20, 200, 0.50, 100) , silent=TRUE)
try(save(result_20t_200c_050na, file="../Data/ace_test/result_20t_200c_050na.Rda"), silent=TRUE)

try(result_50t_100c_010na<-testing.ace(50, 100, 0.10, 100) , silent=TRUE)
try(save(result_50t_100c_010na, file="../Data/ace_test/result_50t_100c_010na.Rda"), silent=TRUE)

try(result_50t_100c_025na<-testing.ace(50, 100, 0.25, 100) , silent=TRUE)
try(save(result_50t_100c_025na, file="../Data/ace_test/result_50t_100c_025na.Rda"), silent=TRUE)

try(result_50t_100c_050na<-testing.ace(50, 100, 0.50, 100) , silent=TRUE)
try(save(result_50t_100c_050na, file="../Data/ace_test/result_50t_100c_050na.Rda"), silent=TRUE)

try(result_50t_200c_010na<-testing.ace(50, 200, 0.10, 100) , silent=TRUE)
try(save(result_50t_200c_010na, file="../Data/ace_test/result_50t_200c_010na.Rda"), silent=TRUE)

try(result_50t_200c_025na<-testing.ace(50, 200, 0.25, 100) , silent=TRUE)
try(save(result_50t_200c_025na, file="../Data/ace_test/result_50t_200c_025na.Rda"), silent=TRUE)

try(result_50t_200c_050na<-testing.ace(20, 200, 0.50, 100) , silent=TRUE)
try(save(result_50t_200c_050na, file="../Data/ace_test/result_50t_200c_050na.Rda"), silent=TRUE)

try(result_100t_200c_010na<-testing.ace(100, 200, 0.10, 100) , silent=TRUE)
try(save(result_100t_200c_010na, file="../Data/ace_test/result_100t_200c_010na.Rda"), silent=TRUE)

try(result_100t_200c_025na<-testing.ace(100, 200, 0.25, 100) , silent=TRUE)
try(save(result_100t_200c_025na, file="../Data/ace_test/result_100t_200c_025na.Rda"), silent=TRUE)

try(result_100t_200c_050na<-testing.ace(100, 200, 0.50, 100) , silent=TRUE)
try(save(result_100t_200c_050na, file="../Data/ace_test/result_100t_200c_050na.Rda"), silent=TRUE)

#Plotting the results
load("../Data/ace_test/result_100t_200c_050na.Rda")
results_list<-result_100t_200c_050na.Rda
boxplot(results_list$ape$correct, results_list$ape$error, results_list$ape$na,
    results_list$ape95$correct, results_list$ape95$error, results_list$ape95$na,
    results_list$claddis$correct, results_list$claddis$error, results_list$claddis$na,
    results_list$claddis95$correct, results_list$claddis95$error, results_list$claddis95$na,
    names=c("ape-co", "ape-er", "ape-na", "a95-co", "a95-er", "a95-na", "cla-co", "cla-er", "cla-na", "c95-co", "c95-er", "c95-na"),
    las=2, col=rep(c("green", "red", "grey"), 4))