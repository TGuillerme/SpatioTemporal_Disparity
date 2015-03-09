#Script for testing the difference between my anc.state function and Graemes AncStateEstMatrix function

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
max.vals<-apply(table, 2, max, na.rm=TRUE)
try(if(grep("&", max.vals)) {
    max.vals[grep("&", max.vals)]<-strsplit(max.vals[grep("&", max.vals)], split="&")[[1]][2]
}, silent=TRUE)
table.nex$max.vals<-as.numeric(max.vals)
#min values
min.vals<-apply(table, 2, min, na.rm=TRUE)
try(if(grep("&", min.vals)) {
    min.vals[grep("&", min.vals)]<-strsplit(min.vals[grep("&", min.vals)], split="&")[[1]][2]
}, silent=TRUE)
table.nex$min.vals<-as.numeric(min.vals)
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