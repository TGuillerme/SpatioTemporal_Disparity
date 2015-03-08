#Script for testing the difference between my anc.state function and Graemes AncStateEstMatrix function

#Dummy data
set.seed(1)
table<-matrix(data=sample(c(0,1), 100, replace=TRUE), nrow=10)
set.seed(1)
tree<-rtree(10)
rownames(table)<-tree$tip.label

#Creating the nexus data table (as in ReadMorphNexus output)
table.nex<-list()
#header
table.nex$header<-"test"
#matrix
table.nex$matrix<-table
#ordering
table.nex$ordering<-rep("unord", 10)
#weighting
table.nex$weights<-rep(1, 10)
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



#anc.state (ace{ape}) method
set.seed(1)
test_anc.state<-anc.state(tree, table, method='ML', verbose=TRUE)
nodes_ape<-test_anc.state$state[-c(1:10),]
#AncStateEstMatrix (rerootingMethods{phytools}) method
set.seed(1)
nodes_phy<-AncStateEstMatrix(table.nex, tree, estimate.allchars=TRUE, estimate.tips=FALSE)



