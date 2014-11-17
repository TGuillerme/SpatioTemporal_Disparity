#Data input (test)

#This is a pruned version of Beck & Lee tree and matrix for testing the STD pipeline analysis (based on Claddis package).

#Reading in the non modified data
source("Beck.data.R")

#Selecting only the euarchontoglires (node 168)
euarch.tree<-extract.clade(tree, 168)
euarch.tree$root.time<-max(tree.age(euarch.tree)[,1])

#Pruning the other taxa from the table
euarch.table<-Beck.table[match(euarch.tree$tip.label, rownames(Beck.table)),]

#Removing constant characters
euarch.table.levels<-apply(euarch.table, 2, function(X) length(levels(as.factor(X))))
euarch.table<-euarch.table[,-c(which(euarch.table.levels < 2))]

#Randomly selecting only 50 characters
set.seed(1)
euarch.table<-euarch.table[,c(sample(1:ncol(euarch.table), 50))]

#Creating the nexus data table (as in ReadMorphNexus output)
euarch.nex<-Beck.nex
#matrix
euarch.nex$matrix<-euarch.table
#ordering
euarch.nex$ordering<-rep("unord", 50)
#weighting
euarch.nex$weights<-rep(1, 50)
#max valsues
max.vals<-apply(euarch.table, 2, max, na.rm=TRUE)
try(if(grep("&", max.vals)) {
    max.vals[grep("&", max.vals)]<-strsplit(max.vals[grep("&", max.vals)], split="&")[[1]][2]
}, silent=TRUE)
euarch.nex$max.vals<-as.numeric(max.vals)
#min valsues
min.vals<-apply(euarch.table, 2, min, na.rm=TRUE)
try(if(grep("&", min.vals)) {
    min.vals[grep("&", min.vals)]<-strsplit(min.vals[grep("&", min.vals)], split="&")[[1]][2]
}, silent=TRUE)
euarch.nex$min.vals<-as.numeric(min.vals)

#Ancestral states
#Using STD anc.state function
#Replacing NAs by "?" in the matrix
matrix<-ifelse(is.na(euarch.table), '?', euarch.table)
STD.ace<-anc.state(euarch.tree, matrix, method='ML', verbose=TRUE)
STD.ace<-anc.unc(STD.ace, 0.5) #WARNING MINIMUM CONFIDENCE
STD.ace.table<-ifelse(STD.ace$state == '?', NA, STD.ace$state)

#Remove node labels otherwise reroot{phytools} freaks out
tmp.tree<-euarch.tree
tmp.tree$node.label<-NULL
CLADDIS.ace<-AncStateEstMatrix(euarch.nex, tmp.tree, estimate.allchars=TRUE, estimate.tips=FALSE)

#Comment on AncStateEstMatrix:
#Doesn't allow any uncertainty in the ancestral states reconstruction
#Maybe update AncStateEstMatrix allowing more options + saving likelihood + verbose

#STD matrix
STD.nex<-euarch.nex
STD.nex$matrix<-STD.ace.table

#Claddis matrix
CLADDIS.nex<-euarch.nex
CLADDIS.nex$matrix<-rbind(CLADDIS.nex$matrix, CLADDIS.ace)
#Changing the tips names
rownames(CLADDIS.nex$matrix)[(Ntip(euarch.tree)+1):nrow(CLADDIS.nex$matrix)]<-euarch.tree$node.label



#Reading some not so good data (diet and body mass)
table<-read.table("../Data/TestingData.csv", row.names=1, header=TRUE, sep=",")
#logBM
table$log_bm<-log(table$mass_g)
#adding mean logBM for NAs
table$log_bm[which(is.na(table$log_bm))]<-mean(table$log_bm, na.rm=TRUE) #STOOPID
#Estimate ancestral body masses and ancestral diet
#Diet
anc.diet<-ace(as.factor(table[,5]), euarch.tree, type="d")$lik.anc
rownames(anc.diet)<-euarch.tree$node.label
diets<-NULL
for (node in 1:nrow(anc.diet)) {
    diets[[node]]<-names(which((anc.diet[node,]) == max(anc.diet[node,])))
}

#BM
anc.bm<-ace(as.numeric(table[,6]), euarch.tree, type="c")$ace
names(anc.bm)<-euarch.tree$node.label

table.anc<-table[-1,]
rownames(table.anc)<-euarch.tree$node.label
table.anc[,c(1:4)]<-NA
table.anc$Diet<-diets
table.anc$log_bm<-anc.bm

table<-rbind(table, table.anc)
