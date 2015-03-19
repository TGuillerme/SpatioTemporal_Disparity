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

#After running the Data_setup script, load the different results

######################
#Select the data
######################

#Selecting the file
#SLATER 2013 MEE
#chain_name<-"Slater2013"
#data_path<-"../Data/"
#file_matrix<-"../Data/2013-Slater-MEE-matrix-morpho.nex"
#file_tree<-"../Data/2013-Slater-MEE-TEM.tre"
#int_breaks<-rev(seq(from=0, to=250, by=32.5))
#slices<-rev(seq(from=0, to=250, by=20))+5
#slices[length(slices)]<-0
#KT_bin=5.5
#KT_sli=10

#BECK 2014 ProcB
chain_name<-"Beck2014"
data_path<-"../Data/"
file_matrix<-"../Data/2014-Beck-ProcB-matrix-morpho.nex"
file_tree<-"../Data/2014-Beck-ProcB-TEM.tre"
int_breaks<-rev(seq(from=0, to=150, by=20))+5
int_breaks[length(int_breaks)]<-0
slices<-rev(seq(from=0, to=150, by=10))
KT_bin=4.5
KT_sli=9.5

######################
#Tree and matrix
######################

#matrix
Nexus_data<-ReadMorphNexus(file_matrix)
Nexus_matrix<-Nexus_data$matrix
#tree
Tree_data<-read.nexus(file_tree)

#Cleaning the matrices and the trees
#Remove species with only missing data before hand
if (any(apply(as.matrix(Nexus_matrix), 1, function(x) levels(as.factor((x)))) == "?")) {
    Nexus_matrix<-Nexus_matrix[-c(as.vector(which(apply(as.matrix(Nexus_matrix), 1, function(x) levels(as.factor(x))) == "?"))),]
}

#Cleaning the tree and the table
#making the saving folder
tree<-clean.tree(Tree_data, Nexus_matrix)
table<-clean.table(Nexus_matrix, Tree_data)
Nexus_data$matrix<-table

#Forcing the tree to be binary
tree<-bin.tree(tree)

#Adding node labels to the tree
tree$node.label<-paste("n",seq(1:Nnode(tree)), sep="")
#Setting the tree root age
ages_data<-tree.age(tree)
tree$root.time<-max(ages_data[,1])

######################
#FADLAD file
######################

#Load the F/LAD for Beck
FADLAD<-read.csv(paste(data_path, chain_name, "_FADLAD.csv", sep=""), row.names=1)

######################
#Ancestral states reconstruction files
######################

load(paste(data_path, chain_name, "/", chain_name, "_ancestral_states-claddis.Rda", sep="")) #anc_states

######################
#Distance matrices
######################

load(paste(data_path, chain_name, "/", chain_name, "_distance-tips.Rda", sep="")) #dist_tips
load(paste(data_path, chain_name, "/", chain_name, "_distance-nodes.Rda", sep="")) #dist_nodes
load(paste(data_path, chain_name, "/", chain_name, "_distance-nodes95.Rda", sep="")) #dist_nodes95

#Remove the inapplicable characters
trimmed_max_data_tips<-TrimMorphDistMatrix(dist_tips$max.dist.matrix)
trimmed_max_data_nodes<-TrimMorphDistMatrix(dist_nodes$max.dist.matrix)
trimmed_max_data_nodes95<-TrimMorphDistMatrix(dist_nodes95$max.dist.matrix)
#Remove the dropped taxa from the tree
tree_tips<-drop.tip(tree, trimmed_max_data_tips$removed.taxa)
tree_nodes<-drop.tip(tree, trimmed_max_data_nodes$removed.taxa)
tree_nodes95<-drop.tip(tree, trimmed_max_data_nodes95$removed.taxa)
#Remove the eventual inapplicable nodes
trimmed_max_data_nodes$dist.matrix<-trimmed_max_data_nodes$dist.matrix[c(tree_nodes$tip.label, tree_nodes$node.label),c(tree_nodes$tip.label, tree_nodes$node.label)]
trimmed_max_data_nodes95$dist.matrix<-trimmed_max_data_nodes95$dist.matrix[c(tree_nodes$tip.label, tree_nodes$node.label),c(tree_nodes95$tip.label, tree_nodes95$node.label)]

#List of trees
trees<-list("tips"=tree_tips, "nodes"=tree_nodes, "nodes95"=tree_nodes95)

######################
#PCO
######################

pco_data_tips<-cmdscale(trimmed_max_data_tips$dist.matrix, k=nrow(trimmed_max_data_tips$dist.matrix) - 1, add=T)$points
pco_data_nodes<-cmdscale(trimmed_max_data_nodes$dist.matrix, k=nrow(trimmed_max_data_nodes$dist.matrix) - 1, add=T)$points
pco_data_nodes95<-cmdscale(trimmed_max_data_nodes95$dist.matrix, k=nrow(trimmed_max_data_nodes95$dist.matrix) - 1, add=T)$points

#Storing as a list
pco_data<-list("tips"=pco_data_tips, "nodes"=pco_data_nodes, "nodes95"=pco_data_nodes95)

######################
#Disparity
######################

#Calculating the rarefaction
#rarefaction_median<-disparity(pco_data, rarefaction=TRUE, verbose=TRUE, central_tendency=median)

#Generating the different intervals PCOs
pco_int_tips<-int.pco(pco_data_tips, tree_tips, int_breaks, include.nodes=FALSE, FAD_LAD=FADLAD, diversity=TRUE)
int_tips_div<-pco_int_tips[[2]] ; pco_int_tips<-pco_int_tips[[1]] 
pco_int_nodes<-int.pco(pco_data_nodes, tree_nodes, int_breaks, include.nodes=TRUE, FAD_LAD=FADLAD, diversity=TRUE)
int_nodes_div<-pco_int_nodes[[2]] ; pco_int_nodes<-pco_int_nodes[[1]]
pco_int_nodes95<-int.pco(pco_data_nodes95, tree_nodes95, int_breaks, include.nodes=TRUE, FAD_LAD=FADLAD, diversity=TRUE)
int_nodes95_div<-pco_int_nodes95[[2]] ; pco_int_nodes95<-pco_int_nodes95[[1]]

#Calculating the disparity per intervals
disp_int_tips<-time.disparity(pco_int_tips, verbose=TRUE)
#disp_int_tips_95axis<-time.disparity(pco_int_tips, verbose=TRUE, rm.last.axis=TRUE)
disp_int_nodes<-time.disparity(pco_int_nodes, verbose=TRUE)
#disp_int_nodes_95axis<-time.disparity(pco_int_nodes, verbose=TRUE, rm.last.axis=TRUE)
disp_int_nodes95<-time.disparity(pco_int_nodes95, verbose=TRUE)
#disp_int_nodes95_95axis<-time.disparity(pco_int_nodes95, verbose=TRUE, rm.last.axis=TRUE)

#Generating the different PCO slices
#methods list
methods=c("random", "acctran", "deltran", "proximity")
#nodes
pco_slices_nodes<-list()
slices_nodes_div<-list()
for (type in 1:length(methods)) {
    pco_slices_nodes[[type]]<-slice.pco(pco_data_nodes, tree_nodes, slices, method=methods[[type]], FAD_LAD=FADLAD, verbose=TRUE, diversity=TRUE)
    slices_nodes_div[[type]]<-pco_slices_nodes[[type]][[2]] ; pco_slices_nodes[[type]]<-pco_slices_nodes[[type]][[1]]
}
names(pco_slices_nodes)<-names(slices_nodes_div)<-methods

#nodes95
pco_slices_nodes95<-list()
slices_nodes95_div<-list()
for (type in 1:length(methods)) {
    pco_slices_nodes95[[type]]<-slice.pco(pco_data_nodes95, tree_nodes95, slices, method=methods[[type]], FAD_LAD=FADLAD, verbose=TRUE, diversity=TRUE)
    slices_nodes95_div[[type]]<-pco_slices_nodes95[[type]][[2]] ; pco_slices_nodes95[[type]]<-pco_slices_nodes95[[type]][[1]]
}
names(pco_slices_nodes95)<-names(slices_nodes95_div)<-paste(methods, "95", sep="")

#Calculating the disparity per interval per list
#nodes
disp_slices_nodes<-list()
for (type in 1:length(methods)) {
    disp_slices_nodes[[type]]<-time.disparity(pco_slices_nodes[[type]], verbose=TRUE)
}
names(disp_slices_nodes)<-names(pco_slices_nodes)

#nodes95
disp_slices_nodes95<-list()
for (type in 1:length(methods)) {
    disp_slices_nodes95[[type]]<-time.disparity(pco_slices_nodes95[[type]], verbose=TRUE)
}
names(disp_slices_nodes95)<-names(pco_slices_nodes95)

######################
#Plot the disparity
######################

#The following is a "BIG" plot of 
quartz(width = 22.4, height = 15.6) #A5 landscape
#Windows dimensions
op<-par(mfrow=c(5, 11), bty="l")# oma=c(bottom, left, top, right)
#Centroid
plot.disparity(disp_int_tips, rarefaction=FALSE, xlab="", ylab="Distance from centroid", measure="Cent.dist", main="Intervals: tips", diversity=int_tips_div)
abline(v= KT_bin, col="red")
plot.disparity(disp_int_nodes, rarefaction=FALSE, xlab="", ylab="", measure="Cent.dist", main="Intervals: nodes", diversity=int_nodes_div)
abline(v= KT_bin, col="red")
plot.disparity(disp_int_nodes95, rarefaction=FALSE, xlab="", ylab="", measure="Cent.dist", main="Intervals: nodes(95)", diversity=int_nodes95_div)
abline(v= KT_bin, col="red")
plot.disparity(disp_slices_nodes$random, rarefaction=FALSE, xlab="", ylab="", measure="Cent.dist", main="Slices: random", diversity=slices_nodes_div$random)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes$acctran, rarefaction=FALSE, xlab="", ylab="", measure="Cent.dist", main="Slices: acctran", diversity=slices_nodes_div$acctran)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes$deltran, rarefaction=FALSE, xlab="", ylab="", measure="Cent.dist", main="Slices: deltran", diversity=slices_nodes_div$deltran)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes$proximity, rarefaction=FALSE, xlab="", ylab="", measure="Cent.dist", main="Slices: proximity", diversity=slices_nodes_div$proximity)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes95$random, rarefaction=FALSE, xlab="", ylab="", measure="Cent.dist", main="Slices: random (95)", diversity=slices_nodes95_div$random)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes95$acctran, rarefaction=FALSE, xlab="", ylab="", measure="Cent.dist", main="Slices: acctran (95)", diversity=slices_nodes95_div$acctran)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes95$deltran, rarefaction=FALSE, xlab="", ylab="", measure="Cent.dist", main="Slices: deltran (95)", diversity=slices_nodes95_div$deltran)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes95$proximity, rarefaction=FALSE, xlab="", ylab="", measure="Cent.dist", main="Slices: proximity (95)", diversity=slices_nodes95_div$proximity)
abline(v= KT_sli, col="red")

#Sum of ranges
plot.disparity(disp_int_tips, rarefaction=FALSE, xlab="", ylab="Sum of ranges", measure="Sum.range", diversity=int_tips_div)
abline(v= KT_bin, col="red")
plot.disparity(disp_int_nodes, rarefaction=FALSE, xlab="", ylab="", measure="Sum.range", diversity=int_nodes_div)
abline(v= KT_bin, col="red")
plot.disparity(disp_int_nodes95, rarefaction=FALSE, xlab="", ylab="", measure="Sum.range", diversity=int_nodes95_div)
abline(v= KT_bin, col="red")
plot.disparity(disp_slices_nodes$random, rarefaction=FALSE, xlab="", ylab="", measure="Sum.range", diversity=slices_nodes_div$random)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes$acctran, rarefaction=FALSE, xlab="", ylab="", measure="Sum.range", diversity=slices_nodes_div$acctran)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes$deltran, rarefaction=FALSE, xlab="", ylab="", measure="Sum.range", diversity=slices_nodes_div$deltran)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes$proximity, rarefaction=FALSE, xlab="", ylab="", measure="Sum.range", diversity=slices_nodes_div$proximity)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes95$random, rarefaction=FALSE, xlab="", ylab="", measure="Sum.range", diversity=slices_nodes95_div$random)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes95$acctran, rarefaction=FALSE, xlab="", ylab="", measure="Sum.range", diversity=slices_nodes95_div$acctran)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes95$deltran, rarefaction=FALSE, xlab="", ylab="", measure="Sum.range", diversity=slices_nodes95_div$deltran)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes95$proximity, rarefaction=FALSE, xlab="", ylab="", measure="Sum.range", diversity=slices_nodes95_div$proximity)
abline(v= KT_sli, col="red")

#Sum of variance
plot.disparity(disp_int_tips, rarefaction=FALSE, xlab="", ylab="Sum of variance", measure="Sum.var", diversity=int_tips_div)
abline(v= KT_bin, col="red")
plot.disparity(disp_int_nodes, rarefaction=FALSE, xlab="", ylab="", measure="Sum.var", diversity=int_nodes_div)
abline(v= KT_bin, col="red")
plot.disparity(disp_int_nodes95, rarefaction=FALSE, xlab="", ylab="", measure="Sum.var", diversity=int_nodes95_div)
abline(v= KT_bin, col="red")
plot.disparity(disp_slices_nodes$random, rarefaction=FALSE, xlab="", ylab="", measure="Sum.var", diversity=slices_nodes_div$random)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes$acctran, rarefaction=FALSE, xlab="", ylab="", measure="Sum.var", diversity=slices_nodes_div$acctran)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes$deltran, rarefaction=FALSE, xlab="", ylab="", measure="Sum.var", diversity=slices_nodes_div$deltran)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes$proximity, rarefaction=FALSE, xlab="", ylab="", measure="Sum.var", diversity=slices_nodes_div$proximity)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes95$random, rarefaction=FALSE, xlab="", ylab="", measure="Sum.var", diversity=slices_nodes95_div$random)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes95$acctran, rarefaction=FALSE, xlab="", ylab="", measure="Sum.var", diversity=slices_nodes95_div$acctran)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes95$deltran, rarefaction=FALSE, xlab="", ylab="", measure="Sum.var", diversity=slices_nodes95_div$deltran)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes95$proximity, rarefaction=FALSE, xlab="", ylab="", measure="Sum.var", diversity=slices_nodes95_div$proximity)
abline(v= KT_sli, col="red")

#Product of ranges
plot.disparity(disp_int_tips, rarefaction=FALSE, xlab="", ylab="Product of range", measure="Prod.range", diversity=int_tips_div)
abline(v= KT_bin, col="red")
plot.disparity(disp_int_nodes, rarefaction=FALSE, xlab="", ylab="", measure="Prod.range", diversity=int_nodes_div)
abline(v= KT_bin, col="red")
plot.disparity(disp_int_nodes95, rarefaction=FALSE, xlab="", ylab="", measure="Prod.range", diversity=int_nodes95_div)
abline(v= KT_bin, col="red")
plot.disparity(disp_slices_nodes$random, rarefaction=FALSE, xlab="", ylab="", measure="Prod.range", diversity=slices_nodes_div$random)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes$acctran, rarefaction=FALSE, xlab="", ylab="", measure="Prod.range", diversity=slices_nodes_div$acctran)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes$deltran, rarefaction=FALSE, xlab="", ylab="", measure="Prod.range", diversity=slices_nodes_div$deltran)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes$proximity, rarefaction=FALSE, xlab="", ylab="", measure="Prod.range", diversity=slices_nodes_div$proximity)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes95$random, rarefaction=FALSE, xlab="", ylab="", measure="Prod.range", diversity=slices_nodes95_div$random)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes95$acctran, rarefaction=FALSE, xlab="", ylab="", measure="Prod.range", diversity=slices_nodes95_div$acctran)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes95$deltran, rarefaction=FALSE, xlab="", ylab="", measure="Prod.range", diversity=slices_nodes95_div$deltran)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes95$proximity, rarefaction=FALSE, xlab="", ylab="", measure="Prod.range", diversity=slices_nodes95_div$proximity)
abline(v= KT_sli, col="red")

#Product of variance
plot.disparity(disp_int_tips, rarefaction=FALSE, xlab="Mya", ylab="Product of variance", measure="Prod.var", diversity=int_tips_div)
abline(v= KT_bin, col="red")
plot.disparity(disp_int_nodes, rarefaction=FALSE, xlab="Mya", ylab="", measure="Prod.var", diversity=int_nodes_div)
abline(v= KT_bin, col="red")
plot.disparity(disp_int_nodes95, rarefaction=FALSE, xlab="Mya", ylab="", measure="Prod.var", diversity=int_nodes95_div)
abline(v= KT_bin, col="red")
plot.disparity(disp_slices_nodes$random, rarefaction=FALSE, xlab="Mya", ylab="", measure="Prod.var", diversity=slices_nodes_div$random)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes$acctran, rarefaction=FALSE, xlab="Mya", ylab="", measure="Prod.var", diversity=slices_nodes_div$acctran)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes$deltran, rarefaction=FALSE, xlab="Mya", ylab="", measure="Prod.var", diversity=slices_nodes_div$deltran)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes$proximity, rarefaction=FALSE, xlab="Mya", ylab="", measure="Prod.var", diversity=slices_nodes_div$proximity)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes95$random, rarefaction=FALSE, xlab="Mya", ylab="", measure="Prod.var", diversity=slices_nodes95_div$random)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes95$acctran, rarefaction=FALSE, xlab="Mya", ylab="", measure="Prod.var", diversity=slices_nodes95_div$acctran)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes95$deltran, rarefaction=FALSE, xlab="Mya", ylab="", measure="Prod.var", diversity=slices_nodes95_div$deltran)
abline(v= KT_sli, col="red")
plot.disparity(disp_slices_nodes95$proximity, rarefaction=FALSE, xlab="Mya", ylab="", measure="Prod.var", diversity=slices_nodes95_div$proximity)
abline(v= KT_sli, col="red")

par(op)
