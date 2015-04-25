#Script for doing ancestral states reconstruction and distance matrix from an input nexus file

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
library(disparity)

###################################
#
# BECK ANALYSIS
#
###################################

#Selecting the file
chain_name<-"Beck2014"
data_path<-"../Data/"
file_matrix<-"../Data/2014-Beck-ProcB-matrix-morpho.nex"
file_tree<-"../Data/2014-Beck-ProcB-TEM.tre"

#matrix
Nexus_data<-ReadMorphNexus(file_matrix)
Nexus_matrix<-Nexus_data$matrix
#tree
Tree_data<-read.nexus(file_tree)

######################################
# Cleaning the matrices and the trees
######################################

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

#Adding node labels and the age to the tree
tree<-lapply.root(tree, max(tree.age(tree)$age))

######################################
# Tree visualisation
######################################

#Plot the tree
geoscalePhylo(ladderize(tree), cex.age=0.6, cex.ts=0.8, cex.tip=0.5)

######################################
# Disparity visualisation
######################################

#Load the disparity data
name<-load(paste(data_path, chain_name, "/", chain_name, "_disparity-intervalsdist_tips-gower.dist.matrix.Rda", sep=""))
dis_tmp<-get(name)
dis_ran<-dis_tmp$disparity
div_ran<-dis_tmp$diversity

#Isolated the non-rarefaction disparity results
dis_ran_max<-extract.disp(dis_ran, rarefaction="max")
dis_ran_min<-extract.disp(dis_ran, rarefaction="min")
dis_ran_med<-extract.disp(dis_ran, rarefaction=10)

op<-par(mfrow=c(3,1), bty="l")
plot.disparity(dis_ran_max, diversity=log(dis_ran_max$rarefaction), main="max")
plot.disparity(dis_ran_min, diversity=log(dis_ran_min$rarefaction), main="min")
plot.disparity(dis_ran_med, diversity=log(dis_ran_med$rarefaction), main="10")
par(op)

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
