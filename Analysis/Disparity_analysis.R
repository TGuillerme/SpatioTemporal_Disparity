#Script for doing ancestral states reconstruction and distance matrix from an input nexus file

library(devtools)
install_github("TGuillerme/SpatioTemporal_Disparity/Functions/disparity")
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
name<-load(paste(data_path, chain_name, "/", chain_name, "_disparity-sli_prodist_nodes95-gower.dist.matrix.Rda", sep=""))
dis_tmp<-get(name)
#Isolating the disparity
disparity_full<-dis_tmp$disparity
#Isolating the diversity
diversity_full<-dis_tmp$diversity

#Isolated the non-rarefaction disparity results
dis_ran_max<-extract.disp(disparity_full$quantiles, rarefaction="max")
dis_ran_min<-extract.disp(disparity_full$quantiles, rarefaction="min")
dis_ran_med<-extract.disp(disparity_full$quantiles, rarefaction=10)

op<-par(mfrow=c(3,1), bty="l")
plot.disparity(dis_ran_max, diversity=log(dis_ran_max$rarefaction), main="max")
plot.disparity(dis_ran_min, diversity=log(dis_ran_min$rarefaction), main="min")
plot.disparity(dis_ran_med, diversity=log(dis_ran_med$rarefaction), main="10")
par(op)

#Pretty plot with the tree and the disparity
op<-par(mfrow=c(2,1), bty="l")
geoscalePhylo(ladderize(tree), cex.age=0.6, cex.ts=0.8, cex.tip=0.5)
plot.disparity(dis_ran_max, diversity=log(dis_ran_max$rarefaction))
par(op)