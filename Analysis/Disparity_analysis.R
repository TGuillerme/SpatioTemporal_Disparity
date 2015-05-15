#Loading the package
library(devtools)
install_github("TGuillerme/SpatioTemporal_Disparity/Functions/disparity")
library(disparity)

###################################
#
# MAIN FIGURE
#
###################################

######################################
# Isolating the data
######################################

#Selecting the slater data
chain_name<-"Slater2013"
data_path<-"../Data/"
file_matrix<-"../Data/2013-Slater-MEE-matrix-morpho.nex"
file_tree<-"../Data/2013-Slater-MEE-TEM.tre"
disparity_data<-"_disparity-sli_prodist_nodes95-gower.dist.matrix.Rda"
#Extracting all the data
tmp<-read.data(chain_name, data_path, file_matrix, file_tree, disparity_data)

#Isolating the slater data
disparity_full_slater<-tmp[[3]]$disparity
#Isolating the slater diversity
diversity_full_slater<-tmp[[3]]$diversity
#Isolating the slater tree
tree_slater<-tmp[[2]]

#Selecting the Beck data
chain_name<-"Beck2014"
data_path<-"../Data/"
file_matrix<-"../Data/2014-Beck-ProcB-matrix-morpho.nex"
file_tree<-"../Data/2014-Beck-ProcB-TEM.tre"
#Extracting all the data
tmp<-read.data(chain_name, data_path, file_matrix, file_tree, disparity_data)

#Isolating the slater data
disparity_full_beck<-tmp[[3]]$disparity
#Isolating the slater diversity
diversity_full_beck<-tmp[[3]]$diversity
#Isolating the slater tree
tree_beck<-tmp[[2]]

######################################
# Tree visualisation
######################################

#Plot the tree (slater)
geoscalePhylo(ladderize(tree_slater), cex.age=0.6, cex.ts=0.8, cex.tip=0.5)
dev.new()
geoscalePhylo(ladderize(tree_beck), cex.age=0.6, cex.ts=0.8, cex.tip=0.5)
dev.new()

######################################
# Disparity visualisation
######################################

dis_ran_max_beck<-extract.disp(disparity_full_beck$quantiles, rarefaction="max")
dis_ran_max_slater<-extract.disp(disparity_full_slater$quantiles, rarefaction="max")

#Pretty plot with the tree and the disparity
op<-par(mfrow=c(2,1), bty="n", mar=c(4,4,4,4))
plot.disparity(dis_ran_max_beck, diversity=log(dis_ran_max_beck$rarefaction), main="Eutherian", xlab="Time (Mya)", y2lab="Diversity (log)")
abline(v=22, col="red")
plot.disparity(dis_ran_max_slater, diversity=log(dis_ran_max_slater$rarefaction), main="Mammaliformes", xlab="Time (Mya)", y2lab="Diversity (log)")
abline(v=22, col="red")
par(op)
