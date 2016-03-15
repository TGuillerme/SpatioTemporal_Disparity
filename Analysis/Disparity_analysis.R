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

#Setting the variables
#Constant
data_path<-"../Data/"
disaprity_pro<-"-disp_sli_nodes95_pro.Rda"
disparity_ran<-"-disp_sli_nodes95_ran.Rda"
diversity_ful<-"-slices_nodes95_div.Rda"
distance_gowr<-"_distance-nodes95.Rda"

#Data
chain_name<-c("Slater2013", "Beck2014")
file_matrix<-c("../Data/2013-Slater-MEE-matrix-morpho.nex", "../Data/2014-Beck-ProcB-matrix-morpho.nex")
file_tree<-c("../Data/2013-Slater-MEE-TEM.tre","../Data/2014-Beck-ProcB-TEM.tre")

#Loading the data
#Disparity + tree
slat_tmp1<-read.data(chain_name[1], data_path, file_matrix[1], file_tree[1], disaprity_pro)
slat_tmp2<-read.data(chain_name[1], data_path, file_matrix[1], file_tree[1], disparity_ran)
beck_tmp1<-read.data(chain_name[2], data_path, file_matrix[2], file_tree[2], disaprity_pro)
beck_tmp2<-read.data(chain_name[2], data_path, file_matrix[2], file_tree[2], disparity_ran)
#Diversity
slat_div<-load(paste(data_path, chain_name[1], "/", chain_name[1], diversity_ful, sep=""))
beck_div<-load(paste(data_path, chain_name[2], "/", chain_name[2], diversity_ful, sep=""))

#Extracting the data
#Tree
tree_slater<-slat_tmp1[[2]]
tree_beck  <-beck_tmp1[[2]]
#Disparity
disparity_full_pro_slater<-slat_tmp1[[3]]
disparity_full_ran_slater<-slat_tmp2[[3]]
disparity_full_pro_beck  <-beck_tmp1[[3]]
disparity_full_ran_beck  <-beck_tmp2[[3]]

#Diversity
diversity_full_slater<-get(slat_div)
diversity_full_beck  <-get(beck_div)

######################################
# Tree visualisation
######################################

#Plot the tree (slater)
geoscalePhylo(ladderize(tree_slater), cex.age=0.6, cex.ts=0.7, cex.tip=0.5, units=c("Period","Epoch"), boxes="Epoch")
abline(v=240, col="red")
dev.new()
geoscalePhylo(ladderize(tree_beck), cex.age=0.6, cex.ts=0.7, cex.tip=0.5, units=c("Period","Epoch"), boxes="Epoch")
abline(v=106, col="red")
dev.new()

######################################
# Disparity visualisation
######################################
dis_pro_max_beck<-extract.disp(disparity_full_ran_beck$quantiles, rarefaction="max")
dis_pro_max_slater<-extract.disp(disparity_full_ran_slater$quantiles, rarefaction="max")
dis_ran_max_beck<-extract.disp(disparity_full_pro_beck$quantiles, rarefaction="max")
dis_ran_max_slater<-extract.disp(disparity_full_pro_slater$quantiles, rarefaction="max")

#Pretty plot with the tree and the disparity
op<-par(mfrow=c(2,2), bty="n", mar=c(4,4,4,4))
plot.disparity(dis_ran_max_beck, diversity=dis_ran_max_beck$rarefaction, main="Eutheria (punctuated)", xlab="Time (Ma)", y2lab="", ylab="Tips and nodes median distance from centroid")
abline(v=22, col="red")
plot.disparity(dis_ran_max_slater, diversity=dis_ran_max_slater$rarefaction, main="Mammaliaformes (punctuated)", xlab="Time (Ma)", y2lab="Number of tips and nodes", ylab="")
abline(v=22, col="red")
plot.disparity(dis_pro_max_beck, diversity=dis_pro_max_beck$rarefaction, main="Eutheria (gradual)", xlab="Time (Ma)", y2lab="", ylab="Tips and nodes median distance from centroid")
abline(v=22, col="red")
plot.disparity(dis_pro_max_slater, diversity=dis_pro_max_slater$rarefaction, main="Mammaliaformes (gradual)", xlab="Time (Ma)", y2lab="Number of tips and nodes", ylab="")
abline(v=22, col="red")
par(op)

######################################
# Rarefaction
######################################
dis_pro_max_beck<-extract.disp(disparity_full_ran_beck$quantiles, rarefaction=8)
dis_pro_max_slater<-extract.disp(disparity_full_ran_slater$quantiles, rarefaction="min")
dis_ran_max_beck<-extract.disp(disparity_full_pro_beck$quantiles, rarefaction=8)
dis_ran_max_slater<-extract.disp(disparity_full_pro_slater$quantiles, rarefaction="min")

#Pretty plot with the tree and the disparity
op<-par(mfrow=c(2,2), bty="n", mar=c(4,4,4,4))
plot.disparity(dis_ran_max_beck, diversity=dis_ran_max_beck$rarefaction, main="Eutheria (punctuated)", xlab="Time (Ma)", y2lab="", ylab="Tips and nodes median distance from centroid")
abline(v=22, col="red")
plot.disparity(dis_ran_max_slater, diversity=dis_ran_max_slater$rarefaction, main="Mammaliaformes (punctuated)", xlab="Time (Ma)", y2lab="Number of tips and nodes", ylab="")
abline(v=22, col="red")
plot.disparity(dis_pro_max_beck, diversity=dis_pro_max_beck$rarefaction, main="Eutheria (gradual)", xlab="Time (Ma)", y2lab="", ylab="Tips and nodes median distance from centroid")
abline(v=22, col="red")
plot.disparity(dis_pro_max_slater, diversity=dis_pro_max_slater$rarefaction, main="Mammaliaformes (gradual)", xlab="Time (Ma)", y2lab="Number of tips and nodes", ylab="")
abline(v=22, col="red")
par(op)
