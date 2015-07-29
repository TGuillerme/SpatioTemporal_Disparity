
#Load the functions and the packages
library(disparity)

###################
#Reading the files
###################

#Selecting the file
chain_name='Slater2013'
data_path='../Data/'
file_matrix='../Data/2013-Slater-MEE-matrix-morpho.nex'
file_tree='../Data/2013-Slater-MEE-TEM.tre'
slices=rev(seq(from=0, to=295, by=5))
FADLAD='../Data/Slater2013_FADLAD.csv'

#matrix
Nexus_data<-ReadMorphNexus(file_matrix)
Nexus_matrix<-Nexus_data$matrix
#tree
Tree_data<-read.nexus(file_tree)

#FAD/LAD
FADLAD<-read.csv(FADLAD, row.names=1)

######################################
#Cleaning the matrices and the trees
######################################

#Remove species with only missing data before hand
if(any(apply(is.na(Nexus_matrix), 1, all))) {
    Nexus_matrix<-Nexus_matrix[-c(which(apply(is.na(Nexus_matrix), 1, all))),]
}

#Cleaning the tree and the table
#making the saving folder
tree<-clean.tree(Tree_data, Nexus_matrix)
table<-clean.table(Nexus_matrix, Tree_data)
Nexus_data$matrix<-table

#Forcing the tree to be binary
tree<-bin.tree(tree)

#Adding node labels to the tree
tree<-lapply.root(tree, max(tree.age(tree)$age)) 

#load the distance matrix
load(paste(data_path, chain_name, '/', chain_name, '_distance-nodes95.Rda', sep='')) #dist_nodes95
trimmed_max_data_nodes95<-TrimMorphDistMatrix(dist_nodes95$gower.dist.matrix)
tree_nodes95<-drop.tip(tree, trimmed_max_data_nodes95$removed.taxa) ; tree_nodes95$root.time<-max(tree.age(tree_nodes95)[,1])
trimmed_max_data_nodes95$dist.matrix<-trimmed_max_data_nodes95$dist.matrix[c(tree_nodes95$tip.label, tree_nodes95$node.label),c(tree_nodes95$tip.label, tree_nodes95$node.label)]
#pco
pco_data_nodes95<-cmdscale(trimmed_max_data_nodes95$dist.matrix, k=nrow(trimmed_max_data_nodes95$dist.matrix) - 2, add=T)$points


#slices
#pco_slices_nodes95_ran<-slice.pco(pco_data_nodes95, tree_nodes95, slices, method='random', FAD_LAD=FADLAD, verbose=TRUE, diversity=TRUE)
#pco_slices_nodes95_pro<-slice.pco(pco_data_nodes95, tree_nodes95, slices, method='proximity', FAD_LAD=FADLAD, verbose=TRUE)
#slices_nodes95_div<-pco_slices_nodes95_ran[[2]] ; pco_slices_nodes95_ran<-pco_slices_nodes95_ran[[1]]
#save(slices_nodes95_div,file=paste(data_path, chain_name, '/',chain_name,'-Full-slices_nodes95_div.Rda', sep=''))
load(paste(data_path, chain_name, '/',chain_name,'-Full-slices_nodes95_div.Rda', sep=''))


#Disparity
#disp_sli_nodes95_ran<-time.disparity(pco_slices_nodes95_ran, method='centroid', verbose=TRUE, rarefaction=TRUE, save.all=TRUE)
#save(disp_sli_nodes95_ran, file=paste(data_path, chain_name, '/',chain_name,'-Full-disp_sli_nodes95_ran.Rda', sep=''))
load(paste(data_path, chain_name, '/',chain_name,'-Full-disp_sli_nodes95_ran.Rda', sep=''))

#disp_sli_nodes95_pro<-time.disparity(pco_slices_nodes95_pro, method='centroid', verbose=TRUE, rarefaction=TRUE, save.all=TRUE)
#save(disp_sli_nodes95_pro, file=paste(data_path, chain_name, '/',chain_name,'-Full-disp_sli_nodes95_pro.Rda', sep=''))
load(paste(data_path, chain_name, '/',chain_name,'-Full-disp_sli_nodes95_pro.Rda', sep=''))

######################################
#Plotting the disparity
######################################
#load the relevant files generated before!

#corrected rar div
rar_diversity<-extract.disp(disp_sli_nodes95_pro[[1]], rarefaction=8)$rarefaction

op<-par(mfrow=c(3,2), bty="n", mar=c(4,4,4,4))
plot.disparity(extract.disp(disp_sli_nodes95_ran[[1]], rarefaction="max"), diversity=slices_nodes95_div, main="Mammaliaformes (punctuated)", xlab="Time (Mya)", y2lab="", ylab="Median distance from centroid")
abline(v=47, col="red")
plot.disparity(extract.disp(disp_sli_nodes95_pro[[1]], rarefaction="max"), diversity=slices_nodes95_div, main="Mammaliaformes (gradual)", xlab="Time (Mya)", y2lab="Species richness", ylab="")
abline(v=47, col="red")
plot.disparity(extract.disp(disp_sli_nodes95_ran[[1]], rarefaction=8), diversity=rar_diversity, main="Mammaliaformes (punctuated - 8 taxa)", xlab="Time (Mya)", y2lab="", ylab="Median distance from centroid")
abline(v=47, col="red")
plot.disparity(extract.disp(disp_sli_nodes95_pro[[1]], rarefaction=8), diversity=rar_diversity, main="Mammaliaformes (gradual - 8 taxa)", xlab="Time (Mya)", y2lab="Species richness", ylab="")
abline(v=47, col="red")
plot.disparity(extract.disp(disp_sli_nodes95_ran[[1]], rarefaction="min"), diversity=rep(3,60), main="Mammaliaformes (punctuated - 3 taxa)", xlab="Time (Mya)", y2lab="", ylab="Median distance from centroid")
abline(v=47, col="red")
plot.disparity(extract.disp(disp_sli_nodes95_pro[[1]], rarefaction="min"), diversity=rep(3,60), main="Mammaliaformes (gradual - 3 taxa)", xlab="Time (Mya)", y2lab="Species richness", ylab="")
abline(v=47, col="red")
par(op)
