
#Load the functions and the packages
library(disparity)

###################
#Reading the files
###################

#Selecting the file
chain_name='Beck2014'
data_path='../../Data/'
file_matrix='../../Data/2014-Beck-ProcB-matrix-morpho.nex'
file_tree='../../Data/2014-Beck-ProcB-TEM.tre'
intervals=as.numeric(strsplit(c(noquote('170.300,168.300,166.100,163.500,157.300,152.100,145.000,139.800,132.900,129.400,125.000,113.000,100.500,93.900,89.800,86.300,83.600,72.100,66.000,61.600,59.200,56.000,47.800,41.300,38.000,33.900,28.100,23.030,23.030,20.440,15.970,13.820,11.620,7.246,5.333,0.000')), split=',')[[1]])
slices=as.numeric(strsplit(c(noquote('170,165,160,155,150,145,140,135,130,125,120,115,110,105,100,95,90,85,80,75,70,65,60,55,50,45,40,35,30,25,20,15,10,5,0')), split=',')[[1]])
FADLAD='../../Data/Beck2014_FADLAD.csv'

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
trimmed_max_data_nodes95<-TrimMorphDistMatrix(dist_nodes95$max.dist.matrix)
tree_nodes95<-drop.tip(tree, trimmed_max_data_nodes95$removed.taxa) ; tree_nodes95$root.time<-max(tree.age(tree_nodes95)[,1])
trimmed_max_data_nodes95$dist.matrix<-trimmed_max_data_nodes95$dist.matrix[c(tree_nodes95$tip.label, tree_nodes95$node.label),c(tree_nodes95$tip.label, tree_nodes95$node.label)]
#pco
pco_data_nodes95<-cmdscale(trimmed_max_data_nodes95$dist.matrix, k=nrow(trimmed_max_data_nodes95$dist.matrix) - 2, add=T)$points


#slices
pco_slices_nodes95_acc<-slice.pco(pco_data_nodes95, tree_nodes95, slices, method='acctran', FAD_LAD=FADLAD, verbose=TRUE, diversity=TRUE)
slices_nodes95_div<-pco_slices_nodes95_acc[[2]] ; pco_slices_nodes95_acc<-pco_slices_nodes95_acc[[1]]


#Disparity
disp_sli_nodes95_acc<-time.disparity(pco_slices_nodes95_acc, verbose=TRUE, rarefaction=TRUE, save.all=TRUE)
save(disp_sli_nodes95_acc, file=paste(data_path, chain_name, '/',chain_name,'-disp_sli_nodes95_acc.Rda', sep=''))

