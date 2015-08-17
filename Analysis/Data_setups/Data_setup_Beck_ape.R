#Script for doing ancestral states reconstruction and distance matrix from an input nexus file

#Loading the package
library(devtools)
install_github("TGuillerme/SpatioTemporal_Disparity/Functions/disparity")
library(disparity)


###################
#Reading the files
###################

#Selecting the file
chain_name<-"Beck2014" #"Beck2014" #"Slater2013"
data_path<-"../Data/"
file_matrix<-"../Data/2014-Beck-ProcB-matrix-morpho.nex" #2014-Beck-ProcB-matrix-morpho.nex #"../Data/2013-Slater-MEE-matrix-morpho.nex"
file_tree<-"../Data/2014-Beck-ProcB-TEM.tre" #2014-Beck-ProcB-TEM.tre #"../Data/2013-Slater-MEE-TEM.tre"

#matrix
Nexus_data<-ReadMorphNexus(file_matrix)
Nexus_matrix<-Nexus_data$matrix
#tree
Tree_data<-read.nexus(file_tree)

######################################
#Cleaning the matrices and the trees
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

#Adding node labels to the tree
tree$node.label<-paste("n",seq(1:Nnode(tree)), sep="")


#Set folder for saving
system(paste("mkdir ",data_path, chain_name, sep=""))

####################################
#Ancestral states reconstruction
####################################

anc_states<-anc.state(tree, Nexus_data, method='ML-ape', verbose=TRUE)
save(anc_states, file=paste("../Data/",chain_name,"ape_ancestral_states.Rda", sep=""))

#anc_states<-anc.state(tree, Nexus_data, method='ML-claddis', verbose=TRUE)
#save(anc_states, file=paste(data_path, chain_name, "/",chain_name,"Claddis_ancestral_states.Rda", sep=""))

####################################
#Distance matrix
####################################

#Distance matrix using tips only
matrix_tips<-Nexus_data
message("\nCalculating the distance matrix for the tips only...", appendLF=FALSE)
dist_tips<-MorphDistMatrix.verbose(matrix_tips, verbose=TRUE)
message("Done.\n", appendLF=FALSE)
save(dist_tips, file=paste(data_path, chain_name, "/",chain_name,"ape_distance-tips.Rda", sep="")) #dist_tips

#Distance matrix using also nodes
matrix_nodes<-Nexus_data
matrix_nodes$matrix<-anc_states$state
message("\nCalculating the distance matrix for the tips and the nodes...", appendLF=FALSE)
dist_nodes<-MorphDistMatrix.verbose(matrix_nodes, verbose=TRUE)
message("Done.\n", appendLF=FALSE)
save(dist_nodes, file=paste(data_path, chain_name, "/",chain_name,"ape_distance-nodes.Rda", sep="")) #dist_nodes

#Distance matrix using nodes with 95 CI
matrix_nodes95<-Nexus_data
matrix_nodes95$matrix<-anc.unc(anc_states, 0.95, missing=NA)$state
message("\nCalculating the distance matrix for the tips and the nodes with a 95 CI...", appendLF=FALSE)
dist_nodes95<-MorphDistMatrix.verbose(matrix_nodes95, verbose=TRUE)
message("Done.\n", appendLF=FALSE)
save(dist_nodes95, file=paste(data_path, chain_name, "/",chain_name,"ape_distance-nodes95.Rda", sep="")) #dist_nodes95


