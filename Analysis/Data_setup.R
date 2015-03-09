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
source("functions.R")

###################
#Reading the files
###################

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

#WARNING: problem with NAs and ML-ape

#if(all(Nexus_data$ordering == "unord")) {
    #If no ordered characters
#    anc_states<-anc.state(tree, Nexus_data, method='ML-ape', verbose=TRUE)
#    save(anc_states, file=paste("../Data/",chain_name,"_ancestral_states-ape.Rda", sep=""))
    #load(paste("../Data/",chain_name,"_ancestral_states-ape.Rda", sep="")) #anc_states
#} else {
    anc_states<-anc.state(tree, Nexus_data, method='ML-claddis', verbose=TRUE)
    save(anc_states, file=paste(data_path, chain_name, "/",chain_name,"_ancestral_states-claddis.Rda", sep=""))
    #load(paste("../Data/",chain_name,"_ancestral_states-claddis.Rda", sep="")) #anc_states
#}

####################################
#Distance matrix
####################################

#Distance matrix using tips only
matrix_tips<-Nexus_data
message("Calculating the distance matrix for the tips only...", appendLF=FALSE)
dist_tips<-MorphDistMatrix(matrix_tips)
message("Done.\n", appendLF=FALSE)
save(dist_tips, file=paste(data_path, chain_name, "/",chain_name,"_distance-tips.Rda", sep="")) #dist_tips

#Distance matrix using also nodes
matrix_nodes<-Nexus_data
matrix_nodes$matrix<-anc_states$states
message("Calculating the distance matrix for the tips and the nodes...", appendLF=FALSE)
dist_nodes<-MorphDistMatrix(matrix_nodes)
message("Done.\n", appendLF=FALSE)
save(dist_nodes, file=paste(data_path, chain_name, "/",chain_name,"_distance-nodes.Rda", sep="")) #dist_nodes



#dist.data <- MorphDistMatrix(nexus.data)
