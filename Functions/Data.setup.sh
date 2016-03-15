#!/bin/sh
##########################
#Shell script for setting up the data for disparity analysis.
##########################
#SYNTAX:
#sh Data.setup.sh <chain> <path> <matrix> <tree> <type>
#with:
#<chain> the name of the chain to generate task files for
#<path> the path where the data will be stored under the chain name
#<matrix> the path to the morphological matrix (can be relative)
#<tree> the path to the phylogenetic tree (can be relative)
#<type> either "ape" or "Claddis" to calculate the distance using one or the other method
#########################
#version 0.1
#----
#guillert(at)tcd.ie - 08/04/2015
###########################

#INPUT
chain=$1
path=$2
matrix=$3
tree=$4
ace=$5

#Set the saving path folder
mkdir ${path}/${chain}

#R CODE TEMPLATE
echo "
#Load the functions and the packages
library(disparity)

###################
#Reading the files
###################

#Selecting the file
chain_name='${chain}'
data_path='${path}'
file_matrix='${matrix}'
file_tree='${tree}'

#matrix
Nexus_data<-ReadMorphNexus(file_matrix)
Nexus_matrix<-Nexus_data\$matrix
#tree
Tree_data<-read.nexus(file_tree)

######################################
#Cleaning the matrices and the trees
######################################

#Remove species with only missing data before hand
if (any(apply(as.matrix(Nexus_matrix), 1, function(x) levels(as.factor((x)))) == '?')) {
    Nexus_matrix<-Nexus_matrix[-c(as.vector(which(apply(as.matrix(Nexus_matrix), 1, function(x) levels(as.factor(x))) == '?'))),]
}

#Cleaning the tree and the table
#making the saving folder
tree<-clean.tree(Tree_data, Nexus_matrix)
table<-clean.table(Nexus_matrix, Tree_data)
Nexus_data\$matrix<-table

#Forcing the tree to be binary
tree<-bin.tree(tree)

#Adding node labels to the tree
tree\$node.label<-paste('n',seq(1:Nnode(tree)), sep='')
" > ${chain}-data_setup.R

#ADD THE OPTIONS (WHAT TO CALCULATE)

#Ancestral states
if echo $ace | grep 'ape' > /dev/null
then
    #Add the APE ace
    echo "
    ####################################
    #Ancestral states reconstruction
    ####################################
    anc_states<-anc.state(tree, Nexus_data, method='ML-ape', verbose=TRUE)
    save(anc_states, file=paste('../Data/',chain_name,'_ancestral_states-ape.Rda', sep=''))
    " >> ${chain}-data_setup.R

else
    #Add the Claddis ace
    echo "
    ####################################
    #Ancestral states reconstruction
    ####################################
    anc_states<-anc.state(tree, Nexus_data, method='ML-claddis', verbose=TRUE)
    save(anc_states, file=paste(data_path, chain_name, '/',chain_name,'_ancestral_states-claddis.Rda', sep=''))
    " >> ${chain}-data_setup.R
fi

#Add the distance script
echo "
####################################
#Distance matrix
####################################
#Distance with tips only
matrix_tips<-Nexus_data
message('\\nCalculating the distance matrix for the tips only...', appendLF=FALSE)
dist_tips<-MorphDistMatrix.verbose(matrix_tips, verbose=TRUE)
message('Done.\\n', appendLF=FALSE)
save(dist_tips, file=paste(data_path, chain_name, '/',chain_name,'_distance-tips.Rda', sep='')) #dist_tips

#Distance matrix using also nodes
matrix_nodes<-Nexus_data
matrix_nodes\$matrix<-anc_states\$state
message('\\nCalculating the distance matrix for the tips and the nodes...', appendLF=FALSE)
dist_nodes<-MorphDistMatrix.verbose(matrix_nodes, verbose=TRUE)
message('Done.\\n', appendLF=FALSE)
save(dist_nodes, file=paste(data_path, chain_name, '/',chain_name,'_distance-nodes.Rda', sep='')) #dist_nodes

#Distance matrix using also nodes (safe)
matrix_nodes95<-Nexus_data
matrix_nodes95\$matrix<-anc.unc(anc_states, 0.95, missing=NA)\$state
message('\\nCalculating the distance matrix for the tips and the nodes with a 95 CI...', appendLF=FALSE)
dist_nodes95<-MorphDistMatrix.verbose(matrix_nodes95, verbose=TRUE)
message('Done.\\n', appendLF=FALSE)
save(dist_nodes95, file=paste(data_path, chain_name, '/',chain_name,'_distance-nodes95.Rda', sep='')) #dist_nodes95
" >> ${chain}-data_setup.R

#Remove template
rm ${chain}-data_setup.R