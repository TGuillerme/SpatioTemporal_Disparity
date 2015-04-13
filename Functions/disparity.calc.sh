#!/bin/sh
##########################
#Shell script for Calculating the disparity from a distance matrix
##########################
#SYNTAX:
#sh disparity.calc.sh <chain> <path> <matrix> <tree> <file_dist> <distance> <intervals> <slices> <FADLAD>
#with:
#<chain> the name of the chain to generate task files for
#<path> the path where the data will be stored under the chain name
#<matrix> the path to the morphological matrix (can be relative)
#<file_dist> the path to a distance matrix (can be relative)
#<distance> which distance to use
#<tree> the path to the phylogenetic tree (can be relative)
#<intervals> a series of time intervals
#<slices> a series of time slices
#<FADLAD> a csv file with FADLAD data
#########################
#version 0.2
#Update: matching trimmed data with tree
#----
#guillert(at)tcd.ie - 13/04/2015
###########################

#INPUT
chain=$1
path=$2
matrix=$3
tree=$4
file_dist=$5
distance=$6
intervals=$7
slices=$8
FADLAD=$9

#Distance file prefix (for saving)
dist_pre=$(echo $file_dist | sed 's/.*\/'"${chain}"'//g' | sed 's/.Rda//')

######################################
#R template
######################################
echo "
#Load the functions and the packages
library(disparity)

######################
#Reading the files
######################

#Selecting the file
chain_name='${chain}'
data_path='${path}'
file_matrix='${matrix}'
file_tree='${tree}'
file_dist='${file_dist}'
distance='${distance}'
intervals=as.numeric(strsplit(c(noquote('${intervals}')), split=',')[[1]])
slices=as.numeric(strsplit(c(noquote('${slices}')), split=',')[[1]])
FADLAD='${FADLAD}'

#matrix
Nexus_data<-ReadMorphNexus(file_matrix)
Nexus_matrix<-Nexus_data\$matrix
#tree
Tree_data<-read.nexus(file_tree)

######################
#Cleaning the matrices and the trees
######################

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
tree<-lapply.root(tree, max(tree.age(tree)\$age)) 

#Load the FAD/LAD file
FADLAD<-read.csv(paste(data_path, chain_name, '_FADLAD.csv', sep=''), row.names=1)

######################
#Loading the distance matrix
######################
mat_name<-load(file_dist)
dist_mat<-extract.dist(get(mat_name), distance)

######################
#Trim the matrix and the tree (if necessary)
######################
trimmed_data<-TrimMorphDistMatrix(dist_mat)
tree<-drop.tip(tree, trimmed_data\$removed.taxa) ; tree\$root.time<-max(tree.age(tree)[,1])
#drop nodes from the distance matrix
if(length(trimmed_data\$removed.data) != 0) {
    trimmed_data\$dist.matrix<-trimmed_data\$dist.matrix[-which(is.na(match(rownames(trimmed_data\$dist.matrix), c(tree\$tip.label, tree\$node.label)))),-which(is.na(match(rownames(trimmed_data\$dist.matrix), c(tree\$tip.label, tree\$node.label))))]
}

######################
#pco
######################
pco_data<-cmdscale(trimmed_data\$dist.matrix, k=nrow(trimmed_data\$dist.matrix) - 1, add=T)\$points

######################
#Disparity
######################

#Intervals

#Test whether to include the nodes or not
if(nrow(trimmed_data\$dist.matrix) == (Ntip(tree) + Nnode(tree))) {
    incl_nodes<-TRUE
} else {
    incl_nodes<-FALSE
}

#Generating the different intervals PCOs
pco_int<-int.pco(pco_data, tree, intervals, include.nodes=incl_nodes, FAD_LAD=FADLAD, diversity=TRUE)
diversity_int<-pco_int[[2]] ; pco_int<-pco_int[[1]] 

#Calculating the disparity per intervals
disp_int<-time.disparity(pco_int, verbose=TRUE, rarefaction=TRUE)

#Saving
disparity_intervals<-list('disparity'=disp_int, 'diversity'=diversity_int)
save(disparity_intervals, file=paste(data_path, chain_name, '/',chain_name,'_disparity-intervals', mat_name,'-',distance ,'.Rda', sep=''))

#Slices (if incl_nodes == TRUE)
if(incl_nodes==TRUE) {
    #Generating the slices for each method
    pco_sli_ran<-slice.pco(pco_data, tree, slices, method='random', FAD_LAD=FADLAD, verbose=TRUE, diversity=TRUE)
    diversity_sli<-pco_sli_ran[[2]] ; pco_sli_ran<-pco_sli_ran[[1]]
    pco_sli_del<-slice.pco(pco_data, tree, slices, method='deltran', FAD_LAD=FADLAD, verbose=TRUE)
    pco_sli_acc<-slice.pco(pco_data, tree, slices, method='acctran', FAD_LAD=FADLAD, verbose=TRUE)
    pco_sli_pro<-slice.pco(pco_data, tree, slices, method='proximity', FAD_LAD=FADLAD, verbose=TRUE)

    #Calculating the disparity per interval
    disp_sli_ran<-time.disparity(pco_sli_ran, verbose=TRUE, rarefaction=TRUE)
    disparity_sli_ran<-list('disparity'=disp_sli_ran, 'diversity'=diversity_sli)
    save(disparity_sli_ran, file=paste(data_path, chain_name, '/',chain_name,'_disparity-sli_ran', mat_name,'-',distance ,'.Rda', sep=''))

    disp_sli_del<-time.disparity(pco_sli_del, verbose=TRUE, rarefaction=TRUE)
    disparity_sli_del<-list('disparity'=disp_sli_del, 'diversity'=diversity_sli)
    save(disparity_sli_del, file=paste(data_path, chain_name, '/',chain_name,'_disparity-sli_del', mat_name,'-',distance ,'.Rda', sep=''))

    disp_sli_acc<-time.disparity(pco_sli_acc, verbose=TRUE, rarefaction=TRUE)
    disparity_sli_acc<-list('disparity'=disp_sli_acc, 'diversity'=diversity_sli)
    save(disparity_sli_acc, file=paste(data_path, chain_name, '/',chain_name,'_disparity-sli_acc', mat_name,'-',distance ,'.Rda', sep=''))

    disp_sli_pro<-time.disparity(pco_sli_pro, verbose=TRUE, rarefaction=TRUE)
    disparity_sli_pro<-list('disparity'=disp_sli_pro, 'diversity'=diversity_sli)
    save(disparity_sli_pro, file=paste(data_path, chain_name, '/',chain_name,'_disparity-sli_pro', mat_name,'-',distance ,'.Rda', sep=''))
}
" > ${chain}-disparity-${dist_pre}-${distance}.R