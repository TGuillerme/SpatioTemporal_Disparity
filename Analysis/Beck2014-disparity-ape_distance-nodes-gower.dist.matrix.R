
#Load the functions and the packages
library(disparity)

######################
#Reading the files
######################

#Selecting the file
chain_name='Beck2014'
data_path='../Data/'
file_matrix='../Data/2014-Beck-ProcB-matrix-morpho.nex'
file_tree='../Data/2014-Beck-ProcB-TEM.tre'
file_dist='../Data/Beck2014/Beck2014ape_distance-nodes.Rda'
distance='gower.dist.matrix'
intervals='c(150,140,130,120,110,100,90,80,70,60,50,40,30,20,10,0)'
slices='c(150,148,146,144,142,140,138,136,134,132,130,128,126)'
FADLAD='../Data/Beck2014_FADLAD.csv'

#matrix
Nexus_data<-ReadMorphNexus(file_matrix)
Nexus_matrix<-Nexus_data$matrix
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
Nexus_data$matrix<-table

#Forcing the tree to be binary
tree<-bin.tree(tree)

#Adding node labels to the tree
tree<-lapply.root(tree, max(tree.age(tree)$age)) 

#Load the FAD/LAD file
FADLAD<-read.csv(paste(data_path, chain_name, '_FADLAD.csv', sep=''), row.names=1)

######################
#Loading the distance matrix
######################
mat_name<-load(file_dist)
distance<-extract.dist(get(mat_name), distance)

######################
#Trim the matrix and the tree (if necessary)
######################
trimmed_data<-TrimMorphDistMatrix(distance)
tree<-(tree, trimmed_data$removed.taxa) ; tree$root.time<-max(tree.age(tree_tips)[,1])

######################
#pco
######################
pco_data<-cmdscale(trimmed_data$dist.matrix, k=nrow(trimmed_data$dist.matrix) - 1, add=T)$points

######################
#Disparity
######################

#Intervals

#Test whether to include the nodes or not
if(nrow(trimmed_data$dist.matrix) == (Ntip(tree) + Nnode(tree))) {
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

