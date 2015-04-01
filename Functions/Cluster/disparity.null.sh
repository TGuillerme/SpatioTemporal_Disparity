#!/bin/sh
##########################
#Shell script for Calculating the disparity null models (using 8 cores)
##########################
#SYNTAX:
#sh disparity.null.sh <chain> <path> <matrix> <tree> <disparity> <type>
#with:
#<chain> the name of the chain to generate task files for
#<path> the path where the data will be stored under the chain name
#<matrix> the path to the morphological matrix (can be relative)
#<tree> the path to the phylogenetic tree (can be relative)
#<disparity> the disparity metric
#<type> either random or sim.char
#<iterations> a number of iterations
#########################
#version 0.1
#----
#guillert(at)tcd.ie - 31/03/2015
###########################

chain=$1
path=$2
matrix=$3
tree=$4
disparity=$5
type=$6

######################################
#R template
######################################
echo "
#Load the functions and the packages
source('functions.R')

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
tree<-lapply.root(tree, max(tree.age(tree)\$age))

######################################
#Calculating null disparity
######################################

rand_data<-null.data(tree=step_tree, matrix=apply(Nexus_data$matrix, 2, states.count), matrix.model='${type}', replicates=<REPLICATES>, verbose=TRUE, include.nodes=TRUE)
ran_matrix<-lapply(rand_data, make.nexus)
#distance
dist_ran<-lapply(ran_matrix, MorphDistMatrix.verbose, verbose=TRUE)
dist_ran<-lapply(dist_ran, extract.dist, distance='max.dist.matrix')
#pco
pco_ran<-lapply(dist_ran, function(X) cmdscale(X, k=nrow(X)-1, add=TRUE)\$points)
#intervals
pco_ran_int<-lapply(pco_ran, function(X) int.pco(X, step_tree, intervals, include.nodes=TRUE))
#diversity
div_ran_int<-int.pco(pco_ran[[1]], step_tree, intervals, include.nodes=TRUE, diversity=TRUE)[[2]]
#disparity
disp_ran_int<-lapply(pco_ran_int, time.disparity, verbose=TRUE, method='${disparity}', save.all=TRUE)
#Combine results
disp_ran_int<-combine.disp(disp_ran_int)
save(disp_ran_int, file=paste(data_path, chain_name, '/',chain_name,'-diparity_null_${type}_${disparity}<NAME>.Rda', sep=''))
" > disparity.null.template.tmp

######################################
#Tasker
######################################

for ((rep=1 ;  rep<=8 ;  rep+=1))
do  
    #make the replicate files
    sed 's/<REPLICATES>/13/g' disparity.null.template.tmp | sed 's/<NAME>/'"${rep}"'/g' > ${chain}-${type}_${disparity}.${rep}.R
    #make the shell files
    echo "R --no-save < ${chain}-${type}_${disparity}.${rep}.R" > ${chain}-${type}_${disparity}.${rep}.sh
    #making the config file
    n=$(( $rep - 1 )) 
    echo "$n sh ${chain}-${type}_${disparity}.${rep}.sh" >> ${chain}-${type}_${disparity}.config
done

#Remove the template
rm disparity.null.template.tmp

echo "#!/bin/sh
#SBATCH -n 8
#SBATCH -t 4-00:00:00
#SBATCH -p compute
#SBATCH -J D-${chain}
srun --multi-prog ${chain}-null_${disparity}.config" > ${chain}-${type}_${disparity}-launcher.sh

echo 'Running the R tasks:'
echo "sbatch ${chain}-${type}_${disparity}-launcher.sh"