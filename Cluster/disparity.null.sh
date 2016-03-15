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
#<intervals> a series of time intervals
#<slices> a series of time slices
#########################
#version 0.2.1
#Update: RAM friendly version
#Update: Faster cleaning part
#----
#guillert(at)tcd.ie - 18/05/2015
###########################

chain=$1
path=$2
matrix=$3
tree=$4
disparity=$5
type=$6
intervals=$7
slices=$8

######################################
#R template
######################################
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
intervals=as.numeric(strsplit(c(noquote('${intervals}')), split=',')[[1]])
slices=as.numeric(strsplit(c(noquote('${slices}')), split=',')[[1]])

#matrix
Nexus_data<-ReadMorphNexus(file_matrix)
Nexus_matrix<-Nexus_data\$matrix
#tree
Tree_data<-read.nexus(file_tree)

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
Nexus_data\$matrix<-table

#Forcing the tree to be binary
tree<-bin.tree(tree)

#Adding node labels to the tree
tree<-lapply.root(tree, max(tree.age(tree)\$age))

######################################
#Calculating null disparity 
######################################

rand_data<-null.data(tree=tree, matrix=apply(Nexus_data\$matrix, 2, states.count), matrix.model='${type}', replicates=<REPLICATES>, verbose=TRUE, include.nodes=TRUE)
ran_matrix<-lapply(rand_data, make.nexus)
#distance
dist_ran<-lapply(ran_matrix, MorphDistMatrix.verbose, verbose=TRUE)
dist_ran<-lapply(dist_ran, extract.dist, distance='gower.dist.matrix')
#pco
pco_data<-lapply(dist_ran, function(X) cmdscale(X, k=nrow(X)-1, add=TRUE)\$points)

######################################
#intervals
######################################

#intervals
pco_time<-lapply(pco_data, function(X) int.pco(X, tree, intervals, include.nodes=TRUE))
#diversity
div_data<-int.pco(pco_data[[1]], tree, intervals, include.nodes=TRUE, diversity=TRUE)[[2]]
#disparity
dis_data<-lapply(pco_time, time.disparity, verbose=TRUE, method='${disparity}', save.all=TRUE)
#Combine results
dis_data<-combine.disp(dis_data)
#list for output
disp_null_int<-list('diversity'=div_data, 'disparity'=dis_data)
save(disp_null_int, file=paste(data_path, chain_name, '/',chain_name,'-diparity_null-int-${type}_${disparity}<NAME>.Rda', sep=''))

#Saving RAM!
pco_time<-NULL
div_data<-NULL
dis_data<-NULL
disp_null_int<-NULL

######################################
#Slices
######################################

#slices - random
pco_time<-lapply(pco_data, function(X) slice.pco(X, tree, slices, verbose=TRUE, method='random'))
#diversity
div_data<-slice.pco(pco_data[[1]], tree, slices, verbose=TRUE, method='random')[[2]]
#disparity
dis_data<-lapply(pco_time, time.disparity, verbose=TRUE, method='${disparity}', save.all=TRUE)
#list for output
disp_null_int<-list('diversity'=div_data, 'disparity'=dis_data)
save(disp_null_int, file=paste(data_path, chain_name, '/',chain_name,'-diparity_null-sli_ran-${type}_${disparity}<NAME>.Rda', sep=''))

#slices - acctran
pco_time<-lapply(pco_data, function(X) slice.pco(X, tree, slices, verbose=TRUE, method='acctran'))
#disparity
dis_data<-lapply(pco_time, time.disparity, verbose=TRUE, method='${disparity}', save.all=TRUE)
#list for output
disp_null_int<-list('diversity'=div_data, 'disparity'=dis_data)
save(disp_null_int, file=paste(data_path, chain_name, '/',chain_name,'-diparity_null-sli_acc-${type}_${disparity}<NAME>.Rda', sep=''))

#slices - deltran
pco_time<-lapply(pco_data, function(X) slice.pco(X, tree, slices, verbose=TRUE, method='deltran'))
#disparity
dis_data<-lapply(pco_time, time.disparity, verbose=TRUE, method='${disparity}', save.all=TRUE)
#list for output
disp_null_int<-list('diversity'=div_data, 'disparity'=dis_data)
save(disp_null_int, file=paste(data_path, chain_name, '/',chain_name,'-diparity_null-sli_del-${type}_${disparity}<NAME>.Rda', sep=''))

#slices - proximity
pco_time<-lapply(pco_data, function(X) slice.pco(X, tree, slices, verbose=TRUE, method='proximity'))
#disparity
dis_data<-lapply(pco_time, time.disparity, verbose=TRUE, method='${disparity}', save.all=TRUE)
#list for output
disp_null_int<-list('diversity'=div_data, 'disparity'=dis_data)
save(disp_null_int, file=paste(data_path, chain_name, '/',chain_name,'-diparity_null-sli_pro-${type}_${disparity}<NAME>.Rda', sep=''))
" > disparity.null.template.tmp

######################################
#Tasker
######################################

for ((rep=1 ;  rep<=24 ;  rep+=1))
do  
    #make the replicate files
    sed 's/<REPLICATES>/5/g' disparity.null.template.tmp | sed 's/<NAME>/'"${rep}"'/g' > ${chain}-${type}_${disparity}.${rep}.R
    #make the shell files
    echo "R --no-save < ${chain}-${type}_${disparity}.${rep}.R > /dev/null" > ${chain}-${type}_${disparity}.${rep}.sh
    #making the config file
    n=$(( $rep - 1 )) 
    echo "$n sh ${chain}-${type}_${disparity}.${rep}.sh" >> ${chain}-${type}_${disparity}.config
done

#Remove the template
rm disparity.null.template.tmp

echo "#!/bin/sh
#SBATCH -n 24
#SBATCH -t 4-00:00:00
#SBATCH -p compute
#SBATCH -J D-${chain}
source /etc/profile.d/modules.sh
export http_proxy=http://proxy.tchpc.tcd.ie:8080
srun --multi-prog ${chain}-${type}_${disparity}.config" > ${chain}-${type}_${disparity}-launcher.sh

echo 'Running the R tasks:'
echo "sbatch ${chain}-${type}_${disparity}-launcher.sh"