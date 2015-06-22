#!/bin/sh
##########################
#Shell script for Calculating the disparity per chain for a cluster
##########################
#SYNTAX:
#sh disparity.tasker.sh <chain> <path> <matrix> <tree> <intervals> <slices> <FADLAD>
#with:
#<chain> the name of the chain to generate task files for
#<path> the path where the data will be stored under the chain name
#<matrix> the path to the morphological matrix (can be relative)
#<tree> the path to the phylogenetic tree (can be relative)
#<intervals> a series of time intervals
#<slices> a series of time slices
#<FADLAD> a csv file with FADLAD data
#########################
#version 0.2.2
#Update: Faster cleaning part
#Update: fixed config file
#Update: fixed error in slices files
#Update: added observed data calculation
#Update: fixed diversity bug + added diversity count saving
#Update: add the three centroid disparity scenarios
#----
#guillert(at)tcd.ie - 22/06/2015
###########################

#INPUT
chain=$1
path=$2
matrix=$3
tree=$4
intervals=$5
slices=$6
FADLAD=$7

#Script storing folder
rm -R R_scripts_${chain}
mkdir R_scripts_${chain}

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
FADLAD='${FADLAD}'

#matrix
Nexus_data<-ReadMorphNexus(file_matrix)
Nexus_matrix<-Nexus_data\$matrix
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
Nexus_data\$matrix<-table

#Forcing the tree to be binary
tree<-bin.tree(tree)

#Adding node labels to the tree
tree<-lapply.root(tree, max(tree.age(tree)\$age)) " > R_scripts_${chain}/R_script_template.tmp

######################################
#Distance matrices
######################################
cp R_scripts_${chain}/R_script_template.tmp R_scripts_${chain}/R_disparity_tips.tmp
cp R_scripts_${chain}/R_script_template.tmp R_scripts_${chain}/R_disparity_nodes.tmp
cp R_scripts_${chain}/R_script_template.tmp R_scripts_${chain}/R_disparity_nodes95.tmp
echo "
#load the distance matrix
load(paste(data_path, chain_name, '/', chain_name, '_distance-tips.Rda', sep='')) #dist_tips
trimmed_max_data_tips<-TrimMorphDistMatrix(dist_tips\$max.dist.matrix)
tree_tips<-drop.tip(tree, trimmed_max_data_tips\$removed.taxa) ; tree_tips\$root.time<-max(tree.age(tree_tips)[,1])
#pco
pco_data_tips<-cmdscale(trimmed_max_data_tips\$dist.matrix, k=nrow(trimmed_max_data_tips\$dist.matrix) - 2, add=T)\$points
" >> R_scripts_${chain}/R_disparity_tips.tmp
echo "
#load the distance matrix
load(paste(data_path, chain_name, '/', chain_name, '_distance-nodes.Rda', sep='')) #dist_nodes
trimmed_max_data_nodes<-TrimMorphDistMatrix(dist_nodes\$max.dist.matrix)
tree_nodes<-drop.tip(tree, trimmed_max_data_nodes\$removed.taxa) ; tree_nodes\$root.time<-max(tree.age(tree_nodes)[,1])
trimmed_max_data_nodes\$dist.matrix<-trimmed_max_data_nodes\$dist.matrix[c(tree_nodes\$tip.label, tree_nodes\$node.label),c(tree_nodes\$tip.label, tree_nodes\$node.label)]
#pco
pco_data_nodes<-cmdscale(trimmed_max_data_nodes\$dist.matrix, k=nrow(trimmed_max_data_nodes\$dist.matrix) - 2, add=T)\$points
" >> R_scripts_${chain}/R_disparity_nodes.tmp
echo "
#load the distance matrix
load(paste(data_path, chain_name, '/', chain_name, '_distance-nodes95.Rda', sep='')) #dist_nodes95
trimmed_max_data_nodes95<-TrimMorphDistMatrix(dist_nodes95\$max.dist.matrix)
tree_nodes95<-drop.tip(tree, trimmed_max_data_nodes95\$removed.taxa) ; tree_nodes95\$root.time<-max(tree.age(tree_nodes95)[,1])
trimmed_max_data_nodes95\$dist.matrix<-trimmed_max_data_nodes95\$dist.matrix[c(tree_nodes95\$tip.label, tree_nodes95\$node.label),c(tree_nodes95\$tip.label, tree_nodes95\$node.label)]
#pco
pco_data_nodes95<-cmdscale(trimmed_max_data_nodes95\$dist.matrix, k=nrow(trimmed_max_data_nodes95\$dist.matrix) - 2, add=T)\$points
" >> R_scripts_${chain}/R_disparity_nodes95.tmp

######################################
#Intervals
######################################
cp R_scripts_${chain}/R_disparity_tips.tmp R_scripts_${chain}/R_disparity_tips_int.R
cp R_scripts_${chain}/R_disparity_nodes.tmp R_scripts_${chain}/R_disparity_nodes_int.R
cp R_scripts_${chain}/R_disparity_nodes95.tmp R_scripts_${chain}/R_disparity_nodes95_int.R

echo "
#Intervals
pco_int_tips<-int.pco(pco_data_tips, tree_tips, intervals, include.nodes=FALSE, FAD_LAD=FADLAD, diversity=TRUE)
int_tips_div<-pco_int_tips[[2]] ; pco_int_tips<-pco_int_tips[[1]]
save(int_tips_div,file=paste(data_path, chain_name, '/',chain_name,'-int_tips_div.Rda', sep=''))
" >> R_scripts_${chain}/R_disparity_tips_int.R

echo "
#Intervals
pco_int_nodes<-int.pco(pco_data_nodes, tree_nodes, intervals, include.nodes=TRUE, FAD_LAD=FADLAD, diversity=TRUE)
int_nodes_div<-pco_int_nodes[[2]] ; pco_int_nodes<-pco_int_nodes[[1]]
save(int_nodes_div,file=paste(data_path, chain_name, '/',chain_name,'-int_nodes_div.Rda', sep=''))
" >> R_scripts_${chain}/R_disparity_nodes_int.R

echo "
#Intervals
pco_int_nodes95<-int.pco(pco_data_nodes95, tree_nodes95, int_breaks, include.nodes=TRUE, FAD_LAD=FADLAD, diversity=TRUE)
int_nodes95_div<-pco_int_nodes95[[2]] ; pco_int_nodes95<-pco_int_nodes95[[1]]
save(int_nodes95_div,file=paste(data_path, chain_name, '/',chain_name,'-int_nodes95_div.Rda', sep=''))
" >> R_scripts_${chain}/R_disparity_nodes95_int.R

######################################
#Slices
######################################
cp R_scripts_${chain}/R_disparity_nodes.tmp R_scripts_${chain}/R_disparity_nodes_sli_ran.R
cp R_scripts_${chain}/R_disparity_nodes.tmp R_scripts_${chain}/R_disparity_nodes_sli_acc.R
cp R_scripts_${chain}/R_disparity_nodes.tmp R_scripts_${chain}/R_disparity_nodes_sli_del.R
cp R_scripts_${chain}/R_disparity_nodes.tmp R_scripts_${chain}/R_disparity_nodes_sli_pro.R
cp R_scripts_${chain}/R_disparity_nodes95.tmp R_scripts_${chain}/R_disparity_nodes95_sli_ran.R
cp R_scripts_${chain}/R_disparity_nodes95.tmp R_scripts_${chain}/R_disparity_nodes95_sli_acc.R
cp R_scripts_${chain}/R_disparity_nodes95.tmp R_scripts_${chain}/R_disparity_nodes95_sli_del.R
cp R_scripts_${chain}/R_disparity_nodes95.tmp R_scripts_${chain}/R_disparity_nodes95_sli_pro.R

echo "
#slices
pco_slices_nodes_ran<-slice.pco(pco_data_nodes, tree_nodes, slices, method='random', FAD_LAD=FADLAD, verbose=TRUE, diversity=TRUE)
slices_nodes_div<-pco_slices_nodes_ran[[2]] ; pco_slices_nodes_ran<-pco_slices_nodes_ran[[1]]
save(slices_nodes_div,file=paste(data_path, chain_name, '/',chain_name,'-slices_nodes_div.Rda', sep=''))
" >> R_scripts_${chain}/R_disparity_nodes_sli_ran.R

echo "
#slices
pco_slices_nodes_acc<-slice.pco(pco_data_nodes, tree_nodes, slices, method='acctran', FAD_LAD=FADLAD, verbose=TRUE, diversity=TRUE)
slices_nodes_div<-pco_slices_nodes_acc[[2]] ; pco_slices_nodes_acc<-pco_slices_nodes_acc[[1]]
" >> R_scripts_${chain}/R_disparity_nodes_sli_acc.R

echo "
#slices
pco_slices_nodes_del<-slice.pco(pco_data_nodes, tree_nodes, slices, method='deltran', FAD_LAD=FADLAD, verbose=TRUE, diversity=TRUE)
slices_nodes_div<-pco_slices_nodes_del[[2]] ; pco_slices_nodes_del<-pco_slices_nodes_del[[1]]
" >> R_scripts_${chain}/R_disparity_nodes_sli_del.R

echo "
#slices
pco_slices_nodes_pro<-slice.pco(pco_data_nodes, tree_nodes, slices, method='proximity', FAD_LAD=FADLAD, verbose=TRUE, diversity=TRUE)
slices_nodes_div<-pco_slices_nodes_pro[[2]] ; pco_slices_nodes_pro<-pco_slices_nodes_pro[[1]]
" >> R_scripts_${chain}/R_disparity_nodes_sli_pro.R

echo "
#slices
pco_slices_nodes95_ran<-slice.pco(pco_data_nodes95, tree_nodes95, slices, method='random', FAD_LAD=FADLAD, verbose=TRUE, diversity=TRUE)
slices_nodes95_div<-pco_slices_nodes95_ran[[2]] ; pco_slices_nodes95_ran<-pco_slices_nodes95_ran[[1]]
save(slices_nodes95_div,file=paste(data_path, chain_name, '/',chain_name,'-slices_nodes95_div.Rda', sep=''))
" >> R_scripts_${chain}/R_disparity_nodes95_sli_ran.R

echo "
#slices
pco_slices_nodes95_acc<-slice.pco(pco_data_nodes95, tree_nodes95, slices, method='acctran', FAD_LAD=FADLAD, verbose=TRUE, diversity=TRUE)
slices_nodes95_div<-pco_slices_nodes95_acc[[2]] ; pco_slices_nodes95_acc<-pco_slices_nodes95_acc[[1]]
" >> R_scripts_${chain}/R_disparity_nodes95_sli_acc.R

echo "
#slices
pco_slices_nodes95_del<-slice.pco(pco_data_nodes95, tree_nodes95, slices, method='deltran', FAD_LAD=FADLAD, verbose=TRUE, diversity=TRUE)
slices_nodes95_div<-pco_slices_nodes95_del[[2]] ; pco_slices_nodes95_del<-pco_slices_nodes95_del[[1]]
" >> R_scripts_${chain}/R_disparity_nodes95_sli_del.R

echo "
#slices
pco_slices_nodes95_pro<-slice.pco(pco_data_nodes95, tree_nodes95, slices, method='proximity', FAD_LAD=FADLAD, verbose=TRUE, diversity=TRUE)
slices_nodes95_div<-pco_slices_nodes95_pro[[2]] ; pco_slices_nodes95_pro<-pco_slices_nodes95_pro[[1]]
" >> R_scripts_${chain}/R_disparity_nodes95_sli_pro.R

######################################
#Disparity
######################################
echo "
#Disparity
disp_int_tips<-time.disparity(pco_int_tips, verbose=TRUE, rarefaction=TRUE, save.all=TRUE)
save(disp_int_tips, file=paste(data_path, chain_name, '/',chain_name,'-disp_int_tips.Rda', sep=''))
#Observed disparity
disp_int_tips_obs<-time.disparity(pco_int_tips, method='centroid', bootstraps=0, verbose=TRUE, rarefaction=TRUE, save.all=TRUE, centroid.type='full')
save(disp_int_tips_obs,file=paste(data_path, chain_name, '/',chain_name,'-disp_int_tips_obs.Rda', sep=''))
#Observed disparity (BS)
disp_int_tips_obs_BS<-time.disparity(pco_int_tips, method='centroid', bootstraps=1000, verbose=TRUE, rarefaction=TRUE, save.all=TRUE, centroid.type='full')
save(disp_int_tips_obs_BS,file=paste(data_path, chain_name, '/',chain_name,'-disp_int_tips_obs_BS.Rda', sep=''))

" >> R_scripts_${chain}/R_disparity_tips_int.R

echo "
#Disparity
disp_int_nodes<-time.disparity(pco_int_nodes, verbose=TRUE, rarefaction=TRUE, save.all=TRUE)
save(disp_int_nodes, file=paste(data_path, chain_name, '/',chain_name,'-disp_int_nodes.Rda', sep=''))
#Observed disparity
disp_int_nodes_obs<-time.disparity(pco_int_nodes, method='centroid', bootstraps=0, verbose=TRUE, rarefaction=TRUE, save.all=TRUE, centroid.type='full')
save(disp_int_nodes_obs,file=paste(data_path, chain_name, '/',chain_name,'-disp_int_nodes_obs.Rda', sep=''))
#Observed disparity (BS)
disp_int_nodes_obs_BS<-time.disparity(pco_int_nodes, method='centroid', bootstraps=1000, verbose=TRUE, rarefaction=TRUE, save.all=TRUE, centroid.type='full')
save(disp_int_nodes_obs_BS,file=paste(data_path, chain_name, '/',chain_name,'-disp_int_nodes_obs_BS.Rda', sep=''))
" >> R_scripts_${chain}/R_disparity_nodes_int.R

echo "
#Disparity
disp_int_nodes95<-time.disparity(pco_int_nodes95, verbose=TRUE, rarefaction=TRUE, save.all=TRUE)
save(disp_int_nodes95, file=paste(data_path, chain_name, '/',chain_name,'-disp_int_nodes95.Rda', sep=''))
#Observed disparity
disp_int_nodes95_obs<-time.disparity(pco_int_nodes95, method='centroid', bootstraps=0, verbose=TRUE, rarefaction=TRUE, save.all=TRUE, centroid.type='full')
save(disp_int_nodes95_obs,file=paste(data_path, chain_name, '/',chain_name,'-disp_int_nodes95_obs.Rda', sep=''))
disp_int_nodes95_obs_BS<-time.disparity(pco_int_nodes95, method='centroid', bootstraps=1000, verbose=TRUE, rarefaction=TRUE, save.all=TRUE, centroid.type='full')
save(disp_int_nodes95_obs_BS,file=paste(data_path, chain_name, '/',chain_name,'-disp_int_nodes95_obs_BS.Rda', sep=''))
" >> R_scripts_${chain}/R_disparity_nodes95_int.R

echo "
#Disparity
disp_sli_nodes_ran<-time.disparity(pco_slices_nodes_ran, verbose=TRUE, rarefaction=TRUE, save.all=TRUE)
save(disp_sli_nodes_ran, file=paste(data_path, chain_name, '/',chain_name,'-disp_sli_nodes_ran.Rda', sep=''))
#Observed disparity
disp_sli_nodes_ran_obs<-time.disparity(pco_slices_nodes_ran, method='centroid', bootstraps=0, verbose=TRUE, rarefaction=TRUE, save.all=TRUE, centroid.type='full')
disp_sli_nodes_ran_obs,file=paste(data_path, chain_name, '/',chain_name,'-disp_sli_nodes_ran_obs.Rda', sep=''))
#Observed disparity (BS)
disp_sli_nodes_ran_obs_BS<-time.disparity(pco_slices_nodes_ran, method='centroid', bootstraps=1000, verbose=TRUE, rarefaction=TRUE, save.all=TRUE, centroid.type='full')
save(disp_sli_nodes_ran_obs_BS,file=paste(data_path, chain_name, '/',chain_name,'-disp_sli_nodes_ran_obs_BS.Rda', sep=''))
" >> R_scripts_${chain}/R_disparity_nodes_sli_ran.R

echo "
#Disparity
disp_sli_nodes_del<-time.disparity(pco_slices_nodes_del, verbose=TRUE, rarefaction=TRUE, save.all=TRUE)
save(disp_sli_nodes_del, file=paste(data_path, chain_name, '/',chain_name,'-disp_sli_nodes_del.Rda', sep=''))
#Observed disparity
disp_sli_nodes_del_obs<-time.disparity(pco_slices_nodes_del, method='centroid', bootstraps=0, verbose=TRUE, rarefaction=TRUE, save.all=TRUE, centroid.type='full')
disp_sli_nodes_del_obs,file=paste(data_path, chain_name, '/',chain_name,'-disp_sli_nodes_del_obs.Rda', sep=''))
#Observed disparity (BS)
disp_sli_nodes_del_obs_BS<-time.disparity(pco_slices_nodes_del, method='centroid', bootstraps=1000, verbose=TRUE, rarefaction=TRUE, save.all=TRUE, centroid.type='full')
save(disp_sli_nodes_del_obs_BS,file=paste(data_path, chain_name, '/',chain_name,'-disp_sli_nodes_del_obs_BS.Rda', sep=''))
" >> R_scripts_${chain}/R_disparity_nodes_sli_del.R

echo "
#Disparity
disp_sli_nodes_acc<-time.disparity(pco_slices_nodes_acc, verbose=TRUE, rarefaction=TRUE, save.all=TRUE)
save(disp_sli_nodes_acc, file=paste(data_path, chain_name, '/',chain_name,'-disp_sli_nodes_acc.Rda', sep=''))
#Observed disparity
disp_sli_nodes_acc_obs<-time.disparity(pco_slices_nodes_acc, method='centroid', bootstraps=0, verbose=TRUE, rarefaction=TRUE, save.all=TRUE, centroid.type='full')
disp_sli_nodes_acc_obs,file=paste(data_path, chain_name, '/',chain_name,'-disp_sli_nodes_acc_obs.Rda', sep=''))
#Observed disparity (BS)
disp_sli_nodes_acc_obs_BS<-time.disparity(pco_slices_nodes_acc, method='centroid', bootstraps=1000, verbose=TRUE, rarefaction=TRUE, save.all=TRUE, centroid.type='full')
save(disp_sli_nodes_acc_obs_BS,file=paste(data_path, chain_name, '/',chain_name,'-disp_sli_nodes_acc_obs_BS.Rda', sep=''))
" >> R_scripts_${chain}/R_disparity_nodes_sli_acc.R

echo "
#Disparity
disp_sli_nodes_pro<-time.disparity(pco_slices_nodes_pro, verbose=TRUE, rarefaction=TRUE, save.all=TRUE)
save(disp_sli_nodes_pro, file=paste(data_path, chain_name, '/',chain_name,'-disp_sli_nodes_pro.Rda', sep=''))
#Observed disparity
disp_sli_nodes_pro_obs<-time.disparity(pco_slices_nodes_pro, method='centroid', bootstraps=0, verbose=TRUE, rarefaction=TRUE, save.all=TRUE, centroid.type='full')
disp_sli_nodes_pro_obs,file=paste(data_path, chain_name, '/',chain_name,'-disp_sli_nodes_pro_obs.Rda', sep=''))
#Observed disparity (BS)
disp_sli_nodes_pro_obs_BS<-time.disparity(pco_slices_nodes_pro, method='centroid', bootstraps=1000, verbose=TRUE, rarefaction=TRUE, save.all=TRUE, centroid.type='full')
save(disp_sli_nodes_pro_obs_BS,file=paste(data_path, chain_name, '/',chain_name,'-disp_sli_nodes_pro_obs_BS.Rda', sep=''))
" >> R_scripts_${chain}/R_disparity_nodes_sli_pro.R

echo "
#Disparity
disp_sli_nodes95_ran<-time.disparity(pco_slices_nodes95_ran, verbose=TRUE, rarefaction=TRUE, save.all=TRUE)
save(disp_sli_nodes95_ran, file=paste(data_path, chain_name, '/',chain_name,'-disp_sli_nodes95_ran.Rda', sep=''))
#Observed disparity
disp_sli_nodes95_ran_obs<-time.disparity(pco_slices_nodes95_ran, method='centroid', bootstraps=0, verbose=TRUE, rarefaction=TRUE, save.all=TRUE, centroid.type='full')
disp_sli_nodes95_ran_obs,file=paste(data_path, chain_name, '/',chain_name,'-disp_sli_nodes95_ran_obs.Rda', sep=''))
#Observed disparity (BS)
disp_sli_nodes95_ran_obs_BS<-time.disparity(pco_slices_nodes95_ran, method='centroid', bootstraps=1000, verbose=TRUE, rarefaction=TRUE, save.all=TRUE, centroid.type='full')
save(disp_sli_nodes95_ran_obs_BS,file=paste(data_path, chain_name, '/',chain_name,'-disp_sli_nodes95_ran_obs_BS.Rda', sep=''))
" >> R_scripts_${chain}/R_disparity_nodes95_sli_ran.R

echo "
#Disparity
disp_sli_nodes95_del<-time.disparity(pco_slices_nodes95_del, verbose=TRUE, rarefaction=TRUE, save.all=TRUE)
save(disp_sli_nodes95_del, file=paste(data_path, chain_name, '/',chain_name,'-disp_sli_nodes95_del.Rda', sep=''))
#Observed disparity
disp_sli_nodes95_del_obs<-time.disparity(pco_slices_nodes95_del, method='centroid', bootstraps=0, verbose=TRUE, rarefaction=TRUE, save.all=TRUE, centroid.type='full')
disp_sli_nodes95_del_obs,file=paste(data_path, chain_name, '/',chain_name,'-disp_sli_nodes95_del_obs.Rda', sep=''))
#Observed disparity (BS)
disp_sli_nodes95_del_obs_BS<-time.disparity(pco_slices_nodes95_del, method='centroid', bootstraps=1000, verbose=TRUE, rarefaction=TRUE, save.all=TRUE, centroid.type='full')
save(disp_sli_nodes95_del_obs_BS,file=paste(data_path, chain_name, '/',chain_name,'-disp_sli_nodes95_del_obs_BS.Rda', sep=''))
" >> R_scripts_${chain}/R_disparity_nodes95_sli_del.R

echo "
#Disparity
disp_sli_nodes95_acc<-time.disparity(pco_slices_nodes95_acc, verbose=TRUE, rarefaction=TRUE, save.all=TRUE)
save(disp_sli_nodes95_acc, file=paste(data_path, chain_name, '/',chain_name,'-disp_sli_nodes95_acc.Rda', sep=''))
#Observed disparity
disp_sli_nodes95_acc_obs<-time.disparity(pco_slices_nodes95_acc, method='centroid', bootstraps=0, verbose=TRUE, rarefaction=TRUE, save.all=TRUE, centroid.type='full')
disp_sli_nodes95_acc_obs,file=paste(data_path, chain_name, '/',chain_name,'-disp_sli_nodes95_acc_obs.Rda', sep=''))
#Observed disparity (BS)
disp_sli_nodes95_acc_obs_BS<-time.disparity(pco_slices_nodes95_acc, method='centroid', bootstraps=1000, verbose=TRUE, rarefaction=TRUE, save.all=TRUE, centroid.type='full')
save(disp_sli_nodes95_acc_obs_BS,file=paste(data_path, chain_name, '/',chain_name,'-disp_sli_nodes95_acc_obs_BS.Rda', sep=''))
" >> R_scripts_${chain}/R_disparity_nodes95_sli_acc.R

echo "
#Disparity
disp_sli_nodes95_pro<-time.disparity(pco_slices_nodes95_pro, verbose=TRUE, rarefaction=TRUE, save.all=TRUE)
save(disp_sli_nodes95_pro, file=paste(data_path, chain_name, '/',chain_name,'-disp_sli_nodes95_pro.Rda', sep=''))
#Observed disparity
disp_sli_nodes95_pro_obs<-time.disparity(pco_slices_nodes95_pro, method='centroid', bootstraps=0, verbose=TRUE, rarefaction=TRUE, save.all=TRUE, centroid.type='full')
disp_sli_nodes95_pro_obs,file=paste(data_path, chain_name, '/',chain_name,'-disp_sli_nodes95_pro_obs.Rda', sep=''))
#Observed disparity (BS)
disp_sli_nodes95_pro_obs_BS<-time.disparity(pco_slices_nodes95_pro, method='centroid', bootstraps=1000, verbose=TRUE, rarefaction=TRUE, save.all=TRUE, centroid.type='full')
save(disp_sli_nodes95_pro_obs_BS,file=paste(data_path, chain_name, '/',chain_name,'-disp_sli_nodes95_pro_obs_BS.Rda', sep=''))
" >> R_scripts_${chain}/R_disparity_nodes95_sli_pro.R


#Remove temporary scripts
rm R_scripts_${chain}/*.tmp

######################################
#Tasker
######################################

#Shell R jobs
n=0
for f in R_scripts_${chain}/*.R
do
    script=$(basename ${f} .R)
    echo "R --no-save < ${script}.R" > R_scripts_${chain}/${script}.sh
    echo "$n sh ${script}.sh" >> R_scripts_${chain}/Rdisparity${chain}.config
    n=$(( $n + 1 )) 
done


echo "#!/bin/sh
#SBATCH -n 11
#SBATCH -t 1-00:00:00
#SBATCH -p compute
#SBATCH -J D-${chain}
srun --multi-prog Rdisparity${chain}.config" > R_scripts_${chain}/Rdisparity${chain}.sh

echo 'Running the R tasks:'
echo "sbatch R_scripts_${chain}/Rdisparity${chain}.sh"


