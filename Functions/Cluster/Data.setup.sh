#!/bin/sh
##########################
#Shell script for setting up the data for disparity analysis.
##########################
#SYNTAX:
#sh Data.setup.sh <chain> <method> <moduleslist> <split>
#with:
#<chain> the name of the chain to generate task files for
#<path> the path where the data will be stored under the chain name
#<matrix> the path to the morphological matrix (can be relative)
#<tree> the path to the phylogenetic tree (can be relative)
#<ace> either "TRUE" to calculate it (two tasks) or the path to an ancestral state matrix. WARNING, if distance is set to "TRUE" only the tip distance can be calculated simultaneously.
#########################
#version 0.1
#----
#guillert(at)tcd.ie - 31/03/2015
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
tree\$node.label<-paste('n',seq(1:Nnode(tree)), sep='')
" > R_template.R

#ADD THE OPTIONS (WHAT TO CALCULATE)

#Ancestral states
if echo $ace | grep 'TRUE' > /dev/null
then
    #Use the R template
    cp R_template.R R_script_ace_ape.R
    cp R_template.R R_script_ace_claddis.R

    #Add the APE ace
    echo "
    ####################################
    #Ancestral states reconstruction
    ####################################
    anc_states<-anc.state(tree, Nexus_data, method='ML-ape', verbose=TRUE)
    save(anc_states, file=paste('../Data/',chain_name,'_ancestral_states-ape.Rda', sep=''))
    " >> R_script_ace_ape.R

    #Add the Claddis ace
    echo "
    ####################################
    #Ancestral states reconstruction
    ####################################
    anc_states<-anc.state(tree, Nexus_data, method='ML-claddis', verbose=TRUE)
    save(anc_states, file=paste(data_path, chain_name, '/',chain_name,'_ancestral_states-claddis.Rda', sep=''))
    " >> R_script_ace_claddis.R
else
    #Just add the ace path to the template
    echo"
    ####################################
    #Ancestral states reconstruction
    ####################################
    load('${ace}')
    " >> R_template.R
fi

#Distance matrices
if echo $ace | grep 'TRUE' > /dev/null
then
    cp R_template.R R_script_dist_tips.R
    cp R_template.R R_script_dist_nodes_ape.R
    cp R_template.R R_script_dist_nodes95_ape.R
    cp R_template.R R_script_dist_nodes_claddis.R
    cp R_template.R R_script_dist_nodes95_claddis.R
    
    #Add the distance script (tips)
    echo "
    ####################################
    #Distance matrix
    ####################################
    matrix_tips<-Nexus_data
    message('\\nCalculating the distance matrix for the tips only...', appendLF=FALSE)
    dist_tips<-MorphDistMatrix.verbose(matrix_tips, verbose=TRUE)
    message('Done.\\n', appendLF=FALSE)
    save(dist_tips, file=paste(data_path, chain_name, '/',chain_name,'_distance-tips.Rda', sep='')) #dist_tips
    " >> R_script_dist_tips.R

    #Add the distance script (nodes-ape)
    echo "
    ####################################
    #Distance matrix
    ####################################
    #Ancestral states
    load(paste(data_path, chain_name, '/',chain_name,'_ancestral_states-ape.Rda', sep=''))
    #Distance matrix using also nodes
    matrix_nodes<-Nexus_data
    matrix_nodes\$matrix<-anc_states\$state
    message('\\nCalculating the distance matrix for the tips and the nodes...', appendLF=FALSE)
    dist_nodes<-MorphDistMatrix.verbose(matrix_nodes, verbose=TRUE)
    message('Done.\\n', appendLF=FALSE)
    save(dist_nodes, file=paste(data_path, chain_name, '/',chain_name,'_distance-nodes-ape.Rda', sep='')) #dist_nodes
    " >> R_script_dist_nodes_ape.R

    #Add the distance script (nodes95-ape)
    echo "
    ####################################
    #Distance matrix
    ####################################
    #Ancestral states
    load(paste(data_path, chain_name, '/',chain_name,'_ancestral_states-ape.Rda', sep=''))
    #Distance matrix using also nodes
    matrix_nodes95<-Nexus_data
    matrix_nodes95\$matrix<-anc.unc(anc_states, 0.95, missing=NA)\$state
    message('\\nCalculating the distance matrix for the tips and the nodes with a 95 CI...', appendLF=FALSE)
    dist_nodes95<-MorphDistMatrix.verbose(matrix_nodes95, verbose=TRUE)
    message('Done.\\n', appendLF=FALSE)
    save(dist_nodes95, file=paste(data_path, chain_name, '/',chain_name,'_distance-nodes95-ape.Rda', sep='')) #dist_nodes95
    " >> R_script_dist_nodes95_ape.R

    #Add the distance script (nodes-claddis)
    echo "
    ####################################
    #Distance matrix
    ####################################
    #Ancestral states
    load(paste(data_path, chain_name, '/',chain_name,'_ancestral_states-claddis.Rda', sep=''))
    #Distance matrix using also nodes
    matrix_nodes<-Nexus_data
    matrix_nodes\$matrix<-anc_states\$state
    message('\\nCalculating the distance matrix for the tips and the nodes...', appendLF=FALSE)
    dist_nodes<-MorphDistMatrix.verbose(matrix_nodes, verbose=TRUE)
    message('Done.\\n', appendLF=FALSE)
    save(dist_nodes, file=paste(data_path, chain_name, '/',chain_name,'_distance-nodes-claddis.Rda', sep='')) #dist_nodes
    " >> R_script_dist_nodes_claddis.R

    #Add the distance script (nodes95-claddis)
    echo "
    ####################################
    #Distance matrix
    ####################################
    #Ancestral states
    load(paste(data_path, chain_name, '/',chain_name,'_ancestral_states-claddis.Rda', sep=''))
    #Distance matrix using also nodes
    matrix_nodes95<-Nexus_data
    matrix_nodes95\$matrix<-anc.unc(anc_states, 0.95, missing=NA)\$state
    message('\\nCalculating the distance matrix for the tips and the nodes with a 95 CI...', appendLF=FALSE)
    dist_nodes95<-MorphDistMatrix.verbose(matrix_nodes95, verbose=TRUE)
    message('Done.\\n', appendLF=FALSE)
    save(dist_nodes95, file=paste(data_path, chain_name, '/',chain_name,'_distance-nodes95-claddis.Rda', sep='')) #dist_nodes95
    " >> R_script_dist_nodes95_claddis.R


else
    cp R_template.R R_script_dist_tips.R
    cp R_template.R R_script_dist_nodes.R
    cp R_template.R R_script_dist_nodes95.R

    #Add the distance script (tips)
    echo "
    ####################################
    #Distance matrix
    ####################################
    matrix_tips<-Nexus_data
    message('\\nCalculating the distance matrix for the tips only...', appendLF=FALSE)
    dist_tips<-MorphDistMatrix.verbose(matrix_tips, verbose=TRUE)
    message('Done.\\n', appendLF=FALSE)
    save(dist_tips, file=paste(data_path, chain_name, '/',chain_name,'_distance-tips.Rda', sep='')) #dist_tips
    " >> R_script_dist_tips.R

    #Add the distance script (nodes)
    echo "
    ####################################
    #Distance matrix
    ####################################
    #Distance matrix using also nodes
    matrix_nodes<-Nexus_data
    matrix_nodes\$matrix<-anc_states\$state
    message('\\nCalculating the distance matrix for the tips and the nodes...', appendLF=FALSE)
    dist_nodes<-MorphDistMatrix.verbose(matrix_nodes, verbose=TRUE)
    message('Done.\\n', appendLF=FALSE)
    save(dist_nodes, file=paste(data_path, chain_name, '/',chain_name,'_distance-nodes.Rda', sep='')) #dist_nodes
    " >> R_script_dist_nodes_ape.R

    #Add the distance script (nodes95)
    echo "
    ####################################
    #Distance matrix
    ####################################
    #Distance matrix using also nodes
    matrix_nodes95<-Nexus_data
    matrix_nodes95\$matrix<-anc.unc(anc_states, 0.95, missing=NA)\$state
    message('\\nCalculating the distance matrix for the tips and the nodes with a 95 CI...', appendLF=FALSE)
    dist_nodes95<-MorphDistMatrix.verbose(matrix_nodes95, verbose=TRUE)
    message('Done.\\n', appendLF=FALSE)
    save(dist_nodes95, file=paste(data_path, chain_name, '/',chain_name,'_distance-nodes95.Rda', sep='')) #dist_nodes95
    " >> R_script_dist_nodes95_ape.R

fi

#CREATE THE CLUSTER R SCRIPT TASKER
#Preparing the R-codes jobs
if echo $ace | grep 'TRUE' > /dev/null
then
    #Shell R jobs
    echo "R --no-save < R_script_ace_ape" > R_script_ace_ape.sh
    echo "R --no-save < R_script_ace_claddis" > R_script_ace_claddis.sh
    echo "R --no-save < R_script_dist_tips.R" > R_script_dist_tips.sh
    echo "R --no-save < R_script_dist_nodes_ape.R" > R_script_dist_nodes_ape.sh
    echo "R --no-save < R_script_dist_nodes95_ape.R" > R_script_dist_nodes95_ape.sh
    echo "R --no-save < R_script_dist_nodes_claddis.R" > R_script_dist_nodes_claddis.sh
    echo "R --no-save < R_script_dist_nodes_claddis.R" > R_script_dist_nodes_claddis.sh
    #Shell R jobs config (1)
    echo "0 sh R_script_ace_ape.sh" >> Rjob${chain}-1.config
    echo "1 sh R_script_ace_claddis.sh" >> Rjob${chain}-1.config
    echo "2 sh R_script_dist_tips.sh" >> Rjob${chain}-1.config
    #Shell R jobs config (2)
    echo "0 sh R_script_dist_nodes_ape.sh" >> Rjob${chain}-2.config
    echo "1 sh R_script_dist_nodes95_ape.sh" >> Rjob${chain}-2.config
    echo "2 sh R_script_dist_nodes_claddis.sh" >> Rjob${chain}-2.config
    echo "3 sh R_script_dist_nodes95_claddis.sh" >> Rjob${chain}-2.config
    #Preparing the batch file (1)
    echo "#!/bin/sh
    #SBATCH -n 8
    #SBATCH -t 4-00:00:00
    #SBATCH -p compute
    #SBATCH -J ${chain}-1
    srun --multi-prog Rjob${chain}-1.config" > Rjob${chain}-1.sh
    #Preparing the batch file (2)
    echo "#!/bin/sh
    #SBATCH -n 8
    #SBATCH -t 4-00:00:00
    #SBATCH -p compute
    #SBATCH -J ${chain}-2
    srun --multi-prog Rjob${chain}-2.config" > Rjob${chain}-2.sh
    echo 'Running the R tasks:'
    echo 'sbatch Rjob${chain}-1.sh ; sbatch Rjob${chain}-2.sh'

else
    #Shell R jobs
    echo "R --no-save < R_script_dist_tips.R" > R_script_dist_tips.sh
    echo "R --no-save < R_script_dist_nodes.R" > R_script_dist_nodes.sh
    echo "R --no-save < R_script_dist_nodes95.R" > R_script_dist_nodes95.sh
    #Shell R jobs config
    echo "0 sh R_script_dist_tips.sh" >> Rjob${chain}.config
    echo "1 sh R_script_dist_nodes.sh" >> Rjob${chain}.config
    echo "2 sh R_script_dist_nodes95.sh" >> Rjob${chain}.config
    #Preparing the batch file
    echo "#!/bin/sh
    #SBATCH -n 8
    #SBATCH -t 4-00:00:00
    #SBATCH -p compute
    #SBATCH -J ${chain}
    srun --multi-prog Rjob${chain}.config" > Rjob${chain}.sh
    echo 'Running the R tasks:'
    echo 'sbatch Rjob${chain}.sh'
fi

#Remove template
rm R_template.R