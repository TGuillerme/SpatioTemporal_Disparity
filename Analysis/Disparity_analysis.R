#Loading the package
library(devtools)
install_github("TGuillerme/SpatioTemporal_Disparity/Functions/disparity")
library(disparity)

###################################
#
# MAIN FIGURE
#
###################################

#Data extraction function
read.data<-function(chain_name, data_path, file_matrix, file_tree, disparity_data) {
    #matrix
    Nexus_data<-ReadMorphNexus(file_matrix)
    Nexus_matrix<-Nexus_data$matrix
    #tree
    Tree_data<-read.nexus(file_tree)

    ######################################
    # Loading the data
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

    #Adding node labels and the age to the tree
    tree<<-lapply.root(tree, max(tree.age(tree)$age))

    #Load the disparity data
    name<-load(paste(data_path, chain_name, "/", chain_name, disparity_data, sep=""))
    dis_tmp<<-get(name)
}

######################################
# Isolating the data
######################################

#Selecting the slater data
chain_name<-"Slater2013"
data_path<-"../Data/"
file_matrix<-"../Data/2013-Slater-MEE-matrix-morpho.nex"
file_tree<-"../Data/2013-Slater-MEE-TEM.tre"
disparity_data<-"_disparity-sli_prodist_nodes95-gower.dist.matrix.Rda"
#Extracting all the data
read.data(chain_name, data_path, file_matrix, file_tree, disparity_data)

#Isolating the slater data
disparity_full_slater<-dis_tmp$disparity
#Isolating the slater diversity
diversity_full_slater<-dis_tmp$diversity
#Isolating the slater tree
tree_slater<-tree

#Selecting the Beck data
chain_name<-"Beck2014"
data_path<-"../Data/"
file_matrix<-"../Data/2014-Beck-ProcB-matrix-morpho.nex"
file_tree<-"../Data/2014-Beck-ProcB-TEM.tre"
#Extracting all the data
read.data(chain_name, data_path, file_matrix, file_tree, disparity_data)

#Isolating the slater data
disparity_full_beck<-dis_tmp$disparity
#Isolating the slater diversity
diversity_full_beck<-dis_tmp$diversity
#Isolating the slater tree
tree_beck<-tree

######################################
# Tree visualisation
######################################

#Plot the tree (slater)
    geoscalePhylo(ladderize(tree_slater), cex.age=0.6, cex.ts=0.8, cex.tip=0.5)
    dev.new()
    geoscalePhylo(ladderize(tree_beck), cex.age=0.6, cex.ts=0.8, cex.tip=0.5)

######################################
# Disparity visualisation
######################################

dis_ran_max_beck<-extract.disp(disparity_full_beck$quantiles, rarefaction="max")
dis_ran_max_slater<-extract.disp(disparity_full_slater$quantiles, rarefaction="max")

#Pretty plot with the tree and the disparity
op<-par(mfrow=c(2,1), bty="n", mar=c(4,4,4,4))
plot.disparity(dis_ran_max_beck, diversity=log(dis_ran_max_beck$rarefaction), main="Eutherian", xlab="Time (Mya)", y2lab="Diversity (log)")
abline(v=22, col="red")
plot.disparity(dis_ran_max_slater, diversity=log(dis_ran_max_slater$rarefaction), main="Mammaliformes", xlab="Time (Mya)", y2lab="Diversity (log)")
abline(v=22, col="red")
par(op)





#Isolated the non-rarefaction disparity results
dis_ran_max_slater<-extract.disp(disparity_full_slater$quantiles, rarefaction="max")
dis_ran_min_slater<-extract.disp(disparity_full_slater$quantiles, rarefaction="min")
dis_ran_med_slater<-extract.disp(disparity_full_slater$quantiles, rarefaction=10)

op<-par(mfrow=c(3,1), bty="l")
plot.disparity(dis_ran_max, diversity=log(dis_ran_max$rarefaction), main="max")
abline(v=22, col="red")
plot.disparity(dis_ran_min, diversity=log(dis_ran_min$rarefaction), main="min")
abline(v=22, col="red")
plot.disparity(dis_ran_med, diversity=log(dis_ran_med$rarefaction), main="10")
abline(v=22, col="red")
par(op)

#Pretty plot with the tree and the disparity
op<-par(mfrow=c(2,1), bty="l")
geoscalePhylo(ladderize(tree), cex.age=0.6, cex.ts=0.8, cex.tip=0.5)
plot.disparity(dis_ran_max, diversity=log(dis_ran_max$rarefaction))
par(op)