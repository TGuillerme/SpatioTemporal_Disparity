#Load the functions and the packages
library(disparity)

##############################
#Build the Eutheria morphospace
##############################

chain_name='Beck2014'
data_path='../Data/'
file_matrix='../Data/2014-Beck-ProcB-matrix-morpho.nex'
file_tree='../Data/2014-Beck-ProcB-TEM.tre'
intervals=as.numeric(strsplit(c(noquote('170.300,168.300,166.100,163.500,157.300,152.100,145.000,139.800,132.900,129.400,125.000,113.000,100.500,93.900,89.800,86.300,83.600,72.100,66.000,61.600,59.200,56.000,47.800,41.300,38.000,33.900,28.100,23.030,23.030,20.440,15.970,13.820,11.620,7.246,5.333,0.000')), split=',')[[1]])
slices=as.numeric(strsplit(c(noquote('170,165,160,155,150,145,140,135,130,125,120,115,110,105,100,95,90,85,80,75,70,65,60,55,50,45,40,35,30,25,20,15,10,5,0')), split=',')[[1]])
FADLAD='../Data/Beck2014_FADLAD.csv'

#matrix
Nexus_data<-ReadMorphNexus(file_matrix)
Nexus_matrix<-Nexus_data$matrix
#tree
Tree_data<-read.nexus(file_tree)

#FAD/LAD
FADLAD<-read.csv(FADLAD, row.names=1)

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

#Adding node labels to the tree
tree<-lapply.root(tree, max(tree.age(tree)$age)) 

#load the distance matrix
load(paste(data_path, chain_name, '/', chain_name, '_distance-nodes95.Rda', sep='')) #dist_nodes95
trimmed_max_data_nodes95<-TrimMorphDistMatrix(dist_nodes95$gower.dist.matrix)
tree_nodes95<-drop.tip(tree, trimmed_max_data_nodes95$removed.taxa) ; tree_nodes95$root.time<-max(tree.age(tree_nodes95)[,1])
trimmed_max_data_nodes95$dist.matrix<-trimmed_max_data_nodes95$dist.matrix[c(tree_nodes95$tip.label, tree_nodes95$node.label),c(tree_nodes95$tip.label, tree_nodes95$node.label)]
#pco
Eutheria_morphospace<-cmdscale(trimmed_max_data_nodes95$dist.matrix, k=nrow(trimmed_max_data_nodes95$dist.matrix) - 2, add=T)$points


##############################
#Build the Mammaliaformes morphospace
##############################
chain_name='Slater2013'
data_path='../Data/'
file_matrix='../Data/2013-Slater-MEE-matrix-morpho.nex'
file_tree='../Data/2013-Slater-MEE-TEM.tre'
intervals=as.numeric(strsplit(c(noquote('170.300,168.300,166.100,163.500,157.300,152.100,145.000,139.800,132.900,129.400,125.000,113.000,100.500,93.900,89.800,86.300,83.600,72.100,66.000,61.600,59.200,56.000,47.800,41.300,38.000,33.900,28.100,23.030,23.030,20.440,15.970,13.820,11.620,7.246,5.333,0.000')), split=',')[[1]])
slices=as.numeric(strsplit(c(noquote('170,165,160,155,150,145,140,135,130,125,120,115,110,105,100,95,90,85,80,75,70,65,60,55,50,45,40,35,30,25,20,15,10,5,0')), split=',')[[1]])
FADLAD='../Data/Slater2013_FADLAD.csv'

#matrix
Nexus_data<-ReadMorphNexus(file_matrix)
Nexus_matrix<-Nexus_data$matrix
#tree
Tree_data<-read.nexus(file_tree)

#FAD/LAD
FADLAD<-read.csv(FADLAD, row.names=1)

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

#Adding node labels to the tree
tree<-lapply.root(tree, max(tree.age(tree)$age)) 

#load the distance matrix
load(paste(data_path, chain_name, '/', chain_name, '_distance-nodes95.Rda', sep='')) #dist_nodes95
trimmed_max_data_nodes95<-TrimMorphDistMatrix(dist_nodes95$gower.dist.matrix)
tree_nodes95<-drop.tip(tree, trimmed_max_data_nodes95$removed.taxa) ; tree_nodes95$root.time<-max(tree.age(tree_nodes95)[,1])
trimmed_max_data_nodes95$dist.matrix<-trimmed_max_data_nodes95$dist.matrix[c(tree_nodes95$tip.label, tree_nodes95$node.label),c(tree_nodes95$tip.label, tree_nodes95$node.label)]
#pco
Mammaliaformes_morphospace<-cmdscale(trimmed_max_data_nodes95$dist.matrix, k=nrow(trimmed_max_data_nodes95$dist.matrix) - 2, add=T)$points


#######################
# Morphospaces dimensions
#######################

#Eutheria
ncol(Eutheria_morphospace)

#Mammaliaformes
ncol(Mammaliaformes_morphospace)

#######################
# Morphospaces loadings
#######################
op<-par(mfrow=c(2,1), bty="n")
barplot(cumsum(apply(Eutheria_morphospace, 2, var) / sum(apply(Eutheria_morphospace, 2, var))),
    main="Eutheria", ylab="Cumulative relative variance", xlab=paste("Dimensions (", ncol(Eutheria_morphospace), ")", sep=""))
barplot(cumsum(apply(Mammaliaformes_morphospace, 2, var) / sum(apply(Mammaliaformes_morphospace, 2, var))),
    main="Mammaliaformes", ylab="Cumulative relative variance", xlab=paste("Dimensions (", ncol(Mammaliaformes_morphospace), ")", sep=""))
par(op)

#######################
# Morphospaces plots
#######################
fun.ord.plot<-function(morphospace) {
    op<-par(mfrow=c(3,3), bty="n", mar=c(4,4,2,2))

    sequence <- matrix(data=c(seq(from=1, to=17, by=2), seq(from=2, to=18, by=2)), nrow=2, byrow=TRUE)
    sequence <- unlist(apply(sequence, 2, list), recursive=FALSE)

    plot.fun <- function(sequence_element, morphospace) {
        plot(morphospace[,sequence_element],
            xlab=paste("Axis ", sequence_element[[1]], " (",
                round((apply(morphospace, 2, var) / sum(apply(morphospace, 2, var)))[sequence_element[[1]]]*100, digit=2)
            ,"%)", sep=""),
            ylab=paste("Axis ", sequence_element[[2]], " (",
                round((apply(morphospace, 2, var) / sum(apply(morphospace, 2, var)))[sequence_element[[2]]]*100, digit=2)
            ,"%)", sep="")
        )
    }

    lapply(sequence, plot.fun, morphospace)

    par(op)
}

#Eutheria
fun.ord.plot(Eutheria_morphospace)
#variance explained:
sum(as.vector(apply(Eutheria_morphospace, 2, var) / sum(apply(Eutheria_morphospace, 2, var)))[1:18])

#Mammaliaformes
fun.ord.plot(Mammaliaformes_morphospace)
#variance explained:
sum(as.vector(apply(Mammaliaformes_morphospace, 2, var) / sum(apply(Mammaliaformes_morphospace, 2, var)))[1:18])
