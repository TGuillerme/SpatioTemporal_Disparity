#Paleo style (use CLADIS) per bin

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

#Load the data
source("Test.data.R")


######################
#Creating the matrix with ACE
######################


#Renaming the matrix to match with Graeme's workshop
#Choose one of the following matrices:
#nexus.data<-euarch.nex #Tips and nodes
#nexus.data<-STD.nex #Tips and nodes STD method
nexus.data<-CLADDIS.nex #Tips and nodes CLADDIS method
#Problem with CLADDIS method: no account for uncertainty


######################
#Running the PCO
######################


#Safe taxonomic reduction
#safe.data <- SafeTaxonomicReduction(nexus.data) #Removes nodes

#Distance matrix
dist.data <- MorphDistMatrix(nexus.data)

#Trim data
trimmed.max.data <-TrimMorphDistMatrix(dist.data$max.dist.matrix)

#Run the PCO
pco.data <- cmdscale(trimmed.max.data$dist.matrix, k=nrow(trimmed.max.data$dist.matrix) - 1, add=T)$points


######################
#Plotting the tree
######################


#Renaming the tree
tree.data<-euarch.tree

#Tree ages (useless?)
ages.data<-tree.age(tree.data)
tree.data$root.time<-max(ages.data[,1])
#FAD/LAD
ages.data<-data.frame("FAD"=tree.age(tree.data)[1:Ntip(tree),1], "LAD"=tree.age(tree.data)[1:Ntip(tree.data),1], row.names=tree.age(tree.data)[1:Ntip(tree.data),2])

#Plot the tree
geoscalePhylo(ladderize(tree.data), cex.age=0.6, cex.ts=0.8, cex.tip=1)


######################
#Disparity
######################

#Selecting ALL the pco axis
pco.data

#Calculating the rarefaction
rarefaction_median<-disparity(pco.data, rarefaction=TRUE, verbose=TRUE, central_tendency=median)
rarefaction_mean<-disparity(pco.data, rarefaction=TRUE, verbose=TRUE, central_tendency=mean)
rarefaction_median_Nth<-disparity(pco.data, rarefaction=TRUE, verbose=TRUE, central_tendency=median)
rarefaction_mean_Nth<-disparity(pco.data, rarefaction=TRUE, verbose=TRUE, central_tendency=mean)


dev.new()
Disparity.Output
rarefaction<-rarefaction_mean_Nth
op<-par(mfrow=c(3,2))
plot.disparity(rarefaction, measure="Cent.dist", rarefaction=TRUE, xlab="Taxa", ylab="Distance from centroid")
plot.disparity(rarefaction, measure="Sum.range", rarefaction=TRUE, xlab="Taxa", ylab="Sum of ranges")
plot.disparity(rarefaction, measure="Prod.range", rarefaction=TRUE, xlab="Taxa", ylab="Product of ranges")
plot.disparity(rarefaction, measure="Sum.var", rarefaction=TRUE, xlab="Taxa", ylab="Sum of variance")
plot.disparity(rarefaction, measure="Prod.var", rarefaction=TRUE, xlab="Taxa", ylab="Product of variance")
par(op)

#Making bins
bins_breaks<-c(80,60,40,20,0)
pco_binned<-bin.pco(pco.data, tree.data, bins_breaks, include.nodes=FALSE)

disparity_binned<-lapply(pco_binned, disparity, verbose=TRUE)



ChangesInBins   
EdgeLengthsInBins

#Slicing the tree
#set the slices
slices<-seq(from=40, to=80, by=5)
#set the data as pco.scores
pco.scores<-list(pco.data)
names(pco.scores)<-"scores" ; class(pco.scores) <- "pco.scores"
#slicing
std.slice_pro<-std.slice(tree.data, pco.scores, slices=slices, method="PROXIMITY")
#Plot (visual)
#plot.std(std.slice_pro, legend=FALSE, pars=c(3,4))


#Calculate the disparity per time slice
disparity_slices<-list()
for (slice in 1:length(slices)){
    disparity_slices[[slice]]<-disparity(std.slice_pro[[slice]]$sub_scores, rarefaction=FALSE)
}
names(disparity_slices)<-slices
#Transform the results as a matrix
disparity_slices_table<-matrix(ncol=ncol(disparity_slices[[1]]), data=unlist(disparity_slices), byrow=TRUE)
disparity_slices_table<-as.data.frame(disparity_slices_table)
colnames(disparity_slices_table)<-names(disparity_slices[[1]])

#Plotting the disparity through time
plot.disparity(disparity_slices_table, measure="Prod.var", rarefaction=FALSE, xlab="Mya", ylab="Disparity (Distance from centroid)")
op<-par(mfrow=c(3,2))
plot.disparity(disparity_slices_table, measure="Cent.dist", rarefaction=FALSE, xlab="Time since present", ylab="Disparity (Distance from centroid)")
plot.disparity(disparity_slices_table, measure="Sum.range", rarefaction=FALSE, xlab="Time since present", ylab="Disparity (Sum of ranges)")
plot.disparity(disparity_slices_table, measure="Prod.range", rarefaction=FALSE, xlab="Time since present", ylab="Disparity (Product of ranges)")
plot.disparity(disparity_slices_table, measure="Sum.var", rarefaction=FALSE, xlab="Time since present", ylab="Disparity (Sum of variance)")
plot.disparity(disparity_slices_table, measure="Prod.var", rarefaction=FALSE, xlab="Time since present", ylab="Disparity (Product of variance)")
par(op)




#Calculating the PCO/MDS with a scaled euclidean distance matrix and removing the NAs
pco<-pco.std(submatrix, distance="euclidean", scale=TRUE, center=FALSE, na.rm=TRUE, correction="none")

#Visualising the axis variance load
plot.std(pco, legend=TRUE)

#Creating the pco.scores object (containing the axis and the taxonomy)
#Change the taxonomy list into group list containing not only taxonomic data (BM, etc...).
pco.scores<-as.pco.scores(tree, pco, n.axis=2, taxonomy.list)

#Full pco plot
plot.std(pco.scores, legend=TRUE, main="Full character-space")

#Creating the slice list
std.slice_acc<-std.slice(tree.data, pco.scores, slices, method="ACCTRAN")
std.slice_del<-std.slice(tree, pco.scores, slices, method="DELTRAN")
std.slice_pro<-std.slice(tree, pco.scores, slices, method="PROXIMITY")

#Plot the two series of slices
plot.std(std.slice_acc, legend=TRUE, pars=c(3,3))
plot.std(std.slice_del, legend=TRUE, pars=c(3,3))
plot.std(std.slice_pro, legend=TRUE, pars=c(3,3))