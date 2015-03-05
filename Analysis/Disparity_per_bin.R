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

######################
#Testing with Beck data
######################

#Reading in the data
source("Beck.data.R")

#Renaming the Beck nexus file to match with the rest of the workshop
nexus.data <- Beck.nex

#Include ancestral states
#matrix<-anc.matrix.save$state
#matrix<-ifelse(matrix == '?', NA, matrix)
#nexus.data$matrix<-matrix

#Safe Taxonomic Reduction (Wilkinson 1995; Systematic Biology).
#Removes taxa we know (under parsimony) can only fall out in particular place(s) in the tree.
#safe.data <- SafeTaxonomicReduction(nexus.data)
#save(safe.data, file="../Data/2014-Beck-reduced_tax_matrix2.Rda")
#load("../Data/2014-Beck-reduced_tax_matrix2.Rda")

#Distance matrix
#Detail on the metrics: http://www.slideshare.net/graemelloyd/new-methodologies-for-the-use-of-cladistictype-matrices-to-measure-morphological-disparity-and-evolutionary-rate
#dist.data <- MorphDistMatrix(nexus.data)
#save(dist.data, file="../Data/2014-Beck-dist_matrices2.Rda")
load("../Data/2014-Beck-dist_matrices2.Rda") #dist.data

#Performs MDS on the GED matrix.
#warning("MDS on the GED matrix which 'cheats' by filling in those missing character distances.")
#cmdscale(dist.data$GED.dist.matrix)

#Remove the unaplicable characters
trimmed.max.data <-TrimMorphDistMatrix(dist.data$max.dist.matrix)

# We can see what taxa have been removed by typing:
trimmed.max.data$removed.taxa
# Remove the droped taxa from the tree
tree<-drop.tip(tree, trimmed.max.data$removed.taxa)

#Check gaps in the matrix
any(is.na(trimmed.max.data$dist.matrix))

#PCO

#Performs MDS on the MOD matrix
#cmdscale(trimmed.max.data$dist.matrix)

# We can maximise our axes by upping the value "k" (an option in the function) to N - 1 (the maximum number of axes for N objects, i.e., N taxa).
# In addition we want to use another option in the function (add) which gets around the negative eigenvalue problem that can cause downstream problems (e.g., a scree plot with negative values).
# We can specify these options fairly easily and store our answer in a new variable (pco.data) and this time we will just part of the output ($points) which are the values for our taxa on every ordination axis:
pco.data <- cmdscale(trimmed.max.data$dist.matrix, k=nrow(trimmed.max.data$dist.matrix) - 1, add=T)$points


######################
#Plot the tree
######################

#Renaming the tree
tree.data<-tree

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

#Calculating the rarefaction
#rarefaction_median<-disparity(pco.data, rarefaction=TRUE, verbose=TRUE, central_tendency=median)

#Making bins
bins_breaks<-rev(hist(ages.data[,1])$breaks)+5
bins_breaks[9]<-0
pco_binned<-bin.pco(pco.data, tree.data, bins_breaks, include.nodes=TRUE)
#Calculating the disparity per bins
disparity_binned_table<-bin.disparity(pco_binned, verbose=TRUE)

op<-par(mfrow=c(3,2))
plot.disparity(disparity_binned_table, rarefaction=FALSE, xlab="bins (Mya)", ylab="Distance from centroid", measure="Cent.dist")
abline(v=c(5,6), col="red")
plot.disparity(disparity_binned_table, rarefaction=FALSE, xlab="bins (Mya)", ylab="Sum of ranges", measure="Sum.range")
abline(v=c(5,6), col="red")
plot.disparity(disparity_binned_table, rarefaction=FALSE, xlab="bins (Mya)", ylab="Sum of variance", measure="Sum.var")
abline(v=c(5,6), col="red")
plot.disparity(disparity_binned_table, rarefaction=FALSE, xlab="bins (Mya)", ylab="Product of ranges", measure="Prod.range")
abline(v=c(5,6), col="red")
plot.disparity(disparity_binned_table, rarefaction=FALSE, xlab="bins (Mya)", ylab="Product of variance", measure="Prod.var")
abline(v=c(5,6), col="red")
par(op)