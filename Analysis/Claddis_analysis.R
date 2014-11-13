#Header
if(grep("TGuillerme", getwd())) {
    setwd('~/PhD/Projects/SpatioTemporal_Disparity/Analysis')
} else {
    warning("You might have to change the directory!")
}
if(!grep("SpatioTemporal_Disparity/Analysis", getwd())) {
    stop("Wrong directory!\nThe current directory must be:\nSpatioTemporal_Disparity/Analysis/\nYou can clone the whole repository from:\nhttps://github.com/TGuillerme/SpatioTemporal_Disparity")
}

source("functions.R")


#This script is based on Graeme Lloyd's workshop from the Linnean Society Meeting on the 11/11/2014.
#Workshop info: http://www.linnean.org/Meetings-and-Events/Events/Radiation+and+Extinction+-+Investigating+Clade+Dynamics+in+Deep+Time
#Workshop script: http://www.graemetlloyd.com/teaching/RE2014/disparity_and_rates.r

#Reading in the data
source("Beck.data.R")

#Renaming the Beck nexus file to match with the rest of the workshop
nexus.data <- Beck.nex

#Safe Taxonomic Reduction (Wilkinson 1995; Systematic Biology).
#Removes taxa we know (under parsimony) can only fall out in particular place(s) in the tree.
#safe.data <- SafeTaxonomicReduction(nexus.data)
#save(safe.data, file="../Data/2014-Beck-reduced_tax_matrix.Rda")
load("../Data/2014-Beck-reduced_tax_matrix.Rda")

#Distance matrix
#Detail on the metrics: http://www.slideshare.net/graemelloyd/new-methodologies-for-the-use-of-cladistictype-matrices-to-measure-morphological-disparity-and-evolutionary-rate
#dist.data <- MorphDistMatrix(nexus.data)
#save(dist.data, file="../Data/2014-Beck-dist_matrices.Rda")
load("../Data/2014-Beck-dist_matrices.Rda")

#Performs MDS on the GED matrix.
#warning("MDS on the GED matrix which 'cheats' by filling in those missing character distances.")
#cmdscale(dist.data$GED.dist.matrix)

#Remove the unaplicable characters
trimmed.max.data <-TrimMorphDistMatrix(dist.data$max.dist.matrix)

# We can see what taxa have been removed by typing:
trimmed.max.data$removed.taxa

#Check gaps in the matrix
any(is.na(trimmed.max.data$dist.matrix))

#Performs MDS on the MOD matrix
#cmdscale(trimmed.max.data$dist.matrix)

# We can maximise our axes by upping the value "k" (an option in the function) to N - 1 (the maximum number of axes for N objects, i.e., N taxa).
# In addition we want to use another option in the function (add) which gets around the negative eigenvalue problem that can cause downstream problems (e.g., a scree plot with negative values).
# We can specify these options fairly easily and store our answer in a new variable (pco.data) and this time we will just part of the output ($points) which are the values for our taxa on every ordination axis:
pco.data <- cmdscale(trimmed.max.data$dist.matrix, k=nrow(trimmed.max.data$dist.matrix) - 1, add=T)$points

# Before we plot this lets get the data we need to make a scree plot:
scree.data <- apply(pco.data, 2, var) / sum(apply(pco.data, 2, var)) * 100

# We can make a simple plot of this:
plot(scree.data, type="l", xlab="Ordination axis", ylab="Percentage variance")

# Before we start plotting a useful thing to do is define our plotting axes first so we can edit these later to easily plot different axes:
PCOx <- 1
PCOy <- 2

# We can use this line to plot our data:
plot(pco.data[, PCOx], pco.data[, PCOy], xlab=paste("PCO ", PCOx, " (", round(scree.data[PCOx], 2), "% variance)", sep=""), ylab=paste("PCO ", PCOy, " (", round(scree.data[PCOy], 2), "% variance)", sep=""), pch=19)

# Add labels
text(pco.data[, PCOx], pco.data[, PCOy], rownames(pco.data))

# And read it in and store in a new variable (tree.data) using the ape function read.tree:
tree.data <- tree


# There is a slight problem here:
plot(tree.data)

# Building the taxa ages data frame (TO FIX)
age.tree<-tree.age(tree)
ages.data<-data.frame("FAD"=tree.age(tree)[1:Ntip(tree),1], "LAD"=tree.age(tree)[1:Ntip(tree),1], row.names=tree.age(tree)[1:Ntip(tree),2])

#Replaces 0s by 1 (for now) and LAD by -1
ages.data[which(ages.data[,1]==0),1:2]<-1
ages.data[,2]<-ages.data[,1]-1


# I could give an entire workshop on how to time-scale a tree of fossil taxa, but for now we will just use the paleotree package and the function timePaleoPhy.
# Here we are timescaling our tree treating our tip dates (FAD-LAD) as true ranges and forcing a minimum branch-length of 2 million years.
# Note that many approaches require all branch-lengths to be positive (no zeroes) so this at least deals with that problem.
# When we do this we may as well overwrite our tree.data variable again (R will not warn you when you do this so this is can be a bad habit!):
tree.data <- timePaleoPhy(tree.data, ages.data, type="mbl", vartime=1)

# This can be important when plotting your tree against time as most phylogenetic software assumes the youngest tip terminates at the present.

# You can plot your tree in a much prettier way using Mark Bells awesome geoscalePhylo function in the strap package:
tree.data$root.time<-171.786
geoscalePhylo(ladderize(tree.data), cex.age=0.6, cex.ts=0.8, cex.tip=0.5)

# Now we have a time-scaled tree we can calculate rates.
# As I (briefly) mentioned this is hard to do for time series so please disregard that output from the function (hopefully I will get a good approach working in future though).
# At the moment the funciton still wants some time bins so in the below we will simply set six equally spaced time bins from the root to the youngest tip.
# Another option in the function is to set the alpha value for the significance tests.
# Note that the function already accounts for multiple comparisons using the Benjamini and Hochberg 1995 JRSSB False Discovery Rate (FDR).
# We can run the rate tests using DiscreteCharacterRate and store the results in a new variable (rate.data):
rate.data <- DiscreteCharacterRate(tree.data, nexus.data, seq(tree.data$root.time, tree.data$root.time - max(diag(vcv(tree.data))), length.out=6), alpha=0.01)

# Again this is a list and the output is quite verbose, e.g....:
rate.data$branch.results

# ...but we can just look at the branches that are significantly high or low:
rate.data$branch.results[, c("ml.signif.hi", "ml.signif.lo")]

# In this case we just have one significantly high rate and no significantly low rates.

# However, a better way to look at the data is visually by colouring the branches.
# We can start by creating a vector of colours for our branches:
edge.color <- rep("black", nrow(tree.data$edge))
edge.color[which(rate.data$branch.results[, "ml.signif.hi"] == 1)] <- "red"
edge.color[which(rate.data$branch.results[, "ml.signif.lo"] == 1)] <- "blue"

# We can now plot our tree with branches coloured by rate (black = non-significant rates, red = significantly high rates, blue = significiantly low rates):
geoscalePhylo(ladderize(tree.data), cex.age=0.6, cex.ts=0.8, cex.tip=1, edge.color=edge.color[match(ladderize(tree.data)$edge[, 2], tree.data$edge[,2])])

# We might also want to do this for clades...:
rate.data$node.results

# ...but this time there are no signficant rates:
rate.data$node.results[, c("ml.signif.hi", "ml.signif.lo")]

# However, as mentioned in my slides we can also compare results for just internal or just terminal branches.
# Here we do have one significant result for the internal branches:
rate.data$node.results[, c("ml.signif.hi.ib", "ml.signif.lo.ib")]

# Note that there are NAs in here for clades that are "cherries", i.e., nodes that lead to just two tips.
# (These contain no internal branches and hence cannot have a value.)
# We can 
node.color <- rep("white", nrow(rate.data$node.results))
node.color[which(rate.data$node.results[, "ml.signif.hi.ib"] == 1)] <- "red"
node.color[which(rate.data$node.results[, "ml.signif.lo.ib"] == 1)] <- "blue"
node.color[which(is.na(rate.data$node.results[, "ml.signif.lo.ib"]))] <- NA

# Now wec an plot our tree...:
geoscalePhylo(tree.data, cex.age=0.6, cex.ts=0.8, cex.tip=0.5)

# ...and plot our node results on top:
nodelabels(node=rate.data$node.results[, "node"][!is.na(node.color)], pch=21, col="black", bg=node.color[!is.na(node.color)])

# Another thing we can do with our tree is plot it into our ordination space from earlier.
# As before we will define our axes first to make it easier later to change what we plot:
PCOx <- 1
PCOy <- 2

# Because we trimmed our matrix down we also need to remove taxa from our tree that are not in our PCO data.
# We can do this using drop.tip in ace and setdiff (in base R) and store our answer in a new variable (plot.tree):
plot.tree <- drop.tip(tree.data, setdiff(tree.data$tip.label, rownames(pco.data)))

# Now we have a tree we can use to get values for our ancestors.
# Here we will use the ace function in ape:
PCOx.anc <- ace(pco.data[, PCOx], plot.tree, type="continuous")$ace
PCOy.anc <- ace(pco.data[, PCOy], plot.tree, type="continuous")$ace

# Now we an add in the values for our tips to gte vectors for both our x and our y:
all.PCOx <- c(pco.data[match(plot.tree$tip.label, rownames(pco.data)), PCOx], PCOx.anc)
all.PCOy <- c(pco.data[match(plot.tree$tip.label, rownames(pco.data)), PCOy], PCOy.anc)

# These can be modified further into two matrices that give the x and y coordinates needed to plot each branch:
branch.xs <- cbind(all.PCOx[plot.tree$edge[, 1]], all.PCOx[plot.tree$edge[, 2]])
branch.ys <- cbind(all.PCOy[plot.tree$edge[, 1]], all.PCOy[plot.tree$edge[, 2]])

# Make an empty plot (this allows us to plot our branches before our tips so the former do not appear on top of the latter):
plot(pco.data[, PCOx], pco.data[, PCOy], xlab=paste("PCO ", PCOx, " (", round(scree.data[PCOx], 2), "% variance)", sep=""), ylab=paste("PCO ", PCOy, " (", round(scree.data[PCOy], 2), "% variance)", sep=""), type="n")

# Now we can go through our branches one by one and plot them in grey...:
for(i in 1:nrow(branch.xs)) lines(x=branch.xs[i,], y=branch.ys[i,], col="grey", lwd=2)

# ...and finally add our tips as black circles...:
points(pco.data[, PCOx], pco.data[, PCOy], pch=19)

# ...and add in our taxon names:
text(pco.data[, PCOx], pco.data[, PCOy], rownames(pco.data), cex=0.6)

# That is as far as we are going to go today, but there are some other functions in the package (and more are on the way):
# 
# AncStateEstMatrix - another way to get ancestral values that can be used in ordiantion spaces (this is also used by the discrete rates function)
# FindAncestor - a way to get the internal node for the MRCA of two or more tips
# GetNodeAges - a way to get all the node (tip and internal) dates for a time-scaled tree
# 
# Please feel free to post pull requests on GitHub (https://github.com/graemetlloyd/Claddis) or email me with queries: graemetlloyd@gmail.com
