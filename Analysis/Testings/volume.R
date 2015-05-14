# Claddis suggestions

#############################
# Function for loading raw scripts from GitHub 
#############################
# Just loading the whole bunch of functions. I you want to look at them in detail, I've added the GitHub links for each individual function as I call them through this script
library(devtools)
install_github("TGuillerme/SpatioTemporal_Disparity/Functions/disparity")
library(disparity)

#############################
# Data loading
#############################
# This is just exactly your workshop script up until the distance matrix calculations (without the "install.packages")

# Load the morphological matrix
nexus.data <- ReadMorphNexus("http://www.graemetlloyd.com/nexus/Cullen_etal_2013a.nex")

# Load the tree
tree.data <- read.tree("http://www.graemetlloyd.com/firstmpt/Cullen_etal_2013a.tre")

# Fix polytomies.
tree.data <- multi2di(tree.data)

# Read the age data
ages.data <- read.table("http://www.graemetlloyd.com/teaching/RE2014/Cullenages.txt", row.names=1, sep="\t", header=T)

# Building the tree
tree.data <- timePaleoPhy(tree.data, ages.data, type="mbl", vartime=2)

# Safe taxonomic reduction
safe.data <- SafeTaxonomicReduction(nexus.data)

#############################
# Distance matrix
#############################

# Distance matrix
dist.data <- MorphDistMatrix(nexus.data)

# For small datasets this is instantly but I'm playing around with reasonably big ones (100 taxa * 400 characters) and some can take a day or so to compute.
# To make sure the script is not frozen, I've just added a verbose function.
# I clumsily renamed it MorphDistMatrix.verbose to avoid any confusion.
browseURL("https://github.com/TGuillerme/SpatioTemporal_Disparity/blob/master/Functions/disparity/R/MorphDistMatrix.verbose.R")

# The function is exactly the same but can take an optional additional argument: verbose (TRUE (default) or FALSE)
# With verbose = FALSE, the function is exactly the same
dist.data.verbose <- MorphDistMatrix.verbose(nexus.data, verbose=FALSE)

#The verbose option just adds some dots and some text (probably worth modifying but that's the general idea)
dist.data.verbose <- MorphDistMatrix.verbose(nexus.data, verbose=TRUE)

# Just for the sake of testing the differences
for (different_distances in 1:length(dist.data)) {
    print( all( (dist.data.verbose[[different_distances]] == dist.data[[different_distances]]), na.rm=TRUE) )
}

# Let's just continue through your workshop to get to the disparity part

# Trim inapplicable data
trimmed.max.data <-TrimMorphDistMatrix(dist.data$max.dist.matrix)

# Remove the trimmed taxa from the tree for later
tree.data<-drop.tip(tree.data, trimmed.max.data$removed.taxa)

# PCO on MOD distance
pco <- cmdscale(trimmed.max.data$dist.matrix, k=nrow(trimmed.max.data$dist.matrix) - 1, add=T, eig=TRUE)

#########################
# ABOUT THE PC
#########################

# Difference between pcoa and cmdscale??
cmd <- cmdscale(trimmed.max.data$dist.matrix, k=9, add=F, eig=TRUE)
pcoa<- pcoa(trimmed.max.data$dist.matrix)

round(cmd$eig, digit=10) == round(pcoa$values$Eigenvalues, digit=10)
as.matrix(round(cmd$points, digit=10)) == as.matrix(round(pcoa$vectors, digit=10)) #weirdly the 7th dimension is reversed???

#In cmdscale the add=T uses a Cailliez (1983) correction correcting for negative eigenvalues. 
#To maximise the variance in our analysis, we create a space of n dimensions where n cam be represented by exactly in at most k-1 dimensions. Here we use k-2 dimensions because the two last eigenvalues are null (following Cailliez correction) 

cmd <- cmdscale(trimmed.max.data$dist.matrix, k=nrow(trimmed.max.data$dist.matrix) - 2, add=T, eig=TRUE)
pcoa<- pcoa(trimmed.max.data$dist.matrix, correction="cailliez")
#Cailliez correction applied to negative eigenvalues: D' = -0.5*(D + 1.64469512246817 )^2, except diagonal elements

all(round(cmd$eig, digit=10) == round(pcoa$values$Corr_eig, digit=10))
as.matrix(round(cmd$points[,-13], digit=10)) == as.matrix(round(pcoa$vectors.cor, digit=10)) #weirdly the 8th and 9th dimension are reversed???

#########################
# VOLUME
#########################

#Dummy test
clad.mat<-matrix(data=c(0,0,0,1,0,0,1,1,0,1,1,1), nrow=3) ; rownames(clad.mat) <- letters[1:3]
dist.mat<-as.matrix(dist(clad.mat, diag=TRUE))
pco<-cmdscale(dist.mat, eig=TRUE)
pco.data<-pco$points
pco.eig<-pco$eig
plot(pco.data)

round(pco.eig[-3], digit=10) == round(apply(cov(pco.data),2,sum)*2, digit=10)

#volume.fast(pco.data, pco.eig)
#volume(pco.data)
disparity(pco.data, "volume")

#ok for this simple case

#Dummy test 2
clad.mat<-matrix(data=sample(0:1,1000, replace=TRUE), nrow=20)
dist.mat<-as.matrix(dist(clad.mat, diag=TRUE))
pco<-cmdscale(dist.mat, k=nrow(dist.mat)-1, eig=TRUE, add=F)
pco.data<-pco$points
pco.eig<-pco$eig
plot(pco.data)

round(pco.eig[-20], digit=10) == round(apply(cov(pco.data),2,sum)*19, digit=10)

#volume.fast(pco.data, pco.eig)
#volume(pco.data)
disparity(pco.data, "volume")

#volume increases with the number of characters (increased distances)

# Do the proper Multidimensional Scaling (corrected for negative eigenvalues)
pco <- cmdscale(trimmed.max.data$dist.matrix, k=nrow(trimmed.max.data$dist.matrix) - 2, add=T, eig=TRUE)
pco.data <- pco$points #generates k-2 eigenvectors (=n cladistic space dimensions)
pco.eigen<- pco$eig #generates k eigenvalues (including two null eigenvalues)

#Calculate the hyper-dimensional volume (i.e. size of the cladistic-space)

#volume.fast(pco.data, pco.eigen)
#volume(pco.data)
disparity(pco.data, "volume")

#Relation between eigen value and eigen vectors

set.seed(0)
X<-abs(matrix(rnorm(9),3,3))
diag(X)<-0
X[upper.tri(X)]<-X[lower.tri(X)]
first<-eigen(cov(X))
second<- apply(cov(first$vectors),2, sum)
#eigen(cov(X))$values != cmdscale(X, eig=TRUE)$eig #problem here?

eig.val<-cmdscale(X, eig=TRUE)$eig
eig.vec<-cmdscale(X, eig=TRUE)$points

round(eig.val[-3], digit=10) == round(apply(cov(eig.vec),2,sum)*2, digit=10)

#Works only for a distance matrix and when vectors are reported as in the cmdscale/pcoa algorithm?


#############################
# Volume through time
############################# 

# Let's first make a vector of intervals boundaries let's make it very simple: 3 bins with at least three taxa.
bins<-rev(seq(from=0, to=100, by=20))
# Note that the bins must be in reverse chronological order (time since the present)

# Now we can separate the pco.data in different bins
int.pco(pco.data, tree.data, bins)

# We can also make it more accurate by adding the FAD-LAD data
int.pco(pco.data, tree.data, bins, FAD_LAD=ages.data)

# We can also count the diversity
pco_in_bins<-int.pco(pco.data, tree.data, bins, FAD_LAD=ages.data, diversity=TRUE)
# The new object is a list of the binned pco data
pco_in_bins$pco_intervals
# And of the taxonomic counts per bins
pco_in_bins$diversity


# This function just takes the binned pco data from int.pco.
# Let's just calculate an easy (the default options)
disparity_per_bin<-time.disparity(pco_in_bins$pco_intervals, relative=FALSE, verbose=TRUE, bootstraps=1000)

# Et voilà!
disparity_per_bin

#############################
# Plotting the results
############################# 

# Finally I'm working on a last series of functions on plotting the results, it's still in a more "drafty" part than the rest but this can give you an idea.
browseURL("https://github.com/TGuillerme/SpatioTemporal_Disparity/blob/master/Functions/disparity/R/plot.disparity.R")
# The idea is to make a plot S3 method that recognise the data you input and gives a nice plot with the default plotting option (xlim, main, col, etc...)
# The function takes the data to plot (one of the disparity tables), which disparity metric to plot and various graphical options.
# The disparity measure to plot can be set to default which will be the first column of the table ("Cent.dist" by default)
# Here's for plotting the rarefaction curve of the rarefaction analysis ran above for the centroid distance
plot.disparity(rarefaction_results, xlab="Number of taxa")

# Or for the for classic disparity metrics
op<-par(mfrow=c(2,2))
plot.disparity(rarefaction_results, measure="Sum.range", xlab="Number of taxa")
plot.disparity(rarefaction_results, measure="Prod.range", xlab="Number of taxa")
plot.disparity(rarefaction_results, measure="Sum.var", xlab="Number of taxa")
plot.disparity(rarefaction_results, measure="Prod.var", xlab="Number of taxa")
par(op)

# We can also plot the disparity through time (more exciting!) by giving the disparity_per_bin data.
# This should take the exact same options as before with the default measurement to be plot being the first one in the table.
plot.disparity(disparity_per_bin)

# Finally we can also add the diversity curve calculated by the int.pco function using the diversity option.
plot.disparity(disparity_per_bin, diversity=pco_in_bins$diversity)

# Or for the other disparity metrics
op<-par(mfrow=c(2,3))
plot.disparity(disparity_per_bin, measure="Volume", diversity=pco_in_bins$diversity, ylim=c(0, 500))
plot.disparity(disparity_per_bin, measure="Cent.dist", diversity=pco_in_bins$diversity)
plot.disparity(disparity_per_bin, measure="Sum.range", diversity=pco_in_bins$diversity)
plot.disparity(disparity_per_bin, measure="Prod.range", diversity=pco_in_bins$diversity)
plot.disparity(disparity_per_bin, measure="Sum.var", diversity=pco_in_bins$diversity)
plot.disparity(disparity_per_bin, measure="Prod.var", diversity=pco_in_bins$diversity)
par(op)

# I developed more function on how to look at the tree through time (using ancestral states, etc...) but here's the global idead