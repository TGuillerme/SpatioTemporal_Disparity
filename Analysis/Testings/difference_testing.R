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
#browseURL("https://github.com/TGuillerme/SpatioTemporal_Disparity/blob/master/Functions/disparity/R/MorphDistMatrix.verbose.R")

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
pco <- cmdscale(trimmed.max.data$dist.matrix, k=nrow(trimmed.max.data$dist.matrix) - 2, add=T, eig=TRUE)
pco <- pco$points

#bins
bins<-c(160, 140, 120, 100, 80, 60)

# We can also make it more accurate by adding the FAD-LAD data
pco_int<-int.pco(pco, tree.data, bins, FAD_LAD=ages.data)
pco_div<-cor.diversity(pco_int$diversity)

#Disparity
obs<-disparity(pco, bootstraps=5, save.all=TRUE, rarefaction=TRUE)
bs<-disparity(pco, method="centroid", bootstraps=1000, save.all=TRUE, centroid.type="full")

apply(bs$centroid, 2, median)
obs$centroid

obs$table[1,2]
bs$table[1,2]


disparity_int<-time.disparity(pco_int, method="centroid", verbose=TRUE, rarefaction=TRUE, save.all=TRUE, bootstraps=10, centroid.type="full")
#PROBLEM WITH COMBINING SLICES (TAXA DUPLICATION)
#PROBLEM WITH save.all=TRUE
#PROBLEM WITH rarefaction=TRUE



#Test comparing the distribution of medians of the centroid distances (BOOSTRAPED)
#Test comparing the distribution the centroid distances (NON-BOOSTRAPED)