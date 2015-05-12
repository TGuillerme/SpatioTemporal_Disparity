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
pco.data <- cmdscale(trimmed.max.data$dist.matrix, k=nrow(trimmed.max.data$dist.matrix) - 1, add=T)$points

#############################
# Disparity
#############################

# So here is a series of function for calculating different measures of disparity. I think it might be nice to add that to the package.
browseURL("https://github.com/TGuillerme/SpatioTemporal_Disparity/blob/master/Functions/disparity/R/disparity.R")

# The function intakes multiple arguments that allows to calculate various metrics of disparity (sum of ranges, prod of ranges, etc...), bootstrap, do rarefaction curves and remove some pco axis
# The first argument must be a pco matrix (output from cmdscale).
# The second argument, method, are which disparity metrics to calculate, it can calculate any metric among these: "centroid", "sum.range", "product.range", "sum.variance", "product.variance"
# More than one can be calculate at the same time. (default = all of them)
# The third argument is the number of bootstraps to do run (default = 1000)
# The forth argument is the size of the confidence intervals to report for the bootstraps (default = 50 and 95)
# The fifth argument is the central tendency to report. It can be any function like median or mean or anything else the user can come with (default = median)
# Then there is a series of logical arguments:
# Whether to run a rarefaction analysis (i.e. run the bootstraps n-1 times with each time removing randomly one taxa)
# Whether to be verbose (I like it for big datasets ;))
# Whether to remove the last pco axis (TRUE removes every axis beyond 95% of the variance but any value can be given as a threshold by the user(e.g. rm.last.axis=75 for 75%))
# However, as you said, this is just a non justifiable way to remove outliers from the dataset.
# Finaly the last argument (save.all) is whether to save all the values of the bootstraps/rarefaction draws. 

# For example, to calculate the 4 "classic" (but not excellent statistically) measures of disparity with all the other default options:
disparity(pco.data, method=c("sum.range", "product.range", "sum.variance", "product.variance"))

# We can also change the confidence interval and the central tendency for all the disparity metrics:
disparity(pco.data, CI=75, central_tendency=mean)

# Or including rarefaction and being verbose:
rarefaction_results<-disparity(pco.data, rarefaction=TRUE, verbose=TRUE)
rarefaction_results

# As you can see, the output is a table of the central tendency and the confidence intervals.
# But you can also save all the values using save.all
all_the_values<-disparity(pco.data, method="centroid", bootstraps=10, save.all=TRUE)

# This outputs a list of one table as before...
all_the_values$table

# As well as a matrix of all the bootstraps for the disparity metric (here 10 bootstraps for 14 taxa)
all_the_values$centroid

#############################
# Disparity through time
############################# 

# Another series of functions is to look at the disparity through time.
# I came up with various functions on how to look at diversity through time but here's the simplest one.
browseURL("https://github.com/TGuillerme/SpatioTemporal_Disparity/blob/master/Functions/disparity/R/int.pco.R")

# It just classically put species in bins (time intervals) and makes a list of them using the pco.data
# It takes as first argument the pco.data calculated by the cmdscale function.
# Then it needs a tree, a list of intervals and an optional FAD-LAD table.
# Note that if some taxa (or none) are not present in the FAD-LAD table, it will consider that their age is just the branch length and that they don't span through time.
# Then two options: whether to include nodes (let's leave that one for now)
# And whether to calculate the taxonomic diversity per interval (counting the taxa)

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

# We can now calculate the disparity for each slice using a lapply wrapping function for the disparity function that uses the same arguments as the disparity function
browseURL("https://github.com/TGuillerme/SpatioTemporal_Disparity/blob/master/Functions/disparity/R/time.disparity.R")

# This function just takes the binned pco data from int.pco.
# Let's just calculate an easy (the default options)
disparity_per_bin<-time.disparity(pco_in_bins$pco_intervals)

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
op<-par(mfrow=c(2,2))
plot.disparity(disparity_per_bin, measure="Sum.range", diversity=pco_in_bins$diversity)
plot.disparity(disparity_per_bin, measure="Prod.range", diversity=pco_in_bins$diversity)
plot.disparity(disparity_per_bin, measure="Sum.var", diversity=pco_in_bins$diversity)
plot.disparity(disparity_per_bin, measure="Prod.var", diversity=pco_in_bins$diversity)
par(op)

# I developed more function on how to look at the tree through time (using ancestral states, etc...) but here's the global idead