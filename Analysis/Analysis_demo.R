# This is a tidied up script going through each step of the disparity analysis on a subset of data for speeding up computational
# time.
# Please email me any questions or query (guillert@tcd.ie).

# This demo is divided in six parts:
# A - Preparing the data 
# B - Estimating the ancestral characters states
# C - Computing the distance matrix
# D - Measuring disparity through time
# E - Visualising disparity through time
# F - Measuring the effect of the K-Pg boundary
# The data necessary for all this analysis is located in the data folder from this repository. If you are running this script from
# the Analysis folder, and that you downloaded all the data you should have no problems running this script. Otherwise it is
# necessary to download both the tree from Beck & Lee 2014 and the the cladistic matrix. Both can be found at the following links
# @@@@@@@@@ ADD LINKS @@@@@@@@@
# The lines you might to modify the path will contain the following tag
    # ~~ PATH ~~

# Before the analysis, we need to download all the functions and packages used in this demo. Most functions are based on Graeme
# Lloyd's Claddis package:
browseURL("http://cran.r-project.org/web/packages/Claddis/")
# and on a draft version of a future package that I will develop along with Graeme Lloyd. For now this package is not properly tested
# and neither is it documented so I will describe the few used functions as we use them. This package, temporarily called 'disparity'
# can be downloaded from my gitHub repository. Because the developement of this package is not finished, there is no manuals
# associated to each function. However, you can find details on each functions in each function script here:
browseURL("https://github.com/TGuillerme/SpatioTemporal_Disparity/tree/master/Functions/disparity/R")

# We can use devtools to install it.
library(devtools)
install_github("TGuillerme/SpatioTemporal_Disparity/Functions/disparity")
library(disparity)

# To speed up the analysis, we are going to use only Beck & Lee (2014; Proceedings B) data.

#########################
#
#   A - Preparing the data
#
#########################

# The tree
Tree_data<-read.nexus("../Data/2014-Beck-ProcB-TEM.tre") # ~~ PATH ~~

# The cladistic matrix
Nexus_data<-ReadMorphNexus("../Data/2014-Beck-ProcB-matrix-morpho.nex") # ~~ PATH ~~
# Isolating the matrix
Nexus_matrix<-Nexus_data$matrix

# We then need to clean both the tree and the matrix using the clean.tree and clean.table function so the data matches exactly the
#tree.
# Cleaning the tree
Tree_data<-clean.tree(Tree_data, Nexus_matrix)
# Cleaning the matrix
Nexus_data$matrix<-clean.table(Nexus_matrix, Tree_data)
# Making sure they match
nrow(Nexus_data$matrix) == Ntip(Tree_data)

# A last step for the tree is to add the node labels and to make sure the tree has a root time information. We can do that using the
# lapply.root function
tree<-lapply.root(Tree_data)

# We can now visualise the tree using the geoscalePhylo function from the 'strap' package
geoscalePhylo(ladderize(tree), cex.age=0.6, cex.ts=0.7, cex.tip=0.5, units=c("Period","Epoch"), boxes="Epoch")
# We can also add the K-Pg boundary
abline(v=106.3, col="red")


#########################
#
#   B - Estimating the ancestral states
#
#########################

# We can estimate the ancestral states using the anc.state function. This function intakes:
# -A tree
# -A nexus data object
# -A specification on the method used: here we are going to use the ML-claddis method that uses the rerooting method to reconstruct
# ancestral states and Lloyd's method for dealing with missing data (see manuscript for more details).

# Warning, this task can be fairly long because it's estimating the states of 421 characters for 101 nodes. Skip this step if you are
# not patient.
#anc_states<-anc.state(tree, Nexus_data, method='ML-claddis', verbose=TRUE)

# For speeding up this demonstration, we can directly load the previously calculated states.
anc_states<-get(load("../Data/Beck2014/Beck2014_ancestral_states.Rda")) # ~~ PATH ~~

# This object contains both the ancestral states for each nodes added to the cladistic matrix but also the marginal probability
# associated to each state.
tail(anc_states$state[,1:10])
tail(anc_states$prob[,1:10])

# Because the problems linked to estimating ancestral states, we can use the anc.unc function to set a threshold for the marginal
# likelihood under which estimations are replaced by missing data.
save_anc_states<-anc.unc(anc_states, 0.95, missing=NA)$state

# The cladistic matrix now contains missing data
tail(save_anc_states[,1:10])

# We can now create the estimated cladistic matrix (i.e. observed states for tips + estimated states for nodes).
estimated_cladistic_matrix<-Nexus_data
estimated_cladistic_matrix$matrix<-save_anc_states

#########################
#
#   C - Computing the distance matrix
#
#########################

# We can then use the full estimated cladistic matrix to compute the distance matrix between every tips and every nodes. To do so, we
# can use a slightly modified version from the 'MorphDistMatrix' function from 'Claddis' that as been modified to be verbose and is
# now called MorphDistMatrix.verbose:

# Warning this step also takes some times (calculating 203 pairwise distances) so skip this step if you are not patient.
#distance_matrix<-MorphDistMatrix.verbose(estimated_cladistic_matrix, verbose=TRUE)

# For speeding up this demonstration, we can directly load the previously calculated states.
distance_matrix<-get(load("../Data/Beck2014/Beck2014_distance-nodes95.Rda")) # ~~ PATH ~~

# The MorphDistMatrix produces 5 different matrices using 5 different distances. In this manuscript we are just going to use the gower
# distance matrix (i.e. the euclidean distance divided by the number of shared characters).
gower_distance_matrix<-distance_matrix$gower.dist.matrix


#########################
#
#   D - Measuring disparity through time
#
#########################

# Before measuring disparity we need to set up two time variables: the time of the slices (every 5 million years from 170 Mya to the
# present) and getting information on the fossil occurrence (FAD/LAD). The occurrence data has already been compiled in a .csv file.
slices=rev(seq(from=0, to=170, by=5))
FADLAD<-read.csv('../Data/Beck2014_FADLAD.csv', row.names=1) # ~~ PATH ~~

# The first step is to transform the distance matrix into a cladisto-space (i.e. ordination, e.g. pco) using the cmdscale function.
# Because our distance are not Euclidean we are using the Cailliez correction (option 'add=TRUE') and setting the dimensions to the
# number of tips + nodes -2
cladisto_space<-cmdscale(gower_distance_matrix, k=nrow(gower_distance_matrix) - 2, add=TRUE)$points

# We time-sample this cladisto-space using the 'slice.pco' function. This function needs:
# -The full cladisto-space
# -The tree
# -The slices
# -The slicing method (the evolutionary model)
# -The FAD/LAD data
# It can also take various options such as being verbose (verbose=TRUE) and counting the taxonomic diversity at each sub-sample
# (diversity=TRUE). Here we are going to use the gradual model (method="proximity") but other methods are available (see the function
# script header for info).
cladisto_space_sliced<-slice.pco(cladisto_space, tree, slices, method='proximity', FAD_LAD=FADLAD, verbose=TRUE, diversity=TRUE)

# We can then isolate the sliced cladisto_space and the taxonomic diversity
taxonomic_diversity<-cladisto_space_sliced$diversity
cladisto_space_sliced<-cladisto_space_sliced$pco_slices

# We can already visualise the taxonomic diversity
plot(taxonomic_diversity, type="l", xlab="sub-samples")

# We can now measure disparity using the time.disparity function that intakes the the sliced cladisto-space and the metric we want to
# use for calculating disparity. Note that this function also intakes many other options.
disparity_through_time<-time.disparity(cladisto_space_sliced, method="centroid", verbose=TRUE)

# This function takes a bit of time because it bootstraps each sub-sample of the data set.

#########################
#
#   E - Visualising disparity through time
#
#########################

# The time.disparity function outputs disparity in a table format with the median distance from centroid and the confidence intervals
# limits for each slice.
head(disparity_through_time)

# We can plot these results using the plot.disparity function.
op<-par(bty="n", mar=c(4,4,4,4))
plot.disparity(disparity_through_time, diversity=taxonomic_diversity, xlab="Time (Mya)", y2lab="Diversity",
    ylab="Median distance from centroid")
abline(v=22, col="red")
par(op)

###################################
#
#   F - Measuring the effect of the K-Pg boundary
#
###################################

# We can then measure whether the changes in disparity observed through time are significantly different than 0. We can use the
# permanova function (adonis) embeded in the disparity.test.time function that intakes the sliced cladisto-space
disparity.test.time(cladisto_space_sliced, method="euclidean", permutations=1000)

# Finally, if a positive effect of time is detected, we can then run a series of post-hoc t-tests to measure the difference between
# the last sub-sample of the Cretaceous and all the sub-samples of the Cenozoic by using the disparity.test function.
reftest_pro_slater<-disparity.test(cladisto_space_sliced[22:35], method="centroid", test="reference", bootstraps=1000)
reftest_pro_slater