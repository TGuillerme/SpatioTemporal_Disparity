# This script was written for the workshop at: http://www.linnean.org/Meetings-and-Events/Events/Radiation+and+Extinction+-+Investigating+Clade+Dynamics+in+Deep+Time
# This script is available to download from: http://www.graemetlloyd.com/teaching/RE2014/disparity_and_rates.r

# Install the packages on which Claddis depends from CRAN (it will not work without these!):
install.packages(c("phytools", "strap", "ape", "gdata"), dependencies=T)

# Install the paleotree package from CRAN:
install.packages("paleotree", dependencies=T)

# Install the strap package from CRAN:
install.packages("strap", dependencies=T)

# Install the devtools package from CRAN:
install.packages("devtools")

# Load the devtools package into R:
library(devtools)

# Install the Claddis package from GitHub:
install_github("graemetlloyd/Claddis")

# The lines above here only need to be run the first time you use this script

# Load the Claddis package into R:
library(Claddis)

# Load the paleotree library into R (we will use this for time-scaling our tree later):
library(paleotree)

# Load the strap library into R (we will use this for plotting a time-scaled tree):
library(strap)

# By typing a question mark and the package name we can bring up some basic details:
?Claddis

# Try clicking on the 'index' link for a list of functions then on a function name for a basic help file (these are in various states of development!)

# Before we can use any of the other functions we really need to read in some data.
# To do this we will use ReadMorphNexus, a function designed to read in morphological data in the common #NEXUS format.
# Note that for today we will read in a fairly small data set as many of the methods will run too slowly on a larger file to be useful in a 45-minute session!
# HOWEVER, by all means try to at least load your own data with this function as this is usually the hardest step.
# Here we are going to use Cullen et al (2013):
browseURL("http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0058853")

# This has a few advantages for us: 1. it is small, 2. it is open access, 3. it is on dinosaurs (where age data is easily available).
# First of all we can simply look at the data on line:
browseURL("http://www.graemetlloyd.com/nexus/Cullen_etal_2013a.nex")

# You can download this file and read it in to R from your own hard drive, but for now we will read it in direct from the web as this means a single address (should!) work for everyone:
nexus.data <- ReadMorphNexus("http://www.graemetlloyd.com/nexus/Cullen_etal_2013a.nex")

# In the above line we have stored (<-) the data in a variable (nexus.data).
# This data is stored in an R format known as a "list", which is (one) way to store different types of data together.
# The components of this list (in this case) have names that can be seen by typing:
names(nexus.data)

# We can view individual components of this list by using the dollar symbol and a name from above, e.g.:
nexus.data$matrix

# This shows just the matrix itself.
# Note that all types of missing data (including gaps/inapplicables) are replaced with NA and the data are stored as strings (programming-speak for text).
# This might seem odd given that the data are coded as numbers, but is necessary due to polymorphisms (I am not going to talk about these).
# This is someting you should bear in mind if you want to run code on the data (i.e., mathematical operators will not work unless you alter the data first!).
# We can quickly look at an example of a polymorphism: character 22 for the taxon "Sinornithomimus_dongi":
nexus.data$matrix["Sinornithomimus_dongi", 22]

# Now we have data can try using some of the other functions.
# One function that we can run on just a morphological matrix (with no additional data) is Safe Taxonomic Reduction (Wilkinson 1995; Systematic Biology).
# This is a simple way to remove taxa we know (under parsimony) can only fall out in particular place(s) in the tree.
# We can thus "safely" remove them prior to inference and save ourselves some computation time.
# Here we will run it and store the results in a new variable (safe.data):
safe.data <- SafeTaxonomicReduction(nexus.data)

# This is also a list. For now we can just look at the first part of it:
safe.data$str.list

# You should see a matrix with three columns (Junior, Senior, and Rule).
# Here "Junior" is a taxon that can be safely removed, "Senior" is the taxon (or taxa) it will be expected to fall out next to and "Rule" is the rule (see Willinson 1995) under which it can be safely removed.
# In this case there is only one taxon ("Dry_Island_bonebed_specimens") that can be removed.
# Note that if you run this on your own data you may not be able to remove any taxa.
# If this occurs instead of a list you will simply get a string (text) telling you "No taxa can be safely removed".

# Another function we can use on just a matrix is MorphDistMatrix, which converts a cladistic matrix into a distance matrix using the various diferent metrics I mentioned in my talk: http://www.slideshare.net/graemelloyd/new-methodologies-for-the-use-of-cladistictype-matrices-to-measure-morphological-disparity-and-evolutionary-rate
# Again we will use a new variable (dist.data) to store the output:
dist.data <- MorphDistMatrix(nexus.data)

# This is (again!) a list. We can see the names of each part again using names():
names(dist.data)

# These are the four metrics I discussed in my talk (max = MOD), plus an additional matrix giving the number of characters that are scored in BOTH taxa for each pairwise comparison.
# Lets check this last one quickly to see if there are any zeroes (i.e., incalculable distances):
any(dist.data$comp.char.matrix == 0)

# You should see it return a logical (TRUE or FALSE).
# In this case it is TRUE, which is a problem if we wanted to ordinate our data using anything other than GED.
# We can show this very easily using the cmdscale (Classic Multi-dimensional Scaling; i.e., principal coordinates) function and the raw distance matrix:
cmdscale(dist.data$raw.dist.matrix)

# You should see we get the following message: "Error in cmdscale(dist.data$raw.dist.matrix) : NA values not allowed in 'd'".
# Here "d" is our distance matrix.
# You should get the same error with both "$gower.dist.matrix" and "$max.dist.matrix", but not with "$GED.dist.matrix" (which "cheats" by filling in those missing character distances).
# We can take a quick look at it to see if we can spot where the problem is (this time using the MOD metric):
dist.data$max.dist.matrix

# You should see two NA values here, where we compare "Qiupalong_henanensis" with "Pelecanimimus_polyodon" and "Dry_Island_bonebed_specimens" with "Pelecanimimus_polyodon".
# We could thus remove both "Qiupalong_henanensis" and "Dry_Island_bonebed_specimens", or just "Pelecanimimus_polyodon".
# If we want to go down the latter route (which means we retain as many taxa as possible) we can use another function (TrimMorphDistMatrix) and another new variable (trimmed.max.data):
trimmed.max.data <-TrimMorphDistMatrix(dist.data$max.dist.matrix)

# This is (again!) a list.
# We can see what taxa have been removed by typing:
trimmed.max.data$removed.taxa

# This should just show us the name "Pelecanimimus_polyodon", which is good as that is what we were expecting!
# We now also have a distance matrix without any gaps (NAs):
any(is.na(trimmed.max.data$dist.matrix))

# This should be FALSE, which means we can hand it to cmdscale and not get an error:
cmdscale(trimmed.max.data$dist.matrix)

# This should give a two-column matrix, but this is not what we really want as forcing all the variance on to two axes will confer a cost in terms of fidelity to the true distance.
# We can maximise our axes by upping the value "k" (an option in the function) to N - 1 (the maximum number of axes for N objects, i.e., N taxa).
# In addition we want to use another option in the function (add) which gets around the negative eigenvalue problem that can cause downstream problems (e.g., a scree plot with negative values).
# We can specify these options fairly easily and store our answer in a new variable (pco.data) and this time we will just use part of the output ($points) which are the values for our taxa on every ordination axis:
pco.data <- cmdscale(trimmed.max.data$dist.matrix, k=nrow(trimmed.max.data$dist.matrix) - 1, add=T)$points

# Before we plot this lets get the data we need to make a scree plot:
scree.data <- apply(pco.data, 2, var) / sum(apply(pco.data, 2, var)) * 100

# We can make a simple plot of this:
plot(scree.data, type="l", xlab="Ordination axis", ylab="Percentage variance")

# There are a couple of things you should notice here:
# 1. The variance on the first two axes are pretty low.
# This should lead us to be cautious in drawing conclusions from any single bivariate plot of the data.
# 2. There are two "doglegs" in the data.
# The second one (bending downwards at the right of the plot) is related to that negative eigenvalue problem and our associated correction.
# In this case this means the only axis we can really ignore is the 13th one (everything else contains some variance we should probably care about).

# Before we start plotting a useful thing to do is define our plotting axes first so we can edit these later to easily plot different axes:
PCOx <- 1
PCOy <- 2

# We can use this line to plot our data:
plot(pco.data[, PCOx], pco.data[, PCOy], xlab=paste("PCO ", PCOx, " (", round(scree.data[PCOx], 2), "% variance)", sep=""), ylab=paste("PCO ", PCOy, " (", round(scree.data[PCOy], 2), "% variance)", sep=""), pch=19)

# This can be a bit tricky to interpret, not least of all because we do not know which point is which taxon.
# We can fix this by adding the taxon names:
text(pco.data[, PCOx], pco.data[, PCOy], rownames(pco.data))

# Unfortunately this is not the nicest plot in the world, but for now we will move on.
# So far we have done everything with just a single input (a cladistic matrix), but adding additional data (a topology plus ages for our tips) opens up more analyses.
# To keep things simple we will use just a single tree (the first MPT) for our data.
# We can view the Newick string here...:
browseURL("http://www.graemetlloyd.com/firstmpt/Cullen_etal_2013a.tre")

# ...and read it in and store in a new variable (tree.data) using the ape function read.tree:
tree.data <- read.tree("http://www.graemetlloyd.com/firstmpt/Cullen_etal_2013a.tre")

# Note: if your own tree data is in #NEXUS format you should instead use read.nexus.

# There is a slight problem here:
plot(tree.data)

# There are a couple of polytomies in our tree, but many phylogenetic approaches expect a fully bifurcating tree.
# We could break our polytomies in some sensible way, but for now lets just do this at random using the ape function multi2di and overwrite our tree.data variable:
tree.data <- multi2di(tree.data)

# Double-checking we can now see there are no polytomies:
plot(tree.data)

# Our tree is useless without branch-lengths (time) and as this is not a time-scaled tree we need to import some age data before we can do this.
# We can take a look at how this data should be formatted here...:
browseURL("http://www.graemetlloyd.com/teaching/RE2014/Cullenages.txt")

# ...and read it in here (storing in a new variable, ages.data):
ages.data <- read.table("http://www.graemetlloyd.com/teaching/RE2014/Cullenages.txt", row.names=1, sep="\t", header=T)

# I could give an entire workshop on how to time-scale a tree of fossil taxa, but for now we will just use the paleotree package and the function timePaleoPhy.
# Here we are timescaling our tree treating our tip dates (FAD-LAD) as true ranges and forcing a minimum branch-length of 2 million years.
# Note that many approaches require all branch-lengths to be positive (no zeroes) so this at least deals with that problem.
# When we do this we may as well overwrite our tree.data variable again (R will not warn you when you do this so this is can be a bad habit to develop!):
tree.data <- timePaleoPhy(tree.data, ages.data, type="mbl", vartime=2)

# Note: an important additional variable for time-scaled trees of fossil taxa is the root age:
tree.data$root.time

# This can be important when plotting your tree against time as most phylogenetic software assumes the youngest tip terminates at the present.

# You can plot your tree in a much prettier way using Mark Bells awesome geoscalePhylo function in the strap package:
geoscalePhylo(ladderize(tree.data), cex.age=0.6, cex.ts=0.8, cex.tip=1)

# Now we have a time-scaled tree we can calculate rates.
# As I (briefly) mentioned this is hard to do for time series so please disregard that output from the function (hopefully I will get a good approach working in future though).
# At the moment the function still wants some time bins so in the below we will simply set six equally spaced time bins from the root to the youngest tip.
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

# We can now plot our tree with branches coloured by rate (black = non-significant rates, red = significantly high rates, blue = significantly low rates):
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
# We can create a vector of colours for plotting these nodes:
node.color <- rep("white", nrow(rate.data$node.results))
node.color[which(rate.data$node.results[, "ml.signif.hi.ib"] == 1)] <- "red"
node.color[which(rate.data$node.results[, "ml.signif.lo.ib"] == 1)] <- "blue"
node.color[which(is.na(rate.data$node.results[, "ml.signif.lo.ib"]))] <- NA

# Now we can plot our tree...:
geoscalePhylo(tree.data, cex.age=0.6, cex.ts=0.8, cex.tip=1)

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

# Now we can add in the values for our tips to get vectors for both our x and our y:
all.PCOx <- c(pco.data[match(plot.tree$tip.label, rownames(pco.data)), PCOx], PCOx.anc)
all.PCOy <- c(pco.data[match(plot.tree$tip.label, rownames(pco.data)), PCOy], PCOy.anc)

# These can be modified further into two matrices that give the x and y coordinates needed to plot each branch:
branch.xs <- cbind(all.PCOx[plot.tree$edge[, 1]], all.PCOx[plot.tree$edge[, 2]])
branch.ys <- cbind(all.PCOy[plot.tree$edge[, 1]], all.PCOy[plot.tree$edge[, 2]])

# Make an empty plot (this allows us to plot our branches before our tips so the former do not appear on top of the latter):
plot(pco.data[, PCOx], pco.data[, PCOy], xlab=paste("PCO ", PCOx, " (", round(scree.data[PCOx], 2), "% variance)", sep=""), ylab=paste("PCO ", PCOy, " (", round(scree.data[PCOy], 2), "% variance)", sep=""), type="n")

# Now we can go through our branches one by one and plot them in grey...:
for(i in 1:nrow(branch.xs)) lines(x=branch.xs[i,], y=branch.ys[i,], col="grey", lwd=2)

# ...add our tips as black circles...:
points(pco.data[, PCOx], pco.data[, PCOy], pch=19)

# ...and finally add in our taxon names:
text(pco.data[, PCOx], pco.data[, PCOy], rownames(pco.data), cex=0.6)

# That is as far as we are going to go today, but there are some other functions in the package (and more are on the way):
# 
# AncStateEstMatrix - another way to get ancestral values that can be used in ordiantion spaces (this is also used by the discrete rates function)
# FindAncestor - a way to get the internal node for the MRCA of two or more tips
# GetNodeAges - a way to get all the node (tip and internal) dates for a time-scaled tree
# 
# Please feel free to post pull requests on GitHub (https://github.com/graemetlloyd/Claddis) or email me with queries: graemetlloyd@gmail.com
