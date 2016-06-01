#####
# Disparity calculation using the dispRity package
####

######
# Before starting
#####

##### Package

# Loading the package
# The latest version of devtools and ape
if(!require(devtools)) install.packages("devtools")
if(!require(ape)) install.packages("ape")

## The latest version of Claddis (on development on my own branch)
install_github("TGuillerme/Claddis")

## The latest version of dispRity 
install_github("TGuillerme/dispRity")


##### Data

## Loading the packages
library(Claddis) ; library(dispRity)

## Loading the discrete morphological matrix
matrix_Beck <- ReadMorphNexus("../Data/2014-Beck-ProcB-matrix-morpho.nex")
matrix_Hall <- ReadMorphNexus("../Data/2015-Halliday-LinSoc-matrix-morpho.nex")

## Loading the tree
tree_Beck <- read.nexus("../Data/2014-Beck-ProcB-TEM.tre")
tree_Hall <- read.nexus("../Data/2015-Halliday-LinSoc-tree.nex")

#### TO FIX PROPERLY!
#
## Arbitrarily solving 0 branch lengths in Halliday's tree by setting them to
## 1% of the minimum branch length
min_brlen <- min(tree_Hall$edge.length[-which(tree_Hall$edge.length == 0)])*0.1
tree_Hall$edge.length[which(tree_Hall$edge.length == 0)] <- min_brlen
#
#### TO FIX PROPERLY!

## Making sure the tree and the matrix are matching
clean.MorphNexus <- function(matrix, tree) {
    cleaned_data <- clean.data(matrix$matrix, tree)
    tree <- cleaned_data$tree
    matrix$matrix <- cleaned_data$data
    return(list(tree, matrix))
}

tmp <- clean.MorphNexus(matrix_Beck, tree_Beck)
tree_Beck <- tmp[[1]] ; matrix_Beck <- tmp[[2]]

tmp <- clean.MorphNexus(matrix_Hall, tree_Hall[[1]])
tree_Hall <- tmp[[1]] ; matrix_Hall <- tmp[[2]]


## Adding node labels to the tree
tree_Beck <- makeNodeLabel(tree_Beck, method = "number", prefix = "n")
tree_Hall <- makeNodeLabel(tree_Hall, method = "number", prefix = "n")

## Adding root.time to the tree
tree_Beck$root.time <- max(tree.age(tree_Beck)[,1])
tree_Hall$root.time <- max(tree.age(tree_Hall)[,1])

## loading the first/last occurrence data
FADLAD_Beck <- read.csv("../Data/Beck2014_FADLAD.csv", row.names = 1)
FADLAD_Hall <- read.csv("../Data/Halliday2015_FADLAD.csv", row.names = 1)



##### Ancestral states

## Estimating the ancestral states (with a uncertainty threshold of 0.95)
matrix_nodes_Beck <- AncStateEstMatrix(matrix_Beck, tree_Beck,
    estimate.allchars = TRUE, uncertainty.threshold = 0.95)
## Saving the matrix
save(matrix_nodes_Beck, file = "../Data/ancestral_states_beck2014.Rda")

## Estimating the ancestral states (with a uncertainty threshold of 0.95)
matrix_nodes_Hall <- AncStateEstMatrix(matrix_Hall, tree_Hall,
    estimate.allchars = TRUE, uncertainty.threshold = 0.95)
## Saving the matrix
save(matrix_nodes_Hall, file = "../Data/ancestral_states_halliday2015.Rda")


###### Loading ancestral nodes

## Loading the matrix
load("../Data/ancestral_states_beck2014.Rda")
load("../Data/ancestral_states_halliday2015.Rda")
## The matrix was assigned to the object matrix_nodes (see above).

## Adding the node names as row names
row.names(matrix_nodes_Beck) <- tree_Beck$node.label
row.names(matrix_nodes_Hall) <- tree_Hall$node.label

## Combining the tips and the nodes in the same matrix
matrix_Beck$matrix <- rbind(matrix_Beck$matrix, matrix_nodes_Beck) 
matrix_Hall$matrix <- rbind(matrix_Hall$matrix, matrix_nodes_Hall) 


##### Calculating the distance matrix

## Calculating the Gower distance
matrix_dist_Beck <- MorphDistMatrix.fast(matrix_Beck, distance = "Gower")
save(matrix_dist_Beck, file = "../Data/distmatrix_beck2014.Rda")
matrix_dist_Hall <- MorphDistMatrix.fast(matrix_Hall, distance = "Gower")
save(matrix_dist_Hall, file = "../Data/distmatrix_halliday2015.Rda")



####### Ordinating the matrices

## Loading the distance matrices
load("../Data/distmatrix_beck2014.Rda")
load("../Data/distmatrix_halliday2015.Rda")

#### TO FIX PROPERLY!
#
# Replacing NAs (incomparable taxa) by median distance (not excellent practice)
# But there's only 8 NAs in 120409 distances (0.006%) so it won't influence the
# results (I guess)
matrix_dist_Hall[which(is.na(matrix_dist_Hall))] <- mean(matrix_dist_Hall, na.rm = TRUE)
#
#### TO FIX PROPERLY!

## Ordinating the matrices
matrix_ord_Beck <- cmdscale(matrix_dist_Beck,
                            k = nrow(matrix_dist_Beck) - 2, add = T)$points
matrix_ord_Hall <- cmdscale(matrix_dist_Hall,
                            k = nrow(matrix_dist_Hall) - 2, add = T)$points



####### Slicing the data

## Creating the time slices
time_slices <- rev(seq(from = 0, to = 120, by = 5))

## Creating the time series
time_series_Beck <- time.series(matrix_ord_Beck, tree_Beck,
    method = "continuous", time = time_slices, model = "gradual",
    FADLAD = FADLAD_Beck, verbose = TRUE)
time_series_Hall <- time.series(matrix_ord_Hall, tree_Hall,
    method = "continuous", time = time_slices, model = "gradual",
    FADLAD = FADLAD_Hall, verbose = TRUE)



##### Bootstrapping the data

## Bootstrapping the time series
time_series_Beck <- boot.matrix(time_series_Beck, bootstrap = 1000)
time_series_Hall <- boot.matrix(time_series_Hall, bootstrap = 1000)




##### Calculating Disparity

## Calculating the distance from the centroids
dist_centroids_Beck <- dispRity(time_series_Beck, metric = centroids)
dist_centroids_Hall <- dispRity(time_series_Hall, metric = centroids)

## Calculating the distance from the centre of the space
dist_centre_Beck <- dispRity(time_series_Beck, metric = centroids,
                               centroid = rep(0, ncol(matrix_ord_Beck)))
dist_centre_Hall <- dispRity(time_series_Hall, metric = centroids,
                               centroid = rep(0, ncol(matrix_ord_Hall)))

## Calculating the medians of the distances from centroids/centre
med_dist_centroids_Beck <- dispRity(dist_centroids_Beck, metric = median)
med_dist_centroids_Hall <- dispRity(dist_centroids_Hall, metric = median)
med_dist_centre_Beck <- dispRity(dist_centre_Beck, metric = median)
med_dist_centre_Hall <- dispRity(dist_centre_Hall, metric = median)




###### Plotting disparity

## Graphical parameters
op <- par(mfrow = c(2,2), bty = "n")

## Plotting the distance from centroids
plot(med_dist_centroids_Beck, main = "Beck 2014",
    ylab = "Median distance from centroids", xlab = "")
abline(v = 12, col = "red")
plot(med_dist_centroids_Hall, main = "Halliday 2015", ylab = "", xlab = "")
abline(v = 12, col = "red")

## Plotting the distances from centre
plot(med_dist_centre_Beck, ylab = "Median distance from center",
    xlab = "Age (Mya)")
abline(v = 12, col = "red")
plot(med_dist_centre_Hall, ylab = "", xlab = "Age (Mya)")
abline(v = 12, col = "red")



###### Effect of K-Pg?

## Testing any change in group size (centroids) after K-Pg
test.dispRity(get.dispRity(dist_centroids_Beck, what = c(19:33)), test = t.test,
    concatenate = FALSE, comparisons = "referential", correction = "bonferroni")
test.dispRity(get.dispRity(dist_centroids_Hall, what = c(19:33)), test = t.test,
    concatenate = FALSE, comparisons = "referential", correction = "bonferroni")
## Testing any change in group position (centre) after K-Pg
test.dispRity(get.dispRity(dist_centre_Beck, what = c(19:33)), test = t.test,
    concatenate = FALSE, comparisons = "referential", correction = "bonferroni")
test.dispRity(get.dispRity(dist_centre_Hall, what = c(19:33)), test = t.test,
    concatenate = FALSE, comparisons = "referential", correction = "bonferroni")