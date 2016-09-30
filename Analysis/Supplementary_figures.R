###################################
## 1- Setting up the data
###################################

set.seed(123)
## Loading the packages
library(Claddis) ; library(dispRity)

## Loading the data
###################################

## Loading the discrete morphological matrix
matrix_Beck <- ReadMorphNexus("../Data/2014-Beck-ProcB-matrix-morpho.nex")
matrix_Hall <- ReadMorphNexus("../Data/2015-Halliday-LinSoc-matrix-morpho.nex")

## Loading the tree
tree_Beck <- read.nexus("../Data/2014-Beck-ProcB-TEM.tre")
tree_Hall <- read.nexus("../Data/2015-Halliday-LinSoc-tree.nex")[[1]]

## Arbitrarily solving 0 branch lengths in Halliday's tree by setting them to
## 1% of the minimum branch length
zero_brlen <- which(tree_Hall$edge.length == 0)
min_brlen <- min(tree_Hall$edge.length[-zero_brlen])*0.1
tree_Hall$edge.length[which(tree_Hall$edge.length == 0)] <- min_brlen

## Making sure the tree and the matrix are matching
clean.MorphNexus <- function(matrix, tree) {
    cleaned_data <- clean.data(matrix$matrix, tree)
    tree <- cleaned_data$tree
    matrix$matrix <- cleaned_data$data
    return(list(tree, matrix))
}

tmp <- clean.MorphNexus(matrix_Beck, tree_Beck)
tree_Beck <- tmp[[1]] ; matrix_Beck <- tmp[[2]]

tmp <- clean.MorphNexus(matrix_Hall, tree_Hall)
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

## loading the ancestral states reconstructions
###################################

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

## Ordinating the data
###################################

## Loading the distance matrices
load("../Data/distmatrix_beck2014.Rda")
load("../Data/distmatrix_halliday2015.Rda")

#### TO FIX PROPERLY!
#
# Replacing NAs (incomparable taxa) by mean distance (not excellent practice)
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



###################################
## 2- Calculating disparity
###################################

## Creating the time slices/bins
###################################

## Creating equal time bins (5 Mya)
time_binsEQ <- rev(seq(from = 0, to = 120, by = 5))
## Creating stratigraphic bins (mean = 5 Mya)
time_binsST <- c(113.000,100.500,93.900,89.800,86.300,83.600,72.100,66.000,
    61.600,59.200,56.000,47.800,41.300,38.000,33.900,28.100,23.030,20.440,
    15.970,13.820,11.620,7.246,5.333,0.000)
## Creating the time slices
time_slices <- rev(seq(from = 0, to = 120, by = 5))

## Creating the time series
time_binsEQ_Beck <- time.series(matrix_ord_Beck, tree_Beck, method = "discrete", time = time_binsEQ, FADLAD = FADLAD_Beck, verbose = TRUE, inc.nodes = TRUE)
time_binsST_Beck <- time.series(matrix_ord_Beck, tree_Beck, method = "discrete", time = time_binsST, FADLAD = FADLAD_Beck, verbose = TRUE, inc.nodes = TRUE)
time_sliGRA_Beck <- time.series(matrix_ord_Beck, tree_Beck, method = "continuous", time = time_slices, FADLAD = FADLAD_Beck, model = "gradual", verbose = TRUE)
time_sliPUN_Beck <- time.series(matrix_ord_Beck, tree_Beck, method = "continuous", time = time_slices, FADLAD = FADLAD_Beck, model = "punctuated", verbose = TRUE)
time_sliACC_Beck <- time.series(matrix_ord_Beck, tree_Beck, method = "continuous", time = time_slices, FADLAD = FADLAD_Beck, model = "acctran", verbose = TRUE)
time_sliDEL_Beck <- time.series(matrix_ord_Beck, tree_Beck, method = "continuous", time = time_slices, FADLAD = FADLAD_Beck, model = "deltran", verbose = TRUE)

time_binsEQ_Hall <- time.series(matrix_ord_Hall, tree_Hall, method = "discrete", time = time_binsEQ, FADLAD = FADLAD_Hall, verbose = TRUE, inc.nodes = TRUE)
time_binsST_Hall <- time.series(matrix_ord_Hall, tree_Hall, method = "discrete", time = time_binsST, FADLAD = FADLAD_Hall, verbose = TRUE, inc.nodes = TRUE)
time_sliGRA_Hall <- time.series(matrix_ord_Hall, tree_Hall, method = "continuous", time = time_slices, FADLAD = FADLAD_Hall, model = "gradual", verbose = TRUE)
time_sliPUN_Hall <- time.series(matrix_ord_Hall, tree_Hall, method = "continuous", time = time_slices, FADLAD = FADLAD_Hall, model = "punctuated", verbose = TRUE)
time_sliACC_Hall <- time.series(matrix_ord_Hall, tree_Hall, method = "continuous", time = time_slices, FADLAD = FADLAD_Hall, model = "acctran", verbose = TRUE)
time_sliDEL_Hall <- time.series(matrix_ord_Hall, tree_Hall, method = "continuous", time = time_slices, FADLAD = FADLAD_Hall, model = "deltran", verbose = TRUE)

## Removing small intervals
time_binsEQ_Beck <- merge.time.series(time_binsEQ_Beck)
time_binsST_Beck <- merge.time.series(time_binsST_Beck)
time_binsEQ_Hall <- merge.time.series(time_binsEQ_Hall)
time_binsST_Hall <- merge.time.series(time_binsST_Hall)


## Bootstrapping the matrices
###################################
BS_time_binsEQ_Beck <- boot.matrix(time_binsEQ_Beck, 500)
BS_time_binsST_Beck <- boot.matrix(time_binsST_Beck, 500)
BS_time_sliGRA_Beck <- boot.matrix(time_sliGRA_Beck, 500)
BS_time_sliPUN_Beck <- boot.matrix(time_sliPUN_Beck, 500)
BS_time_sliACC_Beck <- boot.matrix(time_sliACC_Beck, 500)
BS_time_sliDEL_Beck <- boot.matrix(time_sliDEL_Beck, 500)
BS_time_binsEQ_Hall <- boot.matrix(time_binsEQ_Hall, 500)
BS_time_binsST_Hall <- boot.matrix(time_binsST_Hall, 500)
BS_time_sliGRA_Hall <- boot.matrix(time_sliGRA_Hall, 500)
BS_time_sliPUN_Hall <- boot.matrix(time_sliPUN_Hall, 500)
BS_time_sliACC_Hall <- boot.matrix(time_sliACC_Hall, 500)
BS_time_sliDEL_Hall <- boot.matrix(time_sliDEL_Hall, 500)

## Calculating disparity
###################################
# Median distance from centroids
med_dist_centroids_Beck_binsEQ <- dispRity(BS_time_binsEQ_Beck, metric = c(median, centroids))
med_dist_centroids_Beck_binsST <- dispRity(BS_time_binsST_Beck, metric = c(median, centroids))
med_dist_centroids_Beck_sliGRA <- dispRity(BS_time_sliGRA_Beck, metric = c(median, centroids))
med_dist_centroids_Beck_sliPUN <- dispRity(BS_time_sliPUN_Beck, metric = c(median, centroids))
med_dist_centroids_Beck_sliACC <- dispRity(BS_time_sliACC_Beck, metric = c(median, centroids))
med_dist_centroids_Beck_sliDEL <- dispRity(BS_time_sliDEL_Beck, metric = c(median, centroids))

# Median distance from centre
centre = rep(0, ncol(matrix_ord_Beck))
med_dist_centre_Beck_binsEQ <- dispRity(dispRity(BS_time_binsEQ_Beck, metric = centroids, centroid = centre), median)
med_dist_centre_Beck_binsST <- dispRity(dispRity(BS_time_binsST_Beck, metric = centroids, centroid = centre), median)
med_dist_centre_Beck_sliGRA <- dispRity(dispRity(BS_time_sliGRA_Beck, metric = centroids, centroid = centre), median)
med_dist_centre_Beck_sliPUN <- dispRity(dispRity(BS_time_sliPUN_Beck, metric = centroids, centroid = centre), median)
med_dist_centre_Beck_sliACC <- dispRity(dispRity(BS_time_sliACC_Beck, metric = centroids, centroid = centre), median)
med_dist_centre_Beck_sliDEL <- dispRity(dispRity(BS_time_sliDEL_Beck, metric = centroids, centroid = centre), median)

# sum ranges
sum_ranges_Beck_binsEQ <- dispRity(BS_time_binsEQ_Beck, metric = c(sum, ranges))
sum_ranges_Beck_binsST <- dispRity(BS_time_binsST_Beck, metric = c(sum, ranges))
sum_ranges_Beck_sliGRA <- dispRity(BS_time_sliGRA_Beck, metric = c(sum, ranges))
sum_ranges_Beck_sliPUN <- dispRity(BS_time_sliPUN_Beck, metric = c(sum, ranges))
sum_ranges_Beck_sliACC <- dispRity(BS_time_sliACC_Beck, metric = c(sum, ranges))
sum_ranges_Beck_sliDEL <- dispRity(BS_time_sliDEL_Beck, metric = c(sum, ranges))

# prod ranges
prod_ranges_Beck_binsEQ <- dispRity(BS_time_binsEQ_Beck, metric = c(prod, ranges))
prod_ranges_Beck_binsST <- dispRity(BS_time_binsST_Beck, metric = c(prod, ranges))
prod_ranges_Beck_sliGRA <- dispRity(BS_time_sliGRA_Beck, metric = c(prod, ranges))
prod_ranges_Beck_sliPUN <- dispRity(BS_time_sliPUN_Beck, metric = c(prod, ranges))
prod_ranges_Beck_sliACC <- dispRity(BS_time_sliACC_Beck, metric = c(prod, ranges))
prod_ranges_Beck_sliDEL <- dispRity(BS_time_sliDEL_Beck, metric = c(prod, ranges))

# sum variances
sum_variances_Beck_binsEQ <- dispRity(BS_time_binsEQ_Beck, metric = c(sum, variances))
sum_variances_Beck_binsST <- dispRity(BS_time_binsST_Beck, metric = c(sum, variances))
sum_variances_Beck_sliGRA <- dispRity(BS_time_sliGRA_Beck, metric = c(sum, variances))
sum_variances_Beck_sliPUN <- dispRity(BS_time_sliPUN_Beck, metric = c(sum, variances))
sum_variances_Beck_sliACC <- dispRity(BS_time_sliACC_Beck, metric = c(sum, variances))
sum_variances_Beck_sliDEL <- dispRity(BS_time_sliDEL_Beck, metric = c(sum, variances))

# prod variances
prod_variances_Beck_binsEQ <- dispRity(BS_time_binsEQ_Beck, metric = c(prod, variances))
prod_variances_Beck_binsST <- dispRity(BS_time_binsST_Beck, metric = c(prod, variances))
prod_variances_Beck_sliGRA <- dispRity(BS_time_sliGRA_Beck, metric = c(prod, variances))
prod_variances_Beck_sliPUN <- dispRity(BS_time_sliPUN_Beck, metric = c(prod, variances))
prod_variances_Beck_sliACC <- dispRity(BS_time_sliACC_Beck, metric = c(prod, variances))
prod_variances_Beck_sliDEL <- dispRity(BS_time_sliDEL_Beck, metric = c(prod, variances))


# Median distance from centroids
med_dist_centroids_Hall_binsEQ <- dispRity(BS_time_binsEQ_Hall, metric = c(median, centroids))
med_dist_centroids_Hall_binsST <- dispRity(BS_time_binsST_Hall, metric = c(median, centroids))
med_dist_centroids_Hall_sliGRA <- dispRity(BS_time_sliGRA_Hall, metric = c(median, centroids))
med_dist_centroids_Hall_sliPUN <- dispRity(BS_time_sliPUN_Hall, metric = c(median, centroids))
med_dist_centroids_Hall_sliACC <- dispRity(BS_time_sliACC_Hall, metric = c(median, centroids))
med_dist_centroids_Hall_sliDEL <- dispRity(BS_time_sliDEL_Hall, metric = c(median, centroids))

# Median distance from centre
centre = rep(0, ncol(matrix_ord_Hall))
med_dist_centre_Hall_binsEQ <- dispRity(dispRity(BS_time_binsEQ_Hall, metric = centroids, centroid = centre), median)
med_dist_centre_Hall_binsST <- dispRity(dispRity(BS_time_binsST_Hall, metric = centroids, centroid = centre), median)
med_dist_centre_Hall_sliGRA <- dispRity(dispRity(BS_time_sliGRA_Hall, metric = centroids, centroid = centre), median)
med_dist_centre_Hall_sliPUN <- dispRity(dispRity(BS_time_sliPUN_Hall, metric = centroids, centroid = centre), median)
med_dist_centre_Hall_sliACC <- dispRity(dispRity(BS_time_sliACC_Hall, metric = centroids, centroid = centre), median)
med_dist_centre_Hall_sliDEL <- dispRity(dispRity(BS_time_sliDEL_Hall, metric = centroids, centroid = centre), median)

# sum ranges
sum_ranges_Hall_binsEQ <- dispRity(BS_time_binsEQ_Hall, metric = c(sum, ranges))
sum_ranges_Hall_binsST <- dispRity(BS_time_binsST_Hall, metric = c(sum, ranges))
sum_ranges_Hall_sliGRA <- dispRity(BS_time_sliGRA_Hall, metric = c(sum, ranges))
sum_ranges_Hall_sliPUN <- dispRity(BS_time_sliPUN_Hall, metric = c(sum, ranges))
sum_ranges_Hall_sliACC <- dispRity(BS_time_sliACC_Hall, metric = c(sum, ranges))
sum_ranges_Hall_sliDEL <- dispRity(BS_time_sliDEL_Hall, metric = c(sum, ranges))

# prod ranges
prod_ranges_Hall_binsEQ <- dispRity(BS_time_binsEQ_Hall, metric = c(prod, ranges))
prod_ranges_Hall_binsST <- dispRity(BS_time_binsST_Hall, metric = c(prod, ranges))
prod_ranges_Hall_sliGRA <- dispRity(BS_time_sliGRA_Hall, metric = c(prod, ranges))
prod_ranges_Hall_sliPUN <- dispRity(BS_time_sliPUN_Hall, metric = c(prod, ranges))
prod_ranges_Hall_sliACC <- dispRity(BS_time_sliACC_Hall, metric = c(prod, ranges))
prod_ranges_Hall_sliDEL <- dispRity(BS_time_sliDEL_Hall, metric = c(prod, ranges))

# sum variances
sum_variances_Hall_binsEQ <- dispRity(BS_time_binsEQ_Hall, metric = c(sum, variances))
sum_variances_Hall_binsST <- dispRity(BS_time_binsST_Hall, metric = c(sum, variances))
sum_variances_Hall_sliGRA <- dispRity(BS_time_sliGRA_Hall, metric = c(sum, variances))
sum_variances_Hall_sliPUN <- dispRity(BS_time_sliPUN_Hall, metric = c(sum, variances))
sum_variances_Hall_sliACC <- dispRity(BS_time_sliACC_Hall, metric = c(sum, variances))
sum_variances_Hall_sliDEL <- dispRity(BS_time_sliDEL_Hall, metric = c(sum, variances))

# prod variances
prod_variances_Hall_binsEQ <- dispRity(BS_time_binsEQ_Hall, metric = c(prod, variances))
prod_variances_Hall_binsST <- dispRity(BS_time_binsST_Hall, metric = c(prod, variances))
prod_variances_Hall_sliGRA <- dispRity(BS_time_sliGRA_Hall, metric = c(prod, variances))
prod_variances_Hall_sliPUN <- dispRity(BS_time_sliPUN_Hall, metric = c(prod, variances))
prod_variances_Hall_sliACC <- dispRity(BS_time_sliACC_Hall, metric = c(prod, variances))
prod_variances_Hall_sliDEL <- dispRity(BS_time_sliDEL_Hall, metric = c(prod, variances))


###################################
## 3- Plotting!
###################################

#The following is a "BIG" plot comparing all methods/metrics for Slater
quartz(width = 15.6, height = 11.2) #A5 landscape
#Windows dimensions
op<-par(mfrow=c(6, 6), bty="n", mar=c(4,4,4,4))# oma=c(bottom, left, top, right)
#Centroid
plot.dispRity(med_dist_centroids_Beck_binsEQ, xlab="", ylab="Median distance from centroid", measure="Cent.dist", main="Intervals (equal)", elements=TRUE, cex.xaxis=0.7)
abline(v= 17.5, col="red")
plot.disparity(med_dist_centroids_Beck_binsST, xlab="", ylab="", measure="Cent.dist", main="Intervals (stratigraphy)", diversity=dis_nodes_slater$rarefaction, y2lab="", cex.xaxis=0.7)
abline(v= 17.5, col="red")
plot.disparity(med_dist_centroids_Beck_sliGRA, xlab="", ylab="", measure="Cent.dist", main="Slices (punctuated)", diversity=slat_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(med_dist_centroids_Beck_sliPUN, xlab="", ylab="", measure="Cent.dist", main="Slices (punctuated:acctran)", diversity=slat_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(med_dist_centroids_Beck_sliACC, xlab="", ylab="", measure="Cent.dist", main="Slices (punctuated:deltran)", diversity=slat_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(med_dist_centroids_Beck_sliDEL, xlab="", ylab="", measure="Cent.dist", main="Slices (gradual)", diversity=slat_div_nod, y2lab="Number of tips and nodes")
abline(v= 22, col="red")

#Sum of ranges
plot.disparity(dis_tips_slater, xlab="", ylab="Sum of ranges", measure="Sum.range", main="Intervals (tips only)", diversity=dis_tips_slater$rarefaction, y2lab="", cex.xaxis=0.7)
abline(v= 17.5, col="red")
plot.disparity(dis_nodes_slater, xlab="", ylab="", measure="Sum.range", main="Intervals (tips and nodes)", diversity=dis_nodes_slater$rarefaction, y2lab="", cex.xaxis=0.7)
abline(v= 17.5, col="red")
plot.disparity(dis_ran_slater, xlab="", ylab="", measure="Sum.range", main="Slices (punctuated)", diversity=slat_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_acc_slater, xlab="", ylab="", measure="Sum.range", main="Slices (punctuated:acctran)", diversity=slat_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_del_slater, xlab="", ylab="", measure="Sum.range", main="Slices (punctuated:deltran)", diversity=slat_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_pro_slater, xlab="", ylab="", measure="Sum.range", main="Slices (gradual)", diversity=slat_div_nod, y2lab="Number of tips and nodes")
abline(v= 22, col="red")

#Sum of variance
plot.disparity(dis_tips_slater, xlab="", ylab="Sum of variances", measure="Sum.var", main="Intervals (tips only)", diversity=dis_tips_slater$rarefaction, y2lab="", cex.xaxis=0.7)
abline(v= 17.5, col="red")
plot.disparity(dis_nodes_slater, xlab="", ylab="", measure="Sum.var", main="Intervals (tips and nodes)", diversity=dis_nodes_slater$rarefaction, y2lab="", cex.xaxis=0.7)
abline(v= 17.5, col="red")
plot.disparity(dis_ran_slater, xlab="", ylab="", measure="Sum.var", main="Slices (punctuated)", diversity=slat_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_acc_slater, xlab="", ylab="", measure="Sum.var", main="Slices (punctuated:acctran)", diversity=slat_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_del_slater, xlab="", ylab="", measure="Sum.var", main="Slices (punctuated:deltran)", diversity=slat_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_pro_slater, xlab="", ylab="", measure="Sum.var", main="Slices (gradual)", diversity=slat_div_nod, y2lab="Number of tips and nodes")
abline(v= 22, col="red")

#Product of ranges
plot.disparity(dis_tips_slater, xlab="", ylab="Product of ranges", measure="Prod.range", main="Intervals (tips only)", diversity=dis_tips_slater$rarefaction, y2lab="", cex.xaxis=0.7)
abline(v= 17.5, col="red")
plot.disparity(dis_nodes_slater, xlab="", ylab="", measure="Prod.range", main="Intervals (tips and nodes)", diversity=dis_nodes_slater$rarefaction, y2lab="", cex.xaxis=0.7)
abline(v= 17.5, col="red")
plot.disparity(dis_ran_slater, xlab="", ylab="", measure="Prod.range", main="Slices (punctuated)", diversity=slat_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_acc_slater, xlab="", ylab="", measure="Prod.range", main="Slices (punctuated:acctran)", diversity=slat_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_del_slater, xlab="", ylab="", measure="Prod.range", main="Slices (punctuated:deltran)", diversity=slat_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_pro_slater, xlab="", ylab="", measure="Prod.range", main="Slices (gradual)", diversity=slat_div_nod, y2lab="Number of tips and nodes")
abline(v= 22, col="red")

#Product of variance
plot.disparity(dis_tips_slater, xlab="", ylab="Product of variances", measure="Prod.var", main="Intervals (tips only)", diversity=dis_tips_slater$rarefaction, y2lab="", cex.xaxis=0.7)
abline(v= 17.5, col="red")
plot.disparity(dis_nodes_slater, xlab="", ylab="", measure="Prod.var", main="Intervals (tips and nodes)", diversity=dis_nodes_slater$rarefaction, y2lab="", cex.xaxis=0.7)
abline(v= 17.5, col="red")
plot.disparity(dis_ran_slater, xlab="Time (Ma)", ylab="", measure="Prod.var", main="Slices (punctuated)", diversity=slat_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_acc_slater, xlab="Time (Ma)", ylab="", measure="Prod.var", main="Slices (punctuated:acctran)", diversity=slat_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_del_slater, xlab="Time (Ma)", ylab="", measure="Prod.var", main="Slices (punctuated:deltran)", diversity=slat_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_pro_slater, xlab="Time (Ma)", ylab="", measure="Prod.var", main="Slices (gradual)", diversity=slat_div_nod, y2lab="Number of tips and nodes")
abline(v= 22, col="red")

par(op)