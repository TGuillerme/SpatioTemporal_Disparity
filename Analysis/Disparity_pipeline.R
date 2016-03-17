#####
# Disparity calculation using the dispRity package
####

######
# Before starting
#####

# Loading the package
library(Claddis)
library(dispRity)

# Loading the data
tree <- read.nexus("../Data/2014-Beck-ProcB-TEM.tre")

# Loading the matrix
matrix <- ReadMorphNexus("../Data/2014-Beck-ProcB-matrix-morpho.nex")

# loading the first/last occurrence data
FADLAD <- read.csv("../Data/Beck2014_FADLAD.csv", row.names=1)


#####
# Calculating the ancestral nodes states
#####

# Adding node labels to the tree
tree <- makeNodeLabel(tree, method = "number", prefix = "n")

# Adding root.time to the tree
tree$root.time <- max(tree.age(tree)[,1])

# Ancestral states reconstruction
#matrix_nodes <- AncStateEstMatrix(matrix, tree, estimate.allchars = TRUE)
#TG: Add the 0.95 lower limit for ancestral reconstruction using the anc.unc function.

# Load the data (shortcut)
load("../Data/Beck2014/Beck2014_ancestral_states.Rda")
source("../Functions/disparity/R/sanitizing.R")
source("../Functions/disparity/R/anc.unc.R")
matrix_nodes <- anc.unc(anc_states, 0.95, NA)$state

#####
# Ordinating the data
#####

# Calculating the distance matrix
#dist_matrix <- MorphDistMatrix(matrix)
#TG: Calculate only the gower distance + use the fast algorithm.

# Load the data (shortcut)
load("../Data/Beck2014/Beck2014_distance-nodes95.Rda")
dist_matrix <- dist_nodes95$gower

# Ordinate the data
ord_matrix <- cmdscale(dist_matrix, k=nrow(dist_matrix) - 2, add=T)$points

#####
# Creating the series
#####

# Custom series (pre/after K-Pg mammals)
# Creating the empty data frame
categories <- as.data.frame(matrix(data = rep(NA, nrow(ord_matrix)), nrow = nrow(ord_matrix), ncol = 2, dimnames = list(rownames(ord_matrix))))
# Adding column names
colnames(categories) <- c("KPg", "")
# Adding the factors for all mammals
categories[,2] <- rep("all", nrow(ord_matrix))
# Adding the pre/post factors for the K-Pg limit
categories[which(tree.age(tree)$age <= 66), 1] <- "pre"
categories[which(tree.age(tree)$age > 66), 1] <- "post"
# Creating the custom series
cust_series_matrix <- cust.series(ord_matrix, categories)

# Time series (every 5 Ma)
# Creating the time slices
time <- rev(seq(from = 0, to = 160, by = 5))
# Creating the time series
time_series_matrix <- time.series(ord_matrix, tree, method = "continuous", time = time, model = "gradual", FADLAD = FADLAD, verbose = TRUE)

#####
# Bootstrapping the data
#####

# Bootstrapping the custom series
cust_series_matrix <- boot.matrix(cust_series_matrix, bootstrap = 100) #TG: set to 1000

# Bootstrapping the time series
time_series_matrix <- boot.matrix(time_series_matrix, bootstrap = 100) #TG: set to 1000

#####
# Calculating disparity
#####

# Calculating disparity in the custom series
disparity_cust <- dispRity(cust_series_matrix, metric = c(median, centroids))
#TG: median centroids = average radius of the morphspace. Maybe test sum of centroids: total dispersion of the morphospace (probably sensitive to morphospace).

# Calculating disparity in the time series
disparity_time <- dispRity(time_series_matrix, metric = c(median, centroids))

#####
# Displaying the results
#####

# Summarizing the disparity
summary(disparity_cust)
summary(disparity_time)

# Plotting the results
plot(disparity_cust, type = "polygon")
plot(disparity_time)

#####
# Calculating the effect of disparity
#####

# Is the disparity before and after K-Pg siginificatively different?
test.dispRity(disparity_cust, test = wilcox.test, correction = "bonferroni")

# How much do the distributions overlap?
test.dispRity(disparity_cust, test = bhatt.coeff)
#TG: disparity after K-Pg is significantly different than before?

# Difference between each series?
test.dispRity(disparity_time, test = t.test, correction = "bonferroni", comparisons = "sequential")
test.dispRity(disparity_time, test = bhatt.coeff, comparisons = "sequential")

# Sequential test
plot(disparity_time)
seq_disparity <- test.dispRity(disparity_time, test = sequential.test, family = gaussian, add = TRUE, lines.args = list(col = "red", lty = 3, lwd = 1.5))

# Difference after K-Pg
test.dispRity(get.dispRity(disparity_time, what = c(19:33)), test = t.test, comparisons = "referential")



#####
# Calculating effect of sample size
#####


summary(disparity_cust)
summary(disparity_time)

# Plotting the results
plot(disparity_cust, type = "polygon")
plot(disparity_time)
