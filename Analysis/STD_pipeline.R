#Graeme's tips
#Use MOD distance
#Use ACE only on nodes, not on tips
#Time series on the distance matrix, not on the ordination



##################################################
#------------------------------------------------
#STD pipelined analysis based on Claddis package
#------------------------------------------------
##################################################


######################
#Data input
######################

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

#Load the data
source("Test.data.R")


######################
#Creating the matrix with ACE
######################


#Renaming the matrix to match with Graeme's workshop
#Choose one of the following matrices:
#nexus.data<-euarch.nex #Tips and nodes
#nexus.data<-STD.nex #Tips and nodes STD method
nexus.data<-CLADDIS.nex #Tips and nodes CLADDIS method
#Problem with CLADDIS method: no account for uncertainty


######################
#Running the PCO
######################


#Safe taxonomic reduction
#safe.data <- SafeTaxonomicReduction(nexus.data) #Removes nodes

#Distance matrix
dist.data <- MorphDistMatrix(nexus.data)

#Trim data
trimmed.max.data <-TrimMorphDistMatrix(dist.data$max.dist.matrix)

#Run the PCO
pco.data <- cmdscale(trimmed.max.data$dist.matrix, k=nrow(trimmed.max.data$dist.matrix) - 1, add=T)$points


######################
#Selecting the axis
######################


#Axis loading
axis.load <- apply(pco.data, 2, var) / sum(apply(pco.data, 2, var)) * 100

#Selecting the axis with a 95% threshold
threshold=95
selected.axis<-which(cumsum(axis.load) < threshold)

#Plot the cumsum
op<-par(mfrow=c(1,2))
plot(axis.load, type="l", xlab="Ordination axis", ylab="Percentage variance")
barplot(cumsum(axis.load),xlab="Ordination axis", ylab="Proportional cumulative variance")
abline(h=threshold)
text(1, 97, paste(threshold, "% of cumulative variance", sep=""), pos=4, cex=0.5)
text(1, 90, paste(length(which(cumsum(axis.load) < threshold)), "selected axis"), pos=4, cex=0.5)
par(op)


######################
#Plotting the tree
######################


#Renaming the tree
tree.data<-euarch.tree

#Tree ages (useless?)
ages.data<-tree.age(tree.data)
tree.data$root.time<-max(ages.data[,1])
#FAD/LAD
ages.data<-data.frame("FAD"=tree.age(tree)[1:Ntip(tree),1], "LAD"=tree.age(tree)[1:Ntip(tree),1], row.names=tree.age(tree)[1:Ntip(tree),2])

#Plot the tree
geoscalePhylo(ladderize(tree.data), cex.age=0.6, cex.ts=0.8, cex.tip=1)

######################
#Rates analysis
######################

do.rate.analysis=FALSE
if(do.rate.analysis==TRUE) {
    tree.tmp<-tree.data
    tree.tmp$node.label<-NULL
    rate.data <- DiscreteCharacterRate(tree.tmp, nexus.data, seq(tree.tmp$root.time, tree.tmp$root.time - max(diag(vcv(tree.tmp))), length.out=6), alpha=0.01) #remove tip labels
    edge.color <- rep("black", nrow(tree.data$edge))
    edge.color[which(rate.data$branch.results[, "ml.signif.hi"] == 1)] <- "red"
    edge.color[which(rate.data$branch.results[, "ml.signif.lo"] == 1)] <- "blue"
    geoscalePhylo(ladderize(tree.data), cex.age=0.6, cex.ts=0.8, cex.tip=1, edge.color=edge.color[match(ladderize(tree.data)$edge[, 2], tree.data$edge[,2])])

    node.color <- rep("white", nrow(rate.data$node.results))
    node.color[which(rate.data$node.results[, "ml.signif.hi.ib"] == 1)] <- "red"
    node.color[which(rate.data$node.results[, "ml.signif.lo.ib"] == 1)] <- "blue"
    node.color[which(is.na(rate.data$node.results[, "ml.signif.lo.ib"]))] <- NA
    geoscalePhylo(tree.data, cex.age=0.6, cex.ts=0.8, cex.tip=1)
    nodelabels(node=rate.data$node.results[, "node"][!is.na(node.color)], pch=21, col="black", bg=node.color[!is.na(node.color)])
}


######################
#Plotting the results
######################


#tree
plot.tree<-tree.data

#tips + nodes PCO
PCOx<-1
PCOy<-2
all.PCOx<-pco.data[,1]
all.PCOy<-pco.data[,2]

#These can be modified further into two matrices that give the x and y coordinates needed to plot each branch:
branch.xs <- cbind(all.PCOx[plot.tree$edge[, 1]], all.PCOx[plot.tree$edge[, 2]])
branch.ys <- cbind(all.PCOy[plot.tree$edge[, 1]], all.PCOy[plot.tree$edge[, 2]])

#Plot
main="Euarchontoglires through time - TEST"
plot(pco.data[, PCOx], pco.data[, PCOy], xlab=paste("PCO ", PCOx, " (", round(scree.data[PCOx], 2), "% variance)", sep=""), ylab=paste("PCO ", PCOy, " (", round(scree.data[PCOy], 2), "% variance)", sep=""), type="n", main=main)
#Setting phylogenetic groups
Glires<-1:9
Euarchonta<-10:20
#Plotting the Glires tree
for(i in 1:nrow(branch.xs[Glires,])){
    lines(x=branch.xs[i,], y=branch.ys[i,], col="lightgreen", lwd=2)
}
#Plotting the Euarchonta tree
for(i in 10:nrow(branch.xs)){
    lines(x=branch.xs[i,], y=branch.ys[i,], col="lightblue", lwd=2)
}

#Silly ordering for plotting the gradient color
order.table<-table
order.table$extant<-table$log_bm
order.table<-order.table[ do.call(order, order.table), ]
BM_gradient<-colorRampPalette(c("yellow", "red"))
order.table$colors<-BM_gradient(length(order.table$extant))[as.numeric(cut(order.table$extant, breaks = nrow(order.table)))]
#Adding the tips and nodes with their BM values
points(pco.data[, PCOx], pco.data[, PCOy], pch=19, col=order.table$colors)
points(pco.data[, PCOx], pco.data[, PCOy], pch=1, col="black")
#Adding tip/node names
text(pco.data[, PCOx], pco.data[, PCOy], rownames(pco.data), cex=0.6, pos=3)

#Plotting the Convex Hull per groups
groups=length(levels(as.factor(table[,5])))

#Plotting the polygons
for(group in 1:groups) {
    n<-which(table[,5] == levels(as.factor(table[,5]))[group])
    chull.group<-pco.data[n,c(1:2)]
    polygon(chull.group[chull(chull.group),], border="gray", lty=(group+1))
}

#Global legend
legend(-1, -0.6, c("Glires", "Euarchontes", "Generalist", "Specialist"), col=c("lightgreen","lightblue","grey", "grey"), lty=c(1,1,2,3), bty="n", cex=0.8)
legend(-0.5, -0.6, c("Low body mass", "High body mass"), col=c("yellow","red"), pch=19, bty="n", cex=0.8)


######################
#Disparity analysis
######################

#NPMANOVA on diet (check if categories occupy a different position in morphospace (overlap or not)) (e.g. Stayton 2005 and Ruta 2013)
cat.col=5
PC.man<-adonis(pco.data[,selected.axis]~table[,cat.col], data=table, permutation=1000, method="euclidean") #vegan - think permutation is doing jackknife

#Distance from centroid
#Function to calculate mean euclidean distances from the centroid - Sive Finlay 06/05/2014
euc.dist.cent<-function(PCdata, method="euclidean"){
    #Centroid (mean score of each PC axis)
    centroid<-apply(PCdata, 2, mean)

    #Euclidean distances to the centroid
    cent.dist<-NULL
    for (j in 1:nrow(PCdata)){
        cent.dist[j] <- dist(rbind(PCdata[j,], centroid), method=method)
    }
    return(cent.dist)
}
  
gene.dist<-euc.dist.cent(pco.data[which(table[,cat.col]=='generalist'), selected.axis])
spec.dist<-euc.dist.cent(pco.data[which(table[,cat.col]=='specialist'), selected.axis])


#Compare the distance from centroid between the categories
cent.dist<-matrix(nrow=nrow(pco.data[,selected.axis]), ncol=length(levels(as.factor(table[,cat.col]))))
colnames(cent.dist)<-c("dist", "category")
cent.dist[,1]<-c(gene.dist, spec.dist)
cent.dist[,2]<-c(rep("generalist", length(gene.dist)), rep("specialist", length(spec.dist)))

#t-test
comp.cent2 <- t.test(gene.dist, spec.dist)

#Mean and standard error from distance from centroid
mean.se<-matrix(nrow=length(levels(as.factor(table[,cat.col]))), ncol=2)
colnames(mean.se) <- c("mean", "se")
rownames(mean.se) <- c("generalist", "specialist")
mean.se[,1] <- c(mean(gene.dist), mean(spec.dist))
mean.se[,2] <- c(std.error(gmole.cent.dist), std.error(tenrec.cent.dist))

sive=FALSE
if(sive==TRUE){
    ##########################
    #Compare positions in morphospace
    ########################

    #NPMANOVA of the PC axes  (e.g. Stayton 2005 and Ruta 2013)
      PC.man <- adonis(PC95axes~sp.fam$Family, data=sp.fam, permutations=999, method="euclidean")

    #################
    #Diversity of families based on distances to centroid
    #################

    #Distance from each species to that family's centroid
      gmole.cent.dist <- euc.dist.cent (gmolePC)
      tenrec.cent.dist <- euc.dist.cent (tenrecPC)
      
    #Compare the distances to centroids in tenrecs and golden moles
      cent.dist <- matrix(nrow=nrow(PC95axes), ncol=2)
        colnames(cent.dist) <- c("dist", "group")
        cent.dist[,1] <- c(gmole.cent.dist, tenrec.cent.dist)
        cent.dist[,2] <- c(rep("gmole", length(gmole.cent.dist)), rep("tenrec", length(tenrec.cent.dist)))


    #Compare the two groups with a t test
      comp.cent <- t.test(as.numeric(cent.dist[,1]) ~ cent.dist[,2])

      
    #Mean and standard error of those distances from the centroid
      mean.se <- matrix(nrow=2, ncol=2)
        colnames(mean.se) <- c("mean", "se")
        rownames(mean.se) <- c("gmole", "tenrec")
      mean.se[,1] <- c(mean(gmole.cent.dist), mean(tenrec.cent.dist))
      mean.se[,2] <- c(std.error(gmole.cent.dist), std.error(tenrec.cent.dist))
}


#Jackknife permutation test (sample size effect)

sive=FALSE
if(sive==TRUE){
    ################################
    #Pairwise permutation tests
    ##############################
    #Check that the significant differences are not just an artefact of differences in sample size

    #Observed differences in the mean
      obs.mean.diff <- mean.se[2,1] - mean.se[1,1]  

    #Permutation test for significant difference in the mean
      perm.mean <- group.diff(1000, sp.fam$Family, PC95axes, mean)
    #Test for significant difference
      perm.mean.pvalue <- pvalue.dist(perm.mean, obs.mean.diff)


    #Summary table of the results
          perm.res.summary <- matrix(NA, nrow=1, ncol=6)
            colnames(perm.res.summary) <- c("obs.tenrec", "obs.gmole", "obs.diff", "perm.min", "perm.max", "pvalue") 
            perm.res.summary[1,] <- c(mean.se[2,1], mean.se[1,1], obs.mean.diff, min(perm.mean), max(perm.mean), perm.mean.pvalue)
}