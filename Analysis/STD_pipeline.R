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

#STD pipelined analysis based on Claddis package

#Data input
source("Test.data.R")

#Renaming the matrix to match with Graeme's workshop
#Choose one of the following matrices:
#nexus.data<-euarch.nex #Tips and nodes
#nexus.data<-STD.nex #Tips and nodes STD method
nexus.data<-CLADDIS.nex #Tips and nodes CLADDIS method

#Safe taxonomic reduction
#safe.data <- SafeTaxonomicReduction(nexus.data) #Removes nodes

#Distance matrix
dist.data <- MorphDistMatrix(nexus.data)

#Trim data
trimmed.max.data <-TrimMorphDistMatrix(dist.data$max.dist.matrix)

#Run the PCO
pco.data <- cmdscale(trimmed.max.data$dist.matrix, k=nrow(trimmed.max.data$dist.matrix) - 1, add=T)$points

#Cumsum
scree.data <- apply(pco.data, 2, var) / sum(apply(pco.data, 2, var)) * 100

#Plot the cumsum
plot(scree.data, type="l", xlab="Ordination axis", ylab="Percentage variance")

#Plot the PCO (maybe useless)
PCOx <- 1 ; PCOy <- 2
plot(pco.data[, PCOx], pco.data[, PCOy], xlab=paste("PCO ", PCOx, " (", round(scree.data[PCOx], 2), "% variance)", sep=""), ylab=paste("PCO ", PCOy, " (", round(scree.data[PCOy], 2), "% variance)", sep=""), pch=19)
text(pco.data[, PCOx], pco.data[, PCOy], rownames(pco.data), cex=0.5, pos=3)

#Renaming the tree
tree.data<-euarch.tree

#Tree ages (useless?)
ages.data<-tree.age(tree.data)
tree.data$root.time<-max(ages.data[,1])
#FAD/LAD
ages.data<-data.frame("FAD"=tree.age(tree)[1:Ntip(tree),1], "LAD"=tree.age(tree)[1:Ntip(tree),1], row.names=tree.age(tree)[1:Ntip(tree),2])

#Plot the tree
geoscalePhylo(ladderize(tree.data), cex.age=0.6, cex.ts=0.8, cex.tip=1)

#Rates analysis
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

#Plot PCO with tree
#tree
plot.tree<-tree.data

#ace PCO
#PCOx.anc<-pco.data[-(1:Ntip(tree.data)),1]
#PCOy.anc<-pco.data[-(1:Ntip(tree.data)),2]

#tips + nodes PCO
all.PCOx<-pco.data[,1]
all.PCOy<-pco.data[,2]

#These can be modified further into two matrices that give the x and y coordinates needed to plot each branch:
branch.xs <- cbind(all.PCOx[plot.tree$edge[, 1]], all.PCOx[plot.tree$edge[, 2]])
branch.ys <- cbind(all.PCOy[plot.tree$edge[, 1]], all.PCOy[plot.tree$edge[, 2]])

#Plot
plot(pco.data[, PCOx], pco.data[, PCOy], xlab=paste("PCO ", PCOx, " (", round(scree.data[PCOx], 2), "% variance)", sep=""), ylab=paste("PCO ", PCOy, " (", round(scree.data[PCOy], 2), "% variance)", sep=""), type="n", main="Euarchontoglires through time - TEST")
Glires<-1:9
Euarchonta<-10:20
for(i in 1:nrow(branch.xs[Glires,])){
    lines(x=branch.xs[i,], y=branch.ys[i,], col="lightgreen", lwd=2)
}
for(i in 10:nrow(branch.xs)){
    lines(x=branch.xs[i,], y=branch.ys[i,], col="lightblue", lwd=2)
}

order.table<-table
order.table$extant<-table$log_bm
order.table<-order.table[ do.call(order, order.table), ]
BM_gradient<-colorRampPalette(c("yellow", "red"))
order.table$colors<-BM_gradient(length(order.table$extant))[as.numeric(cut(order.table$extant, breaks = nrow(order.table)))]
points(pco.data[, PCOx], pco.data[, PCOy], pch=19, col=order.table$colors)
points(pco.data[, PCOx], pco.data[, PCOy], pch=1, col="black")
text(pco.data[, PCOx], pco.data[, PCOy], rownames(pco.data), cex=0.6, pos=3)

groups=length(levels(as.factor(table[,5])))

for(group in 1:groups) {
    n<-which(table[,5] == levels(as.factor(table[,5]))[group])
    chull.group<-pco.data[n,c(1:2)]
    polygon(chull.group[chull(chull.group),], border="gray", lty=(group+1))
}
legend(-1, -0.6, c("Glires", "Euarchontes", "Generalist", "Specialist"), col=c("lightgreen","lightblue","grey", "grey"), lty=c(1,1,2,3), bty="n", cex=0.8)
legend(-0.5, -0.6, c("Low body mass", "High body mass"), col=c("yellow","red"), pch=19, bty="n", cex=0.8)

