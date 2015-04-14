#Loading the data
if(grep("TGuillerme", getwd())) {
    setwd('~/PhD/Projects/SpatioTemporal_Disparity/Analysis')
} else {
    warning("You might have to change the directory!")
}
if(!grep("SpatioTemporal_Disparity/Analysis", getwd())) {
    stop("Wrong directory!\nThe current directory must be:\nSpatioTemporal_Disparity/Analysis/\nYou can clone the whole repository from:\nhttps://github.com/TGuillerme/SpatioTemporal_Disparity")
}

source("functions.R")
source("set.data.R")
set.data("Beck", with.anc.matrix=TRUE)
anc.matrix<-anc.unc(anc.matrix.save, 0.95)
submatrix<-anc.matrix



#Checking for the differences with/without scaling
pco.scores_sTcT<-as.pco.scores(tree, pco.std(submatrix, distance="euclidean", scale=TRUE, center=TRUE, na.rm=TRUE, correction="none"), n.axis=2, taxonomy.list)
pco.scores_sFcT<-as.pco.scores(tree, pco.std(submatrix, distance="euclidean", scale=FALSE, center=TRUE, na.rm=TRUE, correction="none"), n.axis=2, taxonomy.list)
pco.scores_sTcF<-as.pco.scores(tree, pco.std(submatrix, distance="euclidean", scale=TRUE, center=FALSE, na.rm=TRUE, correction="none"), n.axis=2, taxonomy.list)
pco.scores_sFcF<-as.pco.scores(tree, pco.std(submatrix, distance="euclidean", scale=FALSE, center=FALSE, na.rm=TRUE, correction="none"), n.axis=2, taxonomy.list)
op<-par(mfrow=c(2,2)) 
plot.std(pco.scores_sTcT, main="Euclidean, scale=T, center=T")
plot.std(pco.scores_sFcT, main="Euclidean, scale=F, center=T")
plot.std(pco.scores_sTcF, main="Euclidean, scale=T, center=F")
plot.std(pco.scores_sFcF, main="Euclidean, scale=F, center=F")
par(op)

#Checking for the different distances methods
pco.scores_bray<-as.pco.scores(tree, pco.std(submatrix, distance="bray", scale=FALSE, center=FALSE, na.rm=TRUE, correction="none"), n.axis=2, taxonomy.list) #results may be meaningless because data have negative entries in method “bray”
pco.scores_jacc<-as.pco.scores(tree, pco.std(submatrix, distance="jaccard", scale=FALSE, center=FALSE, na.rm=TRUE, correction="none"), n.axis=2, taxonomy.list) #results may be meaningless because data have negative entries in method “jaccard”
pco.scores_eucl<-as.pco.scores(tree, pco.std(submatrix, distance="euclidean", scale=TRUE, center=FALSE, na.rm=TRUE, correction="none"), n.axis=2, taxonomy.list) #
pco.scores_manh<-as.pco.scores(tree, pco.std(submatrix, distance="manhattan", scale=TRUE, center=FALSE, na.rm=TRUE, correction="none"), n.axis=2, taxonomy.list)
op<-par(mfrow=c(2,2)) 
plot.std(pco.scores_bray, main="Bray, scale=F, center=F")
plot.std(pco.scores_jacc, main="Jaccard, scale=F, center=F")
plot.std(pco.scores_eucl, main="Euclidean, scale=T, center=F")
plot.std(pco.scores_manh, main="Manhattan, scale=T, center=F")
par(op)

#Checking for the different distances methods
pco.scores_eucl_nonCor<-as.pco.scores(tree, pco.std(submatrix, distance="euclidean", scale=TRUE, center=FALSE, na.rm=TRUE, correction="none"), n.axis=2, taxonomy.list) #results may be meaningless because data have negative entries in method “bray”
pco.scores_eucl_corLin<-as.pco.scores(tree, pco.std(submatrix, distance="euclidean", scale=TRUE, center=FALSE, na.rm=TRUE, correction="lingoes"), n.axis=2, taxonomy.list) #results may be meaningless because data have negative entries in method “jaccard”
pco.scores_eucl_corCai<-as.pco.scores(tree, pco.std(submatrix, distance="euclidean", scale=TRUE, center=FALSE, na.rm=TRUE, correction="cailliez"), n.axis=2, taxonomy.list) #
op<-par(mfrow=c(2,2)) 
plot.std(pco.scores_eucl_nonCor, main="No correction, scale=F, center=F")
plot.std(pco.scores_eucl_corLin, main="Lingoes correction, scale=F, center=F")
plot.std(pco.scores_eucl_corCai, main="Cailliez correction, scale=T, center=F")
par(op)
#CORRECTION DOESN'T SEEM TO AFFECT THE PCO


#Checking for different parameters spaces
#Full
submatrix.full<-anc.matrix ; submatrix.full$state<-submatrix.full$state
pco.full<-pco.std(submatrix.full, distance="euclidean", scale=TRUE, center=FALSE, na.rm=TRUE, correction="none")
pco.scores.full<-as.pco.scores(tree, pco.full, n.axis=2, taxonomy.list)
#Dental
submatrix.dental<-anc.matrix ; submatrix.dental$state<-submatrix.dental$state[, Dental]
pco.dental<-pco.std(submatrix.dental, distance="euclidean", scale=TRUE, center=FALSE, na.rm=TRUE, correction="none")
pco.scores.dental<-as.pco.scores(tree, pco.dental, n.axis=2, taxonomy.list)
#cranial
submatrix.cranial<-anc.matrix ; submatrix.cranial$state<-submatrix.cranial$state[, Cranial]
pco.cranial<-pco.std(submatrix.cranial, distance="euclidean", scale=TRUE, center=FALSE, na.rm=TRUE, correction="none")
pco.scores.cranial<-as.pco.scores(tree, pco.cranial, n.axis=2, taxonomy.list)
#postcranial
submatrix.PostCranial<-anc.matrix ; submatrix.PostCranial$state<-submatrix.PostCranial$state[, PostCranial]
pco.PostCranial<-pco.std(submatrix.PostCranial, distance="euclidean", scale=TRUE, center=FALSE, na.rm=TRUE, correction="none")
pco.scores.PostCranial<-as.pco.scores(tree, pco.PostCranial, n.axis=2, taxonomy.list)


op<-par(mfrow=c(2,2)) 
plot.std(pco.scores.full, legend=TRUE, pos.leg=c(-10, 10), xlim=c(-10,10), ylim=c(-10,10), main="Full character-space")
plot.std(pco.scores.dental, xlim=c(-10,10), ylim=c(-10,10), main="Dental")
plot.std(pco.scores.cranial, xlim=c(-10,10), ylim=c(-10,10), main="Cranial")
plot.std(pco.scores.PostCranial, xlim=c(-10,10), ylim=c(-10,10), main="PostCranial")
par(op)


#Checking various axis combination
pco.12<-pco.std(submatrix, distance="euclidean", scale=TRUE, center=FALSE, na.rm=TRUE, correction="none")
pco.12$vectors<-pco.12$vectors[,]
pcoscores.12<-as.pco.scores(tree, pco.12, n.axis=2, taxonomy.list)

pco.34<-pco.std(submatrix, distance="euclidean", scale=TRUE, center=FALSE, na.rm=TRUE, correction="none")
pco.34$vectors<-pco.34$vectors[,-c(1:2)]
pcoscores.34<-as.pco.scores(tree, pco.34, n.axis=2, taxonomy.list)

pco.56<-pco.std(submatrix, distance="euclidean", scale=TRUE, center=FALSE, na.rm=TRUE, correction="none")
pco.56$vectors<-pco.56$vectors[,-c(1:4)]
pcoscores.56<-as.pco.scores(tree, pco.56, n.axis=2, taxonomy.list)

pco.78<-pco.std(submatrix, distance="euclidean", scale=TRUE, center=FALSE, na.rm=TRUE, correction="none")
pco.78$vectors<-pco.78$vectors[,-c(1:6)]
pcoscores.78<-as.pco.scores(tree, pco.78, n.axis=2, taxonomy.list)

pco.910<-pco.std(submatrix, distance="euclidean", scale=TRUE, center=FALSE, na.rm=TRUE, correction="none")
pco.910$vectors<-pco.910$vectors[,-c(1:8)]
pcoscores.910<-as.pco.scores(tree, pco.910, n.axis=2, taxonomy.list)

pco.1112<-pco.std(submatrix, distance="euclidean", scale=TRUE, center=FALSE, na.rm=TRUE, correction="none")
pco.1112$vectors<-pco.1112$vectors[,-c(1:10)]
pcoscores.1112<-as.pco.scores(tree, pco.1112, n.axis=2, taxonomy.list)

pco.1314<-pco.std(submatrix, distance="euclidean", scale=TRUE, center=FALSE, na.rm=TRUE, correction="none")
pco.1314$vectors<-pco.1314$vectors[,-c(1:12)]
pcoscores.1314<-as.pco.scores(tree, pco.1314, n.axis=2, taxonomy.list)

pco.1516<-pco.std(submatrix, distance="euclidean", scale=TRUE, center=FALSE, na.rm=TRUE, correction="none")
pco.1516$vectors<-pco.1516$vectors[,-c(1:14)]
pcoscores.1516<-as.pco.scores(tree, pco.1516, n.axis=2, taxonomy.list)

pco.101102<-pco.std(submatrix, distance="euclidean", scale=TRUE, center=FALSE, na.rm=TRUE, correction="none")
pco.101102$vectors<-pco.101102$vectors[,-c(1:100)]
pcoscores.101102<-as.pco.scores(tree, pco.101102, n.axis=2, taxonomy.list)

op<-par(mfrow=c(3,3)) 
plot.std(pcoscores.12, xlim=c(-8,8), ylim=c(-8,8), main="Axis 1 2")
plot.std(pcoscores.34, xlim=c(-8,8), ylim=c(-8,8), main="Axis 3 4")
plot.std(pcoscores.56, xlim=c(-8,8), ylim=c(-8,8), main="Axis 5 6")
plot.std(pcoscores.78, xlim=c(-8,8), ylim=c(-8,8), main="Axis 7 8")
plot.std(pcoscores.910, xlim=c(-8,8), ylim=c(-8,8), main="Axis 9 10")
plot.std(pcoscores.1112, xlim=c(-8,8), ylim=c(-8,8), main="Axis 11 12")
plot.std(pcoscores.1314, xlim=c(-8,8), ylim=c(-8,8), main="Axis 13 14")
plot.std(pcoscores.1516, xlim=c(-8,8), ylim=c(-8,8), main="Axis 15 16")
plot.std(pcoscores.101102, xlim=c(-8,8), ylim=c(-8,8), main="Axis 101 102")
par(op)
