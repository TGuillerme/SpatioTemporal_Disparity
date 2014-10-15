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
source("set.data.R")

#Choose the data set to load:
##"Beck" for the Beck & Lee 2014 Proc.B containing fossil and living placental mammals
##"Slater" for the Slater 2013 MEE containing fossil and living mammaliforms
#Use the with.anc.matrix option to skip the ancestral states estimation part (long)
set.data("Beck", with.anc.matrix=TRUE)

#Ancestral state matrix
##Skiped for now
#anc.matrix.save<-anc.state(tree, table, model='ML', verbose=TRUE)

#Recalculating the matrix with a 0.95 probability lower limit
##Using the load anc.matrix.save data
anc.matrix<-anc.unc(anc.matrix.save, 0.95)

#Submatrix
submatrix<-anc.matrix
##Using the loaded character list
submatrix$state<-submatrix$state[, Dental]

#Calculating the PCO/MDS with a scaled euclidean distance matrix and removing the NAs
pco<-pco.std(submatrix, distance="euclidean", scale=TRUE, center=FALSE, na.rm=TRUE, correction="none")

#Visualising the axis variance load
plot.std(pco, legend=TRUE)

#Creating the pco.scores object (containing the axis and the taxonomy)
pco.scores<-as.pco.scores(tree, pco, n.axis=2, taxonomy.list)

#Full pco plot
plot.std(pco.scores, legend=TRUE, main="Full character-space")

#Creating the slice list
std.slice_acc<-std.slice(tree, pco.scores, slices, method="ACCTRAN")
std.slice_del<-std.slice(tree, pco.scores, slices, method="DELTRAN")

#Plot the two series of slices
plot.std(std.slice_acc, legend=TRUE, pars=c(3,3), pos.leg=c(-5,6))
plot.std(std.slice_del, legend=TRUE, pars=c(3,3), pos.leg=c(-5,6))









#Check the for different parameters spaces
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
