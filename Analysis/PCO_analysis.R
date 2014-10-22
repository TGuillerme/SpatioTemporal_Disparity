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

#Choose the data set to load:
##"Beck" for the Beck & Lee 2014 Proc.B containing fossil and living placental mammals
##"Slater" for the Slater 2013 MEE containing fossil and living mammaliforms
#Use the with.anc.matrix option to skip the ancestral states estimation part (long)
#source("Slater.data.R")
#source("Beck.data.R")


#Ancestral state matrix
##Skiped for now
#anc.matrix.save<-anc.state(tree, table, method='ML', verbose=TRUE)

#Recalculating the matrix with a 0.95 probability lower limit
##Using the load anc.matrix.save data
anc.matrix<-anc.unc(anc.matrix.save, 0.95)

#Submatrix
submatrix<-anc.matrix
##Using the loaded character list
submatrix$state<-submatrix$state[,]

#Calculating the PCO/MDS with a scaled euclidean distance matrix and removing the NAs
pco<-pco.std(submatrix, distance="euclidean", scale=TRUE, center=FALSE, na.rm=TRUE, correction="none")

#Visualising the axis variance load
plot.std(pco, legend=TRUE)

#Creating the pco.scores object (containing the axis and the taxonomy)
#Change the taxonomy list into group list containing not only taxonomic data (BM, etc...).
pco.scores<-as.pco.scores(tree, pco, n.axis=2, taxonomy.list)

#Full pco plot
plot.std(pco.scores, legend=TRUE, main="Full character-space")

#Creating the slice list
std.slice_acc<-std.slice(tree, pco.scores, slices, method="ACCTRAN")
std.slice_del<-std.slice(tree, pco.scores, slices, method="DELTRAN")

#Plot the two series of slices
plot.std(std.slice_acc, legend=TRUE, pars=c(3,3))
plot.std(std.slice_del, legend=TRUE, pars=c(3,3))
