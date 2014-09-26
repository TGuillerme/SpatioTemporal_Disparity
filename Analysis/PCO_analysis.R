#Header
if(grep("TGuillerme", getwd())) {
    setwd('~/PhD/Projects/SpatioTemporal_Disparity/Analysis')
} else {
    warning("You might have to change the directory!")
}
library(ape)
source("functions.R")

#Data input
Slater.table<-read.table("../Data/2013-Slater-MEE-morpho.table", header=F, sep=" ", row.names=1) 
Slater.tree<-read.nexus('../Data/2013-Slater-MEE-TEM.tre')

#Remove species with only missing data before hand
Slater.table<-Slater.table[-c(as.vector(which(apply(as.matrix(Slater.table), 1, function(x) levels(as.factor(x))) == "?"))),]
Slater.table<-Slater.table[-c(grep("Aegialodon", row.names(Slater.table)), grep("Murtoilestes", row.names(Slater.table))),] #Aegialodon and Murtoilestes are bugged

#Cleaning the tree and the table
tree<-clean.tree(Slater.tree, Slater.table)
table<-clean.table(Slater.table, Slater.tree)
#Making the tree binary
tree<-bin.tree(tree)
#adding node names
tree$node.label<-paste("n",seq(1:Nnode(tree)), sep="")


#test
expect_equal(Ntip(tree), nrow(table))
suppressWarnings({eucl.table<-dist(table, method = "euclidean")})
pco<-pcoa(eucl.table)
expect_is(pco, "pcoa")


#Ancestral state matrix
anc.matrix<-anc.state(tree, table, model='ML', verbose=TRUE)

#PCO
dist.matrix<-dist(anc.matrix$state, method="euclidian")
pco<-pcoa(dist.matrix)
pco.scores<-pco$vectors

#Clades
#Australosphenids (monotremes and relatives) node 160
pco.scores<-set.group(tree, pco.scores, type='clade', node=160, name='Australosphenids')
#Marsupialiomorphs (marsupials and relatives) node 122
pco.scores<-set.group(tree, pco.scores, tax.col="taxonomy", type='clade', node=122, name='Marsupialiomorphs')
#Placentaliomorphs (placentals and relatives) node 106
pco.scores<-set.group(tree, pco.scores, tax.col="taxonomy", type='clade', node=106, name='Placentaliomorphs')
#Stem mammaliforms node 91 to 101
pco.scores<-set.group(tree, pco.scores, tax.col="taxonomy", type='grade', node=c(91, 101), name='Stem_mammaliforms')
#Stem theriforms node 102 105
pco.scores<-set.group(tree, pco.scores, tax.col="taxonomy", type='grade', node=c(102, 105), name='Stem_theriforms')

#Full pco plot
plot.pco(pco.scores, "taxonomy", main="Entire \"morphospace\"", legend=TRUE)

#Time slices pco plots
pco.slice(tree, pco.scores, 9, 'ACCTRAN', tax.col="taxonomy", legend=FALSE)