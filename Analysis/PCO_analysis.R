#Header
if(grep("TGuillerme", getwd())) {
    setwd('~/PhD/Projects/SpatioTemporal_Disparity/Analysis')
} else {
    warning("You might have to change the directory!")
}
source("functions.R")

#data='Slater'
data='Beck'

if (data == 'Slater') {
    #Data input (Slater)
    Slater.table<-read.table("../Data/2013-Slater-MEE-morpho.table", header=F, sep=" ", row.names=1) 
    Slater.tree<-read.nexus('../Data/2013-Slater-MEE-TEM.tre')

    #Remove species with only missing data before hand
    if (any(apply(as.matrix(Beck.table), 1, function(x) levels(as.factor((x)))) == "?")) {
        Slater.table<-Slater.table[-c(as.vector(which(apply(as.matrix(Slater.table), 1, function(x) levels(as.factor(x))) == "?"))),]
    }
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

    #Characters list
    Mandible<-seq(1:39)
    Premolars<-seq(40:56)
    Molar_morpho<-seq(57:126)
    Molar_wear<-seq(127:141)
    Dental_other<-seq(142:173)
    Vertebrae<-seq(174:185)
    Shoulder<-seq(186:208)
    Forelimb<-seq(209:224)
    Pelvic<-seq(225:237)
    Hindlimb<-seq(238:287)
    Postcranial_other<-seq(288:292)
    Basicranium<-seq(293:365)
    Middle_ear<-seq(366:387)
    Cranial_other<-seq(388:436)
    Cranial_vault<-seq(437:443)
    Soft_tissue<-seq(444:445)
    #Groups
    Dental<-c(Premolars, Molar_morpho, Molar_wear, Dental_other)
    Cranial<-c(Mandible, Basicranium, Middle_ear, Cranial_other, Cranial_vault)
    PostCranial<-c(Vertebrae, Shoulder, Forelimb, Pelvic, Hindlimb, Postcranial_other)
}

if (data == 'Beck') {
    #Data input (Beck)
    Beck.table<-read.table("../Data/2014-Beck-ProcB-morpho.table", header=F, sep=" ", row.names=1) 
    Beck.tree<-read.nexus('../Data/2014-Beck-ProcB-TEM.tre')

    #Remove species with only missing data before hand
    if (any(apply(as.matrix(Beck.table), 1, function(x) levels(as.factor((x)))) == "?")) {
        Beck.table<-Beck.table[-c(as.vector(which(apply(as.matrix(Beck.table), 1, function(x) levels(as.factor(x))) == "?"))),]
    }

    #Cleaning the tree and the table
    tree<-clean.tree(Beck.tree, Beck.table)
    table<-clean.table(Beck.table, Beck.tree)
    #Making the tree binary
    tree<-bin.tree(tree)
    #adding node names
    tree$node.label<-paste("n",seq(1:Nnode(tree)), sep="")   

    #test
    expect_equal(Ntip(tree), nrow(table))
    suppressWarnings({eucl.table<-dist(table, method = "euclidean")})
    pco<-pcoa(eucl.table)
    expect_is(pco, "pcoa")

    #Characters list
    Dentition_general<-c(seq(1:5), 420)
    Incisors<-c(seq(6:22), 417)
    Canine<-seq(23:28)
    Premolars<-seq(29:60)
    Molars<-c(seq(61:127), 421)
    Mandible<-seq(128:158)
    Rostrum<-seq(159:182)
    Palate<-seq(183:195)
    Zygoma<-c(seq(196:202), 418)
    Orbit<-seq(203:225)
    Braincase<-seq(226:233)
    Mesocranium<-seq(234:249)
    Basicranium<-c(seq(250:332), 419)
    Occiput<-seq(333:338)
    Vertebrae<-seq(339:351)
    Forelimb<-seq(352:369)
    Hindlimb<-seq(370:416)
    Dental<-c(Dentition_general, Incisors, Canine, Premolars, Molars)
    Cranial<-c(Mandible, Rostrum, Palate, Zygoma, Braincase, Orbit, Mesocranium, Basicranium, Occiput)
    PostCranial<-c(Vertebrae, Forelimb, Hindlimb)
}

#Ancestral state matrix
anc.matrix.save<-anc.state(tree, table, model='ML', verbose=TRUE)

#Recalculating the matrix with a 0.95 probability lower limit
anc.matrix<-anc.unc(anc.matrix.save, 0.95)

#Submatrix
submatrix<-anc.matrix
submatrix$state<-submatrix$state

#Calculating the PCO/MDS with a scaled euclidean distance matrix and removing the NAs
pco<-pco.std(submatrix, distance="euclidean", scale=TRUE, center=FALSE, na.rm=TRUE, correction="none")

#Visualising the axis variance load
plot.std(pco, legend=TRUE)

taxonomy.list<-list("Placental"=153, "Stem_placental"=c(103,153)) #Beck
#taxonomy.list<-list("Australosphenids"=160, "Marsupialiomorphs"=122, "Placentaliomorphs"=106, "Stem_mammaliforms"=c(91, 101), "Stem_theriforms"=c(102,105)) #Slater

#Creating the pco.scores object (containing the axis and the taxonomy)
pco.scores<-as.pco.scores(tree, pco, n.axis=2, taxonomy.list)

#Full pco plot
plot.std(pco.scores, legend=TRUE)









#lim
xlim=c(min(pco.scores[,1]), max(pco.scores[,1]))
ylim=c(min(pco.scores[,2]), max(pco.scores[,2]))

#Time slices pco plots Beck
pco.slice(tree, pco.scores, c(0, 40, 50, 60, 70, 80, 90, 100, 110), 'ACCTRAN', tax.col="taxonomy", legend=FALSE, pars=c(3,3), xlim=xlim, ylim=ylim)
dev.new()
pco.slice(tree, pco.scores, c(0, 40, 50, 60, 70, 80, 90, 100, 110), 'DELTRAN', tax.col="taxonomy", legend=FALSE, pars=c(3,3), xlim=xlim, ylim=ylim)

#Time slices pco plots Slate
pco.slice(tree, pco.scores, c(0,25,50,75,100,125,150,175,200), 'ACCTRAN', tax.col="taxonomy", legend=FALSE, pars=c(3,3), xlim=xlim, ylim=ylim)
pco.slice(tree, pco.scores, c(0,25,50,75,100,125,150,175,200), 'DELTRAN', tax.col="taxonomy", legend=FALSE, pars=c(3,3), xlim=xlim, ylim=ylim)
