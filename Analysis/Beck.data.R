#Data input (Beck)

#Data cleaning repeatability
if(length(grep("Windows", Sys.info()["sysname"]) == 1) {
    stop("UNIX language based machines only.")
} else {
    #Isolate the morphological data from the raw Beck&Lee matrix (UNIX!)
    system("
        #Changing the matrix dimension
        sed 's/Dimensions ntax=106 nchar=8959;/Dimensions ntax=106 nchar=421;/' ../Data/2014-Beck-ProcB-matrix-raw.nex |
        #Changing the matrix type
        sed 's/datatype=mixed(standard:1-421,DNA:422-8959)/datatype=standard/' |
        #removing the phylogenetic analysis
        sed '237,1408d' |
        #removing the molecular data
        sed '126,231d' |
        #removing the comments in the header
        sed '3,13d' > ../Data/2014-Beck-ProcB-matrix-morpho.nex
    ")
}

#Read nexus table
Beck.nex<-ReadMorphNexus("../Data/2014-Beck-ProcB-matrix-morpho.nex")
Beck.table<-Beck.nex$matrix

#Read tree
Beck.tree<-read.nexus('../Data/2014-Beck-ProcB-TEM.tre')

#Remove species with only missing data before hand
if (any(apply(as.matrix(Beck.table), 1, function(x) levels(as.factor((x)))) == "?")) {
    Beck.table<-Beck.table[-c(as.vector(which(apply(as.matrix(Beck.table), 1, function(x) levels(as.factor(x))) == "?"))),]
}

#Cleaning the tree and the table
tree<-clean.tree(Beck.tree, Beck.table)
table<-clean.table(Beck.table, Beck.tree)
cat("Created the morphological characters table as:\ntable")
#Making the tree binary
tree<-bin.tree(tree)
#adding node names
tree$node.label<-paste("n",seq(1:Nnode(tree)), sep="")
cat("Created the tree as:\ntree")   

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
character.details<-c("Dentition_general","Incisors","Canine","Premolars","Molars","Mandible","Rostrum","Palate","Zygoma","Orbit","Braincase","Mesocranium","Basicranium","Occiput","Vertebrae","Forelimb","Hindlimb")

#Groups
Dental<-c(Dentition_general, Incisors, Canine, Premolars, Molars)
Cranial<-c(Mandible, Rostrum, Palate, Zygoma, Braincase, Orbit, Mesocranium, Basicranium, Occiput)
PostCranial<-c(Vertebrae, Forelimb, Hindlimb)
cat("Created the different morphological characters categories as:\nDental\nCranial\nPostCranial\nCharacters details list is available in the object:\ncharacter.details\n")

#taxonomy.list
taxonomy.list<-list("Placental"=153, "Stem_placental"=c(103,153)) #Beck
cat("Created the taxonomical list as:\ntaxonomy.list\n")

#Age slices
slices<-c(0, 40, 50, 60, 70, 80, 90, 100, 110)
cat("Created the age slices list as:\nslices\n")

#load("https://raw.githubusercontent.com/TGuillerme/SpatioTemporal_Disparity/master/Data/beck.mat.Rda")
load("../Data/beck.ML.mat.rda")
anc.matrix.save<-anc.matrix.beck.ML
cat("Created the ancestral matrix in \'ML\' as:\nanc.matrix.save\n")