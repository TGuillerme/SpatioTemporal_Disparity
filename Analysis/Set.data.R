#'Funtion' for setting the data for the pco analysis.

set.data<-function(data, with.anc.matrix=FALSE)
{
    #SANITIZING
    #data
    check.class(data, "character", ' must be \"Slater\" or \"Beck\".')
    check.length(data, 1, ' must be \"Slater\" or \"Beck\".')
    if(data != "Slater") {
        if(data != "Beck") {
            stop('Data must be \"Slater\" or \"Beck\".')
        }
    }

    #with.anc.matrix
    check.class(with.anc.matrix, "logical", ' must be logical.')

    #Reading the data
    if (data == 'Slater') {
        #Data input (Slater)
        #Slater.table<-read.table("https://raw.githubusercontent.com/TGuillerme/SpatioTemporal_Disparity/master/Data/2013-Slater-MEE-morpho.table", header=F, sep=" ", row.names=1) 
        Slater.table<-read.table("../Data/2013-Slater-MEE-morpho.table", header=F, sep=" ", row.names=1) 
        #Slater.tree<-read.nexus('https://raw.githubusercontent.com/TGuillerme/SpatioTemporal_Disparity/master/Data/2013-Slater-MEE-TEM.tre')
        Slater.tree<-read.nexus('../Data/2013-Slater-MEE-TEM.tre')

        #Remove species with only missing data before hand
        Slater.table<-Slater.table[-c(as.vector(which(apply(as.matrix(Slater.table), 1, function(x) levels(as.factor(x))) == "?"))),]
        Slater.table<-Slater.table[-c(grep("Aegialodon", row.names(Slater.table)), grep("Murtoilestes", row.names(Slater.table))),] #Aegialodon and Murtoilestes are bugged

        #Cleaning the tree and the table
        tree<-clean.tree(Slater.tree, Slater.table)
        table<<-clean.table(Slater.table, Slater.tree)
        cat("\nCreated the morphological characters table as:\ntable")

        #Making the tree binary
        tree<<-bin.tree(tree)
        #adding node names
        tree$node.label<<-paste("n",seq(1:Nnode(tree)), sep="")
        cat("\nCreated the tree as:\ntree")

        #test
        expect_equal(Ntip(tree), nrow(table))
        suppressWarnings({eucl.table<-dist(table, method = "euclidean")})
        pco<-pcoa(eucl.table)
        expect_is(pco, "pcoa")

        #Characters list
        Mandible<<-seq(1:39)
        Premolars<<-seq(40:56)
        Molar_morpho<<-seq(57:126)
        Molar_wear<<-seq(127:141)
        Dental_other<<-seq(142:173)
        Vertebrae<<-seq(174:185)
        Shoulder<<-seq(186:208)
        Forelimb<<-seq(209:224)
        Pelvic<<-seq(225:237)
        Hindlimb<<-seq(238:287)
        Postcranial_other<<-seq(288:292)
        Basicranium<<-seq(293:365)
        Middle_ear<<-seq(366:387)
        Cranial_other<<-seq(388:436)
        Cranial_vault<<-seq(437:443)
        Soft_tissue<<-seq(444:445)
        character.details<<-c("Mandible","Premolars","Molar_morpho","Molar_wear","Dental_other","Vertebrae","Shoulder","Forelimb","Pelvic","Hindlimb","Postcranial_other","Basicranium","Middle_ear","Cranial_other","Cranial_vault","Soft_tissue")
        
        #Groups
        Dental<<-c(Premolars, Molar_morpho, Molar_wear, Dental_other)
        Cranial<<-c(Mandible, Basicranium, Middle_ear, Cranial_other, Cranial_vault)
        PostCranial<<-c(Vertebrae, Shoulder, Forelimb, Pelvic, Hindlimb, Postcranial_other)
        cat("\nCreated the different morphological characters categories as:\nDental\nCranial\nPostCranial\nCharacters details list is available in the object:\ncharacter.details")

        #taxonomy.list
        taxonomy.list<<-list("Australosphenids"=160, "Marsupialiomorphs"=122, "Placentaliomorphs"=106, "Stem_mammaliforms"=c(91, 101), "Stem_theriforms"=c(102,105)) #Slater
        cat("\nCreated the taxonomical list as:\ntaxonomy.list")

        #Age slices
        slices<<-c(0,25,50,75,100,125,150,175,200)
        cat("\nCreated the age slices list as:\nslices")

        if(with.anc.matrix==TRUE) {
            #load("https://raw.githubusercontent.com/TGuillerme/SpatioTemporal_Disparity/master/Data/slater.mat.Rda")
            load("../Data/slater.mat.Rda")
            anc.matrix.save<<-anc.matrix.save.slater
            cat("\nCreated the ancestral matrix in \'ML\' as:\nanc.matrix.save")
        }
    }

    if (data == 'Beck') {
        #Data input (Beck)
        #Beck.table<-read.table("https://raw.githubusercontent.com/TGuillerme/SpatioTemporal_Disparity/master/Data/2014-Beck-ProcB-morpho.table", header=F, sep=" ", row.names=1) 
        Beck.table<-read.table("../Data/2014-Beck-ProcB-morpho.table", header=F, sep=" ", row.names=1) 
        #Beck.tree<-read.nexus('https://raw.githubusercontent.com/TGuillerme/SpatioTemporal_Disparity/master/Data/2014-Beck-ProcB-TEM.tre')
        Beck.tree<-read.nexus('../Data/2014-Beck-ProcB-TEM.tre')

        #Remove species with only missing data before hand
        if (any(apply(as.matrix(Beck.table), 1, function(x) levels(as.factor((x)))) == "?")) {
            Beck.table<-Beck.table[-c(as.vector(which(apply(as.matrix(Beck.table), 1, function(x) levels(as.factor(x))) == "?"))),]
        }

        #Cleaning the tree and the table
        tree<-clean.tree(Beck.tree, Beck.table)
        table<<-clean.table(Beck.table, Beck.tree)
        cat("\nCreated the morphological characters table as:\ntable")
        #Making the tree binary
        tree<<-bin.tree(tree)
        #adding node names
        tree$node.label<<-paste("n",seq(1:Nnode(tree)), sep="")
        cat("\nCreated the tree as:\ntree")   

        #test
        expect_equal(Ntip(tree), nrow(table))
        suppressWarnings({eucl.table<-dist(table, method = "euclidean")})
        pco<-pcoa(eucl.table)
        expect_is(pco, "pcoa")

        #Characters list
        Dentition_general<<-c(seq(1:5), 420)
        Incisors<<-c(seq(6:22), 417)
        Canine<<-seq(23:28)
        Premolars<<-seq(29:60)
        Molars<<-c(seq(61:127), 421)
        Mandible<<-seq(128:158)
        Rostrum<<-seq(159:182)
        Palate<<-seq(183:195)
        Zygoma<<-c(seq(196:202), 418)
        Orbit<<-seq(203:225)
        Braincase<<-seq(226:233)
        Mesocranium<<-seq(234:249)
        Basicranium<<-c(seq(250:332), 419)
        Occiput<<-seq(333:338)
        Vertebrae<<-seq(339:351)
        Forelimb<<-seq(352:369)
        Hindlimb<<-seq(370:416)
        character.details<<-c("Dentition_general","Incisors","Canine","Premolars","Molars","Mandible","Rostrum","Palate","Zygoma","Orbit","Braincase","Mesocranium","Basicranium","Occiput","Vertebrae","Forelimb","Hindlimb")

        #Groups
        Dental<<-c(Dentition_general, Incisors, Canine, Premolars, Molars)
        Cranial<<-c(Mandible, Rostrum, Palate, Zygoma, Braincase, Orbit, Mesocranium, Basicranium, Occiput)
        PostCranial<<-c(Vertebrae, Forelimb, Hindlimb)
        cat("\nCreated the different morphological characters categories as:\nDental\nCranial\nPostCranial\nCharacters details list is available in the object:\ncharacter.details")

        #taxonomy.list
        taxonomy.list<<-list("Placental"=153, "Stem_placental"=c(103,153)) #Beck
        cat("\nCreated the taxonomical list as:\ntaxonomy.list")

        #Age slices
        slices<<-c(0, 40, 50, 60, 70, 80, 90, 100, 110) #Beck
        cat("\nCreated the age slices list as:\nslices")

        if(with.anc.matrix==TRUE) {
            #load("https://raw.githubusercontent.com/TGuillerme/SpatioTemporal_Disparity/master/Data/beck.mat.Rda")
            load("../Data/beck.mat.Rda")
            anc.matrix.save<<-anc.matrix.save.beck
            cat("\nCreated the ancestral matrix in \'ML\' as:\nanc.matrix.save")
        }
    }
}
