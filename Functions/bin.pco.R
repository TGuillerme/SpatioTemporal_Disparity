##########################
#bin.pco
##########################
#Select the number of taxa per bin.
#v0.2
##########################
#SYNTAX :
#<pco_data> the pco data to split in bins.
#<tree> a 'phylo' object. The tree must be dated.
#<bins> a series of bins breaks limits.
#<FAD_LAD> a data.frame containing the first and last apparition datums. If none is provided, or if taxa are missing, taxa are assumed to have the same FAD and LAD.
#<include.nodes> logical, whether to include nodes or not in the bins. default = FALSE. If TRUE, the nodes must be the same name in the pco_data and in the tree.
##########################
#Update: fixed FAD_LAD to me more plastic: if input FAD_LAD contains extra taxa, they are now being discarded from the analysis.
#----
#guillert(at)tcd.ie 12/03/2014
##########################

bin.pco<-function(pco_data, tree, bins, include.nodes=FALSE, FAD_LAD) {

    #SANITIZING
    #pco
    check.class(pco_data, 'matrix', ' must be a pco scores matrix.')

    #tree
    check.class(tree, 'phylo', ' must be a phylo object.')
    #the tree must be dated
    if(length(tree$root.time)==0){
        stop("Tree must be a dated tree with $root.time.")
    }

    #bins
    check.class(bins, 'numeric', ' must be numeric.')
    #length must be greater than one
    if(length(bins) < 2) {
        stop("At least two breaks should be specified for the bins.")
    }

    #FAD_LAD
    if(missing(FAD_LAD)) {
        #Create the FAD_LAD table
        FAD_LAD<-data.frame("FAD"=tree.age(tree)[1:Ntip(tree),1], "LAD"=tree.age(tree)[1:Ntip(tree),1], row.names=tree.age(tree)[1:Ntip(tree),2])
        message("No FAD_LAD table has been provided so every tip is assumed to bin single points in time.")
    } else {
        #Check if the FAD_LAD contains all taxa
        if(any(tree$tip.label %in% rownames(FAD_LAD) == FALSE)) {
            message("Some tips have FAD/LAD and are assumed to bin single points in time.")
            #If not generate the FAD_LAD for the missing taxa
            missing_FADLAD<-which(is.na(match(tree$tip.label, rownames(FAD_LAD))))
            add_FAD_LAD<-data.frame(tree.age(tree)[missing_FADLAD,1], tree.age(tree)[missing_FADLAD,1], row.names=tree.age(tree)[missing_FADLAD,2])
            colnames(add_FAD_LAD)<-colnames(FAD_LAD)
            FAD_LAD<-rbind(FAD_LAD, add_FAD_LAD)
        }
        #Remove FAD_LAD taxa not present in the tree
        if(nrow(FAD_LAD) != Ntip(tree)) {
            FAD_LAD<-FAD_LAD[-c(which(is.na(match(rownames(FAD_LAD), tree$tip.label)))),]
        }

    }

    #include.nodes
    check.class(include.nodes, 'logical', " must be logical.")
    #Check if nodes are present in the pco_data object and in the tree
    if(include.nodes == TRUE) {
        #Check if node labels are present
        if(length(tree$node.label) == 0) {
            stop("Provided tree has no nodes labels.")
        } else {
            for (node in 1:length(tree$node.label)) {
                if(length(grep(tree$node.label[node], rownames(pco_data))) == 0) {
                    stop(paste("node", tree$node.label[node], "not found."))
                }
            }
        }

        #Check if the pco_data contains more nodes
        if(nrow(pco_data) > (Ntip(tree)+Nnode(tree))) {
            message("Some rows in pco_data are not present in the tree!")
        }

    } else {
        #Check if the pco_data contains more nodes
        if(nrow(pco_data) > (Ntip(tree)+Nnode(tree))) {
            message("Some rows in pco_data are not present in the tree!")
        }
    }

    #BINING THE PCO
    #ages of tips/nodes + FAD/LAD
    ages_tree_FAD<-tree.age(tree)
    ages_tree_LAD<-tree.age(tree)
    #Change the age if FAD or LAD are higher/lower than the age of the tip
    for(tip in 1:nrow(FAD_LAD)) {
        #Replace age of the tip if FAD is higher
        if(FAD_LAD[tip,1] > ages_tree_FAD$ages[which(ages_tree_FAD$edges == rownames(FAD_LAD)[tip])]) {
            ages_tree_FAD$ages[which(ages_tree_FAD$edges == rownames(FAD_LAD)[tip])]<-FAD_LAD[tip,1]
        }
        #Replace age of the tip if LAD is lower
        if(FAD_LAD[tip,2] < ages_tree_LAD$ages[which(ages_tree_LAD$edges == rownames(FAD_LAD)[tip])]) {
            ages_tree_LAD$ages[which(ages_tree_LAD$edges == rownames(FAD_LAD)[tip])]<-FAD_LAD[tip,2]
        }
    }

    #Empty list element per bin
    bin_elements<-NULL
    bin_elements<-list()

    #Attribute each taxa/node to it's bin
    for (bin in 1:(length(bins)-1)) {
        #Select the elements of one bin
        bin_elements[[bin]]<-ages_tree_FAD$edges[which(ages_tree_FAD$ages >= bins[bin+1] & ages_tree_LAD$ages <= bins[bin])]
    }
    
    #Remove the nodes (if necessary)
    if(include.nodes==FALSE) {
        for (bin in 1:length(bin_elements)) {
        #Remove nomatch with tree$tip.label
            bin_elements[[bin]]<-bin_elements[[bin]][match(tree$tip.label, bin_elements[[bin]])[-which(is.na(match(tree$tip.label, bin_elements[[bin]])))]]
        }
    }

    #Making the pco binned list
    pco_bins<-NULL
    pco_bins<-list()

    for (bin in 1:length(bin_elements)) {
        #Matching list
        matching<-match(as.character(bin_elements[[bin]]),as.character(rownames(pco_data)))
        #If only one taxa is matching, make sure it's not a vector
        if(length(matching) == 1) {
            pco_bins[[bin]]<-matrix(data=pco_data[matching,], nrow=1)
            rownames(pco_bins[[bin]])<-rownames(pco_data)[matching]
        } else {
            pco_bins[[bin]]<-pco_data[matching,]
        }
    }

    #Naming the bins
    name_list<-NULL
    for(bin in 1:length(bin_elements)) {
        name_list[bin]<-paste(bins[bin], bins[bin+1], sep="-")
    }
    
    #If bin is empty, send warning and delete the bin
    #list of empty bins (empty)
    empty_bins<-NULL
    for (bin in 1:length(pco_bins)) {
        if(length(pco_bins[[bin]]) == 0) {
            #Remove the bin
            empty_bins[bin]<-bin
            #Select the empty bin
            empty_bin<-paste(bins[bin], bins[bin+1], sep="-")
            message("The following bin is empty: ", empty_bin, ".")
        }
    }

    #If any empty bins
    if(!is.null(empty_bins)) {
        #NA removal from empty_bins vector (if any)
        if(any(is.na(empty_bins))) {
            empty_bins<-empty_bins[-which(is.na(empty_bins))]
        }
        #Removing the empty bins
        pco_bins<-pco_bins[c(-empty_bins)]
        #Removing the empty bins names
        name_list<-name_list[-empty_bins]
    }
    
    names(pco_bins)<-name_list

    return(pco_bins)
}