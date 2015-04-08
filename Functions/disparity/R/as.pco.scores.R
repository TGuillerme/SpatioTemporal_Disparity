##########################
#Creates a pco.scores object
##########################
#Transforms a pcoa object in a pco.scores object: extracts the PC coordinates and adds a taxonomy data
#v0.1
##########################
#SYNTAX :
#<tree> a 'phylo' object
#<pco> a 'pcoa' object
#<n.axis> the number of axis to extract from pco (default=2)
#<taxonomy.list> an optional list of taxonomic nodes links to taxonomic names present in the given tree
#Example" list("<clade_name>"=<node.number>, ...); grades can also be specified by using list("<crade_name>"=c(<first_node>,<last_node>), ...)
##########################
#----
#guillert(at)tcd.ie 13/10/2014
##########################

as.pco.scores<-function(tree, pco, n.axis=2, taxonomy.list) {

    #SANITIZING
    #tree
    check.class(tree, 'phylo', ' must be a \"phylo\" object.')
    #must have node labels
    if(is.null(tree$node.label)) {
        stop('The tree must have node label names.')
    }

    #pco
    check.class(pco, 'pcoa', ' must be a \"pcoa\" object.')

    #n.axis
    check.class(n.axis, 'numeric', ' must be a single numeric value.')
    check.length(n.axis, 1, ' must be a single numeric value.')
    if(n.axis < 2) {
        stop("At least two pco axis must be selected.")
    }

    #taxonomy.list
    if(missing(taxonomy.list)) {
        do.taxonomy<-FALSE
    } else {
        do.taxonomy<-TRUE
        check.class(taxonomy.list, 'list', ' must be a list of taxonomic names and associated nodes.')
        #Check if nodes of taxonomy.list are present in the tree
        if(!all(match(as.vector(unlist(taxonomy.list)), tree$edge))) {
            stop("Nodes specified in \"taxonomy.list\" are not found in thee given tree.")
        }
    }

    #CREATING THE PCO.SCORES OBJECT
    #Extracting the pco axis
    pco.scores<-pco$vectors[,1:n.axis]

    #Creating the taxonomic table
    if(do.taxonomy==TRUE) {
        taxonomy<-data.frame(row.names=c(tree$tip.label, tree$node.label), 'group'=rep(NA, length(c(tree$tip.label, tree$node.label))))
        #Attributing a group name if absent
        if(length(names(taxonomy.list[1])) == 0) {
            group.name<-"Group1"
        } else {
            group.name<-names(taxonomy.list[1])
        }
        #Checking if one (=clade) or two (=grade) nodes numbers are associated to Group1
        if(length(taxonomy.list[[1]]) == 1) {
            group.type<-"clade"
        } else {
            group.type<-"grade"
        }

        #Attributing "Group1" to the taxa names of group1
        if(group.type == 'clade') {
            taxa<-set.group_clade(tree, taxonomy.list[[1]])
        } else {
            taxa<-set.group_grade(tree, taxonomy.list[[1]])
        }
        taxonomy[taxa,1]<-group.name

        if(length(taxonomy.list) > 1) {
            for (tax.group in 2:length(taxonomy.list)) {
                #Attributing a group name if absent
                if(is.null(names(taxonomy.list[tax.group]))) {
                    group.name<-paste("Group", tax.group, sep="")
                } else {
                    group.name<-names(taxonomy.list[tax.group])
                }
                #Checking if one (=clade) or two (=grade) nodes numbers are associated to Group1
                if(length(taxonomy.list[[tax.group]]) == 1) {
                    group.type<-"clade"
                } else {
                    group.type<-"grade"
                }

                #Attributing "Groupn" to the taxa names of groupn
                if(group.type == 'clade') {
                    taxa<-set.group_clade(tree, taxonomy.list[[tax.group]])
                } else {
                    taxa<-set.group_grade(tree, taxonomy.list[[tax.group]])
                }
                taxonomy[taxa,1]<-group.name

            }
        }
    }

    #OUTPUT (pco.scores object)
    if(do.taxonomy==TRUE) {
        out<-list("scores"=pco.scores, "taxonomy"=taxonomy)
    } else {
        out<-list("scores"=pco.scores)
    }
    class(out)<-"pco.scores"

    return(out)
#End
}