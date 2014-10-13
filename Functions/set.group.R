#OBSOLETE FUNCTION see pco.scores()
##########################
#Set a group of taxa in pco plot
##########################
#Extract all the species from a clade or a grade
#v0.2
#Update: can now do multiple groups in a single line
#Update: produces a table of the species taxonomy for pco.scores
#TO DO!
##########################
#SYNTAX :
#<tree> a 'phylo' object
#<pco.names> a list of species from the pco object
#<type> the type of group ('clade' or 'grade')
#<node> the node number of the origin of the clade or two nodes for a grade
#<name> the name of the group
##########################
#----
#guillert(at)tcd.ie 25/09/2014
##########################


set.group<-function(tree, pco.scores, tax.col, type, node, name) {

#SANITYZING

    #tree
    check.class(tree, 'phylo', ' must be a phylo object.')
    #must have node labels
    if(is.null(tree$node.label)) {
        stop('The tree must have node label names.')
    }

    #pco.scores
    #check.class(pco.scores, 'matrix', ' must be a pco axis matrix.')
    #if(!grep('Axis', names(pco.scores))) {
    #    stop('pco.scores must be a pco axis matrix.')
    #}

    #tax.col
    if(missing(tax.col)) {
        tax.col<-'taxonomy'
    } else {
        check.class(tax.col, 'character', ' must be a character string.')
        if(!grep(tax.col, names(pco.scores))) {
            stop('tax.col not found in pco.scores.')
        }
    }

    #type
    check.class(type, 'character', ' must be a \'clade\' or \'grade\'.')
    if(type != 'clade') {
        if (type != 'grade') {
            stop('type must be a \'clade\' or \'grade\'.')
        }
    }

    #node
    check.class(node, 'numeric', ' must be numeric.')
    if(type == 'clade') {
        check.length(node, 1, ' must be a single node.')
    }
    if(type == 'grade') {
        check.length(node, 2, ' must be a two nodes.')
    }

    #name
    check.class(name, 'character', ' must be a single character string.')
    check.length(name, 1, ' must be a single character string.')

#SETING GROUP

    if(type == 'clade') {
        taxa<-set.group_clade(tree, node)
    } else {
        taxa<-set.group_grade(tree, node)
    }

    pco.scores<-set.group_add.tax(pco.scores, tax.col, taxa, name)
}