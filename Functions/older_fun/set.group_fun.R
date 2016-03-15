#FUNCTIONS FOR set.group

#Extracting the tips and the nodes names from a clade
set.group_clade<-function(tree, node) {
    clade<-extract.clade(tree, node)
    names<-c(clade$tip.label, clade$node.label)
    return(names)
}

#Extracting the tips and the nodes names from a grade
set.group_grade<-function(tree, node) {
    clade<-extract.clade(tree, node[1])
    drop<-extract.clade(tree, node[2])$tip.label
    grade<-drop.tip(clade, drop)
    names<-c(grade$tip.label, grade$node.label)
}

#Adding a taxonomic group to a pco.scores
set.group_add.tax<-function(pco.scores, tax.col, taxa, name) {
    if(length(grep(tax.col, colnames(pco.scores))) == 0) {
        tax.column<-ncol(pco.scores)+1
        pco.scores<-cbind(pco.scores, rep(NA, nrow(pco.scores)))
        colnames(pco.scores)[tax.column]<-tax.col
    }
    pco.table<-as.data.frame(pco.scores)
    tax.column<-which(colnames(pco.table) == tax.col)
    pco.table[match(taxa, rownames(pco.table)), tax.column]<-name
    return(pco.table)
}
