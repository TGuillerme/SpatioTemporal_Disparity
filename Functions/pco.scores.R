##########################
#Creates a pco.scores object
##########################
#Transforms a pcoa object in a pco.scores object: extracts the PC coordinates and adds a taxonomy column
#v0.1
##########################
#SYNTAX :
#<tree> a 'phylo' object
#<pco> a 'pcoa' object
#<axis> the number of axis to extract from pco
#<tax.col> the name of the column containing the taxonomy (default is 'taxonomy')
#<tax.data> the taxonomic data as a set.group object
##########################
#----
#guillert(at)tcd.ie 09/10/2014
##########################