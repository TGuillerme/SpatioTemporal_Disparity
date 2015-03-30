#!/bin/sh
##########################
#Shell script for setting up the data for disparity analysis.
##########################
#SYNTAX:
#sh Data.setup.sh <chain> <method> <moduleslist> <split>
#with:
#<chain> the name of the chain to generate task files for
#<path> the path where the data will be stored under the chain name
#<matrix> the path to the morphological matrix (can be relative)
#<tree> the path to the phylogenetic tree (can be relative)
#<ace> either "TRUE" to calculate it (two tasks) or the path to an ancestral state matrix. WARNING, if distance is set to "TRUE" only the tip distance can be calculated simultaneously.
#########################
#version 0.1
Data.setup_version="Data.setup v0.1"
#----
#guillert(at)tcd.ie - 30/03/2015
###########################

#INPUT
chain=$1
path=$2
moduleslist=$3
split=$4
