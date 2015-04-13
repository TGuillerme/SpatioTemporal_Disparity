# Mass extinction and niche replacement around the K-Pg boundary in mammals.
[Thomas Guillerme](http://tguillerme.github.io) and [Natalie Cooper](https://http://nhcooper123.github.io/).

This repository contains all the code and data used in the manuscript.
###### Manuscript in prep.
Part of the code used in the analysis is based on [Graeme T Lloyd](http://graemetlloyd.com/)'s [Claddis](https://github.com/graemetlloyd/Claddis) package. Make sure to have a look, it's great!

-------

To load all the functions:
------------------
```r
#install.packages("devtools")
library(devtools)
install_github("TGuillerme/SpatioTemporal_Disparity/Functions/disparity")
library(disparity)
```
Note that the following "package" is not a real package, just a set of functions to be use to reproduce this paper's analysis.

## Data
All the data will be available once every analysis is run (hopefully that means 'soon').

## Analysis
Here are the following steps for the running the full analysis (`shell` version only - a thorough `R` version should be coming soon.)
####Estimating ancestral characters and calculating the distance matrices (`shell`)
The [`Data.setup.sh`](https://github.com/TGuillerme/SpatioTemporal_Disparity/blob/master/Functions/Data.setup.sh) script allows to generate the proper `R` script to run. The analysis may take a bit long on certain machines so I advice to run it as a back ground task on a terminal using the following script (as an example):

```
# Generating the R script
sh Data.setup.sh "Beck2014" "../Data/" "../Data/2014-Beck-ProcB-matrix-morpho.nex" "../Data/2014-Beck-ProcB-TEM.tre" "Claddis"

# Launching the R script (as a verbose background task)
R --no-save < Data_setup_Beck_Claddis.R > /dev/null
```

####Calculating the disparity (`shell` version only - a thorough `R` version should be coming soon.)
The [`disparity.calc.sh`](https://github.com/TGuillerme/SpatioTemporal_Disparity/blob/master/Functions/disparity.calc.sh) script generates the `R` script. The analysis also takes some time so same advice as above. Here's a following up example:

```
# Setting the intervals and the slices values
intervals="170,155,140,125,110,95,80,65,50,35,20,0"
slices="170,165,160,155,150,145,140,135,130,125,120,115,110,105,100,95,90,85,80,75,70,65,60,55,50,45,40,35,30,25,20,15,10,5,0"

# Generating the R script
sh disparity.calc.sh "Beck2014" "../Data/" "../Data/2014-Beck-ProcB-matrix-morpho.nex" "../Data/2014-Beck-ProcB-TEM.tre" "../Data/Beck2014/Beck2014Claddis_distance-nodes95.Rda" "gower.dist.matrix" $intervals $slices "../Data/Beck2014_FADLAD.csv"

# Launching the R script (as a verbose background task)
R --no-save < Beck2014-disparity-Claddis_distance-nodes95-gower.dist.matrix.R > /dev/null
```

####Null data (`shell` version only)
In the mean time you can calculate the expected null disparity using the [`disparity.null.sh`](https://github.com/TGuillerme/SpatioTemporal_Disparity/blob/master/Functions/disparity.null.sh) script. Note that this is a single core version, to make faster (up to 24 cores at the moment), use the [cluster version](https://github.com/TGuillerme/SpatioTemporal_Disparity/blob/master/Functions/Cluster/disparity.null.sh) using [`slurm`](https://en.wikipedia.org/wiki/Slurm_Workload_Manager). Here is an example for calcualting the expected null centroid distance disparity for the previous examples under a Mkn model (Brownian):

```
# Setting the intervals and the slices values
intervals="170,155,140,125,110,95,80,65,50,35,20,0"
slices="170,165,160,155,150,145,140,135,130,125,120,115,110,105,100,95,90,85,80,75,70,65,60,55,50,45,40,35,30,25,20,15,10,5,0"

#Generating the R script
sh disparity.null.sh "Beck2014" "../Data/" "../Data/2014-Beck-ProcB-matrix-morpho.nex" "../Data/2014-Beck-ProcB-TEM.tre" "centroid" "sim.char" $intervals $slices 

# Launching the R script (as a verbose background task)
R --no-save < Beck2014-sim.char-centroid.R > /dev/null
```

###### More details and update are about to come as we complete the study.
