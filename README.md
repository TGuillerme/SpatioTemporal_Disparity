# Mammalian morphological diversity does not increase in response to the Cretaceous-Paleogene mass extinction.
[Thomas Guillerme](http://tguillerme.github.io) and [Natalie Cooper](http://nhcooper123.github.io/).

This repository contains all the code and data used in the manuscript.
###### Manuscript in review in Current Biology.
Part of the code used in the analysis is based on [Graeme T Lloyd](http://graemetlloyd.com/)'s [Claddis](https://github.com/graemetlloyd/Claddis) package. Make sure to have a look, it's great!

Also note that all the code used in this analysis is now properly implemented in the [dispRity](https://github.com/TGuillerme/dispRity) package.

-------

To load all the functions:
------------------
```r
#install.packages("devtools")
library(devtools)
install_github("TGuillerme/SpatioTemporal_Disparity/Functions/disparity")
library(disparity)
```
Note that the following "package" is not a *real* package, just a set of functions to be used to reproduce this paper's analysis. I kept the code as it was for the analysis for reproducibility standards but please do check out the [dispRity](https://github.com/TGuillerme/dispRity) package for the proper code implementatation if you want to use some of the features of this analysis on your own data!

## Data <a href="http://figshare.com/articles/Mammalian_morphological_diversity_does_not_increase_in_response_to_the_Cretaceous_Paleogene_mass_extinction_and_the_extinction_of_the_non_avian_dinosaurs_/1539545"><img src="http://tguillerme.github.io/images/logo-FS.png" height="26" widht="26"/></a> 
All the raw data (trees, cladistic matrices and FAD/LAD) are available in the [data folder](https://github.com/TGuillerme/SpatioTemporal_Disparity/tree/master/Data) on this repository . Both the raw data and the processed data (distance matrices and disparity measurements) are also avaialble on [figshare](http://figshare.com/articles/Mammalian_morphological_diversity_does_not_increase_in_response_to_the_Cretaceous_Paleogene_mass_extinction_and_the_extinction_of_the_non_avian_dinosaurs_/1539545). 

## Analysis

The analysis is divided into four steps
* 1.Estimating the ancestral characters states and the distance matrix for each cladistic matrix (based on [Claddis](https://cran.r-project.org/web/packages/Claddis/) `R` packages).
* 2.Calculating disparity
* 3.Plotting the results
* 4.Testing the differences in disparity around the K-Pg boundary

Note that the two first steps take some time (several hours on your usual desktop machine) so the results from these two steps is made directly available [here](http://figshare.com/articles/Mammalian_morphological_diversity_does_not_increase_in_response_to_the_Cretaceous_Paleogene_mass_extinction_and_the_extinction_of_the_non_avian_dinosaurs_/1539545).
Note also that the three last steps are generalised in the [dispRity](https://github.com/TGuillerme/dispRity) package (still in development though).

#### 1-Estimating the ancestral characters states and the distance matrix
The code is available separately for the [Mammaliaformes](https://github.com/TGuillerme/SpatioTemporal_Disparity/blob/master/Analysis/Data_setups/Data_setup_Slater_claddis.R) and the [Eutheria](https://github.com/TGuillerme/SpatioTemporal_Disparity/blob/master/Analysis/Data_setups/Data_setup_Beck_Claddis.R) datasets. The scripts are divided in four sections:
* **Loading the data.**
* **Cleaning the matrices and the trees:** removing the species with only missing data (i.e. the living species with only molecular data) and making sure the cladistic matrix and the tree matches exactly (label wise).
* **Ancestral states reconstruction:** using a modified version of the `Claddis::AncStateEstMatrix` function (details [here](https://github.com/TGuillerme/SpatioTemporal_Disparity/blob/master/Functions/disparity/R/anc.state.R) and [here (internal)](https://github.com/TGuillerme/SpatioTemporal_Disparity/blob/master/Functions/disparity/R/anc.state_fun.R) that saves the scaled likelihood values for each reconstruction. This allows to get rid of all the characters < 0.95 scaled likelihood (using [anc.unc](https://github.com/TGuillerme/SpatioTemporal_Disparity/blob/master/Functions/disparity/R/anc.unc.R)).
* **Distance matrix:** calculating the Gower distance using the `Claddis::MorphDistMatrix` function.

#### 2-Calculating disparity
All the different scripts for calculating disparity (i.e. different time-slicing methods, time binning method, different disparity metrics) are available in [this folder](https://github.com/TGuillerme/SpatioTemporal_Disparity/tree/master/Analysis/Disparity_calculations). Files are named as follows:
```
Dataset_disparity_nodes_timeSeries_timeSlicing.R
```
with `Dataset` being either `Slater2013` (Mammaliaformes) or `Beck2014` (Eutheria); `nodes` being either `nodes95` for matrices including only nodes with > 0.95 scaled likelihood ancestral reconstructions or just `nodes` for the matrices including all ancestral reconstructions; `timeSeries` being either `sli` for time-slicing method or `int` for time-bining method; and `timeSlicing` (if `sli`) for the time-slicing method (`acc` = ACCTRAN, `del` = DELTRAN, `ran` = punctuated (random) and `pro` = gradual (proximity)).

Each file does:
* **Loading the data.**
* **Cleaning the distance matrix:** using the `Claddis::TrimMorphDistMatrix` function for removing taxa with no overlaping data
* **Ordinating the matrix:** using the `stats::cmdscale` function for ordinating the data
* **Creating the time series:** (or intervals) for generating the time series using the [slice.pco](https://github.com/TGuillerme/SpatioTemporal_Disparity/blob/master/Functions/disparity/R/slice.pco.R) or [int.pco](https://github.com/TGuillerme/SpatioTemporal_Disparity/blob/master/Functions/disparity/R/int.pco.R) functions
* **Calculating disparity:** using the [time.disparity](https://github.com/TGuillerme/SpatioTemporal_Disparity/blob/master/Functions/disparity/R/time.disparity.R) function.

#### 3-Plotting the results
Following the two previous steps, the data can be plotted via [this script](https://github.com/TGuillerme/SpatioTemporal_Disparity/blob/master/Analysis/Disparity_analysis.R).

#### 4-Testing the differences in disparity around the K-Pg boundary
Independently, the differences between the last slice of the Cretaceous (70 Ma) and the slices in the Cenozoic up to 35 Ma can be analysed via [this script](https://github.com/TGuillerme/SpatioTemporal_Disparity/blob/master/Analysis/Disparity_statistics.R).

#### Supplementary analysis
Three supplementary scripts allow to generate the tables and the figures from the appendices:
* [Full analysis for the Mammaliaformes dataset](https://github.com/TGuillerme/SpatioTemporal_Disparity/blob/master/Analysis/Supplementary_Slater_full.R)
* [Plotting the results using all the disparity metrics](https://github.com/TGuillerme/SpatioTemporal_Disparity/blob/master/Analysis/Supplementary_results.R)
* [Testing the differences around the K-Pg boundary using all the disparity metrics](https://github.com/TGuillerme/SpatioTemporal_Disparity/blob/master/Analysis/Supplementary_results.R)
