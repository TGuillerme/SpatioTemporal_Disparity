##########################
#Reading landmark data (TPS style)
##########################
#Reading landmark coordinates output by TPS
#v.0.1
##########################
#SYNTAX :
#<data> the path to a csv file.
#
##########################
#----
#guillert(at)tcd.ie - 14/07/2016
##########################
#Requirements:
#-R 3


#This amounted to 47 total points collected: 7 ‘internal’ landmarks (figure 2a), 20 points around the talonid, and 20 points around the trigonid

#######
# Prepare the data
#######

# Reading in raw data
raw_coordinates <- read.csv("../Data/molar_LM_coords_GrossnickleNewham2016.csv", header = F, stringsAsFactors = F)

# Removing the taxa with no coordinates
raw_coordinates <- raw_coordinates[-which(raw_coordinates[,1] == 0),]
# Removing the associated ID and LMs
LMs <- which(raw_coordinates[,1] == "LM=47")
# Which of these have an NA in second position
removed_taxa <- which(is.na(raw_coordinates[LMs+1,2]))
raw_coordinates <- raw_coordinates[-c(LMs[removed_taxa], LMs[removed_taxa]+1),]

# Getting the species names
species_names_raws <- grep("ID=", raw_coordinates[,1])
species_names_list <- as.list(as.character(raw_coordinates[species_names_raws, 1]))
species_names <- unlist(lapply(species_names_list, strsplit, split = "ID="))
species_names <- species_names[-which(species_names == "")]

# Getting the number of specimens
number_specimen <- length(species_names)

# Getting landmarks function
get.landmarks <- function(raw_data, number_of_specimens, shift_size, landmarks) {
    list_out <- list()
    for(spec in 1:number_specimen) {
        #specimen shifter
        shifter <- (shift_size)*(spec-1)
        #storing the internal landmarks in the list
        list_out[[spec]] <- data.matrix(raw_data[(landmarks[1]+shifter):(landmarks[2]+shifter),])
    }
    return(list_out)
}

# Getting the internal landmarks
internal_landmarks <- get.landmarks(raw_coordinates, number_specimen, 49, c(2,8))
talonid_slides <- get.landmarks(raw_coordinates, number_specimen, 49, c(9,28))
trigonid_slides <- get.landmarks(raw_coordinates, number_specimen, 49, c(29,48))

# Adding the specimens names
names(internal_landmarks) <- names(talonid_slides) <- names(trigonid_slides) <- species_names

# Remove the coordinates with only 0s




# Plotting the teeth!
plot.teeth <- function(which_outline, talonid_slides, trigonid_slides, internal_landmarks, legend = TRUE, axis = TRUE) {
    op <- par(bty = "n")
    # Plots the talonid
    if(axis) {
        plot(rbind(talonid_slides[[which_outline]], trigonid_slides[[which_outline]], internal_landmarks[[which_outline]]), col = "white", xlab = "", ylab = "")
    } else {
        plot(rbind(talonid_slides[[which_outline]], trigonid_slides[[which_outline]], internal_landmarks[[which_outline]]), col = "white", xlab = "", ylab = "", xaxt = 'n', yaxt = 'n')
    }
    lines(talonid_slides[[which_outline]], col = "blue")
    # Plots the trigonid
    lines(trigonid_slides[[which_outline]], col = "red")
    # Plots the cusp (landmarks)
    points(internal_landmarks[[which_outline]], pch = 19)
    # Legend
    if(legend) {
        limits <- apply(rbind(talonid_slides[[which_outline]], trigonid_slides[[which_outline]], internal_landmarks[[which_outline]]), 2, range)
        legend(limits[1,1], limits[2,2], legend = c("Talonid", "Trigonid", "Cusps"), lty = 1, pch = 19, col = c("blue", "red", "black"), cex = 0.75, bty = "n")
    }
    par(op)
}

par(mfrow = c(5,7))
for(i in 1:35) {
    plot.teeth(i, talonid_slides, trigonid_slides, internal_landmarks, legend = FALSE, axis = FALSE)
}

rotate.coordinates <- function(coordinates, rotation, mirror) {

    warning("proportions not preserved!")
    rotated_coords <- coordinates

    if(!missing(rotation) && rotation != 0) {
        # rotation_matrix <- matrix(c(cos(rotation), -sin(rotation), sin(rotation), cos(rotation)), 2, 2)
        # rotated_coords <- coordinates %*% rotation_matrix

        rotated_coords[,1] = coordinates[,1] * cos(rotation) - coordinates[,2] * sin(rotation)
        rotated_coords[,2] = coordinates[,2] * cos(rotation) + coordinates[,1] * sin(rotation)
        # x' = x cos f - y sin f
        # y' = y cos f + x sin f

    }
    if(!missing(mirror)) {
        if(mirror == "vertical") {
            rotated_coords[,1] <- rotated_coords[,1]*-1
        }
        if(mirror == "horizontal") {
            rotated_coords[,2] <- rotated_coords[,2]*-1
        }
    }

    return(rotated_coords)
}

plot(rotate.coordinates(bla, mirror = "vertical"), type = "l")

#######
# 2D Geomorph
#######


#######
# 2D Fourier
#######


