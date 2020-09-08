elimCellDups <- function(pointsFrame, baseRaster, longLatFields=c('LONG_WGS84', 'LAT_WGS84'), priority=NULL) {
# elimCellDups For a set of geographic points, this script eliminates all but one point in each cell. An optional priority vector can be passed to the function so that points with higher priority are kept.
#
# ARGUMENTS
#
# pointsFrame
# data frame of points with at least the fields named in longLatFields
#
# baseRaster
# Raster* object against which to compare points data frame
#
# longLatFields
# character list of the names of the variables that have longitude and latitude... longitude must be listed first
#
# priority
# optional, numeric or character vector same length as pointsFrame has rows, with values corresponding to priority in which records in the same cell will be kept (though only one point per cell will be kept)... these MUST be specified in alphanumeric order, so that, for example, a record with a priority value of '1' will be preferred over one with a value of '2', which will be preferred over a value of 'alpha' which is preferred over 'beta'
#
# VALUES
#
# cleanedData
# data frame with just one point per cell in baseRaster
#
# BAUHAUS
# - for all points get cell numbers in which they fall in baseRaster
# - adjoin pointsFrame and cell number vector
# - sort by priority if extant
# - eliminate all but one point with duplicate cell numbers
#
# EXAMPLE
# theseRecords <- elimCellDups(pointsFrame=data.frame(LAT_WGS84=c(38, 38.001, 39, 40, 40.01, 38), LONG_WGS84=c(-90, -90, -90, -90, -90, -90), Extra=letters[1:6] ), longLatFields=c('LONG_WGS84', 'LAT_WGS84'), baseRaster=raster('C:/ecology/DropBox/r/areaRasters/areaRaster_easternUsgsProvinces_terrestrialCellsEqual1.tif'), Priority=c(2, 2, 2, 2, 2, 1) )
#
# SOURCE	source('C:/ecology/Drive/Workshops/SDM from Start to Finish (KSU, 2016-02)/Scripts/Eliminate Points in Same Cell of a Raster.r')
#
# LICENSE
# This document is copyright ©2011 by Adam B. Smith.  This document is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3, or (at your option) any later version.  This document is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. Copies of the GNU General Public License versions are available at http://www.R-project.org/Licenses/
#
# AUTHOR	Adam B. Smith | Missouri Botanical Garden | adamDOTsmithATmobotDOTorg
# DATE		2012-01 Adapted from "Utility - Eliminate Species Records in Same Cell of a Raster.r"
# REVISIONS 2012-01 Now returns empty data frame if input frame had no rows

#############################
## libraries and functions ##
#############################

require(raster)

##########################
## initialize variables ##
##########################

##########
## MAIN ##
##########

if (nrow(pointsFrame) > 0) { # frame has at least one row

	# get cell numbers for each point and adjoin with data frame
	pointsFrame$cellNoTEMP <- cellFromXY(object=baseRaster, xy=cbind(pointsFrame[ , longLatFields[1] ], pointsFrame[ , longLatFields[2] ]) )

	# remember original row names
	pointsFrame$origRowNamesTEMP <- rownames(pointsFrame)

	# sort by priority
	if (!is.null(priority)) pointsFrame <- pointsFrame[ order(priority), ]

	# sort by cell number
	pointsFrame <- pointsFrame[ order(pointsFrame$cellNoTEMP), ]

	# reassign rownames so that the "unique" function keeps top-most row
	rownames(pointsFrame) <- 1:nrow(pointsFrame)

	# get top-most point in each cell
	cleanFrame <- pointsFrame[ rownames( unique( as.data.frame(pointsFrame$cellNoTEMP) ) ), ]

	# re-assign original row names
	rownames(cleanFrame) <- cleanFrame$origRowNamesTEMP

	# remove column with original row names
	cleanFrame$origRowNamesTEMP <- NULL
	
	# remove cell number column
	cleanFrame$cellNoTEMP <- NULL

	# original order
	cleanFrame <- cleanFrame[ order( rownames(cleanFrame) ), ]

} else { # frame was empty!

	cleanFrame <- data.frame()
	
}

return(cleanFrame)

}

