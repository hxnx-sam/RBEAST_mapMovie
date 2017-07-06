# function to make succession of maps images
# extract of earlier code, h5n8_global_consortium_movie_only.R (15 Apr 2016 - 8 May 2017)
# S J Lycett
# 4 July 2017
# 6 July 2017

# load required libraries
library(maps)
library(mapdata)
library(mapproj)
library(RgoogleMaps)
library(geosphere)

# load utility functions (assumes these are in same directory as mapMovie.R and this is the current working directory)
#Rpath <- ""
#source(paste(Rpath,"get_BEAST_cols.R",sep=""))
#source(paste(Rpath,"calcDecimalDate.R",sep=""))
#source(paste(Rpath,"read_BEAST_tree_latlon.R",sep=""))
#source(paste(Rpath,"makeGiantTbl.R",sep=""))

########################################################################################
# MAIN FUNCTION
########################################################################################

# plot a sample of the phylogeographic construction on a centre north map
# points are coloured by trait (Host)
# requires 'giant_tbl' input data
# assumes that the coordinates in giant_tbl are already converted for plotting
# assumes that you want to plot point colours according to a discrete trait, and that this trait is also in the giant_tbl (here Host)
# generates a series of png images, one for each time point
northMapMovieImages <- function( tblName=tblName, trait="Host",  ustateNames=c(-1),
						fillCol="grey90", npts=1500, res=150, w=1920/res, h=1080/res,
						cex=1, inc=7/365, startTime=2013.5, endTime=2015.3, 
						showNorth=TRUE, showLegend=TRUE, 
						showTitle=TRUE, titleTxt="Infection route from virus sequence data",
						useAllReps=FALSE,
						imagePath=imagePath, imName=paste(trait,"_latlon_movie_",sep="") ) {
											
	projection 		<- "azequalarea"
	orientation		<- c(90,0,0)
	ylim	     		<- c(-40,89)
	
	# load giant table
	print(paste("Loading giant table:",tblName))
	load( tblName )
	
	cn	   <- colnames(giant_tbl)
	to_j   <- which( cn == paste("to",trait,sep="") )
	from_j <- which( cn == paste("from",trait,sep="") )
	
	uhosts <- levels( giant_tbl[,to_j] )	#levels(giant_tbl$toHost)
	nhosts <- length(uhosts)
	
	# option to use user supplied state names instead of the names in the giant_tbl
	# e.g. these, instead of Dom-ans, Dom-gal, Wild-Long, Wild-Short etc
	# ustateNames = c("Domestic Ducks & Geese","Chickens & Turkeys","Wild Long Range Migrators","Other Wild Birds")
	if (ustateNames[1]==-1) {
		ustateNames=uhosts
	}
	
	hCols 	<- get_BEAST_cols(nhosts, bright=0.8, sat=0.8)
	
	timePts    	<- seq(startTime,endTime,inc)
	midTime    	<- (giant_tbl$toTime + giant_tbl$fromTime)/2
	
	# these are called Lat and Lon but are actually not, these are really x and y
	midLats		<- (giant_tbl$toLat + giant_tbl$fromLat)/2
	midLons		<- (giant_tbl$toLon + giant_tbl$fromLon)/2

	toHostCols	<- hCols[match( giant_tbl[,to_j], 	uhosts)]
	fromHostCols<- hCols[match( giant_tbl[,from_j], uhosts)]
	
	count		<- 100
	
	for (i in 1:length(timePts)) {
		
			# lineages which exist at this time
			is	<- timePts[i]-inc
			ie 	<- timePts[i]+inc

			if (useAllReps) {
				tinds1	<- which( ( midTime < ie ) & ( midTime >= is) & 
						    (giant_tbl$fromTime <= timePts[i]) & (giant_tbl$toTime >= timePts[i])  )

			} else {
				tinds1	<- which( ( midTime < ie ) & ( midTime >= is) & (giant_tbl$repInd==1) &
						    (giant_tbl$fromTime <= timePts[i]) & (giant_tbl$toTime >= timePts[i])  )
			}

			if (length(tinds1) > npts) {
				tinds1  <- sample(tinds1, npts)
			}
			fract <- (giant_tbl$toTime[tinds1]-timePts[i])/( giant_tbl$toTime[tinds1]-giant_tbl$fromTime[tinds1] )

			tempHostCol <- fromHostCols[tinds1]
			kk		<- which(fract > 0.5)
			if (length(kk) > 0) {
				tempHostCol[kk] <- toHostCols[tinds1[kk]]
			}
			 			
			imageName <- paste(imagePath,imName,count,".png",sep="")
			png(file=imageName, width=w*res, height=h*res, res=res)

			op	<- par(mar=c(1,1,1,1))
				map("worldHires", projection=projection,orientation=orientation, 
					ylim=ylim, fill=TRUE, col="white", border=NA, bg=fillCol)
			
				points( midLats[tinds1], midLons[tinds1], col=tempHostCol, cex=cex )

				# add north
				if (showNorth) {
					points(0, 0, pch=3, bg="black", cex=cex)
				}
			
				if (showLegend) {
					legend("bottomright",paste(ustateNames),pch=21,pt.bg=hCols,bty="n")
					legend("bottomleft",invertDecimalDate(timePts[i],ddmmyy=TRUE),pch=NA,bty="n")
				}

				if (showTitle) {
					legend("topleft",titleTxt,pch=NA, bty="n")
				}

			par(op)

			dev.off()
			

			count <- count + 1
		}	

		print(paste("Done",(count-101),"images"))

}



# plot a sample of the phylogeographic construction on a map
# points are coloured by trait (Host)
# requires 'giant_tbl' input data
# generates a series of png images, one for each time point
mapMovieImages <- function( tblName=tblName, trait="Host",  ustateNames=c(-1),
					fillCol="grey90", npts=1500, res=150, w=1920/res, h=1080/res,
					cex=1, inc=7/365, startTime=2013.5, endTime=2015.3, useAllReps=FALSE,
					showNorth=TRUE, showLegend=TRUE, 
					showTitle=TRUE, titleTxt="Infection route from virus sequence data",
					convertCoords=TRUE,
					projection="azequalarea", orientation=c(90,0,0), ylim=c(-40,89), xlim=c(-180,180),
					imagePath=imagePath, imName=paste(trait,"_latlon_movie_",sep="") ) {
											
	
	# load giant table
	print(paste("Loading giant table:",tblName))
	load( tblName )
	
	cn	   <- colnames(giant_tbl)
	to_j   <- which( cn == paste("to",trait,sep="") )
	from_j <- which( cn == paste("from",trait,sep="") )
	
	uhosts <- levels( giant_tbl[,to_j] )	#levels(giant_tbl$toHost)
	nhosts <- length(uhosts)
	
	# option to use user supplied state names instead of the names in the giant_tbl
	# e.g. these, instead of Dom-ans, Dom-gal, Wild-Long, Wild-Short etc
	# ustateNames = c("Domestic Ducks & Geese","Chickens & Turkeys","Wild Long Range Migrators","Other Wild Birds")
	if (ustateNames[1]==-1) {
		ustateNames=uhosts
	}
	
	hCols	 	<- get_BEAST_cols(nhosts, bright=0.8, sat=0.8)
	
	timePts    	<- seq(startTime,endTime,inc)
	midTime    	<- (giant_tbl$toTime + giant_tbl$fromTime)/2
	
	# if convertCoords = TRUE, then these are really Lat and Lons and require converting to the map projection coordinates for plotting
	midLats		<- (giant_tbl$toLat + giant_tbl$fromLat)/2
	midLons		<- (giant_tbl$toLon + giant_tbl$fromLon)/2

	toHostCols	<- hCols[match( giant_tbl[,to_j], 	uhosts)]
	fromHostCols<- hCols[match( giant_tbl[,from_j], uhosts)]
	
	count		<- 100
	
	for (i in 1:length(timePts)) {
		
			# lineages which exist at this time
			is	<- timePts[i]-inc
			ie 	<- timePts[i]+inc

			if (useAllReps) {
				tinds1	<- which( ( midTime < ie ) & ( midTime >= is) & 
						    (giant_tbl$fromTime <= timePts[i]) & (giant_tbl$toTime >= timePts[i])  )

			} else {
				tinds1	<- which( ( midTime < ie ) & ( midTime >= is) & (giant_tbl$repInd==1) &
						    (giant_tbl$fromTime <= timePts[i]) & (giant_tbl$toTime >= timePts[i])  )
			}

			if (length(tinds1) > npts) {
				tinds1  <- sample(tinds1, npts)
			}
			fract <- (giant_tbl$toTime[tinds1]-timePts[i])/( giant_tbl$toTime[tinds1]-giant_tbl$fromTime[tinds1] )

			tempHostCol <- fromHostCols[tinds1]
			kk		<- which(fract > 0.5)
			if (length(kk) > 0) {
				tempHostCol[kk] <- toHostCols[tinds1[kk]]
			}
			 			
			imageName <- paste(imagePath,imName,count,".png",sep="")
			png(file=imageName, width=w*res, height=h*res, res=res)

			op	<- par(mar=c(1,1,1,1))
				map("worldHires", projection=projection,orientation=orientation, 
					ylim=ylim, xlim=xlim, fill=TRUE, col="white", border=NA, bg=fillCol)
					
				if (convertCoords) {
					# convert lat lons to local coordinates for plotting on the map
					coords <- mapproject(midLons[tinds1], midLats[tinds1], projection=projection, orientation=orientation)
					points(coords, col=tempHostCol, cex=cex)
				} else {
					# coordinates were already converted before running BEAST and making giant_tbl
					points( midLats[tinds1], midLons[tinds1], col=tempHostCol, cex=cex )
				}

				# add north
				if (showNorth) {
					points(0, 0, pch=3, bg="black", cex=cex)
				}
			
				if (showLegend) {
					legend("bottomright",paste(ustateNames),pch=21,pt.bg=hCols,bty="n")
					legend("bottomleft",invertDecimalDate(timePts[i],ddmmyy=TRUE),pch=NA,bty="n")
				}

				if (showTitle) {
					legend("topleft",titleTxt,pch=NA, bty="n")
				}

			par(op)

			dev.off()
			

			count <- count + 1
		}


		print(paste("Done",(count-101),"images"))	

}

