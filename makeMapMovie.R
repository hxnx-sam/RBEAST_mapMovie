# Full process to make points movie from BEAST runs
# e.g. for H5N8 data
# S J Lycett
# 6 July 2017

# load utility functions (assumes these are in same directory as makeMapMovie.R)

# set this to the correct directory, or leave as Rpath="" if already in correct directory
Rpath <- "Rcode//mapMovie//"

source(paste(Rpath,"getEl.R",sep=""))
source(paste(Rpath,"get_BEAST_cols.R",sep=""))
source(paste(Rpath,"calcDecimalDate.R",sep=""))
source(paste(Rpath,"read_BEAST_tree_latlon.R",sep=""))
source(paste(Rpath,"makeGiantTbl.R",sep=""))
source(paste(Rpath,"mapMovie.R",sep=""))


#####################################################################################################
# ACTUAL EXAMPLE for H5N8 North Pole Movie
#####################################################################################################


# STEP 1 convert the beast runs to *.trs files
# this will take a long time (e.g. 30 mins)
	path 				<- "D://Movie_Examples//H5N8//"
	rootname 			<- "H5N8_seg4_130_empirical"
	discreteTreesName 	<- paste(path,rootname,"_Host_asym_11.trees.txt",sep="")
	latlonTreesName 		<- paste(path,rootname,"_latlon_hbr_1.trees.txt",sep="")
	trait				<- "Host"
	burnin			<- 1000

	# small test
	nsamples			<- -1		# -1 = do all samples
	step				<- 1		# 1 = no further thinning
	
	trsNames	<- convert_BEAST_runs(	discreteTreesName=discreteTreesName, trait=trait,
					latlonTreesName=latlonTreesName,
					burnin=burnin, nsamples=nsamples, step=step)
						
# STEP 2 make the giant_tbl
	tblName	<- makeGiantTbl(rootname = paste(path,rootname,sep=""),
						r_discreteTreesName = trsNames$r_discreteTreesName,
						r_latlonTreesName = trsNames$r_latlonTreesName,
						nper=1)
	
# STEP 3 make the movie png images
	imagePath 	<- "D://Movie_Examples//H5N8//png//"
	res 		<- 150
	w		<- 1920
	h		<- 1080
	npts		<- 1500
	startTime	<- 2014.0
	endTime	<- 2015.25
	inc		<- 7/365
	fillCol	<- "grey90"
	northMapMovieImages( tblName=tblName, trait=trait, 
					fillCol=fillCol, npts=npts, res=res, w=w/res, h=h/res,
					cex=1, inc=inc, startTime=startTime, endTime=endTime, 
					imagePath=imagePath )

# STEP 4 use an external tool to combine the pngs to a movie file, e.g. gif, mpeg4, avi etc
# for example: ImageMagick
# https://www.r-bloggers.com/animated-plots-with-r/
# From command line on PC, cd to directory, then use this command:
# convert Host_latlon_movie_1*.png -set delay 25 -loop 0 Host_latlon.gif


