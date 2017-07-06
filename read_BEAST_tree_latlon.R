# convert BEAST trees to R object
# S J Lycett
# 10 July 2015 - see h5n8_movie.R
# 28 Oct 2015 - processLatLonTree_oct2015
# 30 Aug 2016 - option to remove rate from tree lines (required if not empirical trees)

library(ape)
library(geosphere)


##############################################################################################################
# START INTERNAL CUSTOM FUNCTIONS - DO NOT MODIFY

####################################################################
# see readBeastLatLonTrees.R

# used to get the numbers -> names translation for the taxa
getTranslation	<- function( lines ) {
	ts 		<- grep("Translate", lines)
	te 		<- grep(";", lines)
	te 		<- te[which( te > ts )][1]
	taxaLines 	<- lines[ (ts+1) : (te-1) ]

	taxaLines   <- gsub("'", "", taxaLines)

	#taxa		<- apply(as.matrix(taxaLines), 1, getEl, sep="'", ind=2)
	#tipNumbers	<- apply(as.matrix(taxaLines), 1, getEl, sep="'", ind=1)

	taxa		<- apply(as.matrix(taxaLines), 1, getEl, sep=" ", ex=1, reconstruct=TRUE)
	tipNumbers  <- apply(as.matrix(taxaLines), 1, getEl, sep=" ", ind=1)

	taxa		<- gsub(",", "", taxa)
	tipNumbers  <- gsub("\t", "", tipNumbers)
	tipNumbers  <- gsub(" ", "", tipNumbers)

	translateTbl <- cbind(tipNumbers, taxa)

	return( translateTbl )
}

# internal function to get latitude and longitude from node name
getLatLon	<- function( txt ) {
	#latlonRegex <- "latlon=\\{[\\-]?[0-9]+\\.[0-9]+\\|[\\-]?[0-9]+\\.[0-9]+\\}"
	latlonRegex <- "latlon=\\{[\\-]?[0-9]+\\.[0-9Ee\\-]+\\|[\\-]?[0-9]+\\.[0-9Ee\\-]+\\}"
	llmatch	<- gregexpr(latlonRegex,txt)
	if (llmatch[[1]][1]==-1) {
		# presumed NaN's
		latlon <- c(0,0)

		latlonRegex <- "latlon=\\{NaN\\|NaN\\}"
		llmatch	<- gregexpr(latlonRegex,txt)

		is		<- llmatch[[1]]
		ie		<- llmatch[[1]] + attributes(llmatch[[1]])$match.length - 1
		old_txt	<- substring(txt, is, ie)
		for (k in 1:length(old_txt)) {
			txt <- gsub(old_txt[k],"",txt,fixed=TRUE)
		}
		old_txt	<- gsub("latlon=\\{","",old_txt)
		old_txt	<- gsub("\\}","",old_txt)

	} else {
		is		<- llmatch[[1]]
		ie		<- llmatch[[1]] + attributes(llmatch[[1]])$match.length - 1
		old_txt	<- substring(txt, is, ie)
		for (k in 1:length(old_txt)) {
			txt <- gsub(old_txt[k],"",txt,fixed=TRUE)
		}
		old_txt	<- gsub("latlon=\\{","",old_txt)
		old_txt	<- gsub("\\}","",old_txt)
		latlon	<- as.numeric(strsplit(old_txt,"\\|")[[1]])
	}
	
	return( list( txt=txt, latlon=latlon) )
}

# internal function to extract properties from node name
extractNodeProperties <- function( txt ) {

	edgeLen	<- as.numeric(strsplit(txt,"\\]")[[1]][2])
	txt		<- strsplit(txt,"\\]")[[1]][1]
	res 		<- getLatLon(txt)
	latlon 	<- res$latlon	
	remTxt	<- res$txt
	els		<- strsplit(remTxt,"\\|")[[1]]
	#inds		<- which( apply(as.matrix(els), 1, nchar) > 0 )
	inds		<- which( apply(as.matrix(els), 1, nchar) > 1 )

	if (length(inds) > 0) {
		props		<- els[inds]
		props		<- gsub("\\&","",props)
		props		<- gsub("\"","",props)
		props		<- gsub("--","",props)

		props		<- t(matrix(unlist(strsplit(props,"=")),2,length(props)))
		ii		<- grep("rate",props[,1])
		edgeProps	<- props[ii,]
		props		<- props[ setdiff(1:length(props[,1]), ii), ]
	} else {
		edgeProps	<- c()
		props		<- c()
	}

	return( list(latlon=latlon, props=props, edgeLen=edgeLen, edgeProps=edgeProps) )
}

# internal function
getEdgeLen_val	<- function( prop ) {
	return (prop$edgeLen )
}

# internal function
getLatLon_val 	<- function( prop ) {
	return( prop$latlon )
}

# internal function
getPropNames_val <- function( prop ) {
	return( prop[,1] )
}

# internal function
getHost_val	<- function( prop, ind=2 ) {
	return( prop$props[ind] )
}

# internal functions, used in readBeastRateTrees - originally from readBeastMultiRateTree.R
# from population_fitness_features.R
# internal function
distFromRoot	<- function( tr ) {

	rootNode	<- length(tr$tip.label)+1
	nodeDists	<- array(0, max(tr$edge))
	
	toProcess	<- c(rootNode)

	while ( length(toProcess) > 0 ) {
		einds 			<- which(tr$edge[,1]==toProcess[1])

		if (length(einds) > 0) {
			children 			<- tr$edge[einds,2]
			nodeDists[children] 	<- nodeDists[toProcess[1]] + tr$edge.length[einds]

			toProcess			<- c(toProcess, children)
		}
		toProcess			<- setdiff(toProcess, toProcess[1])
	}

	tr$nodeDists <- nodeDists

	return ( tr )
}

# internal function
nodeTimes	<- function(tr, youngestTip=2011.027) {

	if ( !any(attributes(tr)$names == "nodeDists") ) {
		tr <- distFromRoot(tr)
	}

	tr$rootHeight<- max(tr$nodeDists)
	tr$nodeTimes <- youngestTip - tr$rootHeight + tr$nodeDists

	return ( tr )
}



# function to process a single BEAST tree line
processLatLonTree <- function( trLine, youngestTip=youngestTip, translateTbl=translateTbl ) {

		#trLine   <- lines[i]
		els	   <- strsplit(trLine, "\\[\\&R\\] ")[[1]]
		tempTr   <- els[2]

		# remove edge properties - attaching to previous child node
		tempTr   <- gsub("\\]\\:\\[\\&","--\\&",tempTr)

		# remove \\. from rate names
		tempTr   <- gsub("\\.rate","-rate",tempTr)

		# replace all , in [] with |
		b_regex	<- "\\[\\&"
		b_is		<- gregexpr(b_regex,tempTr)[[1]]
		#c_regex	<- "[0-9]\\]" 17 may 2015, need to allow }] for hbr empirical
		c_regex	<- "\\]"
		c_ie		<- gregexpr(c_regex,tempTr)[[1]]
		c_ie		<- c(c_ie, nchar(tempTr))
		old_txt	<- substring(tempTr, b_is, c_ie)
		new_txt	<- gsub("\\,","\\|", old_txt)
		for (k in 1:length(old_txt)) {
			tempTr <- gsub(old_txt[k],new_txt[k],tempTr,fixed=TRUE)
		}

		# read into tree object
		tr		<- read.tree(text=tempTr)

		# process tip and node names
		tipnames	<- apply(as.matrix(tr$tip.label), 1, getEl, sep="\\[")
		tipNos	<- tipnames[1,]
		tipProps	<- apply( as.matrix(tipnames[2,]), 1, extractNodeProperties )
				


		nodenames	<- apply(as.matrix(tr$node.label), 1, getEl, sep="\\[")
		nodeProps	<- apply( as.matrix(nodenames[2,]), 1, extractNodeProperties )
		
	
		all_props	<- c(tipProps, nodeProps)
		n		<- length(all_props)
		edgeLen	<- array(0,n)
		latlon	<- matrix(0, n, 2)
		for (k in 1:n) {
			edgeLen[k] <- getEdgeLen_val( all_props[[k]] )
			latlon[k,] <- as.matrix(getLatLon_val( all_props[[k]] ))
		}

		tr$latlon	   <- latlon

		#if ( all(is.finite(edgeLen)) & all(edgeLen != 0) ) {
			toNodes	   <- tr$edge[,2]
			edgeLen_toNodes<- edgeLen[toNodes]
			tr$edge.length <- edgeLen_toNodes
		#}

		tr$tip.label   <- tipNos
		tinds		   <- match(tipNos, translateTbl[,1])
		all( tipNos==translateTbl[tinds,1] )
		tr$tip.label   <- translateTbl[tinds,2]
		tr$node.label  <- NULL

		tr		   <- nodeTimes(tr, youngestTip=youngestTip)

		return( tr )

}

# function to process a single tree line
# for ape 3.3 seems that format has changed..
processLatLonTree_oct2015 <- function(trLine, youngestTip=youngestTip, translateTbl=translateTbl, removeRate=FALSE ) {
	#trLine   <- lines[i]
	els	   <- strsplit(trLine, "\\[\\&R\\] ")[[1]]
	tempTr   <- els[2]

	if (removeRate) {
		tempTr <- gsub("\\[\\&rate=[-]?[0-9]+\\.[0-9Ee\\-]+]","",tempTr)
	}

	latlonregex <- "\\[\\&latlon=\\{[-]?[0-9]+\\.[0-9E\\-]+,[-]?[0-9]+\\.[0-9E\\-]+\\},latlon=\\{[-]?[0-9]+\\.[0-9E\\-]+,[-]?[0-9]+\\.[0-9E\\-]+\\}\\]"
	pos		<-gregexpr(latlonregex,tempTr)

	llvals	<- substring(tempTr, pos[[1]], pos[[1]]+attributes(pos[[1]])$match.length-1)
	
	nn		<- length(translateTbl[,1])
	nn		<- 2*nn-1
   if (length(llvals) != nn) {
		#stop("processLatLonTree_oct2015 - didnt recognise all the lat-lon labels")
		print("processLatLonTree_oct2015 - didnt recognise all the lat-lon labels")
		return( NULL )
   } else {

	lats		<- apply(as.matrix(llvals), 1, getEl, ind=1, sep=",")
	lats		<- apply(as.matrix(lats), 1, getEl, ind=2, sep="=\\{")
	lons		<- apply(as.matrix(llvals), 1, getEl, ind=2, sep=",")
	lons		<- gsub("\\}", "", lons)
	newllstr    <- paste("lat_",lats,"-lon_",lons,sep="")

	tempTr2	<- tempTr
	for (i in 1:length(llvals)) {
		tempTr2 <- gsub(llvals[i],  newllstr[i], tempTr2, fixed=TRUE)
	}
	
	# read into tree object
	tr		<- read.tree(text=tempTr2)

	# get lats and lons for tips
	tipNames	<- apply(as.matrix(tr$tip.label), 1, getEl, sep="lat_")	
	tipProps	<- tipNames[2,]
	tipNames	<- tipNames[1,]

	tiplats	<- apply(as.matrix(tipProps), 1, getEl, sep="-lon_")
	tiplons	<- tiplats[2,]
	tiplats	<- tiplats[1,]

	# get lats and lons for nodes
	nodelats	<- apply(as.matrix(tr$node.label), 1, getEl, sep="-lon_")
	nodelons	<- nodelats[2,]
	nodelats	<- gsub("lat_","",nodelats[1,])

	# re-attach lats and lons as attributes in tree, ordered same as tips and nodes
	lat		<- as.numeric(c(tiplats, nodelats))
	lon		<- as.numeric(c(tiplons, nodelons))
	latlon	<- cbind(lat,lon)
	tr$latlon	<- latlon

	tinds		   <- match(tipNames, translateTbl[,1])
	all( tipNames==translateTbl[tinds,1] )
	tr$tip.label   <- translateTbl[tinds,2]
	tr$node.label  <- NULL

	tr		   <- nodeTimes(tr, youngestTip=youngestTip)

	return( tr )

    }

}

# for one trait only, for empirical trees
processDiscreteTree_oct2015 <- function(trLine, youngestTip=youngestTip, 
						translateTbl=translateTbl, trait="Host" ) {
	#trLine   <- lines[i]
	els	   <- strsplit(trLine, "\\[\\&R\\] ")[[1]]
	tempTr   <- els[2]

	#remove Host rate
	#hr_regex <- "\\[\\&Host\\.rate=[-]?[0-9]+\\.[0-9]+\\]"
	hr_regex <- "\\[\\&Host\\.rate=[-]?[0-9]+\\.[0-9]+[Ee]?[-]?[0-9]?\\]"
	if (trait != "Host") {
		hr_regex <- gsub("Host",trait, hr_regex)
	}
	tempTr   <- gsub(hr_regex, "", tempTr)

	#make node annotations part of node
	tempTr   <- gsub("\\[\\&","__",tempTr)
	tempTr   <- gsub("\\]","",tempTr)

	tr	   <- read.tree(text=tempTr)
	tipNames <- tr$tip.label
	tipProps <- apply(as.matrix(tipNames), 1, getEl, sep=paste("__",trait,"=",sep=""))
	tipNames <- tipProps[1,]
	tipProps <- gsub("\"","",tipProps[2,])

	nodeProps<- gsub(paste("__",trait,"=",sep=""), "", tr$node.label)
	nodeProps<- gsub("\"", "", nodeProps)

	tr$state 		<- c(tipProps, nodeProps)
	tr$tip.state 	<- tipProps
	tr$node.state	<- nodeProps

	tinds		   <- match(tipNames, translateTbl[,1])
	all( tipNames==translateTbl[tinds,1] )
	tr$tip.label   <- translateTbl[tinds,2]
	tr$node.label  <- NULL

	tr		   <- nodeTimes(tr, youngestTip=youngestTip)

	return( tr )

	
}



# function to process a single BEAST tree line
processLatLonHostTree <- function( trLine, youngestTip=youngestTip, translateTbl=translateTbl ) {

		#trLine   <- lines[i]
		els	   <- strsplit(trLine, "\\[\\&R\\] ")[[1]]
		tempTr   <- els[2]

		# remove edge properties - attaching to previous child node
		tempTr   <- gsub("\\]\\:\\[\\&","--\\&",tempTr)

		# remove \\. from rate names
		tempTr   <- gsub("\\.rate","-rate",tempTr)

		# replace all , in [] with |
		b_regex	<- "\\[\\&"
		b_is		<- gregexpr(b_regex,tempTr)[[1]]
		c_regex	<- "[0-9]\\]"
		c_ie		<- gregexpr(c_regex,tempTr)[[1]]
		c_ie		<- c(c_ie, nchar(tempTr))
		old_txt	<- substring(tempTr, b_is, c_ie)
		new_txt	<- gsub("\\,","\\|", old_txt)
		for (k in 1:length(old_txt)) {
			tempTr <- gsub(old_txt[k],new_txt[k],tempTr,fixed=TRUE)
		}

		# read into tree object
		tr		<- read.tree(text=tempTr)

		# process tip and node names
		tipnames	<- apply(as.matrix(tr$tip.label), 1, getEl, sep="\\[")
		tipNos	<- tipnames[1,]
		tipProps	<- apply( as.matrix(tipnames[2,]), 1, extractNodeProperties )
				
		nodenames	<- apply(as.matrix(tr$node.label), 1, getEl, sep="\\[")
		nodeProps	<- apply( as.matrix(nodenames[2,]), 1, extractNodeProperties )
	
		all_props	<- c(tipProps, nodeProps)
		n		<- length(all_props)
		edgeLen	<- array(0,n)
		latlon	<- matrix(0, n, 2)
		host		<- array(0,n)
		for (k in 1:n) {
			edgeLen[k] <- getEdgeLen_val( all_props[[k]] )
			latlon[k,] <- as.matrix(getLatLon_val( all_props[[k]] ))
			host[k]    <- getHost_val( all_props[[k]] )
		}

		tr$state 	   <- host
		tr$latlon	   <- latlon

		toNodes	   <- tr$edge[,2]
		edgeLen_toNodes<- edgeLen[toNodes]
		tr$edge.length <- edgeLen_toNodes

		tr$tip.label   <- tipNos
		tinds		   <- match(tipNos, translateTbl[,1])
		all( tipNos==translateTbl[tinds,1] )
		tr$tip.label   <- translateTbl[tinds,2]
		tr$node.label  <- NULL

		tr		   <- nodeTimes(tr, youngestTip=youngestTip)

		return( tr )

}




# function to process a single BEAST tree line
processLatLonHostSubtypeTree <- function( trLine, youngestTip=youngestTip, translateTbl=translateTbl ) {

		#trLine   <- lines[i]
		els	   <- strsplit(trLine, "\\[\\&R\\] ")[[1]]
		tempTr   <- els[2]

		# remove edge properties - attaching to previous child node
		tempTr   <- gsub("\\]\\:\\[\\&","--\\&",tempTr)

		# remove \\. from rate names
		tempTr   <- gsub("\\.rate","-rate",tempTr)

		# replace all , in [] with |
		b_regex	<- "\\[\\&"
		b_is		<- gregexpr(b_regex,tempTr)[[1]]
		c_regex	<- "[0-9]\\]"
		c_ie		<- gregexpr(c_regex,tempTr)[[1]]
		c_ie		<- c(c_ie, nchar(tempTr))
		old_txt	<- substring(tempTr, b_is, c_ie)
		new_txt	<- gsub("\\,","\\|", old_txt)
		for (k in 1:length(old_txt)) {
			tempTr <- gsub(old_txt[k],new_txt[k],tempTr,fixed=TRUE)
		}

		# read into tree object
		tr		<- read.tree(text=tempTr)

		# process tip and node names
		tipnames	<- apply(as.matrix(tr$tip.label), 1, getEl, sep="\\[")
		tipNos	<- tipnames[1,]
		tipProps	<- apply( as.matrix(tipnames[2,]), 1, extractNodeProperties )
				
		nodenames	<- apply(as.matrix(tr$node.label), 1, getEl, sep="\\[")
		nodeProps	<- apply( as.matrix(nodenames[2,]), 1, extractNodeProperties )
	
		all_props	<- c(tipProps, nodeProps)
		n		<- length(all_props)
		edgeLen	<- array(0,n)
		latlon	<- matrix(0, n, 2)
		state		<- matrix(0, n, 2)
		for (k in 1:n) {
			edgeLen[k] <- getEdgeLen_val( all_props[[k]] )
			latlon[k,] <- as.matrix(getLatLon_val( all_props[[k]] ))
			state[k,]  <- all_props[[k]]$props[,2]
		}

		tr$state 	   <- state
		tr$latlon	   <- latlon

		toNodes	   <- tr$edge[,2]
		edgeLen_toNodes<- edgeLen[toNodes]
		tr$edge.length <- edgeLen_toNodes

		tr$tip.label   <- tipNos
		tinds		   <- match(tipNos, translateTbl[,1])
		all( tipNos==translateTbl[tinds,1] )
		tr$tip.label   <- translateTbl[tinds,2]
		tr$node.label  <- NULL

		tr		   <- nodeTimes(tr, youngestTip=youngestTip)

		return( tr )

}


	# interpolate lat lon on a great circle - using geoSphere
	interpol_ll <- function( fractT=fractT, llpair=llpair, n=n, tol = 0.0000001) {

		fromPt 		<- llpair[c(2,1)]
		toPt   		<- llpair[c(4,3)]
		if ( sum(abs(fromPt-toPt)) > tol ) {
			midLine		<- gcIntermediate( fromPt, toPt, n=n, addStartEnd=TRUE)
			nearestN		<- round(n*fractT)+1
			yy			<- midLine[nearestN,]
			return(yy)
		} else {
			return( toPt )
		}
	}




# END INTERNAL CUSTOM FUNCTIONS
##############################################################################################################

# DECPRECATED FUNCTIONS, DUE TO APE FORMAT CHANGE - THESE ARE UNLIKELY TO WORK NOW (2017)
useDECPRECATED <- FALSE
if ( useDECPRECATED ) {

# DEPRECATED
# function to process a BEAST trees file set (many lines)
# uses processLatLonHostTree
# returns trs object (a list of processed trees)
read_BEAST_trees_LatLonHost	<- function( treesName, nsamples=-1, burnin=0, step=1 ) {

	#includeNodeTimes 	<- TRUE
	#burnin		<- 0
	#step			<- 1
	#nsamples		<- -1

	lines			<- readLines(treesName)
	translateTbl 	<- getTranslation( lines )

	taxa		<- translateTbl[,2]
	decDates	<- as.numeric(apply(as.matrix(taxa), 1, getEl, final=TRUE, sep="\\|"))
	youngestTip <- max(decDates)
	
	trInds	<- grep("tree STATE",lines)

	# remove bad trees
	exInds	<- grep("NaN",lines)
	if (length(exInds) > 0) {
		print( paste("bad trees in",length(exInds),"lines out of",length(trInds),"trees") )
		trInds	<- setdiff(trInds,exInds)
	}
	ii		<- seq( (burnin+1), length(trInds), step )
	
	if ( (nsamples > 0) & (nsamples <= length(ii)) ) {
		ii	<- ii[1:nsamples]
	}
	trInds	<- trInds[ii]

	lines		<- lines[trInds]
	nsamples	<- length(lines)

	trs		<- vector("list", nsamples)
	trs		<- lapply(lines, processLatLonHostTree, youngestTip=youngestTip, translateTbl=translateTbl)

	return ( trs )

}





# DEPRECATED
# function to process a BEAST trees file set (many lines)
# uses processLatLonHostSubtypeTree
# returns trs object (a list of processed trees)
read_BEAST_trees_LatLonHostSubtype	<- function( treesName, nsamples=-1, burnin=0, step=1 ) {

	#includeNodeTimes 	<- TRUE
	#burnin		<- 0
	#step			<- 1
	#nsamples		<- -1

	lines			<- readLines(treesName)
	translateTbl 	<- getTranslation( lines )

	taxa		<- translateTbl[,2]
	decDates	<- as.numeric(apply(as.matrix(taxa), 1, getEl, final=TRUE, sep="\\|"))
	youngestTip <- max(decDates)
	
	trInds	<- grep("tree STATE",lines)

	# remove bad trees
	exInds	<- grep("NaN",lines)
	if (length(exInds) > 0) {
		print( paste("bad trees in",length(exInds),"lines out of",length(trInds),"trees") )
		trInds	<- setdiff(trInds,exInds)
	}
	ii		<- seq( (burnin+1), length(trInds), step )
	
	if ( (nsamples > 0) & (nsamples <= length(ii)) ) {
		ii	<- ii[1:nsamples]
	}
	trInds	<- trInds[ii]

	lines		<- lines[trInds]
	nsamples	<- length(lines)

	trs		<- vector("list", nsamples)
	trs		<- lapply(lines, processLatLonHostSubtypeTree, youngestTip=youngestTip, translateTbl=translateTbl)

	return ( trs )

}


}
# END DEPRECATED


##############################################################################################################
# MAIN FUNCTION TO USE IN OWN CODE
# this calls the internal custom functions above
# DO NOT MODIFY
##############################################################################################################


# function to process a BEAST trees file set (many lines)
# returns trs object (a list of processed trees)
# 28 oct 2015 - using processLatLonTree_oct2015
# 30 Aug 2016 - option to remove rate (rate is included if didnt do empirical trees)
read_BEAST_trees_LatLon	<- function( treesName, nsamples=-1, burnin=0, step=1, removeRate=FALSE ) {

	#includeNodeTimes 	<- TRUE
	#burnin		<- 0
	#step			<- 1
	#nsamples		<- -1

	lines			<- readLines(treesName)
	translateTbl 	<- getTranslation( lines )

	taxa		<- translateTbl[,2]
	decDates	<- as.numeric(apply(as.matrix(taxa), 1, getEl, final=TRUE, sep="\\|"))
	youngestTip <- max(decDates)
	
	trInds	<- grep("tree STATE",lines)

	# remove bad trees
	exInds	<- grep("NaN",lines)
	if (length(exInds) > 0) {
		print( paste("bad trees in",length(exInds),"lines out of",length(trInds),"trees") )
		trInds	<- setdiff(trInds,exInds)
	}
	ii		<- seq( (burnin+1), length(trInds), step )
	
	if ( (nsamples > 0) & (nsamples <= length(ii)) ) {
		ii	<- ii[1:nsamples]
	}
	trInds	<- trInds[ii]

	lines		<- lines[trInds]
	nsamples	<- length(trInds)

	trs		<- vector("list", nsamples)
	trs		<- lapply(lines, processLatLonTree_oct2015, youngestTip=youngestTip, translateTbl=translateTbl, removeRate=removeRate)

	return ( trs )

}

# function to process a BEAST trees file set (many lines); from empirical trees
# returns trs object (a list of processed trees)
# 29 oct 2015 - using processDiscreteTree_oct2015
read_BEAST_trees_Discrete	<- function( treesName, nsamples=-1, burnin=0, step=1, trait="Host", 
								youngestTip=-1 ) {

	#includeNodeTimes 	<- TRUE
	#burnin		<- 0
	#step			<- 1
	#nsamples		<- -1

	lines			<- readLines(treesName)
	translateTbl 	<- getTranslation( lines )

	taxa		<- translateTbl[,2]

	if (youngestTip==-1) {
		decDates	<- as.numeric(apply(as.matrix(taxa), 1, getEl, final=TRUE, sep="\\|"))
		youngestTip <- max(decDates)
	}
	
	trInds	<- grep("tree STATE",lines)

	# remove bad trees
	exInds	<- grep("NaN",lines)
	if (length(exInds) > 0) {
		print( paste("bad trees in",length(exInds),"lines out of",length(trInds),"trees") )
		trInds	<- setdiff(trInds,exInds)
	}
	ii		<- seq( (burnin+1), length(trInds), step )
	
	if ( (nsamples > 0) & (nsamples <= length(ii)) ) {
		ii	<- ii[1:nsamples]
	}
	trInds	<- trInds[ii]

	lines		<- lines[trInds]
	nsamples	<- length(trInds)

	trs		<- vector("list", nsamples)
	trs		<- lapply(lines, processDiscreteTree_oct2015, youngestTip=youngestTip, translateTbl=translateTbl, trait=trait)

	return ( trs )

}

#####################################################################################################
# 6 July 2017
# all in one function to do Host (discrete trait) and latlon (continous trait) trees

convert_BEAST_runs <- function(	discreteTreesName=discreteTreesName, trait="Host",
						latlonTreesName=latlonTreesName,
						burnin=1000, nsamples=-1, step=1) {
	
	print( paste("Formatting",discreteTreesName) )
	outname1	<- gsub(".trees.txt",paste(".",trait,".trs.Rdata",sep=""),discreteTreesName,fixed=TRUE)
	trs 		<- read_BEAST_trees_Discrete( discreteTreesName, trait=trait, burnin=burnin, step=step, nsamples=nsamples)	
	save(trs, file=outname1)
	
	print( paste("Formatting",latlonTreesName) )
	outname2	<- gsub(".trees.txt",".latlon.trs.Rdata",latlonTreesName,fixed=TRUE)
	trs 		<- read_BEAST_trees_LatLon( latlonTreesName, burnin=burnin, step=step, nsamples=nsamples)
	save(trs, file=outname2)

	return( list(r_discreteTreesName=outname1, r_latlonTreesName=outname2) )
	
}


