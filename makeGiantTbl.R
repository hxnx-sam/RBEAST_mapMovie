# make giant table - collect the host, subtype and latlon information per node from the trees
# assumes used empirical trees
# S J Lycett
# 4 July 2017 - extract from 2015-2017 code
# 6 July 2017

calcTreeLength <- function( tr ) {
	return( sum(tr$edge.length) )
}


# make giant table from converted Host + latlon trees (trs.Rdata) files
# actually can use any discrete trait in place of host
makeGiantTbl <- function(rootname=rootname, 
						 r_discreteTreesName = paste(rootname,".Host.trs.Rdata",sep=""),
						 r_latlonTreesName = paste(rootname,".latlon.trs.Rdata",sep=""),
						 nper=3, ntrees=-1) {
	# load host trs
	load(r_discreteTreesName)
	host_trs <- trs
	hlens	 <- lapply(host_trs, calcTreeLength)
	uhlens <- unique(hlens)
	print(paste("There are",length(uhlens),"unique discrete trs"))
	
	# load latlon trs
	load(r_latlonTreesName)
	latlon_trs <- trs
	llens	 <- lapply(latlon_trs, calcTreeLength)	
	ullens <- unique(llens)
	print(paste("There are",length(ullens),"unique latlon trs"))
	
	# all use the sample empirical tree set, this should be 1000 trees (or so)
	ulens	 <- unique(c(hlens,llens))
	print(paste("There are",length(ulens),"unique trs"))
	
	if (ntrees < 1) {
		ntrees	 <- length(ulens)
	}
	if (ntrees > length(ulens)) {
		ntrees <- length(ulens)
	}
	
	tol		<- 1e-8
	first		<- TRUE
	for (k in 1:ntrees) {

		skip <- FALSE

		h_kk <- which( abs(unlist(hlens)-ulens[[k]]) <= tol )
		if (length(h_kk) > 1) {
			h_kk <- sample(h_kk,nper,replace=length(h_kk) < nper)
		} else {
			if (length(h_kk)==1) {
				h_kk <- array(h_kk, nper)
			} else {
				skip <- TRUE
				print(paste("Warning tree",k,"not found in discrete trs"))
			}
		}
	
		l_kk <- which( abs(unlist(llens)-ulens[[k]]) <= tol )
		if (length(l_kk) > 1) {
			l_kk <- sample(l_kk,nper,replace=length(l_kk) < nper)
		} else {
			if (length(l_kk)==1) {
				l_kk <- array(l_kk, nper)
			} else {
				skip <- TRUE
				print(paste("Warning tree",k,"not found in latlon trs"))
			}
		}
	
		if (!skip) {
			k_tbl		<- c()
			for (p in 1:nper) {
			
				tr   			<- host_trs[[h_kk[p]]]
				fromHost 		<- tr$state[tr$edge[,1]]
				toHost   		<- tr$state[tr$edge[,2]]
				fromTime		<- tr$nodeTimes[tr$edge[,1]]
				toTime		<- tr$nodeTimes[tr$edge[,2]]

				tr			<- latlon_trs[[l_kk[p]]]
				fromLat		<- tr$latlon[tr$edge[,1],1]
				fromLon		<- tr$latlon[tr$edge[,1],2]
				toLat			<- tr$latlon[tr$edge[,2],1]
				toLon			<- tr$latlon[tr$edge[,2],2]

				treeNumber		<- array(k, length(tr$edge.length))
				repInd		<- array(p, length(tr$edge.length))
				tempTbl			<- cbind(treeNumber,repInd,
											fromTime,toTime,
											fromHost,toHost,
											fromLat,fromLon,toLat,toLon)

				if (p==1) {
					k_tbl <- tempTbl
				} else {
					k_tbl <- rbind(k_tbl, tempTbl)
				}
			
			}
			# end p
			
			if (first) {
				giant_tbl 	<- k_tbl
				first 	<- FALSE
			} else {
				giant_tbl 	<- rbind(giant_tbl, k_tbl)
			}

		} 
		# end !skip
	}	

	gtname <- paste(rootname,".giant_tbl.Rdata",sep="")	
	if (!first) {
		write.table(giant_tbl, file=paste(rootname,"_giant_tbl.txt",sep=""), sep="\t", 
				col.names=TRUE, row.names=FALSE, quote=FALSE)
		giant_tbl <- read.table( paste(rootname,"_giant_tbl.txt",sep=""), sep="\t",
					header=TRUE)
		save(giant_tbl, file=gtname )
	} else {
		print("Warning couldnt match the trees, the output is empty")
	}
	
	return( gtname )
	
}



# make giant table from converted Host, Subtype and latlon trees (trs.Rdata) files
# not as many checks in here as in the above
makeGiantTbl_Host_Subtype_LatLon <- function(rootname=rootname, 
						 r_hostTreesName = paste(rootname,".Host.trs.Rdata",sep=""),
						 r_subtypeTreesName = paste(rootname,".Subtype.trs.Rdata",sep=""),
						 r_latlonTreesName = paste(rootname,".latlon.trs.Rdata",sep=""),
						 nper=3, ntrees=-1) {
	# load host trs
	load(r_hostTreesName)
	host_trs <- trs
	hlens	 <- lapply(host_trs, calcTreeLength)
	
	# load subtype trs
	load(r_subtypeTreesName)
	subtype_trs <- trs
	slens	 <- lapply(subtype_trs, calcTreeLength)
	
	# load latlon trs
	load(r_latlonTreesName)
	latlon_trs <- trs
	llens	 <- lapply(latlon_trs, calcTreeLength)
	
	# all use the sample empirical tree set, this should be 1000 trees (or so)
	ulens	 <- unique(hlens)
	
	if (ntrees < 1) {
		ntrees	 <- length(ulens)
	}
	
	first		<- TRUE
	tol		<- 1e-8
	for (k in 1:ntrees) {
		h_kk <- which( abs(unlist(hlens)-ulens[[k]]) <= tol )
		h_kk <- sample(h_kk,nper,replace=length(h_kk) < nper)

		s_kk <- which( abs(unlist(slens)-ulens[[k]]) <= tol )
		s_kk <- sample(s_kk,nper,replace=length(s_kk) < nper)
	
		l_kk <- which( abs(unlist(llens)-ulens[[k]]) <= tol )
		l_kk <- sample(l_kk,nper,replace=length(l_kk) < nper)
	
		k_tbl		<- c()
			for (p in 1:nper) {
			
				tr   			<- host_trs[[h_kk[p]]]
				fromHost 		<- tr$state[tr$edge[,1]]
				toHost   		<- tr$state[tr$edge[,2]]
				fromTime		<- tr$nodeTimes[tr$edge[,1]]
				toTime			<- tr$nodeTimes[tr$edge[,2]]

				tr   			<- subtype_trs[[s_kk[p]]]
				fromSubtype 	<- tr$state[tr$edge[,1]]
				toSubtype   	<- tr$state[tr$edge[,2]]
				
				tr				<- latlon_trs[[l_kk[p]]]
				fromLat			<- tr$latlon[tr$edge[,1],1]
				fromLon			<- tr$latlon[tr$edge[,1],2]
				toLat			<- tr$latlon[tr$edge[,2],1]
				toLon			<- tr$latlon[tr$edge[,2],2]

				treeNumber		<- array(k, length(tr$edge.length))
				repInd		<- array(p, length(tr$edge.length))
				tempTbl			<- cbind(treeNumber,repInd,
											fromTime,toTime,
											fromHost,toHost,
											fromSubtype,toSubtype,
											fromLat,fromLon,toLat,toLon)

				if (p==1) {
					k_tbl <- tempTbl
				} else {
					k_tbl <- rbind(k_tbl, tempTbl)
				}
			
			}
			
		if (first) {
			giant_tbl <- k_tbl
			first <- FALSE
		} else {
			giant_tbl <- rbind(giant_tbl, k_tbl)
		}

	}
	
	gtname <- paste(rootname,".giant_tbl.Rdata",sep="")	
	if (!first) {
		write.table(giant_tbl, file=paste(rootname,"_giant_tbl.txt",sep=""), sep="\t", 
				col.names=TRUE, row.names=FALSE, quote=FALSE)
		giant_tbl <- read.table( paste(rootname,"_giant_tbl.txt",sep=""), sep="\t",
					header=TRUE)
		save(giant_tbl, file=gtname )
	} else {
		print("Warning couldnt match the trees, the output is empty")
	}
	
	return( gtname )
	
}