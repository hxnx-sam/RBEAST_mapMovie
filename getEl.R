# useful function to split a sequence name
# S. J. Lycett
# 19 May 2011
# 6  Oct 2011

getEl	<- function( line, sep=",", ind=-1, final=FALSE, reconstruct=FALSE, ex=-1, fromEnd=FALSE ) {
	els	<- strsplit(line, sep)[[1]]

	if (ind[1] != -1) {
		if (fromEnd) {
			ind <- length(els)-(ind-1)
		}
	}

	if (final) {
		return( els[length(els)] )
	} else {

		if (reconstruct) {
			if (ex[1] > 0) {
				if (fromEnd) {
					ex <- length(els)-(ex-1)
				}
				ind <- setdiff((1:length(els)),ex)
			}

			newLine <- els[ind[1]]
			if (length(ind) > 1) {
				for (i in 2:length(ind)) {
					newLine <- paste(newLine, els[ind[i]], sep=sep)
				}
			}
			return ( newLine )
		} else {
			if ( ind[1] == -1 ) {
				return( els )
			} else {
				return( els[ind] )
			}
		}
	}
}