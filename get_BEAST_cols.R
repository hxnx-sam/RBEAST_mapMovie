# functions to get BEAST discrete trait colours
# S J Lycett
# 22 Dec 2015
# see H5NX work
# 11 March 2017 - added transparency

########################################################################################
# UTILITY FUNCTION
########################################################################################

get_BEAST_cols <- function( nstates, sat=0.7, bright=0.9, transparency=1 ) {
	if (transparency < 1) {
		ucols		<- hsv( (0:(nstates-1))/nstates, sat, bright, transparency)
	} else {
		ucols		<- hsv( (0:(nstates-1))/nstates, sat, bright)
	}
	return( ucols )
}