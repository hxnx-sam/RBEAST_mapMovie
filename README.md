# RBEAST_mapMovie
R code to plot BEAST trees output on a map and write a series of png images to make into a movie file

This repository contains custom R functions for reading and post-processing BEAST output trees.
Assuming that the BEAST trees were generated with discrete trait and continous trait spatial information, they can be used to create a movie of infection spread.

See the main file makeMapMovie.R for how to do this using H5N8 2014/2015 outbreak as an example

A version of these scripts were used to make the movie in the supplementary of:
The Global Consortium for H5N8 and Related Influenza Viruses: Lycett S, Bodewes R et al 
"Role for migratory wild birds in the global spread of avian influenza H5N8"
Science 14 Oct 2016: Vol 354, Issue 6309, pp 213-217
http://science.sciencemag.org/content/354/6309/213  
http://science.sciencemag.org/content/sci/suppl/2016/10/18/354.6309.213.DC1/aaf8852-MovieS1.mp4

