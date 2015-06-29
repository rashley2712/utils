#!/usr/bin/env python
import numpy, math
import matplotlib.pyplot
import matplotlib.image as mpimg
import argparse, sys
import astropy.io.fits
import astropy.stats
import loadingSavingUtils, statsUtils
import scipy.optimize
import time
import random
import timeClasses
import trm.sla

	

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Test an ephemeris.')
	parser.add_argument('-e', '--filename', type=str, help='Ephemeris .dat file.')
	parser.add_argument('csvfile', type=str, help='CSV file containing MJDs.')
	arg = parser.parse_args()
	print arg
	if arg.filename == None:
		print "Need an ephemeris file."
		sys.exit()

	ephemeris = timeClasses.ephemerisObject()
	
	ephemeris.loadFromFile(arg.filename)
	
	print ephemeris

	# Now load some light-curve data
	columnNames, photometry = loadingSavingUtils.loadNewCSV(arg.csvfile)
	xColumn = columnNames[0]
	
	obsLong = 98.48
	obsLat = 18.57
	obsAlt = 2457.2
	ra = ephemeris.ra /15.
	dec = ephemeris.dec
	print "Input ra", ra, "dec", dec, 
	for MJD in photometry[xColumn]:
		result = trm.sla.utc2tdb(MJD, obsLong, obsLat, obsAlt, ra, dec)
		bmjd = result[2]
		print "%6.8f --> %6.8f"%(MJD, bmjd)
	
