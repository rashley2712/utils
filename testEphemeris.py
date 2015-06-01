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

	

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Test an ephemeris.')
	parser.add_argument('-e', '--filename', type=str, help='Ephemeris .dat file.')
	arg = parser.parse_args()
	print arg
	if arg.filename == None:
		print "Need an ephemeris file."
		sys.exit()

	ephemeris = timeClasses.ephemerisObject()
	
	ephemeris.loadFromFile(arg.filename)
	
	print ephemeris
	
	print "======================================================="
	print "Red"
	MJD = 57163.17148243
	
	print "Phase of %f is %f"%(MJD, ephemeris.getPhase(MJD)) 
	print "N Orbits of %f is %d"%(MJD, ephemeris.getOrbits(MJD)) 
	
	phase = ephemeris.getPhase(MJD)
	phaseDifference = 1 - phase
	
	ominusc = phaseDifference * ephemeris.Period
	
	print "O-C is %f days, or %f seconds."%(ominusc, ominusc*86400.)
	
	
	
	print "======================================================="
	print "Green"
	MJD = 57163.17148299
	
	print "Phase of %f is %f"%(MJD, ephemeris.getPhase(MJD)) 
	print "N Orbits of %f is %d"%(MJD, ephemeris.getOrbits(MJD)) 
	
	phase = ephemeris.getPhase(MJD)
	phaseDifference = 1 - phase
	
	ominusc = phaseDifference * ephemeris.Period
	
	print "O-C is %f days, or %f seconds."%(ominusc, ominusc*86400.)

	print "======================================================="
	print "Blue"
	MJD = 57163.17148816
	
	print "Phase of %f is %f"%(MJD, ephemeris.getPhase(MJD)) 
	print "N Orbits of %f is %d"%(MJD, ephemeris.getOrbits(MJD)) 
	
	phase = ephemeris.getPhase(MJD)
	phaseDifference = 1 - phase
	
	ominusc = phaseDifference * ephemeris.Period
	
	print "O-C is %f days, or %f seconds."%(ominusc, ominusc*86400.)


	print "======================================================="
	print "LT Rise"
	MJD = 56944.9054808
	
	print "Phase of %f is %f"%(MJD, ephemeris.getPhase(MJD)) 
	print "N Orbits of %f is %d"%(MJD, ephemeris.getOrbits(MJD)) 
	
	phase = ephemeris.getPhase(MJD)
	phaseDifference = 1 - phase
	
	ominusc = phaseDifference * ephemeris.Period
	
	print "O-C is %f days, or %f seconds."%(ominusc, ominusc*86400.)


	print "======================================================="
	print "trm BMJD"
	MJD = 57163.17069173
	
	print "Phase of %f is %f"%(MJD, ephemeris.getPhase(MJD)) 
	print "N Orbits of %f is %d"%(MJD, ephemeris.getOrbits(MJD)) 
	
	phase = ephemeris.getPhase(MJD)
	phaseDifference = 1 - phase
	
	ominusc = phaseDifference * ephemeris.Period
	
	print "O-C is %f days, or %f seconds."%(ominusc, ominusc*86400.)
