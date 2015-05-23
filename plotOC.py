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
import csv
	

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Plot an O-C.')
	parser.add_argument('filename', type=str, help='CSV file of eclipse times.')
	parser.add_argument('-e', '--ephemerisfilename', type=str, help='Ephemeris .dat file.')
	arg = parser.parse_args()
	print arg
	if arg.filename == None:
		print "Need an ephemeris file."
		sys.exit()

	ephemeris = timeClasses.ephemerisObject()
	
	ephemeris.loadFromFile(arg.ephemerisfilename)
	
	print ephemeris
	
	csvfile = open(arg.filename, 'r')
	reader = csv.reader(csvfile, delimiter=',')
	headings = reader.next()
	columns = []
	for h in headings:
		columnName = h.strip()
		columns.append(columnName)
	
	MJDs = []
	errors = []	
	for line in reader:
		values = [v.strip() for v in line]
		MJD = float(values[0])
		error = float(values[1])
		MJDs.append(MJD)
		errors.append(error)
	
	ocs = []
	cycles = []
	for MJD, error in zip(MJDs, errors):
		phase = ephemeris.getPhase(MJD)
		cycle = ephemeris.getOrbits(MJD)
		if phase>0.5: phase-=1
		phaseDifference = phase 
		print phaseDifference
		#if phaseDifference < 0: phaseDifference-=1
		ominusc = phaseDifference * ephemeris.Period
		print MJD, phase, phaseDifference, ominusc, ominusc*86400., cycle
		ocs.append(ominusc*86400.)
		cycles.append(cycle)
		
	matplotlib.pyplot.figure(figsize=(12, 8))
	matplotlib.pyplot.scatter(cycles, ocs)
	fig = matplotlib.pyplot.gcf()
	
	matplotlib.pyplot.show()
	fig.savefig('huaqr_oc.eps',dpi=100, format='eps')
	fig.savefig('huaqr_oc.png',dpi=100, format='png')
	sys.exit()
	
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
