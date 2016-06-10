#!/usr/bin/env python
import numpy, math
import argparse, sys
import astropy.io.fits
import astropy.stats
import loadingSavingUtils, statsUtils
import timeClasses
import generalUtils
import slbarycentric

if __name__ == "__main__":
	telescopes = []
	# MacDonald
	telescope = {}
	telescope['name'] = "MacDonald Observatory, Texas"
	telescope['longitude'] = 104.0225
	telescope['latitude'] = 30.67139
	telescope['altitude'] = 2070.
	telescopes.append(telescope)
	
	# VLT
	telescope = {}
	telescope['name'] = "Very Large Telescope, UT3 (Melipal)"
	telescope['longitude'] = 289.5972
	telescope['latitude'] = -24.6253
	telescope['altitude'] = 2635.
	telescopes.append(telescope)
	
	# TNT observatory
	telescope = {}
	telescope['name'] = "Thai National Observatory"
	telescope['longitude'] = 98.48
	telescope['latitude'] = 18.57
	telescope['altitude'] = 2457.
	telescopes.append(telescope)
	
	# WHT observatory
	telescope = {}
	telescope['name'] = "William Herschel Telescope"
	telescope['longitude'] = 342.1184
	telescope['latitude'] = 28.7606
	telescope['altitude'] = 2326.
	telescopes.append(telescope)
	
	# SAFT observatory
	telescope = {}
	telescope['name'] = "Warwick One Metre"
	telescope['longitude'] = 342.1184
	telescope['latitude'] = 28.7606
	telescope['altitude'] = 2326.
	telescopes.append(telescope)
	
	# NTT observatory
	telescope = {}
	telescope['name'] = "New Technology Telescope"
	telescope['longitude'] = 289.27
	telescope['latitude'] = 29.2567
	telescope['altitude'] = 2347.
	telescopes.append(telescope)

	parser = argparse.ArgumentParser(description='Computes the BMJD from the MJD.')
	parser.add_argument('inputfiles', type=str, nargs='+', help='Input data in CSV format')
	parser.add_argument('-e', '--ephemeris', type=str, help='Ephemeris file containing the RA and DEC of the target.')
	parser.add_argument('--list', action='store_true', help='Specify this option if the input file is actually a list of input files.')
	
	arg = parser.parse_args()
	print arg
	
	if arg.ephemeris!=None:
		# Load the ephemeris file
		hasEphemeris = True
		ephemeris = timeClasses.ephemerisObject()
		ephemeris.loadFromFile(arg.ephemeris)
		print ephemeris
	else:
		hasEphemeris = False
		print "We need to know the coordinates in order to compute the BMJD. Please place them in the ephemeris file."
		sys.exit()
		
	filenames = []
	if arg.list:
		# Load the list of files.
		if len(arg.inputfiles)>1: 
			print "You can only give me one list of filenames."
			sys.exit()
		filename = arg.inputfiles[0]
		fileList = open(filename, 'r')
		for line in fileList:
			filenames.append(str(line.strip()))
	else:
		filenames = arg.inputfiles
	
	
	allData = []
	
	for filename in filenames:
		columnNames, photometry = loadingSavingUtils.loadNewCSV(filename)
		photometry = generalUtils.filterOutNaNs(photometry)
		photometry['runName'] = filename
		allData.append(photometry)

	xColumn = columnNames[5]
	yColumn = columnNames[1]
	yErrors = columnNames[2]
	
	""" Data is now loaded 
	"""
	# Sort the data by mjd
	sortedData = sorted(allData, key=lambda object: object[xColumn][0], reverse = False)
	allData = sortedData
	
	print allData[0]
	times = []
	for d in allData[0][xColumn]:
		times.append(d)
	
	print "RA, DEC of the target:", ephemeris.ra, ephemeris.dec	
	
	print "Calculating barycentric MJD or BMJD"	
	
	telescope = telescopes[0]
	obsLong = telescope['longitude']
	obsLat = telescope['latitude']
	obsAlt = telescope['altitude']
	obsName = telescope['name'] 
	
	print "Using the following telescope information:\nName: %s\nLongitude: %.2f [deg]\tLatitude: %.2f [deg]\t Altitude: %.1f [m]"%(obsName, obsLong, obsLat, obsAlt)

	obsLocation = astropy.coordinates.EarthLocation(lon = obsLong, lat = obsLat, height=obsAlt)

	ra, dec = ephemeris.ra, ephemeris.dec
	targetRADEC = generalUtils.toSexagesimal((ra,dec))
	print "Target position: %s (%f, %f)"%(targetRADEC, ra, dec)

	targetCoords = astropy.coordinates.SkyCoord(ra, dec, unit='deg')
	BMJD = []
	t = times[0]
	observationTime = slbarycentric.Time(t, format='mjd', location = obsLocation)
	for index, t in enumerate(times):
		observationTime.__init__(t, format='mjd', location = obsLocation)
		delta, bcor = observationTime.bcor(targetCoords)
		bmjd = float(bcor.mjd)
		BMJD.append(bmjd)
		sys.stdout.write("\r[%d/%d]  MJD %5.8f ---> BMJD %5.8f  = %f seconds   "%(index, len(times)-1, t, bmjd, delta))
		sys.stdout.flush()
	print
	
	for t,b in zip(times, BMJD):
		print t, b
	
	output = "newdata.csv"	
	outputFile = open(output, 'wt')
	outputFile.write("BMJD, MJD, Mag, MagError\n")
	for bmjd, mjd, y, yerr in zip(BMJD, times, allData[0][yColumn], allData[0][yErrors]):
		print bmjd, mjd, y, yerr
		outputFile.write("%f, %f, %f, %f\n"%(bmjd, mjd, y, yerr))
	outputFile.close()
	
		
