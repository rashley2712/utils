#!/usr/bin/env python
import sys
try: 
	sys.path.remove('/home/astro/phsgan/Python64/lib/python/site-packages/astropy-0.3.2-py2.6-linux-x86_64.egg')
except (ValueError):
	print "No need to fix sys.path"
import numpy, math
import matplotlib.pyplot
import argparse
import astropy.io.fits
import astropy.stats
import loadingSavingUtils, statsUtils
import trm.dnl.molly
import trm.sla
import spectrumClasses, timeClasses
import scipy.optimize
import copy

def quad(x, a1, a2, a3):
	y = a1 * x * x + a2 *x + a3
	return y

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Loads data that has been downloaded from the Catalina archive. Plots them using MatPlotLib')
	parser.add_argument('inputFiles', type=str, nargs='+', help='Downloaded Catalina CSV file.')
	parser.add_argument('-e', type=str, help='Optional ephemeris file')
	 
	arg = parser.parse_args()
	print arg
	
	print "Astropy version:", astropy.__version__
	
	if arg.e!=None:
		# Load the ephemeris file
		hasEphemeris = True
		ephemeris = timeClasses.ephemerisObject()
		ephemeris.loadFromFile(arg.e)
		print ephemeris
	else:
		hasEphemeris = False

	
	for fileIndex, f in enumerate(arg.inputFiles):
	
		catalinaFile = f
	
		columnNames, data = loadingSavingUtils.loadNewCSV(f)

	print columnNames

	if hasEphemeris:
		BMJDs = []
		longitude = -110.73167
		latitude  = 32.41667
		elev      = 2510.
		print ephemeris.ra, ephemeris.dec
		for MJD in data['MJD']:
			results = trm.sla.utc2tdb(MJD, longitude, latitude, elev, ephemeris.ra/15., ephemeris.dec)
			print "%5.8f --> %5.8f"%(MJD, results[2])
			BMJDValue = results[2]
			BMJDs.append(BMJDValue)

		BJDs = [b + 2400000.5 for b in BMJDs]
		data['BJD'] = BJDs
		# print data['BJD']

		phases = [ephemeris.getPhase(h) for h in BJDs]
		offsetPhases = []
		for p in phases:
			if p<0.5:
				offsetPhase = p + 1.0
			else:
				offsetPhase = p
			offsetPhases.append(offsetPhase)

		# print phases

		data['phase'] = offsetPhases
	
	matplotlib.pyplot.figure(figsize=(14, 8))

	xValues = data['phase']
	yValues = data['Mag']
	yErrors = data['Magerr']
	matplotlib.pyplot.errorbar(xValues, yValues, color='b', yerr=yErrors, fmt = '.', ecolor='b', capsize=0)
	matplotlib.pyplot.gca().invert_yaxis()
	matplotlib.pyplot.show(block = False)
	fig = matplotlib.pyplot.gcf()
	fig.savefig('catalina_folded.eps',dpi=100, format='eps')
	fig.savefig('catalina_folded.png',dpi=100, format='png')


	left = 0.90
	right = 1.10
	trimmedX = []
	trimmedY = []
	trimmedError = []
	for x, y, err in zip(xValues, yValues, yErrors):
		if x<right and x>left:
			trimmedX.append(x)
			trimmedY.append(y)
			trimmedError.append(err)


	xValues = trimmedX
	yValues = trimmedY
	yErrors = trimmedError
	matplotlib.pyplot.figure(figsize=(14, 8))
	matplotlib.pyplot.errorbar(xValues, yValues, color='b', yerr=yErrors, fmt = '.', ecolor='b', capsize=0)
	matplotlib.pyplot.gca().invert_yaxis()
	fig = matplotlib.pyplot.gcf()
	matplotlib.pyplot.show()
	fig.savefig('catalina_folded_zoom.eps',dpi=100, format='eps')
	fig.savefig('catalina_folded_zoom.png',dpi=100, format='png')

