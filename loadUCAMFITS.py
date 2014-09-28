#!/usr/bin/env python
import numpy, math
import matplotlib.pyplot
import matplotlib.image as mpimg
import argparse
import astropy.io.fits
import astropy.stats


	
class observation:
	def __init__(self, MJD): 
		self.MJD = MJD
	
def getDates(data):
	returnArray = []
	for d in data:
		date = d['MJD']
		returnArray.append( date )
	
	return returnArray
		
	
def getColumn(data, columnName):
	returnArray = []
	for d in data:
		value = d[columnName]
		returnArray.append(value)
	
	return returnArray

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Loads a FITS file produced by ULTRACAM and plots the data.')
	parser.add_argument('datafiles', nargs='+', type=str, help='Input data file(s)')
	parser.add_argument('--outputFilename', '-o', type=str, default = "default_out.fits", help = "Output filename")
	
	arg = parser.parse_args()
	print arg
	
	datafile = arg.datafiles[0]
	
	inputFile = astropy.io.fits.open(datafile)
	
	CCDs = { 'r': 'CCD 1', 'g': 'CCD 2', 'b': 'CCD 3'}
	colours = ['r', 'g', 'b']
	colourDescriptions = {'r' : 'red', 'g': 'green', 'b': 'blue'}
	allPhotometry = { 'r': [], 'g': [], 'b': []}
	
	mainHeaders = inputFile[0].header
	print "File:", datafile	
	
	for c in colours:
		headers = inputFile[CCDs[c]].header
		data = inputFile[CCDs[c]].data
		columns = inputFile[CCDs[c]].columns

		print "Loading the fits card %s"%(CCDs[c]), "...   "
		photometry = []

		for d in data:
			dataObject = {}
			for col in columns: 
				dataObject[col.name] = d[columns.names.index(col.name)]

			photometry.append(dataObject)

		print "         ... Loaded %d datapoints for the colour %s"%(len(photometry), colourDescriptions[c])
		allPhotometry[c] = photometry
	
	inputFile.close()
	
	redPoints = allPhotometry['r']
	print len(redPoints)
	
	countsRed = getColumn(redPoints, 'Counts_1')
	comparisonRed = getColumn(redPoints, 'Counts_2')
	ratiosRed = numpy.array(countsRed) / numpy.array(comparisonRed)
	MJDs = getDates(redPoints)
	print countsRed
	
	matplotlib.pyplot.plot(MJDs, ratiosRed)
	
	matplotlib.pyplot.show()
