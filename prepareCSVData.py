#!/usr/bin/env python
import numpy, math
import matplotlib.pyplot
import matplotlib.image as mpimg
import argparse, sys
import astropy.io.fits
import astropy.stats
import loadingSavingUtils

	
def appendPhotometry(existing, new):
	print "Appending photometry"
	colours = ['r', 'g', 'b']
	for c in colours:
		photometry = existing[c]
		newPhotometry = new[c]
		for n in newPhotometry:
			photometry.append(n)
			
	for c in colours:
		numDataPoints = len(existing[c])
		print c, "now has", numDataPoints, "data points"

	return existing
		

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Reads CSV and FITS files and turns them into standard JSON files for loading into my various plotting programs,')
	parser.add_argument('datafile', nargs='+', type=str, help='Input data file(s)')
	parser.add_argument('-o', '--outfile', type=str, default='default', help='Name of the output file.')
	parser.add_argument('--channels', nargs='+', default = ['r', 'g', 'b'], type=str, help='Channels to plot')
	 
	arg = parser.parse_args()
	print arg
	colours = arg.channels
		
	# Find last '.' in the filename and grab the file extension
	filename = str(arg.datafile[0])
	index = filename.rfind('.')
	extension = filename[index+1:].lower()
	filenameNoExtension = filename[:index]
	filetype = "unknown"
	if extension=="fits":
		filetype = 'fits'
	elif extension=="csv":
		filetype = "csv"
	
	if filetype == 'unknown':
		print "Sorry that filetype is unknown. Try CSV or FITS."
		
	
	if filetype == 'fits':
		print "Loading", filename
		photometry = loadingSavingUtils.loadFITSFile(filename, colours = colours)
		
		if (len(arg.datafile)>1):
			additionalFiles = arg.datafile[1:]
			for newFilename in additionalFiles:
				print "...also loading", newFilename
				additionalPhotometry = loadingSavingUtils.loadFITSFile(newFilename)
				
				photometry = appendPhotometry(photometry, additionalPhotometry)
	
	
	for c in colours:
		photometry[c], numRemoved = loadingSavingUtils.removeNegativeValues(photometry[c])
		if numRemoved>0: print "..removed %d negative values from %s."%(numRemoved, c)
	
	for c in colours:
		photometry[c], numRemoved = loadingSavingUtils.removeZeroValues(photometry[c])
		if numRemoved>0: print "..removed %d zero values from %s."%(numRemoved, c)
	
	
	if arg.outfile == 'default':
		outputFilename = filenameNoExtension + ".csv"
	else:
		outputFilename = arg.outfile
	loadingSavingUtils.writeCSV(outputFilename, photometry, colours = colours)
	
	
	
	sys.exit()
	
	
