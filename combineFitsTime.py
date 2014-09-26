#!/usr/bin/env python
import numpy, math
import matplotlib.pyplot
import matplotlib.image as mpimg
import argparse
import astropy.io.fits
import astropy.stats

def findMatchingTime(data, target):
	reading = (0, -1)
	
	for d in data:
		time = d[0]
		if time==target:
			reading = (time, d[1])
			continue
	
	return reading

def findClosestTime(data, target):
	distance = 1000.
	reading = (0, 0)
	for d in data:
		time = d[0]
		gap = time - target
		if abs(gap) < distance:
			distance = abs(gap)
			reading = (time, d[1])
	return reading
	
def filterData(dataArray):
	cleanData = []
	
	# This step removes the negative values
	for d in dataArray:
		items = len(d)
		if d[1]<0: continue
		if items>2:
			if d[2]<0: continue
		cleanData.append(d)
		
	return cleanData
	
def calcStats(data):
	np = numpy.array(data)
	values = np[:, 1]
	return (numpy.mean(values), numpy.std(values))
	
class observation:
	def __init__(self, MJD): 
		self.MJD = MJD
		

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Combines multiple photometry FITS files in the time domain.')
	parser.add_argument('datafiles', nargs='+', type=str, help='Input data file(s)')
	parser.add_argument('--outputFilename', '-o', type=str, default = "default_out.fits", help = "Output filename")
	
	arg = parser.parse_args()
	print arg
	
	datafile = arg.datafiles[0]
	
	inputFile = astropy.io.fits.open(datafile)
	
	CCDs = { 'r': 'CCD 1', 'g': 'CCD 2', 'b': 'CCD 3'}
	colours = ['r', 'g', 'b']
	
	c = colours[0]
	
	mainHeaders = inputFile[0].header
	print inputFile.info()	
	
	headers = inputFile[CCDs[c]].header
	data = inputFile[CCDs[c]].data
	columns = inputFile[CCDs[c]].columns
	
		
	for f in arg.datafiles[1:]:
		print "Opening:", f
		newFile = astropy.io.fits.open(f)
		print newFile.info()
		newHeaders = inputFile[CCDs[c]].header
		newData = inputFile[CCDs[c]].data
		newColumns = inputFile[CCDs[c]].columns
		for d in newData:
			print 
		newFile.close()
	
	print data
	print type(data)
		
	
	print "Writing to: ", arg.outputFilename	
	inputFile.writeto(arg.outputFilename, clobber = 'yes')
	
	for f in arg.datafiles[1:]:
		newFile = astropy.io.fits.open(f)
		newData = inputFile[CCDs['r']].data
		newFile.append(arg.outputFilename, newData)
