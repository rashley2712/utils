#!/usr/bin/env python
import numpy, math
import argparse
import astropy.io.fits
import astropy.wcs
import csv
import astropy.coordinates.angles
import astropy.units as u
import subprocess

import matplotlib.pyplot
import matplotlib.image

class wcssolution:
	def __init__(self, filename):
		self.filename = filename[:-5]
		self.catalogObjects = []
		self.solvedObjects = []
		
	def setCatalogObjects(self, objects):
		self.catalogObjects = objects
		
	def setSolvedObjects(self, objects):
		self.solvedObjects = objects
		
	def getCatalogObjects(self):
		return self.catalogObjects
		
	def __str__(self):
		return self.filename + " containing " + str(len(self.solvedObjects)) + ' total objects, with ' + str(len(self.catalogObjects)) + ' found in the catalog.'
		
	def findNearestObject(self, targetRADEC):
		lowestSeparation = 360
		matchedObject = {}
		
		for index, i in enumerate(self.solvedObjects):
			solvedRA, solvedDEC = i
			c1 = astropy.coordinates.ICRS(solvedRA.degree, solvedDEC.degree, unit=(u.degree, u.degree))
			c2 = astropy.coordinates.ICRS(targetRADEC[0].degree, targetRADEC[1].degree, unit=(u.degree, u.degree))
			separation = c1.separation(c2)
			if separation.degree>1:
				# Exit here... there is not going to be a match to this field
				return -1, separation.degree
			#print i, index, separation.arcsecond
			if separation.degree<(lowestSeparation):
				lowestSeparation = separation.degree
				matchedObject = {}
				matchedObject['fieldRA'] = solvedRA
				matchedObject['fieldDEC']=  solvedDEC
			print index, separation.arcsecond 
				
		print matchedObject
		return matchedObject, lowestSeparation

def readCSVCatalog(filename):
	objectList = []
	inputFile = open(arg.csvfile, 'r')
	csvReader = csv.reader(inputFile)
	for line in csvReader:
		if line[0][0] == '#':
			continue
		object = {}
		object['name'] = line[0]
		object['ra'] = astropy.coordinates.angles.Angle(line[1], 'hour')
		object['dec'] = astropy.coordinates.angles.Angle(line[2], 'degree')
		coords = [ line[1], line[2] ]
		objectList.append(object)
	
	return objectList

def getRefCoords(filename):	
	wcsFile = astropy.io.fits.open(filename)
	header = wcsFile[0].header
		
	equinox = float(header['EQUINOX'])
	referenceCoord = astropy.coordinates.angles.Angle(float(header['CRVAL1']), 'degree'), astropy.coordinates.angles.Angle(float(header['CRVAL2']), 'degree')
	referencePixel = float(header['CRPIX1']), float(header['CRPIX2'])
	
	return (referenceCoord)		
	
def loadOriginalFITSheaders(filename):
	inputFile = astropy.io.fits.open(filename)
	
def loadWCSsolution(filename):
	solution = wcssolution(filename)
	
	inputFile = astropy.io.fits.open(filename)
	
	headers = inputFile[1].header
	data = inputFile[1].data
	columns = inputFile[1].columns
	catalogObjects = []
	for item in data:
		catalogObject = {}
		for column in columns.names: 
			catalogObject[column] = item[columns.names.index(column)]
		catalogObjects.append(catalogObject)
	solution.setCatalogObjects(catalogObjects)
	
	# Also load the rdls file...
	rdlsFilename = filename[:-4] + "rdls"
	inputFile = astropy.io.fits.open(rdlsFilename)
	headers = inputFile[1].header
	data = inputFile[1].data
	columns = inputFile[1].columns
	solvedObjects = []
	for item in data:
		ra =  item[columns.names.index('RA')]
		dec = item[columns.names.index('DEC')]
		solvedRA = astropy.coordinates.angles.Angle(ra, 'degree')
		solvedDEC = astropy.coordinates.angles.Angle(dec, 'degree')
		solvedObjects.append((solvedRA, solvedDEC))
	solution.setSolvedObjects(solvedObjects)
	
	print solution
	
	return solution

def percentiles(data, lo, hi):
    """ Returns a normalised array where lo percent of the pixels are 0 and hi percent of the pixels are 255
    """
    max = data.max()
    dataArray = data.flatten()
    pHi = numpy.percentile(dataArray, hi)
    pLo = numpy.percentile(dataArray, lo)
    range = pHi - pLo
    scale = range/255
    data = numpy.clip(data, pLo, pHi)
    data-= pLo
    data/=scale
    return data


if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Loads a WCS solution and previews in matplotlib.')
	parser.add_argument('coreFilename', type=str, help='Filename (without extension) of the Astrometry.net solution files')
	
	arg = parser.parse_args()
	
	# Open the image and display it in MATPLOTLIB
	newFile = astropy.io.fits.open(arg.coreFilename + '.new')
	header = newFile[0].header
	
	imageData = newFile[0].data
	w = astropy.wcs.WCS(newFile[0].header)
	newFile.close()
		
	# Open the extracted sources
	axyFile = astropy.io.fits.open(arg.coreFilename + '.axy')
	header = axyFile[0].header
	axyData = axyFile[1].data
	axyFile.close()
	
	# Open the solved WCS solutions
	rdlsFile = astropy.io.fits.open(arg.coreFilename + '.rdls')
	header = rdlsFile[0].header
	rdlsData = rdlsFile[1].data
	rdlsFile.close()
	
	for r in rdlsData:
		pixelPosition = w.all_world2pix(r[0], r[1], 1)
		x, y = pixelPosition[0], pixelPosition[1]
		print r, w.all_world2pix(r[0], r[1], 1), x, y
		

	
	enhancedImage = percentiles(imageData, 50, 98)
	#enhancedImage = imageData
	matplotlib.pyplot.figure(figsize=(10, 10))
	matplotlib.pyplot.title("Original image")
	image = matplotlib.pyplot.imshow(enhancedImage, cmap='gray', interpolation='none')
	matplotlib.pyplot.show(block = False)
	
	matplotlib.pyplot.figure(figsize=(10, 10))
	matplotlib.pyplot.title("Extracted sources")
	image = matplotlib.pyplot.imshow(enhancedImage, cmap='gray', interpolation='none')
	for a in axyData:
		x = a[0]
		y = a[1]
		matplotlib.pyplot.gca().add_artist(matplotlib.pyplot.Circle((x,y), 25, color='red', fill=False, linewidth=1.0))
	matplotlib.pyplot.show(block = False)
	
	matplotlib.pyplot.figure(figsize=(10, 10))
	matplotlib.pyplot.title("Solved positions")
	image = matplotlib.pyplot.imshow(enhancedImage, cmap='gray', interpolation='none')
	
	for r in rdlsData:
		pixelPosition = w.all_world2pix(r[0], r[1], 1)
		x, y = pixelPosition[0], pixelPosition[1]
		matplotlib.pyplot.gca().add_artist(matplotlib.pyplot.Circle((x,y), 25, color='red', fill=False, linewidth=1.0))
	
		
	matplotlib.pyplot.show()

