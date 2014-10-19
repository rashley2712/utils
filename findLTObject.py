#!/usr/bin/env python
import numpy, math
import argparse
import astropy.io.fits
import astropy.wcs
import csv
import astropy.coordinates.angles
import astropy.units as u
import subprocess
import glob

class wcssolution:
	def __init__(self, filename):
		self.filename = filename[:-5]
		self.catalogObjects = []
		self.solvedObjects = []
		self.extractedObjects = []
		self.wcs = None
		
	def setWCS(self, w):
		self.wcs = w
		
	def setCatalogObjects(self, objects):
		self.catalogObjects = objects
		
	def setSolvedObjects(self, objects):
		self.solvedObjects = objects
		
	def setExtractedObjects(self, objects):
		self.extractedObjects = objects
		
	def __str__(self):
		return self.filename + " containing " + str(len(self.extractedObjects)) + ' total objects, with ' + str(len(self.catalogObjects)) + ' found in the catalog.'
		
	def findNearestObject(self, targetRADEC):
		lowestSeparation = 360
		matchedObject = {}
		
		for index, i in enumerate(self.extractedObjects):
			w = self.wcs
			(x, y) = i['x'], i['y']
			world = w.all_pix2world(x, y, 1)
			ra = world[0]
			dec = world[1]
			c1 = astropy.coordinates.ICRS(ra, dec, unit=(u.degree, u.degree))
			c2 = astropy.coordinates.ICRS(targetRADEC[0].degree, targetRADEC[1].degree, unit=(u.degree, u.degree))
			separation = c1.separation(c2)
			if separation.degree>1:
				# Exit here... there is not going to be a match to this field
				return -1, separation.degree
			#print i, index, separation.arcsecond
			if separation.degree<(lowestSeparation):
				lowestSeparation = separation.degree
				matchedObject = {}
				matchedObject['fieldRA'] = ra
				matchedObject['fieldDEC']=  dec
			#print index, separation.arcsecond 
				
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

	# Also load the .wcs file in order to read the full WCS solution parameters
	wcsFilename = filename[:-4] + "wcs"
	wcsFile = astropy.io.fits.open(wcsFilename)
	header = wcsFile[0].header
	w = astropy.wcs.WCS(wcsFile[0].header)
	wcsFile.close()
	solution.setWCS(w)
	
	# Also load the .axy file in order to ensure we have *all* of the objects (not just the matched ones)
	axyFilename = filename[:-4] + "axy"
	axyFile = astropy.io.fits.open(axyFilename)
	headers = axyFile[1].header
	data = axyFile[1].data
	columns = axyFile[1].columns
	extractedObjects = []
	for item in data:
		extractedObject = {}
		extractedObject['x'] =  item[columns.names.index('X')]
		extractedObject['y'] =  item[columns.names.index('Y')]
		extractedObject['flux'] =  item[columns.names.index('FLUX')]
		extractedObject['background'] =  item[columns.names.index('BACKGROUND')]
		extractedObjects.append(extractedObject)
	solution.setExtractedObjects(extractedObjects)

	print solution
	
	return solution

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Find LT object. Looks in a CSV catalog for an object and matches it to an object in an LT image.')
	parser.add_argument('csvfile', type=str, help='Input catalog in CSV format')
	parser.add_argument('-n', '--nomatch', action='store_true', help='Suppress reporting of non-matches')
	
	arg = parser.parse_args()
	
	objectList = readCSVCatalog(arg.csvfile)
	
	fileList = glob.glob("*.corr")
	
	wcsSolutions = []
	for f in fileList:
		print "Loading Astrometry.net solution file: ",  
		wcs_solution = loadWCSsolution(f)
		wcsSolutions.append(wcs_solution)
		# Also load the WCS FITS headers to get a coord somewhere in the field
		wcsFilename = f[:-4] + "wcs"
		refCoords = getRefCoords(wcsFilename)
		print "... reference coord is:", refCoords[0].to_string(unit=u.hour, sep=':'), refCoords[1].to_string(unit=u.degree, sep=':')
	print "Loaded a total of ", len(wcsSolutions), " wcs solutions."
	
	for i in range(len(objectList)):
	#for i in range(18,19):
		target = objectList[i]
		#print "Looking for a match for [%d]:"%i, target['name'], " { RA:", target['ra'], "DEC:", target['dec'], "}"
		matches = []
		for w in wcsSolutions:
			fieldmatch, distance = w.findNearestObject((target['ra'], target['dec']))
			if (fieldmatch==-1):
				if (arg.nomatch==False):
					print "No match to this field [%s]"%w.filename[:-4]
			elif distance<(0.16): 
				match = {}
				#print "Found in field", distance, w.filename
				match['object'] = fieldmatch
				match['field'] = w
				match['distance'] = distance
				matches.append(match)
		
		if len(matches)>0:
			print "For:", target['name'], "[%d]"%i
			
			for m in matches:						
				matchFilename = m['field'].filename
				matchedObject = m['object']
				distance = m['distance']
				
				print "  Found a match in file: ", matchFilename
				
				field_ra = astropy.coordinates.angles.Angle(matchedObject['fieldRA'], 'degree')
				field_dec = astropy.coordinates.angles.Angle(matchedObject['fieldDEC'], 'degree')
				print "    Position given in the csv file :", target['ra'].to_string(unit=u.hour, sep=':'), target['dec'].to_string(unit=u.degree, sep=':')
				print "    Astrometry computed position   :", field_ra.to_string(unit=u.hour, sep=':'), field_dec.to_string(unit=u.degree, sep=':')
				c1 = astropy.coordinates.ICRS(field_ra.degree, field_dec.degree, unit=(u.degree, u.degree))
				c2 = astropy.coordinates.ICRS(target['ra'].degree, target['dec'].degree, unit=(u.degree, u.degree))
				separation = c1.separation(c2)
				print "    Separation between input csv file and astrometry computed position: %3.3f"%separation.arcsecond, "arcseconds"
				
			print
				#if matchedObject['catalogMatch'] == 1:
				#	print "  Found in catalog"
				#	catRA  = matchedObject['catRA']
				#	catDEC = matchedObject['catDEC']
				#	print "    Catalog position           :", catRA.to_string(unit=u.hour, sep=':'), catDEC.to_string(unit=u.degree, sep=':')
				#else:
				#	print "    Not found in Astrometry.net's source catalog"
				#	print "    Object with pixel coordinates of (%3.2f, %3.2f)"%(matchedObject['pixel_x'], matchedObject['pixel_y']),
				#print "is %3.3f arcseconds away"%(distance * 60 * 60)
				#c1 = astropy.coordinates.ICRS(index_ra.degree, index_dec.degree, unit=(u.degree, u.degree))
				#c2 = astropy.coordinates.ICRS(field_ra.degree, field_dec.degree, unit=(u.degree, u.degree))
				#separation = c1.separation(c2)
				#print "    Separation between catalog and astrometry computed position: %3.3f"%separation.arcsecond, "arcseconds"
				#print 
				
			
		
				
	
