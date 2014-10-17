#!/usr/bin/env python
import numpy, math
import argparse
import astropy.io.fits
import csv
import astropy.coordinates.angles
import astropy.units as u
import subprocess
import glob

class wcssolution:
	def __init__(self, filename):
		self.filename = filename
		self.objects = []
		
	def setObjects(self, objects):
		self.objects = objects
		
	def getObjects(self):
		return self.objects
		
	def __str__(self):
		return self.filename + " containing " + str(len(self.objects)) + ' objects.'
		
	def findNearestObject(self, targetRADEC):
		closestDistance = 360
		matchedObject = {}
		for o in self.objects:
			catRA = astropy.coordinates.angles.Angle(o['field_ra'], 'degree')
			catDEC= astropy.coordinates.angles.Angle(o['field_dec'], 'degree')
			c1 = astropy.coordinates.ICRS(catRA.degree, catDEC.degree, unit=(u.degree, u.degree))
			c2 = astropy.coordinates.ICRS(targetRADEC[0].degree, targetRADEC[1].degree, unit=(u.degree, u.degree))
			separation = c1.separation(c2)
			totalDistance = separation.degree
			if totalDistance > 1:
				# More than one degree away... this ain't gonna be a match... don't bother checking the rest of the objects.
				return o, totalDistance
			if totalDistance<closestDistance:
				closestDistance = totalDistance
				matchedObject = o
		return matchedObject, closestDistance

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
	solution.setObjects(catalogObjects)
	
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
		print "Loading Astrometry.net solution file: ", f, 
		wcs_solution = loadWCSsolution(f)
		wcsSolutions.append(wcs_solution)
		# Also load the WCS FITS headers to get a coord somewhere in the field
		wcsFilename = f[:-4] + "wcs"
		refCoords = getRefCoords(wcsFilename)
		print "... reference coord is:", refCoords[0].to_string(unit=u.hour, sep=':'), refCoords[1].to_string(unit=u.degree, sep=':')
	print "Loaded a total of ", len(wcsSolutions), " wcs solutions."
	
	for i in range(len(objectList)):
	#for i in range(70,71):
		target = objectList[i]
		#print "Looking for a match for [%d]:"%i, target['name'], " { RA:", target['ra'], "DEC:", target['dec'], "}"
		matches = []
		for w in wcsSolutions:
			fieldmatch, distance = w.findNearestObject((target['ra'], target['dec']))
			if distance<(0.16): 
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
				print "    Object with pixel coordinates of (%3.2f, %3.2f)"%(matchedObject['field_x'], matchedObject['field_y']),
				print "is %3.3f arcseconds away"%(distance * 60 * 60)
				field_ra = astropy.coordinates.angles.Angle(matchedObject['field_ra'], 'degree')
				field_dec = astropy.coordinates.angles.Angle(matchedObject['field_dec'], 'degree')
				index_ra = astropy.coordinates.angles.Angle(matchedObject['index_ra'], 'degree')
				index_dec = astropy.coordinates.angles.Angle(matchedObject['index_dec'], 'degree')
				print "    Position given in the csv file :", target['ra'].to_string(unit=u.hour, sep=':'), target['dec'].to_string(unit=u.degree, sep=':')
				print "    Catalog position               :", index_ra.to_string(unit=u.hour, sep=':'), index_dec.to_string(unit=u.degree, sep=':')
				print "    Astrometry computed position   :", field_ra.to_string(unit=u.hour, sep=':'), field_dec.to_string(unit=u.degree, sep=':')
				c1 = astropy.coordinates.ICRS(index_ra.degree, index_dec.degree, unit=(u.degree, u.degree))
				c2 = astropy.coordinates.ICRS(field_ra.degree, field_dec.degree, unit=(u.degree, u.degree))
				separation = c1.separation(c2)
				print "    Separation between catalog and astrometry computed position: %3.3f"%separation.arcsecond, "arcseconds"
				print 
				
			
		
				
	
