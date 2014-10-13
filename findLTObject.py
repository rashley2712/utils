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
			#raDistance = catRA - targetRADEC[0]
			#decDistance = catDEC - targetRADEC[1]
			#totalDistance = math.sqrt( raDistance.degree * raDistance.degree + decDistance.degree*decDistance.degree)
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
		#print matchedObject, "at a distance of", closestDistance
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
		#c = SkyCoord(coords, FK4, unit=(u.deg, u.hourangle), obstime="J1992.21")
		objectList.append(object)
	
	return objectList
	

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
	
	arg = parser.parse_args()
	print arg
	
	objectList = readCSVCatalog(arg.csvfile)
	
	fileList = glob.glob("*.corr")
	
	wcsSolutions = []
	for f in fileList:
		print "Loading Astrometry.net solution file: ", f
		wcs_solution = loadWCSsolution(f)
		wcsSolutions.append(wcs_solution)	
	print "Loaded a total of:", len(wcsSolutions), " wcs solutions."
	
	for target in objectList:
	
		print "Looking for a match for:", target['name'], " { RA:", target['ra'], "DEC:", target['dec'], "}"
	
		
		
		closestDistance = 360
		for w in wcsSolutions:
			match, distance = w.findNearestObject((target['ra'], target['dec']))
			if (distance<closestDistance):
				closestDistance = distance
				matchedObject = match
				matchedField = w.filename[:-4] + "fits"
		
		if closestDistance > 0.16:
			print "    No match: (closest was %3.2f degrees away)"%(closestDistance)
		else:
			print "    Found a match in file: ", matchedField, " Object with pixel coordinates of (%3.2f, %3.2f)"%(matchedObject['field_x'], matchedObject['field_y']),
			print "is %3.2f arcseconds away"%(closestDistance * 60 * 60)
			
		print
		
				
	
