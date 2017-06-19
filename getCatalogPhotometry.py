#!/usr/bin/env python
import sys, os
import numpy
import argparse
import astropy.io.fits
import astroquery
from astropy import coordinates
from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
	
def getCatalogList(ra, dec):
	Vizier.ROW_LIMIT = 50
	c = coordinates.SkyCoord(ra,dec,unit=('deg','deg'),frame='icrs')
	results = Vizier.query_region(coordinates = c, radius= 5.0/60. * u.deg)
	print results
			
def getGalexResults(ra, dec):
	c = coordinates.SkyCoord(ra,dec,unit=('deg','deg'),frame='icrs')
	print "going to Astroquery for GALEX data"
	v = Vizier(columns=['_RAJ2000', '_DEJ2000','FUV', 'e_FUV', 'NUV', 'e_NUV', 'objid'])
	result = v.query_region(coordinates = c, radius="0d6m0s", catalog ='II/312/ais', verbose=True)
	print result[0]
	

def get2MASSResults(ra, dec, debug=False):
	radiusString = "0d6m0s"
	MASSkeys = ["J", "Jerr", "H", "Herr", "K", "Kerr"]
	Vizierkeys = ["Jmag", "e_Jmag", "Hmag", "e_Hmag", "Kmag", "e_Kmag"]
	c = coordinates.SkyCoord(ra,dec,unit=('deg','deg'),frame='icrs')
	if debug: print "going to Astroquery for 2MASS data"
	allColumns = ['_RAJ2000', '_DEJ2000', "_2MASS"] + Vizierkeys
	v = Vizier(columns=allColumns)
	result = v.query_region(coordinates = c, radius="0d6m0s", catalog ='II/246', verbose=debug)
	MASStable = result[0]
	print MASStable
	if len(MASStable) == 0:
		print "No results found within %s of %s"%(radiusString, c.to_string('hmsdms'))
		return None
	(catalog_RAs, catalog_DECs) = (MASStable['_RAJ2000'], MASStable['_DEJ2000'])
	catalog = SkyCoord(ra = catalog_RAs, dec = catalog_DECs)
	results = c.separation(catalog).degree * 3600
	MASStable['r'] = results
	MASStable.sort('r')
	if debug: print MASStable
	distance = MASStable['r'][0]
	result = {'distance': distance }
	for index, m in enumerate(MASSkeys):
		result[m] = MASStable[Vizierkeys[index]][0]	
	
	return result
			
	
class targetObject:
	
	def __init__(self):
		self.name = "unknown"
		self.ra = None
		self.dec = None
		self.photometry = []
		
	def __str__(self):
		retStr = "%s: %f %f"%(self.name, self.ra, self.dec) 
		return retStr
	
if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Gets catalog photometry for a list of targets.')
	parser.add_argument('filename', type=str, help='Target list ... tsv with "name\tra(deg)\tdec(deg)" columns.')	
	
	arg = parser.parse_args()
	objects = []
	inputFile = open(arg.filename, 'rt')
	for index, line in enumerate(inputFile):
		if line[0]=='#': continue
		fields = line.strip().split("\t")
		try:
			name = str(fields[0])
			ra= float(fields[1])
			dec= float(fields[2])
			target = targetObject()
			target.name = name
			target.dec = dec
			target.ra = ra
			objects.append(target)
		except ValueError as e:
			print "WARNING: Could not parse the data on line %d."%(index+1)
	
	print "%d targets loaded"%len(objects)

	for o in objects:
		print o
	
	objects = [objects[0]]
	for o in objects:
		o.MASS = get2MASSResults(o.ra, o.dec, debug=True)
		print o.name, o.MASS
	
	getCatalogList(o.ra, o.dec)
	sys.exit()
	
	
