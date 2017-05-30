#!/usr/bin/env python
import sys, os
import numpy, math
import argparse
import loadingSavingUtils, generalUtils
import spectrumClasses, timeClasses
import scipy.optimize
import copy
import ppgplot
import gSheets

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


class object:
	def __init__(self, id):
		self.id = id
		self.MJD = []
		self.mag = []
		self.err = []
		self.ephemeris = None
		
	def appendData(self, data):
		self.MJD.append(data['MJD'])
		self.mag.append(data['mag'])
		self.err.append(data['err'])
		
	def loadEphemeris(self):
		# Look in the local directory for a file called 'id'-ephem.dat and load it
		filename = self.id + "-ephem.dat"
		if os.path.exists(filename):
			self.ephemeris = timeClasses.ephemerisObject()
			self.ephemeris.loadFromFile(filename)
			return True
		
		return False
		
	def setHJD(self, HJD):
		self.HJD = HJD
		
	    

if __name__ == "__main__":  
	parser = argparse.ArgumentParser(description='Gets the observation dates for a set of targets and dumps them to disk.')
	parser.add_argument('objects', type=str, help='Object name list')
	arg = parser.parse_args()

	# Load the list objects we are going to plot
	objects = []
	objectFile = open(arg.objects, 'rt')
	for f in objectFile:
		o = object(f.strip())
		objects.append(o)
    
	
	# Load the fitted wavelength data from the Google Doc
	docInstance = gSheets.gSheetObject()
	docInstance.initCredentials()
	docInstance.setDocID('11fsbzSII1u1-O6qQUB8P0RzvJ8MzC5VHIASsZTYplXc')
	
	for o in objects:
		docInstance.setObjectName(o.id)
		docInstance.loadAllReadings()
    
		data = docInstance.readings
		print "Loaded data for %s from the Google Doc."%o.id
		print "%d data points loaded."%len(data)
		dates = []
		velocities = []
		velErrors = []
		good = []
		print "HJD\t\tVelocity (km/s)\tVel error"
		for index, d in enumerate(data):
			good = d['good']
			dates.append(d['HJD'])
			velocities.append(d['RV'])
			velErrors.append(d['RV error'])
    		print "%f\t%f\t%f"%(d['HJD'], d['RV'], d['RV error'])
    		
		o.HJD = dates
		o.velocities = velocities
		o.velErrors = velErrors
    	
	objectNames = [o.id for o in objects]
	print objectNames, "loaded."
	
	outputFile = open('observations.txt', 'wt')
	for o in objects:
		for date in o.HJD:
			outputFile.write("%s\t%f\n"%(o.id, date))
	outputFile.close()

