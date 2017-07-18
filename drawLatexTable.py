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
import matplotlib

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
	parser = argparse.ArgumentParser(description='Loads RV values from the Google sheet and writes it as a LaTex table.')
    # parser.add_argument('inputfile', type=str, help='Filename of the CSV file containing the RVs.')
	parser.add_argument('objects', type=str, help='Object name list')
	parser.add_argument('-o', '--output', type=str, default = 'table.tex', help='Output filename. (Default: table.tex)')
	arg = parser.parse_args()
	outputFilename = arg.output	

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
			if good==0: 
				print bcolors.WARNING + "%f\t%f\t%f"%(d['HJD'], d['RV'], d['RV error']) + bcolors.ENDC 
			else: 
				print "%f\t%f\t%f"%(d['HJD'], d['RV'], d['RV error'])
				dates.append(d['HJD'])
				velocities.append(d['RV'])
				velErrors.append(d['RV error'])
    
		o.HJD = dates
		o.velocities = velocities
		o.velErrors = velErrors
    	
	
	print "Writing latex file to:", outputFilename
	outputFile = open(outputFilename, 'wt')
	heading = """\\begin{table}
\\begin{tabular}{ l  l  l  l }
	\\hline
	WD & HJD & Radial velocity & Error (km/s) \\\\
	   &     & (km/s)          & (km/s) \\\\
"""
	outputFile.write(heading)
	for o in objects:
		outputFile.write("\t\\hline\n")
		for hjd, v, verr in zip(o.HJD, o.velocities, o.velErrors):
			print o.id, hjd, v, verr
			outputFile.write("\t%s & %f & %4.1f & %2.1f \\\\\n"%(o.id, hjd, v, verr))
	footer = """  \hline
  \end{tabular}
  \label{tab:rvdata}
\end{table}
"""
	outputFile.write(footer)	
	outputFile.close()
	
	print "Writing raw ascii file too"
	outputFile = open("dump.txt", 'wt')
	for o in objects:
		for hjd, v, verr in zip(o.HJD, o.velocities, o.velErrors):
			outputFile.write("%s\t%f\t%4.1f\t%2.1f\n"%(o.id, hjd, v, verr))
	outputFile.close()

