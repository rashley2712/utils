#!/usr/bin/env python

import argparse
import datetime
import datetime
import sys
import urllib2
import os

if __name__ == "__main__":
	

	parser = argparse.ArgumentParser(description='Looks for a run_log file from IDS and then uses it to seperate data files by central wavelength of observation.')
 
	arg = parser.parse_args()
	print arg
		
	files = os.listdir(".")
	
	print files
	
	for f in files:
		if ".int" in f:
			print "Found a file that appears to be an INT run log file:", f
			runlogfile = f
	
	filenames = []
	cenWavelengths = []		
	readLogfile = open(runlogfile, 'rt')
	for line in readLogfile:
		# print line.strip().split()
		parameters = line.strip().split()
		# print "Parameters:", len(parameters)
		cenLambda = line[112:118]
		try:
			cenWavelength = float(cenLambda)
		except ValueError:
			cenWavelength = -1
		if (len(parameters) > 1) and (cenWavelength!=-1): 
			filename = "r" + parameters[0] + ".fit"
			# print "Filename: %s CenLambda: %4.1f"%(filename, cenWavelength)
			filenames.append(filename)
			cenWavelengths.append(int(cenWavelength))
	
	uniqueWavelengths = []
	for w in cenWavelengths:
		if w not in uniqueWavelengths:
			uniqueWavelengths.append(w)
							
		
	readLogfile.close()
	
	print uniqueWavelengths
	
	currentDir = os.getcwd()
	print "Current Dir:", currentDir
	
	# Create the required sub-directories
	for d in uniqueWavelengths:
		newDirectory ="%dlambda"%d
		if not os.path.exists(newDirectory):
			print "Creating the directory: " + newDirectory
			os.mkdir(newDirectory)

	# Move the files
	for f, w in zip(filenames, cenWavelengths):
		g = "%dlambda/%s"%(w, f)
		print f, w, g
		fromFile = "./" + f
		toFile = "./" + g 
		print fromFile, toFile
		os.rename(fromFile, toFile)
	
	sys.exit()
	
	"""
		crosswordFilenames = []
		dates = []
		for f in files:
			if f.find('gdn.quick.')!=-1:
				crosswordFilenames.append(f)
				dateString = f[10:18]
				dates.append(dateString)
		print "Crosswords found in root folder..."
		print crosswordFilenames
		daysOld = []
		for d in dates:
			date = datetime.datetime.strptime(d, '%Y%m%d')
			days = (todaysDate - date).days
			daysOld.append(days)
			
		oldCrosswords = []	
		oldCrosswordDates = []
		for index, f in enumerate(crosswordFilenames):
			if daysOld[index] > 7:
				oldCrosswords.append(f)
				oldCrosswordDates.append(dates[index])
				
		print "Crosswords older than 7 days..."
		print oldCrosswords
		
		for index, filename in enumerate(oldCrosswords):
			date = datetime.datetime.strptime(oldCrosswordDates[index], '%Y%m%d')
			print filename,date
			month = date.month
			monthString = months[month-1]
			year = date.year
			print year, monthString
			directory = str(year) + "-" + monthString
			fullDirectory = dropboxPath + "/" + directory
			if not os.path.exists(fullDirectory):
				print "Creating the directory: " + fullDirectory
				os.mkdir(fullDirectory)
			oldFilename = dropboxPath + "/" + filename
			newFilename = dropboxPath + "/" + directory + "/" + filename
			print oldFilename, newFilename 
			os.rename(oldFilename, newFilename)
				
			
	print 'Completed successfully'
	
	"""
