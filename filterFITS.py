#!/usr/bin/env python 

import argparse
import sys, subprocess, os
import json, re, math
import datetime
from astropy.io import fits

def getUserHome():
	homeDir = os.path.expanduser('~')
	return str(homeDir)

def getUsername():
	username = os.getlogin()
	return str(username)

if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description='Filters a list of files to second list of files that only contain matching FITS criteria.')
	parser.add_argument('--path', type=str, help="Path to use to search for files.")
	parser.add_argument('--list', type=str, help="List of files to be checked.")
	parser.add_argument('-fh', '--FITSheader', type=str, help="The FITS header to look for.")
	parser.add_argument('-fv', '--FITSvalue', type=str, help="Value to match in the FITS header.")
	parser.add_argument('--cone', action='store_true', help='Specify this option if you want a cone search.')
	parser.add_argument('--radius', type=float, help="The cone search radius in degrees")
	parser.add_argument('--ra', type=float, help="The cone search RA value in degrees")
	parser.add_argument('--dec', type=float, help="The cone search DEC value in degrees")
	parser.add_argument('-o', '--output', type=str, help="Name for the output list.")
	
	args = parser.parse_args()
	now = datetime.datetime.now()
	print "filterFITS.py running at:", now
	print 
	
	if args.output is None:
		outputFilename = "output.list"
	else:
		outputFilename = args.output
	
	if args.list is not None:
		print "List specified:", args.list
		FITSFilenames = []
		# Load the list of files.
		filename = args.list
		fileList = open(filename, 'r')
		for line in fileList:
			FITSFilenames.append(str(line.strip()))
	else:
		print "Generating list from input path."
		# Get a list of files in the folder
		# First, check if the source data is there
		if not os.path.exists(args.path):
			print "The folder for the source data %s could not be found. Exiting."%dataPath
			sys.exit()
	
		searchString = "r.*.fits.fz"
		search_re = re.compile(searchString)
		 
		FITSFilenames = []
		for root, dirs, filenames in os.walk(args.path):
			for f in filenames:
				m = search_re.match(f)
				if (m): 
					FITSFilenames.append(os.path.join(root, f))

		print "%d files found."%len(FITSFilenames)
	
		
	if args.cone:
		coneSearch = True
		raDegrees = args.ra
		decDegrees = args.dec
		radiusDegrees = args.radius
		if (args.ra is None):
			print "Please specific an RA for the cone search."
			sys.exit()
		if (args.dec is None):
			print "Please specific a DEC for the cone search."
			sys.exit()
		if (args.radius is None):
			print "Please specific a radius for the cone search."
			sys.exit()
	else:
		coneSearch = False
		desiredParameter = args.FITSheader
		desiredValue = args.FITSvalue
		if desiredParameter is None: 
			print "Please specify a FITS header to match"
			sys.exit()
		if desiredValue is None: 
			print "Please specify a value for the FITS header %s to match"%desiredParameter
			sys.exit()

	raHeader = ["CRVAL1"]
	decHeader = ["CRVAL2"]
    
	if coneSearch:
		outputFileList = []
		
		for f in FITSFilenames:
			hdulist = fits.open(f)
			skip = False
			raValue = -1
			decValue = -1
			status = str(f)
			for card in hdulist:
				for r in raHeader:
					if r in card.header.keys():
						raValue = card.header[r]
						status+= " RA: %f"%raValue
				for d in decHeader:
					if d in card.header.keys():
						decValue = card.header[d]
						status+= " DEC: %f"%decValue
			if raValue == -1:
				status+= "  WARNING: Could not find an RA value in the FITS headers."
				skip = True
			if decValue == -1: 
				status+=  "  WARNING: Could not find a DEC value in the FITS headers."
				skip = True				
			hdulist.close()    

			if skip: 
				print status
				continue
			# Calculate if RA and DEC is within the cone
			raDistance = abs(raValue - raDegrees)
			decDistance = abs(decValue - decDegrees)
			totalDistance = math.sqrt(raDistance**2 + decDistance**2)
			status+= "  Distance: %f degrees"%totalDistance
			distanceMinutes = totalDistance * 60.
			if totalDistance <= radiusDegrees:
				outputFileList.append(f)
			
			print status
			
		print outputFileList
		
        
	                
	else:
		outputFileList = []
		for f in FITSFilenames:
			hdulist = fits.open(f)
			headerFound = False
			for card in hdulist:
				# Search for valuable metadata in the FITS headers
				try:
					value = card.header[desiredParameter]
					headerFound = True
				except KeyError:
					# print("WARNING: Could not find the FITS header you were looking for in this card: %s"%(desiredParameter))
					value = None
			if headerFound:
				print f, value, "...", 
				if value==desiredValue:
					outputFileList.append(f)
					print "Matches the criterion [%s = %s]"%(desiredParameter, desiredValue)
				else:
					print "No match."
	
			hdulist.close()
	
	outputfile = open(outputFilename, 'wt')
	for f in outputFileList:
		outputfile.write(f + "\n")
	outputfile.close()
			
	
	sys.exit()
	
	
	
