#!/usr/bin/env python 

import argparse
import sys, subprocess, os
import json, re
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
	parser.add_argument('list', type=str, help="List of files to be checked.")
	parser.add_argument('-fh', '--FITSheader', type=str, help="The FITS header to look for.")
	parser.add_argument('-fv', '--FITSvalue', type=str, help="Value to match in the FITS header.")
	parser.add_argument('--cone', action='store_true', help='Specify this option if you want a cone search.')
	parser.add_argument('--radius', type=float, help="The cone search radius in minutes of arc")
	parser.add_argument('--ra', type=float, help="The cone search RA value in degrees")
	parser.add_argument('--dec', type=float, help="The cone search DEC value in degrees")
	parser.add_argument('-o', '--output', type=str, help="Name for the output list.")
	
	args = parser.parse_args()
	now = datetime.datetime.now()
	print "filterFITS.py running at:", now
	print 
	
	if args.list is not None:
		print "List specified:", args.list
		FITSFilenames = []
		# Load the list of files.
		filename = args.list
		fileList = open(filename, 'r')
		for line in fileList:
			FITSFilenames.append(str(line.strip()))
	else:
		print "Please specify a list of files."
		sys.exit()
	
	# print "Found the following FITS files:", FITSFilenames
	
	
	
	if args.cone==True:
		raDegrees = args.ra
		decDegree = args.dec
		coneSearch = True
	else:
		desiredParameter = args.FITSheader
		desiredValue = args.FITSvalue
		if desiredParameter is None: 
			print "Please specify a FITS header to match"
			sys.exit()
		if desiredValue is None: 
			print "Please specify a value for the FITS header %s to match"%desiredParameter
			sys.exit()
			

	print "Searching for:", desiredParameter, " = ", desiredValue
	
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
				outputFileListfilesToCopy.append(f)
				print "Matches the criterion [%s = %s]"%(desiredParameter, desiredValue)
			else:
				print "No match."
	
		hdulist.close()
	
			
	
	sys.exit()
	
	
	
