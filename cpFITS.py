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
	
	parser = argparse.ArgumentParser(description='Copies files from local folder (or an input list of files) to destination folder based on matching FITS headers.')
	parser.add_argument('destination', type=str, help="The destination folder.")
	parser.add_argument('--list', type=str, help="List of files to be checked.")
	parser.add_argument('-fh', '--FITSheader', type=str, help="The FITS header to look for.")
	parser.add_argument('-fv', '--FITSvalue', type=str, help="Value to match in the FITS header.")
	parser.add_argument('--ext', type=str, default='.*.(fits|fits.gz|fits.fz|fit)', help="FITS file descriptor: Pattern to match in filenames. Default is '.*.(fits|fits.gz|fits.fz|fit)'")
	parser.add_argument('--cone', action='store_true', help='Specify this option if you want a cone search.')
	parser.add_argument('--radius', type=float, help="The cone search radius in minutes of arc")
	parser.add_argument('--ra', type=float, help="The cone search RA value in degrees")
	parser.add_argument('--dec', type=float, help="The cone search DEC value in degrees")
	
	
	
	args = parser.parse_args()
	now = datetime.datetime.now()
	print "cpFITS.py running at:", now
	print 
	
	dirName = "."
	
	if args.list is not None:
		print "List specified:", args.list
		FITSFilenames = []
		# Load the list of files.
		filename = args.list
		fileList = open(filename, 'r')
		for line in fileList:
			FITSFilenames.append(str(line.strip()))
	else:
		
		folders = os.walk(dirName)
	
		fileList = []
		for root, dirs, files in folders:
			if root == dirName:
				fileList = files
			
	
		FITSFilenames = []
		search_re = re.compile(args.ext)
		for file in fileList:
			m = search_re.match(file)
			if (m): 
				FITSFilenames.append(file)
	
	if not os.path.exists(args.destination):
		print "Output folder, %s, does not exist... Creating it."%args.destination
		os.makedirs(args.destination)
	
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
	
	filesToCopy = []
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
				filesToCopy.append(f)
				print "Will copy to:", args.destination
			else:
				print "Won't copy."
	
		hdulist.close()
	
	for f in filesToCopy:
		copyCommand = ["cp -v " + f + " " + args.destination + "/."]
		subprocess.call(copyCommand, shell=True)
			
	
	sys.exit()
	
	
	
	try:
		objectFile = open(filename, 'rt')
	except:
		print "Could not find a list of candidate objects. Please create one and place it in '~/objects.theli'."
		
	objectList = []
	for line in objectFile:
		if line[0] == '#': continue
		parts = line.strip().split()
		target = {}
		if len(parts)<1: continue
		target['name'] = parts[1]
		if parts[0] == '*':
			target['date'] = lastnightString
		elif parts[0] == '.':
			target['date'] = todayString
		else:
			target['date'] = parts[0]
		objectList.append(target)
		
	for target in objectList:
		print "Date: %s, Target: %s"%(target['date'], target['name'])
	print
		
	for target in objectList:
		searchPath = "/obsdata/inta/" + target['date'] + "/r*.fit"
		print "Looking at FITS headers of files in %s for OBJECT: %s."%(searchPath, target['name'])
		
		getheadCommand = ["gethead " + searchPath + " OBJECT"]
		# getheadCommand.append(searchPath)
		# getheadCommand.append("OBJECT")
		print getheadCommand
		response = subprocess.check_output(getheadCommand, shell=True)
		lines = response.split('\n')
		wantedFiles = []
		for line in lines:
			parts = line.strip().split()
			if len(parts)<2: continue
			filename = parts[0]
			targetname = parts[1]
			if targetname == target['name']:
				wantedFiles.append(filename)
				
		if len(wantedFiles)==0:
			print "Did not find any files containing this object"
		else:
			print "The following files contain observations of this object: ", wantedFiles
			targetPath = args.outputdir + "/" + target['date'] + "/" + target['name']
			if not os.path.exists(targetPath):
				os.makedirs(targetPath)
			print "Copying these files to %s"%targetPath
			for w in wantedFiles:
				copyCommand = ["cp -vu" + " /obsdata/inta/" + target['date'] + "/" + w + " " + targetPath + "/."]
				subprocess.call(copyCommand, shell=True)
			
		print 
	
