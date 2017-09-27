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
	
	parser = argparse.ArgumentParser(description='Lists FITS headers.')
	parser.add_argument('--list', type=str, help="List of files to be checked.")
	parser.add_argument('-fh', '--FITSheader', nargs='+', type=str, help="The FITS header to look for.")
	parser.add_argument('--ext', type=str, default='.*.(fits|fits.gz|fits.fz|fit)', help="FITS file descriptor: Pattern to match in filenames. Default is '.*.(fits|fits.gz|fits.fz|fit)'")
	parser.add_argument('--debug', action="store_true", help="Output debug information")
	
	args = parser.parse_args()
	if args.debug: 
		debug = True
	else: debug = False
	search_re = re.compile(args.ext)
	if debug: print args	

	now = datetime.datetime.now()
	if debug: print "lsFITS.py running at:", now, "\n"
	
	dirName = "."
	
	if args.list is not None:
		if debug: print "List specified:", args.list
		FITSFilenames = []
		# Load the list of files.
		filename = args.list
		fileList = open(filename, 'r')
		for line in fileList:
			m = search_re.match(line)
			if (m): 
				FITSFilenames.append(str(line.strip()))
	else:
		
		folders = os.walk(dirName)
	
		fileList = []
		for root, dirs, files in folders:
			if root == dirName:
				fileList = files
			
	
		FITSFilenames = []
		for file in fileList:
			m = search_re.match(file)
			if (m): 
				FITSFilenames.append(file)
	
	# print "Found the following FITS files:", FITSFilenames
	
	FITSFilenames = sorted(FITSFilenames)

	if args.FITSheader==None:
		FITSheaders= 'all'
	else:
		FITSheaders = args.FITSheader
			
	
	if debug: print "Searching for the headers :", FITSheaders
	
	headerList = []

	for f in FITSFilenames:
		if debug: print "In file:", f
		hdulist = fits.open(f)
		for card in hdulist:
			header = card.header
			if FITSheaders=='all':		
				for key in header.keys():
					if debug: print "%s\t%s"%(key, header[key])
					headerListObject = {}
					headerListObject['header'] = key
					headerListObject['value'] = header[key]
					headerListObject['filename'] = f
					headerList.append(headerListObject)

			else:
				for key in header.keys():
					for match in FITSheaders:
						if match==key:
							if debug: print "%s\t%s"%(key, header[key])
							headerListObject = {}
							headerListObject['header'] = key
							headerListObject['value'] = header[key]
							headerListObject['filename'] = f
							headerList.append(headerListObject)


	for results in headerList:
		print "%s\t%s\t%s\t"%(results['filename'], results['header'], results['value'])













			
