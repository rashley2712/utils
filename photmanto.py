#!/usr/bin/env python
import sys
import cmd
import commandModule
import astropy.io.fits
import argparse

def loadFromFITSFile(filename):
	print "Loading:", filename
	
	inputFile = astropy.io.fits.open(filename)
	
	fileInfo = inputFile.info()
	
	colours = ['r', 'g', 'b']
	CCDs = { 'r': 'CCD 1', 'g': 'CCD 2', 'b': 'CCD 3'}
	
	c = colours[0]
	
	headers = inputFile[CCDs[c]].header
	data = inputFile[CCDs[c]].data
	columns = inputFile[CCDs[c]].columns
		
	inputFile.close()




if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='General purpose tools for loading, manipulating and plotting UCAM photometry data.')
	parser.add_argument('-s', '--script', type=str, help='The name of a script file containing commands to execute.')
	arg = parser.parse_args()
	print arg
	
	commands = commandModule.photCommands
	commands().cmdloop()

	sys.exit()
