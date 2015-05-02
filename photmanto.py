#!/usr/bin/env python
import sys
import commandModule
import astropy.io.fits
import argparse

def loadFromFITSFile(filename):
	
	inputFile = astropy.io.fits.open(filename)
	
	fileInfo = inputFile.info()
	
	print fileInfo
	
	CCD = 'CCD 1'
	fitsColumns = ["Counts_1"]
	
	headers = inputFile[CCD].header
	
	data = inputFile[CCD].data
	columns = inputFile[CCD].columns
	
	allData = []
	
	for index, item in enumerate(data):
		reading = {}
		reading['frame'] = index
		for col in columns.names:
			value = item[columns.names.index(col)]
			reading[col] = (value)
		allData.append(reading)
		sys.stdout.write("\rReading line %d"%index)
		sys.stdout.flush()
	
	inputFile.close()

	rows = len(allData)
	sys.stdout.write("\rRead %d lines with the following columns, %s\n"%(rows, str(columns.names)))
	sys.stdout.flush()
	
	
	



if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='General purpose tools for loading, manipulating and plotting ULTRACAM and ULTRASPEC photometry data.')
	parser.add_argument('script', type=str, nargs='?', help='The name of a script file containing commands to execute.')
	arg = parser.parse_args()
	print arg
	
	commands = commandModule.photCommands
	
	if arg.script != None:
		input = open(arg.script, 'rt')
		print "Running the commands found in :", arg.script
		try:
			commands(stdin=input).cmdloop()
		finally:
			input.close()
	
	commands.prompt = 'photmanto> '
	commands.use_rawinput = True
	commands().cmdloop()

	sys.exit()
