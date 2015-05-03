#!/usr/bin/env python
import sys
import commandModule
import astropy.io.fits
import argparse
import photometryClasses

slots = photometryClasses.slotCollection()

def loadFromFITSFile(filename):
	
	inputFile = astropy.io.fits.open(filename)
	
	fileInfo = inputFile.info()
	
	print fileInfo
	
	CCD = 'CCD 1'
	
	print inputFile[0].header
	headers = inputFile[CCD].header
	
	data = inputFile[CCD].data
	columns = inputFile[CCD].columns
	
	allData = []
	
	for index, item in enumerate(data):
		allData.append(item)
		sys.stdout.write("\rReading line %d"%index)
		sys.stdout.flush()
		# if index>100: break
	
	inputFile.close()
	rows = len(allData)
	sys.stdout.write("\rRead %d lines with the following columns, %s\n"%(rows, str(columns.names)))
	sys.stdout.flush()
	
	# Count the number of apertures in this data (using this method, the max is 9!)
	maxApertureIndex = 0
	for column in columns.names:
		try:
			apertureIndex = int(column[-1])
		except ValueError:
			apertureIndex = 0
		if apertureIndex > maxApertureIndex:
			maxApertureIndex = apertureIndex
	print "This data file has %d apertures."%(maxApertureIndex)
	
	MJDIndex = columns.names.index('MJD')
	MJD = [ data[MJDIndex] for data in allData]
	for aperture in range(1, maxApertureIndex+1):
		print "Processing aperture data:", aperture
		photometry = photometryClasses.photometryObject()
		photometry.time = MJD
		photometry.timeDescription = 'MJD'
		photometry.addValueDescription('MJD')
		exposureIndex = columns.names.index('Expose')
		photometry.addValueDescription('Expose')
		FWHMIndex = columns.names.index('FWHM')
		photometry.addValueDescription('FWHM')
		betaIndex = columns.names.index('beta')
		photometry.addValueDescription('beta')
		xIndex = columns.names.index('X_'+str(aperture))
		photometry.addValueDescription('X')
		yIndex = columns.names.index('Y_'+str(aperture))
		photometry.addValueDescription('Y')
		countsIndex = columns.names.index('Counts_' + str(aperture))
		photometry.addValueDescription('Counts')
		countsErrorIndex = columns.names.index('Sigma_' + str(aperture))
		photometry.addValueDescription('Sigma')
		skyCountsIndex = columns.names.index('Sky_' + str(aperture))
		photometry.addValueDescription('Sky')
		errorFlagIndex = columns.names.index('Eflag_' + str(aperture))
		photometry.addValueDescription('ErrorFlag')
		
		measurementArray = [(data[MJDIndex], data[exposureIndex], data[FWHMIndex], data[betaIndex], \
		                     data[xIndex], data[yIndex], data[countsIndex], data[countsErrorIndex],  \
							 data[skyCountsIndex], int(data[errorFlagIndex])) for data in allData]
		photometry.addData(measurementArray)
		slot = photometryClasses.slotObject()
		slot.channels = ['ULTRASPEC']
		slot.target = "Test"
		numSlots = slots.addSlot(slot)
		print photometry.valueDescriptions, measurementArray[0:3]
		print "Added the data to a new slot. Total number of slots is now: %d"%(numSlots)
		
	return

def listAllSlots(options):
	slotInfo = slots.getSlotInfo()
	print slotInfo
	return
	
		
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='General purpose tools for loading, manipulating and plotting ULTRACAM and ULTRASPEC photometry data.')
	parser.add_argument('script', type=str, nargs='?', help='The name of a script file containing commands to execute.')
	arg = parser.parse_args()
	# print arg
	
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
