#!/usr/bin/env python
import sys
import commandModule
import astropy.io.fits
import argparse
import photometryClasses
import generalUtils
import photmantoPlotting

slots = photometryClasses.slotCollection()
state = {'plotter':'pgplot'}
plotterHandle = 0

def loadFromFITSFile(filename, maxRows=0):
	""" Loads photometry data from a FITS file generated by ulog2fits. If maxRows!=0 then it loads up to maxRows."""
	inputFile = astropy.io.fits.open(filename)
	#fileInfo = inputFile.info()
	
	CCD = 'CCD 1'
	
	headerBlock = str(inputFile[0].header)
	
	# Get some header info
	targetName = generalUtils.getKeyValueFromFITSHeader('target', headerBlock)
	filterName = generalUtils.getKeyValueFromFITSHeader('filter', headerBlock)
	runName = generalUtils.getKeyValueFromFITSHeader('Data file name', headerBlock, terminator=' ')
	headers = inputFile[CCD].header
	
	data = inputFile[CCD].data
	columns = inputFile[CCD].columns
	
	allData = []
	
	for index, item in enumerate(data):
		allData.append(item)
		# sys.stdout.write("\rReading line %d"%index)
		# sys.stdout.flush()
		if maxRows!=0 and index>=maxRows-1: break
	
	inputFile.close()
	rows = len(allData)
	# sys.stdout.write("\rRead %d lines with the following columns, %s\n"%(rows, str(columns.names)))
	# sys.stdout.flush()
	
	# Count the number of apertures in this data (using this method, the max is 9!)
	maxApertureIndex = 0
	for column in columns.names:
		try:
			apertureIndex = int(column[-1])
		except ValueError:
			apertureIndex = 0
		if apertureIndex > maxApertureIndex:
			maxApertureIndex = apertureIndex
	# print "This data file has %d apertures."%(maxApertureIndex)
	
	MJDIndex = columns.names.index('MJD')
	MJD = [ data[MJDIndex] for data in allData]
	for aperture in range(1, maxApertureIndex+1):
		print "Processing aperture data:", aperture
		photometry = photometryClasses.photometryObject()
		photometry.times = MJD
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
		slot.target = targetName
		slot.filter = filterName
		slot.aperture = aperture
		slot.headers = headerBlock
		slot.runName = runName
		slot.photometry = photometry
		numSlots = slots.addSlot(slot)
		# print "Added the data to a new slot. Total number of slots is now: %d"%(numSlots)
		print slot
		
	return

def listAllSlots(options):
	slotInfo = slots.getSlotInfo()
	print slotInfo
	return
	
def plot(slotnumber):
	slot = slots.getSlot(slotnumber)
	global plotterHandle
	if (state['plotter'] == 'pgplot'):
		plotterHandle = photmantoPlotting.pgplot(slot, plotterHandle)
	else:
		plotterHandle = photmantoPlotting.matplot(slot, plotterHandle)
	return 
	
def showState():
	print "Current state variables:"
	print state
	return
	
def setState(variable, value):
	state[variable] = value
	print "New state value: %s = %s "%(variable, value)
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
