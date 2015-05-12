#!/usr/bin/env python
import sys
import commandModule
import astropy.io.fits
import argparse
import photometryClasses
import generalUtils
import photmantoPlotting
import json, numpy
import copy
import slbarycentric
# import astropy.coordinates.EarthLocation

slots = photometryClasses.slotCollection()
state = {	'plotter'	: 'matplot', \
			'overplot'	: False, \
			'yerrors'	: False, \
			'plotlimits': 'auto', \
			'plotcolour': 'r', \
		 	'plotterhandle': None}

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
	for aperture in range(1, maxApertureIndex+1):
		print "Loading data for aperture #", aperture
		
		""" Try a different approach to loading this stuff """
		photometry = {}
		photometry['MJD'] = 		data.field('MJD')
		photometry['exposure'] = 	data.field('Expose')
		photometry['FWHM'] = 		data.field('FWHM')
		photometry['beta'] = 		data.field('beta')
		photometry['x'] = 			data.field('X_' + str(aperture))
		photometry['y'] = 			data.field('Y_' + str(aperture))
		photometry['counts'] = 		data.field('Counts_' + str(aperture))
		photometry['sigma'] = 		data.field('Sigma_' + str(aperture))
		photometry['sky'] = 		data.field('Sky_' + str(aperture))
		photometry['sigma'] = 		data.field('Sigma_' + str(aperture))
		photometry['error'] = 		data.field('Eflag_' + str(aperture))
		
		id = slots.getNextSlotID()
		print "new ID:", id
		slot = photometryClasses.slotObject(id)
		slot.channels = ['ULTRASPEC']
		slot.target = targetName
		slot.filter = filterName
		slot.aperture = aperture
		slot.headers = headerBlock
		slot.runName = runName
		slot.setPhotometry(photometry)
		slot.setTimeColumn('MJD')
		numSlots = slots.addSlot(slot)
		# print "Added the data to a new slot. Total number of slots is now: %d"%(numSlots)
		print slot
		
	return

def listAllSlots(options):
	if slots.getNumSlots()==0: 
		print "No slots"
		return
	if options=="-l":
		slotInfo = slots.getSlotInfo(long = True)
	else:
		slotInfo = slots.getSlotInfo()
	print slotInfo
	return
	
def plot(slotID):
	global state, slot
	
	slot = slots.getSlotByID(slotID)
	if (state['plotter'] == 'pgplot'):
		state = photmantoPlotting.pgplot(slot, plotterHandle)
	else:
		state = photmantoPlotting.matplot(slot, state)
	return 
	
def showState():
	print "Current state variables:"
	for key in state.keys():
		print "%s \t: \t%s \t\t\t[%s]"%(key, state[key], type(state[key]).__name__)
	return
	
def setState(variable, value):
	state[variable] = value
	if value == 'on':
		state[variable] = True
	if value == 'off':
		state[variable] = False
	print "New state value: %s = %s "%(variable, value)
	return
	
def saveSession(filename):
	if filename==None or filename=="":
		sessionFilename = 'session.ptm'
		dataFilename = 'data.ptm'
	else:
		sessionFilename+= '.session.ptm'
		dataFilename+= '.data.ptm'
		
	print "Saving session to file:", sessionFilename
	outputfile = open(sessionFilename, "w")
	json.dump(state, outputfile)
	outputfile.close()
	
	print "Saving slots to file:", dataFilename
	outputfile = open(dataFilename, "w")
	
	slotData = []
	for s in slots.slotList:
		dataObject =  s.toJSON()
		slotData.append(dataObject)
	json.dump(slotData, outputfile)
	outputfile.close()
	
	return

def restoreSession(filename):
	global state
	if filename==None or filename=="":
		sessionFilename = 'session.ptm'
		dataFilename = 'data.ptm'
	else:
		sessionFilename+= '.session.ptm'
		dataFilename+= '.data.ptm'
		
	print "Loading session from file:", sessionFilename
	inputfile = open(sessionFilename, "r")
	stateObject = json.load(inputfile)
	for key in stateObject.keys():
		keyString = str(key)
		value = stateObject[key]
		if type(value) is unicode: 
			value = str(value)
		state[keyString] = value
	inputfile.close()
	state['plotterhandle'] = None
	showState()
	
	print "Loading slots from file:", dataFilename
	inputfile = open(dataFilename, "r")
	
	slotData = []
	allData = json.load(inputfile)
	for index, s in enumerate(allData):
		loadedSlot = json.loads(s)
		loadedSlotID = loadedSlot["id"]
		newSlot = photometryClasses.slotObject(loadedSlotID)
		newSlot.initFromJSON(loadedSlot)
		photometry = {}
		photometryColumns = loadedSlot['columns']
		for c in photometryColumns:
			photometry[str(c)] = numpy.asarray(loadedSlot[c])
		newSlot.setPhotometry(photometry)
		if (slots.exists(loadedSlotID)): 
			print "Warning: We are about to replace the slot with slot ID:%d"%loadedSlotID
			slots.replace(newSlot)
		else:
			print "New slot created: %d - for: %s"%(loadedSlotID, newSlot)
			slots.addSlot(newSlot)
	inputfile.close()
	
	return
	
def showColumns(slotID):
	slot = slots.getSlotByID(slotID)
	columns = slot.getPhotometryColumnList()
	print "Available columns:", columns
	print "Time [%s]"%slot.timeColumn,
	if slot.yColumn!="": print " y-axis [%s]"%slot.yColumn,
	if slot.yError!="": print " y-errors [%s]"%slot.yError,
	print 
	return

def calculateBMJD(slotID):
	slot = slots.getSlotByID(slotID)
	times = slot.getPhotometryColumn('MJD')
	obsLong = 98.48
	obsLat = 18.57
	obsAlt = 2457.
	obsLocation = astropy.coordinates.EarthLocation(lon = obsLong, lat = obsLat, height=obsAlt)
	targetRASex = "07 11 26"
	targetDecSex = "+44 04 05"
	ra, dec = generalUtils.fromSexagesimal(targetRASex, targetDecSex)
	print "Calculating barycentric MJD or BMJD"
	print "Observatory location: Lat: %f [deg] Long: %f [deg] Height: %f [m]"%(obsLat, obsLong, obsAlt)
	print "Target position: %s, %s (%f, %f)"%(targetRASex, targetDecSex, ra, dec)
	targetCoords = astropy.coordinates.SkyCoord(ra, dec, unit='deg')
	BMJD = []
	t = times[0]
	observationTime = slbarycentric.Time(t, format='mjd', location = obsLocation)
	for index, t in enumerate(times):
		observationTime.__init__(t, format='mjd', location = obsLocation)
		delta, bcor = observationTime.bcor(targetCoords)
		bmjd = float(bcor.mjd)
		BMJD.append(bmjd)
		sys.stdout.write("\r[%d/%d]  MJD %5.8f ---> BMJD %5.8f  = %f seconds"%(index, len(times)-1, t, bmjd, delta))
		sys.stdout.flush()
	print
	slot.addColumn("BMJD", numpy.array(BMJD))
	return
	
def writeCSV(filename, slotID):
	slot = slots.getSlotByID(slotID)
	outputfile = open(filename, 'w')
	timeLabel = slot.timeColumn
	yLabel = slot.yColumn
	errorLabel = slot.yError
	print timeLabel, yLabel, errorLabel
	outputfile.write("%s, %s, %s\n"%(timeLabel, yLabel, errorLabel))
	xValues = slot.getPhotometryColumn(timeLabel)
	yValues = slot.getPhotometryColumn(yLabel)
	errorValues = slot.getPhotometryColumn(errorLabel)
	
	csvValues = zip(xValues, yValues, errorValues)
	for c in csvValues:
		outputfile.write("%5.8f, %s, %s\n"%(c[0], str(c[1]), str(c[2])))
	
	outputfile.close()
	
	return
	
def calculateMinutes(slotID):
	slot = slots.getSlotByID(slotID)
	times = slot.getPhotometryColumn(slot.timeColumn)
	beginTime = min(times)
	minutesArray = []
	for t in times:
		minutes = (t - beginTime)*24*60
		print t, "converted to", minutes
		minutesArray.append(minutes)
	slot.addColumn('minutes', minutesArray)
	return

def showHeader(slotID):
	slot = slots.getSlot(slotID)
	header = slot.headers
	print header	
	return
	
def setSlotProperty(slotID, property, value):
	slot = slots.getSlotByID(slotID)
	if property=='xaxis':
		slot.setTimeColumn(value)
		print "setting xaxis"
	if property=='yaxis':
		if (slot.setYColumn(value)): print "setting yaxis to " + value
		else: print value, " not found."
	setattr(slot, property, value)
	return
	
def copySlot(fromID, toID):
	if (not slots.exists(fromID)):
		print "Aborting opteration: No slot known with source ID: %d"%fromID
		return
	if (slots.exists(toID)):
		print "Aborting opteration: A slot already exists with destination ID: %d"%toID
		return
	oldSlot = slots.getSlot(fromID)
	newSlot = copy.copy(oldSlot)
	newSlot.id = toID
	slots.addSlot(newSlot)
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
