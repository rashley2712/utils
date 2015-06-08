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
import math, time

# import astropy.coordinates.EarthLocation

slots = photometryClasses.slotCollection()
state = {	'plotter'	: 'matplot', \
			'overplot'	: False,  \
			'yerrors'	: False,  \
			'plotlimits': 'auto', \
			'plotcolour': 'r',    \
			'xlabel'	: 'auto', \
			'ylabel'	: 'auto', \
			'ps'		: False,  \
		 	'plotterhandle': None}

telescopes = []
telescope = {}
telescope['name'] = "Very Large Telescope, UT3 (Melipal)"
telescope['longitude'] = 289.5972
telescope['latitude'] = -24.6253
telescope['altitude'] = 2635.
telescopes.append(telescope)

print "Telescope information:"
print telescopes

def loadFromLogFile(filename, maxRows = 0):
	""" Loads the photometry data from a log file created by Tom Marsh ULTRACAM pipeline. """
	ultracam = False
	ultraspec = True
	inputFile = open(filename, 'r')
	
	xValues = []
	yValues = []
	frameList = []
	headerBlock = ""
	runName = "--unknown--"
	telescope = "--unknown--"
	targetName = "--unknown--"
	filterName = "--unknown--"
	PI = "--unknown--"
	columnCount = 0
	for line in inputFile:
		if line[0] == '#':
			headerBlock+=line
			if ("target" in line) and ("estimated" not in line):
				targetName = generalUtils.getBetweenChars(line, '=', '/').strip()
				print "Target: %s"%targetName
			if ("filters" in line):
				filterName = generalUtils.getBetweenChars(line, '=', '/').strip()
				print "Filters: %s"%filterName
			if ("Telescope" in line) and ("observing" not in line):
				telescope = generalUtils.getBetweenChars(line, '=', '/').strip()
				print "Telescope: %s"%telescope
			if (" pi " in line):
				PI = generalUtils.getBetweenChars(line, '=', '/').strip()
				print "PI: %s"%PI
			if (" Data file name " in line):
				runName = generalUtils.getBetweenChars(line, '=', '\n').strip()
				print "run data file: %s"%runName
			if (" Server file name " in line):
				runName = generalUtils.getBetweenChars(line, '=', '\n').strip()
				print "run data file: %s"%runName
				
		if line[0] != '#':
			params = line.split()
			# print params
			frameIndex = int(params[0])
			frameList.append(frameIndex)
			columnCount = len(params)
	firstFrame = frameList[0]
	countRecurrence = 0
	for f in frameList:
		if f == firstFrame: countRecurrence+=1
	
	numApertures = int( ((columnCount-7)/14) )
	print "ColumnCount: ", columnCount, "which means %d apertures."%numApertures
	# frameList = generalUtils.removeDuplicatesFromList(frameList)
	print "The run in file %s contains %d frames. Start frame: %d End frame: %d"%(filename, len(frameList), min(frameList), max(frameList))
	if countRecurrence == 3:
		print "This file has 3 CCDs. It is an ULTRACAM file."
		ultracam = True
		ultraspec = False
	if countRecurrence == 1: 
		print "This file has 1 CCD. It is an ULTRASPEC file."
		ultracam = False
		ultraspec = True

	if (ultracam): CCDs = [1, 2, 3]
	else: CCDs = [1]
	for CCD in CCDs: 
		for aperture in range(1, numApertures+1):
			apertureIndex = 14*(aperture-1) + 7
			print "Reading data for aperture %d, CCD %d"%(aperture, CCD)
			inputFile.seek(0)
			MJDs = []
			counts = []
			skys = []
			sigmas = []
			errors = []
			timeFlags = []
			exposures = []
			FWHMs = []
			betas = []
			xs = []
			ys = []
			lineCounter = 0
			for line in inputFile:
				lineCounter+= 1
				sys.stdout.write("\rLine number: %d    "%(lineCounter))
				sys.stdout.flush()
				if line[0] != '#':
					params = line.split()
					# print params
					CCDValue = int(params[4])
					apertureValue = int(params[apertureIndex])
					if CCDValue == CCD: 
						frameIndex = int(params[0])
						MJDs.append(float(params[1]))
						timeFlags.append(int(params[2]))
						exposures.append(float(params[3]))
						FWHMs.append(float(params[5]))
						betas.append(float(params[6]))
						xs.append(float(params[apertureIndex + 1]))
						ys.append(float(params[apertureIndex + 2]))
						counts.append(float(params[apertureIndex + 7]))
						sigmas.append(float(params[apertureIndex + 8]))
						skys.append(float(params[apertureIndex + 9]))
						errors.append(int(params[apertureIndex + 13]))
					
			photometry = {}
			
			photometry['MJD'] = numpy.array(MJDs)
			photometry['exposure'] = numpy.array(exposures)
			photometry['FWHM'] = numpy.array(FWHMs)
			photometry['beta'] = numpy.array(betas)
			photometry['x'] = numpy.array(xs)
			photometry['y'] = numpy.array(ys)
			photometry['counts'] = numpy.array(counts)
			photometry['sigma'] = numpy.array(sigmas)
			photometry['sky'] = numpy.array(skys)
			photometry['error'] = numpy.array(errors)	
		
			id = slots.getNextSlotID()
			print "new ID:", id
			slot = photometryClasses.slotObject(id)
			slot.setPhotometry(photometry)
			slot.setTimeColumn('MJD')
			slot.setYColumn('counts')
			slot.setYError('sigma')
			slot.target = targetName
			slot.filter = filterName
			slot.aperture = aperture
			slot.headers = headerBlock
			slot.runName = runName
			slot.telescope = telescope
			slot.CCD = "CCD %d"%CCD
			numSlots = slots.addSlot(slot)
			print "Added the data to a new slot. Total number of slots is now: %d"%(numSlots)
			print slot
	
	inputFile.close()
	return
	

def loadFromFITSFile(filename, maxRows=0):
	""" Loads photometry data from a FITS file generated by ulog2fits. If maxRows!=0 then it loads up to maxRows."""
	ultracam = False
	ultraspec = False
	inputFile = astropy.io.fits.open(filename)
	fileInfo = inputFile.info()
	
	print fileInfo
	print len(inputFile)
	if len(inputFile)==4:
		print "We have an ULTRACAM file..."
		ultracam = True
	if len(inputFile)==2:
		print "We have an ULTRASPEC file..."
		ultraspec = True
	
	if ultraspec:
		CCDs = ['CCD 1']
	if ultracam:
		CCDs = ['CCD 1', 'CCD 2', 'CCD 3']
	
	
	headerBlock = str(inputFile[0].header)
	
	# Get some header info
	targetName = generalUtils.getKeyValueFromFITSHeader('target', headerBlock)
	filterName = generalUtils.getKeyValueFromFITSHeader('filter', headerBlock)
	runName = generalUtils.getKeyValueFromFITSHeader('Data file name', headerBlock, terminator=' ')
	
	for CCD in CCDs:
		headers = inputFile[CCD].header
		data = inputFile[CCD].data
		columns = inputFile[CCD].columns
	
		allData = []
	
		for index, item in enumerate(data):
			allData.append(item)
			if maxRows!=0 and index>=maxRows-1: break
	
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
		for aperture in range(1, maxApertureIndex+1):
			print "Loading data for aperture #", aperture
		
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
			slot.setYColumn('counts')
			slot.setYError('sigma')
			slot.CCD = CCD
			numSlots = slots.addSlot(slot)
			# print "Added the data to a new slot. Total number of slots is now: %d"%(numSlots)
			print slot
	
	inputFile.close()
		
	return

def listSlots(slotIDs, long=False):
	for s in slotIDs:
		slot = slots.getSlotByID(s)
		if slot:
			if long: print slot.longString() 
			else: print str(slot)
	return

def listAllSlots(long=False):
	if slots.getNumSlots()==0: 
		print "No slots"
		return
	slotInfo = slots.getSlotInfo(long)
	print slotInfo
	return
	
def plot(slotIDs):
	global state, slot
	
	for s in slotIDs:
		slot = slots.getSlotByID(s)
		if slot:
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
	if not slot:
		print "No slot with slot id %d found."%slotID
		return
	columns = slot.getPhotometryColumnList()
	print "Available columns:", columns
	print "Time [%s]"%slot.timeColumn,
	if slot.yColumn!="": print " y-axis [%s]"%slot.yColumn,
	if slot.yError!="": print " y-errors [%s]"%slot.yError,
	print 
	return

def catSlot(slotID):
	slot = slots.getSlotByID(slotID)
	if not slot:
		print "No slot with slot id %d found."%slotID
		return
	slot.catData()
	return
	
def removeZeros(slotID):
	slot = slots.getSlotByID(slotID)
	if not slot:
		print "No slot with slot id %d found."%slotID
		return
	yaxis = slot.yColumn
	yValues = slot.getPhotometryColumn(yaxis)
	print yValues
	mask = numpy.zeros(len(yValues))
	for index, y in enumerate(yValues):
		if y == 0:
			mask[index] = 1
	slot.applyMask(mask)
			
	return
	

def calculateBMJD(slotID):
	slot = slots.getSlotByID(slotID)
	times = slot.getPhotometryColumn('MJD')
	# TNT observatory
	obsLong = 98.48
	obsLat = 18.57
	obsAlt = 2457. 
	# WHT observatory
	# obsLong = 342.1184
	# obsLat = 28.7606
	# obsAlt = 2326. 
	# VLT (Paranal)
	obsLong = 289.5972
	obsLat = -24.6253
	obsAlt = 2635. 
	
	obsLocation = astropy.coordinates.EarthLocation(lon = obsLong, lat = obsLat, height=obsAlt)
	
	targetRASex = "07 11 26"			# CSS081231
	targetDecSex = "+44 04 05"
	# targetRASex = "21 07 58.188"		# HU Aqr  21 07 58.188 -05 17 40.47
	# targetDecSex = "-05 17 40.47"
	targetRASex = "14 09 07.46"			# V834 Cen
	targetDecSex = "-45 17 17.1"
	
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
		sys.stdout.write("\r[%d/%d]  MJD %5.8f ---> BMJD %5.8f  = %f seconds   "%(index, len(times)-1, t, bmjd, delta))
		sys.stdout.flush()
	print
	slot.addColumn("BMJD", numpy.array(BMJD), clobber=True)
	return
	
def divide(slotA_ID, slotB_ID, slotD_ID):
	slotA = slots.getSlotByID(slotA_ID)
	slotB = slots.getSlotByID(slotB_ID)
	slotD = slots.getSlotByID(slotD_ID)
	if not slotA: 
		print "No slot found with slotID:", slotA_ID
		return	
	if not slotB: 
		print "No slot found with slotID:", slotB_ID
		return
	if slotD: 
		print "Destination slot ID: %d is not empty. Aborting."%slotD_ID	
		return
	timesA = slotA.getPhotometryColumn(slotA.timeColumn)
	timesB = slotB.getPhotometryColumn(slotB.timeColumn)
	print "Time axis for slot A (%d) is: %s with a length of %d data points."%(slotA_ID, slotA.timeColumn, len(timesA))
	print "Time axis for slot B (%d) is: %s with a length of %d data points."%(slotB_ID, slotB.timeColumn, len(timesB))
	if len(timesA) != len(timesB):
		print "Lengths don't match. Aborting."
		return
	valuesA = slotA.getPhotometryColumn(slotA.yColumn)
	valuesB = slotB.getPhotometryColumn(slotB.yColumn)
	errorsA = slotA.getPhotometryColumn(slotA.yError)
	errorsB = slotB.getPhotometryColumn(slotB.yError)
	if len(valuesA) != len(valuesB):
		print "Y values are different lengths. Aborting."
		return
	if len(errorsA) != len(errorsB):
		print "Y values are different lengths. Aborting."
		return
	timesD = timesA
	valuesD = []
	errorsD = []
	for (A, B, aA, bB) in zip(valuesA, valuesB, errorsA, errorsB):
		D = A/B
		dD = D * math.sqrt((aA/A)**2 + (bB/B)**2)
		print "A: %f [%f], B: %f [%f] = D: %f [%f]"%(A, aA, B, bB, D, dD)
		valuesD.append(D)
		errorsD.append(dD)
	valuesD = numpy.array(valuesD)
	errorsD = numpy.array(errorsD) 
	
	copySlot(slotA_ID, slotD_ID)
	destinationSlot = slots.getSlotByID(slotD_ID)
	print destinationSlot
	destinationSlot.addColumn(destinationSlot.timeColumn, timesD, clobber=True)
	destinationSlot.addColumn(destinationSlot.yColumn, valuesD, clobber=True)
	destinationSlot.addColumn(destinationSlot.yError, errorsD, clobber=True)
	return
	
def sigmaclip(slotID, factor, groupSize):
	slot = slots.getSlotByID(slotID)
	if not slot:
		print "No slot found with slot ID:", slotID
		return
	values = numpy.array(slot.getPhotometryColumn(slot.yColumn))
	values = numpy.ma.masked_array(values)
	errorFlags = slot.getPhotometryColumn("error")
	for index, e in enumerate(errorFlags):
		if e!=0: 
			print "Error flag at:%d, %d"%(index, e)
			values[index] = numpy.ma.masked
	
	means = []
	sigmas = []
	halfGroup = groupSize/2
	print "Using group size: %d,half group: %d"%(groupSize, halfGroup)
	# time.sleep(2)
	if groupSize!=1:
		for index in range(halfGroup):
			leftIndex = 0
			rightIndex = groupSize
			subSet = values[leftIndex:rightIndex+1]
			# print "index: %d  range [%d-%d] -> %s"%(index, leftIndex, rightIndex, str(subSet))
			mean = numpy.mean(subSet)
			sigma = numpy.std(subSet)
			means.append(mean)
			sigmas.append(sigma)
			
		# time.sleep(5)	
		
		for index in range(halfGroup, len(values) - halfGroup):
			leftIndex = index-halfGroup
			rightIndex = index+halfGroup
			subSet = values[leftIndex:rightIndex+1]
			# print "index: %d  range [%d-%d] -> %s"%(index, leftIndex, rightIndex, str(subSet))
			mean = numpy.mean(subSet)
			sigma = numpy.std(subSet)
			means.append(mean)
			sigmas.append(sigma)
		# time.sleep(5)
	
		for index in range(len(values)-halfGroup, len(values)):
			leftIndex = len(values) - groupSize
			rightIndex = len(values)
			subSet = values[leftIndex:rightIndex+1]
			# print "index: %d  range [%d-%d] -> %s"%(index, leftIndex, rightIndex, str(subSet))
			mean = numpy.mean(subSet)
			sigma = numpy.std(subSet)
			means.append(mean)
			sigmas.append(sigma)
	
	pointsMask = numpy.zeros(len(values))
	indices = range(len(values))	
	for (index, mean, sigma, value) in zip(indices, means, sigmas, values):
		print index, "Mean:", mean, "Sigma:", sigma, "Value:", value
		distanceFromMean = abs(value - mean)
		if distanceFromMean > (sigma * factor):
			print "Rejecting this point"
			pointsMask[index] = 1
			
	print pointsMask
	
	slot.applyMask(pointsMask)
	
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

def calculateMean(slotID):
	slot = slots.getSlotByID(slotID)
	yValues = slot.getPhotometryColumn(slot.yColumn)
	yErrors = slot.getPhotometryColumn(slot.yError)
	weights = 1.0/yErrors
	mean = numpy.average(yValues, weights = weights)
	print "Weighted mean:", mean
	print "Std dev:", numpy.std(yValues)
	return


def calculateCountRate(slotID):
	slot = slots.getSlotByID(slotID)
	
	counts = slot.getPhotometryColumn('counts')
	sigmas = slot.getPhotometryColumn('sigma')
	exposures = slot.getPhotometryColumn('exposure')
	if len(counts) == 0:
		print "Could not locate the 'counts' column. Aborting operation."
		return
	if len(exposures) == 0:
		print "Could not locate the 'exposure' column. Aborting operation."
		return
	if len(exposures) != len(counts):
		print "For some reason the lengths of the 'exposure' and the 'counts' column are different lengths. Need to abort the operation."
		return
	countrates = []
	sigma_crs = []
	for s, c, e in zip(sigmas, counts, exposures):
		cr = c / e
		sigma_cr = s / e
		print "Count: %f [%f], Exposure: %f seconds, CountRate: %f [%f] c/s"%(c, s, e, cr, sigma_cr)
		countrates.append(cr)
		sigma_crs.append(sigma_cr)
		
	slot.addColumn('countrate', countrates)
	slot.addColumn('sigma_cr', sigma_crs)
	return


def showHeader(slotID):
	slot = slots.getSlot(slotID)
	header = slot.headers
	print header	
	return
	
def setSlotProperty(slotID, property, value):
	slot = slots.getSlotByID(slotID)
	if property=='xaxis':
		if (slot.setTimeColumn(value)): print "setting xaxis to " + value
		else: print value, "is not a valid column in slot ", slotID
	if property=='yaxis':
		if (slot.setYColumn(value)): print "setting yaxis to " + value
		else: print value, "is not a valid column in slot ", slotID
	if property=='yerrors':
		if (slot.setYError(value)): print "setting yerrors to " + value
		else: print value, "is not a valid column in slot ", slotID
	setattr(slot, property, value)
	return
	
def copySlot(fromID, toID):
	if (not slots.exists(fromID)):
		print "Aborting operation: No slot known with source ID: %d"%fromID
		return
	if (slots.exists(toID)):
		print "Aborting operation: A slot already exists with destination ID: %d"%toID
		return
	oldSlot = slots.getSlot(fromID)
	newSlot = copy.deepcopy(oldSlot)
	newSlot.id = toID
	slots.addSlot(newSlot)
	return

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='General purpose tool for loading, manipulating and plotting ULTRACAM and ULTRASPEC photometry data.')
	parser.add_argument('script', type=str, nargs='?', help='The name of a script file containing commands to execute.')
	arg = parser.parse_args()
	# print arg
	
	commands = commandModule.photCommands
	
	if arg.script != None:
		if arg.script=="restore":
			commands().do_restore("")
		else: 
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
