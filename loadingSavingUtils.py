import astropy.io.fits
import json, csv

def loadFITSFile(filename, *arguments, **keywords):

	inputFile = astropy.io.fits.open(filename)

	requestedColumns = ['MJD', 'Counts_1', 'Counts_2', 'Sigma_1', 'Sigma_2']
	colours = ['r', 'g', 'b']
	
	for kw in keywords.keys(): 
		if kw=='colours':
			colours = keywords[kw]
		if kw=='columns':
			requestedColumns = keywords[kw]

	CCDs = { 'r': 'CCD 1', 'g': 'CCD 2', 'b': 'CCD 3'}
	
	
	allPhotometry = {}
	
	for c in colours:
		headers = inputFile[CCDs[c]].header
		data = inputFile[CCDs[c]].data
		columns = inputFile[CCDs[c]].columns
		
		channelPhotometry = []
		for item in data:
			record = {}
			for col in requestedColumns:
				value = item[columns.names.index(col)]
				record[col] = value
			
			channelPhotometry.append(record)
		allPhotometry[c] = channelPhotometry
	
	print "Loaded data from fits file"
	print "Columns:", requestedColumns
	
	for c in colours:
		numDataPoints = len(allPhotometry[c])
		print c, "has", numDataPoints, "data points"
	
	return allPhotometry

def removeNegativeValues(photometry):
	newPhotometry = []
	for d in photometry:
		if d['Counts_1']<0: continue
		if d['Counts_2']<0: continue
		newPhotometry.append(d)
	numRemoved = len(photometry) - len(newPhotometry)
	return newPhotometry, numRemoved

def removeZeroValues(photometry):
	newPhotometry = []
	for d in photometry:
		if d['Counts_1']==0: continue
		if d['Counts_2']==0: continue
		newPhotometry.append(d)
	numRemoved = len(photometry) - len(newPhotometry)
	return newPhotometry, numRemoved
		

def loadCSV(filename):
	photometry = {}
	colours = ['r', 'g', 'b']
	photometry['r'] = []
	photometry['g'] = []
	photometry['b'] = []
	
	csvfile = open(filename, 'r')
	reader = csv.reader(csvfile, delimiter=',')
	headings = reader.next()
	print headings
	columns = []
	for h in headings:
		columnName = h.strip()
		columns.append(columnName)
		
	for line in reader:
		values = [v.strip() for v in line]
		print values
		colour = values[0]
		record = {}
		for index in range(len(columns)-1):
			record[columns[index+1]] = float(line[index+1])
		print record
		channelPhotometry = photometry[colour]
		channelPhotometry.append(record)
		photometry[colour] = channelPhotometry
		 
	for c in colours:
		numDataPoints = len(photometry[c])
		print c, "has", numDataPoints, "data points"
		
	return photometry
	

def writeCSV(filename, allData, **keywords):
	outputfile = open(filename, "w")
	
	colours = ['r', 'g', 'b']
	
	for kw in keywords.keys(): 
		if kw=='colours':
			colours = keywords[kw]
		if kw=='columns':
			requestedColumns = keywords[kw]
	
	redPhotometry = allData['r']
	
	headerString = 'channel'
	keys = []
	for key, value in redPhotometry[0].items():
		headerString+= ', ' + key
		keys.append(key)
		
	outputfile.write(headerString)
	outputfile.write('\n')

	for c in colours:
		photometry = allData[c]
		for p in photometry:
			outStr = c
			for k in keys:
				outStr+= ", " + str(p[k])
			outputfile.write(outStr)
			outputfile.write('\n')
		
	
	outputfile.close()
		
	

def findMatchingTime(data, target):
	reading = (0, -1)
	
	for d in data:
		time = d[0]
		if time==target:
			reading = (time, d[1])
			continue
	
	return reading

def findClosestTime(data, target):
	distance = 1000.
	reading = (0, 0)
	for d in data:
		time = d[0]
		gap = time - target
		if abs(gap) < distance:
			distance = abs(gap)
			reading = (time, d[1])
	return reading
	
def binData(x, y, bin):
	newX = []
	newY = []
	
	newLength = len(x)/bin
	
	for i in range(newLength):
		xVal = 0
		yVal = 0
		for j in range(bin):
			xVal+= x[i*bin +j]
			yVal+= y[i*bin +j]
		xVal = xVal/bin
		yVal = yVal/bin
		newX.append(xVal)
		newY.append(yVal)
	
	return newX, newY
	
def filterData(dataArray):
	cleanData = []
	
	# This step removes the negative values
	for d in dataArray:
		items = len(d)
		if d[1]<0: continue
		if items>2:
			if d[2]<0: continue
		cleanData.append(d)
		
	return cleanData
	
def calcStats(data):
	np = numpy.array(data)
	values = np[:, 1]
	return (numpy.mean(values), numpy.std(values))
		
