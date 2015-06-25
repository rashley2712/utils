import math

def toSexagesimal(world):
	raDeg = world[0]
	ra = raDeg/15.
	hours = int(ra)
	minutes = (ra - int(ra)) * 60
	seconds = (minutes - int(minutes)) * 60
				
	dec = world[1]
	decDegrees = int(dec)
	decMinutes = (dec - int(dec)) * 60
	decSeconds = (decMinutes - int(decMinutes)) * 60
		
	outString = "RA: %02d:%02d:%02.1f"%(hours, minutes, seconds)
	outString+= " DEC: %02d:%02d:%02.3f"%(dec, decMinutes, decSeconds)
	return outString

def filterOutNaNs(data):
	""" Filter out NaN entries from a dictionary containing any number of arrays (all of the same length)"""
	newData = {}
	for key in data.keys():
		newData[key] = []
	for index, d in enumerate(data[key]):
		for key in data.keys():
			value = data[key][index]
			if not math.isnan(value):
				newData[key].append(value)
	return newData			
	
def parseIntegerList(line):
	""" Takes a string containing a list directive and turns it into a list of valid integers. Removes duplicates.
	eg 1,2,4,5  --->  [1, 2, 4, 5]
	   1-5,8    --->  [1, 2, 3, 4, 5, 8]"""
	   
	intList = []
	# Remove any spaces
	line = line.replace(' ', '')
	# Take each comma separated part...
	parts = line.split(',')
	for p in parts:
		# If a range a-b is specified...
		if '-' in p:
			try:
				limits = p.split('-')
				start = int(limits[0])
				end = int(limits[1])
				numbers = range(start, end+1)
				for number in numbers:
					intList.append(number)
			except ValueError:
				print "Didn't understand the range you entered."
		# If just a single number is specified
		else:
			try:
				number = int(p)
				intList.append(number)
			except ValueError:
				print "Didn't understand the number you entered."
	intList = removeDuplicatesFromList(intList)
	return intList
	
def removeDuplicatesFromList(list):
	newList = []
	for index, l in enumerate(list):
		if l not in list[:index]:
			newList.append(l)
	return newList
	
def getBetweenChars(inputString, startChar, endChar):
	""" Gets a portion of a string between two characters """
	returnString = ""
	if (startChar not in inputString):
		return returnString
	if (endChar not in inputString):
		return returnString
	startPosition = inputString.find(startChar)
	endPosition = inputString[startPosition:].find(endChar)
	returnString = inputString[startPosition+1:startPosition+endPosition]
	return returnString

def getBetweenStrings(inputString, startString, endString):
	""" Gets a portion of a string between two strings """
	returnString = ""
	if (startString not in inputString):
		return returnString
	if (endString not in inputString):
		return returnString
	startPosition = inputString.find(startString)
	endPosition = inputString[startPosition:].find(endString)
	returnString = inputString[startPosition+len(startString):startPosition+endPosition]
	return returnString
		
	
def fromSexagesimal(raStr, decStr):
	""" Format for input ra and dec are 'HH:MM:SS.dd' and 'nDD:MM:SS.dd'
									or 	'HH MM SS.dd' and 'nDD MM SS.dd'
	"""
	separator = ':'
	if raStr.find(separator)==-1:
		separator = ' '
	raPieces = raStr.split(separator)
	raHours = int(raPieces[0])
	raMinutes = int(raPieces[1])
	raSeconds = float(raPieces[2])
	ra = 15 * (raHours + raMinutes/60.0 + raSeconds / 3600.0)
	
	decPieces = decStr.split(separator)
	if decPieces[0][0]=='-':
		south = True
	else:
		south = False
		
	decHours = int(decPieces[0])
	decMinutes = int(decPieces[1])
	decSeconds = float(decPieces[2])
	
	if south:
		dec = decHours - decMinutes/60.0 - decSeconds / 3600.0
	else:
		dec = decHours + decMinutes/60.0 + decSeconds / 3600.0
			
	return (ra, dec)

def convertMJDtoHJD(MJD, coords):
	ra = coords[0]
	dec = coords[1]
	
def getKeyValueFromFITSHeader(key, headerBlock, terminator='/'):
	""" Get the value from a 'key' = 'value' comment in a FITS header. Returns 'Unknown' if not found"""
	
	valueIndex = headerBlock.find(key)
	if valueIndex != -1:
		valueIndex = headerBlock.find('=', valueIndex) + 2
		endIndex = headerBlock.find(terminator, valueIndex)
		value = headerBlock[valueIndex:endIndex]
		value = value.rstrip()
	else:
		value = "Unknown"
		
	return value
	
