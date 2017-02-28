import math, numpy, os, sys
from PIL import Image,ImageDraw,ImageFont

def percentiles(data, lo, hi):
    """ Returns a normalised array where lo percent of the pixels are 0 and hi percent of the pixels are 255
    """
    max = data.max()
    dataArray = data.flatten()
    pHi = numpy.percentile(dataArray, hi)
    pLo = numpy.percentile(dataArray, lo)
    range = pHi - pLo
    scale = range/255
    data = numpy.clip(data, pLo, pHi)
    data-= pLo
    data/=scale
    return data

def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")


def formatValueError(value, error, significant = 1):
	""" Returns two neatly formatted numbers (as strings) taking into account an appropriate number of significant figures based on the error. """
	if error==0: decimalPlaces = 1
	else: decimalPlaces = round(math.log10(error))
	if decimalPlaces < 1:
		decimalPlaces = abs(decimalPlaces) + (significant)
		value = round(value * 10 ** decimalPlaces) / 10**decimalPlaces
		error = round(error * 10 ** decimalPlaces) / 10**decimalPlaces
	
	valueStr = str(value)
	errorStr = str(error)
	
	return valueStr, errorStr

	
def changeExtension(filename, extension):
	return os.path.splitext(filename)[0] + "." + extension 
	
def writePNG(imageArray, filename, caption = ""):
	""" Writes to a PNG file using the PIL library. Adds a caption if sent in the parameters. Also adds a .png extension if it isn't already there in 'filename' """
	imgData = numpy.rot90(imageArray, 3)
	imgSize = numpy.shape(imgData)
	imgLength = imgSize[0] * imgSize[1]
	testData = numpy.reshape(imgData, imgLength, order="F")
	img = Image.new("L", imgSize)
	palette = []
	for i in range(256):
		palette.extend((i, i, i)) # grey scale
		img.putpalette(palette)
	img.putdata(testData)
	outputFilename = changeExtension(filename, "png")
	print ("Writing PNG file: " + outputFilename) 
	img.save(outputFilename, "PNG", clobber=True)
	
	"""image = Image.fromarray(imageArray)
	imageCopy = image.copy()    # We need to copy the image so we don't alter the original when adding the caption.
	if (caption!=""): 
		font = ImageFont.truetype(config.FONT, 25) 
		draw = ImageDraw.Draw(imageCopy)
		if (imageCopy.mode == "L"):
			draw.text((0, 0), caption, 255, font = font)
		else: 
			draw.text((0, 0), caption, (255, 255, 255), font = font)
	
	if (filename[-4:]!=".png"): filename+= ".png"
	imageCopy.save(filename, "PNG")"""
	
def toSexagesimal(world):
	raDeg = world[0]
	ra = raDeg/15.
	hours = int(ra)
	minutes = (ra - int(ra)) * 60
	seconds = (minutes - int(minutes)) * 60
				
	dec = world[1]
	decDegrees = int(dec)
	if dec>0: decSign = "+"
	else: decSign = "-"
	decMinutes = (dec - int(dec)) * 60
	decSeconds = (decMinutes - int(decMinutes)) * 60
		
	outString = "RA: %02d:%02d:%04.1f"%(hours, minutes, seconds)
	outString+= " DEC: %s%02d:%02d:%04.3f"%(decSign, dec, decMinutes, decSeconds)
	return outString
	
def writeFriendlyTimeSeconds(seconds):
	""" Writes a friendly time (in hours and/or minutes) based on an input of seconds
	"""
	minutes = seconds / 60.
	timeStr = str(int(minutes)) + " minutes"
	if minutes > 60: 
		hours = minutes/60.
		minutes = minutes - int(hours)*60
		timeStr= "%d hour"%(int(hours))
		if int(hours)!=1: timeStr+="s";
		timeStr+= ", %d minute"%(int(minutes))
	if int(minutes)!=1: timeStr+="s";
			
	return timeStr

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
	
