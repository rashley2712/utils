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
	dec = decHours + decMinutes/60.0 + decSeconds / 3600.0
	if south:
		dec = decHours - decMinutes/60.0 - decSeconds / 3600.0
		
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
	
