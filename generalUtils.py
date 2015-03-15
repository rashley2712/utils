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
	
def fromSexagesimal(coords):
	raStr = split(coords[0])
	
	return (0,0)

def convertMJDtoHJD(MJD, coords):
	ra = coords[0]
	dec = coords[1]
	
	
