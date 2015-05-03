import numpy
import json


class slotCollection:
	def __init__(self):
		self.slotList = []
		
	def addSlot(self, slot):
		self.slotList.append(slot)
		return len(self.slotList)
		
	def getSlotInfo(self):
		retString = "Slot info...\n"
		for index, s in enumerate(self.slotList):
			retString+= "%d : %s"%(index, s.target) + "\n"
		return retString
			
		

class slotObject:
	""" A class containing time-series photometry for an object... similar in concept to Tom Marsh's slot in his Molly software
	"""
	def __init__(self):
		self.channels = []   # A list of channel descriptions (eg 'red', 'green', 'blue')
		self.target = "None" # Name of the target object
		self.object = None   # An object of class targetObject containing meta-data about the target
		self.photometry = {}
		
	def addPhotometry(self, photometry, channel):
		self.photometry['channel'] = photometry
		self.channels.append(channel)
		
class photometryObject:
	def __init__(self):
		self.times = []
		self.values = []
		self.timeDescription = "Unknown"
		self.valueDescriptions = []
		
	def addValueDescription(self, description):
		self.valueDescriptions.append(description)
		return len(self.valueDescriptions)
		
	def setPhotometry(self, times, values, timeDescription, valueDescription):
		self.times = times
		self.values = values
		self.timeDescription = timeDescription
		self.valueDescription = valueDescription
		
	def addData(self, data):
		self.values.append(data)
		

class ephemerisObject:
	def __init__(self):
		self.T0 = 0
		self.T0_error = 0
		self.Period = 0 
		self.Period_error = 0
		self.ra = 0
		self.dec = 0
		
	def setCoords(self, ra, dec):
		self.ra = ra
		self.dec = dec
		
		
	def getPhase(self, HJD):
		HJD_difference = HJD - self.T0
		#norbits = int( HJD_difference / self.Period)
		phase = (HJD_difference % self.Period) / self.Period
		return phase
		
	def getOrbits(self, HJD):
		HJD_difference = HJD - self.T0
		norbits = int( HJD_difference / self.Period)
		return norbits
		
	def loadFromFile(self, filename):
		file = open(filename, 'r')
		for line in file:
			tokens = line.split()
			if len(tokens)>1: 
				if tokens[0] == 'T0': self.T0 = float(tokens[1])
				if tokens[0] == 'T0_error': self.T0_error = float(tokens[1])
				if tokens[0] == 'E': self.Period = float(tokens[1])
				if tokens[0] == 'E_error': self.Period_error = float(tokens[1])
				if tokens[0] == 'J2000': 
					coords = tokens[1:]
					self.ra, self.dec = self.parseCoords(coords)
					
			
	def parseCoords(self, coords):
		print "Given coords:", coords
		raHours = int(coords[0])
		raMinutes = int(coords[1])
		raSeconds = int(coords[2])
		raFraction = float(raSeconds)/3600. + float(raMinutes)/60.
		raTotal = float(raHours) + raFraction
		raDegrees = raTotal * 15.0
		decDegreesString = str(coords[3])
		sign = decDegreesString[0]
		decDegrees = abs(float(coords[3]))
		decMinutes = float(coords[4])
		decSeconds = float(coords[5])
		decCalc = decDegrees + decMinutes/60. + decSeconds/3600.
		if sign == '-':
			decCalc = -1.0 * decCalc
		#print "RA:", raDegrees, "DEC:", decCalc
		return raDegrees, decCalc
			
	def __str__(self):
		outString = "T0: %7.8f [%7.8f] + E X %7.10f [%7.10f]"%(self.T0, self.T0_error, self.Period, self.Period_error)
		return outString
			
