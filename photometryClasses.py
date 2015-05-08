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
		#idList = [s['id'] for s in self.slotList]
		#print "IDs"
		orderedSlotList = sorted(self.slotList, key=lambda object: object.id)
		for index, s in enumerate(orderedSlotList):
			retString+= "%s"%str(s) + "\n"
		return retString
		
	def getSlot(self, number):
		return self.slotList[number]
		
	def getSlotByID(self, slotID):
		slotIDList = [s.id for s in self.slotList]
		try:
			index = slotIDList.index(slotID)
			return self.slotList[index]
		except ValueError:
			return False
		
	def exists(self, number):
		""" Checks to see if the slotID already exists """
		slotIDList = [s.id for s in self.slotList]
		try:
			index = slotIDList.index(number)
			print number," id exists in slot index:", index
			return True
		except ValueError:
			return False
		return False
		
	def replace(self, slot):
		slotIDList = [s.id for s in self.slotList]
		index = slotIDList.index(slot.id)
		self.slotList[index] = slot	
		
	def getNextSlotID(self):
		if len(self.slotList)==0: return 0
		slotIDList = [s.id for s in self.slotList]
		nextID = max(slotIDList) + 1
		return nextID
			

class slotObject:
	""" A class containing time-series photometry for an object... similar in concept to Tom Marsh's slot in his Molly software
	"""
	def __init__(self, id):
		self.id = id
		self.channels = []   # A list of channel descriptions (eg 'red', 'green', 'blue')
		self.target = "None" # Name of the target object
		self.object = None   # An object of class targetObject containing meta-data about the target
		self.photometry = {}
		self.headers = "No header information loaded."
		self.filter = "n/a"
		self.aperture = 0
		self.runName = ""	
		self.photometryColumns = []
		self.timeColumn = ""
		self.yColumn = ""
		self.yError = ""
		self.times = []
		
	def initFromJSON(self, jsonObject):
		self.target = jsonObject['target']
		self.headers = jsonObject['headers']
		self.filter = jsonObject['filter']
		self.aperture = jsonObject['aperture']
		self.runName = jsonObject['runName'] 
		self.columns = jsonObject['columns']
		self.timeColumn = jsonObject['timecolumn']
		try:
			self.yColumn = jsonObject['ycolumn']
		except:
			self.yColumn = ""
		try:
			self.yError = jsonObject['yerror']
		except:
			self.yError = ""
			
		
	def getColumn(self, columnName):
		valueDescriptions = self.photometry.valueDescriptions
		index = valueDescriptions.index(columnName)
		data = [v[index] for v in self.photometry.values]
		return data
		
	def setPhotometry(self, photometry):
		self.photometry = photometry
		self.getPhotometryColumnList()
		
	def getPhotometryColumn(self, columnName):
		return self.photometry[columnName]
		
	def setTimeColumn(self, columnName):
		self.times = self.photometry[columnName]
		self.timeColumn = columnName
		return True
		
	def getPhotometryColumnList(self):
		columns = []
		for key in self.photometry.keys():
			columns.append(key)
		self.photometryColumns = columns
		return columns
		
	def getColumnLengths(self):
		columnLengths = {}
		for c in self.photometryColumns:
			length = len(self.photometry[c])
			columnLengths[c] = length
		return columnLengths
		
	def toJSON(self):
		me = {}
		me['id'] = self.id
		me['target'] = self.target
		me['headers'] = self.headers
		me['filter'] = self.filter
		me['aperture'] = self.aperture
		me['runName'] = self.runName
		self.getPhotometryColumnList()
		me['columns'] = self.photometryColumns
		me['timecolumn'] = self.timeColumn
		me['ycolumn'] = self.yColumn
		me['yerrorcolumn'] = self.yError
		
		
		for c in self.photometryColumns:
			me[c] = self.photometry[c].tolist()
		return json.dumps(me)	
	
	def __str__(self):
		retStr = "ID: %d Run file: %s \tTarget: %s \tFilter: %s \tAperture: %d \t Length: %d"%(self.id, self.runName, self.target, self.filter, self.aperture, len(self.photometry['MJD']))
		return retStr

	
class photometryObject:
	def __init__(self):
		self.times = []
		self.values = []
		self.timeDescription = "Unknown"
		self.valueDescriptions = []
		self.dataLength = 0
		
	def addValueDescription(self, description):
		self.valueDescriptions.append(description)
		return len(self.valueDescriptions)
		
	def setPhotometry(self, times, values, timeDescription, valueDescription):
		self.times = times
		self.values = values
		self.timeDescription = timeDescription
		self.valueDescription = valueDescription
		
	def addData(self, data):
		self.values+= data
		self.dataLength = len(self.values)
		return self.dataLength
		

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
			
