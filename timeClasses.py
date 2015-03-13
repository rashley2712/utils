import numpy
import json

class ephemerisObject:
	def __init__(self):
		self.T0 = 0
		self.T0_error = 0
		self.Period = 0 
		self.Period_error = 0
		
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
			
			
	def __str__(self):
		outString = "T0: %7.8f [%7.8f] + E X %7.10f [%7.10f]"%(self.T0, self.T0_error, self.Period, self.Period_error)
		return outString
			
