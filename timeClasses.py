import numpy
import json
import generalUtils
import astropy
import slbarycentric
import sys

class heliocentric:
	telescopes = [
		{ 'name': 'CSS', 'longitude': -110.73167 , 'latitude': 32.41667, 'elev': 2510. },
		{ 'name': 'SSS', 'longitude': -149.1 , 'latitude': -31.3, 'elev': 1150. },
		{ 'name': 'MLS', 'longitude': -110.7889 , 'latitude': 32.4433, 'elev': 2790. }
		]
	
	def __init__(self):
		self.telescope = None
		self.target = None
	
	def setTelescope(self, telescopeName):
		for t in heliocentric.telescopes:
			if t['name'] == telescopeName:
				self.telescope = t
				return True
		return False
		
	def setTarget(self, ra, dec):
		self.target = { 'ra': ra, 'dec': dec}
		
	def convertMJD(self, MJD):
		if self.telescope is None:
			print "We don't know the location of the telescope. Exiting"
			return
			
		obsLocation = astropy.coordinates.EarthLocation(lon = self.telescope['longitude'], lat = self.telescope['latitude'], height = self.telescope['elev'])

		if self.target is None:
			print "We don't know the coordinates of the target. Exiting."

		targetRADEC = generalUtils.toSexagesimal((self.target['ra'], self.target['dec']))
		print "Target position: %s (%f, %f)"%(targetRADEC, self.target['ra'], self.target['dec'])

		targetCoords = astropy.coordinates.SkyCoord(self.target['ra'], self.target['dec'], unit='deg')
		BMJD = []
		t = MJD[0]
		observationTime = slbarycentric.Time(t, format='mjd', location = obsLocation)
		for index, t in enumerate(MJD):
			# observationTime.__init__(t, format='mjd', location = obsLocation)
			observationTime = slbarycentric.Time(t, format='mjd', location = obsLocation)
			delta, bcor = observationTime.bcor(targetCoords)
			bmjd = float(bcor.mjd)
			BMJD.append(bmjd)
			sys.stdout.write("\r[%d/%d]  MJD %5.8f ---> BMJD %5.8f  = %f seconds   "%(index, len(MJD)-1, t, bmjd, delta))
			sys.stdout.flush()
			#print("[%d/%d]  MJD %5.8f ---> BMJD %5.8f  = %f seconds   "%(index, len(MJD)-1, t, bmjd, delta))
			
		print
		
		return BMJD
	

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
		#if phase>0.5: phase = phase - 1.0
		return phase
		
	def getOrbits(self, HJD):
		HJD_difference = HJD - self.T0
		norbits = int(HJD_difference / self.Period)
		upperorbit = int(round(HJD_difference / self.Period))
		return norbits, upperorbit
		
	def getOffsetOrbits(self, HJD):
		HJD_difference = HJD - self.T0 - self.Period/2.0
		norbits = int( HJD_difference / self.Period)
		return norbits
		
	def loadFromFile(self, filename):
		file = open(filename, 'r')
		for line in file:
			if line[0] == '#': continue
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
		raSeconds = float(coords[2])
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
		print "RA:", raDegrees, "DEC:", decCalc
		return raDegrees, decCalc
			
	def __str__(self):
		outString = "T0: %7.8f [%7.8f] + E x %7.10f [%7.10f]"%(self.T0, self.T0_error, self.Period, self.Period_error)
		return outString
			
