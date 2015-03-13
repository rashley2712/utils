import numpy
import json

class spectrumObject:
	def __init__(self):
		self.wavelengths = []
		self.flux = []
		self.length = 0 
		self.wavelengthRange = (0, 0)
		
	def setData(self, wavelengths, flux):
		if len(wavelengths) != len(flux):
			return -1
		self.wavelengths = []
		self.flux = []
		for w, f in zip(wavelengths, flux):
			self.wavelengths.append(w)
			self.flux.append(f)
		self.length = len(wavelengths)
		self.wavelengthRange = (min(wavelengths), max(wavelengths))
		
		return self.length
		
	def trimWavelengthRange(self, lower, upper):
		newWavelengths = []
		newFlux = []
	
		if lower>=upper: return self.length
		
		for w, f in zip(self.wavelengths, self.flux):
			if (w>lower) and (w<upper):
				newWavelengths.append(w)
				newFlux.append(f)
		self.length = len(newWavelengths)
		self.wavelengths = newWavelengths
		self.flux = newFlux
		self.wavelengthRange = (min(wavelengths), max(wavelengths))
			
		return self.length		
		
	def getWavelengths(self):
		return self.wavelengths
		
	def getFlux(self):
		return self.flux
		
	def getNearestFlux(self, wavelength):
		minDistance = self.wavelengthRange[1]
		for w, f in zip (self.wavelengths, self.flux):
			distance = abs(w - wavelength)
			if distance < minDistance:
				minDistance = distance
				result = f
		return result
		
	def getSubsetByWavelength(self, lower, upper):
		newWavelengths = []
		newFlux = []
		for w, f in zip(self.wavelengths, self.flux):
			if (w>=lower) and (w<=upper):
				newWavelengths.append(w)
				newFlux.append(f)
		return (newWavelengths, newFlux)
		
	def parseHeaderInfo(self, headers):
		self.objectName = headers['Object']
		self.telescope = headers['telescope']
		self.ra = headers['ra']
		self.dec = headers['dec']
		self.UTC = headers['UTC']
		self.RJD = headers['RJD']
		self.equinox = headers['Equinox']
		self.Vearth = headers['Vearth']
		self.hourangle = headers['Hour angle']
		self.longitude = headers['Longitude']
		self.latitude = headers['Latitude']
		self.siderial = headers['Sidereal time']
		self.site = headers['Site']
		self.day = headers['day']
		self.month = headers['month']
		self.year = headers['year']
		self.dwell = headers['dwell']
		self.HJD = headers['hjd']
		self.airmass = headers['Airmass']
		self.galLatitude = headers['Gal latitude']
		self.galLongitude = headers['Gal longitude']
		self.extractPosition = headers['Extract position']
		
		return self.objectName

	def appendDataAtNewWavelengths(self, wavelengths, flux):
		currentWavelengthRange = self.wavelengthRange
		newWavelengthsRange = (min(wavelengths), max(wavelengths))
		print "Current wavelength range:", currentWavelengthRange
		print "New data wavelength range:", newWavelengthsRange
		if currentWavelengthRange[1] > newWavelengthsRange[0]:
			print "There *IS* an overlap in wavelengths"
			print "First data loaded takes preference"
			startWavelength = currentWavelengthRange[1]
			print "Starting to add new data from above: ", startWavelength
			for w, f in zip(wavelengths, flux):
				if w>startWavelength:
					self.wavelengths.append(w)
					self.flux.append(f)
					
		self.length = len(self.wavelengths)
		self.wavelengthRange = (min(wavelengths), max(wavelengths))
		
		return self.length
		
	def appendNewData(self, newSpectrum):
		print "HJD of existing spectrum: %7.7f\nHJD of spectrum to be added: %7.7f"%(self.HJD, newSpectrum.HJD)
		timeDifferenceSeconds = 86400. * (self.HJD - newSpectrum.HJD)
		print "Time difference is %f seconds."%timeDifferenceSeconds
		
		currentWavelengthRange = self.wavelengthRange
		wavelengths = newSpectrum.getWavelengths()
		flux = newSpectrum.getFlux()
		newWavelengthsRange = (min(wavelengths), max(wavelengths))
		print "Current wavelength range:", currentWavelengthRange
		print "New data wavelength range:", newWavelengthsRange
		
		# Append data at the end of the current data
		if currentWavelengthRange[0] < newWavelengthsRange[0]:
			if currentWavelengthRange[1] > newWavelengthsRange[0]:
				print "There *IS* an overlap in wavelengths"
				print "First data loaded takes preference"
				startWavelength = currentWavelengthRange[1]
				print "Starting to add new data from above: ", startWavelength
				for w, f in zip(wavelengths, flux):
					if w>startWavelength:
						self.wavelengths.append(w)
						self.flux.append(f)
						
		# Add new data to the beginning of the old data
		if currentWavelengthRange[0] > newWavelengthsRange[0]:
			if currentWavelengthRange[1] > newWavelengthsRange[0]:
				print "There *IS* an overlap in wavelengths"
				print "First data loaded takes preference"
				endWavelength = currentWavelengthRange[0]
				newFlux = []
				newWavelengths = []
				print "Starting to add new data until: ", endWavelength
				for w, f in zip(wavelengths, flux):
					if w<endWavelength:
						newWavelengths.append(w)
						newFlux.append(f)
				for w, f in zip(self.wavelengths, self.flux):
					newWavelengths.append(w)
					newFlux.append(f)
			self.wavelengths = newWavelengths
			self.flux = newFlux
		
						
		self.length = len(self.wavelengths)
		self.wavelengthRange = (min(wavelengths), max(wavelengths))
		
		return self.wavelengthRange

