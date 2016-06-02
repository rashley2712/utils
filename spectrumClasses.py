import numpy
import json
import scipy.interpolate
import scipy.integrate

class spectrumObject:
	def __init__(self):
		self.wavelengths = []
		self.flux = []
		self.length = 0 
		self.wavelengthRange = (0, 0)
		self.name = 'unknown'
		self.loadedFromFilename = 'unknown'
		self.wavelengthUnits = 'unknown'
		self.fluxUnits = 'unknown'
		
	def getProperty(self, property):
		try:
			data = getattr(self, property)
			return data
		except AttributeError:
			return None
		  
		
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
		
	def sortData(self):
		wavelengths = self.wavelengths
		fluxes = self.flux
		list1, list2 = zip(*sorted(zip(wavelengths, fluxes)))
		self.wavelengths = list1
		self.flux = list2

	def resample(self, sampleWavelengths):
		startWavelength = min(sampleWavelengths)
		endWavelength = max(sampleWavelengths)
		self.trimWavelengthRange(startWavelength, endWavelength)
		print "num points:", len(self.wavelengths)
		spline = scipy.interpolate.splrep(self.wavelengths, self.flux, s=0)
		sampleFlux = scipy.interpolate.splev(sampleWavelengths, spline, der=0)
		self.wavelengths = sampleWavelengths
		self.flux = sampleFlux
		return sampleFlux

	def writeCSV(self, filename):
		outputfile = open(filename, 'w')
		outputfile.write("wavelength, flux\n")
		for w, f in zip(self.wavelengths, self.flux):
			outputfile.write("%f, %f\n"%(w, f))
		outputfile.close()
		
	def snipWavelengthRange(self, lower, upper):
		""" Removes a section from the spectrum """
		newWavelengths = []
		newFlux = []
		if lower>=upper: return self.length
			
		for w, f in zip(self.wavelengths, self.flux):
			if (w<lower) or (w>upper):
				newWavelengths.append(w)
				newFlux.append(f)
		self.length = len(newWavelengths)
		self.wavelengths = newWavelengths
		self.flux = newFlux
		self.wavelengthRange = (min(self.wavelengths), max(self.wavelengths))
			
		return self.length		
	
	def integrate(self, wavelengthrange = (-1, -1)):
		""" Integrates under the spectrum between two wavelength limits. Defaults to all of the spectrum """
		if wavelengthrange[0]!=-1:
			wavelengths, fluxes = self.getSubsetByWavelength(wavelengthrange[0], wavelengthrange[1])
		else: 
			wavelengths, fluxes = (self.wavelengths, self.flux)
		
		total = scipy.integrate.simps(fluxes, wavelengths)
		return total

	def divide(self, constant):
		""" Divides the spectrum by a constant value """
		newFlux = []
		for w, f in zip(self.wavelengths, self.flux):
			newFlux.append(f / constant)
		self.flux = newFlux
		return 

	def subtractSpectrum(self, subtractSpectrum):
		if len(self.wavelengths)!=len(subtractSpectrum.wavelengths):
			print "Can't subtract spectra of different lengths."
			return
		newFlux = []
		for (aw, af, bw, bf) in zip(self.wavelengths, self.flux, subtractSpectrum.wavelengths, subtractSpectrum.flux):
			f = af - bf
			#print af, bf, f
			newFlux.append(f)
		self.flux = newFlux
		return
		
	def trimWavelengthRange(self, lower, upper):
		""" Trims out the lower and upper portions of the spectrum """ 
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
		self.wavelengthRange = (min(self.wavelengths), max(self.wavelengths))
			
		return self.length		
		
	def getWavelengths(self):
		return self.wavelengths
		
	def convertFluxes(self):
		print "Current fluxUnits are:", self.fluxUnits
		
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

	def writeToJSON(self, filename):
		object = {}
		for key in self.__dict__.keys():
			data = getattr(self, key)
			# print key, type(data)
			if type(data)==numpy.float32:
				data = float(data)
			if type(data)==numpy.ndarray:
				data = numpy.array(data).tolist()
			if type(data)==list:
				data = numpy.array(data).tolist()
			object[key] = data
			
		outputfile = open(filename, 'w')
		json.dump(object, outputfile)
		outputfile.close()
	
	def loadFromJSON(self, filename):
		inputfile = open(filename, "r")
		jsonObject = json.load(inputfile)
		for key in jsonObject.keys():
			keyString = str(key)
			value = jsonObject[key]
			if type(value) is unicode: 
				value = str(value)
			if type(value) is list:
				value = numpy.array(value)
			setattr(self, key, value)
		inputfile.close()
		self.loadedFromFilename = filename
		
		
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

