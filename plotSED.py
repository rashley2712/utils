#!/usr/bin/env python
import sys, os
import numpy, math
import argparse
import astropy.io.fits
import astropy.stats
import loadingSavingUtils, statsUtils
import spectrumClasses, timeClasses
import scipy.optimize
import copy
import ppgplot
import scipy.signal as signal

magOffsets = { 'fuv': 0, 'nuv': 0, 'u': 0, 'g': 0, 'r': 0, 'i' : 0, 'z': 0, 'j': 0.89, 'h': 1.37, 'k': 1.84 }


class object:
	def __init__(self, id="unknown"):
		self.objectID = id
		self.magData = []
		self.bandNames = None
		self.wavelength = None
		self.magnitudes = None
		
	def addData(self, band, mag, magerr):
		dataPoint = {}
		dataPoint['band'] = band
		try:
			dataPoint['mag'] = float(mag)
			dataPoint['magerr'] = float(magerr)
			self.magData.append(dataPoint)
		except ValueError:
			print "No magnitude for %s band"%band
			
	def calculateOffsets(self):
		for d in self.magData:
			d['mag']+= magOffsets[d['band']]
			
		
	def calculateFluxes(self):
		for d in self.magData:
			flux = 3631000 * 10**(d['mag']/-2.5) 	
			# fluxerr =  -2 * numpy.log(10) * 3631000 * d['magerr'] / (5 * 10**(2*d['mag']/5))
			fluxerr = flux * numpy.log(10) * d['magerr'] / (-2.5)
  
			d['flux'] = flux
			d['fluxerr'] = fluxerr
			print d['flux'], d['fluxerr'], d['magerr']
			# print d
			
	def addWavelengths(self, lookup):
		for d in self.magData:
			for w in lookup:
				if d['band'] == w['band']: d['wavelength'] = w['wavelength']
	
	def getFluxData(self):
		return [d['wavelength'] for d in self.magData] , [d['flux'] for d in self.magData], [d['fluxerr'] for d in self.magData], [d['band'] for d in self.magData]
			
	def __str__(self):
		retStr = self.objectID + ": "
		for d in self.magData:
			retStr +=str(d)
		return retStr

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Loads CSV file and plots Spectral Energy distribution.')
	parser.add_argument('filename', type=str, nargs='+', help='Filename of CSV file..')	
	parser.add_argument('--device', type=str, default='/xs', help = "Dump plot to the following device.")
	
	
	arg = parser.parse_args()
	
	files = arg.filename
	bandNames = []
	wavelengths = []
	objects = []
	state = "headers"
	for f in files:
		inputFile = open(f, 'rt')
		for l in inputFile:
			line = l.strip()
			if len(line.strip()) == 0: continue
			# print line
			
			if state == "bands":
				bandNames = line.split(',')
				bandNames = [b.strip() for b in bandNames]
				print "Filter bands:", bandNames
				state = "headers"
				
			if state == "wavelength":
				wavelengths = line.split(',')
				wavelengths = [float(w) for w in wavelengths]
				print "Effective wavelengths:", wavelengths
				state = "headers"
				
			if state == "magnitudes":
				objectID = line.split(',')[0].strip()
				newObject = object(objectID)
				newObject.bandNames = bandNames
				newObject.wavelengths = wavelengths
				magnitudes = line.split(',')[1:]
				print magnitudes
				if len(magnitudes)!=2 * len(bandNames):
					print "Wrong number of data points for %s. Fill in missing data with '-'"
					continue
				else:
					for index, b in enumerate(bandNames):
						print index, b, magnitudes[index*2], magnitudes[index*2 + 1]
						newObject.addData(b,  magnitudes[index*2], magnitudes[index*2 + 1])
				objects.append(newObject)
				print newObject
				
			if "bands" in line:
				state = "bands"
			if "effective" in line:
				state = "wavelength"
			if "abmagnitudes" in line:
				state = "magnitudes"
				
		inputFile.close()
		
	if len(wavelengths)!=len(bandNames):
		print "You didn't specify the same number of effective wavelengths and names for the passbands."
		sys.exit()
	wavelengthLookup = []
	for w, b in zip(wavelengths, bandNames):
		wl = {'band': b, 'wavelength': w}
		wavelengthLookup.append(wl)
	print wavelengthLookup
	
		
	print "%d targets loaded"%len(objects)

	for o in objects:
		print "Calculating AB magnitude offsets for %s"%o.objectID
		o.calculateOffsets()
		o.calculateFluxes()
		
		o.addWavelengths(wavelengthLookup)

	PGPlotWindow = ppgplot.pgopen(arg.device) 
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
	ppgplot.pgslct(PGPlotWindow)   
	ppgplot.pgsci(1)
	ppgplot.pgask(True)
	
	for o in objects:
		ppgplot.pgsch(1.6)
		wavelengths, fluxes, fluxerrors, bands = o.getFluxData()
		fluxMax = max(fluxes)
		fluxMin = min(fluxes)
		wavelengthMin = min(wavelengths)
		wavelengthMax = max(wavelengths)
		ppgplot.pgenv(500,25000 , 0, fluxMax*1.2, 0)
		ppgplot.pglab("wavelength [\A]", "f\d\gn\u [mJy]", o.objectID)
		ppgplot.pgsch(1.0)
		ppgplot.pgpt(wavelengths, fluxes)
		ppgplot.pgerrb(2, wavelengths, fluxes, fluxerrors, 0)
		ppgplot.pgerrb(4, wavelengths, fluxes, fluxerrors, 0)
	
	ppgplot.pgclos()	
	
	sys.exit()
	
	
