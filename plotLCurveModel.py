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


def quad(x, a1, a2, a3):
	y = a1 * x * x + a2 *x + a3
	return y


class object:
	def __init__(self, id):
		self.id = id
		self.MJD = []
		self.mag = []
		self.err = []
		self.ephemeris = None
		self.data = []
		self.hasEphemeris = False
		
	def appendData(self, dataDict):
		# self.MJD.append(data['MJD'])
		# self.mag.append(data['mag'])
		# self.err.append(data['err'])
		self.data.append(dataDict)
	
	def getColumn(self, columnName):
		return [d[columnName] for d in self.data]	
		
	def loadEphemeris(self):
		# Look in the local directory for a file called 'id'-ephem.dat and load it
		filename = self.id + "-ephem.dat"
		if os.path.exists(filename):
			self.ephemeris = timeClasses.ephemerisObject()
			self.ephemeris.loadFromFile(filename)
			self.hasEphemeris = True
			return True
		
		return False
		
	def computePhases(self):
		if not self.hasEphemeris: return
		for d in self.data:
			d['phase'] = self.ephemeris.getPhase(d['HJD'])
			print d['HJD'], d['phase']
			
		
	def setHJDs(self, MJD, HJD):
		keys = [d['MJD'] for d in self.data]
		dates = zip(MJD, HJD)
		for index, d in enumerate(dates):
			self.data[index]['HJD'] = d[1]
		
	def convertFluxMagnitude(self):
	 
		fudge = 10
		for d in self.data:
			d['mag'] = -2.5 * numpy.log10(d['flux'])
			d['err'] = -2.5 * d['flux_err'] / numpy.log(10) / d['flux'] /fudge
			# print d['mag'], d['err']
			
	def filterData(self, columnName, limit1, limit2):
		newData = []
		if limit1>limit2:
			temp = limit1
			limit1 = limit2
			limit2 = temp
		for d in self.data:
			value = d[columnName]
			if value>limit1 and value<limit2:
				newData.append(d)
			else:
				print "Filtering out a point: %s %f not between %f and %f"%(columnName, value, limit1, limit2)
		self.data = newData 
	
	def writeData(self, filename):
		outputFile = open(filename, 'wt')
		for d in self.data:
			date = d['HJD']
			phase = d['phase']
			flux = 10**(d['mag'] / -2.5)
			flux_error = flux * numpy.log(10) * d['err']
			print date, phase, d['mag'], flux, flux_error
			outputFile.write("%f 0 %E %E 1 1\n"%(phase, flux, flux_error))
		
		outputFile.close()
		
	def sortData(self, columnName):
		print "Sorting data by:", columnName
		sortedData = self.data.sort(key=lambda x:x[columnName])
		
	def computeHJDs(self):
		if self.hasEphemeris:
			print o.id, o.ephemeris
			MJD = o.getColumn('MJD')
			correctHelio = timeClasses.heliocentric()
			correctHelio.setTelescope('CSS') 
			correctHelio.setTarget(o.ephemeris.ra, o.ephemeris.dec)
			BMJD = correctHelio.convertMJD(MJD)
			HJD = [b + 2400000.5 for b in BMJD]
			self.setHJDs(MJD, HJD)
		
	
if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Plots lcurve model alongside observed data.')
	parser.add_argument('observed', type=str, help='Observed data.')	
	parser.add_argument('model', type=str, help='Modelled data.')	
	parser.add_argument('name', type=str, help='Object name.')	
	parser.add_argument('--ps', action='store_true', help = "Dump plots to ps files instead of the screen.")
	arg = parser.parse_args()

	print "Astropy version:", astropy.__version__

	objects = []
	observedData = object('observed')
	dataFile = open(arg.observed, 'rt')
	for line in dataFile:
		fields = line.strip().split(' ')
		data = {}
		data['phase'] = float(fields[0])
		data['observedtime'] = float(fields[1])
		data['flux'] = float(fields[2])
		data['flux_err'] = float(fields[3])
		# print data
		observedData.appendData(data)
	dataFile.close()
	objects.append(observedData)
	modelledData = object('modelled')
	dataFile = open(arg.model, 'rt')
	for line in dataFile:
		fields = line.strip().split(' ')
		data = {}
		data['phase'] = float(fields[0])
		data['observedtime'] = float(fields[1])
		data['flux'] = float(fields[2])
		data['flux_err'] = float(fields[3])
		print data
		modelledData.appendData(data)
	dataFile.close()
	
	print "%d targets loaded"%len(objects)

	observedData.convertFluxMagnitude()
	modelledData.convertFluxMagnitude()
	
	if arg.ps: device = arg.name + "_lcurve.ps/ps"
	else: device = "/xs"
	PGPlotWindow = ppgplot.pgopen(device) 
	ppgplot.pgask(True)
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
	ppgplot.pgslct(PGPlotWindow)   
	ppgplot.pgsci(1)
	ppgplot.pgask(False)
	phases = observedData.getColumn('phase')
	flux = observedData.getColumn('flux')
	flux_err = observedData.getColumn('flux_err')
	mag = observedData.getColumn('mag')
	err = observedData.getColumn('err')
	maxFlux = 0
	for f, fe in zip(flux, flux_err):
		if f + fe > maxFlux: maxFlux = f + fe
	
	
	meanError = numpy.mean(flux_err)
	# Duplicate data out to phase 2.0
	extendedPhases = copy.deepcopy(phases)
	for p in phases:
		extendedPhases.append(p + 1.0)
	phases = extendedPhases
	flux.extend(flux)
	flux_err.extend(flux_err)
	mag.extend(mag)
	err.extend(err)
			
	"""
	ppgplot.pgsch(1.6)
	ppgplot.pgenv(0, 2, numpy.max(mag), numpy.min(mag), 0, 0)
	ppgplot.pglab("Phase", "PTF magnitude", "%s"%(arg.name))
	ppgplot.pgsch(1.0)
	ppgplot.pgpt(phases, mag)
	ppgplot.pgerrb(2, phases, mag, err, 0)
	ppgplot.pgerrb(4, phases, mag, err, 0)
	
	ppgplot.pgsci(2)
	model = modelledData.getColumn('mag')
	model.extend(model)
	ppgplot.pgsls(2)
	ppgplot.pgslw(7)
	ppgplot.pgline(phases, model)
	"""
	flux_err = [f/10 for f in flux_err]
	ppgplot.pgsch(1.6)
	ppgplot.pgenv(0, 2, 0, maxFlux, 0, 0)
	ppgplot.pglab("Phase", "PTF flux", "%s"%(arg.name))
	ppgplot.pgsch(1.0)
	ppgplot.pgpt(phases, flux)
	ppgplot.pgerrb(2, phases, flux, flux_err, 0)
	ppgplot.pgerrb(4, phases, flux, flux_err, 0)
	
	ppgplot.pgsci(2)
	ppgplot.pgsls(2)
	ppgplot.pgslw(7)
	modelPhases = modelledData.getColumn('phase')
	temp = copy.deepcopy(modelPhases)
	for p in modelPhases: temp.append(p + 1.0)
	modelPhases = temp
	model = modelledData.getColumn('flux')
	model.extend(model)
	ppgplot.pgline(modelPhases, model) 	
	ppgplot.pgsci(1)
	ppgplot.pgsls(1)
	ppgplot.pgslw(1)

	# Plot the inset subplot showing the eclipse
	subFlux = []
	subFluxErr = []
	subPhases = []
	subModelPhases = []
	subModel = []
	startPhase = 0.96
	endPhase = 1.04
	print "maxFlux:", maxFlux
	for index, phase in enumerate(phases):
		if phase > startPhase and phase < endPhase:
			subFlux.append(flux[index])
			subPhases.append(phase)
			subFluxErr.append(flux_err[index])
			
	for index, phase in enumerate(modelPhases):
		if phase > startPhase and phase < endPhase:
			subModel.append(model[index])
			subModelPhases.append(phase)
	
	print ppgplot.pgqvp(0)
	(x1, x2, y1, y2) = ppgplot.pgqvp(0)
	print ppgplot.pgqwin(0)
	xlower = (x1 + x2) / 2 + 0.05
	xupper = x2 - 0.05
	yupper = (y1 + y2) / 2 - 0.05
	ylower = 0.20
	ppgplot.pgsvp(xlower, xupper, ylower, yupper)
	ppgplot.pgswin(startPhase, endPhase, 0, numpy.max(subFlux))
	print ppgplot.pgqvp(0)
	print ppgplot.pgqwin(0)
	# ppgplot.pgsubp(4, 3)
	# ppgplot.pgpanl(2, 2)
	# ppgplot.pgenv(startPhase, endPhase, 0, numpy.max(flux), 0, -1)
	# ppgplot.pglab("Phase", "PTF flux", "%s"%(arg.name))
	ppgplot.pgbox("GBC", 0, 0, "GBC", 0, 0)
	ppgplot.pgpt(subPhases, subFlux)
	ppgplot.pgerrb(2, subPhases, subFlux, subFluxErr, 0)
	ppgplot.pgerrb(4, subPhases, subFlux, subFluxErr, 0)
	
	ppgplot.pgline(subModelPhases, subModel) 	
	
	ppgplot.pgclos()	
	

	sys.exit()
	

