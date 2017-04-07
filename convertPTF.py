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

	parser = argparse.ArgumentParser(description='Converts data from Paul''s PTF extraction into something more useable.')
	parser.add_argument('object', type=str, nargs='+', help='Object name.')	
	parser.add_argument('--ps', action='store_true', help = "Dump plots to ps files instead of the screen.")
	parser.add_argument('-l', '--list', action='store_true', help="Object name is a list of objects.")
	
	
	arg = parser.parse_args()
	print arg
	if arg.list: objectList = True
	else: objectList = False
	
	print "Astropy version:", astropy.__version__

	objectNames = []
	if objectList:
		nameFile = open(arg.object[0], 'rt')
		for l in nameFile:
			objectnames.append(l.strip())
		nameFile.close()
	else:
		objectNames = arg.object
        
	objects = []
	for o in objectNames:
		target = object(o)
		dataFilename = o + '_ptf.dat'
		ptfFile = open(dataFilename, 'rt')
		
		for line in ptfFile:
			line = line.strip()
			if line[0]=='#': 
				print "COMMENT: ", line[2:]
				if 'Reference' in line:
					startDate = float(line.split(':')[-1].strip(' '))
			else:
				fields = line.split(' ')
				data = {}
				data['HJD'] = startDate + float(fields[0])
				data['mag'] = float(fields[1])
				data['err'] = float(fields[2])
				target.appendData(data)
		ptfFile.close()
		objects.append(target)


	

	
	print "%d targets loaded"%len(objects)

	#for o in objects:
	#	o.filterData('mag', 10, 16.5)

	if arg.ps: device = "lightcurves.ps/ps"
	else: device = "/xs"
	PGPlotWindow = ppgplot.pgopen(device) 
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
	ppgplot.pgslct(PGPlotWindow)   
	ppgplot.pgsci(1)
	ppgplot.pgask(False)
	for index, o in enumerate(objects):
		HJD = o.getColumn('HJD')
		mag = o.getColumn('mag')
		err = o.getColumn('err')
		startDate = numpy.min(HJD)
		dates = [ d - startDate for d in HJD]
		endDate = numpy.max(HJD)
		magMax = numpy.max(mag) + err[numpy.argmax(mag)]
		magMin = numpy.min(mag) - err[numpy.argmin(mag)]
		meanError = numpy.mean(err)
		print "%s Start date: %f, End date: %f"%(o.id, startDate, endDate)
		ppgplot.pgenv(0, numpy.max(dates), magMax + meanError*2, magMin - meanError*2, 0, 0)
		ppgplot.pgpt(dates, mag)
		ppgplot.pgerrb(2, dates, mag, err, 0)
		ppgplot.pgerrb(4, dates, mag, err, 0)
		ppgplot.pglab("Days since %f"%startDate, "PTF mag", "%s [%d]"%(o.id, len(HJD)))
	
	ppgplot.pgclos()	
	
	# Load the ephemerides
	for o in objects:
		hasEphemeris = o.loadEphemeris()
		print o.ephemeris
	
	
	# Write the object data to a textfile
	for o in objects:
		o.computePhases()
		o.sortData('phase')
		o.writeData(o.id + "_phases.dat")
		
		
	##########################################################################################################################
	# Periodograms 
	##########################################################################################################################
	
	from astropy.stats import LombScargle

	if arg.ps: device = "periodograms.ps/ps"
	else: device = "/xs"
	pgramPlotWindow = ppgplot.pgopen(device)  
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
	ppgplot.pgslct(pgramPlotWindow)   
	ppgplot.pgask(True)

	for o in objects:
		HJD = o.getColumn('HJD')
		mag = o.getColumn('mag')
		err = o.getColumn('err')
		frequency, power = LombScargle(HJD, mag, err).autopower(nyquist_factor=20)
		bestFrequency = frequency[numpy.argmax(power)]
		period = 1/bestFrequency
		ppgplot.pgenv(numpy.min(frequency) ,numpy.max(frequency) , numpy.min(power), numpy.max(power), 0, 0)
		ppgplot.pglab("Frequency (days\u-1\d)", "Power", "Periodogram: %s period: %f"%(o.id, period) )
		ppgplot.pgline(frequency, power)
		print "Best frequency for %s is %f cycles/day or a period of %f days."%(o.id, bestFrequency, period) 
	ppgplot.pgclos()
	
	##########################################################################################################################
	# Phase Plots 
	##########################################################################################################################
	if arg.ps: device = "phaseplots.ps/ps"
	else: device = "/xs"
	phasePlotWindow = ppgplot.pgopen(device)  
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
	ppgplot.pgslct(phasePlotWindow)   
	ppgplot.pgsci(1)
	ppgplot.pgpap(10, 0.618)
	ppgplot.pgask(True)
	
	for o in objects:
		if o.hasEphemeris:
			HJD = o.getColumn('HJD')
			mag = o.getColumn('mag')
			err = o.getColumn('err')
			phases = [o.ephemeris.getPhase(h) for h in HJD]
			extendPhases = copy.deepcopy(phases)
			for p in phases:
				extendPhases.append(p + 1.0)
			phases = extendPhases
			magMax = numpy.max(mag) + err[numpy.argmax(mag)]
			magMin = numpy.min(mag) - err[numpy.argmin(mag)]
			meanError = numpy.mean(err)
			mag.extend(mag)
			err.extend(err)
			# print phases
			ppgplot.pgenv(0. ,2.0 , magMax + meanError*2, magMin - meanError*2, 0, 0)
			ppgplot.pglab("Phase", "PTF mag", "Phase plot: %s [%d]"%(o.id, len(phases)/2) )
			ppgplot.pgpt(phases, mag)
			ppgplot.pgerrb(2, phases, mag, err, 0)
			ppgplot.pgerrb(4, phases, mag, err, 0)
			
	ppgplot.pgclos()
	sys.exit()
	

