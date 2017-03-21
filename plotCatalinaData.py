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
		
	def setHJDs(self, MJD, HJD):
		keys = [d['MJD'] for d in self.data]
		print keys
		dates = zip(MJD, HJD)
		for index, d in enumerate(dates):
			self.data[index]['HJD'] = d[1]
			
		
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

	parser = argparse.ArgumentParser(description='Loads data that has been downloaded from the Catalina archive. Plots them using MatPlotLib')
	parser.add_argument('inputFiles', type=str, nargs='+', help='Downloaded Catalina CSV file.')
	parser.add_argument('-e', type=str, help='Optional ephemeris file')
	parser.add_argument('--device', type=str, default="/xs", help="Output device. Default is '/xs'")
	parser.add_argument('--ps', action='store_true', help = "Dump plots to ps files instead of the screen.")
	 
	arg = parser.parse_args()
	print arg
	
	print "Astropy version:", astropy.__version__
	
	if arg.e!=None:
		# Load the ephemeris file
		hasEphemeris = True
		ephemeris = timeClasses.ephemerisObject()
		ephemeris.loadFromFile(arg.e)
		print ephemeris
	else:
		hasEphemeris = False

	
	columns = [ {'name': 'ID',        'type':'str'}, 
				{'name': 'CRTSID',    'type':'str'},
                {'name': 'mag',       'type':'float'},
                {'name': 'err',       'type':'float'},
                {'name': 'ra',        'type':'float'},
                {'name': 'dec',       'type':'float'},
                {'name': 'MJD',       'type':'float'},
                {'name': 'blend',     'type':'int'} ]
	columnsLongForm = [ {'name': 'ID',        	'type':'str'}, 
						{'name': 'MasterFrame', 'type':'str'},
                		{'name': 'CRTSID',      'type':'float'},
                		{'name': 'mag',       	'type':'float'},
                		{'name': 'err',      	'type':'float'},
                		{'name': 'ra',        	'type':'float'},
                		{'name': 'dec',       	'type':'float'},
                		{'name': 'FWHM',       	'type':'float'},
                		{'name': 'var',       	'type':'float'},
                		{'name': 'FrameID',     'type':'float'},
                		{'name': 'MJD',       	'type':'float'},
                		{'name': 'airmass',     'type':'float'},
                		{'name': 'exposure',    'type':'float'},
                		{'name': 'X',       	'type':'float'},
                		{'name': 'Y',       	'type':'float'},
                		{'name': 'Flux',       	'type':'float'},
                		{'name': 'Area',       	'type':'float'},
                		{'name': 'Flags',       'type':'int'},
                		{'name': 'Theta',       'type':'float'},
                		{'name': 'Elong',       'type':'float'},
                		{'name': 'NuMax',       'type':'float'},
                		{'name': 'blend',     	'type':'int'} ]
				
	columnNames = ['ID', 'CRTSID', 'mag', 'err', 'ra', 'dec', 'MJD', 'blend' ]
	data = []
	for fileIndex, f in enumerate(arg.inputFiles):
		catalinaFile = open(f, 'rt')
		headings = catalinaFile.readline()
		if len(headings.split(',')) == 8:
			longForm = False
			print "Short form of Catalina data. Available columns are:", [c['name'] for c in columns]
		elif len(headings.split(',')) == 22:
			longForm = True
			columns = columnsLongForm
			print "Long form Catalina data. Available columns are:", [c['name'] for c in columns]
		else:
			print "Something is wrong with the input. CRTS data should have 8 cols (short form) or 22 cols (long form)."
			sys.exit()
		for line in catalinaFile:
			fields = line.split(',')
			d = {}
			for index, column in enumerate(columns):
				value = fields[index].strip(' \t\n\r')
				if column['type']=='str': value = str(value)
				if column['type']=='float': value = float(value)
				if column['type']=='int': value = int(value)
				d[column['name']] = value
			print d
			data.append(d)
			

	
	# Separate different objects
	objects = []
	ids = []
	for d in data:
		id = d['ID']
		if id not in ids:
			ids.append(id)
	print ids
	
	for id in ids:
		o = object(id)
		for d in data:
			if d['ID'] == o.id:
				o.appendData(d)
		objects.append(o)
	
	print "%d targets loaded"%len(objects)

	if arg.ps: device = "lightcurves.ps/ps"
	else: device = "/xs"
	PGPlotWindow = ppgplot.pgopen(device) 
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
	ppgplot.pgslct(PGPlotWindow)   
	ppgplot.pgsci(1)
	ppgplot.pgask(False)
	for index, o in enumerate(objects):
		MJD = o.getColumn('MJD')
		mag = o.getColumn('mag')
		err = o.getColumn('err')
		startDate = numpy.min(MJD)
		endDate = numpy.max(MJD)
		magMax = numpy.max(mag) + err[numpy.argmax(mag)]
		magMin = numpy.min(mag) - err[numpy.argmin(mag)]
		meanError = numpy.mean(err)
		print "%s Start date: %f, End date: %f"%(o.id, startDate, endDate)
		ppgplot.pgenv(startDate, endDate, magMax + meanError*2, magMin - meanError*2, 0, 0)
		ppgplot.pgpt(MJD, mag)
		ppgplot.pgerrb(2, MJD, mag, err, 0)
		ppgplot.pgerrb(4, MJD, mag, err, 0)
		ppgplot.pglab("MJD", "CRTS mag", "%s [%d]"%(o.id, len(MJD)))
	
	ppgplot.pgclos()	
	
	
	# Compute HJDs for the observations
	for o in objects:
		hasEphemeris = o.loadEphemeris()
		if hasEphemeris: o.computeHJDs()
	
			

	
	##########################################################################################################################
	# Periodograms 
	##########################################################################################################################
	plo = 0.01
	phi = 1.00
	if arg.ps: device = "pgrams.ps/ps"
	else: device = "/xs"
	pgramPGPlotWindow = ppgplot.pgopen(device)  
	pgPlotTransform = [0, 1, 0, 0, 0, 1]	
	ppgplot.pgslct(pgramPGPlotWindow)
	for o in objects:
		if o.hasEphemeris:
			HJD = o.getColumn('HJD')
			mag = o.getColumn('mag')
			err = o.getColumn('err')
			x = numpy.array(HJD)
			y = numpy.array(mag)
			# Subtract the mean from the y-data
			y_mean = numpy.mean(y)
			y = y - y_mean
			periods = numpy.linspace(plo, phi, 1000)
			ang_freqs = 2 * numpy.pi / periods
			power = signal.lombscargle(x, y, ang_freqs)
			# normalize the power
			N = len(x)
			power *= 2 / (N * y.std() ** 2)
	
			ppgplot.pgenv(min(periods), max(periods), 0, max(power), 0, 0)
			ppgplot.pgline(periods, power)
			ppgplot.pglab("Period (d)", "Amplitude", "Lomb-Scargle: " + o.id)
			bestPeriod = periods[numpy.argmax(power)]
			lc = ppgplot.pgqci()
			ls = ppgplot.pgqls()
			ppgplot.pgsci(3)
			ppgplot.pgsls(2)
			ppgplot.pgline([bestPeriod, bestPeriod], [0, max(power)])
			ppgplot.pgsci(lc)
			ppgplot.pgsls(ls)
			print "Best period: %f days or %f hours"%(bestPeriod, bestPeriod * 24.)

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
	ppgplot.pgask(True)
	
	for o in objects:
		if o.hasEphemeris:
			HJD = o.getColumn('HJD')
			mag = o.getColumn('mag')
			err = o.getColumn('err')
			phases = [o.ephemeris.getPhase(h) for h in HJD]
			offsetPhases = []
			for p in phases:
				if p<0.5: offsetPhases.append(p + 1.0)
				else: offsetPhases.append(p)
			phases = offsetPhases
			# print phases
			magMax = numpy.max(mag) + err[numpy.argmax(mag)]
			magMin = numpy.min(mag) - err[numpy.argmin(mag)]
			meanError = numpy.mean(err)
			ppgplot.pgenv(0.5 ,1.5 , magMax + meanError*2, magMin - meanError*2, 0, 0)
			ppgplot.pgpt(phases, mag)
			ppgplot.pgerrb(2, phases, mag, err, 0)
			ppgplot.pgerrb(4, phases, mag, err, 0)
			ppgplot.pglab("Phase", "CRTS mag", "Phase plot: %s [%d]"%(o.id, len(phases)) )

			extraColumn = 'FWHM'
			yHeight = (magMax + meanError*2) - (magMin - meanError*2)
			additionalData = o.getColumn(extraColumn)
			FWHMrange = numpy.max(additionalData) - numpy.min(additionalData)
			ppgplot.pgpt(phases, additionalData)
			print additionalData
			print yHeight, FWHMrange
			scaledData = [yHeight/FWHMrange * (fwhm - numpy.min(additionalData)) + numpy.min(mag) for fwhm in additionalData]
			print scaledData
			ppgplot.pgsci(2)
			ppgplot.pgpt(phases, scaledData)
			print "%s min: %f max: %f"%(extraColumn, numpy.min(additionalData), numpy.max(additionalData))
	ppgplot.pgclos()
	sys.exit()
	

