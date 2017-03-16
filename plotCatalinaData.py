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
		
	def appendData(self, data):
		self.MJD.append(data['MJD'])
		self.mag.append(data['mag'])
		self.err.append(data['err'])
		
	def loadEphemeris(self):
		# Look in the local directory for a file called 'id'-ephem.dat and load it
		filename = self.id + "-ephem.dat"
		if os.path.exists(filename):
			self.ephemeris = timeClasses.ephemerisObject()
			self.ephemeris.loadFromFile(filename)
			return True
		
		return False
		
	def setHJD(self, HJD):
		self.HJD = HJD
		
	
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

	
	columnNames = ['ID', 'CRTSID', 'mag', 'err', 'ra', 'dec', 'MJD', 'blend' ]
	data = []
	for fileIndex, f in enumerate(arg.inputFiles):
		catalinaFile = open(f, 'rt')
		headings = catalinaFile.readline()
		for line in catalinaFile:
			fields = line.split(',')
			d = {}
			d['ID'] = fields[0].strip(' \t\n\r')
			d['CRTSID'] = fields[1].strip(' \t\n\r')
			d['mag'] = float(fields[2].strip(' \t\n\r'))
			d['err'] = float(fields[3].strip(' \t\n\r'))
			d['ra'] = float(fields[4].strip(' \t\n\r'))
			d['dec'] = float(fields[5].strip(' \t\n\r'))
			d['MJD'] = float(fields[6].strip(' \t\n\r'))
			d['blend'] = fields[7].strip(' \t\n\r')
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
		startDate = numpy.min(o.MJD)
		endDate = numpy.max(o.MJD)
		magMax = numpy.max(o.mag) + o.err[numpy.argmax(o.mag)]
		magMin = numpy.min(o.mag) - o.err[numpy.argmin(o.mag)]
		meanError = numpy.mean(o.err)
		print "%s Start date: %f, End date: %f"%(o.id, startDate, endDate)
		ppgplot.pgenv(startDate, endDate, magMax + meanError*2, magMin - meanError*2, 0, 0)
		ppgplot.pgpt(o.MJD, o.mag)
		ppgplot.pgerrb(2, o.MJD, o.mag, o.err, 0)
		ppgplot.pgerrb(4, o.MJD, o.mag, o.err, 0)
		ppgplot.pglab("MJD", "CRTS mag", "%s [%d]"%(o.id, len(o.MJD)))
	
	ppgplot.pgclos()	
	

	# Compute HJDs for the observations
	for o in objects:
		hasEphemeris = o.loadEphemeris()
		if hasEphemeris:
			print o.id, o.ephemeris
			correctHelio = timeClasses.heliocentric()
			correctHelio.setTelescope('CSS') 
			correctHelio.setTarget(o.ephemeris.ra, o.ephemeris.dec)
			BMJD = correctHelio.convertMJD(o.MJD)
			HJD = [b + 2400000.5 for b in BMJD]
			o.setHJD(HJD)
			

	
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
		hasEphemeris = o.loadEphemeris()
		if hasEphemeris:
			x = numpy.array(o.HJD)
			y = numpy.array(o.mag)
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
		hasEphemeris = o.loadEphemeris()
		if hasEphemeris:
			"""print o.id, o.ephemeris
			JD = [b + 2400000.5 for b in o.MJD]
			correctHelio = timeClasses.heliocentric()
			correctHelio.setTelescope('CSS') 
			correctHelio.setTarget(o.ephemeris.ra, o.ephemeris.dec)
			BMJD = correctHelio.convertMJD(o.MJD)
			HJD = [b + 2400000.5 for b in BMJD]
			"""
			phases = [o.ephemeris.getPhase(h) for h in o.HJD]
			offsetPhases = []
			for p in phases:
				if p<0.5: offsetPhases.append(p + 1.0)
				else: offsetPhases.append(p)
			phases = offsetPhases
			# print phases
			magMax = numpy.max(o.mag) + o.err[numpy.argmax(o.mag)]
			magMin = numpy.min(o.mag) - o.err[numpy.argmin(o.mag)]
			meanError = numpy.mean(o.err)
			ppgplot.pgenv(0.5 ,1.5 , magMax + meanError*2, magMin - meanError*2, 0, 0)
			ppgplot.pgpt(phases, o.mag)
			ppgplot.pgerrb(2, phases, o.mag, o.err, 0)
			ppgplot.pgerrb(4, phases, o.mag, o.err, 0)
			ppgplot.pglab("Phase", "CRTS mag", "Phase plot: %s [%d]"%(o.id, len(phases)) )

	ppgplot.pgclos()
	sys.exit()
	

