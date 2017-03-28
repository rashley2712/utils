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

	parser = argparse.ArgumentParser(description='Converts data from Paul''s PTF extraction into something more useable.')
	parser.add_argument('object', type=str, nargs='+', help='Object name.')	
	parser.add_argument('--ps', action='store_true', help = "Dump plots to ps files instead of the screen.")
	
	
	arg = parser.parse_args()
	print arg
	
	print "Astropy version:", astropy.__version__

	objects = []
	
	for o in arg.object:
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

		objects.append(target)

	

	
	print "%d targets loaded"%len(objects)


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
		endDate = numpy.max(HJD)
		magMax = numpy.max(mag) + err[numpy.argmax(mag)]
		magMin = numpy.min(mag) - err[numpy.argmin(mag)]
		meanError = numpy.mean(err)
		print "%s Start date: %f, End date: %f"%(o.id, startDate, endDate)
		ppgplot.pgenv(startDate, endDate, magMax + meanError*2, magMin - meanError*2, 0, 0)
		ppgplot.pgpt(HJD, mag)
		ppgplot.pgerrb(2, HJD, mag, err, 0)
		ppgplot.pgerrb(4, HJD, mag, err, 0)
		ppgplot.pglab("HJD", "PTF mag", "%s [%d]"%(o.id, len(HJD)))
	
	ppgplot.pgclos()	
	
	# Compute HJDs for the observations
	for o in objects:
		hasEphemeris = o.loadEphemeris()
		print o.ephemeris
		
	##########################################################################################################################
	# Phase Plots 
	##########################################################################################################################
	if arg.ps: device = "phaseplots.ps/ps"
	else: device = "/xs"
	phasePlotWindow = ppgplot.pgopen(device)  
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
	ppgplot.pgslct(phasePlotWindow)   
	ppgplot.pgsci(1)
	ppgplot.pgpap(12, 0.618)
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
			mag.extend(mag)
			err.extend(err)
			# print phases
			magMax = numpy.max(mag) + err[numpy.argmax(mag)]
			magMin = numpy.min(mag) - err[numpy.argmin(mag)]
			meanError = numpy.mean(err)
			ppgplot.pgenv(0. ,2.0 , magMax + meanError*2, magMin - meanError*2, 0, 0)
			ppgplot.pglab("Phase", "PTF mag", "Phase plot: %s [%d]"%(o.id, len(phases)/2) )
			ppgplot.pgsch(1.0)
			ppgplot.pgpt(phases, mag)
			ppgplot.pgerrb(2, phases, mag, err, 0)
			ppgplot.pgerrb(4, phases, mag, err, 0)
			
	ppgplot.pgclos()
	sys.exit()
	

