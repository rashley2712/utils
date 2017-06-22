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
import matplotlib.pyplot
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
		
		
class population:
	def __init__(self, name):
		self.periods = []
		self.objectNames = []
		self.studyname = name
		self.plotColour = 'g'
		
	def addPeriod(self, objectName, period):
		self.periods.append(period)
		self.objectNames.append(objectName)
		
	def getLogPeriods(self):
		return [numpy.log10(p) for p in self.periods]
		
	
if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Loads CSV file and plots histogram.')
	parser.add_argument('filename', type=str, nargs='+', help='Filename of CSV file..')	
	parser.add_argument('--save', type=str, help = "Dump the plot to file. Specify the filename with an extension.")
	
	
	arg = parser.parse_args()
	
	files = arg.filename
	names = []
	periods = []
	study = []
	populations = []
	for f in files:
		p = population(f)
		inputFile = open(f, 'rt')
		headers = inputFile.readline().strip()
		print headers
		for line in inputFile:
			fields = line.strip().split(",")
			try:
				name = str(fields[0])
				period = float(fields[-1])
				periods.append(period)
				names.append(name)
				study.append(f)
				p.addPeriod(name, period)
			except (ValueError, IndexError) as e:
				print "No valid period for ", name
		populations.append(p)
	print zip(names, periods, study)
	logPeriods = [numpy.log10(p) for p in periods]
	# print zip(names, logPeriods)
	
	populations[1].plotColour = 'r'
	
	for p in populations:
		print p.studyname, p.getLogPeriods()
	
	print "%d targets loaded"%len(names)

	# Draw the histogram of the input data
	figure1 = matplotlib.pyplot.figure(figsize=(7,10))
	for p in populations:
		binwidth=0.5
		data = p.getLogPeriods()
		n, bins, patches = matplotlib.pyplot.hist(p.getLogPeriods(), bins=numpy.arange(min(data), max(data) + binwidth, binwidth), facecolor=p.plotColour, alpha=0.75, cumulative=False)
		print n, bins, patches
	print "Total in bins:", sum(n)
	matplotlib.pyplot.xlabel('$log_{10}(P_{orb})$ [d]',fontsize=18)
	matplotlib.pyplot.ylabel('N',fontsize=18)
	matplotlib.pyplot.title('Observed period distribution',fontsize=18)
	matplotlib.pyplot.grid(True)

	matplotlib.pyplot.show()
	if arg.save is not None:
		figure1.savefig(arg.save, bbox_inches='tight')

	"""if arg.ps: device = "histogram.ps/ps"
	else: device = "/xs"
	PGPlotWindow = ppgplot.pgopen(device) 
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
	ppgplot.pgslct(PGPlotWindow)   
	ppgplot.pgsci(1)
	ppgplot.pgask(False)
	ppgplot.pgsch(1.6)
	ppgplot.pgenv(-1.5, 1.5, 0, 7, 0)
	ppgplot.pglab("log\d10\u(P\dorb\u/days)", "N", "")
	ppgplot.pgsch(1.0)
	# ppgplot.pgsvp(0.1, 0.9, 0.1, 0.9)
	ppgplot.pghist(logPeriods, -1.5, 1.5, 9, 5)
	
	ppgplot.pgclos()	
	"""
	sys.exit()
	
	
