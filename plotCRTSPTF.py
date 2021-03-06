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
import matplotlib.pyplot
from matplotlib import rc



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
		self.PTFData = []

	def appendData(self, dataDict):
		# self.MJD.append(data['MJD'])
		# self.mag.append(data['mag'])
		# self.err.append(data['err'])
		self.data.append(dataDict)

	def appendPTFData(self, data):
		print data['HJD'], data['mag']
		self.PTFData.append(data)		
	
	def getColumn(self, columnName):
		return [d[columnName] for d in self.data]	
	
	def getPTFColumn(self, columnName):
		return [d[columnName] for d in self.PTFData]	
	
	def meanError(self):
		CRTSmean = numpy.mean([d['err'] for d in self.data])
		PTFmean = numpy.mean([d['err'] for d in self.PTFData])
		return (CRTSmean + PTFmean)/2.
	
	def magMax(self):
		CRTSmax = max([d['mag'] for d in self.data])
		PTFmax = max([d['mag'] for d in self.PTFData])
		print PTFmax, CRTSmax
		if CRTSmax>PTFmax: return CRTSmax
		else: return PTFmax
	
	def magMin(self):
		CRTSmin = min([d['mag'] for d in self.data])
		PTFmin = min([d['mag'] for d in self.PTFData])
		if CRTSmin<PTFmin: return CRTSmin
		else: return PTFmin	

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

	parser = argparse.ArgumentParser(description='Loads data that has been downloaded from the Catalina archive. Plots them using MatPlotLib')
	parser.add_argument('inputFiles', type=str, nargs='+', help='Downloaded Catalina CSV file.')
	parser.add_argument('-e', type=str, help='Optional ephemeris file')
	parser.add_argument('--device', type=str, default="/xs", help="Output device. Default is '/xs'")
	parser.add_argument('--ps', action='store_true', help = "Dump plots to ps files instead of the screen.")
	parser.add_argument('--column', type=str, help = "Plot an extra column of Catalina data.")
	 
	arg = parser.parse_args()
	print arg
	
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
	
			
	# Load the PTF data
	for o in objects:
		targetName = o.id
		print "Loading PTF data for ", targetName
		dataFilename = targetName + '_ptf.dat'
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
				o.appendPTFData(data)
		ptfFile.close()
	
	
	##########################################################################################################################
	# Phase Plots 
	##########################################################################################################################
	if arg.ps: device = "phaseplots.ps/ps"
	else: device = "/xs"
	phasePlotWindow = ppgplot.pgopen(device)  
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
	ppgplot.pgslct(phasePlotWindow)   
	ppgplot.pgsci(1)
	# ppgplot.pgpap(3, 0.618)
	ppgplot.pgask(True)
	
	figSize = 10
	matplotlib.pyplot.figure(figsize=(figSize, figSize / 1.618))
	rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	## for Palatino and other serif fonts use:
	rc('font',**{'family':'serif','serif':['Palatino']})
	rc('text', usetex=True)


	labelSize = 20
	tickSize = 18
	
	
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
			ppgplot.pgsch(1.6)
			meanError = o.meanError()
			magMin = o.magMin()
			magMax = o.magMax()
			print "Max %f, min %f, error %f"%(magMax, magMin, meanError)
			ppgplot.pgenv(0. ,2.0 , magMax + meanError*2, magMin - meanError*2, 0, 0)
			# ppgplot.pglab("Phase", "CRTS mag", "Phase plot: %s [%d]"%(o.id, len(phases)/2) )
			ppgplot.pglab("Phase", "CRTS/PTF mag", "WD%s"%o.id)
			ppgplot.pgsch(1.0)
			ppgplot.pgpt(phases, mag)
			ppgplot.pgerrb(2, phases, mag, err, 0)
			ppgplot.pgerrb(4, phases, mag, err, 0)
			
			matplotlib.pyplot.xlabel("Phase", size = labelSize)
			matplotlib.pyplot.ylabel('CRTS/PTF magnitude', size = labelSize)
			matplotlib.pyplot.errorbar(phases, mag, color='k', yerr=err, fmt = '.', ecolor='0.75', capsize=0, label="CRTS")
			
			
			# Now plot the PTF data
			HJD = o.getPTFColumn('HJD')
			mag = o.getPTFColumn('mag')
			err = o.getPTFColumn('err')
			phases = [o.ephemeris.getPhase(h) for h in HJD]
			print HJD, phases
			extendPhases = copy.deepcopy(phases)
			for p in phases:
				extendPhases.append(p + 1.0)
			phases = extendPhases
			mag.extend(mag)
			err.extend(err)
			
			magMax = numpy.max(mag) + err[numpy.argmax(mag)]
			magMin = numpy.min(mag) - err[numpy.argmin(mag)]
			meanError = numpy.mean(err)
			# ppgplot.pglab("Phase", "CRTS mag", "Phase plot: %s [%d]"%(o.id, len(phases)/2) )
			# ppgplot.pglab("Phase", "PTF mag", "WD %s"%o.id)
			ppgplot.pgsch(1.0)
			ppgplot.pgsci(2)
			ppgplot.pgpt(phases, mag, 4)
			ppgplot.pgerrb(2, phases, mag, err, 0)
			ppgplot.pgerrb(4, phases, mag, err, 0)
			
			matplotlib.pyplot.errorbar(phases, mag, color='k', yerr=err, fmt = '.', ecolor='0.75', capsize=0, label="PTF", marker='s')
			
			
			axes = matplotlib.pyplot.gca()
			for label in (axes.get_xticklabels() + axes.get_yticklabels()):
			#label.set_fontname('Arial')
				label.set_fontsize(tickSize)
			
			print("Phase folded according to the following ephemeris: %s"%o.ephemeris)
			matplotlib.pyplot.legend()
			matplotlib.pyplot.gca().invert_yaxis()
			fig = matplotlib.pyplot.gcf()
			matplotlib.pyplot.show()
			fig.savefig('lightcurve.pdf',dpi=100, format='pdf')
			print("Phase folded according to the following ephemeris: %s"%o.ephemeris)
	ppgplot.pgclos()
	sys.exit()
	

