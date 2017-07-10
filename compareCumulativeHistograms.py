#!/usr/bin/env python
import sys, os
import numpy, math
import argparse
import astropy.io.fits
import astropy.stats
import loadingSavingUtils, statsUtils
import spectrumClasses, timeClasses, generalUtils
import scipy.optimize
import copy
import scipy.signal as signal

import matplotlib.pyplot

		
class population:
	def __init__(self, name):
		self.periods = []
		self.objectNames = []
		self.studyname = name
		self.plotColour = 'r'
		
	def addPeriod(self, objectName, period):
		self.periods.append(period)
		self.objectNames.append(objectName)
		
	def getLogPeriods(self):
		return [numpy.log10(p) for p in self.periods]
		
	def __str__(self):
		return self.studyname + " : " + str([p for p in self.periods]) 
	
if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Compare the cumulative histograms of two period population samples.')
	parser.add_argument('files', type=str, nargs='+', help= 'File(s) containing the population samples.')
	parser.add_argument('--device', type=str, default='/xs', help = "Dump to the following device. Default is '/xs'.")
	parser.add_argument('--save', type=str, help = "Dump the plot to file. Specify the filename with an extension.")
	parser.add_argument('-p', '--probability', type=str,  help= 'File(s) containing an input population probability.')
	
	
	arg = parser.parse_args()

	populations = []
	for f in arg.files:
		p = population(f)
		p.studyname = f
		inputFile = open(f, 'rt')
		headers = inputFile.readline().strip()
		# print headers
		for line in inputFile:
			fields = line.strip().split("\t")
			try:
				name = str(fields[0])
				period = float(fields[-1])
				p.addPeriod(name, period)
			except (ValueError, IndexError) as e:
				print "No valid period for ", name
		populations.append(p)
	
	for p in populations:
		if p.studyname == 'nebot_periods.tsv':
			fixPeriods = [pdays/24. for pdays in p.periods]
			p.periods = fixPeriods 
			p.plotColour = 'k'
			p.studyname = "Nebot et al (2014)"
			
	if arg.probability is not None:
		print arg.probability 
		inputFile = open(arg.probability, 'rt')
		headers = inputFile.readline().strip()
		bins = []
		probabilities = []
		for line in inputFile:
			fields = line.strip().split('\t')
			print fields
			try:
				bin = float(fields[0])
				probability = float(fields[1])
				bins.append(fields[0])
				probabilities.append(fields[1])
			except (ValueError, IndexError):
				print "Could not parse a line in the file"
				
		probabilityFigure = matplotlib.pyplot.figure()
		matplotlib.pyplot.step(bins, probabilities)
		matplotlib.pyplot.draw()
		matplotlib.pyplot.show(block = False)
			
	figure2 = matplotlib.pyplot.figure()
	
	for p in populations:
		print p
		binwidth=0.2
		data = p.getLogPeriods()
		n, bins, patches = matplotlib.pyplot.hist(p.getLogPeriods(), facecolor=p.plotColour, alpha=0.75, cumulative=False, normed = True, histtype = 'step')
		# n, bins, patches = matplotlib.pyplot.hist(p.getLogPeriods(), bins=numpy.arange(min(data), max(data) + binwidth, binwidth), facecolor=p.plotColour, alpha=0.75, cumulative=True, normed = True, histtype = 'step')
		print n, bins, patches
	
		print "Total in bins:", sum(n)
	matplotlib.pyplot.xlabel('$log_{10}(P_{orb})$ [d]',fontsize=18)
	matplotlib.pyplot.ylabel('N',fontsize=18)
	matplotlib.pyplot.title('Cumulative period distribution',fontsize=18)
	matplotlib.pyplot.grid(True)

	matplotlib.pyplot.show(block = False)
	if arg.save is not None:
		figure1.savefig(arg.save, bbox_inches='tight')

	figure2 = matplotlib.pyplot.figure()
	for p in populations:
		minimum = numpy.min(p.getLogPeriods())
		maximum = numpy.max(p.getLogPeriods())
		sortedLogPeriods = p.getLogPeriods()
		sortedLogPeriods = sorted(sortedLogPeriods)
		xvalues = []
		yvalues = []
		for index, s in enumerate(sortedLogPeriods):
			y = (index+1.) / float(len(sortedLogPeriods))
			print s, y
			xvalues.append(s)
			yvalues.append(y)
	
		matplotlib.pyplot.step(xvalues, yvalues, color=p.plotColour)
	
	
	matplotlib.pyplot.xlabel('$log_{10}(P_{orb})$ [d]',fontsize=16)
	matplotlib.pyplot.ylabel('$N_{>log_{10}(P_{orb})} / N_{total}$',fontsize=16)
	matplotlib.pyplot.title('Cumulative period distribution',fontsize=16)
	# legend = figure2.legend( (line1), 'Nebot',  loc='upper left', shadow=True)
	matplotlib.pyplot.show(block = False)
	
	testPopulation = sorted(populations[0].getLogPeriods())
	modelPopulation = sorted(populations[1].getLogPeriods())
	print "Test population:", testPopulation
	print "Model population:", modelPopulation
	s=0
	N = len(testPopulation)
	for index, y in enumerate(testPopulation):
		i = index+1 
		f = float(sum(1 for x in modelPopulation if float(x) <= y) / float(len(modelPopulation)))
		fn = float(sum(1 for x in modelPopulation if float(x) <= testPopulation[N - i]) / float(len(modelPopulation)))
		#s+= (2*i-1) / N * (numpy.log(f) + numpy.log(1 - fn) )
		s+= ((2*float(i)-1) / float(N)) * (numpy.log(f) + numpy.log(1.0 - fn))
		print i, y, s, f, fn
		
		
	from scipy import stats
	print stats.ks_2samp(testPopulation, modelPopulation)
	print stats.anderson_ksamp( (testPopulation, modelPopulation) )
	
	
	generalUtils.query_yes_no("Continue?")
	
	sys.exit()
	
