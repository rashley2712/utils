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
import ppgplot
import scipy.signal as signal

import matplotlib
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

def rv(K2, period, phase, date):
	return K2 * numpy.sin(2*numpy.pi/period*date + 2*numpy.pi*phase)	

def massFunction(m1, m2, period, i):
	Msol =  1.989e30  		# Solar mass (kg)
	G = 6.678e-11     		# Gravitational constant (m^3/kg/s^2)
	m1 = m1 * Msol
	m2 = m2 * Msol
	p = period * 86400.  	# Period in seconds
	return numpy.power( ((m1 * numpy.sin(i))**3 / (m1 + m2)**2) * 2. * numpy.pi * G / p  , 1./3.) /1000.
	
class sampleObservations:
	def __init__(self):
		self.targets = []
	
	def addTarget(self, name):
		target = {}
		target['name'] = name
		target['HJD'] = []
		self.targets.append(target)
		
	def addDataToTarget(self, targetName, HJD):
		for t in self.targets:
			if t['name'] == targetName: 
				t['HJD'].append(HJD)
				
	def dumpTargets(self):
		for t in self.targets:
			print t['name'], len(t['HJD']), 'observations'
			print t['HJD']	
			
	def getExtremes(self):
		for index, t in enumerate(self.targets):
			earliest = numpy.min(t['HJD'])
			latest = numpy.max(t['HJD'])
			if index==0:
				minTime = earliest
				maxTime = latest
			if earliest<minTime:
				minTime = earliest
			if latest>maxTime:
				maxTime = latest
		return minTime, maxTime
			
	def getAverageObservationsPerTarget(self):
		numObs = [len(t['HJD']) for t in self.targets]
		return numpy.mean(numObs)
		
		
	def getNumber(self):
		return len(self.targets)	
		
	def getObs(self, index):
		return self.targets[index]
		
	def getObsByName(self, name):
		for t in self.targets:
			if t['name'] == name: return t
	

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Attempt to calculate the probability of detecting RV variability based on the times of observations. ')
	#parser.add_argument('filename', type=str, nargs='+', help='Filename of CSV file..')	
	parser.add_argument('-o', '--observations', type=str, help='Observation times')
	parser.add_argument('--fake', type=int, help='Fake [n] observations evenly spaced over baseline of all observations. Specify [n].')
	parser.add_argument('-n', '--simulations', type=int, default=100, help='Number of random periods to generate. Default is 100.')
	parser.add_argument('--baseline', type=float, help='Baseline in years over which to fake the observations.') 
	parser.add_argument('-s', '--sigma', default=15.0,  type=float, help='Measurement error in km/s of a typical spectral line fit.') 
	parser.add_argument('-p', '--logprange', nargs = 2, type = float, default = [-2, 2], help='Period range in log P to run over. Specify lower and upper bounds. Default is -2 2')
	parser.add_argument('--steps', default=400., type=float, help='Number of steps over the log P range. Default is 400.')
	parser.add_argument('--dump', help='Dump plot to a .pdf file. Specify filename (without extension).')
	arg = parser.parse_args()
	# print arg
	
	if arg.observations is None and arg.fake is None:
		print "Either specify a file of observations times or ask to fake the data."
		sys.exit()
	
	logpStart = arg.logprange[0]
	logpStop = arg.logprange[1]
	
	generalUtils.setMatplotlibDefaults()
	
	# Load the sample times from the sample times file
	sampleTargets = []
	observations = sampleObservations()
	if arg.observations is not None:
		sampleTimesFile = open(arg.observations, 'rt')
		targetNames = []
		for line in sampleTimesFile:
			items = line.strip().split('\t')
			name = str(items[0])
			HJD = float(items[1])
			if name not in targetNames:
				observations.addTarget(name)
				targetNames.append(name)
			observations.addDataToTarget(name, HJD)
		sampleTimesFile.close()
	
	if arg.fake is not None:
		earliestObs, latestObs = observations.getExtremes()
		obsLength = latestObs - earliestObs
		if arg.baseline is not None:
			obsLength = arg.baseline * 365.
		numObs = arg.fake
		fakeObs = [earliestObs + obsLength/numObs * o for o in range(numObs)]
		observations.addTarget('fake')
		for f in fakeObs:
			observations.addDataToTarget('fake', f)
	
	# observations.dumpTargets()
	print "Loaded sample observation times for %d sample objects. On average, %.1f observations per target."%(observations.getNumber(), observations.getAverageObservationsPerTarget())
	
	
	nsteps = arg.steps
	stepSize = (logpStop-logpStart)/nsteps
	logps = []
	probabilities = []
	for logp in numpy.arange(logpStart, logpStop + stepSize, stepSize):
		period = 10**logp
		nSimulations = arg.simulations
		detections = 0
		for simulation in range(nSimulations):
			cosi = numpy.random.rand()
			inclination = numpy.arccos(cosi)
			phase = numpy.random.rand()
			K2 = massFunction(0.6, 0.2, period, inclination)
			# print "Period: %f [days], inclination: %f [deg], K2: %f km/s"%(period, inclination/numpy.pi*180, K2) 
			
			rvSTDs = []
			index = int(numpy.random.rand() * observations.getNumber())
			obs = observations.getObs(index)
			measuredRVs = []
			for date in obs['HJD']:
				measuredRV = rv(K2, period, phase, date)
					# print date, measuredRV
				measuredRVs.append(measuredRV)
			
			rvScatter = numpy.std(measuredRVs)
		
			# print "RV scatter is ", rvScatter 
			if rvScatter > arg.sigma: 
				detections+=1 
				# print "Detected"
		probability = float(detections)/float(nSimulations)
		print "logP %f -> detection probability: %f"%(logp, probability)
		probabilities.append(probability)
		logps.append(logp)
	
	
	
	figure = plt.figure()
	plt.plot(logps, probabilities)
	plt.xlabel('$log_{10}(P_{orb})$ [d]')
	plt.ylabel('Probability of detection')
	axes = plt.gca()
	axes.set_xlim(logpStart, logpStop)
	axes.set_ylim(0, 1.1)
	
	plt.show(block=False)
	if arg.dump is not None:
		plt.savefig(arg.dump + '.pdf', format = 'pdf')
	generalUtils.query_yes_no("Continue?")
	
	sys.exit()
	
	
