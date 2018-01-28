#!/usr/bin/env python
import sys, os
import numpy, math
import argparse
import loadingSavingUtils, generalUtils
import spectrumClasses, timeClasses
import scipy.optimize
import copy
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def sine(x, a0, a1):
    y = a0 * numpy.sin(2.*math.pi*(x + a1))
    return y

def sineGammaAmplitude(x, gamma, amplitude):
	y = gamma + amplitude * numpy.sin(2. * numpy.pi * x)
	return y

def sineFixedFreq(x, gamma, amplitude, phase):
    global frequency
    f = frequency
    y = gamma + amplitude * numpy.sin(2. * numpy.pi * (f * x - t0))
    return y

def sineFixedGammaAmplitude(x, period, phase):
    global gamma, amplitude
    w = 2. * numpy.pi / period
    y = gamma + amplitude * numpy.sin(w * x + phase)
    return y


def sineFreqPhase(x, gamma, amplitude, f, p):
    y = gamma + (amplitude * numpy.sin(2.*math.pi*f*(x + p)))
    return y


def sinePhase(x, gamma, amplitude, phase):
    y = gamma + amplitude * numpy.sin(2. * numpy.pi * x + phase)
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

	def addPgram(self, freq, chisq):
		self.pgram = {}
		self.pgram['freq'] = freq
		self.pgram['chisq'] = chisq
		chiMax = numpy.max(chisq)
		self.pgram['power'] = [1.0 - y/chiMax for y in chisq]

	def getPgram(self):
		return self.pgram['freq'], self.pgram['chisq']

	def getSampledPgram(self, n=10):
		sampledFreq = []
		sampledChisq = []
		for index in numpy.arange(0, len(self.pgram['chisq'])-n-1, n):
			# print index, self.pgram['freq'][int(index + n/2)], self.pgram['power'][index:index+n]
			sampledFreq.append(self.pgram['freq'][int(index + n/2)])
			sampledChisq.append(numpy.mean(self.pgram['chisq'][index:index+n]))
		print "New length", len(self.pgram['chisq']), len(sampledChisq)
		return sampledFreq, sampledChisq


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Loads all of the data produced in ''rvanal'' package and plots them.')
	parser.add_argument('objects', type=str, help='Object name list')
	parser.add_argument('--output', type=str, default='none', help='Output final plot to a .pdf file, specify the name. Default is ''none''')
	arg = parser.parse_args()
	# print arg

  	# Set up the matplotlib environment
  	generalUtils.setMatplotlibDefaults()
  	params = {	'axes.labelsize': 'large',
				'xtick.labelsize': 'large',
				'ytick.labelsize': 'large',
			}
	matplotlib.rcParams.update(params)

	# Load the list objects we are going to plot
	objects = []
	objectFile = open(arg.objects, 'rt')
	for f in objectFile:
		if len(f) < 2: continue
		o = object(f.strip())
		objects.append(o)

	for o in objects:
		filename = str(o.id) + ".tsv"
		try:
			rvFile = open(filename, 'rt')
		except:
			print "Could not find radial velocities data:", filename
			sys.exit()


		dates = []
		velocities = []
		velErrors = []

		print "HJD\t\tVelocity (km/s)\tVel error"
		for index, f in enumerate(rvFile):
			params = f.strip().split('\t')
			try:
				date = float(params[0])
				velocity = float(params[1])
				velError = float(params[2])
			except:
				print "Error parsing data:", f
				continue
			dates.append(date)
			velocities.append(velocity)
			velErrors.append(velError)
			print "%f\t%f\t%f"%(date, velocity, velError)
		rvFile.close()

		o.HJD = dates
		o.velocities = velocities
		o.velErrors = velErrors

		hasEphemeris = o.loadEphemeris()
		if hasEphemeris:
			print o.id, o.ephemeris
		else:
			print "Warning. No ephemeris for id: ", o.id
			sys.exit()

		filename = o.id + "_pgram.dat"
		try:
			rvAnalInput = open(filename, 'rt')
		except:
			print "Could not find pgram data:", filename
			sys.exit()

		freq = []
		chisq = []
		for line in rvAnalInput:
			if line[0]=='#': continue
			if len(line) < 3: continue
			params = line.strip().split(' ')
			freq.append(float(params[0]))
			chisq.append(float(params[1]))
		rvAnalInput.close()
		o.addPgram(freq, chisq)
		print "Loaded %d data points in the periodogram."%len(freq)

	phasedFoldedLightCurves = matplotlib.pyplot.figure(figsize=(8, 11))
	# fig, axes = matplotlib.pyplot.subplots(nrows=len(objects), ncols=2, figsize=(8, 11))
	# fig.tight_layout()


	plotsPerPage = 4

	plotIndex = 0
	pageNumber = 1
	for o in objects:
		plotsPending = True
		if plotIndex%2==0:
			ax1 = matplotlib.pyplot.subplot(plotsPerPage, 2, (plotIndex*2) + 1)
		else:
			ax1 = matplotlib.pyplot.subplot(plotsPerPage, 2, (plotIndex*2))

		phases = numpy.array([o.ephemeris.getPhase(h) for h in o.HJD])
		velocities = numpy.array(o.velocities)
		velErrors = numpy.array(o.velErrors)
		gamma = o.ephemeris.gamma
		K2 = o.ephemeris.K2
		extendedPhases = copy.deepcopy(phases)
		extendedVelocities = copy.deepcopy(velocities)
		extendedVelErrors = copy.deepcopy(velErrors)
		for index, p in enumerate(phases):
			extendedPhases = numpy.append(extendedPhases, (p + 1.0))
			extendedVelocities = numpy.append(extendedVelocities, velocities[index])
			extendedVelErrors = numpy.append(extendedVelErrors, velErrors[index])

		maxVel = max(velocities)
		minVel = min(velocities)
		velRange = (maxVel - minVel)/2.0

		yMin = gamma - 1.2*K2
		yMax = gamma + 1.2*K2
		xFit = numpy.arange(0, 2, 0.02)
		yFit = gamma + K2 * numpy.sin(2*numpy.pi*(xFit))

		matplotlib.pyplot.errorbar(extendedPhases, extendedVelocities, color='k', yerr=extendedVelErrors, fmt='.', capsize=0)
		matplotlib.pyplot.plot(xFit, yFit, color='k', linewidth=1.0)
		matplotlib.pyplot.plot([0, 2], [gamma, gamma], color='k', linestyle='--')
		if plotIndex%2==0: matplotlib.pyplot.ylabel('$K_{sec}$ velocity (km/s)')
		matplotlib.pyplot.xlabel('Phase')
		yPosition = ax1.get_ylim()[1] * 0.8
		matplotlib.pyplot.text(0.22, .85, o.id, fontsize='x-large', transform = ax1.transAxes)
		matplotlib.pyplot.text(0.8, 0.85, "%2.2fd"%o.ephemeris.Period, fontsize='x-large', transform = ax1.transAxes)


		extendedResiduals = numpy.array([ vel - (gamma + K2 * numpy.sin(2*numpy.pi*(phase))) for phase, vel in zip(extendedPhases, extendedVelocities)])
		residuals = numpy.array([ vel - (gamma + K2 * numpy.sin(2*numpy.pi*(phase))) for phase, vel in zip(phases, velocities)])
		sigma = numpy.std(residuals)
		chisq = 0
		for (r, e) in zip(residuals, velErrors):
			chisq+= (r/e)**2

		stddev = numpy.sqrt(chisq/len(residuals))
		#print "Stddev: ", stddev, " from numpy:", sigma
		rejected = numpy.array(numpy.abs(residuals)> 4*sigma)
		extendedRejected = numpy.array(numpy.abs(extendedResiduals)> 4*sigma)
		#print "Rejected:", rejected

		# print "Plotting residuals for %s"%(o.id)
		if plotIndex%2==0:
			ax2 = matplotlib.pyplot.subplot(plotsPerPage, 2, (plotIndex*2 + 3))
		else:
			ax2 = matplotlib.pyplot.subplot(plotsPerPage, 2, (plotIndex*2 + 2))

		extendedAccepted = numpy.logical_not(extendedRejected)
		matplotlib.pyplot.errorbar(extendedPhases[extendedAccepted], extendedResiduals[extendedAccepted], color='k', yerr=extendedVelErrors[extendedAccepted], fmt='.', capsize=0)
		matplotlib.pyplot.plot([0, 2], [sigma, sigma], color='grey', linestyle=':')
		matplotlib.pyplot.plot([0, 2], [0, 0], color='grey', linestyle='--')
		matplotlib.pyplot.plot([0, 2], [-sigma, -sigma], color='grey', linestyle=':')
		print residuals, rejected
		matplotlib.pyplot.errorbar(phases[rejected], residuals[rejected], color='k', yerr=velErrors[rejected], fmt='o', capsize=0, markerfacecolor='white')

		if plotIndex%2==0: matplotlib.pyplot.ylabel('km/s')
		matplotlib.pyplot.xlabel('Phase')


		plotIndex+= 1

		if plotIndex%plotsPerPage == 0:
			print "End of page"
			matplotlib.pyplot.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=.45)
			matplotlib.pyplot.show(block=False)
			if arg.output!='none':
				plotFilename = arg.output + "_page_%d"%pageNumber + ".pdf"
				print "Dumping plot to: ", plotFilename
				matplotlib.pyplot.savefig(arg.output + "_page_%d"%pageNumber + ".pdf")
			pageNumber+= 1
			plotIndex = 0
			phasedFoldedLightCurves = matplotlib.pyplot.figure(figsize=(8, 11))
			# matplotlib.pyplot.gcf().clear()
			plotsPending = False


	if plotsPending:
		print "End of page"
		matplotlib.pyplot.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=.45)
		matplotlib.pyplot.show(block=False)
		if arg.output!='none':
			plotFilename = arg.output + "_page_%d"%pageNumber + ".pdf"
			print "Dumping plot to: ", plotFilename
			matplotlib.pyplot.savefig(arg.output + "_page_%d"%pageNumber + ".pdf")

	generalUtils.query_yes_no("Continue?")

	# Write RV table
	rvFilename = "rvdata.tsv"
	rvFile = open(rvFilename, "wt")
	for o in objects:
		name = o.id
		for mjd, rv, rvErr in zip(o.HJD, o.velocities, o.velErrors):
			rvFile.write("%s\t%f\t%f\t%f\n"%(name, mjd, rv, rvErr))
	rvFile.close()


	sys.exit()
