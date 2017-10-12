#!/usr/bin/env python
import sys, os
import numpy, math
import argparse
import loadingSavingUtils, generalUtils
import spectrumClasses, timeClasses
import scipy.optimize
import copy
import matplotlib
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

def getPeriodogram(pgplotHandle, xdata, ydata, plo, phi):
    x = numpy.array(xdata)
    # Subtract the mean from the y-data
    y_mean = numpy.mean(ydata)
    y = numpy.array(ydata - y_mean)
    periods = numpy.linspace(plo, phi, 8000)
    ang_freqs = 2 * numpy.pi / periods
    power = signal.lombscargle(x, y, ang_freqs)
    # normalize the power
    N = len(x)
    power *= 2 / (N * y.std() ** 2)
    ppgplot.pgslct(pgplotHandle)
    ppgplot.pgeras()
    ppgplot.pgenv(min(periods), max(periods), 0, max(power), 0, 0)
    ppgplot.pgline(periods, power)
    ppgplot.pglab("Period (d)", "Amplitude", "Lomb-Scargle: " + arg.objectname)
    bestPeriod = periods[numpy.argmax(power)]
    lc = ppgplot.pgqci()
    ls = ppgplot.pgqls()
    ppgplot.pgsci(3)
    ppgplot.pgsls(2)
    ppgplot.pgline([bestPeriod, bestPeriod], [0, max(power)])
    ppgplot.pgsci(lc)
    ppgplot.pgsls(ls)
    return bestPeriod

def fitPhasePlot(plotHandle, period, xdata, ydata, yerrors):
    print "Performing phase plot"
    print period
    t0 = xdata[0]
    phaseFirst = [((d - t0) % period) / period for d in xdata]
    yPlot = copy.deepcopy(ydata)
    yErrors = copy.deepcopy(yerrors)
    phases = copy.deepcopy(phaseFirst)
    for index, p in enumerate(phaseFirst):
        phases.append(p + 1.0)
        yPlot.append(ydata[index])
        yErrors.append(yerrors[index])
    
    ppgplot.pgslct(plotHandle)
    ppgplot.pgsci(1)
    ppgplot.pgenv(0, 2.0, numpy.min(yPlot)*1.2, numpy.max(yPlot)*1.2, 0, 0)
    ppgplot.pgpt(phases, yPlot)
    ppgplot.pgerrb(2, phases, yPlot, velocityErrorPlot, 0)
    ppgplot.pgerrb(4, phases, yPlot, yErrors, 0)
    ppgplot.pglab("Phase", "Radial velocity km/s", "Phase plot: period: %f days   %f hours"%(period, period*24))
    
    gammaGuess = numpy.mean(ydata)
    guess = [gammaGuess, numpy.max(ydata), 0]
    upperBounds = [100, 200, 1]
    lowerBounds = [-100, 0 , 0]
    bounds = (lowerBounds, upperBounds)
    print "Guess: ", guess
    results, covariance = scipy.optimize.curve_fit(sinePhase, phaseFirst, ydata, p0 = guess, sigma = yerrors, bounds= bounds)
    errors = numpy.sqrt(numpy.diag(covariance))
    print results
    gammaFit = results[0]
    gammaError = errors[0]
    amplitudeFit = results[1]
    amplitudeError = errors[1]
    phaseFit = results[2]
    phaseError = errors[2]
    print "Result of curve fit: "
    print "\tgamma velocity: \t%s[%s] km/s"%generalUtils.formatValueError(gammaFit, gammaError)
    print "\tv sin i amplitude: \t%s[%s] km/s"%generalUtils.formatValueError(amplitudeFit, amplitudeError)
    print "\tphase: \t\t\t%s[%s]"%generalUtils.formatValueError(phaseFit, phaseError)
    
    xFit = numpy.arange(0, 2, 0.02)
    yFit = gammaFit + amplitudeFit * numpy.sin(2*numpy.pi*(xFit + phaseFit))
    lc = ppgplot.pgqci()
    ls = ppgplot.pgqls()
        
    ppgplot.pgsci(2)
    ppgplot.pgline(xFit, yFit)
    ppgplot.pgsci(3)
    ppgplot.pgsls(2)
    ppgplot.pgline([0, 2], [gammaFit, gammaFit])
    ppgplot.pgsci(lc)
    ppgplot.pgsls(ls)
    
    return (gammaFit, amplitudeFit, phaseFit)

def sineFreeParameters(x, period, gamma, amplitude, phase):
    w = 2. * numpy.pi / period
    y = gamma + amplitude * numpy.sin(w * x + phase)
    return y

def finalFit(period, gamma, amplitude, phase, xdata, ydata, yerrors):
    print "Final fit, starting with guesses of:"
    print "\tGamma velocity = %f km/s"%gamma
    print "\tAmplitude = %f km/s"%amplitude
    print "\tPeriod = %f days"%period
    
    guess = [period, gamma, amplitude, 0]
    upperBounds = [arg.phi, 100, 400, 2*numpy.pi]
    lowerBounds = [0, -100, 0 , 0]
    bounds = (lowerBounds, upperBounds)
    results, covariance = scipy.optimize.curve_fit(sineFreeParameters, xdata, ydata, p0 = guess, sigma = yerrors, bounds=bounds)
    errors = numpy.sqrt(numpy.diag(covariance))
    print "Final results", results
    return zip(results, errors)
    
    
def phasePlot(plotData, plotHandle = -1, device='/xs'):
    # Calculate the phases for the xdata
    global xStart
    xdata = plotData['x']
    ydata = plotData['y']
    yerrors = plotData['yerrors']
    t0 = xdata[0]
    phases = [(d % plotData['period'])/plotData['period'] for d in xdata]
    extendedPhases = numpy.asarray(phases)
    extendedVelocities = numpy.asarray(ydata)
    extendedVelocityErrors = numpy.asarray(yerrors)
    for index, p in enumerate(phases):
        extendedPhases = numpy.append(extendedPhases, p+1)
        extendedVelocities = numpy.append(extendedVelocities, ydata[index])
        extendedVelocityErrors = numpy.append(extendedVelocityErrors, yerrors[index])
    
    if plotHandle == -1:
        plotHandle = ppgplot.pgopen(device)
    ppgplot.pgslct(plotHandle)
    ppgplot.pgenv(0, 2.0, numpy.min(ydata)*1.2, numpy.max(ydata)*1.2, 0, 0)
    periodHours = plotData['period'] * 24.
    periodHoursError = plotData['periodError'] * 24.
    ppgplot.pglab("Phase", "Radial velocity km/s", plotData['object']+ " period: %s[%s]d"%generalUtils.formatValueError(plotData['period'], plotData['periodError']) + " or %s[%s] hours."%generalUtils.formatValueError(periodHours, periodHoursError))
    ppgplot.pgpt(extendedPhases, extendedVelocities, 0)
    ppgplot.pgerrb(2, extendedPhases, extendedVelocities, extendedVelocityErrors, 0)
    ppgplot.pgerrb(4, extendedPhases, extendedVelocities, extendedVelocityErrors, 0)
    ppgplot.pgsci(3)
    ppgplot.pgsls(2)
    ppgplot.pgline([0, 2], [plotData['gamma'], plotData['gamma']])
    ppgplot.pgsci(1)
    ppgplot.pgsls(4)
    ppgplot.pgline([0, 2], [0, 0])
  
    xFit = numpy.arange(0, 2, 0.01)
    yFit = plotData['gamma'] + plotData['amplitude'] * numpy.sin(numpy.pi * 2.0 * xFit + plotData['phase'])
    ppgplot.pgsci(2)
    ppgplot.pgsls(1)
    ppgplot.pgline(xFit, yFit)
    return



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
		return self.pgram['freq'], self.pgram['power']

if __name__ == "__main__":  
	parser = argparse.ArgumentParser(description='Loads all of the data produced in ''rvanal'' package and plots them.')
	parser.add_argument('--device', type=str, default = "/xs", help='[Optional] PGPLOT device. Defaults to "/xs".')
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
			continue
			
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
		rvAnalInput = open(filename, 'rt')
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
		if plotIndex%2==0:
			ax1 = matplotlib.pyplot.subplot(plotsPerPage, 2, (plotIndex*2) + 1)
		else:
			ax1 = matplotlib.pyplot.subplot(plotsPerPage, 2, (plotIndex*2))

		phases = [o.ephemeris.getPhase(h) for h in o.HJD]
		velocities = o.velocities
		velErrors = o.velErrors
		gamma = o.ephemeris.gamma
		K2 = o.ephemeris.K2
		extendedPhases = copy.deepcopy(phases)
		extendedVelocities = copy.deepcopy(velocities)
		extendedVelErrors = copy.deepcopy(velErrors)
		for index, p in enumerate(phases):
			extendedPhases.append(p + 1.0)
			extendedVelocities.append(velocities[index])
			extendedVelErrors.append(velErrors[index])
		 
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

		frequency, power = o.getPgram()
		print len(frequency), "points in the periodogram"
		LSFrequency = frequency[numpy.argmax(power)]
		print LSFrequency, ' period:', 1.0/LSFrequency
		if plotIndex%2==0:
			ax2 = matplotlib.pyplot.subplot(plotsPerPage, 2, (plotIndex*2 + 3))
		else:
			ax2 = matplotlib.pyplot.subplot(plotsPerPage, 2, (plotIndex*2 + 2))
		
		matplotlib.pyplot.plot(frequency, power, linewidth=1.0, color='k', alpha=0.75)
		matplotlib.pyplot.ylim([0, 1.2 * max(power)])
		matplotlib.pyplot.plot([1/o.ephemeris.Period, 1/o.ephemeris.Period], [0,  1.2 * max(power)], linewidth=1.5, linestyle='--',  color='r', alpha=1)
		matplotlib.pyplot.plot([LSFrequency, LSFrequency], [0,  1.2 * max(power)], linewidth=1.5, linestyle='--',  color='b', alpha=1)
		
		if plotIndex%2==0: matplotlib.pyplot.ylabel('Power')
		matplotlib.pyplot.xlabel('Frequency [cycles/day]')
		
		zoomFraction = 0.1
		ax3 = inset_axes(ax2, width="30%", height="50%", loc=1)
		zoomedFrequency = []
		zoomedPower = []
		bestFrequency = 1.0/o.ephemeris.Period
		for index, f in enumerate(frequency):
			if (f > bestFrequency*(1-zoomFraction)) and (f < bestFrequency*(1+zoomFraction)):
				zoomedFrequency.append(f)
				zoomedPower.append(power[index])
			
		ax3.plot(zoomedFrequency, zoomedPower, linewidth=1.0, color='k', alpha=1.0)
		ax3.plot([bestFrequency, bestFrequency], [0, max(power)], color='r', linewidth=1.5, linestyle='--')
		ax3.plot([LSFrequency, LSFrequency], [0, max(power)], color='b', linewidth=1.5, linestyle='--')
		matplotlib.pyplot.yticks(visible=False)
		matplotlib.pyplot.xticks(visible=False)

		
		plotIndex+= 1
		
		if plotIndex%plotsPerPage == 0:
			print "End of page"
			matplotlib.pyplot.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=.45)
			matplotlib.pyplot.show(block=False)
			if arg.output!='none': 
				matplotlib.pyplot.savefig(arg.output + "_page_%d"%pageNumber + ".pdf")
				pageNumber+= 1
			plotIndex = 0
			phasedFoldedLightCurves = matplotlib.pyplot.figure(figsize=(8, 11))
			# matplotlib.pyplot.gcf().clear()
			
	
	print "End of page"
	matplotlib.pyplot.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=.45)
	matplotlib.pyplot.show(block=False)
			
	generalUtils.query_yes_no("Continue?")
		
	sys.exit()
	
    
    
