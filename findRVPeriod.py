#!/usr/bin/env python
import sys, os
import numpy, math
import argparse
import loadingSavingUtils
import spectrumClasses, timeClasses
import scipy.optimize
import copy
import ppgplot

def sine(x, a0, a1):
    y = a0 * numpy.sin(2.*math.pi*(x + a1))
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
	y = gamma + amplitude * numpy.sin(2. * numpy.pi * (x + phase)) 
	return y

def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")

class dataLog:
    def __init__(self):
        self.measurements = []
        
    def addMeasurement(self, HJD, velocity, velocityError, fwhm, wavelength):
        measurement = {}
        measurement['HJD'] = HJD
        measurement['velocity'] = velocity
        measurement['velocityError'] = velocityError
        measurement['fwhm'] = fwhm
        measurement['wavelength'] = wavelength
        
        # Check if the measurement is a duplicate (ie same HJD)
        found = -1
        for index, m in enumerate(self.measurements):
            if m['HJD'] == HJD: 
                found = index
                m['velocity'] = velocity
                m['velocityError'] = velocityError
                m['fwhm'] = fwhm
                m['wavelength'] = wavelength
        if found == -1:
            self.measurements.append(measurement)
            
    def writeToFile(self, filename):
        logFile = open(filename, 'wt')
        logFile.write("HJD, velocity, velocity_error, fwhm, wavelength\n")
        for m in self.measurements: 
            outString = "%10.10f, %10.10f, %10.10f, %10.10f, %10.10f\n"%(m['HJD'], m['velocity'], m['velocityError'], m['fwhm'], m['wavelength'])
            logFile.write(outString)
        logFile.close()	
        
    def getSavedValues(self, HJD):
        for m in self.measurements:
			if (float(m['HJD']) == float(HJD)): 
				return (m['wavelength'], m['fwhm'])
        return (-1, -1)
        
    def sortByHJD(self):
        self.measurements =  sorted(self.measurements, key=lambda object: object['HJD'], reverse = False)
        return
        
    def getData(self):
        HJD = []
        velocity = []
        velErr = []
        for m in self.measurements:
            HJD.append(m['HJD'])
            velocity.append(m['velocity'])
            velErr.append(m['velocityError'])
        return (HJD, velocity, velErr)
        
    def loadFromFile(self, filename):
        if not os.path.exists(filename): return
        inputFile = open(filename, 'rt')
        for line in inputFile:
            parts = line.strip().split(',')
            if parts[0].strip(',') == 'HJD': continue 
            try:
                HJD = float(parts[0].strip(','))
                velocity = float(parts[1].strip(','))
                velocityError = float(parts[2].strip(','))
                fwhm = float(parts[3].strip(','))
                wavelength = float(parts[4].strip(','))
            except ValueError:
                print "Warning: Could not convert one of the values in line.. ", parts
                fwhm = 0
                wavelength = 0
            self.addMeasurement(HJD, velocity, velocityError, fwhm, wavelength)
        inputFile.close()
        
    def __str__(self):
        if len(self.measurements)==0:
            return
        retStr = "HJD, velocity, velocityErr, fwhm, wavelength\n"
        for m in self.measurements: 
            retStr+="%10.10f, %f, %f, %f, %f\n"%(m['HJD'], m['velocity'], m['velocityError'], m['fwhm'], m['wavelength'])
        return retStr

if __name__ == "__main__":	
    parser = argparse.ArgumentParser(description='Loads a CSV file containing HJDs and RVs and tries to fit a period and sinusoid to the data.')
    parser.add_argument('inputfile', type=str, help='Filename of the CSV file containing the RVs.')
    parser.add_argument('--device', type=str, default = "/xs", help='[Optional] PGPLOT device. Defaults to "/xs".')
    parser.add_argument('--title', type=str, help='Title for the plot. Otherwise title will be generated from data in the .CSV file.')
	
    arg = parser.parse_args()
	# print arg
		
    # Load the CSV file into a data object of class dataLog()
    recordedData = dataLog()
    logFilename = arg.inputfile
    
    recordedData.loadFromFile(logFilename)
		
    recordedData.sortByHJD()
    print recordedData
    (dates, velocities, velErrors) = recordedData.getData()
	
    mainPGPlotWindow = ppgplot.pgopen(arg.device)	
    ppgplot.pgask(False)
    pgPlotTransform = [0, 1, 0, 0, 0, 1]
    
    ppgplot.pgslct(mainPGPlotWindow)	
    ppgplot.pgsci(1)
    reducedDates = [d - 2457000 for d in dates]
    startDate = min(reducedDates)
    endDate = max(reducedDates)
    maxVel = max(velocities)
    minVel = min(velocities)
    velRange = max([maxVel, abs(minVel)])
    print "Start date: %f, End date: %f, Velocity range: %f km/s"%(startDate, endDate, velRange)
    ppgplot.pgenv(0, len(dates), - 1.2 * velRange, 1.2 * velRange, 0, 0 )
    #ppgplot.pgenv(startDate, endDate, -velRange, velRange, 0, 0)
    
    ppgplot.pgpt(range(len(dates)), velocities)
    ppgplot.pgerrb(2, range(len(dates)), velocities, velErrors, 0)
    ppgplot.pgerrb(4, range(len(dates)), velocities, velErrors, 0)

    ppgplot.pglab("Point number", "Radial velocity km/s", arg.inputfile)
    
    
    x = numpy.array(reducedDates)
    y = numpy.array(velocities)
    # Subtract the mean from the y-data
    y_mean = numpy.mean(y)
    y = y - y_mean
    import scipy.signal as signal
    
    periods = numpy.linspace(0.01, 10.0, 8000)
    ang_freqs = 2 * numpy.pi / periods
    power = signal.lombscargle(x, y, ang_freqs)
    # normalize the power
    N = len(x)
    power *= 2 / (N * y.std() ** 2)

    pgramPGPlotWindow = ppgplot.pgopen(arg.device)	
    ppgplot.pgask(False)
    pgPlotTransform = [0, 1, 0, 0, 0, 1]

    ppgplot.pgslct(pgramPGPlotWindow)
    ppgplot.pgenv(min(periods), max(periods), 0, max(power), 0, 0)
    ppgplot.pgline(periods, power)
    ppgplot.pglab("Period (d)", "Amplitude", "Lomb-Scargle: " + arg.inputfile[:-4])
	
    bestPeriod = periods[numpy.argmax(power)]
    print "Best period: %f days or %f hours"%(bestPeriod, bestPeriod * 24.)
    
    # Now plot a folded RV curve
    t0 = reducedDates[0]
    periodDays = bestPeriod
    phaseFirst = [((d - t0) % periodDays) / periodDays for d in reducedDates]
    velocityPlot = copy.deepcopy(velocities)
    velocityErrorPlot = copy.deepcopy(velErrors)
    phases = copy.deepcopy(phaseFirst)
    for index, p in enumerate(phaseFirst):
        phases.append(p + 1.0)
        velocityPlot.append(velocities[index])
        velocityErrorPlot.append(velErrors[index])
    
    ppgplot.pgslct(mainPGPlotWindow)
    ppgplot.pgsci(1)
    ppgplot.pgenv(0, 2.0, -velRange*1.2, velRange*1.2, 0, 0)
    ppgplot.pgpt(phases, velocityPlot)
    ppgplot.pgerrb(2, phases, velocityPlot, velocityErrorPlot, 0)
    ppgplot.pgerrb(4, phases, velocityPlot, velocityErrorPlot, 0)
    ppgplot.pglab("Phase", "Radial velocity km/s", arg.inputfile[:-4])
    
    # Now fit a sine wave with this period to the data...
    frequency = 1 / periodDays
    phase = 0
    gamma = 0
    amplitude = velRange
    
    x_values = numpy.array(phaseFirst)
    y_values = numpy.array(velocities)
    y_errors = numpy.array(velErrors)
    
    guess = [gamma, amplitude, phase]
    upperBounds = [100, 200, 1]
    lowerBounds = [-100, 0 , 0]
    bounds = (lowerBounds, upperBounds)
    print "Guess: ", guess
    results, covariance = scipy.optimize.curve_fit(sinePhase, x_values, y_values, p0 = guess, sigma = y_errors, bounds = bounds)
    errors = numpy.sqrt(numpy.diag(covariance))
        
    gammaFit = results[0]
    gammaError = errors[0]
    amplitudeFit = results[1]
    amplitudeError = errors[1]
    phaseFit = results[2]
    phaseError = errors[2]
    print "Result of curve fit: "
    print "\tgamma velocity: \t%f[%f] km/s"%(gammaFit, gammaError)
    print "\tv sin i amplitude: \t%f[%f] km/s"%(amplitudeFit, amplitudeError)
    print "\tphase: \t\t\t%f[%f]"%(phaseFit, phaseError)
    
    xFit = numpy.arange(0, 2, 0.02)
    yFit = gammaFit + amplitudeFit * numpy.sin(2*numpy.pi*(xFit + phaseFit))
        
    ppgplot.pgsci(2)
    ppgplot.pgline(xFit, yFit)
    ppgplot.pgsci(3)
    ppgplot.pgsls(2)
    ppgplot.pgline([0, 2], [gammaFit, gammaFit])
    
    # Now try to tweak the period with a fixed amplitude and gamma velocity
    phaseGuess = 0
    guess = [periodDays, phaseGuess]
    gamma = gammaFit
    amplitude = amplitudeFit
    results, covariance = scipy.optimize.curve_fit(sineFixedGammaAmplitude, reducedDates, y_values, p0 = guess, sigma = y_errors)
    errors = numpy.sqrt(numpy.diag(covariance))
    
    tweakedPeriod = results[0]
    tweakedPeriodError = errors[0]
    phase = results[1]
    phaseError = errors[1]
    print "Tweaked Period: %f[%f]"%(tweakedPeriod, tweakedPeriodError)
    print "Tweaked Phase: %f[%f]"%(phase, phaseError)
    
    phases = [(d % tweakedPeriod)/tweakedPeriod for d in reducedDates]
    extendedPhases = numpy.asarray(phases)
    extendedVelocities = numpy.asarray(y_values)
    extendedVelocityErrors = numpy.asarray(y_errors)
    for index, p in enumerate(phases):
        extendedPhases = numpy.append(extendedPhases, p+1)
        extendedVelocities = numpy.append(extendedVelocities, y_values[index])
        extendedVelocityErrors = numpy.append(extendedVelocityErrors, y_errors[index])
    
    tweakedPeriodPlot = ppgplot.pgopen(arg.device)
    ppgplot.pgslct(tweakedPeriodPlot)
    ppgplot.pgask(False)
    ppgplot.pgenv(0, 2.0, -velRange*1.2, velRange*1.2, 0, 0)
    yMin = gamma - 1.2*amplitude
    yMax = gamma + 1.2*amplitude
    ppgplot.pgenv(0, 2.0, yMin, yMax, 0, 0)
    ppgplot.pgpt(extendedPhases, extendedVelocities, 0)
    ppgplot.pgerrb(2, extendedPhases, extendedVelocities, extendedVelocityErrors, 0)
    ppgplot.pgerrb(4, extendedPhases, extendedVelocities, extendedVelocityErrors, 0)
    ppgplot.pgsci(3)
    ppgplot.pgsls(2)
    ppgplot.pgline([0, 2], [gammaFit, gammaFit])
  
    
    
    ppgplot.pglab("Phase", "Radial velocity km/s", arg.inputfile[:-4] + " period: " + str(tweakedPeriod))
	
 
    
    sys.exit()
    
    # Try a bit more of a brute force approach
    frequency = bestFreq
    frequencyRange = numpy.arange(13.0, 17.0, 0.001)
    testedFreq = []
    chiSqMeasures = []
    amplitudes = []
    phaseMeasure = []
    for frequency in frequencyRange:
        period = 1.0/frequency
        amplitude = velRange
        
        t0 = reducedDates[0]
        phases = [((d - t0) % period) / period for d in reducedDates]

    
        phaseOffset = 0.0
        guess = numpy.array([amplitude, phaseOffset])
        x_values = phases
        y_values = velocities
        y_errors = velErrors
        results, covariance = scipy.optimize.curve_fit(sine, x_values, y_values, guess, y_errors)
        errors = numpy.sqrt(numpy.diag(covariance))

        amplitude = results[0]
        phase = results[1]
       
        chiSq = 0
        for i, p in enumerate(phases):
            chiSq = chiSq + (velocities[i] - sine(p, amplitude, phase))**2/(velErrors[i]**2)
        
        testedFreq.append(frequency)
        chiSqMeasures.append(chiSq)
        amplitudes.append(amplitude)
        phaseMeasure.append(phase)
        # print frequency, chiSq, amplitude, phase
    
    
        
    index = numpy.argmin(chiSqMeasures)
    lowestChiSquared = chiSqMeasures[index]
    reducedChiSquared = lowestChiSquared / (len(reducedDates)-2)
    bestFrequency = testedFreq[index]
    bestAmplitude = amplitudes[index]
    bestPhase = phaseMeasure[index]
    print "Least squares result: Amplitude: %f km/s, Frequency: %f cycles/day, Phase offset: %f"%(bestAmplitude, bestFrequency, bestPhase)
    print "Reduced Chi squared: ", reducedChiSquared
    chiSqPGPlotWindow = ppgplot.pgopen(arg.device)	
    ppgplot.pgask(False)
    pgPlotTransform = [0, 1, 0, 0, 0, 1]
    ppgplot.pgslct(chiSqPGPlotWindow)
    ppgplot.pgenv(min(testedFreq), max(testedFreq), 0, max(chiSqMeasures), 0, 0)
    ppgplot.pgline(testedFreq, chiSqMeasures)
    ppgplot.pglab("Frequency cycles/day", "Chi-squared", "Least squares: " + arg.inputfile)
	
    period = 1.0/bestFrequency
    t0 = reducedDates[0]
    phases = [((d - t0) % period) / period for d in reducedDates]
    extraPhases = [p + 1.0 for p in phases]
    

    ppgplot.pgslct(mainPGPlotWindow)
    ppgplot.pgsci(1)
    ppgplot.pgeras()
    ppgplot.pgenv(0, 2.0, -velRange, velRange, 0, 0)
    ppgplot.pgpt(phases, velocityPlot)
    ppgplot.pgerrb(2, phases, velocityPlot, velocityErrorPlot, 0)
    ppgplot.pgerrb(4, phases, velocityPlot, velocityErrorPlot, 0)
    ppgplot.pgpt(extraPhases, velocityPlot)
    ppgplot.pgerrb(2, extraPhases, velocityPlot, velocityErrorPlot, 0)
    ppgplot.pgerrb(4, extraPhases, velocityPlot, velocityErrorPlot, 0)
    ppgplot.pglab("Phase", "Radial velocity km/s", arg.inputfile + " %f hours"%(24.0 * period))
    
    lineX = numpy.arange(0, 2, 0.01)
    lineY = sine(lineX, bestAmplitude, bestPhase)
    
    print "Fit: %f km/s, %f cycles/day, %f days, %f hours"%(bestAmplitude, bestFrequency, period, 24. * period)
    
    ppgplot.pgsci(2)
    ppgplot.pgline(lineX, lineY)
   
    #################################################################################################
    if not query_yes_no("Ready to refine the fit?"):
        sys.exit()
    print "Starting with frequency of %f cycles/day and at a phase of %f."%(bestFrequency, bestPhase)
		
    x_values = reducedDates
    y_values = velocities
    y_errors = velErrors
        
    guess = [0, bestAmplitude, bestFrequency, 0]
    amplitude = bestAmplitude
    results, covariance = scipy.optimize.curve_fit(sineFreqPhase, x_values, y_values, guess, y_errors)
    errors = numpy.sqrt(numpy.diag(covariance))
        
    print results
    gamma = results[0]
    gammaError = errors[0]
    finalAmplitude = results[1]
    amplitudeError = errors[1]
    finalFrequency = results[2]
    frequencyError = errors[2]
    t1 = results[3]
    chiSq = 0
    for index, d in enumerate(reducedDates):
        chiSq+= (velocities[index] - sineFreqPhase(d, gamma, finalAmplitude, finalFrequency, t1))**2 / velErrors[index]**2
    reducedChiSq = chiSq/(len(reducedDates) -4)
    period = 1.0/finalFrequency
    periodError = period * frequencyError/finalFrequency
    t0 = reducedDates[0]
    phases = [(d + t1) % period / period  for d in reducedDates]
    extraPhases = [p + 1.0 for p in phases]

    ppgplot.pgslct(mainPGPlotWindow)
    ppgplot.pgsci(1)
    ppgplot.pgeras()
    ppgplot.pgenv(0, 2.0, -1.2*velRange, 1.2*velRange, 0, 0)
    ppgplot.pgpt(phases, velocityPlot)
    ppgplot.pgerrb(2, phases, velocityPlot, velocityErrorPlot, 0)
    ppgplot.pgerrb(4, phases, velocityPlot, velocityErrorPlot, 0)
    ppgplot.pgpt(extraPhases, velocityPlot)
    ppgplot.pgerrb(2, extraPhases, velocityPlot, velocityErrorPlot, 0)
    ppgplot.pgerrb(4, extraPhases, velocityPlot, velocityErrorPlot, 0)
    ppgplot.pglab("Phase", "Radial velocity km/s", arg.inputfile + " %f hours"%(24.0 * period))
    
    lineX = numpy.arange(0, 2, 0.01)
    lineY = sine(lineX, amplitude, 0)
    
    print "Fit: %f[%f] km/s, %f cycles/day, %f[%f] days, %f[%f] hours"%(finalAmplitude, amplitudeError, finalFrequency, period, periodError, 24. * period, periodError*24.)
    print "gamma velocity: %f[%f] km/s"%(gamma, gammaError)
    print "Reduced ChiSquared:", reducedChiSq
    ppgplot.pgsci(2)
    ppgplot.pgline(lineX, lineY)
    
    
    # Now write it to a PS file.
    psPlot = ppgplot.pgopen("fit.ps/ps")	
    ppgplot.pgask(False)
    pgPlotTransform = [0, 1, 0, 0, 0, 1]
    	
    ppgplot.pgsci(1)
    ppgplot.pgenv(0, 2.0, -1.2*velRange, 1.2*velRange, 0, 0)
    ppgplot.pgpt(phases, velocityPlot)
    ppgplot.pgerrb(2, phases, velocityPlot, velocityErrorPlot, 0)
    ppgplot.pgerrb(4, phases, velocityPlot, velocityErrorPlot, 0)
    ppgplot.pgpt(extraPhases, velocityPlot)
    ppgplot.pgerrb(2, extraPhases, velocityPlot, velocityErrorPlot, 0)
    ppgplot.pgerrb(4, extraPhases, velocityPlot, velocityErrorPlot, 0)
    ppgplot.pglab("Phase", "Radial velocity km/s", arg.inputfile + " %f hours"%(24.0 * period))
    
    lineX = numpy.arange(0, 2, 0.01)
    lineY = sine(lineX, amplitude, 0)
        
    ppgplot.pgsci(2)
    ppgplot.pgline(lineX, lineY)
		
    ppgplot.pgend()
    
        
    
