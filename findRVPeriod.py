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
    y = a0 * numpy.sin(2.*math.pi*x + a1)
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
            HJD = float(parts[0].strip(','))
            velocity = float(parts[1].strip(','))
            velocityError = float(parts[2].strip(','))
            fwhm = float(parts[3].strip(','))
            wavelength = float(parts[4].strip(','))
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
    print dates
    reducedDates = [d - 2457000 for d in dates]
    startDate = min(reducedDates)
    endDate = max(reducedDates)
    maxVel = max(velocities)
    minVel = min(velocities)
    velRange = max([maxVel, abs(minVel)])
    print startDate, endDate, velRange
    ppgplot.pgenv(startDate, endDate, -velRange, velRange, 0, 0)
    ppgplot.pgpt(reducedDates, velocities)
    ppgplot.pglab("HJD - 2457000", "Radial velocity km/s", arg.inputfile)

    x = numpy.array(reducedDates)
    y = numpy.array(velocities)
    f = numpy.arange(0.01, 20, 0.01)
    print f
    import scipy.signal as signal
    pgram = signal.lombscargle(x, y, f)
    print pgram
    
    pgramPGPlotWindow = ppgplot.pgopen(arg.device)	
    ppgplot.pgask(False)
    pgPlotTransform = [0, 1, 0, 0, 0, 1]

    ppgplot.pgslct(pgramPGPlotWindow)
    ppgplot.pgenv(min(f), max(f), 0, max(pgram), 0, 0)
    ppgplot.pgline(f, pgram)
    ppgplot.pglab("Frequency cycles/day", "Amplitude", arg.inputfile)
	
    bestFreq = f[numpy.argmax(pgram)]
    bestPeriod = 24.*1.0/ bestFreq
    periodDays = 1.0/bestFreq
    print "Best frequency: %f cycles/day\nBest period: %f hours"%(bestFreq, bestPeriod)
    
    # Now plot a folded RV curve
    t0 = reducedDates[0]
    phaseFirst = [((d - t0) % bestPeriod) / bestPeriod for d in reducedDates]
    velocityPlot = copy.deepcopy(velocities)
    velocityErrorPlot = copy.deepcopy(velErrors)
    phases = copy.deepcopy(phaseFirst)
    for index, p in enumerate(phaseFirst):
        phases.append(p + 1.0)
        velocityPlot.append(velocities[index])
        velocityErrorPlot.append(velErrors[index])
    print phases

    ppgplot.pgslct(mainPGPlotWindow)
    ppgplot.pgsci(1)
    ppgplot.pgenv(0, 2.0, -velRange, velRange, 0, 0)
    ppgplot.pgpt(phases, velocityPlot)
    ppgplot.pgerrb(2, phases, velocityPlot, velocityErrorPlot, 0)
    ppgplot.pgerrb(4, phases, velocityPlot, velocityErrorPlot, 0)
    ppgplot.pglab("Phase", "Radial velocity km/s", arg.inputfile)
    
    # Try a bit more of a brute force approach
    frequency = bestFreq
    frequencyRange = numpy.arange(0.01, 20.0, 0.01)
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
        print frequency, chiSq, amplitude, phase
    
        
    index = numpy.argmin(chiSqMeasures)
    bestFrequency = testedFreq[index]
    bestAmplitude = amplitudes[index]
    bestPhase = phaseMeasure[index]
    
    
    chiSqPGPlotWindow = ppgplot.pgopen(arg.device)	
    ppgplot.pgask(False)
    pgPlotTransform = [0, 1, 0, 0, 0, 1]

    ppgplot.pgslct(chiSqPGPlotWindow)
    ppgplot.pgenv(min(testedFreq), max(testedFreq), 0, max(chiSqMeasures), 0, 0)
    ppgplot.pgline(testedFreq, chiSqMeasures)
	
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
   
    
        
    
