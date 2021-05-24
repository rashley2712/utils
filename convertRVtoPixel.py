#!/usr/bin/env python2
import sys, argparse, subprocess, re, os, shlex


def getFitskeyvalue(output, key):
	startIndex = output.find(key)
	endIndex = output.find('\n', startIndex + 1)
	lineString = output[startIndex:endIndex]
	tokens = lineString.split()
	value = tokens[1]
	return value
	

if __name__ == "__main__":
			
	parser = argparse.ArgumentParser(description='Converts a molly RV output file to from RVs to Pixel shifts given a dispersion factor in km/s/pixel.')
	parser.add_argument('inputfile', type=str, help='Name of the molly file produced by the "rvel" command.')
	parser.add_argument('dispersion', type=float, help='Dispersion factor in km/s/pixel.')
	arg = parser.parse_args()
	
	
	inputFile = open(arg.inputfile, 'r')
	
	shifts = []
	for line in inputFile:
		items = line.split()
		print items
		rvel = {}
		rvel['jd'] = float(items[0])
		rvel['unknown'] = float(items[1])
		rvel['rv'] = float(items[2])
		rvel['rv_error'] = float(items[3])
		shifts.append(rvel)
	inputFile.close()
	
	
	for rvel in shifts:
		rvel['pixel_shift'] = rvel['rv'] / arg.dispersion
		rvel['pixel_shift_error'] = rvel['rv_error'] / arg.dispersion
		
		
	outputFilename = arg.inputfile + ".pxl"
	outputFile = open(outputFilename, 'w')
	for rvel in shifts:
		outputFile.write("  %1.4f     %1.4f\n"%(rvel['pixel_shift'], rvel['pixel_shift_error']))
		
	outputFile.close()
	
		
	sys.exit()
	
	
