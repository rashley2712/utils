#!/usr/bin/env python
import sys, argparse, subprocess, re, os, shlex


def getFitskeyvalue(output, key):
	startIndex = output.find(key)
	endIndex = output.find('\n', startIndex + 1)
	lineString = output[startIndex:endIndex]
	tokens = lineString.split()
	value = tokens[1]
	return value
	

if __name__ == "__main__":
	fitskeyCommand = "/storage/astro1/phsaap/software/gstar/bin/figaro/fitskeys"
			
	parser = argparse.ArgumentParser(description='Finds all the .sdf files in a list of files specified that have a certain binning configuration.')
	parser.add_argument('filelist', nargs='+', type=str, help='List of files to check.')
	parser.add_argument('--bx', type=int, help='x-Binning factor required')
	parser.add_argument('--by', type=int, help='y-Binning factor required')
	parser.add_argument('--output', type=str, default='out.list', help='Name of list file. (default: out.list)')
	arg = parser.parse_args()
	
	filteredList = []
	
	for filename in arg.filelist:
		status = subprocess.Popen([fitskeyCommand, filename], stdout = subprocess.PIPE)
		output = status.stdout.read()
		xBinningFactor = int(getFitskeyvalue(output, 'CCDXBIN'))
		yBinningFactor = int(getFitskeyvalue(output, 'CCDYBIN'))
		print filename, " - Binning factor: (%dx%d)"%(xBinningFactor, yBinningFactor)
		
		if xBinningFactor == arg.bx:
			if yBinningFactor == arg.by:
				filteredList.append(filename)
	
	print filteredList
	
	if len(filteredList>0):
		outputfile = open(arg.output, "w")
		for f in filteredList:
			outputfile.write(f+ '\n')
		outputfile.close()
		
	sys.exit()
	
	
