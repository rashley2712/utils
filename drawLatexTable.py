#!/usr/bin/env python
import sys, os
import numpy, math
import argparse
import loadingSavingUtils, generalUtils
import spectrumClasses, timeClasses
import scipy.optimize
import copy
import matplotlib

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

if __name__ == "__main__":  
	parser = argparse.ArgumentParser(description='Loads RV values from the Google sheet and writes it as a LaTex table.')
	parser.add_argument('inputfile', type=str, help='Filename of the CSV file containing the RVs.')
	parser.add_argument('-o', '--output', type=str, default = 'table.tex', help='Output filename. (Default: table.tex)')
	parser.add_argument('-n', '--rows', type=int, help='Limit the number of rows to convert into LaTex')
	arg = parser.parse_args()
	outputFilename = arg.output	

	# Load the CSV file objects we are going to plot
	inputFile = open(arg.inputfile, 'rt')

	rvData = []
	for line in inputFile:
		fields = line.strip().split('\t')
		print(fields)
		try:
			dataObject = {}
			dataObject['objectname'] = str(fields[0])
			dataObject['hjd'] = float(fields[1])
			dataObject['rv'] = float(fields[2])
			dataObject['rv_error'] = float(fields[3])
			rvData.append(dataObject)
		except expression as identifier:
			print("Warning. Could not parse line: %s"%line.strip())
		
	inputFile.close()
	print("read %d data points"%len(rvData))

	
	print "Writing latex file to:", outputFilename
	outputFile = open(outputFilename, 'wt')
	heading = """\\begin{table}
\\begin{tabular}{ l  l  l  l }
	\\hline
	WD & HJD & Radial velocity & Error  \\\\
	   &     & (km\,s$^{-1}$)          & (km\,s$^{-1}$) \\\\
"""
	outputFile.write(heading)
	object = ""
	counter = 0
	for r in rvData:
		if object!=r['objectname']:
			outputFile.write("\t\\hline\n")
			object = r['objectname']
		outputFile.write("\t%s & %f & %4.1f & %2.1f \\\\\n"%(r['objectname'], r['hjd'], r['rv'], r['rv_error']))
		counter+=1
		if arg.rows is not None and counter>=arg.rows: break
	ending = """  \hline
  \end{tabular}
  \label{tab:rvdata}
\end{table}
"""
	outputFile.write(ending)	
	outputFile.close()
	
	