#!/usr/bin/env python
import sys, os
import numpy, math
import argparse

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Takes the output from Tom''s rvanal routines and rewrites them into a suitable -ephem file. This saves tedious copy''n''paste.')
	parser.add_argument('nthalias', type=str, help='A text file containing the output from ''nthalias''')
	parser.add_argument('-e', '--ephem', type=str, help='Optional ephem file to write the data to')
	arg = parser.parse_args()
	
	resultsFile = open(arg.nthalias, 'rt')
	for f in resultsFile:
		#print f.strip()
		fields = f.strip().split(' ')
		#print fields
		if "Constant term" in f: 
			gamma = float(fields[3])
			gamma_error = float(fields[5])
		if "Semi-amplitude" in f: 
			k2 = float(fields[2])
			k2_error = float(fields[4])
		if "Phase zero" in f: 
			t0 = float(fields[3])
			t0_error = float(fields[5])
		if "Period" in f: 
			E = float(fields[2])
			E_error = float(fields[4])
		if "Chi**2" in f: 
			Chi2 = float(fields[2])
			N = int(fields[4])
	resultsFile.close()
		
	print "T0 %f  +/- %f"%(t0, t0_error)
	print "E %f  +/- %f"%(E, E_error)
	print "Gamma %f  +/- %f"%(gamma, gamma_error)
	print "K2 %f  +/- %f"%(k2, k2_error)
	Chi2_red = Chi2/float(N-4)
	print "Chi2 %f  %d  Chi2_red: %f"%(Chi2, N, Chi2_red)
	
	if arg.ephem is not None:
		ephemFile = open(arg.ephem, 'rt')
		backupLines = []
		for f in ephemFile:
			if "J2000" in f:
				positionLine = f
			else:
				backupLines.append(f)	
					
		ephemFile.close()
		
		ephemFile = open(arg.ephem, 'wt')
		ephemFile.write("T0 %f\n"%(t0))
		ephemFile.write("T0_error %f\n"%(t0_error))
		ephemFile.write("E %f\n"%(E))
		ephemFile.write("E_error %f\n"%(E_error))
		ephemFile.write(positionLine)
		ephemFile.write("Gamma %f\n"%(gamma))
		ephemFile.write("Gamma_error %f\n"%(gamma_error))
		ephemFile.write("K2 %f\n"%(k2))
		ephemFile.write("K2_error %f\n"%(k2_error))
		ephemFile.write("Chi2 %f\n"%(Chi2))
		ephemFile.write("Chi2_red %f\n"%(Chi2_red))
		ephemFile.write("N %f\n"%(N))
		ephemFile.write("\n")
		
		for b in backupLines:
			ephemFile.write("# " + b)
		
		ephemFile.close()
		
