#!/usr/bin/env python3
import sys, os, argparse, re

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Extracts the periods from the ephemeris files.')
	arg = parser.parse_args()
	print(arg)

	fileList = os.listdir()

	# print(fileList)

	search_re = re.compile(".*.-ephem.dat")

	ephemList = []
	for file in fileList:
		match = search_re.match(file)
		if (match):
			ephemList.append(file)

	ephemList = sorted(ephemList)

	outfile = open('periods.dat', 'wt')
	for file in ephemList:
		object = file[0:8]
		print(object)
		fileHandle = open(file, 'rt')
		for line in fileHandle:
			if line[0]=='#': continue
			line = line.strip()
			if 'E ' in line:
				fields = line.split(' ')
				period = float(fields[1])
				print(period)
				outfile.write("%s, %f\n"%(object, period))
		fileHandle.close()

	outfile.close()
