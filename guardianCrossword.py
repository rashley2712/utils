#!/usr/bin/env python

import argparse
import datetime
import datetime
import sys
import urllib2
import os

if __name__ == "__main__":
	days = ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday']
	testBaseURL = "http://www.devicecloud.co.uk/crosswords/"
	baseURL = "http://static.guim.co.uk/crosswords/pdfs/"
	homeDIR = os.getenv("HOME")
	dropboxPath = homeDIR + "/Dropbox/Crosswords/"
	namePrefix = "gdn.quick."
	nameSuffix = ".pdf"

	parser = argparse.ArgumentParser(description='Downloads the Guardian Quick crosswords and saves (and archives) them to a Dropbox folder.')
	parser.add_argument('--date', default = 'today', type=str, help='Date for the crossword (default: today)')
	parser.add_argument('--test', action='store_true', help='Use the test URL instead of the real Guardian URL.')
 
	arg = parser.parse_args()
	print arg
	if arg.test:
		baseURL = testBaseURL

	todaysDate = datetime.datetime.now()
	requestedDate = todaysDate
	if arg.date!='today':
		try:
			inputDate = datetime.datetime.strptime(arg.date, '%Y-%m-%d')
			requestedDate = inputDate
		except ValueError:
			print "I am not able to understand the date input, please use YYYY-MM-DD"
			sys.exit()

	todayYear = todaysDate.year
	todayMonth = todaysDate.month
	todayDay = todaysDate.day
	todayDayOfWeek = todaysDate.weekday()
	
	requestedYear = requestedDate.year
	requestedDay = requestedDate.day
	requestedMonth = requestedDate.month
	requestedDayOfWeek = requestedDate.weekday()

	dayDifference = todaysDate - requestedDate

	print "Today is: %d-%02d-%02d %s"%(todayYear, todayMonth, todayDay, days[todayDayOfWeek])
	print "You have asked for: %d-%02d-%02d %s"%(requestedYear, requestedMonth, requestedDay, days[requestedDayOfWeek])
	if dayDifference.days<0:
		print "Your requested date is in the future, no crossword yet."
		sys.exit()
	if dayDifference.days>0:
		print 'Your date was %d days ago'%dayDifference.days
    

	if requestedDayOfWeek == 6:
		print "You are requesting a crossword for a Sunday. Try the Observer."
		sys.exit()

	dateString = "%d%02d%02d"%(requestedYear, requestedMonth, requestedDay)

	fullURL = baseURL + namePrefix + dateString + nameSuffix
	print "OK. Now we are ready to fetch: ", fullURL

	try:
		response = urllib2.urlopen(fullURL)
	except  urllib2.HTTPError as e:
		print "We got an error of:", e.code
		sys.exit()
	except urllib2.URLError as e:
		print e.reason
		sys.exit()

	#for property, value in vars(response).iteritems():
	#	print property, "| ", value

	headers = str(response.headers)
	#print "Headers:"
	#print headers

	startIndex = headers.find('Content-Type')
	startIndex+= len("Content-Type: ")
	endIndex = headers.find('\r\n', startIndex)
	contentType = headers[startIndex:endIndex]

	if contentType!='application/pdf':
		print "The server did not return a PDF object."
		sys.exit()

	pdfData = response.read()
	print "Fetched the data ok"

	outputFilename = dropboxPath + namePrefix + dateString + nameSuffix
	outputFile = open(outputFilename, 'w')

	outputFile.write(pdfData)

	outputFile.close()

	print "Written the file to:", outputFilename