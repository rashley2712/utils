import httplib2
import os

from apiclient import discovery
from oauth2client import client
from oauth2client import tools
from oauth2client.file import Storage


# If modifying these scopes, delete your previously saved credentials
# at ~/.credentials/sheets.googleapis.com-python-quickstart.json
SCOPES = 'https://www.googleapis.com/auth/spreadsheets'
CLIENT_SECRET_FILE = 'client_secret.json'
CLIENT_SECRET_FILE = 'sheetsAPI.secret.json'
APPLICATION_NAME = 'Google Sheets API Python Quickstart'


def get_credentials():
    """Gets valid user credentials from storage.

    If nothing has been stored, or if the stored credentials are invalid,
    the OAuth2 flow is completed to obtain the new credentials.

    Returns:
        Credentials, the obtained credential.
    """
    home_dir = os.path.expanduser('~')
    credential_dir = os.path.join(home_dir, '.credentials')
    if not os.path.exists(credential_dir):
        os.makedirs(credential_dir)
    credential_path = os.path.join(credential_dir,
                                   'sheets.googleapis.com-python-quickstart.json')

    store = Storage(credential_path)
    credentials = store.get()
    if not credentials or credentials.invalid:
        flow = client.flow_from_clientsecrets(CLIENT_SECRET_FILE, SCOPES)
        flow.user_agent = APPLICATION_NAME
        if flags:
            credentials = tools.run_flow(flow, store, flags)
        else: # Needed only for compatibility with Python 2.6
            credentials = tools.run(flow, store)
        print('Storing credentials to ' + credential_path)
    return credentials

def getSampleData(spreadsheetId):
	credentials = get_credentials()
	http = credentials.authorize(httplib2.Http())
	discoveryUrl = ('https://sheets.googleapis.com/$discovery/rest?''version=v4')
	service = discovery.build('sheets', 'v4', http=http, discoveryServiceUrl=discoveryUrl)
	
	spreadsheetId = '11fsbzSII1u1-O6qQUB8P0RzvJ8MzC5VHIASsZTYplXc'
	rangeName = 'Sheet1!A1:E10'
	result = service.spreadsheets().values().get(spreadsheetId=spreadsheetId, range=rangeName).execute()
	
	values = result.get('values', [])
	return values

class gSheetObject:
	def __init__(self):
		self.docID = ""
		self.readings = []
		self.sheetName = "object"
		self.columns = ["HJD", "RV", "RV error", "width", "wavelength", "good"]
		
	def setObjectName(self, objectname):
		self.sheetName = objectname
		
	def hasReadingFor(self, HJD):
		for s in self.readings:
			if s['HJD'] == HJD: return True
		return False
		
	def getGoodFlag(self, HJD):
		for s in self.readings:
			if s['HJD'] == HJD and s['good']==1: return True
		return False
	
	def createSheet(self, objectname):
		sheetProperties = { "title": "testsheet"}
		sheetRequest = {
			"properties": {
				sheetProperties
			}
		}
		createRequest = {
			"addSheet": {
				object(sheetRequest)
			}}
			
		result = self.service.spreadsheets().create(body = createRequest).execute()
	

	def getFitByHJD(self, HJD):
		for s in self.readings:
			if s['HJD'] == HJD: return (s['width'], s['wavelength'])
		return False


	def addNewMeasurement(self, HJD, RV, RVerror, width, wavelength, good, overwrite=True):
		reading = {}
		reading['HJD'] = HJD
		reading['RV'] = RV
		reading['RV error'] = RVerror
		reading['width'] = width
		reading['wavelength'] = wavelength
		if good: reading['good'] = 1
		else: reading['good'] = 0
		if overwrite:
			index = -1
			for i, s in enumerate(self.readings):
				if s['HJD'] == HJD: index = i
			if index!=-1:
				self.readings[index] = reading
				print "Overwriting existing reading"
			else:
				self.readings.append(reading)
				print "Adding new reading"
		else:
			self.readings.append(reading)
		
		return len(self.readings)
		
		
	def initCredentials(self):
		self.credentials = get_credentials()
		self.http = self.credentials.authorize(httplib2.Http())
		self.discoveryUrl = ('https://sheets.googleapis.com/$discovery/rest?''version=v4')
		self.service = discovery.build('sheets', 'v4', http=self.http, discoveryServiceUrl=self.discoveryUrl)
		
	def setDocID(self, id):
		self.docID = id

	def writeAllReadings(self):
		values = []
		row = []
		for c in self.columns:
			row.append(c)
		values.append(row)
		for s in self.readings:
			row = []
			for c in self.columns:
				row.append(s[c])
			values.append(row)
		body = { 'values': values }
		rangeName = self.sheetName + "!A:F"
		result = self.service.spreadsheets().values().update( spreadsheetId=self.docID, range=rangeName, valueInputOption="USER_ENTERED", body=body).execute()
		print result
			
	def loadAllReadings(self):
		rangeCells = self.sheetName + "!A:F"
		
		result = self.service.spreadsheets().values().get( spreadsheetId=self.docID, range=rangeCells, valueRenderOption='UNFORMATTED_VALUE').execute()
		values = result.get('values', [])
		
		numRows = len(values)
		columns = values[0]
		print "Columns are:", columns
		print "There are %d rows."%numRows
		newData = []
		self.columns = [c.encode('ascii', 'ignore') for c in columns]
		# self.columns = columns
		for index in range(1, numRows):
			row = values[index]
			readingObject = {}
			for j, columnName in enumerate(self.columns):
				try: data = row[j]
				except IndexError:
					data = 0
				try:
					dataFloat = float(data)
					readingObject[columnName] = dataFloat
				except ValueError:
					readingObject[columnName] = data
			#print readingObject
			newData.append(readingObject)
		self.readings = newData
		return
		
		
	def writeSampleData(self):	
		values =  [
					["Item", "Cost", "Stocked", "Ship Date"],
					["Wheel", "$20.50", "4", "3/1/2016"],
					["Door", "$15", "2", "3/15/2016"],
					["Engine", "$100", "1", "30/20/2016"],
					["Totals", "=SUM(B2:B4)", "=SUM(C2:C4)", "=MAX(D2:D4)"]
				]
		body = { 'values': values }
		result = self.service.spreadsheets().values().update( spreadsheetId=self.docID, range="Sheet1!A1:D5", valueInputOption="USER_ENTERED", body=body).execute()
			
		# result = self.service.spreadsheets().values().batchUpdate(spreadsheetId=self.docID, body=data).execute()
		return None

def oldmain():
    """Shows basic usage of the Sheets API.

    Creates a Sheets API service object and prints the names and majors of
    students in a sample spreadsheet:
    https://docs.google.com/spreadsheets/d/1BxiMVs0XRA5nFMdKvBdBZjgmUUqptlbs74OgvE2upms/edit
    """
    credentials = get_credentials()
    http = credentials.authorize(httplib2.Http())
    discoveryUrl = ('https://sheets.googleapis.com/$discovery/rest?'
                    'version=v4')
    service = discovery.build('sheets', 'v4', http=http,
                              discoveryServiceUrl=discoveryUrl)

    spreadsheetId = '1BxiMVs0XRA5nFMdKvBdBZjgmUUqptlbs74OgvE2upms'
    spreadsheetId = '11fsbzSII1u1-O6qQUB8P0RzvJ8MzC5VHIASsZTYplXc'
    rangeName = 'Class Data!A2:E'
    rangeName = 'Sheet1!A1:E10'
    result = service.spreadsheets().values().get(
        spreadsheetId=spreadsheetId, range=rangeName).execute()
    values = result.get('values', [])

    if not values:
        print('No data found.')
    else:
        print('Name, Major:')
        for row in values:
            # Print columns A and E, which correspond to indices 0 and 4.
            print('%s, %s' % (row[0], row[1]))

