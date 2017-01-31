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
		self.columns = ["HJD", "RV", "RV error"]
		
	def add
		
	def initCredentials(self):
		self.credentials = get_credentials()
		self.http = self.credentials.authorize(httplib2.Http())
		self.discoveryUrl = ('https://sheets.googleapis.com/$discovery/rest?''version=v4')
		self.service = discovery.build('sheets', 'v4', http=self.http, discoveryServiceUrl=self.discoveryUrl)
		
	def setDocID(self, id):
		self.docID = id
		
	def getSampleData(self):
		http = self.credentials.authorize(httplib2.Http())
		discoveryUrl = ('https://sheets.googleapis.com/$discovery/rest?''version=v4')
		service = discovery.build('sheets', 'v4', http=http, discoveryServiceUrl=discoveryUrl)
		spreadsheetId = self.docID
		rangeName = 'Sheet1!A1:E10'
		result = service.spreadsheets().values().get(spreadsheetId=spreadsheetId, range=rangeName).execute()
		values = result.get('values', [])	
		return values
		
	def appendRVReading(self, HJD, RV, RVerror):
		values = [ [ HJD, RV, RVerror] ]
		body = { 'values': values }
		result = self.service.spreadsheets().values().append( spreadsheetId=self.docID, range="Sheet1!A1:A3", valueInputOption="USER_ENTERED", body=body).execute()
		
	def getReadingByHJD(self, HJD):
		rangeCells = self.sheetName + "!A:D"
		
		result = self.service.spreadsheets().values().get( spreadsheetId=self.docID, range=rangeCells).execute()
		print "All values result"
		print result
		values = result.get('values', [])
		for index, v in enumerate(values):
			print index, v
		
		
		
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

