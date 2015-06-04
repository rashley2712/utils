import cmd, sys, os
import photmanto
import os
import readline
import generalUtils
import shlex

class photCommands(cmd.Cmd):
	"""Simple command processor example."""
	prompt = 'photmanto: '
	# Disable rawinput module use
	use_rawinput = False
	# Do not show a prompt after each command read
	prompt = ''
	
	def precmd(self, line):
		if len(line)==0: return line
		if line[0] == '#':
			# This line is a comment, therefore send a 'NOP' to the cmd processor to ignore.
			print "Comment:", line
			return "NOP"
		return line
		
	def do_save(self, line):
		""" save [prefix]
		Save the current session and data. If specified, a prefix will be added to *.session.ptm and *.data.ptm."""
		fileprefix = line.split(' ')[0]
		photmanto.saveSession(fileprefix)
		return
		
	def do_restore(self, line):
		""" restore [prefix]
		Restore the session and data from a previously saved session. All current session information will be lost. If specified, a prefix will be added to *.session.ptm and *.data.ptm."""
		fileprefix = line.split(' ')[0]
		photmanto.restoreSession(fileprefix)
		return
		
	def do_load(self, line):
		""" load [filename] [maxrows]
		Load a file of photometry data.
		Specify maxrows if you want to load less than the full data set. """
		params = line.split(' ')
		filename = params[0]
		if len(params)>1:
			try:
				maxRows = int(params[1])
			except ValueError:
				maxRows = 0
		else:
			maxRows = 0
		if filename == '':
			print "Please specify a filename..."
			return
			
		if not os.path.exists(filename):
			print "Could not find %s. Try '!ls' to list local files."%filename
			return
			
		extension = filename.split('.')[-1]
		if extension=='log':
			photmanto.loadFromLogFile(filename)
		if extension=='fits':
			photmanto.loadFromFITSFile(filename, maxRows)
		return
		
	def do_set(self, line):
		""" set [variable] [state]
		Set a state variable to a specific value  (eg "set plotter pgplot"). 
		To set properties for a specific slot, use 'set slot [slotID] [variable] [value]'"""
		params = line.split(' ')
		params = shlex.split(line)
		if len(params)<2: return
		variable = params[0]
		value = params[1]
		if variable=='slot':
			print "Setting slot value:", params[1:]
			if len(params)<4: 
				print "Too few parameters specified. 'set slot [slotID] [variable] [value]'."
				return
			try:
				slotID = int(params[1])
			except ValueError:
				print "Invalid slot id"
				return
			parameter = params[2]
			value = params[3]
			photmanto.setSlotProperty(slotID, parameter, value)
		else: 
			photmanto.setState(variable, value)
		return 
		
	def do_copy(self, line):
		""" copy [from] [to]
		Copy a slot with id:from to id:to. """
		params = line.split(' ')
		if len(params)<2:
			print "Please specifiy a source and a destination."
			return
		try:
			fromSlot = int(params[0])
			toSlot = int(params[1])
		except ValueError:
			print "Could not recognise source and/or destination. Make sure you use a number for the IDs."
			return
		photmanto.copySlot(fromSlot, toSlot)
		return
		
	def do_env(self, line):
		""" env
		List the current environment variables. """
		photmanto.showState()
		return
	
	def do_calc(self, line):
		""" calc [slotID] [operation]
		Calculate and create a new column in the slot specified by slotID. [operation] can be: 
		- 'bmjd' : Calculates the barycentric MJD from the MJD. Needs an observatory and target location defined. Creates a new column in the slot, called 'BMJD'."""
		params = line.split(' ')
		print "params:", params
		if len(params)<2: 
			print "Specify a slotID and an operation."
			return
		if params[0]=='': 
			print "Please specify a slot ID."
			return
		try: 
			slotID = int(params[0])
		except ValueError:
			print "Did not understand the slot ID:", params[0]
			return
		operation = params[1]
		if operation == 'bmjd': photmanto.calculateBMJD(slotID)
		return
		
	def do_normalise(self, line):
		""" normalise [slotID]
			Calculates the counts per second and puts it in a new column called 'countrate'. Needs a 'counts' column and an 'exposure' column to work from. 
		"""
		params = line.split(' ')
		if params[0]=='': 
			print "Please specify a slot ID."
			return
		try: 
			slotID = int(params[0])
		except ValueError:
			print "Did not understand the slot ID:", params[0]
			return
		photmanto.calculateCountRate(slotID)
		return
		
	def do_minutes(self, line):
		""" minutes [slotID]
			Calculates the minutes from the start of the run and adds a new 'minutes' column. Needs a BMJD or an MJD column to work from. 
		"""
		params = line.split(' ')
		if params[0]=='': 
			print "Please specify a slot ID."
			return
		try: 
			slotID = int(params[0])
		except ValueError:
			print "Did not understand the slot ID:", params[0]
			return
		photmanto.calculateMinutes(slotID)
		return
		
	def do_csv(self, line):
		""" csv [slotID] [filename]
		Write the 3 key columns (xaxis, yaxis, yerrors) to a CSV file."""
		params = line.split(' ')
		if len(params)<2: 
			print "Specify a slotID and a filename."
			return
		if params[0]=='': 
			print "Please specify a slot ID."
			return
		try: 
			slotID = int(params[0])
		except ValueError:
			print "Did not understand the slot ID:", params[0]
			return
		filename = params[1]
		photmanto.writeCSV(filename, slotID)
		return
		
	def do_div(self, line):
		""" div [slot A] [slot B] [slot D]
		Executes a division of [slot A] by [slot B] and places the result in [slot D]."""
		params = line.split(' ')
		try:
			slotA = int(params[0])
			slotB = int(params[1])
			slotD = int(params[2])
		except ValueError:
			print "Could not understand at least one of the slots given. You need to give me 3 slotIDs which are integers."
			return
		photmanto.divide(slotA, slotB, slotD)
		return

	def do_sigmaclip(self, line):
		""" sigmaclip [slot ID] sigma factor [optional step size]
		Performs a sigma-clip of the column in yColumn and creates a mask in the slot. """
		params = line.split(' ')
		try:
			slotID = int(params[0])
		except ValueError:
			print "Could not understand slot ID."
			return
		try:
			factor = float(params[1])
			stepSize = int(params[2])
		except IndexError:
			stepSize = 1
			factor = 4
			print "Using default values: factor: %f, range: %f"%(factor, stepSize)
		
		photmanto.sigmaclip(slotID, factor, stepSize)
		return
		
	def do_cat(self, line):
		""" cat [slot id]
		Prints out all of the data (line by line) in a slot."""
		params = line.split(' ')
		try:
			slotID = int(params[0])
		except ValueError:
			print "Could not understand slot ID."
			return
		photmanto.catSlot(slotID)
		return
		
	def do_removezeros(self, line):
		""" removezeros [slot id]
		Removes all data points from the slot where the y-axis value is equal to zero. """
		params = line.split(' ')
		try:
			slotID = int(params[0])
		except ValueError:
			print "Could not understand slot ID."
			return
		photmanto.removeZeros(slotID)
		return
		
		
		
	def do_plot(self, line):
		""" plot [slot_number]
		Plot the contents of the slot """
		slotIDs = generalUtils.parseIntegerList(line)
		if len(slotIDs)==0:
			print "Nothing to plot."
			return
		photmanto.plot(slotIDs)
		return
		
	def do_times(self, line):
		params = line.split(' ')
		print "params:", params
		if params[0]=='': 
			print "Please specify a slot ID."
			return
		try: 
			slotID = int(params[0])
		except ValueError:
			print "Did not understand the slot ID:", params[0]
			return
		photmanto.showTimes(slotID)
		return
	
	def do_NOP(self, line):
		""" NOP
		Do nothing."""
		return 
		
	def do_test(self, line):
		""" test
		Used for debugging. """
		print "Performing the test command [%s]"%line
		return
	
	def do_lsl(self, line):
		""" lsl [slotIDs]
		Alias for "ls -l"
		Show more info about what's in the slots. """
		if line=="":
			photmanto.listAllSlots(long=True)
			return
		slotIDs = generalUtils.parseIntegerList(line)
		if len(slotIDs)==0: photmanto.listAllSlots(long=True)
		else: photmanto.listSlots(slotIDs, long = True)
		
		return
		
	def do_ls(self, line):
		""" ls [slotIDs]
		Show info about what's in the slots. 
		'*' - all slots
		'1-5' - slots 1 through to 5 (inclusive)
		'1,3,4' - comma separated list of slots"""
		
		if line=="":
			photmanto.listAllSlots()
			return
		slotIDs = generalUtils.parseIntegerList(line)
		if len(slotIDs)==0: photmanto.listAllSlots()
		else: photmanto.listSlots(slotIDs)
		
		return
		
	def do_show(self, line):
		""" show [columns] [slotid]
		Show info about a slot (eg list columns available). """
		params = line.split(' ')
		if len(params)<2: return
		parameter = params[0]
		if parameter == 'columns':
			slotID = int(params[1])
			photmanto.showColumns(slotID)
		if parameter == 'header':
			slotID = int(params[1])
			photmanto.showHeader(slotID)
		
		return
	
	def do_quit(self, line):
		""" quit 
		Leave photmanto and exit to the shell. """
		print "Leaving photmanto. Goodbye."
		sys.exit()
		return True
	
	def do_shell(self, line):
		"Run a shell command"
		print "running shell command:", line
		output = os.popen(line).read()
		print output
		self.last_output = output
		
	def emptyline(self):
		return
	
	def do_EOF(self, line):
		return True
	
	def postloop(self):
		return True

	def do_shell(self, line):
		"Run a shell command"
		print "running shell command:", line
		output = os.popen(line).read()
		print output
		self.last_output = output
