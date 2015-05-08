import cmd, sys, os
import photmanto
import readline

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
			
		photmanto.loadFromFITSFile(filename, maxRows)
		return
		
	def do_set(self, line):
		""" set [variable] [state]
		Set a state variable to a specific value  (eg "set plotter pgplot"). To set values for a specific slot, use 'set slot [slotID] [variable] [value]'"""
		params = line.split(' ')
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
		
	def do_env(self, line):
		""" env
		List the current environment variables. """
		photmanto.showState()
		return
		
		
	def do_plot(self, line):
		""" plot [slot_number]
		Plot the contents of the slot """
		slotNumber = int(line)
		photmanto.plot(slotNumber)
		return
	
	def do_NOP(self, line):
		"""NOP
		Do nothing."""
		return 
		
	def do_test(self, line):
		""" test
		Used for debugging. """
		print "Performing the test command [%s]"%line
		return
		
	def do_ls(self, line):
		""" ls
		Show info about what's in the slots. """
		
		photmanto.listAllSlots(line)
		
		self.do_show(line)
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
