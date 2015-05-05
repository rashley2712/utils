import cmd, sys, os
import photmanto

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
		
	def do_load(self, filename):
		""" load [filename]
		Load a file of photometry data """
		if filename == '':
			print "Please specify a filename..."
			return
		photmanto.loadFromFITSFile(filename)
		return
		
	def do_plot(self, line):
		""" plot [slot_number]
		Plot the contents of the slot """
		slotNumber = int(line)
		photmanto.plot(slotNumber)
		return
	
	def do_NOP(self, line):
		return 
		
	def do_test(self, line):
		print "Performing the test command [%s]"%line
		return
		
	def do_show(self, line):
		""" list [slot numbers]
		Show info about what's in the slots """
		photmanto.listAllSlots(line)
		return
	
	def do_quit(self, line):
		""" quit 
		Leave photmanto and exit to the shell """
		print "Leaving photmanto. Goodbye."
		sys.exit()
		return True
	
	def do_shell(self, line):
		"Run a shell command"
		print "running shell command:", line
		output = os.popen(line).read()
		print output
		self.last_output = output
		
	def empty
	
	def do_EOF(self, line):
		return True
	
	def postloop(self):
		return True
