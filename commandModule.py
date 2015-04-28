import cmd
import photmanto

class photCommands(cmd.Cmd):
	"""Simple command processor example."""

	def do_greet(self, person):
		"""greet [person]
		Greet the named person"""
		if person:
			print "hi,", person
		else:
			print 'hi'
		
	def do_load(self, filename):
		""" load [filename]
		Load a file of photometry data """
		if filename == None:
			print "Please specify a filename..."
			return
		photmanto.loadFromFITSFile(filename)
	
	def do_quit(self, line):
		return True
	
	#def do_EOF(self, line):
	#	return True
	
	def postloop(self):
		print "Leaving the command line processor. Goodbye."
