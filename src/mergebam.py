import os, sys, gzip
import getopt
import time
import subprocess
import re

class baseinfo:
	def __init__(self):
		self.configdist = {}

	def get_config(self, configFile):
		o_config = open (configFile, 'r')
		for config_line in o_config:
			config_line = config_line.strip()
			if len(config_line) == 0:
				pass
			elif config_line.startswith('#'):
				pass
			else:
				(name, path) = re.split(" = ", config_line)
				name = name.strip()
				path = path.strip()
				path = path.replace("'", "")
				path = path.replace('"', '')

				self.configdist[name] = path
#				print (name, path)
		o_config.close()

	def get_path(self, name_variable):
		if name_variable in self.configdist:
			ABS_PATH = self.configdist[name_variable]
			if ABS_PATH == "None":
				ABS_PATH = " "
			return ABS_PATH
		else:
			sys.stderr.write("ERROR - %s - %s is empty, please check the config file you inputted!\n" % (time.asctime(), name_variable))
			sys.exit(-1)

	def Samtools(self):
		abs_path = self.get_path('samtools')
		return abs_path

	def Java(self):
		abs_path = self.get_path('java')
		return abs_path

	def Picard(self):
		abs_path = self.get_path('picard')
		return abs_path

	def PicardMergeparameter(self):
		abs_path = self.get_path('picard_merge_parameter')
		return abs_path

	def PicardMarkparameter(self):
		abs_path = self.get_path('picard_mark_parameter')
		return abs_path

class OUTERSOFT:
	def picard_merge(self, bamlist, mergedbam, java_path, picard_path, picard_paramter):
		bamdir = os.path.dirname(mergedbam)
		multiple_bam = 0
		bamfiles = 0
		bamnum = 0
		rbamlist = open(bamlist, 'r')
		for bamfile in rbamlist:
			bamnum += 1
			bamfile = bamfile.strip()
			if multiple_bam == 0:
				multiple_bam = "I=" + bamfile
				bamfiles = bamfile
			else:
				multiple_bam = multiple_bam + " I=" + bamfile
				bamfiles = bamfiles + ", " + bamfile

		shelldir = bamdir + "/shell"
		if os.path.isdir(shelldir):
			pass
		else:
			os.mkdir(shelldir)
		tmpshell = os.path.join(shelldir, "merge_bam.sh")
		logfile = tmpshell + ".log"
		wtmpshell = open(tmpshell, 'w')
		if bamnum > 1:
			shell_line = " ".join([java_path, "-Xmx5g -jar", picard_path, "MergeSamFiles", multiple_bam, "O=" + mergedbam, picard_paramter, "2>", logfile, "\n"])
		else:
			shell_line = " ".join(["ln -sf", bamfiles, mergedbam, "\n"])
		wtmpshell.write(shell_line)
		shell_line = " ".join([java_path, "-jar", picard_path, "BuildBamIndex", "I=" + mergedbam, "O=" + mergedbam + ".bai", "2>>", logfile, "\n"])
		wtmpshell.write(shell_line)
		wtmpshell.close()
		sys.stderr.write("[ %s ] merge bam files: %s \n\n" % (time.asctime(), bamfiles))
		subprocess.call(["sh", tmpshell])
		sys.stderr.write("[ %s ] bam files has been merged into one file: %s \n\n" % (time.asctime(), mergedbam))
		return(mergedbam)

def usage():
	merge_mark_usage = \
	'''
	merge multiple BAMs from the same individual
	Version: 1.0.0
	Dependents: Python (>=3.0), BWA, SAMtools, Picard (>=2.0)
	Last Updated Date: 2017-11-14
	

	Usage: python mergeMark_bam.py <options>

	Basic options:
		-i --input, file including the path of all the BAMs from the same individual
		-o --output, path of the output BAMs
		-c --config, the path of configuration file [default: ./config/Basic.config]
		-h --help, print help info

	'''
	print(merge_mark_usage)

if __name__ == '__main__':
	if len(sys.argv) < 5:
		usage()
		sys.exit(-1)

	inputbamlist = None
	outputbam = None
	ConfigFile = None
	opts, args = getopt.gnu_getopt(sys.argv[1:], 'i:o:c:h:', ['input', 'output', 'config', 'help'])
	for o, a in opts:
		if o == '-i' or o == '--input':
			inputbamlist = a
		if o == '-o' or o == '--output':
			outputbam = a
		if o == '-c' or o == '--config':
			ConfigFile = a
		if o == '-h' or o == '--help':
			usage()
			sys.exit(-1)

	if ConfigFile == None or os.path.exists(ConfigFile) == False:
		sys.stderr.write("Configuration file does not exist, you can create it using 'python LRTK-seq.py Config'\n")
		sys.exit(-1)

	O = OUTERSOFT()

	G = baseinfo()
	G.get_config(ConfigFile)
	javapath = G.Java()
	picardpath = G.Picard()
	picard_merge_parameter = G.PicardMergeparameter()
	picard_mark_parameter = G.PicardMarkparameter()

	MergedBam = None
	bamdir = os.path.dirname(outputbam)
	MergedBam = O.picard_merge(inputbamlist, outputbam, javapath, picardpath, picard_merge_parameter)
	if os.path.exists(MergedBam) == False:
		sys.stderr.write("\nERROR: %s has not been created, please check the warning info and try again\n\n" % MergedBam)
		sys.exit(-1)
