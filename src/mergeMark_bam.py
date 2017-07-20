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
		rbamlist = open(bamlist, 'r')
		for bamfile in rbamlist:
			bamfile = bamfile.strip()
			if multiple_bam == 0:
				multiple_bam = "I=" + bamfile
				bamfiles = bamfile
			else:
				multiple_bam = multiple_bam + " I=" + bamfile
				bamfiles = bamfiles + ", " + bamfile

		tmpshell = os.path.join(bamdir, "merge_bam.sh")
		wtmpshell = open(tmpshell, 'w')
		shell_line = " ".join([java_path, "-Xmx5g -jar", picard_path, "MergeSamFiles", multiple_bam, "O=" + mergedbam, picard_paramter, "\n"])
		wtmpshell.write(shell_line)
		wtmpshell.close()
		sys.stderr.write("[ %s ] merge bam files: %s \n\n" % (time.asctime(), bamfiles))
		subprocess.call(["sh", tmpshell])
		sys.stderr.write("[ %s ] bam files has been merged into one file: %s \n\n" % (time.asctime(), mergedbam))
		return(mergedbam)
		
	def picard_markdup(self, oribam, markedbam, java_path, picard_path, picard_paramter):
		metrics = markedbam.replace('bam', 'metrics.txt')

		bai = oribam.replace('bam', 'bai')
		isbai = 0
		if os.path.exists(bai):
			isbai = 1
		else:
			bai = oribam + '.bai'
			if os.path.exists(bai):
				isbai = 1
		if isbai:
			pass
		else:
			sys.stderr.write("bai not found, will produce bai file for bam automatically: %s\n" % bai)
			tmpshell = bai + '.sh'
			wtmpshell = open(tmpshell, 'w')
			shell_line = " ".join([java_path, '-jar', picard_path, 'BuildBamIndex', "I=" + oribam, "O=" + bai, "\n"])
			wtmpshell.write(shell_line)
			print(shell_line)
			wtmpshell.close()
			subprocess.call(["sh", tmpshell])
			subprocess.call(["rm", tmpshell])

		sys.stderr.write("[ %s ] marking duplication ... \n" % (time.asctime()))
		markedbamdir = os.path.dirname(markedbam)
		tmpshell = os.path.join(markedbamdir, 'mark_duplicate.sh')
		wtmpshell = open(tmpshell, 'w')
		shell_line = " ".join([java_path, '-jar', picard_path, 'MarkDuplicates', "I=" + oribam, "O=" + markedbam, "M=" + metrics, picard_paramter, "\n"])
		wtmpshell.write(shell_line)
		wtmpshell.close()
		subprocess.call(["sh", tmpshell])
		sys.stderr.write("[ %s ] marked bam file has been written to %s\n" % (time.asctime(), markedbam))
		return(markedbam)

def usage():
	merge_mark_usage = \
	'''
	merge multiple files belong to the same individual, and mark PCR duplication of the merged bam file
	Version: 1.0.0
	Dependents: Python (>=3.0), BWA, SAMtools, Picard (>=2.0)
	Last Updated Date: 2017-06-01
	Contact: meijp@foxmail.com

	Usage: python mergeMark_bam.py <options>

	Options:
		-i --input, file that include the path all bam files that belong to a single individual
		-o --output, path of the output bam file
		-c --config, the path of configuration file [default: outdir/config/QC.config]
		-h --help, help info

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

	if ConfigFile == None:
		script_abs_path = os.path.abspath(sys.argv[0])
		create_config_py = os.path.join(os.path.dirname(script_abs_path), "create_config.py")
		config_dir = os.path.join(outputdir, "config")
		if os.path.isdir(config_dir):
			pass
		else:
			os.mkdir(config_dir)
		tmpshell = os.path.join(config_dir, "cc.sh")
		wtmpshell = open(tmpshell, 'w')
		shell_line = " ".join(["python", create_config_py, "QC -o", config_dir, "\n"])
		wtmpshell.write(shell_line)
		wtmpshell.close()
		subprocess.call(["sh", tmpshell])
		subprocess.call(["rm", tmpshell])
		ConfigFile = os.path.join(config_dir, "QC.config")

	O = OUTERSOFT()

	G = baseinfo()
	G.get_config(ConfigFile)
	javapath = G.Java()
	picardpath = G.Picard()
	picard_merge_parameter = G.PicardMergeparameter()
	picard_mark_parameter = G.PicardMarkparameter()

	bamline = open(inputbamlist, 'r').readlines()
	MergedBam = None
	if len(bamline) > 1:
		bamdir = os.path.dirname(outputbam)
		MergedBam = os.path.join(bamdir, "multiple.merged.bam")
		MergedBam = O.picard_merge(inputbamlist, MergedBam, javapath, picardpath, picard_merge_parameter)
	else:
		MergedBam = bamline[0].strip()
		sys.stderr.write("[ %s ] there is only one bam file in the list: %s\n\n" % (time.asctime(), inputbamlist))

	Markedbam = O.picard_markdup(MergedBam, outputbam, javapath, picardpath, picard_mark_parameter)
	if os.path.exists(Markedbam):
		sys.stderr.write("[ %s ] duplication reads have been marked, and the output has been written to %s \n" % (time.asctime(), Markedbam))
	else:
		sys.stderr.write("\nERROR: %s has not been created, please check the warning info and try again\n\n" % Markedbam)
		sys.exit(-1)
