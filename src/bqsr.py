import os, sys
import getopt
import time
import re, subprocess

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

	def Ref(self):
		abs_path = self.get_path('ref')
		return abs_path

	def Java(self):
		abs_path = self.get_path('java')
		return abs_path

	def Picard(self):
		abs_path = self.get_path('picard')
		return abs_path

	def Gatk(self):
		abs_path = self.get_path('gatk')
		return abs_path

	def Dbsnp(self):
		abs_path = self.get_path('dbsnp')
		return abs_path

	def Gold_indel(self):
		abs_path = self.get_path('gold_indel')
		return abs_path

	def Intervals(self):
		abs_path = self.get_path('intervals')
		return abs_path

class OUTERSOFT:
	def GATK_BQSR(self, oribam, bqsrbam, java_path, picard_path, gatk_path, ref_path, dbsnp_path, gold_indel_path, intervals_path):
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
		bamdir = os.path.dirname(bqsrbam)
		tmpdir = bamdir + "/tmp"
		if os.path.isdir(tmpdir):
			pass
		else:
			os.mkdir(tmpdir)
		shelldir = bamdir + "/shell"
		if os.path.isdir(shelldir):
			pass
		else:
			os.mkdir(shelldir)
		tmpshell = os.path.join(shelldir, "bqsr.sh")
		wtmpshell = open(tmpshell, 'w')
		shell_line = " ".join([java_path, "-Xmx4g -Djava.io.tmpdir=" + tmpdir, "-jar", gatk_path, "-T RealignerTargetCreator -nt 4 -R", ref_path, "-I", oribam, "-known", dbsnp_path, "-known", gold_indel_path, "-o", bqsrbam + ".realn.intervals -L", intervals_path, "\n"])
		wtmpshell.write(shell_line)
		shell_line = " ".join([java_path, "-Xmx4g -Djava.io.tmpdir=" + tmpdir, "-jar", gatk_path, "-T IndelRealigner -LOD 0.4 -R", ref_path, "-I", oribam, "-known", dbsnp_path, "-known", gold_indel_path, "-o", bqsrbam + ".realn.bam --targetIntervals", bqsrbam + ".realn.intervals\n"])
		wtmpshell.write(shell_line)
		shell_line = " ".join([java_path, "-Xmx4g -Djava.io.tmpdir=" + tmpdir, "-jar", gatk_path, "-T BaseRecalibrator -R", ref_path, "-I", bqsrbam + ".realn.bam", "-knownSites", dbsnp_path, "-knownSites", gold_indel_path, "-o", bqsrbam + ".realn.bam.grp\n"])
		wtmpshell.write(shell_line)
		shell_line = " ".join([java_path, "-Xmx4g -Djava.io.tmpdir=" + tmpdir, "-jar", gatk_path, "-T PrintReads -R", ref_path, "-I", bqsrbam + ".realn.bam", "-BQSR", bqsrbam + ".realn.bam.grp", "-o", bqsrbam, "\n"])
		wtmpshell.write(shell_line)
		shell_line = " ".join(["mv", bqsrbam + ".*", tmpdir, "\n"])
		wtmpshell.write(shell_line)
		wtmpshell.close()
		subprocess.call(["sh", tmpshell])
		sys.stderr.write("[ %s ] BQSR has been done: %s \n\n" % (time.asctime(), bqsrbam))
		return(bqsrbam)

def usage():
	bqsr_usage = \
	'''
	Base Quality Score Recalibration using GATK
	Version: 1.0.0
	Dependents: Python (>=3.0), GATK (>=2.0), Java
	Last Updated Date: 2017-06-01
	Contact: meijp@foxmail.com

	Usage: python bqsr.py <options>

	Options:
		-i --input, path of the input bam file
		-o --output, path of the output bam file
		-c --config, the path of configuration file [default: outdir/config/Basic.config]
		-h --help, help info

	'''
	print(bqsr_usage)

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
			inputbam = a
		if o == '-o' or o == '--output':
			outputbam = a
		if o == '-c' or o == '--config':
			ConfigFile = a
		if o == '-h' or o == '--help':
			usage()
			sys.exit(-1)

	if ConfigFile == None:
		sys.stderr.write("Configuration file does not exist, you can create it using create_config.py\n")
		sys.exit(-1)

	O = OUTERSOFT()

	G = baseinfo()
	G.get_config(ConfigFile)
	javapath = G.Java()
	picardpath = G.Picard()
	gatkpath = G.Gatk()
	refpath = G.Ref()
	dbsnppath = G.Dbsnp()
	gold_indelpath = G.Gold_indel()
	intervalspath = G.Intervals()

	outputbam = O.GATK_BQSR(inputbam, outputbam, javapath, picardpath, gatkpath, refpath, dbsnppath, gold_indelpath, intervalspath)

	if os.path.exists(outputbam):
		sys.stderr.write("[ %s ] duplication reads have been marked, and the output has been written to %s \n" % (time.asctime(), outputbam))
	else:
		sys.stderr.write("\nERROR: %s has not been created, please check the warning info and try again\n\n" % outputbam)
		sys.exit(-1)
