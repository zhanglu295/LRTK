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
			sys.stderr.write("\n\nERROR - %s - %s is empty, please check the config file you inputted!\n\n" % (time.asctime(), name_variable))
			time.sleep(5)
			if name_variable == "dbsnp":
				sys.stderr.write("\ndbsnp is not essential, the pipeline will go on\n")
			else:
				sys.exit(-1)
			return 0
#			sys.exit(-1)

	def Ref(self):
		abs_path = self.get_path('ref')
		return abs_path

	def Java(self):
		abs_path = self.get_path('java')
		return abs_path

	def Dbsnp(self):
		abs_path = self.get_path('dbsnp')
		return abs_path

	def Gatk(self):
		abs_path = self.get_path('gatk')
		return abs_path

	def HaplotypeCaller(self):
		abs_path = self.get_path('HaplotypeCaller_par')
		return abs_path
	
	def GenotypeGVCFs(self):
		abs_path = self.get_path('GenotypeGVCFs_par')
		return abs_path

def usage():
	alignment_usage = \
	'''
	call SNV and InDel
	Version: 1.0.0
	Dependents: Python (>=3.0), GATK (>=3.0)
	Last Updated Date: 2017-06-01
	Contact: meijp@foxmail.com

	Usage: python alignment.py <options>

	Options:
		-i --input, bam file
		-o --outputdir, the path of output directory
		-c --config, the path of configuration file [default: outdir/config/Reseq.config]
		-p --parallel, the number of parallele tasks for variant call [default: 3]
		-L --chrlist, file including chromosomes that would call variant
		-h --help, help info

	'''
	print(alignment_usage)

if __name__ == '__main__':
	if len(sys.argv) < 5:
		usage()
		sys.exit(-1)

	inputbam = None
	outputdir = None
	ConfigFile = None
	Chrlist = None
	ParallelNum = 3

	opts, args = getopt.gnu_getopt(sys.argv[1:], 'i:o:c:p:L:h:', ['input', 'outputdir', 'config', 'parallel', 'chrlist', 'help'])
	for o, a in opts:
		if o == '-i' or o == '--input':
			inputbam = a
		if o == '-o' or o == '--outputdir':
			outputdir = a
		if o == '-c' or o == '--config':
			ConfigFile = a
		if o == '-p' or o == '--parallel':
			ParallelNum = a
		if o == '-L' or o == '--chrlist':
			Chrlist = a
		if o == '-h' or o == '--help':
			usage()
			sys.exit(-1)

	if ConfigFile == None or Chrlist == None:
		script_abs_path = os.path.abspath(sys.argv[0])
		create_config_py = os.path.join(os.path.dirname(script_abs_path), "create_config.py")
		config_dir = os.path.join(outputdir, "config")
		if os.path.isdir(config_dir):
			pass
		else:
			os.mkdir(config_dir)
		tmpshell = os.path.join(config_dir, "cc.sh")
		wtmpshell = open(tmpshell, 'w')
		shell_line = " ".join(["python", create_config_py, "Reseq -o", config_dir, "\n"])
		wtmpshell.write(shell_line)
		wtmpshell.close()
		subprocess.call(["sh", tmpshell])
		subprocess.call(["rm", tmpshell])
		ConfigFile = os.path.join(config_dir, "Reseq.config")
		Chrlist = os.path.join(config_dir, "chrlist.txt")

	G = baseinfo()
	G.get_config(ConfigFile)
	javapath = G.Java()
	gatkpath = G.Gatk()
	ref = G.Ref()
	dbsnp = G.Dbsnp()
	HaplotypeCaller_par = G.HaplotypeCaller()
	GenotypeGVCFs_par = G.GenotypeGVCFs()

	vcfdir = os.path.join(outputdir, "vcf")
	if os.path.exists(vcfdir) and os.path.isdir(vcfdir):
		pass
	else:
		os.mkdir(vcfdir)

	shelldir = os.path.join(vcfdir, "shell")
	if os.path.exists(shelldir) and os.path.isdir(shelldir):
		pass
	else:
		os.mkdir(shelldir)

	Chrshell = list()
	rChrlist = open(Chrlist, 'r')
	ChrGvcflist = list()
	chrpriority = None
	AllChrVcf = None
	bychrvcfdir = os.path.join(vcfdir, "bychr")
	if os.path.isdir(bychrvcfdir):
		pass
	else:
		os.makedirs(bychrvcfdir)
	for ch in rChrlist:
		ch = ch.strip()
		chrshell = os.path.join(shelldir, ch + ".vcf.sh")
		chrgvcf = os.path.join(bychrvcfdir, ch + ".gvcf")
		chrvcf = os.path.join(bychrvcfdir, ch + ".vcf")
		wchrshell = open(chrshell, 'w')
		if os.path.isfile(dbsnp):
			shell_line = " ".join(["set -e\n", javapath, "-Djava.io.tmpdir=" + vcfdir, "-jar", gatkpath, "HaplotypeCaller -R", ref, "-I", inputbam, "--emitRefConfidence GVCF", HaplotypeCaller_par, "--dbsnp", dbsnp, "-L", ch, "-O", chrgvcf, "\n"])
		else:
			shell_line = " ".join(["set -e\n", javapath, "-Djava.io.tmpdir=" + vcfdir, "-jar", gatkpath, "HaplotypeCaller -R", ref, "-I", inputbam, "--emitRefConfidence GVCF", HaplotypeCaller_par, "-L", ch, "-O", chrgvcf, "\n"])
		wchrshell.write(shell_line)
		if os.path.isfile(dbsnp):
			shell_line = " ".join([javapath, "-Xmx2g -jar", gatkpath, "-R", ref, "GenotypeGVCFs --variant", chrgvcf, "-O", chrvcf, "--dbsnp", dbsnp, GenotypeGVCFs_par, "-L", ch, "\n"])
		else:
			shell_line = " ".join([javapath, "-Xmx2g -jar", gatkpath, "-R", ref, "GenotypeGVCFs --variant", chrgvcf, "-O", chrvcf, GenotypeGVCFs_par, "-L", ch, "\n"])
		ChrGvcflist.append(chrvcf)
		wchrshell.write(shell_line)
		wchrshell.close()

		Chrshell.append(chrshell)
		if chrpriority == None:
			chrpriority = str(ch)
			AllChrVcf = "--variant:" + ch + " " + chrgvcf
		else:
			chrpriority = chrpriority + "," + str(ch)
			AllChrVcf = AllChrVcf + " --variant:" + ch + " " + chrgvcf
	ShellNum = len(Chrshell)
	## ParallelNum = 3
	finishNum = 0
	pn = 0
	ts = 0
	shell_line = None
	while finishNum < ShellNum:
		if pn < ParallelNum:
			if pn == 0:
				tmpshell = os.path.join(shelldir, "tmp." + str(ts) + ".sh")
				wtmpshell = open(tmpshell, 'w')
				shell_line = "sh " + Chrshell[finishNum] + " &\n"
			else:
				shell_line = shell_line + "sh " + Chrshell[finishNum] + " &\n"
			pn = pn + 1

		finishNum = finishNum + 1
		if pn == ParallelNum:
			shell_line = shell_line + "wait\necho Done\n"
			wtmpshell.write(shell_line)
			wtmpshell.close()
			subprocess.call(["sh", tmpshell])
			pn = 0
			ts = ts + 1
			sys.stderr.write("Command: %s\n" % shell_line)

	if pn > 0:
		tmpshell = os.path.join(shelldir, "tmp." + str(ts) + ".sh")
		wtmpshell = open(tmpshell, 'w')
		shell_line = shell_line + "wait\necho Done\n"
		wtmpshell.write(shell_line)
		wtmpshell.close()
		subprocess.call(["sh", tmpshell])
		sys.stderr.write("Command: %s\n" % shell_line)

	UnphaseVcf = os.path.join(vcfdir, "all.vcf")
	tmpshell = os.path.join(shelldir, "combine_vcf.sh")
	wtmpshell = open(tmpshell, 'w')
	AllChrVcf = AllChrVcf.replace("--variant:", "-I")
	shell_line = " ".join(["set -e\n", javapath, "-Xmx3g -jar", gatkpath, "MergeVcfs -R", ref, AllChrVcf, "-o", UnphaseVcf, "\n"])
	wtmpshell.write(shell_line)
	wtmpshell.close()
	subprocess.call(["sh", tmpshell])

	sys.stderr.write("[ %s ] vcf file has been writeen to %s\n\n" % (time.asctime(), UnphaseVcf))
