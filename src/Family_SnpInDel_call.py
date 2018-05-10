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

	def Hapmap(self):
		abs_path = self.get_path('hapmap')
		return abs_path
	
	def Omni(self):
		abs_path = self.get_path('omni')
		return abs_path
	
	def High_confidence(self):
		abs_path = self.get_path('high_confidence')
		return abs_path
	
	def Mills(self):
		abs_path = self.get_path('mills')
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
		-i --input, gvcf file list
		-o --outputdir, the path of output directory
		-c --config, the path of configuration file [default: outdir/config/Reseq.config]
		-p --parallel, the number of parallele tasks for variant call [default: 4]
		-L --chrlist, file including chromosomes that would call variant
		-h --help, help info

	'''
	print(alignment_usage)

if __name__ == '__main__':
	if len(sys.argv) < 5:
		usage()
		sys.exit(-1)

	inputgvcf = None
	outputdir = None
	ConfigFile = None
	Chrlist = None
	ParallelNum = 4

	opts, args = getopt.gnu_getopt(sys.argv[1:], 'i:o:c:p:L:h:', ['input', 'outputdir', 'config', 'parallel', 'chrlist', 'help'])
	for o, a in opts:
		if o == '-i' or o == '--input':
			inputgvcf = a
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

	if ConfigFile == None or os.path.exists(ConfigFile) == False:
		sys.stderr.write("Configuration file does not exist, you can create it using 'python LRTK-seq.py Config'\n")
		sys.exit(-1)
	
	if Chrlist == None or os.path.exists(Chrlist) == False:
		Chrlist = os.path.dirname(ConfigFile) + "/chrlist.txt"
		sys.stderr.write("[ %s ] Warnings: chr list was not provided, and %s would be used instead." % (time.asctime(), Chrlist))

	G = baseinfo()
	G.get_config(ConfigFile)
	javapath = G.Java()
	gatkpath = G.Gatk()
	ref = G.Ref()
	dbsnp = G.Dbsnp()
	HaplotypeCaller_par = G.HaplotypeCaller()
	GenotypeGVCFs_par = G.GenotypeGVCFs()
	hapmap = G.Hapmap()
	omni = G.Omni()
	high_confidence = G.High_confidence()
	mills = G.Mills()

	vcfdir = outputdir
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
	ChrVcflist = list()
	chrpriority = None
	AllChrVcf = None
	bychrvcfdir = os.path.join(vcfdir, "bychr")
	tmpdir = vcfdir + "/tmp"
	if os.path.exists(tmpdir) == False or os.path.isdir(tmpdir) == False:
		os.mkdir(tmpdir)

	gvcf_dict = dict()
	rinputgvcf = open(inputgvcf, 'r')
	for inputgvcf_info in rinputgvcf:
		#(SampleId, gvcf_path) = re.split("\t", inputgvcf_info.strip())
		gvcf_path = inputgvcf_info.strip()
		rgvcf_path = open(gvcf_path, 'r')
		for gvcf_path_info in rgvcf_path:
			(chrId, chrGvcf) = re.split("\t", gvcf_path_info.strip())
			if chrId not in gvcf_dict:
				gvcf_dict[chrId] = chrGvcf
			else:
				gvcf_dict[chrId] = gvcf_dict[chrId] + " -V " + chrGvcf
		rgvcf_path.close()
	rinputgvcf.close()

	if os.path.isdir(bychrvcfdir):
		pass
	else:
		os.makedirs(bychrvcfdir)
	for ch in rChrlist:
		ch = ch.strip()
		if len(str(ch)) <= 5:
			chrshell = os.path.join(shelldir, ch + ".vcf.sh")
			logfile = chrshell + ".log"
			chrgvcf = os.path.join(bychrvcfdir, ch + ".gvcf.gz")
			chrvcf = os.path.join(bychrvcfdir, ch + ".vcf.gz")
			wchrshell = open(chrshell, 'w')
			shell_line = " ".join(["set -e\n", javapath, "-Djava.io.tmpdir=" + tmpdir, "-Xmx2g -jar", gatkpath, "CombineGVCFs -R", ref, "-V", gvcf_dict[ch], "-O", chrgvcf, "2>", logfile, "\n"])
			wchrshell.write(shell_line)
			shell_line = " ".join([javapath, "-Djava.io.tmpdir=" + tmpdir, "-Xmx2g -jar", gatkpath, "GenotypeGVCFs -R", ref, " -V", chrgvcf, "-O", chrvcf, "--dbsnp", dbsnp, "-L", ch, "2>>", logfile, "\n"])
			wchrshell.write(shell_line)
			shell_line = "echo Work is completed! > " + chrshell + ".sign"
			wchrshell.write(shell_line)
			wchrshell.close()

			ChrVcflist.append(chrvcf)
			Chrshell.append(chrshell)
			if chrpriority == None:
				chrpriority = str(ch)
				AllChrVcf = "-I " + chrvcf
			else:
				chrpriority = chrpriority + "," + str(ch)
				AllChrVcf = AllChrVcf + " -I " + chrvcf

	UnphaseVcf = os.path.join(vcfdir, "all.vcf.gz")
	tmpshell = os.path.join(shelldir, "combine_vcf.sh")
	logfile = tmpshell + ".log"
	wtmpshell = open(tmpshell, 'w')
	#AllChrVcf = AllChrVcf.replace("--variant:", "-I")
	shell_line = " ".join(["set -e\n", javapath, "-Xmx5g -jar", gatkpath, "MergeVcfs -R", ref, AllChrVcf, "-O", UnphaseVcf, "--CREATE_INDEX 2>", logfile, "\n"])
	wtmpshell.write(shell_line)
	recalFile = vcfdir + "/all.snp.recal"
	tranchesFile = vcfdir + "/all.snp.tranches"
	rscriptFile = vcfdir + "/all.snp.R"
	UnphaseVcf1 = vcfdir + "/all.SNP.recal.vcf.gz"
	shell_line = " ".join([javapath, "-Xmx5g -jar", gatkpath, "VariantRecalibrator -an QD -an MQRankSum -an ReadPosRankSum -an FS -an DP -an MQ -an SOR -mode SNP -R", ref, "-V", UnphaseVcf, "--resource hapmap,known=false,training=true,truth=true,prior=15.0:" + hapmap, "--resource omni,known=false,training=true,truth=true,prior=12.0:" + omni, "--resource 1000G,known=false,training=true,truth=false,prior=10.0:" + high_confidence, "--resource dbsnp,known=true,training=false,truth=false,prior=2.0:" + dbsnp, "-O", recalFile, "--tranches-file", tranchesFile, "--max-gaussians 4\n"])
	wtmpshell.write(shell_line)
	shell_line = " ".join([javapath, "-Xmx5g -jar", gatkpath, "ApplyVQSR -ts-filter-level 99.9 -mode SNP -R", ref, "-V", UnphaseVcf, "--recal-file", recalFile, "--tranches-file", tranchesFile, "-O", UnphaseVcf1]) + "\n"
	wtmpshell.write(shell_line)
	shell_line = " ".join(["mv", recalFile, recalFile + ".idx", tranchesFile, UnphaseVcf, UnphaseVcf + ".tbi", tmpdir]) + "\n"
	wtmpshell.write(shell_line)
	recalFile = vcfdir + "/all.indel.recal"
	tranchesFile = vcfdir + "/all.indel.tranches"
	rscriptFile = vcfdir + "/all.indel.R"
	UnphaseVcf2 = vcfdir + "/all.SNP.INDEL.VQSR.vcf.gz"
	shell_line = " ".join([javapath, "-Xmx5g -jar", gatkpath, "VariantRecalibrator -an QD -an FS -an MQ -an MQRankSum -an ReadPosRankSum -an SOR -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 90.0 -R", ref, "-V", UnphaseVcf1, "--resource mills,known=false,training=true,truth=true,prior=12.0:" + mills, "--resource dbsnp,known=true,training=false,truth=false,prior=2.0:" + dbsnp, "-O", recalFile, "--tranches-file", tranchesFile, "--max-gaussians 4"]) + "\n"
	wtmpshell.write(shell_line)
	shell_line = " ".join([javapath, "-Xmx5g -jar", gatkpath, "ApplyVQSR -ts-filter-level 99.9 -mode INDEL -R", ref, "-V", UnphaseVcf1, "--recal-file", recalFile, "--tranches-file", tranchesFile, "-O", UnphaseVcf2]) + "\n"
	wtmpshell.write(shell_line)
	shell_line = " ".join(["mv", recalFile, recalFile + ".idx", tranchesFile, UnphaseVcf1, UnphaseVcf1 + ".tbi", tmpdir]) + "\n"
	wtmpshell.write(shell_line)
	shell_line = " ".join(["echo Done! >", tmpshell + ".sign\n"])
	wtmpshell.write(shell_line)
	wtmpshell.close()

	qsub_shell_txt = outputdir + "/qsub.shell.txt"
	wqsub_shell_txt = open(qsub_shell_txt, 'w')
	for Chrshell_info in Chrshell:
		ch = os.path.basename(Chrshell_info)
		ch = ch.replace(".vcf.sh", "")
		qinfo = ch + "\t" + Chrshell_info + ":6G\t" + tmpshell + ":6G\n"
		wqsub_shell_txt.write(qinfo)
	wqsub_shell_txt.close()
	#subprocess.call(["sh", tmpshell])
