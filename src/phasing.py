import os, sys, getopt
import time
import re, subprocess
import pysam

class baseinfo:
	def __init__(self):
		self.configdist = {}

	def get_config(self, configFile):
		o_config = open(configFile, 'r')
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
		print(name_variable)
		if name_variable in self.configdist:
			ABS_PATH = self.configdist[name_variable]
			if ABS_PATH == "None":
				ABS_PATH = " "
			return ABS_PATH
		else:
			sys.stderr.write("ERROR - %s - %s is empty, please check the config file you inputted!\n" % (time.asctime(), name_variable))
			sys.exit(-1)

	def Java(self):
		abs_path = self.get_path('java')
		return abs_path

	def ExtractHAIRS(self):
		abs_path = self.get_path('extractHAIRS')
		return abs_path

	def HAPCUT2(self):
		abs_path = self.get_path('HAPCUT2')
		return abs_path

	def Fgbio(self):
		abs_path = self.get_path('fgbio')
		return abs_path

	def HapCut_par(self):
		abs_path = self.get_path('HapCut_par')
		return abs_path

def usage():
	phasing_usage = \
	'''
	phasing
	Version: 1.0.0
	Dependents: Python (>=3.0), HapCut, Fgbio
	Last Updated date: 2017-06-01
	Contact: meijp@foxmail.com

	Usage: python phasing.py <options>

	Options:
		-i --input, bam file
		-v --vcf, unphased vcf, both compressed or uncompressed vcf are allowed
		-o --outputdir, the path of output directory
		-c --config, the path of configuration file [default: outdir/config/Reseq.config]
		-h --help, help info

	'''
	print(phasing_usage)

if __name__ == '__main__':
	if len(sys.argv) < 7:
		usage()
		sys.exit(-1)

	Markedbam = None
	UnphaseVcf = None
	ConfigFile = None
	outputdir = None

	opts, args = getopt.gnu_getopt(sys.argv[1:], 'i:o:c:v:h:', ['input', 'outputdir', 'config', 'vcf', 'help'])
	for o, a in opts:
		if o == '-i' or o == '--input':
			Markedbam = a
		if o == '-v' or o == '--vcf':
			UnphaseVcf = a
		if o == '-o' or o == '--outputdir':
			outputdir = a
		if o == '-c' or o == '--config':
			ConfigFile = a
		if o == '-h' or o == '--help':
			usage()
			sys.exit(-1)
		
	if ConfigFile == None:
		sys.stderr.write("configuration file (-c) does not exist! Please generate it using 'python LRTK_simple.py Config'\n")

	G = baseinfo()
	G.get_config(ConfigFile)

	javapath = G.Java()
	Path_extractHAIRS = G.ExtractHAIRS()
	Path_HAPCUT2 = G.HAPCUT2()
	PathSNAPSHOT = G.Fgbio()

	PhaseVcf = None
	vcfname = os.path.basename(UnphaseVcf)
	if vcfname.endswith('.gz'):
		b = str(vcfname)
		b = b[0:(len(b)-3)]
		newUnphaseVcf = outputdir + '/' + b
		subprocess.call(["gunzip -c", UnphaseVcf, ">", newUnphaseVcf])
		UnphaseVcf = newUnphaseVcf
		vcfname = b

	if vcfname.endswith(".vcf"):
		b = str(vcfname)
		PhaseVcf = outputdir + '/' + b[0:(len(b)-4)] + ".phased.vcf"
	else:
		PhaseVcf = outputdir + '/' + vcfname + ".phased.vcf"

	fragment_file = PhaseVcf.replace("phased.vcf", "phased.fragment")
	haplotype_output_file = PhaseVcf.replace("phased.vcf", "haplotype.txt")

	tmpdir = outputdir
	tmpshell = os.path.join(tmpdir, "phase.sh")
	wtmpshell = open(tmpshell, 'w')
	shell_line = " ".join([Path_extractHAIRS, "--bam", Markedbam, "--VCF", UnphaseVcf, "--out", fragment_file, "\n"])
	wtmpshell.write(shell_line)
	shell_line = " ".join([Path_HAPCUT2, "--fragments", fragment_file, "--vcf", UnphaseVcf, "--output", haplotype_output_file, "\n"])
	wtmpshell.write(shell_line)
	shell_line = " ".join([javapath, "-jar", PathSNAPSHOT, "HapCutToVcf -v", UnphaseVcf, "-i", haplotype_output_file, "-o", PhaseVcf, "\n"])
	wtmpshell.write(shell_line)
	wtmpshell.close()

	subprocess.call(["sh", tmpshell])
	print("[ %s ] phased vcf file has been writeen to %s\n\n" % (time.asctime(), PhaseVcf))

	readdict = dict()
	rfragment = open(fragment_file, 'r')
	for fragmentinfo in rfragment:
		fragmentinfolist = re.split("\s", fragmentinfo.strip())
		readdict[fragmentinfolist[1]] = fragmentinfolist[0]
	rfragment.close()

	PhasedBam = Markedbam.replace("bam", "addPhaseInfo.bam")
	rbam = pysam.AlignmentFile(Markedbam, 'rb')
	wbam = pysam.AlignmentFile(PhasedBam, 'wb', template = rbam)
	for read in rbam:
		newinfo = pysam.AlignedSegment()
		newinfo = read
		newtags = newinfo.tags
#		for eachtag in newinfo.tags:
#			if eachtag[0] != "BD" and eachtag[0] != "BI":
#				newtags.append(eachtag)
		addtag = ("HP", 0)
		if newinfo.query_name in readdict:
			addtag = ("HP", readdict[newinfo.query_name])
		newtags.append(addtag)
		newinfo.tags = newtags
		wbam.write(newinfo)
	rbam.close()
	wbam.close()
	print("[ %s ] phasing info has been supplemented: %s \n" % (time.asctime(), PhasedBam))
