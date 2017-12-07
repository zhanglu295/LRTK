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
	Last Updated date: 2017-11-14

	Usage: python phasing.py [options]

	Basic options:
		-i --input, BAM file
		-v --vcf, unphased vcf, either compressed or uncompressed vcf 
		-o --outputdir, the path of output
		-c --config, the path of configuration file [default: ./config/Reseq.config]
		-h --help, print help info

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
		
	if ConfigFile == None or os.path.exists(ConfigFile) == False:
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
	logfile = tmpshell + ".log"
	wtmpshell = open(tmpshell, 'w')
	shell_line = "set -e\n"
	wtmpshell.write(shell_line)
	shell_line = " ".join([Path_extractHAIRS, "--bam", Markedbam, "--VCF", UnphaseVcf, "--out", fragment_file, "2>", logfile, "\n"])
	wtmpshell.write(shell_line)
	shell_line = " ".join([Path_HAPCUT2, "--fragments", fragment_file, "--vcf", UnphaseVcf, "--output", haplotype_output_file, "2>>", logfile, "\n"])
	wtmpshell.write(shell_line)
	shell_line = " ".join([javapath, "-jar", PathSNAPSHOT, "HapCutToVcf -v", UnphaseVcf, "-i", haplotype_output_file, "-o", PhaseVcf, "2>>", logfile, "\n"])
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
	barcode_phase_file = PhasedBam.replace("addPhaseInfo.bam", "addPhaseInfo.barcode.csv")
	rbam = pysam.AlignmentFile(Markedbam, 'rb')
	wbam = pysam.AlignmentFile(PhasedBam, 'wb', template = rbam)
	wbarcode_phase_file = open(barcode_phase_file, 'w')
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
			phased_id = readdict[newinfo.query_name]
			for eachtag in newinfo.tags:
				if eachtag[0] == "BX":
					barcode_info = (eachtag[1]).replace("-1", "")
					barcode_phase_info = barcode_info + "\t" + str(phased_id) + "\t" + newinfo.query_name + "\n"
					wbarcode_phase_file.write(barcode_phase_info)
		newtags.append(addtag)
		newinfo.tags = newtags
		wbam.write(newinfo)
	rbam.close()
	wbarcode_phase_file.close()
	wbam.close()
	print("[ %s ] phasing info has been supplemented: %s \n" % (time.asctime(), PhasedBam))
