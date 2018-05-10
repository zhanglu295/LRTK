import os, sys, gzip
import getopt
import time
import subprocess
import re
import multiprocessing
import random, string
import pysam

### generate clean fastq files and barcode error correction

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

	def Bwa(self):
		abs_path = self.get_path('bwa')
		return abs_path

	def Ref(self):
		abs_path = self.get_path('ref')
		return abs_path

	def Fq_aln_par(self):
		abs_path = self.get_path('fq_aln_parameter')
		return abs_path

	def Fq_sam_par(self):
		abs_path = self.get_path('fq_sampe_parameter')
		return abs_path
	
	def Samtools(self):
		abs_path = self.get_path('samtools')
		return abs_path

	def Picard(self):
		abs_path = self.get_path('picard')
		return abs_path

	def Picard_merge_parameter(self):
		abs_path = self.get_path('picard_merge_parameter')
		return abs_path

def find_samtools_version(samtools_path, outdir = './'):
	randomstring = ''.join(random.sample(string.ascii_letters + string.digits, 8))
	tmpshell = outdir + "/" + randomstring + ".sh"
	tmplog = outdir + "/" + randomstring + ".log"
	wtmpshell = open(tmpshell, 'w')
	shell_line = " ".join([samtools_path, "2>", tmplog, "\n"])
	wtmpshell.write(shell_line)
	wtmpshell.close()
	subprocess.call(["sh", tmpshell])

	sv = 0
	rlog = open(tmplog, 'r')
	for log in rlog:
		if re.search("Version", log):
			loginfo = re.split('\s', log)
			sv = (re.split('\.', loginfo[1]))[0]
	rlog.close()
	subprocess.call(["rm", tmpshell, tmplog])
	return(int(sv))

class OUTERSOFT:
	def bwa_fq(self, fq1, outsam, bwa_path, bwa_aln_parameter, bwa_sam_parameter, RGinfo, samtools_path):
		print("[ %s ] fq alignment starts: %s \n" % (time.asctime(), fq1))
		fqprefix = (os.path.basename(fq1)).replace(".gz", "")
		outdir = os.path.dirname(outsam)
		shelldir = os.path.join(outdir, "shell")
		check_info(shelldir, "dir")
		shell = os.path.join(shelldir, fqprefix + ".aln.sh")
		logfile = shell + ".log"
		oshell = open(shell, 'w')
		if RGinfo != 0:
			RGinfo = '"' + RGinfo + '"'
			shell_line = " ".join([bwa_path, 'mem -p -C', '-R', RGinfo, ref_path, fq1, ">", outsam, '2>', logfile + "\n"])
		else:
			shell_line = " ".join([bwa_path, 'mem -p -C', ref_path, fq1, ">", outsam, '2>', logfile + "\n"])
		oshell.write(shell_line)
		outbam = outsam.replace(".sam", ".bam")
		shell_line = " ".join([samtools_path, "view -S -h -F 256 -F 2048 -b", outsam, ">", outbam]) + "\n"
		oshell.write(shell_line)
		shell_line = "rm " + outsam + "\n"
		oshell.write(shell_line)
		sv = find_samtools_version(samtools_path, shelldir)
		if sv == 0:
			shell_line = " ".join([samtools_path, "sort -m 1G", outbam, outbam]) + "\n"
		else:
			shell_line = " ".join([samtools_path, "sort -m 1G -o", outbam + ".bam", outbam]) + "\n"
		oshell.write(shell_line)
		shell_line = "mv " + outbam + ".bam " + outbam + "\n"
		oshell.write(shell_line)
		oshell.close()

		subprocess.call(["sh", shell])
		print("[ %s ] fq alignment finish: %s \n" % (time.asctime(), outbam))

def check_info(result, attribute):
	if attribute == "file":
		if os.path.isfile(result):
			pass
		else:
			sys.stderr.write("[ %s ] %s does not exist!\n" % (time.asctime(), result))
			sys.exit(-1)
	elif attribute == "dir":
		if os.path.isdir(result):
			pass
		else:
			os.makedirs(result)

def usage():
	alignment_usage = \
	'''
	modify the barcode sequence, purify and split merged fastq file into two files
	Version: 1.0.0
	Dependents: Python (>=3.0), BWA, SAMtools
	Last Updated Date: 2017-11-14

	Usage: python clean_fq_alignment.py <options>

	Basic options:
		-i --input, the input path of compressed or uncompressed fastqs; would be treated as the 1st file of paired fastq files if -2 also in use.
		-o --outputdir, the path of output
		-c --config, the path of configuration file [default: ./config/Basic.config]
		-r --RG, RG info for alignment
		-h --help, print help info
	Advanced options:
		-n --process_num, the max amout of sub-processing [default: 4]

	'''
	print(alignment_usage)

if __name__ == '__main__':
	if len(sys.argv) < 5:
		usage()
		sys.exit(-1)

	RGInfo = 0
	inputfqlist = None
	outputdir = None
	ConfigFile = None
	process_num = 4
	opts, args = getopt.gnu_getopt(sys.argv[1:], 'i:o:c:n:r:h', ['input', 'outputdir', 'config', 'process_num', 'RG', 'help'])
	for o, a in opts:
		if o == '-i' or o == '--input':
			inputfqlist = a
		if o == '-o' or o == '--outputdir':
			outputdir = a
		if o == '-c' or o == '--config':
			ConfigFile = a
		if o == '-n' or o == '--process_num':
			process_num = int(a)
		if o == '-r' or o == '--RG':
			RGInfo = a
		if o == '-h' or o == '--help':
			usage()
			sys.exit(-1)

	if os.path.isfile(inputfqlist):
		pass
	else:
		sys.stderr.write("[ %s ] ERROR: %s does not exist!\n" % (time.asctime(), inputfq))
		sys.exit(-1)

	if ConfigFile == None or os.path.exists(ConfigFile) == False:
		sys.stderr.write("configuration file has not been provided or does not exist, Please create it using 'python LRTK-SEQ.py Config'\n")
		sys.exit(-1)

	O = OUTERSOFT()
	G = baseinfo()
	G.get_config(ConfigFile)

	ref_path = G.Ref()
	ref_path_index = ref_path + ".ann"
	if os.path.isfile(ref_path_index):
		pass
	else:
		sys.stderr.write("[ %s ] ERROR: index of %s does not exist!\n" % (time.asctime(), ref_path))
		sys.exit(-1)
	bwa_path = G.Bwa()
	bwa_aln_parameter = G.Fq_aln_par()
	bwa_sam_parameter = G.Fq_sam_par()
	samtools = G.Samtools()
	picard = G.Picard()
	picard_merge_parameter = G.Picard_merge_parameter()

	otherfiletmplist = list()
	LaneID = os.path.basename(os.path.dirname(inputfqlist))
	SplitSamFileList = outputdir + '/split_bam.' + LaneID + '.txt'
	aln_sign = SplitSamFileList + ".sign"
	LibraryId = os.path.basename(outputdir)

	tmpfqlist = open(inputfqlist, 'r').readlines()
	if process_num > len(tmpfqlist):
		process_num = len(tmpfqlist)
	if os.path.exists(aln_sign) == False or os.path.getsize(aln_sign) == 0:
		wSplitSamFileList = open(SplitSamFileList, 'w')
		im = 0
		pool = multiprocessing.Pool(processes = process_num)
		SplitCleanFqFileList = inputfqlist
		rSplitCleanFqFileList = open(SplitCleanFqFileList, 'r')
		for SplitCleanFqFileInfo in rSplitCleanFqFileList:
			SplitCleanFqFileInfolist = re.split("\t", SplitCleanFqFileInfo.strip())
			if os.path.exists(SplitCleanFqFileInfolist[0]) == False or os.path.getsize(SplitCleanFqFileInfolist[0]) == 0:
				sys.stderr.write("[ %s ] ERROR: %s dose not exist or is empty!\n" % (time.asctime(), SplitCleanFqFileInfolist[0]))
				sys.exit(-1)
			fqsam = outputdir + "/" + (os.path.basename(SplitCleanFqFileInfolist[0])).replace("fq.gz", "sam")
#			fqsam = (str(SplitCleanFqFileInfolist[0]))[0:(len(SplitCleanFqFileInfolist[0])-6)] + ".sam"
			fqbam = fqsam.replace(".sam", ".bam")
			otherfiletmplist.append(fqsam)
			SplitSamFileinfo = fqbam + "\n"

			wSplitSamFileList.write(SplitSamFileinfo)

#		O.bwa_fq(SplitCleanFqFileInfolist[0], fqsam, bwa_path, bwa_aln_parameter, bwa_sam_parameter, RGInfo, samtools)
			im += 1
			if im <= process_num:
				pool.apply_async(O.bwa_fq, (SplitCleanFqFileInfolist[0], fqsam, bwa_path, bwa_aln_parameter, bwa_sam_parameter, RGInfo, samtools, ))
			if im == process_num:
				pool.close()
				pool.join()

				im = 0
				pool = multiprocessing.Pool(processes = process_num)
	
		if im > 0:
			pool.close()
			pool.join()
		wSplitSamFileList.close()

		waln_sign = open(aln_sign, "w")
		waln_sign.write("done!\n")
		waln_sign.close()
	
	SplitBamFilelist = list()
	otherfiletmplist = list()
	rSplitSamFileList = open(SplitSamFileList, 'r')
	for SplitSamFileListInfo in rSplitSamFileList:
		SplitBamFilelist.append(SplitSamFileListInfo.strip())
		otherfiletmplist.append(SplitSamFileListInfo.strip())
	allSplitBamFile = " I=".join(SplitBamFilelist)
	prefix_name = outputdir + "/" + LaneID
	tmpmergeshell = prefix_name + ".merge.sh"
	wtmpmergeshell = open(tmpmergeshell, 'w')
	MergeBamFile = prefix_name + ".sorted.bam"
	if os.path.isfile(MergeBamFile):
		subprocess.call(["rm", MergeBamFile])
	if len(SplitBamFilelist) > 1:
		shell_line = " ".join(["java -jar", picard, "MergeSamFiles I=" + allSplitBamFile, "O=" +  MergeBamFile, picard_merge_parameter]) + "\n"
	else:
		shell_line = " ".join(["cp", SplitBamFilelist[0], MergeBamFile]) + "\n"
	wtmpmergeshell.write(shell_line)
	shell_line = " ".join([samtools, "index", MergeBamFile + "\n"])
	wtmpmergeshell.write(shell_line)
	wtmpmergeshell.close()
	subprocess.call(["sh", tmpmergeshell])
	otherfiletmplist.append(tmpmergeshell)

	tmpdatadir = outputdir + "/tmp"
	check_info(tmpdatadir, "dir")
	for tmpfile in otherfiletmplist:
		subprocess.call(["mv", tmpfile, tmpdatadir])
	subprocess.call(["mv", SplitSamFileList, tmpdatadir])
	#subprocess.call(["mv", aln_sign, tmpdatadir])
#### encode
