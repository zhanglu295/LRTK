import os, sys, gzip
import getopt
import time
import re
import subprocess
import random
import string
import glob
from collections import defaultdict

def LRTK_usage(_Command1_, _Command2_ = "0"):
	helpinfo = dict()
	subhelpinfo = dict()

	usage_all = \
	'''

	LRTK: Linked Reads ToolKit, a toolkit for 10X genomic data analysis, including Basic data preparation, Resequencing analysis and De novo assembly
	Version: 1.0.0
	Dependents: Python (>=3.0), BWA, Picard (>=2.9), java (>=1.8), SAMtools, GATK (>= 3.0)
	Last Updated Date: 2017-06-01
	Contact: meijp@foxmail.com

	Usage: python LRTK.py <command> [options]

	Command:        
                        Config     Generate configuration file
                        Basicall   Execute the whole pipeline of basic data preparation
                        Reseqall   Execute the whole pipeline of resequencing
                        Denovoall  Execute the whole pipeline of de novo assembly
                        Clean      delete temporary files
                                                
                        Basic      Execute selected steps for basic data preparation 
                        Reseq      Execute selected steps for resequencing
                        Denovo     Execute selected steps for de novo assembly

	
	Note: LRTK is allowed to modify configure file to include customerized parameters and datasets.

	'''
	helpinfo["N"] = usage_all

	Config_options = \
	'''

	-o --outputdir, the path of output directory
	-s --softwarepath, the path of directory where all software were installed [default: dirname(LRTK.py)/bin]
	-d --datasetpath, dataset directory [default: dirname(LRTK.py)/dataset]
	-b --bed, target region of WES/target sequencing data, or non-N region of the whole genome [chr	strt	bed]

	'''
	helpinfo["Config"] = Config_options

	Clean_options = \
	'''

	-i --input, the list of sample, the first column must be sample id
	-o --outputdir, the path of the output directory used in LRTK
	-D --delete, list and delete all temporary files [default: only list not delete]

	'''
	helpinfo["Clean"] = Clean_options

	Basicall_options = \
	'''

	requisite command:
	-i --input, the input file containing fastq information (The input file contains three columns:1.Sample ID;2.Library ID;3. Path to sample fastqs).
	-o --outputdir, the path to output

	alternative command:
	-m --molecule, number. molecule length less than it would be discarded [default: 500]
	-p --parallel, the number of fq pairs that processing parallel. The max amount of invoking CPU would be 4*(-p) [default: 1]
	-N --noBX, generate additional fq file that has BX info or not [default: yes]
	-c --config, configuration file [default: outdir/config/Basic.config]
	-s --softwarepath, the path of directory where all software were installed [default: dirname(LRTK.py)/bin]
	-d --datasetpath, dataset directory [default: dirname(LRTK.py)/dataset]
	-b --bed, target region of WES/target sequencing data, or non-N region of the whole genome [chr strt    bed]

	'''
	helpinfo["Basicall"] = Basicall_options

	Basic_options = \
	'''

	CFQ_ALN generate clean fastq files and correct barcode error
 	MARK    merge all bam files belong to the same library of each sample, and barcode aware PCR duplication removal (must complete ALN)
 	BQSR    recalibrate base quality scores, using GATK.
 	STAT    calculate QC statistics, including Cf, Cr, MuFL, NFP etc. (must complete CFQ and ALN)
	MERGE   merge all bam files belong to the same sample

	'''
	helpinfo["Basic"] = Basic_options

	CFQ_ALN_options = \
	'''

	-i --input, the input file containing fastq information (The input file contains three columns:1.Sample ID;2.Library ID;3. Path to sample fastqs).
	-o --outputdir, the path to output
	-p --parallel, the number of fq pairs that processing parallel. The max amount of invoking CPU would be 4*(-p) [default: 1]
	-N --noBX, generate additional fq file that has BX info or not [default: yes]
	-c --config, configuration file [default: outdir/config/Basic.config]
	-s --softwarepath, the path of directory where all software were installed [default: dirname(LRTK.py)/bin]
	-d --datasetpath, dataset directory [default: dirname(LRTK.py)/dataset]
	-b --bed, target region of WES/target sequencing data, or non-N region of the whole genome [chr strt    bed]

	'''
	subhelpinfo["CFQ"] = CFQ_ALN_options

	MAK_options = \
	'''

	-i --input, the input file containing the information of bam files (The input file contains two columns:1.Sample Id;2.Library Id;3.Path to bam)
	-o --outputdir, the path to output
	-p --parallele, the number of fq pairs that processing parallel. The max amount of invoking CPU would be 4*(-p) [default: 1]
	-c --config, configuration file [default: outdir/config/Basic.config]
	-s --softwarepath, the path of directory where all software were installed [default: dirname(LRTK.py)/bin]
	-d --datasetpath, dataset directory [default: dirname(LRTK.py)/dataset]
	-b --bed, target region of WES/target sequencing data, or non-N region of the whole genome [chr strt    bed]

	'''
	subhelpinfo["MARK"] = MAK_options

	BQSR_options = \
	'''

	-i --input, the input file that contains the SAM/BAM files generated by STAT (The input file contains three columns:1.Sample Id;2.Library Id;3.Path to bam)
	-o --outputdir, the path to output
	-p --parallele, the number of fq pairs that processing parallel. The max amount of invoking CPU would be 4*(-p) [default: 1]
	-c --config, configuration file [default: outdir/config/Basic.config]
	-s --softwarepath, the path of directory where all software were installed [default: dirname(LRTK.py)/bin]
	-d --datasetpath, dataset directory [default: dirname(LRTK.py)/dataset]
	-b --bed, target region of WES/target sequencing data, or non-N region of the whole genome [chr strt    bed]

	'''
	subhelpinfo["BQSR"] = BQSR_options

	STAT_options = \
	'''

	-i --input, the input file that contains the SAM/BAM files generated by MAK (The input file contains two columns:1.Sample Id;2.Library Id;3.Path to bam)
	-o --outputdir, the path to output
	-p --parallele, the number of fq pairs that processing parallel. The max amount of invoking CPU would be 4*(-p) [default: 1]
	-m --molecule, number. molecule length less than it would be discarded [default: 500]
	-c --config, configuration file [default: outdir/config/Basic.config]
	-s --softwarepath, the path of directory where all software were installed [default: dirname(LRTK.py)/bin]
	-d --datasetpath, dataset directory [default: dirname(LRTK.py)/dataset]
	-b --bed, target region of WES/target sequencing data, or non-N region of the whole genome [chr strt    bed]

	'''
	subhelpinfo["STAT"] = STAT_options

	MERGE_options = \
	'''

	-i --input, the input file that contains the SAM/BAM files generated by MAK (The input file contains two columns:1.Sample Id;2.Library Id;3.Path to bam)
	-o --outputdir, the path to output
	-c --config, configuration file [default: outdir/config/Basic.config]
	-s --softwarepath, the path of directory where all software were installed [default: dirname(LRTK.py)/bin]
	-d --datasetpath, dataset directory [default: dirname(LRTK.py)/dataset]
	-b --bed, target region of WES/target sequencing data, or non-N region of the whole genome [chr strt    bed]

	'''
	subhelpinfo["MERGE"] = MERGE_options

	Reseqall_options = \
	'''

	-i --input, the input file that contains the BAM files generated by Basicall or ALN
	-o --outputdir, the path to output
	-L --chrlist, list of chromosome [e.g. chr1, chr2, chrX, default: outdir/config/chrlist.txt]
	-p --parallel, the number of CPU allowed to use [default: 1]
	-c --config, configuration file [default: outdir/config/Reseq.config]
	-s --softwarepath, the path of directory where all software were installed [default: dirname(LRTK.py)/bin]
	-d --datasetpath, dataset directory [default: dirname(LRTK.py)/dataset]
	-b --bed, target region of WES/target sequencing data, or non-N region of the whole genome [chr strt    bed]

	'''
	helpinfo["Reseqall"] = Reseqall_options

	Reseq_options = \
	'''
	Varcall     call SNVs and Indels by GATK
	SVcall      call structure variantion by GROC-SVs
	Phasing     phasing variants by HapCUT2
	'''
	helpinfo["Reseq"] = Reseq_options

	Varcall_options = \
	'''

	-i --input, the input file that contains the BAM files generated by Basicall or ALN
	-o --outputdir, the output directory path
	-L --chrlist, list of chromosomes [e.g.chr1,chr2,chrX, default: outdir/config/chrlist.txt]
	-p --parallel, the number of CPU allowed to use [default: 1]
	-c --config, configuration file [default: outdir/config/Reseq.config]
	-s --softwarepath, the path of directory where all software were installed [default: dirname(LRTK.py)/bin]
	-d --datasetpath, dataset directory [default: dirname(LRTK.py)/dataset]

	'''
	subhelpinfo["Varcall"] = Varcall_options

	SVcall_options = \
	'''

	-i --input, the input path that contains the BAM files generated by Basicall or ALN
	-o --outputdir, the path to output
	-c --config, configuration file [default: outdir/config/Reseq.config]
	-s --softwarepath, the path of directory where all software were installed [default: dirname(LRTK.py)/bin]
	-d --datasetpath, dataset directory [default: dirname(LRTK.py)/dataset]

	'''
	subhelpinfo["SVcall"] = SVcall_options

	Phasing_options = \
	'''

	-i --input, the input path that contains the BAM files generated by Basicall or ALN
	-v --vcf, unphased vcf file generated by Varcall or the other variant callers, both compressed or uncompressed vcf files are allowed
	-o --outputdir, the path to output
	-c --config, configuration file [default: outdir/config/Reseq.config]
	-s --softwarepath, the path of directory where all software were installed [default: dirname(LRTK.py)/bin]
	-d --datasetpath, dataset directory [default: dirname(LRTK.py)/dataset]

	'''
	subhelpinfo["Phasing"] = Phasing_options

	Denovoall_options = \
	'''

	-i --input, the input path
	-o --outputdir, the path to output
	-c --config, configuration file [default: outdir/config/Denovo.config]
	-s --softwarepath, the path of directory where all software were installed [default: dirname(LRTK.py)/bin]
	-d --datasetpath, dataset directory [default: dirname(LRTK.py)/dataset]

	'''
	helpinfo["Denovoall"] = Denovoall_options

	Denovo_options = \
	'''

	'''

	if _Command1_ in helpinfo:
		if _Command2_ in subhelpinfo:
			print (subhelpinfo[_Command2_])
		else:
			print (helpinfo[_Command1_])
	else:
		sys.stderr.write("\n\n\t### Unused command line option: %s\n" % _Command1_)
		print (helpinfo["N"])
	sys.exit(1)

def run_parallel(shell_file, maxnum):
	rshell_file = open(shell_file, 'r')
	shell_list = list()
	for shf in rshell_file:
		shf = shf.strip()
		shell_list.append(shf)
	finishNum = 0
	pn = 0
	ts = 0
	Shell_Num = len(shell_list)
	maxnum = int(maxnum)
	shell_line = None
	shelldir = os.path.dirname(shell_file) + "/tmp"
	if os.path.isdir(shelldir):
		pass
	else:
		os.mkdir(shelldir)
	while finishNum < Shell_Num:
		if pn < maxnum:
			if pn == 0:
				tmpshell = os.path.join(shelldir, "tmp." + str(ts) + ".sh")
				wtmpshell = open(tmpshell, 'w')
				shell_line = "sh " + shell_list[finishNum] + " &\n"
			else:
				shell_line = shell_line + "sh " + shell_list[finishNum] + " &\n"
			pn = pn + 1

		finishNum = finishNum + 1
		if pn == maxnum:
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
	elif attribute == "num":
		if re.search(r'\D', str(result)):
			sys.stderr.write("Error: string was found for parallel number, only number is accepted for -p %s \n" % result)
			sys.exit(-1)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		LRTK_usage("N")
	elif len(sys.argv) == 2:
		LRTK_usage(sys.argv[1])
	elif len(sys.argv) == 3:
		LRTK_usage(sys.argv[1], sys.argv[2])

	InputFqList = None
	OutputDir = None
	ParalleleNum = 1

	FqBamList = None
	Basic_config = None
	Reseq_config = None

	SoftwarePathDir = os.path.dirname(os.path.abspath(sys.argv[0])) + '/bin'
	DatasetPahtDir = os.path.dirname(os.path.abspath(sys.argv[0])) + '/dataset'
	nonN_region_bed = DatasetPahtDir + "/genome/nonN.bed"
###################################################### Config ##############################################################
	runConfig = 0
	if sys.argv[1] == "Config":
		runConfig = 1
		opts, args = getopt.gnu_getopt(sys.argv[runConfig:], 'o:s:d:b:', ['outputdir', 'softwarepath', 'datasetpath', 'bed'])
		for o, a in opts:
			if o == '-o' or o == '--outputdir':
				OutputDir = str(a)
				if OutputDir.endswith("/"):
					OutputDir = OutputDir[0:(len(OutputDir)-1)]
			if o == '-s' or o == '--softwarepath':
				SoftwarePathDir = a
			if o == '-d' or o == '--datasetpath':
				DatasetPahtDir = a
			if o == '-b' or o == '--bed':
				nonN_region_bed = a

		config_dir = OutputDir + "/config"
		check_info(config_dir, "dir")
		check_info(SoftwarePathDir, "dir")
		check_info(nonN_region_bed, "file")

		Basic_config = config_dir + "/Basic.config"
		Reseq_config = config_dir + "/Reseq.config"
		ScriptDir = os.path.abspath(os.path.dirname(sys.argv[0]))
		Create_config_script = ScriptDir + "/src/create_config.py"
		ctmpshell = config_dir + "/tmp.sh"
		wctmpshell = open(ctmpshell, 'w')
		shell_line = " ".join(["python", Create_config_script, "Basic", "-o", config_dir, "-s", SoftwarePathDir, "-d", DatasetPahtDir, "-b", nonN_region_bed + "\n"])
		wctmpshell.write(shell_line)
		shell_line = " ".join(["python", Create_config_script, "Reseq", "-o", config_dir, "-s", SoftwarePathDir, "-d", DatasetPahtDir, "-b", nonN_region_bed + "\n"])
		wctmpshell.write(shell_line)
		wctmpshell.close()
		subprocess.call(["sh", ctmpshell])
		subprocess.call(["rm", ctmpshell])
#		subprocess.call(["python", Create_config_script, "Basic", "-o", config_dir, "-s", SoftwarePathDir, "-d", DatasetPahtDir])
#		subprocess.call(["python", Create_config_script, "Reseq", "-o", config_dir, "-s", SoftwarePathDir, "-d", DatasetPahtDir])
		sys.stderr.write("\n\tconfiguration files for Basic, Reseq have been generated: %s, %s\n\n" % (Basic_config, Reseq_config))

####################################################### Config ##############################################################

####################################################### Clean  ##############################################################
	runClean = 0
	delete_tmp = 0
	Original_sample = None
	sampledict = dict()
	if sys.argv[1] == "Clean":
		runClean = 1
		opts, args = getopt.gnu_getopt(sys.argv[runConfig:], 'i:o:D', ['input' ,'outputdir', 'delete'])
		for o, a in opts:
			if o == '-i' or o == '--input':
				Original_sample = a
			if o == '-o' or o == '--outputdir':
				OutputDir = str(a)
				if OutputDir.endswith("/"):
					OutputDir = OutputDir[0:(len(OutputDir)-1)]
			if o == '-D' or o == '--delete':
				delete_tmp = 1
		check_info(Original_sample, "file")
		rOriginal_sample = open(Original_sample, 'r')
		randomstring = "CLEAN_" + ''.join(random.sample(string.ascii_letters + string.digits, 8))
		tmpshelldir = OutputDir + "/tmp"
		check_info(tmpshelldir, "dir")
		CLEANfile = tmpshelldir + "/" + randomstring + ".tmpfile"
		wCLEANfile = open(CLEANfile, 'w')
		for sample in rOriginal_sample:
			sampleinfo = re.split('\t', sample.strip())
			if sampleinfo[0] not in sampledict:
				sampledir = os.path.join(OutputDir, sampleinfo[0])
				for root, dirs, files in os.walk(sampledir):
					for name in files:
						mn = os.path.join(root, name)
						if re.search("tmp", mn):
							if re.search("SamByChr", mn):
								pass
							else:
								shell_line = mn + "\n"
								wCLEANfile.write(shell_line)
					for name in dirs:
						mn = os.path.join(root, name)
						if re.search("tmp", mn):
							if mn.endswith("tmp"):
								pass
							else:
								shell_line = mn + "\n"
								wCLEANfile.write(shell_line)
				sampledict[sampleinfo[0]] = sampleinfo[0]
		wCLEANfile.close()
		rOriginal_sample.close()

		if delete_tmp == 1:
			rCLEANfile = open(CLEANfile, 'r')
			for deletefile in rCLEANfile:
				subprocess.call(["rm -rf", deletefile.strip()])
				sys.stderr.write("delete %s ...\n" % deletefile.strip())
		else:
			sys.stderr.write("\n\tAll temporary files has been listed in %s, you can delete files that listed in it yourself!\n\n\tawk '{print \"rm -rf\", $1}' %s | sh -x\n\n" % (CLEANfile, CLEANfile))
####################################################### Clean  ##############################################################

####################################################### CFQ_ALN #################################################################
	runCFQ = 0
	addBX = 1
	FqBamList = None
	if sys.argv[1] == "Basicall":
		runCFQ = 1
	elif sys.argv[1] == "Basic" and sys.argv[2] == "CFQ":
		runCFQ = 2
	if runCFQ > 0:
		opts, args = getopt.gnu_getopt(sys.argv[runCFQ:], 'i:o:p:c:s:d:b:N', ['input', 'outputdir', 'parallel', 'config', 'softwarepath', 'datasetpath', 'bed', 'noBX'])
		for o, a in opts:
			if o == '-i' or o == '--input':
				InputFqList = a
			if o == '-o' or o == '--outputdir':
				OutputDir = str(a)
				if OutputDir.endswith("/"):
					OutputDir = OutputDir[0:(len(OutputDir)-1)]
			if o == '-p' or o == '--parallel':
				ParalleleNum = a
			if o == '-c' or o == '--config':
				Basic_config = a
			if o == '-N' or o == '--noBX':
				addBX = 0
			if o == '-s' or o == '--softwarepath':
				SoftwarePathDir = a
			if o == '-d' or o == '--datasetpath':
				DatasetPahtDir = a
			if o == '-b' or o == '--bed':
				nonN_region_bed = a

		check_info(InputFqList, "file")
		check_info(OutputDir, "dir")
		check_info(int(ParalleleNum), "num")
		if int(ParalleleNum) > 1:
			sys.stderr.write("The maximum number of %s CPUs would be invoked at the same time\n" % ParalleleNum)

		ScriptDir = os.path.abspath(os.path.dirname(sys.argv[0]))
		print(ScriptDir)
		CFQ_script = ScriptDir + "/src/clean_fq_alignment.py"
		if os.path.isfile(CFQ_script):
			pass
		else:
			sys.stderr.write("%s does not exist, the software package might not been downloaded perfectly!" % CFQ_script)
			sys.exit(-1)
		if Basic_config != None and os.path.isfile(Basic_config):
			pass
		else:
			config_dir = OutputDir + "/config"
			check_info(config_dir, "dir")
			Basic_config = config_dir + "/Basic.config"
			if os.path.isfile(Basic_config):
				pass
			else:
				sys.stderr.write("configuration file has not been provided or does not exist, it would be generated automatically: %s" % Basic_config)
				Create_config_script = ScriptDir + "/src/create_config.py"
				subprocess.call(["python", Create_config_script, "Basic", "-o", config_dir, "-s", SoftwarePathDir, "-d", DatasetPahtDir, "-b", nonN_region_bed])

		rInputFqList = open(InputFqList, 'r')
		Result_List_Dir = OutputDir + "/Result_list"
		check_info(Result_List_Dir, 'dir')
		FqBamList = Result_List_Dir + "/Basic_MARK_input.txt"
		Basic_CFQ_Result_List = Result_List_Dir + "/Basic_CFQ_ALN_result.txt"
		wFqBamList = open(FqBamList, 'w')
		wBasic_CFQ_Result_List = open(Basic_CFQ_Result_List, 'w')
		randomstring = "CFQ_ALN" + ''.join(random.sample(string.ascii_letters + string.digits, 8))
		tmpshelldir = OutputDir + "/tmp"
		check_info(tmpshelldir, "dir")
		CFQshell = tmpshelldir + "/" + randomstring + ".sh"
		wCFQshell = open(CFQshell, 'w')
		for fqinfo in rInputFqList:
			fqinfo = fqinfo.strip()
			(SampleId, LibraryId, FqPath) = re.split("\t", fqinfo)
			FqPathBasename = os.path.basename(FqPath)
			b = str(FqPathBasename)
			if b[(len(b)-9):] == ".fastq.gz":
				b = b[0:(len(b)-9)]
			elif b[(len(b)-6):] == ".fq.gz":
				b = b[0:(len(b)-6)]
			elif b[(len(b)-6):] == ".fastq":
				b = b[0:(len(b)-6)]
			elif b[(len(b)-3):] == ".fq":
				b = b[0:(len(b)-3)]
			FqPathBasename = b

			clean_FQ_output_dir = OutputDir + "/" + SampleId + "/" + LibraryId + "/" + FqPathBasename
			check_info(clean_FQ_output_dir, "dir")

			shelldir = clean_FQ_output_dir + "/shell"
			check_info(shelldir, "dir")
			FQshell = shelldir + "/ALN." + SampleId + ".sh"
			runFQshell = FQshell + "\n"
			wCFQshell.write(runFQshell)
			wFQshell = open(FQshell, 'w')
			RGinfo = "'@RG\\tID:" + LibraryId + "\\tPL:illumina\\tPU:" + FqPathBasename + "\\tLB:" + LibraryId + "\\tSM:" + SampleId + "'"

			if addBX == 0:
				shell_line = " ".join(["python", CFQ_script, "-i", FqPath, "-o", clean_FQ_output_dir, "-c", Basic_config, "-r", RGinfo, "-N", "\n"])
			else:
				shell_line = " ".join(["python", CFQ_script, "-i", FqPath, "-o", clean_FQ_output_dir, "-c", Basic_config, "-r", RGinfo, "\n"])
			wFQshell.write(shell_line)
			wFQshell.close()
			bam = clean_FQ_output_dir + "/" + FqPathBasename + ".sorted.bam"
			baminfo = "\t".join([SampleId, LibraryId, bam]) + "\n"
			wFqBamList.write(baminfo)

			BX_fq_path = clean_FQ_output_dir + "/" + FqPathBasename + ".BX.fq.gz"
			fq_Summary = clean_FQ_output_dir + "/" + FqPathBasename + ".fq_statistics.txt"
			fq_failed = clean_FQ_output_dir + "/" + FqPathBasename + ".failed.fq.gz"
			FastQC_result_path = clean_FQ_output_dir + "/fastqc"
			Result_log = "\t".join([SampleId, LibraryId, "fq including barcode info:", BX_fq_path + "\n"])
			wBasic_CFQ_Result_List.write(Result_log)
			Result_log = "\t".join([SampleId, LibraryId, "reads with abnormal barcode:", fq_failed + "\n"])
			wBasic_CFQ_Result_List.write(Result_log)
			Result_log = "\t".join([SampleId, LibraryId, "data summary of fq", fq_Summary + "\n"])
			wBasic_CFQ_Result_List.write(Result_log)
			Result_log = "\t".join([SampleId, LibraryId, "fastQC result:", FastQC_result_path + "\n"])
			wBasic_CFQ_Result_List.write(Result_log)
			Result_log = "\t".join([SampleId, LibraryId, "bam file:", bam + "\n"])
			wBasic_CFQ_Result_List.write(Result_log)
		wFqBamList.close()
		wCFQshell.close()
		wBasic_CFQ_Result_List.close()

		run_parallel(CFQshell, ParalleleNum)
		sys.stderr.write("[ %s ] files generated by 'Basic CFQ_ALN' have been listed in \n\n\t %s \n\n" % (time.asctime(), Basic_CFQ_Result_List))
		sys.stderr.write("[ %s ] new fq has been listed in \n\n\t %s \n\n\tAnd it's the input file in the 'Basic MARK' step\n\n" % (time.asctime(), FqBamList))
####################################################### CFQ_ALN #################################################################

####################################################### MARK #################################################################
	runMAK = 0
	MarkBamFileList = None
	if sys.argv[1] == "Basicall":
		runMAK = 1
	elif sys.argv[1] == "Basic" and sys.argv[2] == "MARK":
		runMAK = 2
	if runMAK > 0:
		opts, args = getopt.gnu_getopt(sys.argv[runMAK:], 'i:o:p:c:s:d:b:', ['input', 'outputdir', 'parallel', 'config', 'softwarepath', 'datasetpath', 'bed'])
		for o, a in opts:
			if runMAK == 2:
				if o == '-i' or o == '--input':
					FqBamList = a
			if o == '-o' or o == '--outputdir':
				OutputDir = str(a)
				if OutputDir.endswith("/"):
					OutputDir = OutputDir[0:(len(OutputDir)-1)]
			if o == '-p' or o == '--parallel':
				ParalleleNum = a
			if o == '-c' or o == '--config':
				Basic_config = a
			if o == '-s' or o == '--softwarepath':
				SoftwarePathDir = a
			if o == '-d' or o == '--datasetpath':
				DatasetPahtDir = a
			if o == '-b' or o == '--bed':
				nonN_region_bed = a

		check_info(FqBamList, "file")
		check_info(OutputDir, "dir")
		check_info(ParalleleNum, "num")

		if int(ParalleleNum) > 1:
			sys.stderr.write("The maximum number of %s CPUs would be invoked at the same time\n" % ParalleleNum)

		ScriptDir = os.path.dirname(os.path.abspath(sys.argv[0]))
		MAK_script = ScriptDir + "/src/mergeMark_bam.py"
		if os.path.isfile(MAK_script):
			pass
		else:
			sys.stderr.write("%s does not exist, the software package might not been downloaded perfectly!" % ALN_script)
			sys.exit(-1)
		if Basic_config != None and os.path.isfile(Basic_config):
			pass
		else:
			config_dir = OutputDir + "/config"
			check_info(config_dir, "dir")
			Basic_config = config_dir + "/Basic.config"
			if os.path.isfile(Basic_config):
				pass
			else:
				Create_config_script = ScriptDir + "/src/create_config.py"
				subprocess.call(["python", Create_config_script, "Basic", "-o", config_dir, "-s", SoftwarePathDir, "-d", DatasetPahtDir, "-b", nonN_region_bed])

		bamDict = defaultdict(dict)
		rFqBamList = open(FqBamList, 'r')
		for eachBamfile in rFqBamList:
			eachBamfile = eachBamfile.strip()
			(SampleId, LibraryId, bamfile) = re.split("\t", eachBamfile)
			check_info(bamfile, 'file')
			if SampleId in bamDict:
				if LibraryId in bamDict[SampleId]:
					bamDict[SampleId][LibraryId] = bamDict[SampleId][LibraryId] + "\n" + bamfile
				else:
					bamDict[SampleId][LibraryId] = bamfile
			else:
				bamDict[SampleId][LibraryId] = bamfile
		rFqBamList.close()

		tmpshelldir = OutputDir + "/tmp"
		check_info(tmpshelldir, "dir")
		MAKshell = tmpshelldir + "/MAK_" + ''.join(random.sample(string.ascii_letters + string.digits, 8)) + ".sh"
		wMAKshell = open(MAKshell, 'w')
		Result_List_Dir = OutputDir + "/Result_list"
		check_info(Result_List_Dir, 'dir')
		Basic_MARK_Result_List = Result_List_Dir + "/Basic_MARK_result.txt"
		MarkBamFileList = Result_List_Dir + "/Basic_BQSR_input.txt"
		wBasic_MARK_Result_List = open(Basic_MARK_Result_List, 'w')
		wMarkBamFileList = open(MarkBamFileList, 'w')

		for samplekey in bamDict.keys():
			for librarykey in bamDict[samplekey].keys():
				library_bam_dir = OutputDir + "/" + samplekey + "/" + librarykey
				check_info(library_bam_dir, "dir")
				library_bam = OutputDir + "/" + samplekey + "/" + librarykey + "/bam.txt"
				wbam = open(library_bam, 'w')
				allbam = bamDict[samplekey][librarykey] + "\n"
				wbam.write(allbam)
				wbam.close()

				shelldir = OutputDir + "/" + samplekey + "/" + librarykey + "/shell"
				check_info(shelldir, "dir")
				BAMshell = shelldir + "/merge_mark_bam.sh"
				runBAMshell = BAMshell + "\n"
				wMAKshell.write(runBAMshell)

				wBAMshell = open(BAMshell, 'w')
				markedbam = library_bam_dir + "/" + librarykey + ".sorted.merged.marked.bam"
				shell_line = " ".join(["python", MAK_script, "-i", library_bam, "-o", markedbam, "-c", Basic_config, "\n"])
				wBAMshell.write(shell_line)
				wBAMshell.close()
				baminfo = "\t".join([samplekey, librarykey, markedbam]) + "\n"
				wMarkBamFileList.write(baminfo)

				markedlog = library_bam_dir + "/" + librarykey + ".sorted.merged.marked.metrics.txt"
				Result_log = "\t".join([SampleId, LibraryId, "merged & marked bam:", markedbam + "\n"])
				wBasic_MARK_Result_List.write(Result_log)
				Result_log = "\t".join([SampleId, LibraryId, "duplication info:", markedlog + "\n"])
				wBasic_MARK_Result_List.write(Result_log)
		wMAKshell.close()
		wMarkBamFileList.close()
		wBasic_MARK_Result_List.close()

		run_parallel(MAKshell, ParalleleNum)
		sys.stderr.write("[ %s ] files generated by 'Basic MARK' have been listed in \n\n\t %s \n\n" % (time.asctime(), Basic_MARK_Result_List))
		sys.stderr.write("[ %s ] merged and duplication marked bam files of each library have been listed in \n\n\t %s \n\n\tAnd it's the input file in the 'Basic BQSR' step\n\n" % (time.asctime(), MarkBamFileList))
####################################################### MARK #################################################################

####################################################### BQSR ################################################################
	BqsrBamFileList = None
	runBQSR = 0
	if sys.argv[1] == "Basicall":
		runBQSR = 1
	elif sys.argv[1] == "Basic" and sys.argv[2] == "BQSR":
		runBQSR = 2
	if runBQSR > 0:
		opts, args = getopt.gnu_getopt(sys.argv[runBQSR:], 'i:o:p:c:s:d:b:', ['input', 'outputdir', 'parallel', 'config', 'softwarepath', 'datasetpath', 'bed'])
		for o, a in opts:
			if runBQSR == 2:
				if o == '-i' or o == '--input':
					MarkBamFileList = a
			if o == '-o' or o == '--outputdir':
				OutputDir = str(a)
				if OutputDir.endswith("/"):
					OutputDir = OutputDir[0:(len(OutputDir)-1)]
			if o == '-p' or o == '--parallel':
				ParalleleNum = a
			if o == '-c' or o == '--config':
				Basic_config = a
			if o == '-s' or o == '--softwarepath':
				SoftwarePathDir = a
			if o == '-d' or o == '--datasetpath':
				DatasetPahtDir = a
			if o == '-b' or o == '--bed':
				nonN_region_bed = a

		check_info(MarkBamFileList, "file")
		check_info(OutputDir, "dir")
		check_info(ParalleleNum, "num")
		if int(ParalleleNum) > 1:
			sys.stderr.write("The maximum number of %s CPUs would be invoked at the same time\n" % ParalleleNum)

		ScriptDir = os.path.dirname(os.path.abspath(sys.argv[0]))
		BQSR_script = ScriptDir + "/src/bqsr.py"
		if os.path.isfile(BQSR_script):
			pass
		else:
			sys.stderr.write("%s does not exist, the software package might not been downloaded perfectly!" % BQSR_script)
			sys.exit(-1)
		if Basic_config != None and os.path.isfile(Basic_config):
			pass
		else:
			config_dir = OutputDir + "/config"
			check_info(config_dir, "dir")
			Basic_config = config_dir + "/Basic.config"
			if os.path.isfile(Basic_config):
				pass
			else:
				Create_config_script = ScriptDir + "/src/create_config.py"
				subprocess.call(["python", Create_config_script, "Basic", "-o", config_dir, "-s", SoftwarePathDir, "-d", DatasetPahtDir, "-b", nonN_region_bed])

		rMarkBamFileList = open(MarkBamFileList, 'r')
		tmpshelldir = OutputDir + "/tmp"
		check_info(tmpshelldir, "dir")
		BQSRshell = tmpshelldir + "/BQSR_" + ''.join(random.sample(string.ascii_letters + string.digits, 8)) + ".sh"
		wBQSRshell = open(BQSRshell, 'w')
		Result_List_Dir = OutputDir + "/Result_list"
		check_info(Result_List_Dir, 'dir')
		BqsrResultList = Result_List_Dir + "/Basic_BQSR_result.txt"
		BqsrBamFileList = Result_List_Dir + "/Basic_STAT_input.txt"
		wBqsrResultList = open(BqsrResultList, 'w')
		wBqsrBamFileList = open(BqsrBamFileList, 'w')
		for eachBamfile in rMarkBamFileList:
			eachBamfile = eachBamfile.strip()
			(SampleId, LibraryId, bamfile) = re.split("\t", eachBamfile)
			check_info(bamfile, 'file')

			bamdir = os.path.dirname(bamfile)
			shelldir = bamdir + "/shell"
			check_info(shelldir, "dir")
			Eachshell = shelldir + "/run_bqsr.sh"
			runEachshell = Eachshell + "\n"
			wBQSRshell.write(runEachshell)

			BqsrBamFile = bamfile.replace("merged.marked.bam", "merged.marked.BQSR.bam")
			wEachshell = open(Eachshell, 'w')
			shell_line = " ".join(["python", BQSR_script, "-i", bamfile, "-o", BqsrBamFile, "-c", Basic_config, "\n"])
			wEachshell.write(shell_line)
			wEachshell.close()

			BqsrBamFile_info = "\t".join([SampleId, LibraryId, BqsrBamFile]) + "\n"
			wBqsrBamFileList.write(BqsrBamFile_info)
			Result_log = "\t".join([SampleId, LibraryId, "BQSR:", BqsrBamFile + "\n"])
			wBqsrResultList.write(Result_log)
		rMarkBamFileList.close()
		wBQSRshell.close()
		wBqsrResultList.close()
		wBqsrBamFileList.close()

		run_parallel(BQSRshell, ParalleleNum)
		sys.stderr.write("[ %s ] files generated by 'Basic BQSR' have been listed in \n\n\t %s \n\n" % (time.asctime(), BqsrResultList))
		sys.stderr.write("[ %s ] bam files filtered/unfiltered based on molecule have been listed in \n\n\t %s \n\n\tAnd it's the input file in the 'Basic STAT' step\n\n" % (time.asctime(), BqsrBamFileList))
####################################################### BQSR ################################################################

####################################################### STAT ################################################################
	runSTAT = 0
	Molecule_length = str(500)
	FilteredBamFileList = None
	if sys.argv[1] == "Basicall":
		runSTAT = 1
	elif sys.argv[1] == "Basic" and sys.argv[2] == "STAT":
		runSTAT = 2
	if runSTAT > 0:
		opts, args = getopt.gnu_getopt(sys.argv[runSTAT:], 'i:o:m:p:c:s:d:b:', ['input', 'outputdir', 'molecule', 'parallel', 'config', 'softwarepath', 'datasetpath', 'bed'])
		for o, a in opts:
			if runSTAT == 2:
				if o == '-i' or o == '--input':
					BqsrBamFileList = a
			if o == '-o' or o == '--outputdir':
				OutputDir = str(a)
				if OutputDir.endswith("/"):
					OutputDir = OutputDir[0:(len(OutputDir)-1)]
			if o == '-m' or o == '--molecule':
				Molecule_length = a
			if o == '-p' or o == '--parallel':
				ParalleleNum = a
			if o == '-c' or o == '--config':
				Basic_config = a
			if o == '-s' or o == '--softwarepath':
				SoftwarePathDir = a
			if o == '-d' or o == '--datasetpath':
				DatasetPahtDir = a
			if o == '-b' or o == '--bed':
				nonN_region_bed = a

		check_info(BqsrBamFileList, "file")
		check_info(OutputDir, "dir")
		check_info(ParalleleNum, "num")
		if int(ParalleleNum) > 1:
			sys.stderr.write("The maximum number of %s CPUs would be invoked at the same time\n" % ParalleleNum)

		ScriptDir = os.path.dirname(os.path.abspath(sys.argv[0]))
		STAT_script = ScriptDir + "/src/calculate.py"
		if os.path.isfile(STAT_script):
			pass
		else:
			sys.stderr.write("%s does not exist, the software package might not been downloaded perfectly!" % STAT_script)
			sys.exit(-1)
		if Basic_config != None and os.path.isfile(Basic_config):
			pass
		else:
			config_dir = OutputDir + "/config"
			check_info(config_dir, "dir")
			Basic_config = config_dir + "/Basic.config"
			if os.path.isfile(Basic_config):
				pass
			else:
				Create_config_script = ScriptDir + "/src/create_config.py"
				subprocess.call(["python", Create_config_script, "Basic", "-o", config_dir, "-s", SoftwarePathDir, "-d", DatasetPahtDir, "-b", nonN_region_bed])

		rBqsrBamFileList = open(BqsrBamFileList, 'r')
		tmpshelldir = OutputDir + "/tmp"
		check_info(tmpshelldir, "dir")
		STATshell = tmpshelldir + "/STAT_" + ''.join(random.sample(string.ascii_letters + string.digits, 8)) + ".sh"
		wSTATshell = open(STATshell, 'w')
		Result_List_Dir = OutputDir + "/Result_list"
		check_info(Result_List_Dir, 'dir')
		CFCRFileList = Result_List_Dir + "/Basic_STAT_result.txt"
		FilteredBamFileList = Result_List_Dir + "/Basic_MERGE_input.txt"
		wFilteredBamFileList = open(FilteredBamFileList, 'w')
		wCFCRFileList = open(CFCRFileList, 'w')
		for eachBamfile in rBqsrBamFileList:
			eachBamfile = eachBamfile.strip()
			(SampleId, LibraryId, bamfile) = re.split("\t", eachBamfile)
			check_info(bamfile, 'file')

			bamdir = os.path.dirname(bamfile)
			Result_log = "\t".join([SampleId, LibraryId, "CF&CR info:", bamdir + "/statistics/CFCR.stat\n"])
			wCFCRFileList.write(Result_log)
			Result_log = "\t".join([SampleId, LibraryId, "molecule details:", bamdir + "/statistics/molecule.full.gz\n"])
			wCFCRFileList.write(Result_log)
			Result_log = "\t".join([SampleId, LibraryId, "molecule length distribution:", bamdir + "/statistics/molecule_length_distribution.txt\n"])
			wCFCRFileList.write(Result_log)
			Result_log = "\t".join([SampleId, LibraryId, "coverage of each molecule:", bamdir + "/statistics/each_molecule_coverage.txt\n"])
			wCFCRFileList.write(Result_log)
			Result_log = "\t".join([SampleId, LibraryId, "distribution of molecule coverage:", bamdir + "/statistics/molecule_coverage_distribution.txt\n"])
			wCFCRFileList.write(Result_log)
			Result_log = "\t".join([SampleId, LibraryId, "molecule amounts of each barcode:", bamdir + "/statistics/barcode_molecule_amount.txt\n"])
			wCFCRFileList.write(Result_log)
			shelldir = bamdir + "/shell"
			check_info(shelldir, "dir")
			Eachshell = shelldir + "/stat.sh"
			runEachshell = Eachshell + "\n"
			wSTATshell.write(runEachshell)

			wEachshell = open(Eachshell, 'w')
			shell_line = " ".join(["python", STAT_script, "-i", bamfile, "-o", bamdir, "-c", Basic_config, "-m", Molecule_length, "\n"])
			wEachshell.write(shell_line)
			wEachshell.close()

			FilteredBamFile = bamfile.replace("merged.marked.BQSR.bam", "merged.marked.BQSR.filtered.bam")
			FilteredBamFile_info = "\t".join([SampleId, LibraryId, FilteredBamFile]) + "\n"
			wFilteredBamFileList.write(FilteredBamFile_info)
		rBqsrBamFileList.close()
		wSTATshell.close()
		wCFCRFileList.close()
		wFilteredBamFileList.close()

		run_parallel(STATshell, ParalleleNum)
		sys.stderr.write("[ %s ] files generated by 'Basic STAT' have been listed in \n\n\t %s \n\n" % (time.asctime(), CFCRFileList))
		sys.stderr.write("[ %s ] bam files filtered/unfiltered based on molecule have been listed in \n\n\t %s \n\n\tAnd it's the input file in the 'Basic MERGE' step\n\n" % (time.asctime(), FilteredBamFileList))
####################################################### STAT ################################################################

####################################################### MERGE ###############################################################
	FinalBamFileList = None
	runMERGE = 0
	if sys.argv[1] == "Basicall":
		runMERGE = 1
	elif sys.argv[1] == "Basic" and sys.argv[2] == "MERGE":
		runMERGE = 2
	if runMERGE > 0:
		opts, args = getopt.gnu_getopt(sys.argv[runMERGE:], 'i:o:p:c:s:d:b:', ['input', 'outputdir', 'parallel', 'config', 'softwarepath', 'datasetpath', 'bed'])
		for o, a in opts:
			if runMERGE == 2:
				if o == '-i' or o == '--input':
					FilteredBamFileList = a
			if o == '-o' or o == '--outputdir':
				OutputDir = str(a)
				if OutputDir.endswith("/"):
					OutputDir = OutputDir[0:(len(OutputDir)-1)]
			if o == '-p' or o == '--parallel':
				ParalleleNum = a
			if o == '-c' or o == '--config':
				Basic_config = a
			if o == '-s' or o == '--softwarepath':
				SoftwarePathDir = a
			if o == '-d' or o == '--datasetpath':
				DatasetPahtDir = a
			if o == '-b' or o == '--bed':
				nonN_region_bed = a

		check_info(FilteredBamFileList, "file")
		check_info(OutputDir, "dir")
		check_info(ParalleleNum, "num")
		if int(ParalleleNum) > 1:
			sys.stderr.write("The maximum number of %s CPUs would be invoked at the same time\n" % ParalleleNum)

		ScriptDir = os.path.dirname(os.path.abspath(sys.argv[0]))
		MERGE_script = ScriptDir + "/src/mergebam.py"
		if os.path.isfile(MERGE_script):
			pass
		else:
			sys.stderr.write("%s does not exist, the software package might not been downloaded perfectly!" % MERGE_script)
			sys.exit(-1)
		if Basic_config != None and os.path.isfile(Basic_config):
			pass
		else:
			config_dir = OutputDir + "/config"
			check_info(config_dir, "dir")
			Basic_config = config_dir + "/Basic.config"
			if os.path.isfile(Basic_config):
				pass
			else:
				Create_config_script = ScriptDir + "/src/create_config.py"
				subprocess.call(["python", Create_config_script, "Basic", "-o", config_dir, "-s", SoftwarePathDir, "-d", DatasetPahtDir, "-b", nonN_region_bed])

		bamDict = dict()
		rFilteredBamFileList = open(FilteredBamFileList, 'r')
		for eachBamfile in rFilteredBamFileList:
			eachBamfile = eachBamfile.strip()
			(SampleId, LibraryId, bamfile) = re.split("\t", eachBamfile)
			check_info(bamfile, 'file')
			if SampleId in bamDict:
				bamDict[SampleId] = bamDict[SampleId] + "\n" + bamfile
			else:
				bamDict[SampleId] = bamfile
		rFilteredBamFileList.close()

		tmpshelldir = OutputDir + "/tmp"
		check_info(tmpshelldir, "dir")
		MERGEshell = tmpshelldir + "/MERGE_" + ''.join(random.sample(string.ascii_letters + string.digits, 8)) + ".sh"
		wMERGEshell = open(MERGEshell, 'w')
		Result_List_Dir = OutputDir + "/Result_list"
		check_info(Result_List_Dir, 'dir')
		FinalBamFileList = Result_List_Dir + "/Reseq_Varcall_Phasing_input.txt"
		Basic_MERGE_Result_List = Result_List_Dir + "/Basic_MERGE_result.txt"
		wFinalBamFileList = open(FinalBamFileList, 'w')
		wBasic_MERGE_Result_List = open(Basic_MERGE_Result_List, 'w')

		for samplekey in bamDict.keys():
			sample_bam_dir = OutputDir + "/" + samplekey
			check_info(sample_bam_dir, "dir")
			sample_bam = OutputDir + "/" + samplekey + "/bam.txt"
			wbam = open(sample_bam, 'w')
			allbam = bamDict[samplekey] + "\n"
			wbam.write(allbam)
			wbam.close()

			shelldir = OutputDir + "/" + samplekey + "/shell"
			check_info(shelldir, "dir")
			BAMshell = shelldir + "/merge_bam.sh"
			runBAMshell = BAMshell + "\n"
			wMERGEshell.write(runBAMshell)
			wBAMshell = open(BAMshell, 'w')
			mergedbam = sample_bam_dir + "/" + samplekey + ".sorted.merged.marked.BQSR.filtered.bam"
			shell_line = " ".join(["python", MERGE_script, "-i", sample_bam, "-o", mergedbam, "-c", Basic_config, "\n"])
			wBAMshell.write(shell_line)
			wBAMshell.close()

			baminfo = "\t".join([samplekey, mergedbam]) + "\n"
			wFinalBamFileList.write(baminfo)
			Result_log = "\t".join([SampleId, "final alignment result:", mergedbam + "\n"])
			wBasic_MERGE_Result_List.write(Result_log)
		wMERGEshell.close()
		wFinalBamFileList.close()
		wBasic_MERGE_Result_List.close()

		run_parallel(MERGEshell, ParalleleNum)
		sys.stderr.write("[ %s ] files generated by 'Basic MERGE' have been listed in \n\n\t %s \n\n" % (time.asctime(), Basic_MERGE_Result_List))
		sys.stderr.write("[ %s ] finally merged bam files have been listed in \n\t %s \n\n\tAnd it's the input file in the 'Reseq Varcall' and 'Reseq Phasing' step\n\n" % (time.asctime(), FinalBamFileList))

####################################################### MERGE ###############################################################

####################################################### Varcall #############################################################
	UnphasedVcfList = None
	runVAR = 0
	if sys.argv[1] == "Reseqall":
		runVAR = 1
	elif sys.argv[1] == "Reseq" and sys.argv[2] == "Varcall":
		runVAR = 2
	if runVAR > 0:

		all_command = " ".join(sys.argv)
		sys.stderr.write("Command:\t%s\n\n" % all_command)

		opts, args = getopt.gnu_getopt(sys.argv[runVAR:], 'i:o:L:p:c:s:d:b:', ['input', 'outputdir', 'chrlist', 'parallel', 'config', 'softwarepath', 'datasetpath', 'bed'])
		for o, a in opts:
			if o == '-i' or o == '--input':
				FinalBamFileList = a
			if o == '-o' or o == '--outputdir':
				OutputDir = str(a)
				if OutputDir.endswith("/"):
					OutputDir = OutputDir[0:(len(OutputDir)-1)]
			if o == '-L' or o == '--chrlist':
				Chrlist = a
			if o == '-p' or o == '--parallel':
				ParalleleNum = a
			if o == '-c' or o == '--config':
				Reseq_config = a
			if o == '-s' or o == '--softwarepath':
				SoftwarePathDir = a
			if o == '-d' or o == '--datasetpath':
				DatasetPahtDir = a
			if o == '-b' or o == '--bed':
				nonN_region_bed = a

		check_info(FinalBamFileList, "file")
		check_info(OutputDir, "dir")
		check_info(ParalleleNum, "num")

		if int(ParalleleNum) > 1:
			sys.stderr.write("The maximum number of %s CPUs would be invoked at the same time\n" % ParalleleNum)

		ScriptDir = os.path.dirname(os.path.abspath(sys.argv[0]))
		VAR_script = ScriptDir + "/src/SnpInDel_call.py"
		if os.path.isfile(VAR_script):
			pass
		else:
			sys.stderr.write("%s does not exist, the software package might not been downloaded perfectly!" % VAR_script)
			sys.exit(-1)
		if Reseq_config != None and os.path.isfile(Reseq_config):
			pass
		else:
			config_dir = OutputDir + "/config"
			check_info(config_dir, "dir")
			Reseq_config = config_dir + "/Reseq.config"
			if os.path.isfile(Reseq_config):
				pass
			else:
				Create_config_script = ScriptDir + "/src/create_config.py"
				subprocess.call(["python", Create_config_script, "Reseq", "-o", config_dir, "-s", SoftwarePathDir, "-d", DatasetPahtDir, "-b", nonN_region_bed])

		Result_List_Dir = OutputDir + "/Result_list"
		check_info(Result_List_Dir, 'dir')
		Reseq_Varcall_Result_List = Result_List_Dir + "/Reseq_Varcall_result.txt"
		UnphasedVcfList = Result_List_Dir + "/Reseq_Phasing_input.txt"
		wReseq_Varcall_Result_List = open(Reseq_Varcall_Result_List, 'w')
		wUnphasedVcfList = open(UnphasedVcfList, 'w')
		rFinalBamFileList = open(FinalBamFileList, 'r')
		tmpshelldir = OutputDir + "/tmp"
		check_info(tmpshelldir, "dir")
		VARshell = tmpshelldir + "/VAR_" + ''.join(random.sample(string.ascii_letters + string.digits, 8)) + ".sh"
		wVARshell = open(VARshell, 'w')
		for eachBamfile in rFinalBamFileList:
			eachBamfile = eachBamfile.strip()
			(SampleId, bamfile) = re.split("\t", eachBamfile)
			check_info(bamfile, 'file')

			SampleDir = OutputDir + "/" + SampleId
			check_info(SampleDir, "dir")
			shelldir = SampleDir + "/shell"
			check_info(shelldir, "dir")

			Eachshell = shelldir + "/variant_call.sh"
			runEachshell = Eachshell + "\n"
			wVARshell.write(runEachshell)

			wEachshell = open(Eachshell, 'w')
			shell_line = " ".join(["python", VAR_script, "-i", bamfile, "-o", SampleDir, "-c", Reseq_config, "-L", Chrlist, "\n"])
			wEachshell.write(shell_line)
			wEachshell.close()
			unphased_vcf = SampleDir + "/vcf/all.vcf\n"
			wUnphasedVcfList.write(unphased_vcf)

			Result_log = "\t".join([SampleId, "unphased vcf file:", unphased_vcf])
			wReseq_Varcall_Result_List.write(Result_log)
		rFinalBamFileList.close()
		wVARshell.close()
		wUnphasedVcfList.close()
		wReseq_Varcall_Result_List.close()

		run_parallel(VARshell, ParalleleNum)
		sys.stderr.write("[ %s ] files generated by 'Reseq Varcall' have been listed in \n\n\t %s \n\n" % (time.asctime(), Reseq_Varcall_Result_List))
		sys.stderr.write("[ %s ] Variation call results have been listed in \n\n\t %s \n\n\tAnd it's the input file in the 'Reseq Phasing' step\n\n" % (time.asctime(), UnphasedVcfList))
####################################################### Varcall #############################################################

####################################################### Phasing #############################################################
	runPHASE = 0
	if sys.argv[1] == "Reseqall":
		runPHASE = 1
	elif sys.argv[1] == "Reseq" and sys.argv[2] == "Phasing":
		runPHASE = 2
	if runPHASE > 0:
		opts, args = getopt.gnu_getopt(sys.argv[runPHASE:], 'i:o:p:v:c:s:d:b:', ['input', 'outputdir', 'vcf', 'parallel', 'config', 'softwarepath', 'datasetpath', 'bed'])
		for o, a in opts:
			if runPHASE == 2:
				if o == '-v' or o == '--vcf':
					UnphasedVcfList = a
			if o == '-i' or o == '--input':
				FinalBamFileList = a
			if o == '-o' or o == '--outputdir':
				OutputDir = str(a)
				if OutputDir.endswith("/"):
					OutputDir = OutputDir[0:(len(OutputDir)-1)]
			if o == '-p' or o == '--parallel':
				ParalleleNum = a
			if o == '-c' or o == '--config':
				Basic_config = a
			if o == '-s' or o == '--softwarepath':
				SoftwarePathDir = a
			if o == '-d' or o == '--datasetpath':
				DatasetPahtDir = a
			if o == '-b' or o == '--bed':
				nonN_region_bed = a

		check_info(FinalBamFileList, "file")
		check_info(UnphasedVcfList, "file")
		check_info(OutputDir, "dir")
		check_info(ParalleleNum, "num")

		if int(ParalleleNum) > 1:
			sys.stderr.write("The maximum number of %s CPUs would be invoked at the same time\n" % ParalleleNum)

		ScriptDir = os.path.dirname(os.path.abspath(sys.argv[0]))
		PHA_script = ScriptDir + "/src/phasing.py"
		if os.path.isfile(PHA_script):
			pass
		else:
			sys.stderr.write("%s does not exist, the software package might not been downloaded perfectly!" % PHA_script)
			sys.exit(-1)
		if Reseq_config != None and os.path.isfile(Reseq_config):
			pass
		else:
			config_dir = OutputDir + "/config"
			check_info(config_dir, "dir")
			Reseq_config = config_dir + "/Reseq.config"
			if os.path.isfile(Reseq_config):
				pass
			else:
				Create_config_script = ScriptDir + "/src/create_config.py"
				subprocess.call(["python", Create_config_script, "Reseq", "-o", config_dir, "-s", SoftwarePathDir, "-d", DatasetPahtDir, "-b", nonN_region_bed])

		rFinalBamFileList = open(FinalBamFileList, 'r')
		BamDict = dict()
		rUnphasedVcfList = open(UnphasedVcfList, 'r')
		VcfDict = dict()
		for eachBamfile in rFinalBamFileList:
			eachBamfile = eachBamfile.strip()
			(SampleId, bamfile) = re.split("\t", eachBamfile)
			check_info(bamfile, 'file')
			BamDict[SampleId] = bamfile
		for eachVcffile in rUnphasedVcfList:
			eachVcffile = eachVcffile.strip()
			(SampleId, vcffile) = re.split("\t", eachVcffile)
			check_info(vcffile, 'file')
			VcfDict[SampleId] = vcffile
		rFinalBamFileList.close()
		rUnphasedVcfList.close()

		tmpshelldir = OutputDir + "/tmp"
		check_info(tmpshelldir, "dir")
		PHASEshell = tmpshelldir + "/PHASE_" + ''.join(random.sample(string.ascii_letters + string.digits, 8)) + ".sh"
		wPHASEshell = open(PHASEshell, 'w')
		Result_List_Dir = OutputDir + "/Result_list"
		check_info(Result_List_Dir, 'dir')
		phaseVcfList = Result_List_Dir + "/Reseq_Phasing_result.txt"
		wphaseVcfList = open(phaseVcfList, 'w')
		for SampleId in BamDict.keys():
			if SampleId not in VcfDict:
				sys.stderr.write("vcf file not found: %s !\n" % SampleId)
			else:
				bamfile = BamDict[SampleId]
				unphasevcffile = VcfDict[SampleId]
				sampledir = OutputDir + "/" + SampleId + "/" + "phase"
				check_info(sampledir, "dir")
				shelldir = sampledir + "/shell"
				check_info(shelldir, "dir")
				Eachshell = shelldir + "/phase.sh"
				runEachshell = Eachshell + "\n"
				wPHASEshell.write(runEachshell)

				wEachshell = open(Eachshell, 'w')
				shell_line = " ".join(["python", PHA_script, "-i", bamfile, "-v", unphasevcffile, "-o", sampledir, "-c", Reseq_config, "\n"])
				wEachshell.write(shell_line)
				wEachshell.close()

				phasevcffile = None
				vcfname = os.path.basename(unphasevcffile)
				if vcfname.endswith('.gz'):
					b = str(vcfname)
					vcfname = b[0:(len(b)-3)]
				if vcfname.endswith(".vcf"):
					b = str(vcfname)
					phasevcffile = sampledir + "/" + b[0:(len(b)-4)] + ".phased.vcf"
				else:
					phasevcffile = sampledir + "/" + vcfname + ".phased.vcf"
				phasevcffile = SampleId + "\tphased vcf file:\t" + phasevcffile + "\n"
				wphaseVcfList.write(phasevcffile)
		wPHASEshell.close()
		wphaseVcfList.close()

		run_parallel(PHASEshell, ParalleleNum)
		sys.stderr.write("[ %s ] files generated by 'Reseq Phasing' have been listed in \n\n\t %s \n\n" % (time.asctime(), phaseVcfList))
####################################################### Phasing #############################################################
