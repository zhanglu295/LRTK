import os, sys, gzip
import getopt
import time
import re
import subprocess
import random
import string
import glob
from collections import defaultdict
import numpy
import pysam

def LRTK_usage(_Command1_, _Command2_ = "0"):
	helpinfo = dict()
	subhelpinfo = dict()

	usage_all = \
	'''

	LRTK-SEQ: Linked Reads ToolKit, a toolkit for 10X genomic data analysis, including Basic data preparation and Resequencing analysis
	Version: 1.0.0
	Dependents: Python (>=3.0, numpy, pysam, matplotlib), FastQC, BWA, Picard (>=2.9), java (>=1.8), SAMtools, GATK (>= 4.0), sbt, fgbio, HapCut2, NAIBR
	Last Updated Date: 2017-06-01
	Contact: meijp@foxmail.com

	Usage: python3 LRTK-SEQ.py <command> [options]

	Command:        
                        Config     Generate configuration files
                        Basicall   Run the whole pipeline of data preprocessing and alignment, including PRE, MFQ, ALN, MARK, BQSR, STAT and MERGE 
                        Reseqall   Run the whole pipeline of variant calling and phasing, including Varcall(SNV/Indel) and Phasing
                        Clean      delete temporary files
                                                
                        Basic      Run customized steps in `Basicall`
                        Reseq      Run customized steps in `Reseqall`

	Note: LRTK-SEQ is allowed to modify configure file to include customerized parameters and datasets.

	'''
	helpinfo["N"] = usage_all

	Config_options = \
	'''
	python3 LRTK-SEQ.py Config [options]

	Basic options:
	-o --outputdir, output directory.

	Advanced options:
	-s --softwarepath, string, software directory [default: {abs_path(LRTK-SEQ.py)}/bin]
	-d --datasetpath, string, dataset directory [default: {abs_path(LRTK-SEQ.py)}/dataset]
	-b --bed, string, bed of genome regions without Ns [default: {abs_path(LRTK-SEQ.py)}/dataset/GATK_bundle/nonN.bed]

	'''
	helpinfo["Config"] = Config_options

	Clean_options = \
	'''

	-i --input, string, the list of sample, the first column must be sample id
	-o --outputdir, string, the path of the output directory used in LRTK
	-D --delete, flag, list and delete all temporary files [default: only list not delete]

	'''
	helpinfo["Clean"] = Clean_options

	Basicall_options = \
	'''
	python3 LRTK-SEQ.py Basicall [0-8] [options]

	eg. python3 LRTK-SEQ.py Basicall [options]  (run all steps in Basicall step by step)
	eg. python3 LRTK-SEQ.py Basicall 1-3 [options] (run steps 1, 2, and 3 in Basicall)
	eg. python3 LRTK-SEQ.py Basicall 4- [options] (run steps 4, 5, 6, and 7 in Basicall)
	step 1: PRE     generate clean fastq files and correct barcode error
	step 2: MFQ     merge fq based on library or sample
	step 3: ALN     alignment
 	step 4: MARK    barcode-aware PCR duplicates removal
 	step 5: BQSR    base quality score recalibration
 	step 6: STAT    calculate QC statistics based on sample library, including Cf, Cr, MuFL, NFP etc.
	step 7: MERGE   merge all bam files from the same sample

	Warnings: The input file must change according to the starting steps!
	Warnings: The input file must change according to the starting steps!!
	Warnings: The input file must change according to the starting steps!!!
	For instance, if "Basicall" starts with step 3, the input file must change to be the input file of "ALN"

	Basic options:
	-i --input, string, input file for fastq information (4 or 5 columns:1.Sample ID; 2.Library ID; 3. Serial number of the library; 4.Path to the fastqs [5. Path to the fastq 2]), or others(according to the starting step).
	-o --outputdir, string, the path to output
	-c --config, string, configuration file [default: outputdir/config/Basic.config]

	Advanced options:
	-M --mergefq, int [0, 1 or 2. default: 0]. 0: do nothing, MFQ would not run; 1: merge fq files belong to the same library; 2: merge fq files belong to the same sample. If '-M' > 0, all following analysis would be based on the merged file. Other numbers might cause unpredictable error resuls.
	-m --minlen, int, the minimum length (bp) of molecule to be considered [default: 500]
	-p --parallel, int, the number CPUs can be used in parallel [default: 4]
	-z --splitsize, int, the amount of reads of splited fastq, reads num = split_size / 8 [default: 60000000 lines, 7500000 read-pairs, compressed file size: ~300M]

	'''
	helpinfo["Basicall"] = Basicall_options

	Basic_options = \
	'''
	python3 LRTK-SEQ.py Basic <command> [options]

	PRE     generate clean fastq files and correct barcode error
	MFQ     merge fq based on library or sample
	ALN     alignment
 	MARK    barcode-aware PCR duplicates removal
 	BQSR    base quality score recalibration
 	STAT    calculate QC statistics based on sample library, including Cf, Cr, MuFL, NFP etc.
	MERGE   merge all bam files from the same sample

	'''
	helpinfo["Basic"] = Basic_options

	PRE_options = \
	'''
	python3 LRTK-SEQ.py Basic PRE [options]

	Basic options:
	-i --input, string, input file for fastq information (4 or 5 columns:1.Sample ID; 2.Library ID; 3. Serial number of the library; 4.Path to the fastqs [5. Path to the fastq 2])
	-o --outputdir, string, the path to output
	-c --config, string, configuration file [default: outputdir/config/Basic.config]

	Advanced options:
	-p --parallel, int, the number CPUs can be used in parallel [default: 4] 
	-z --splitsize, int, the amount of reads of splited fastq, reads num = split_size / 8 [default: 60000000 lines, 7500000 read-pairs, compressed file size: ~300M]

	'''
	subhelpinfo["PRE"] = PRE_options

	MFQ_options = \
	'''
	python3 LRTK-SEQ.py Basic MFQ [options]

	Basic options:
	-i --input, string, input file for fastq information (4 columns:1.Sample ID; 2.Library ID; 3. fastq list; 4. prefix of new barcode, including 4 letters made up by "A", "T", "C", or "G"), generated by 'PRE'
	-o --outputdir, string, the path to output
	-c --config, string, configuration file [default: outputdir/config/Basic.config]

	Advanced option:
	-M --mergefq, int [0, 1 or 2. default: 0]. 0: do nothing, MFQ would not run; 1: merge fq files belong to the same library; 2: merge fq files belong to the same sample. If '-M' > 0, all following analysis would be based on the merged file. Other numbers might cause unpredictable error resuls.
	-p --parallele, int, the number CPUs can be used in parallel [default: 4] 

	'''
	subhelpinfo["MFQ"] = MFQ_options

	ALN_options = \
	'''
	python LRTK-SEQ.py Basic ALN [options]

	Basic options:
	-i --input, string, input file for fastq information (4 columns:1.Sample ID; 2.Library ID; 3. Serial number of the library; 4.fastqs file list generated in `PRE` or `MFQ`), generated by 'PRE' or 'MFQ'
	-o --outputdir, string, the path to output
	-c --config, string, configuration file [default: outputdir/config/Basic.config]

	Advanced option:
	-p --parallele, int, the number CPUs can be used in parallel [default: 4] 

	'''
	subhelpinfo["ALN"] = ALN_options

	MAK_options = \
	'''
	python3 LRTK-SEQ.py Basic MARK [options]

	Basic options:
	-i --input, string, input file for BAM information (3 columns: 1.Sample Id; 2.Library Id; 3.Path to BAM, generated by `ALN`)
	-o --outputdir, string, the path to output
	-c --config, string, configuration file [default: outputdir/config/Basic.config]

	'''
	subhelpinfo["MARK"] = MAK_options

	BQSR_options = \
	'''
	python3 LRTK-SEQ.py Basic BQSR [options]

	Basic options:
	-i --input, string, input file for BAM information (3 columns: 1.Sample Id; 2.Library Id; 3.Path to BAM, generated by `MARK`)
	-o --outputdir, string, the path to output
	-c --config, string, configuration file [default: outputdir/config/Basic.config]

	'''
	subhelpinfo["BQSR"] = BQSR_options

	STAT_options = \
	'''
	python3 LRTK-SEQ.py Basic STAT [options]

	Basic options
	-i --input, string, input file for BAM information (3 columns: 1.Sample Id; 2.Library Id; 3.Path to BAM, generated by `BQSR`)
	-o --outputdir, string, the path to output
	-c --config, string, configuration file [default: outputdir/config/Basic.config]

	Advanced options
	-m --minlen, int, the minimum length (bp) of molecule to be considered [default: 500]

	'''
	subhelpinfo["STAT"] = STAT_options

	MERGE_options = \
	'''
	python3 LRTK-SEQ.py Basic MERGE [options]

	Basic options:
	-i --input, string, input file for BAM information (3 columns: 1.Sample Id; 2.Library Id; 3.Path to BAM, generated by `STAT`)
	-o --outputdir, string, the path to output
	-c --config, string, configuration file [default: outputdir/config/Basic.config]

	'''
	subhelpinfo["MERGE"] = MERGE_options

	Reseqall_options = \
	'''
	python3 LRTK-SEQ.py Reseqall [options]

	Basic options:
	-i --input, string, input file for BAM information (2 columns: 1. Sample id 2. Path to BAM, generated by `Basicall`)
	-o --outputdir, string, the path to output
    -c --config, string, configuration file [default: outputdir/config/Reseq.config]
    -f --familyinfo, string, family information, including 2 columns: 1. family id, 2. sample id

	'''
	helpinfo["Reseqall"] = Reseqall_options

	Reseq_options = \
	'''
	python3 LRTK-SEQ.py Reseq <command> [options]

	Varcall     call SNVs and Indels by GATK, based on individual
	Famcall     call SNVs and InDels by GATK, based on family
	Phasing     phasing variants by HapCUT2
	'''
	helpinfo["Reseq"] = Reseq_options

	Varcall_options = \
	'''
	python3 LRTK-SEQ.py Reseq Varcall [options]

	Basic options:
	-i --input, string, input file for BAM information (2 columns: 1. Sample id 2. Path to BAM, generated by `Basicall`)
    -o --outputdir, string, the path to output
    -c --config, string, configuration file [default: outputdir/config/Reseq.config]

	'''
	subhelpinfo["Varcall"] = Varcall_options

#	SVcall_options = \
#	'''
#	python LRTK-SEQ.py Reseq SVcall [options]

#	Basic options:
#	-i --input, the BAM information files from Basicall or ALN
#	-o --outputdir, the path to output
#	-c --config, configuration file [default: ./config/Reseq.config]
#
#	Advanced optionL
#	-p --parallel, the number CPUs can be used in parallel [default: 1]
#	-m --min_mapq: Minimum mapping quality for a read to be included in analysis (default: 40)
#	-s --min_sv: Minimum size of a structural variant to be detected (default: lmax, the 95th percentile of the paired-end read insert size distribution)
#	-k --min_barcode minimum number of barcode overlaps supporting a candidate NA (default = 3)

#	'''
#	subhelpinfo["SVcall"] = SVcall_options

	Famcall_options = \
	'''
	python3 LRTK-SEQ.py Reseq Famcall [options]

	Basic options:
	-i --input, string, input file for variation calling (2 columns: 1. Sample id 2. gvcf file list), generated by `Varcall`
	-o --outputdir, string, the path to output
	-c --config, string, configuration file [default: outputdir/config/Reseq.config]
	-f --familyinfo, string, family information, including 2 columns: 1. family id, 2. sample id
	'''

	Phasing_options = \
	'''
	python3 LRTK-SEQ.py Reseq Phasing [options]

	Basic options:
	-i --input, string, the BAM information and vcf information files generated in `Basicall` and `Varcall`. (three columns: 1. Sample ID; 2. indexed bam file; 3. indexed vcf file)
	-o --outputdir, string, the path to output
	-c --config, string, configuration file [default: outputdir/config/Reseq.config]

	'''
	subhelpinfo["Phasing"] = Phasing_options

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

	SCRIPT = defaultdict(dict)

	InputFqList = None
	OutputDir = None
	ParalleleNum = 4

	Bam_for_MARK_List = None
	Basic_config = None
	Reseq_config = None
	Qsub_parameter = ""

	SoftwarePathDir = os.path.dirname(os.path.abspath(sys.argv[0])) + '/bin'
	DatasetPahtDir = os.path.dirname(os.path.abspath(sys.argv[0])) + '/dataset'
	python3 = SoftwarePathDir + "/python"
	if os.path.isfile(python3) == False:
		python3 = "python"
	nonN_region_bed = DatasetPahtDir + "/GATK_bundle/nonN.bed"
###################################################### Config ##############################################################
	runConfig = 0
	if sys.argv[1] == "Config":
		runConfig = 1
		opts, args = getopt.gnu_getopt(sys.argv[runConfig:], 'o:s:d:b:q:', ['outputdir', 'softwarepath', 'datasetpath', 'bed', 'qsub_parameter'])
		for o, a in opts:
			if o == '-o' or o == '--outputdir':
				OutputDir = os.path.abspath(a)
			if o == '-s' or o == '--softwarepath':
				SoftwarePathDir = os.path.abspath(a)
			if o == '-d' or o == '--datasetpath':
				DatasetPahtDir = os.path.abspath(a)
			if o == '-b' or o == '--bed':
				nonN_region_bed = os.path.abspath(a)
			if o == '-q' or o == '--qsub_parameter':
				Qsub_parameter = str(a)

		config_dir = OutputDir + "/config"
		check_info(config_dir, "dir")
		check_info(SoftwarePathDir, "dir")
		check_info(DatasetPahtDir, "dir")
		if os.path.isfile(nonN_region_bed):
			pass
		else:
			FindNonBed_script = os.path.abspath(os.path.dirname(sys.argv[0])) + "/src/find_nonN_region.py"
			subprocess.call([python3, FindNonBed_script, DatasetPahtDir + "/GATK_bundle/Homo_sapiens_assembly38.fasta", nonN_region_bed])
#		check_info(nonN_region_bed, "file")

		Basic_config = config_dir + "/Basic.config"
		Reseq_config = config_dir + "/Reseq.config"
		ScriptDir = os.path.abspath(os.path.dirname(sys.argv[0]))
		Create_config_script = ScriptDir + "/src/create_config.py"
		ctmpshell = config_dir + "/tmp.sh"
		wctmpshell = open(ctmpshell, 'w')
		shell_line = " ".join([python3, Create_config_script, "Basic", "-o", config_dir, "-s", SoftwarePathDir, "-d", DatasetPahtDir, "-b", nonN_region_bed + "\n"])
		wctmpshell.write(shell_line)
		shell_line = " ".join([python3, Create_config_script, "Reseq", "-o", config_dir, "-s", SoftwarePathDir, "-d", DatasetPahtDir, "-b", nonN_region_bed + "\n"])
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
				Original_sample = os.path.abspath(a)
			if o == '-o' or o == '--outputdir':
				OutputDir = os.path.abspath(a)
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

####################################################### Preprocessing #################################################################
	runPRE = 0
	addBX = 1
	Split_FQ_for_ALN_List = None
	Split_FQ_for_MFQ_List = None
	split_size = 60000000
	STEP = 0
	step_dict = dict()
	if sys.argv[1] == "Basicall":
		if (str(sys.argv[2])).startswith("-"):
			runPRE = 1
			for s in range(1, 8):
				step_dict[s] = 1
		else:
			runPRE = 2
			step_list = re.split("-", str(sys.argv[2]))
			if re.search(r'\d', step_list[1]):
				step_list[0] = int(step_list[0])
				step_list[1] = int(step_list[1]) + 1
				for s in range(step_list[0], step_list[1]):
					step_dict[s] = 1
			else:
				for s in range(int(step_list[0]), 8):
					step_dict[s] = 1
	elif sys.argv[1] == "Basic" and sys.argv[2] == "PRE":
		runPRE = 2
		step_dict[1] = 1
	if runPRE > 0 and 1 in step_dict:
		opts, args = getopt.gnu_getopt(sys.argv[runPRE:], 'i:o:p:c:z:M:q:m:', ['input', 'outputdir', 'parallel', 'config', 'splitsize', 'mergefq', 'qsub_parameter', 'minlen'])
		for o, a in opts:
			if o == '-i' or o == '--input':
				InputFqList = os.path.abspath(a)
			if o == '-o' or o == '--outputdir':
				OutputDir = os.path.abspath(a)
			if o == '-p' or o == '--parallel':
				ParalleleNum = a
			if o == '-c' or o == '--config':
				Basic_config = os.path.abspath(a)
			if o == '-z' or o == '--splitsize':
				split_size = str(a)
			if o == '-q' or o == '--qsub_parameter':
				Qsub_parameter = str(a)

		check_info(InputFqList, "file")
		check_info(OutputDir, "dir")
		check_info(int(ParalleleNum), "num")
#		if int(ParalleleNum) > 1:
#			sys.stderr.write("The maximum number of %s CPUs would be invoked at the same time\n" % ParalleleNum)

		ScriptDir = os.path.abspath(os.path.dirname(sys.argv[0]))
		print(ScriptDir)
		PRE_script = ScriptDir + "/src/clean_fq.py"
		if os.path.isfile(PRE_script):
			pass
		else:
			sys.stderr.write("%s does not exist, the software package might not been downloaded perfectly!" % PRE_script)
			sys.exit(-1)
		if Basic_config != None and os.path.isfile(Basic_config):
			pass
		else:
			config_dir = OutputDir + "/config"
			check_info(config_dir, "dir")
			Basic_config = config_dir + "/Basic.config"
			print("[ %s ] Warnings: configuration file was not provided, and %s would be used!\n" % (time.asctime(), Basic_config))
			if os.path.exists(Basic_config) == False:
				sys.stderr.write("configuration file was not provided or does not exist, Please create it using 'python LRTK-SEQ.py Config'\n")
				sys.exit(-1)

		barcode_prefix_dict = {1:"AAAA", 2:"AAAC", 3:"AAAG", 4:"AAAT", 5:"AACA", 6:"AACC", 7:"AACG", 8:"AACT", 9:"AAGA", 10:"AAGC", 11:"AAGG", 12:"AAGT",}

		rInputFqList = open(InputFqList, 'r')
		Result_List_Dir = OutputDir + "/Result_list"
		check_info(Result_List_Dir, 'dir')
		Split_FQ_for_ALN_List = Result_List_Dir + "/Basic_ALN_input.txt"
		Split_FQ_for_MFQ_List = Result_List_Dir + "/Basic_MFQ_input.txt"
		Basic_PRE_Result_List = Result_List_Dir + "/Basic_PRE_result.txt"
		wSplit_FQ_for_ALN_List = open(Split_FQ_for_ALN_List, 'w')
		wSplit_FQ_for_MFQ_List = open(Split_FQ_for_MFQ_List, 'w')
		wBasic_PRE_Result_List = open(Basic_PRE_Result_List, 'w')

		randomstring = "PRE_" + ''.join(random.sample(string.ascii_letters + string.digits, 8))
		tmpshelldir = OutputDir + "/tmp"
		check_info(tmpshelldir, "dir")
		PREshell = tmpshelldir + "/" + randomstring + ".sh"
		wPREshell = open(PREshell, 'w')
		for fqinfo in rInputFqList:
			InputFqListInfolist = re.split("\t", fqinfo.strip())
#			fqinfo = fqinfo.strip()
#			(SampleId, LibraryId, FqPath) = re.split("\t", fqinfo)
			SampleId = InputFqListInfolist[0]
			LibraryId = InputFqListInfolist[1]
			LibraryNum = int(InputFqListInfolist[2])
			FqPath = InputFqListInfolist[3]
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
			if b.endswith("_1"):
				b = b[0:(len(b)-2)]
			FqPathBasename = b

			clean_FQ_output_dir = OutputDir + "/" + SampleId + "/Clean_data/" + LibraryId + "/" + FqPathBasename
			check_info(clean_FQ_output_dir, "dir")

			shelldir = clean_FQ_output_dir + "/shell"
			check_info(shelldir, "dir")
			FQshell = shelldir + "/run_PRE." + SampleId + ".sh"
			runFQshell = FQshell + "\n"
			wPREshell.write(runFQshell)
			wFQshell = open(FQshell, 'w')

			shell_line = "set -e\n"
			wFQshell.write(shell_line)
			SampleId_LibraryId = SampleId + "_" + LibraryId
			if len(InputFqListInfolist) > 4 and os.path.exists(InputFqListInfolist[4]):
				shell_line = " ".join([python3, PRE_script, "-i", FqPath, "-2", InputFqListInfolist[4], "-o", clean_FQ_output_dir, "-c", Basic_config,  "-z", str(split_size), "-n", str(ParalleleNum), "\n"])
			else:
				shell_line = " ".join([python3, PRE_script, "-i", FqPath, "-o", clean_FQ_output_dir, "-c", Basic_config,  "-z", str(split_size), "-n", str(ParalleleNum), "\n"])

			wFQshell.write(shell_line)
			shell_line = " ".join(["echo Done! >", FQshell + ".sign"]) + "\n"
			wFQshell.write(shell_line)
			wFQshell.close()
			if "PRE" not in SCRIPT or SampleId not in SCRIPT["PRE"]:
				SCRIPT["PRE"][SampleId] = FQshell + ":6G"
			else:
				SCRIPT["PRE"][SampleId] = SCRIPT["PRE"][SampleId] + "\t" + FQshell + ":6G"

			BX_fq_path = clean_FQ_output_dir + "/" + FqPathBasename + ".barcode_id.clean.fq.gz"
			BX_a_fq_path = clean_FQ_output_dir + "/" + FqPathBasename + ".barcode_seq.clean.fq.gz"
			fq_Summary = clean_FQ_output_dir + "/" + FqPathBasename + ".fq_statistics.txt"
			fq_failed = clean_FQ_output_dir + "/" + FqPathBasename + ".failed.fq.gz"
			FastQC_result_path = clean_FQ_output_dir + "/fastqc"
			Result_log = "\t".join([SampleId, LibraryId, "fq file with barcode info in read id:", BX_fq_path + "\n"])
			wBasic_PRE_Result_List.write(Result_log)
			Result_log = "\t".join([SampleId, LibraryId, "fq file with barcode info in read sequence:", BX_a_fq_path + "\n"])
			wBasic_PRE_Result_List.write(Result_log)
			Result_log = "\t".join([SampleId, LibraryId, "reads with abnormal barcode:", fq_failed + "\n"])
			wBasic_PRE_Result_List.write(Result_log)
			Result_log = "\t".join([SampleId, LibraryId, "data summary of fq", fq_Summary + "\n"])
			wBasic_PRE_Result_List.write(Result_log)
			Result_log = "\t".join([SampleId, LibraryId, "fastQC result:", FastQC_result_path + "\n"])
			wBasic_PRE_Result_List.write(Result_log)

			Split_FQ_for_ALN_List_info = "\t".join([SampleId, LibraryId, str(LibraryNum), clean_FQ_output_dir + "/Split_final_fq.txt"]) + "\n"
			wSplit_FQ_for_ALN_List.write(Split_FQ_for_ALN_List_info)
			barcode_prefix_sequence = barcode_prefix_dict[LibraryNum]
			Split_FQ_for_MFQ_List_info = "\t".join([SampleId, LibraryId, clean_FQ_output_dir + "/Split_final_fq.txt", barcode_prefix_sequence]) + "\n"
			wSplit_FQ_for_MFQ_List.write(Split_FQ_for_MFQ_List_info)

		wPREshell.close()
		wBasic_PRE_Result_List.close()
		wSplit_FQ_for_ALN_List.close()
		wSplit_FQ_for_MFQ_List.close()

#		run_parallel(PREshell, ParalleleNum)
		sys.stderr.write("[ %s ] files generated by 'Basic PRE' have been listed in \n\n\t %s \n\n" % (time.asctime(), Basic_PRE_Result_List))
		sys.stderr.write("[ %s ] new fq has been listed in \n\n\t %s and %s \n\n\tAnd it's the input file in the 'Basic ALN' and 'Basic MFQ' step\n\n" % (time.asctime(), Split_FQ_for_ALN_List, Split_FQ_for_MFQ_List))
####################################################### Preprocessing #################################################################

####################################################### MFQ #################################################################
	runMFQ = 0
	doMFQ = 0
	step_dict = dict()
	if sys.argv[1] == "Basicall":
		if (str(sys.argv[2])).startswith("-"):
			runMFQ = 1
			for s in range(1, 8):
				step_dict[s] = 1
		else:
			runMFQ = 2
			step_list = re.split("-", str(sys.argv[2]))
			if re.search(r'\d', step_list[1]):
				step_list[0] = int(step_list[0])
				step_list[1] = int(step_list[1]) + 1
				for s in range(step_list[0], step_list[1]):
					step_dict[s] = 1
			else:
				for s in range(int(step_list[0]), 8):
					step_dict[s] = 1
	elif sys.argv[1] == "Basic" and sys.argv[2] == "MFQ":
		runMFQ = 2
		step_dict[2] = 1

	opts, args = getopt.gnu_getopt(sys.argv[runMFQ:], 'i:o:p:c:z:M:q:f:m:', ['input', 'outputdir', 'parallel', 'config', 'splitsize', 'mergefq', 'qsub_parameter', 'familyinfo', 'minlen'])
	for o, a in opts:
		if o == '-M' or o == 'mergefq':
			doMFQ = int(a)
	
	if doMFQ != 0 and 2 in step_dict:
		for o, a in opts:
			if o == '-i' or o == '--input':
				if runMFQ == 2 and (1 not in step_dict):
					Split_FQ_for_MFQ_List = os.path.abspath(a)
			if o == '-o' or o == '--outputdir':
				OutputDir = os.path.abspath(a)
			if o == '-p' or o == '--parallel':
				ParalleleNum = a
			if o == '-c' or o == '--config':
				Basic_config = os.path.abspath(a)
			if o == '-z' or o == '--splitsize':
				split_size = str(a)
			if o == '-q' or o == '--qsub_parameter':
				Qsub_parameter = str(a)

		ScriptDir = os.path.abspath(os.path.dirname(sys.argv[0]))
		MFQ_script = ScriptDir + "/src/mergeFQ.py"
		if os.path.isfile(MFQ_script) == False:
			sys.stderr.write("%s does not exist, the software package might not been downloaded perfectly!" % MFQ_script)
			sys.exit(-1)

		rSplit_FQ_for_MFQ_List = open(Split_FQ_for_MFQ_List, 'r')
		Split_FQ_for_MFQ_List_dict = defaultdict(list)
		FQ_library = defaultdict(dict)
		for Split_FQ_for_MFQ_List_info in rSplit_FQ_for_MFQ_List:
			(SampleID, LibraryID, split_FQ_List, barcode_) = re.split("\t", Split_FQ_for_MFQ_List_info.strip())
			Split_FQ_for_MFQ_List_dict[SampleID].append(Split_FQ_for_MFQ_List_info.strip())
			bxfq = os.path.dirname(split_FQ_List) + "/" + os.path.basename(os.path.dirname(split_FQ_List)) + ".barcode_id.clean.fq.gz"
			if LibraryID not in FQ_library[SampleID]:
				FQ_library[SampleID][LibraryID] = bxfq
			else:
				FQ_library[SampleID][LibraryID] = FQ_library[SampleID][LibraryID] + "\t" + bxfq
		rSplit_FQ_for_MFQ_List.close()

		Result_List_Dir = OutputDir + "/Result_list"
		check_info(Result_List_Dir, 'dir')
		Basic_MFQ_Result_List = Result_List_Dir + "/Basic_MFQ_result.txt"
		Split_FQ_for_ALN_List = Result_List_Dir + "/Basic_ALN_input.txt"

		wBasic_MFQ_Result_List = open(Basic_MFQ_Result_List, 'w')
		randomstring = "MFQ_" + ''.join(random.sample(string.ascii_letters + string.digits, 8))
		tmpshelldir = OutputDir + "/tmp"
		check_info(tmpshelldir, "dir")
		MFQshell = tmpshelldir + "/" + randomstring + ".sh"
		wMFQshell = open(MFQshell, 'w')
		Sample_Merged_all_FQ_list = list()
		if doMFQ == 1:
			for SampleID in FQ_library.keys():
				L_dict = FQ_library[SampleID]
				for LibraryID in L_dict.keys():
					fqpathlist = re.split("\t", FQ_library[SampleID][LibraryID])
					shelldir = OutputDir + "/Clean_data/" + SampleID + "/" + LibraryID + "/shell"
					check_info(shelldir, "dir")
					shell = shelldir + "/run_MFQ." + SampleID + "." + LibraryID + ".sh"
					shell_sign = shell + ".sign"
					merged_lib_fq = OutputDir + "/" + SampleID + "/Clean_data/" + LibraryID + "/" + LibraryID + ".merged.fq.gz"

					Sample_Merged_all_FQ_list.append(merged_lib_fq)
					if os.path.exists(merged_lib_fq) == False or os.path.exists(shell_sign) == False:
						wMFQshell.write(shell + "\n")
						wshell = open(shell, 'w')
						wshell.write("set -e\n")
						for i in range(0, len(fqpathlist)):
							if i == 0:
								shell_line = " ".join(["cat", fqpathlist[i] + ">" + merged_lib_fq + "\n"])
							else:
								shell_line = " ".join(["cat", fqpathlist[i] + ">>" + merged_lib_fq + "\n"])
							wshell.write(shell_line)
						shell_line = " ".join(["echo Done! >", shell + ".sign\n"])
						wshell.write(shell_line)
						wshell.close()

						if "MFQ" not in SCRIPT or SampleID not in SCRIPT["MFQ"]:
							SCRIPT["MFQ"][SampleID] = shell + ":1G"
						else:
							SCRIPT["MFQ"][SampleID] = SCRIPT["MFQ"][SampleID] + "\t" + shell + ":1G"

					Result_log = "\t".join([SampleID, LibraryID, "merged fq:", merged_lib_fq + "\n"])
					wBasic_MFQ_Result_List.write(Result_log)
		elif doMFQ == 2:
			wSplit_FQ_for_ALN_List = open(Split_FQ_for_ALN_List, 'w')
			for Sample_FQ_info in Split_FQ_for_MFQ_List_dict.keys():
				SampleID = Sample_FQ_info
				fqinfo = "\n".join(Split_FQ_for_MFQ_List_dict[SampleID])
				fqinfolist = OutputDir + "/" + SampleID + "/Clean_data/" + SampleID + ".fq.txt"
				wfqinfolist = open(fqinfolist, 'w')
				wfqinfolist.write(fqinfo + "\n")
				wfqinfolist.close()

				lib_num = len(re.split("\n", fqinfo))

				shelldir = OutputDir + "/" + SampleID + "/shell"
				check_info(shelldir, "dir")
				shell = shelldir + "/run_MFQ." + SampleID + ".sh"
				wMFQshell.write(shell + "\n")

				Sample_Merged_FQ = OutputDir + "/" + SampleID + "/Clean_data/" + SampleID + ".merged.fq.gz"
				Sample_Merged_all_FQ_list.append(Sample_Merged_FQ)
				wshell = open(shell, 'w')
				wshell.write("set -e\n")
				shell_line = " ".join([python3, MFQ_script, "-i", fqinfolist, "-o", Sample_Merged_FQ, "-c", Basic_config, "-z", str(split_size), "-n", str(ParalleleNum) + "\n"])
				wshell.write(shell_line)
				shell_line = " ".join(["echo Done! >", shell + ".sign"]) + "\n"
				wshell.write(shell_line)
				wshell.close()

				if "MFQ" not in SCRIPT or SampleID not in SCRIPT["MFQ"]:
					SCRIPT["MFQ"][SampleID] = shell + ":1G"
				else:
					SCRIPT["MFQ"][SampleID] = SCRIPT["MFQ"][SampleID] + "\t" + shell + ":1G"

				Split_FQ_for_ALN_List_info = "\t".join([SampleID, "Merged_libraries\t1", OutputDir + "/" + SampleID + "/Clean_data/Merged_libraries/Split_Modified_FQ.txt"]) + "\n"
				wSplit_FQ_for_ALN_List.write(Split_FQ_for_ALN_List_info)

				Result_log = "\t".join([SampleID, "merged fq" + "(" + str(lib_num) + " libraries):", Sample_Merged_FQ + "\n"])
				wBasic_MFQ_Result_List.write(Result_log)
			wSplit_FQ_for_ALN_List.close()

		wMFQshell.close()
		wBasic_MFQ_Result_List.close()

		#run_parallel(MFQshell, ParalleleNum)
		if doMFQ != 0:
			for Sample_Merged_all_FQ in Sample_Merged_all_FQ_list:
				if os.path.exists(Sample_Merged_all_FQ) == False or os.path.getsize(Sample_Merged_all_FQ) == 0:
					print("[ %s ] ERROR: %s did not generated sucessfully!\n\n" % (time.asctime(), Sample_Merged_all_FQ))
					sys.exit(-1)
		sys.stderr.write("[ %s ] files generated by 'Basic MFQ' have been listed in \n\n\t %s \n\n" % (time.asctime(), Basic_MFQ_Result_List))
####################################################### MFQ #################################################################

####################################################### ALN #################################################################
	runALN = 0
	Bam_for_MARK_List = None
	step_dict = dict()
	if sys.argv[1] == "Basicall":
		if (str(sys.argv[2])).startswith("-"):
			runALN = 1
			for s in range(1, 8):
				step_dict[s] = 1
		else:
			runALN = 2
			step_list = re.split("-", str(sys.argv[2]))
			if re.search(r'\d', step_list[1]):
				step_list[0] = int(step_list[0])
				step_list[1] = int(step_list[1]) + 1
				for s in range(step_list[0], step_list[1]):
					step_dict[s] = 1
			else:
				for s in range(int(step_list[0]), 8):
					step_dict[s] = 1
	elif sys.argv[1] == "Basic" and sys.argv[2] == "ALN":
		runALN = 2
		step_dict[3] = 1
	if runALN > 0 and 3 in step_dict:
		opts, args = getopt.gnu_getopt(sys.argv[runALN:], 'i:o:p:c:M:q:m:', ['input', 'outputdir', 'parallel', 'config', 'mergefq', 'qsub_parameter', 'minlen'])
		for o, a in opts:
			if o == '-i' or o == '--input':
				if runALN == 2 and (2 not in step_dict):
					Split_FQ_for_ALN_List = os.path.abspath(a)
			if o == '-o' or o == '--outputdir':
				OutputDir = os.path.abspath(a)
			if o == '-p' or o == '--parallel':
				ParalleleNum = a
			if o == '-c' or o == '--config':
				Basic_config = os.path.abspath(a)
			if o == '-q' or o == '--qsub_parameter':
				Qsub_parameter = str(a)

		ScriptDir = os.path.abspath(os.path.dirname(sys.argv[0]))
		ALN_script = ScriptDir + "/src/alignment.py"
		if os.path.isfile(ALN_script) == False:
			sys.stderr.write("%s does not exist, the software package might not been downloaded perfectly!" % ALN_script)
			sys.exit(-1)

		rSplit_FQ_for_ALN_List = open(Split_FQ_for_ALN_List, 'r')
		Result_List_Dir = OutputDir + "/Result_list"
		check_info(Result_List_Dir, 'dir')
		Bam_for_MARK_List = Result_List_Dir + "/Basic_MARK_input.txt"
		Basic_ALN_Result_List = Result_List_Dir + "/Basic_ALN_result.txt"
		wBam_for_MARK_List = open(Bam_for_MARK_List, 'w')
		wBasic_ALN_Result_List = open(Basic_ALN_Result_List, 'w')
		randomstring = "ALN_" + ''.join(random.sample(string.ascii_letters + string.digits, 8))
		tmpshelldir = OutputDir + "/tmp"
		check_info(tmpshelldir, "dir")
		ALNshell = tmpshelldir + "/" + randomstring + ".sh"
		wALNshell = open(ALNshell, 'w')
		for fqinfo in rSplit_FQ_for_ALN_List:
			InputFqListInfolist = re.split("\t", fqinfo.strip())
			SampleId = InputFqListInfolist[0]
			LibraryId = InputFqListInfolist[1]
			LibraryNum = InputFqListInfolist[2]
			FqPath = InputFqListInfolist[3]

#			rFqPathList = open(FqPathList, 'r')
#			for FqPath in rFqPathList:
#				FqPath = FqPath.strip()
#				FqPathBasename = (os.path.basename(FqPath)).replace(".BX.modified.fq.gz", "")
			FqPathBasename = os.path.basename(os.path.dirname(FqPath))
			#print(SampleId + "\t" + LibraryId + "\t" + FqPath + "\t" + FqPathBasename)
			clean_FQ_output_dir = OutputDir + "/" + SampleId + "/Alignment/" + LibraryId + "/" + FqPathBasename
			check_info(clean_FQ_output_dir, "dir")

			shelldir = clean_FQ_output_dir + "/shell"
			check_info(shelldir, "dir")
			FQshell = shelldir + "/run_ALN." + SampleId + ".sh"
			runFQshell = FQshell + "\n"
			wALNshell.write(runFQshell)
			wFQshell = open(FQshell, 'w')
			wFQshell.write("set -e\n")
			RGinfo = "'@RG\\tID:" + LibraryId + "\\tPL:illumina\\tPU:" + FqPathBasename + "\\tLB:" + LibraryId + "\\tSM:" + SampleId + "'"

			shell_line = " ".join([python3, ALN_script, "-i", FqPath, "-o", clean_FQ_output_dir, "-c", Basic_config, "-r", RGinfo, "-n", str(ParalleleNum), "\n"])
			wFQshell.write(shell_line)
			shell_line = " ".join(["echo Done! >", FQshell + ".sign"])
			wFQshell.write(shell_line)
			wFQshell.close()

			if "ALN" not in SCRIPT or SampleId not in SCRIPT["ALN"]:
				SCRIPT["ALN"][SampleId] = FQshell + ":12G"
			else:
				SCRIPT["ALN"][SampleId] = SCRIPT["ALN"][SampleId] + "\t" + FQshell + ":12G"
			bam = clean_FQ_output_dir + "/" + FqPathBasename + ".sorted.bam"
			baminfo = "\t".join([SampleId, LibraryId, bam]) + "\n"
			wBam_for_MARK_List.write(baminfo)
			
			Result_log = "\t".join([SampleId, LibraryId, "bam file:", bam + "\n"])
			wBasic_ALN_Result_List.write(Result_log)
#			rFqPathList.close()
		wBam_for_MARK_List.close()
		wALNshell.close()
		wBasic_ALN_Result_List.close()
		rSplit_FQ_for_ALN_List.close()

		#run_parallel(ALNshell, ParalleleNum)
		sys.stderr.write("[ %s ] files generated by 'Basic ALN' have been listed in \n\n\t %s \n\n" % (time.asctime(), Basic_ALN_Result_List))
		sys.stderr.write("[ %s ] new fq has been listed in \n\n\t %s \n\n\tAnd it's the input file in the 'Basic MARK' step\n\n" % (time.asctime(), Bam_for_MARK_List))
####################################################### ALN #################################################################

####################################################### MARK #################################################################
	runMAK = 0
	Bam_for_BQSR_List = None
	step_dict = dict()
	if sys.argv[1] == "Basicall":
		if (str(sys.argv[2])).startswith("-"):
			runMAK = 1
			for s in range(1, 8):
				step_dict[s] = 1
		else:
			runMAK = 2
			step_list = re.split("-", str(sys.argv[2]))
			if re.search(r'\d', step_list[1]):
				step_list[0] = int(step_list[0])
				step_list[1] = int(step_list[1]) + 1
				for s in range(step_list[0], step_list[1]):
					step_dict[s] = 1
			else:
				for s in range(int(step_list[0]), 8):
					step_dict[s] = 1
	elif sys.argv[1] == "Basic" and sys.argv[2] == "MARK":
		runMAK = 2
		step_dict[4] = 1
	if runMAK > 0 and 4 in step_dict:
		opts, args = getopt.gnu_getopt(sys.argv[runMAK:], 'i:o:p:c:M:q:m:', ['input', 'outputdir', 'parallel', 'config', 'mergefq', 'qsub_parameter', 'minlen'])
		for o, a in opts:
			if runMAK == 2:
				if (o == '-i' or o == '--input') and (3 not in step_dict):
					Bam_for_MARK_List = os.path.abspath(a)
			if o == '-o' or o == '--outputdir':
				OutputDir = os.path.abspath(a)
			if o == '-p' or o == '--parallel':
				ParalleleNum = a
			if o == '-c' or o == '--config':
				Basic_config = os.path.abspath(a)
			if o == '-q' or o == '--qsub_parameter':
				Qsub_parameter = str(a)

		check_info(Bam_for_MARK_List, "file")
		check_info(OutputDir, "dir")
		check_info(ParalleleNum, "num")

#		if int(ParalleleNum) > 1:
#			sys.stderr.write("The maximum number of %s CPUs would be invoked at the same time\n" % ParalleleNum)

		ScriptDir = os.path.dirname(os.path.abspath(sys.argv[0]))
		MAK_script = ScriptDir + "/src/mergeMark_bam.py"
		if os.path.exists(MAK_script) == False:
			sys.stderr.write("%s does not exist, the software package might not been downloaded perfectly!" % ALN_script)
			sys.exit(-1)
		if Basic_config != None and os.path.isfile(Basic_config):
			pass
		else:
			config_dir = OutputDir + "/config"
			check_info(config_dir, "dir")
			Basic_config = config_dir + "/Basic.config"
			print("[ %s ] Warnings: configuration file was not provided, and %s would be used!\n" % (time.asctime(), Basic_config))
			if os.path.exists(Basic_config) == False:
				sys.stderr.write("configuration file was not provided or does not exist, Please create it using 'python LRTK-SEQ.py Config'\n")
				sys.exit(-1)

		bamDict = defaultdict(dict)
		rBam_for_MARK_List = open(Bam_for_MARK_List, 'r')
		for eachBamfile in rBam_for_MARK_List:
			eachBamfile = eachBamfile.strip()
			(SampleId, LibraryId, bamfile) = re.split("\t", eachBamfile)
			#check_info(bamfile, 'file')
			if SampleId in bamDict:
				if LibraryId in bamDict[SampleId]:
					bamDict[SampleId][LibraryId] = bamDict[SampleId][LibraryId] + "\n" + bamfile
				else:
					bamDict[SampleId][LibraryId] = bamfile
			else:
				bamDict[SampleId][LibraryId] = bamfile
		rBam_for_MARK_List.close()

		tmpshelldir = OutputDir + "/tmp"
		check_info(tmpshelldir, "dir")
		MAKshell = tmpshelldir + "/MAK_" + ''.join(random.sample(string.ascii_letters + string.digits, 8)) + ".sh"
		wMAKshell = open(MAKshell, 'w')
		Result_List_Dir = OutputDir + "/Result_list"
		check_info(Result_List_Dir, 'dir')
		Basic_MARK_Result_List = Result_List_Dir + "/Basic_MARK_result.txt"
		Bam_for_BQSR_List = Result_List_Dir + "/Basic_BQSR_input.txt"
		wBasic_MARK_Result_List = open(Basic_MARK_Result_List, 'w')
		wBam_for_BQSR_List = open(Bam_for_BQSR_List, 'w')

		for samplekey in bamDict.keys():
			for librarykey in bamDict[samplekey].keys():
				library_bam_dir = OutputDir + "/" + samplekey + "/MarkDup/" + librarykey
				check_info(library_bam_dir, "dir")
				library_bam = OutputDir + "/" + samplekey + "/MarkDup/" + librarykey + "/bam.txt"
				wbam = open(library_bam, 'w')
				allbam = bamDict[samplekey][librarykey] + "\n"
				wbam.write(allbam)
				wbam.close()

				shelldir = OutputDir + "/" + samplekey + "/MarkDup/" + librarykey + "/shell"
				check_info(shelldir, "dir")
				BAMshell = shelldir + "/run_merge_mark_bam.sh"
				runBAMshell = BAMshell + "\n"
				wMAKshell.write(runBAMshell)

				wBAMshell = open(BAMshell, 'w')
				wBAMshell.write("set -e\n")
				markedbam = library_bam_dir + "/" + librarykey + ".sorted.merged.marked.bam"
				shell_line = " ".join([python3, MAK_script, "-i", library_bam, "-o", markedbam, "-c", Basic_config, "\n"])
				wBAMshell.write(shell_line)
				shell_line = " ".join(["echo Done! >", BAMshell + ".sign"]) + "\n"
				wBAMshell.write(shell_line)
				wBAMshell.close()
				if "MAK" not in SCRIPT or samplekey not in SCRIPT["MAK"]:
					SCRIPT["MAK"][samplekey] = BAMshell + ":6G"
				else:
					SCRIPT["MAK"][samplekey] = SCRIPT["MAK"][samplekey] + "\t" + BAMshell + ":6G"
				baminfo = "\t".join([samplekey, librarykey, markedbam]) + "\n"
				wBam_for_BQSR_List.write(baminfo)

				markedlog = library_bam_dir + "/" + librarykey + ".sorted.merged.marked.metrics.txt"
				Result_log = "\t".join([SampleId, LibraryId, "merged & marked bam:", markedbam + "\n"])
				wBasic_MARK_Result_List.write(Result_log)
				Result_log = "\t".join([SampleId, LibraryId, "duplication info:", markedlog + "\n"])
				wBasic_MARK_Result_List.write(Result_log)
		wMAKshell.close()
		wBam_for_BQSR_List.close()
		wBasic_MARK_Result_List.close()

		#run_parallel(MAKshell, ParalleleNum)
		sys.stderr.write("[ %s ] files generated by 'Basic MARK' have been listed in \n\n\t %s \n\n" % (time.asctime(), Basic_MARK_Result_List))
		sys.stderr.write("[ %s ] merged and duplication marked bam files of each library have been listed in \n\n\t %s \n\n\tAnd it's the input file in the 'Basic BQSR' step\n\n" % (time.asctime(), Bam_for_BQSR_List))
####################################################### MARK #################################################################

####################################################### BQSR ################################################################
	Bam_for_STAT_List = None
	runBQSR = 0
	step_dict = dict()
	if sys.argv[1] == "Basicall":
		if (str(sys.argv[2])).startswith("-"):
			runBQSR = 1
			for s in range(1, 8):
				step_dict[s] = 1
		else:
			runBQSR = 2
			step_list = re.split("-", str(sys.argv[2]))
			if re.search(r'\d', step_list[1]):
				step_list[0] = int(step_list[0])
				step_list[1] = int(step_list[1]) + 1
				for s in range(step_list[0], step_list[1]):
					step_dict[s] = 1
			else:
				for s in range(int(step_list[0]), 8):
					step_dict[s] = 1
	elif sys.argv[1] == "Basic" and sys.argv[2] == "BQSR":
		runBQSR = 2
		step_dict[5] = 1
	if runBQSR > 0 and 5 in step_dict:
		opts, args = getopt.gnu_getopt(sys.argv[runBQSR:], 'i:o:p:c:M:q:m:', ['input', 'outputdir', 'parallel', 'config', 'mergefq', 'qsub_parameter', 'minlen'])
		for o, a in opts:
			if runBQSR == 2 and (4 not in step_dict):
				if o == '-i' or o == '--input':
					Bam_for_BQSR_List = os.path.abspath(a)
			if o == '-o' or o == '--outputdir':
				OutputDir = os.path.abspath(a)
			if o == '-p' or o == '--parallel':
				ParalleleNum = a
			if o == '-c' or o == '--config':
				Basic_config = os.path.abspath(a)
			if o == '-q' or o == '--qsub_parameter':
				Qsub_parameter = str(a)

		check_info(Bam_for_BQSR_List, "file")
		check_info(OutputDir, "dir")
		check_info(ParalleleNum, "num")
#		if int(ParalleleNum) > 1:
#			sys.stderr.write("The maximum number of %s CPUs would be invoked at the same time\n" % ParalleleNum)

		ScriptDir = os.path.dirname(os.path.abspath(sys.argv[0]))
		BQSR_script = ScriptDir + "/src/bqsr.py"
		if os.path.exists(BQSR_script) == False:
			sys.stderr.write("%s does not exist, the software package might not been downloaded perfectly!" % BQSR_script)
			sys.exit(-1)
		if Basic_config == None or os.path.exists(Basic_config) == False:
			config_dir = OutputDir + "/config"
			check_info(config_dir, "dir")
			Basic_config = config_dir + "/Basic.config"
			print("[ %s ] Warnings: configuration file was not provided, and %s would be used!\n" % (time.asctime(), Basic_config))
			if os.path.exists(Basic_config) == False:
				sys.stderr.write("configuration file was not provided or does not exist, Please create it using 'python LRTK-SEQ.py Config'\n")
				sys.exit(-1)

		rBam_for_BQSR_List = open(Bam_for_BQSR_List, 'r')
		tmpshelldir = OutputDir + "/tmp"
		check_info(tmpshelldir, "dir")
		BQSRshell = tmpshelldir + "/BQSR_" + ''.join(random.sample(string.ascii_letters + string.digits, 8)) + ".sh"
		wBQSRshell = open(BQSRshell, 'w')
		Result_List_Dir = OutputDir + "/Result_list"
		check_info(Result_List_Dir, 'dir')
		BqsrResultList = Result_List_Dir + "/Basic_BQSR_result.txt"
		Bam_for_STAT_List = Result_List_Dir + "/Basic_STAT_input.txt"
		wBqsrResultList = open(BqsrResultList, 'w')
		wBam_for_STAT_List = open(Bam_for_STAT_List, 'w')
		for eachBamfile in rBam_for_BQSR_List:
			eachBamfile = eachBamfile.strip()
			(SampleId, LibraryId, bamfile) = re.split("\t", eachBamfile)
			#check_info(bamfile, 'file')

			bamdir = OutputDir + "/" + SampleId + "/BQSR/" + LibraryId
			check_info(bamdir, "dir")
			shelldir = bamdir + "/shell"
			check_info(shelldir, "dir")
			Eachshell = shelldir + "/run_bqsr.sh"
			runEachshell = Eachshell + "\n"
			wBQSRshell.write(runEachshell)

			BqsrBamFile = bamdir + "/" + os.path.basename(bamfile).replace("merged.marked.bam", "merged.marked.BQSR.bam")
			wEachshell = open(Eachshell, 'w')
			wEachshell.write("set -e\n")
			shell_line = " ".join([python3, BQSR_script, "-i", bamfile, "-o", BqsrBamFile, "-c", Basic_config, "\n"])
			wEachshell.write(shell_line)
			shell_line = " ".join(["echo Done! >", Eachshell + ".sign"]) + "\n"
			wEachshell.write(shell_line)
			wEachshell.close()

			if "BQSR" not in SCRIPT or SampleId not in SCRIPT["BQSR"]:
				SCRIPT["BQSR"][SampleId] = Eachshell + ":6G"
			else:
				SCRIPT["BQSR"][SampleId] = SCRIPT["BQSR"][SampleId] + "\t" + Eachshell + ":6G"

			BqsrBamFile_info = "\t".join([SampleId, LibraryId, BqsrBamFile]) + "\n"
			wBam_for_STAT_List.write(BqsrBamFile_info)
			Result_log = "\t".join([SampleId, LibraryId, "BQSR:", BqsrBamFile + "\n"])
			wBqsrResultList.write(Result_log)
		rBam_for_BQSR_List.close()
		wBQSRshell.close()
		wBqsrResultList.close()
		wBam_for_STAT_List.close()

		#run_parallel(BQSRshell, ParalleleNum)
		sys.stderr.write("[ %s ] files generated by 'Basic BQSR' have been listed in \n\n\t %s \n\n" % (time.asctime(), BqsrResultList))
		sys.stderr.write("[ %s ] bam files filtered/unfiltered based on molecule have been listed in \n\n\t %s \n\n\tAnd it's the input file in the 'Basic STAT' step\n\n" % (time.asctime(), Bam_for_STAT_List))
####################################################### BQSR ################################################################

####################################################### STAT ################################################################
	runSTAT = 0
	Molecule_length = str(500)
	Bam_for_MERGE_List = None
	step_dict = dict()
	if sys.argv[1] == "Basicall":
		if (str(sys.argv[2])).startswith("-"):
			runSTAT = 1
			for s in range(1, 8):
				step_dict[s] = 1
		else:
			runSTAT = 2
			step_list = re.split("-", str(sys.argv[2]))
			if re.search(r'\d', step_list[1]):
				step_list[0] = int(step_list[0])
				step_list[1] = int(step_list[1]) + 1
				for s in range(step_list[0], step_list[1]):
					step_dict[s] = 1
			else:
				for s in range(int(step_list[0]), 8):
					step_dict[s] = 1
	elif sys.argv[1] == "Basic" and sys.argv[2] == "STAT":
		runSTAT = 2
		step_dict[6] = 1
	if runSTAT > 0 and 6 in step_dict:
		opts, args = getopt.gnu_getopt(sys.argv[runSTAT:], 'i:o:m:p:c:M:q:', ['input', 'outputdir', 'minlen', 'parallel', 'config', 'mergefq', 'qsub_parameter'])
		for o, a in opts:
			if runSTAT == 2 and (5 not in step_dict):
				if o == '-i' or o == '--input':
					Bam_for_STAT_List = os.path.abspath(a)
			if o == '-o' or o == '--outputdir':
				OutputDir = os.path.abspath(a)
			if o == '-m' or o == '--minlen':
				Molecule_length = a
			if o == '-p' or o == '--parallel':
				ParalleleNum = a
			if o == '-c' or o == '--config':
				Basic_config = os.path.abspath(a)
			if o == '-q' or o == '--qsub_parameter':
				Qsub_parameter = str(a)

		check_info(Bam_for_STAT_List, "file")
		check_info(OutputDir, "dir")
		check_info(ParalleleNum, "num")
#		if int(ParalleleNum) > 1:
#			sys.stderr.write("The maximum number of %s CPUs would be invoked at the same time\n" % ParalleleNum)

		ScriptDir = os.path.dirname(os.path.abspath(sys.argv[0]))
		STAT_script = ScriptDir + "/src/calculate.py"
		if os.path.isfile(STAT_script):
			pass
		else:
			sys.stderr.write("%s does not exist, the software package might not been downloaded perfectly!" % STAT_script)
			sys.exit(-1)
		if Basic_config == None or os.path.exists(Basic_config) == False:
			config_dir = OutputDir + "/config"
			check_info(config_dir, "dir")
			Basic_config = config_dir + "/Basic.config"
			print("[ %s ] Warnings: configuration file was not provided, and %s would be used!\n" % (time.asctime(), Basic_config))
			if os.path.exists(Basic_config) == False:
				sys.stderr.write("configuration file was not provided or does not exist, Please create it using 'python LRTK-SEQ.py Config'\n")
				sys.exit(-1)

		rBam_for_STAT_List = open(Bam_for_STAT_List, 'r')
		tmpshelldir = OutputDir + "/tmp"
		check_info(tmpshelldir, "dir")
		STATshell = tmpshelldir + "/STAT_" + ''.join(random.sample(string.ascii_letters + string.digits, 8)) + ".sh"
		wSTATshell = open(STATshell, 'w')
		Result_List_Dir = OutputDir + "/Result_list"
		check_info(Result_List_Dir, 'dir')
		CFCRFileList = Result_List_Dir + "/Basic_STAT_result.txt"
		Bam_for_MERGE_List = Result_List_Dir + "/Basic_MERGE_input.txt"
		wBam_for_MERGE_List = open(Bam_for_MERGE_List, 'w')
		wCFCRFileList = open(CFCRFileList, 'w')
		for eachBamfile in rBam_for_STAT_List:
			eachBamfile = eachBamfile.strip()
			(SampleId, LibraryId, bamfile) = re.split("\t", eachBamfile)
			#check_info(bamfile, 'file')

			bamdir = OutputDir + "/" + SampleId + "/Statistics/" + LibraryId
			check_info(bamdir, "dir")
#			Result_log = "\t".join([SampleId, LibraryId, "CF&CR info:", bamdir + "/statistics/CFCR.stat\n"])
#			wCFCRFileList.write(Result_log)
#			Result_log = "\t".join([SampleId, LibraryId, "molecule details:", bamdir + "/statistics/molecule.full.gz\n"])
#			wCFCRFileList.write(Result_log)
#			Result_log = "\t".join([SampleId, LibraryId, "molecule length distribution:", bamdir + "/statistics/molecule_length_distribution.txt\n"])
#			wCFCRFileList.write(Result_log)
#			Result_log = "\t".join([SampleId, LibraryId, "coverage of each molecule:", bamdir + "/statistics/each_molecule_coverage.txt\n"])
#			wCFCRFileList.write(Result_log)
#			Result_log = "\t".join([SampleId, LibraryId, "distribution of molecule coverage:", bamdir + "/statistics/molecule_coverage_distribution.txt\n"])
#			wCFCRFileList.write(Result_log)
#			Result_log = "\t".join([SampleId, LibraryId, "molecule amounts of each barcode:", bamdir + "/statistics/barcode_molecule_amount.txt\n"])
#			wCFCRFileList.write(Result_log)

			Result_log = "\t".join([SampleId, LibraryId, "molecule info:", bamdir + "/fragment_info.csv"])
			wCFCRFileList.write(Result_log)
			Result_log = "\t".join([SampleId, LibraryId, "data summary:", bamdir + "/summary.xls"])
			wCFCRFileList.write(Result_log)

			shelldir = bamdir + "/shell"
			check_info(shelldir, "dir")
			Eachshell = shelldir + "/run_stat.sh"
			runEachshell = Eachshell + "\n"
			wSTATshell.write(runEachshell)

			wEachshell = open(Eachshell, 'w')
			wEachshell.write("set -e\n")
			shell_line = " ".join([python3, STAT_script, "-i", bamfile, "-o", bamdir, "-c", Basic_config, "-m", Molecule_length, "\n"])
			wEachshell.write(shell_line)
			shell_line = " ".join(["echo Done! >", Eachshell + ".sign"]) + "\n"
			wEachshell.write(shell_line)
			wEachshell.close()

			if "STAT" not in SCRIPT or SampleId not in SCRIPT["STAT"]:
				SCRIPT["STAT"][SampleId] = Eachshell + ":6G"
			else:
				SCRIPT["STAT"][SampleId] = SCRIPT["STAT"][SampleId] + "\t" + Eachshell + ":6G"

			FilteredBamFile = bamdir + "/molecule_filtered.sorted.bam"
			FilteredBamFile_info = "\t".join([SampleId, LibraryId, FilteredBamFile]) + "\n"
			wBam_for_MERGE_List.write(FilteredBamFile_info)
		rBam_for_STAT_List.close()
		wSTATshell.close()
		wCFCRFileList.close()
		wBam_for_MERGE_List.close()

		#run_parallel(STATshell, ParalleleNum)
		sys.stderr.write("[ %s ] files generated by 'Basic STAT' have been listed in \n\n\t %s \n\n" % (time.asctime(), CFCRFileList))
		sys.stderr.write("[ %s ] bam files filtered/unfiltered based on molecule have been listed in \n\n\t %s \n\n\tAnd it's the input file in the 'Basic MERGE' step\n\n" % (time.asctime(), Bam_for_MERGE_List))
####################################################### STAT ################################################################

####################################################### MERGE ###############################################################
	Bam_for_Varcall_List = None
	runMERGE = 0
	step_dict = dict()
	if sys.argv[1] == "Basicall":
		if (str(sys.argv[2])).startswith("-"):
			runMERGE = 1
			for s in range(1, 8):
				step_dict[s] = 1
		else:
			runMERGE = 2
			step_list = re.split("-", str(sys.argv[2]))
			if re.search(r'\d', step_list[1]):
				step_list[0] = int(step_list[0])
				step_list[1] = int(step_list[1]) + 1
				for s in range(step_list[0], step_list[1]):
					step_dict[s] = 1
			else:
				for s in range(int(step_list[0]), 8):
					step_dict[s] = 1
	elif sys.argv[1] == "Basic" and sys.argv[2] == "MERGE":
		runMERGE = 2
		step_dict[7] = 1
	if runMERGE > 0 and 7 in step_dict:
		opts, args = getopt.gnu_getopt(sys.argv[runMERGE:], 'i:o:p:c:M:q:m:', ['input', 'outputdir', 'parallel', 'config', 'mergefq', 'qsub_parameter', 'minlen'])
		for o, a in opts:
			if runMERGE == 2 and (6 not in step_dict):
				if o == '-i' or o == '--input':
					Bam_for_MERGE_List = os.path.abspath(a)
			if o == '-o' or o == '--outputdir':
				OutputDir = os.path.abspath(a)
			if o == '-p' or o == '--parallel':
				ParalleleNum = a
			if o == '-c' or o == '--config':
				Basic_config = os.path.abspath(a)
			if o == '-q' or o == '--qsub_parameter':
				Qsub_parameter = str(a)

		check_info(Bam_for_MERGE_List, "file")
		check_info(OutputDir, "dir")
		check_info(ParalleleNum, "num")
#		if int(ParalleleNum) > 1:
#			sys.stderr.write("The maximum number of %s CPUs would be invoked at the same time\n" % ParalleleNum)

		ScriptDir = os.path.dirname(os.path.abspath(sys.argv[0]))
		MERGE_script = ScriptDir + "/src/mergebam.py"
		if os.path.isfile(MERGE_script):
			pass
		else:
			sys.stderr.write("%s does not exist, the software package might not been downloaded perfectly!" % MERGE_script)
			sys.exit(-1)
		if Basic_config == None or os.path.exists(Basic_config) == False:
			config_dir = OutputDir + "/config"
			check_info(config_dir, "dir")
			Basic_config = config_dir + "/Basic.config"
			print("[ %s ] Warnings: configuration file was not provided, and %s would be used!\n" % (time.asctime(), Basic_config))
			if os.path.exists(Basic_config) == False:
				sys.stderr.write("configuration file was not provided or does not exist, Please create it using 'python LRTK-SEQ.py Config'\n")
				sys.exit(-1)

		bamDict = dict()
		rBam_for_MERGE_List = open(Bam_for_MERGE_List, 'r')
		for eachBamfile in rBam_for_MERGE_List:
			eachBamfile = eachBamfile.strip()
			(SampleId, LibraryId, bamfile) = re.split("\t", eachBamfile)
			#check_info(bamfile, 'file')
			if SampleId in bamDict:
				bamDict[SampleId] = bamDict[SampleId] + "\n" + bamfile
			else:
				bamDict[SampleId] = bamfile
		rBam_for_MERGE_List.close()

		tmpshelldir = OutputDir + "/tmp"
		check_info(tmpshelldir, "dir")
		MERGEshell = tmpshelldir + "/MERGE_" + ''.join(random.sample(string.ascii_letters + string.digits, 8)) + ".sh"
		wMERGEshell = open(MERGEshell, 'w')
		Result_List_Dir = OutputDir + "/Result_list"
		check_info(Result_List_Dir, 'dir')
		Bam_for_Varcall_List = Result_List_Dir + "/Reseq_Varcall_input.txt"
		Basic_MERGE_Result_List = Result_List_Dir + "/Basic_MERGE_result.txt"
		wBam_for_Varcall_List = open(Bam_for_Varcall_List, 'w')
		wBasic_MERGE_Result_List = open(Basic_MERGE_Result_List, 'w')

		for samplekey in bamDict.keys():
			sample_bam_dir = OutputDir + "/" + samplekey + "/Merge"
			check_info(sample_bam_dir, "dir")
			sample_bam = sample_bam_dir + "/bam.txt"
			wbam = open(sample_bam, 'w')
			allbam = bamDict[samplekey] + "\n"
			wbam.write(allbam)
			wbam.close()

			shelldir = sample_bam_dir + "/shell"
			check_info(shelldir, "dir")
			BAMshell = shelldir + "/run_merge_bam.sh"
			runBAMshell = BAMshell + "\n"
			wMERGEshell.write(runBAMshell)
			wBAMshell = open(BAMshell, 'w')
			wBAMshell.write("set -e\n")
			mergedbam = sample_bam_dir + "/" + samplekey + ".sorted.merged.marked.BQSR.filtered.bam"
			shell_line = " ".join([python3, MERGE_script, "-i", sample_bam, "-o", mergedbam, "-c", Basic_config, "\n"])
			wBAMshell.write(shell_line)
			shell_line = " ".join(["echo Done! >", BAMshell + ".sign"]) + "\n"
			wBAMshell.write(shell_line)
			wBAMshell.close()

			if "MERGE" not in SCRIPT or samplekey not in SCRIPT["MERGE"]:
				SCRIPT["MERGE"][samplekey] = BAMshell + ":3G"
			else:
				SCRIPT["MERGE"][samplekey] = SCRIPT["MERGE"][samplekey] + "\t" + BAMshell + ":3G"

			baminfo = "\t".join([samplekey, mergedbam]) + "\n"
			wBam_for_Varcall_List.write(baminfo)
			Result_log = "\t".join([SampleId, "final alignment result:", mergedbam + "\n"])
			wBasic_MERGE_Result_List.write(Result_log)
		wMERGEshell.close()
		wBam_for_Varcall_List.close()
		wBasic_MERGE_Result_List.close()

		#run_parallel(MERGEshell, ParalleleNum)
		sys.stderr.write("[ %s ] files generated by 'Basic MERGE' have been listed in \n\n\t %s \n\n" % (time.asctime(), Basic_MERGE_Result_List))
		sys.stderr.write("[ %s ] finally merged bam files have been listed in \n\t %s \n\n\tAnd it's the input file in the 'Reseq Varcall' and 'Reseq Phasing' step\n\n" % (time.asctime(), Bam_for_Varcall_List))
####################################################### MERGE ###############################################################

####################################################### Varcall #############################################################
	ScriptDir = os.path.dirname(os.path.abspath(sys.argv[0]))
	Bam_VCF_for_Phasing_List = None
	GVCF_for_Familycall_List = None

	runVAR = 0
	Chrlist = None
	if sys.argv[1] == "Reseqall":
		runVAR = 1
	elif sys.argv[1] == "Reseq" and sys.argv[2] == "Varcall":
		runVAR = 2
	if runVAR > 0:

		all_command = " ".join(sys.argv)
		sys.stderr.write("Command:\t%s\n\n" % all_command)

		opts, args = getopt.gnu_getopt(sys.argv[runVAR:], 'i:o:L:p:c:M:q:f:', ['input', 'outputdir', 'chrlist', 'parallel', 'config', 'qsub_parameter', 'familyinfo'])
		for o, a in opts:
			if o == '-i' or o == '--input':
				Bam_for_Varcall_List = os.path.abspath(a)
			if o == '-o' or o == '--outputdir':
				OutputDir = os.path.abspath(a)
			if o == '-L' or o == '--chrlist':
				Chrlist = os.path.abspath(a)
			if o == '-p' or o == '--parallel':
				ParalleleNum = a
			if o == '-c' or o == '--config':
				Reseq_config = os.path.abspath(a)
			if o == '-q' or o == '--qsub_parameter':
				Qsub_parameter = str(a)

		check_info(Bam_for_Varcall_List, "file")
		check_info(OutputDir, "dir")
		check_info(ParalleleNum, "num")
		config_dir = OutputDir + "/config"
		check_info(config_dir, "dir")

#		if int(ParalleleNum) > 1:
#			sys.stderr.write("The maximum number of %s CPUs would be invoked at the same time\n" % ParalleleNum)

		ScriptDir = os.path.dirname(os.path.abspath(sys.argv[0]))
		VAR_script = ScriptDir + "/src/SnpInDel_call.py"
		if os.path.isfile(VAR_script):
			pass
		else:
			sys.stderr.write("%s does not exist, the software package might not been downloaded perfectly!" % VAR_script)
			sys.exit(-1)
		if Reseq_config == None or os.path.exists(Reseq_config) == False:
			Reseq_config = config_dir + "/Reseq.config"
			print("[ %s ] Warnings: configuration file was not provided, and %s would be used!\n" % (time.asctime(), Reseq_config))
			if os.path.exists(Reseq_config) == False:
				sys.stderr.write("configuration file was not provided or does not exist, Please create it using 'python LRTK-SEQ.py Config'\n")
				sys.exit(-1)

		if Chrlist == None or os.path.exists(Chrlist) == False:
			Chrlist = OutputDir + "/config/chrlist.txt"

		Result_List_Dir = OutputDir + "/Result_list"
		check_info(Result_List_Dir, 'dir')
		Reseq_Varcall_Result_List = Result_List_Dir + "/Reseq_Varcall_result.txt"
		Bam_VCF_for_Phasing_List = Result_List_Dir + "/Reseq_Phasing_input.txt"
		GVCF_for_Familycall_List = Result_List_Dir + "/Reseq_Familycall_input.txt"
		wReseq_Varcall_Result_List = open(Reseq_Varcall_Result_List, 'w')
		wBam_VCF_for_Phasing_List = open(Bam_VCF_for_Phasing_List, 'w')
		wGVCF_for_Familycall_List = open(GVCF_for_Familycall_List, 'w')
		rBam_for_Varcall_List = open(Bam_for_Varcall_List, 'r')
		tmpshelldir = OutputDir + "/tmp"
		check_info(tmpshelldir, "dir")
		VARshell = tmpshelldir + "/VAR_" + ''.join(random.sample(string.ascii_letters + string.digits, 8)) + ".sh"
		wVARshell = open(VARshell, 'w')
		for eachBamfile in rBam_for_Varcall_List:
			eachBamfile = eachBamfile.strip()
			(SampleId, bamfile) = re.split("\t", eachBamfile)
			#check_info(bamfile, 'file')

			SampleDir = OutputDir + "/" + SampleId + "/Variation_Call"
			check_info(SampleDir, "dir")
			shelldir = SampleDir + "/shell"
			check_info(shelldir, "dir")

			Eachshell = shelldir + "/run_variant_call.sh"
			runEachshell = Eachshell + "\n"
			wVARshell.write(runEachshell)

			wEachshell = open(Eachshell, 'w')
			wEachshell.write("set -e\n")
			shell_line = " ".join([python3, VAR_script, "-i", bamfile, "-o", SampleDir, "-c", Reseq_config, "-L", Chrlist, "\n"])
			wEachshell.write(shell_line)
			wEachshell.close()

			subprocess.call(["sh", Eachshell])

			SCRIPT["Varcall"][SampleId] = SampleDir + "/qsub.shell.txt"

			unphased_vcf = SampleId + '\t' + bamfile + "\t" + SampleDir + "/all.SNP.INDEL.VQSR.vcf.gz\n"
			wBam_VCF_for_Phasing_List.write(unphased_vcf)

			gvcf_path = SampleId + '\t' + SampleDir + "/gvcf.list\n"
			wGVCF_for_Familycall_List.write(gvcf_path)

			Result_log = "\t".join([SampleId, "unphased vcf file:", unphased_vcf])
			wReseq_Varcall_Result_List.write(Result_log)
		rBam_for_Varcall_List.close()
		wVARshell.close()
		wBam_VCF_for_Phasing_List.close()
		wReseq_Varcall_Result_List.close()
		wGVCF_for_Familycall_List.close()

		#run_parallel(VARshell, ParalleleNum)
		sys.stderr.write("[ %s ] files generated by 'Reseq Varcall' have been listed in \n\n\t %s \n\n" % (time.asctime(), Reseq_Varcall_Result_List))
		sys.stderr.write("[ %s ] Variation call results have been listed in \n\n\t %s \n\n\tAnd it's the input file in the 'Reseq Phasing' step\n\n" % (time.asctime(), Bam_VCF_for_Phasing_List))
####################################################### Varcall #############################################################

####################################################### family based Varcall ################################################
	runFBV = 0
	Chrlist = None
	Full_Family_Info_List = None
	if sys.argv[1] == "Reseqall":
		runFBV = 1
	elif sys.argv[1] == "Reseq" and sys.argv[2] == "Famcall":
		runFBV = 2
	if runFBV > 0:
		all_command = " ".join(sys.argv)
		sys.stderr.write("Command:\t%s\n\n" % all_command)

		opts, args = getopt.gnu_getopt(sys.argv[runFBV:], 'i:o:L:p:c:f:q:', ['input', 'outputdir', 'chrlist', 'parallel', 'config', 'familyinfo', 'qsub_parameter'])
		for o, a in opts:
			if o == '-i' or o == '--input':
				if runFBV == 2:
					GVCF_for_Familycall_List = os.path.abspath(a)
			if o == '-o' or o == '--outputdir':
				OutputDir = os.path.abspath(a)
			if o == '-L' or o == '--chrlist':
				Chrlist = os.path.abspath(a)
			if o == '-p' or o == '--parallel':
				ParalleleNum = a
			if o == '-c' or o == '--config':
				Reseq_config = os.path.abspath(a)
			if o == '-f' or o == '--familyinfo':
				Full_Family_Info_List = os.path.abspath(a)
			if o == '-q' or o == '--qsub_parameter':
				Qsub_parameter = str(a)

		check_info(OutputDir, "dir")
		check_info(ParalleleNum, "num")
		check_info(Full_Family_Info_List, "file")
		check_info(GVCF_for_Familycall_List, "file")
		config_dir = OutputDir + "/config"
		check_info(config_dir, "dir")

#		if int(ParalleleNum) > 1:
#			sys.stderr.write("The maximum number of %s CPUs would be invoked at the same time\n" % ParalleleNum)

		ScriptDir = os.path.dirname(os.path.abspath(sys.argv[0]))
		FBV_script = ScriptDir + "/src/Family_SnpInDel_call.py"
		if os.path.isfile(FBV_script):
			pass
		else:
			sys.stderr.write("%s does not exist, the software package might not been downloaded perfectly!" % FBV_script)
			sys.exit(-1)
		if Reseq_config == None or os.path.exists(Reseq_config) == False:
			Reseq_config = config_dir + "/Reseq.config"
			print("[ %s ] Warnings: configuration file was not provided, and %s would be used!\n" % (time.asctime(), Reseq_config))
			if os.path.exists(Reseq_config) == False:
				sys.stderr.write("configuration file was not provided or does not exist, Please create it using 'python LRTK-SEQ.py Config'\n")
				sys.exit(-1)

		if Chrlist == None or os.path.exists(Chrlist) == False:
			Chrlist = OutputDir + "/config/chrlist.txt"

		Result_List_Dir = OutputDir + "/Result_list"
		check_info(Result_List_Dir, 'dir')
		Reseq_Famcall_Result_List = Result_List_Dir + "/Reseq_Familycall_result.txt"
		wReseq_Famcall_Result_List = open(Reseq_Famcall_Result_List, 'w')
		tmpshelldir = OutputDir + "/tmp"
		check_info(tmpshelldir, "dir")
		FBVshell = tmpshelldir + "/FBV_" + ''.join(random.sample(string.ascii_letters + string.digits, 8)) + ".sh"
		wFBVshell = open(FBVshell, 'w')
		FamilyDataDir = OutputDir + "/Family"
		check_info(FamilyDataDir, "dir")

		gvcf_dict = dict()
		rGVCF_for_Familycall_List = open(GVCF_for_Familycall_List, 'r')
		for GVCF_for_Familycall_List_info in rGVCF_for_Familycall_List:
			(SampleId, gvcf_path) = re.split("\t", GVCF_for_Familycall_List_info.strip())
			gvcf_dict[SampleId] = gvcf_path
		rGVCF_for_Familycall_List.close()

		fam_dict = defaultdict(list)
		sample_uniqe_dict = dict()
		rFull_Family_Info_List = open(Full_Family_Info_List, 'r')
		for Full_Family_Info_List_info in rFull_Family_Info_List:
			Full_Family_Info_List_info_list = re.split("\t", Full_Family_Info_List_info.strip())
			FamilyId = Full_Family_Info_List_info_list[0]
			SampleId = Full_Family_Info_List_info_list[1]
			if SampleId not in sample_uniqe_dict:
				sample_uniqe_dict[SampleId] = 1
				fam_dict[FamilyId].append(SampleId)
		rFull_Family_Info_List.close()

		for FamilyId in fam_dict.keys():
			FamilyDir = FamilyDataDir + "/" + FamilyId + "/Variation_Call"
			check_info(FamilyDir, "dir")
			shelldir = FamilyDir + "/shell"
			check_info(shelldir, "dir")

			famgvcflist = FamilyDir + "/" + FamilyId + ".gvcf.list"
			wfamgvcflist = open(famgvcflist, 'w')
			for SampleId in fam_dict[FamilyId]:
				wfamgvcflist.write(gvcf_dict[SampleId] + "\n")

				SCRIPT["FBV"][SampleId] = FamilyDir + "/qsub.shell.txt"

			wfamgvcflist.close()

			Eachshell = shelldir + "/run_variant_call.sh"
			runEachshell = Eachshell + "\n"
			wFBVshell.write(runEachshell)

			wEachshell = open(Eachshell, 'w')
			wEachshell.write("set -e\n")
			shell_line = " ".join([python3, FBV_script, "-i", famgvcflist, "-o", FamilyDir, "-c", Reseq_config, "-L", Chrlist, "\n"])
			wEachshell.write(shell_line)
			wEachshell.close()
			family_vcf = FamilyId + "\t" + FamilyDir + "/all.SNP.INDEL.VQSR.vcf.gz" + "\n"
			Result_log = "\t".join([FamilyId, "family-based vcf file:", family_vcf]) + "\n"
			wReseq_Famcall_Result_List.write(family_vcf)

			subprocess.call(["sh", Eachshell])
		wReseq_Famcall_Result_List.close()
		wFBVshell.close()

		#run_parallel(FBVshell, ParalleleNum)
		sys.stderr.write("[ %s ] files generated by 'Reseq Famcall' have been listed in \n\n\t %s \n\n" % (time.asctime(), Reseq_Varcall_Result_List))

####################################################### family based Varcall ################################################

####################################################### Phasing #############################################################
	runPHASE = 0
	if sys.argv[1] == "Reseqall":
		runPHASE = 1
	elif sys.argv[1] == "Reseq" and sys.argv[2] == "Phasing":
		runPHASE = 2
	if runPHASE > 0:
		opts, args = getopt.gnu_getopt(sys.argv[runPHASE:], 'i:o:p:v:c:M:q:f:', ['input', 'outputdir', 'vcf', 'parallel', 'config', 'mergefq', 'qsub_parameter', 'familyinfo'])
		for o, a in opts:
			if o == '-i' or o == '--input':
				if runPHASE == 2:
					Bam_VCF_for_Phasing_List = os.path.abspath(a)
			if o == '-o' or o == '--outputdir':
				OutputDir = os.path.abspath(a)
			if o == '-p' or o == '--parallel':
				ParalleleNum = a
			if o == '-c' or o == '--config':
				Reseq_config = os.path.abspath(a)
			if o == '-q' or o == '--qsub_parameter':
				Qsub_parameter = str(a)

		check_info(Bam_VCF_for_Phasing_List, "file")
		check_info(OutputDir, "dir")
		check_info(ParalleleNum, "num")

#		if int(ParalleleNum) > 1:
#			sys.stderr.write("The maximum number of %s CPUs would be invoked at the same time\n" % ParalleleNum)

		ScriptDir = os.path.dirname(os.path.abspath(sys.argv[0]))
		PHA_script = ScriptDir + "/src/phasing.py"
		if os.path.isfile(PHA_script):
			pass
		else:
			sys.stderr.write("%s does not exist, the software package might not been downloaded perfectly!" % PHA_script)
			sys.exit(-1)
		if Reseq_config == None or os.path.exists(Reseq_config) == False:
			config_dir = OutputDir + "/config"
			check_info(config_dir, "dir")
			Reseq_config = config_dir + "/Reseq.config"
			if os.path.exists(Reseq_config) == False:
				sys.stderr.write("configuration file was not provided or does not exist, Please create it using 'python LRTK-SEQ.py Config'\n")
				sys.exit(-1)

		BamDict = dict()
		rBam_VCF_for_Phasing_List = open(Bam_VCF_for_Phasing_List, 'r')
		VcfDict = dict()
		for eachVcffile in rBam_VCF_for_Phasing_List:
			eachVcffile = eachVcffile.strip()
			tlist = re.split("\t", eachVcffile)
			if len(tlist) != 3:
				sys.stderr.write("[ %s ] ERROR: vcf file has not enough columns, it should be '1. sample id; 2. indexed bam file; 3. indexed vcf file'\n" % time.asctime())
				sys.exit(-1)
			(SampleId, bamfile, vcffile) = re.split("\t", eachVcffile)
			check_info(bamfile, 'file')
			BamDict[SampleId] = bamfile
			#check_info(vcffile, 'file')
			VcfDict[SampleId] = vcffile
		rBam_VCF_for_Phasing_List.close()

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
				sampledir = OutputDir + "/" + SampleId + "/" + "Phasing"
				check_info(sampledir, "dir")
				shelldir = sampledir + "/shell"
				check_info(shelldir, "dir")
				Eachshell = shelldir + "/run_phase.sh"
				runEachshell = Eachshell + "\n"
				wPHASEshell.write(runEachshell)

				wEachshell = open(Eachshell, 'w')
				wEachshell.write("set -e\n")
				shell_line = " ".join([python3, PHA_script, "-i", bamfile, "-v", unphasevcffile, "-o", sampledir, "-c", Reseq_config, "\n"])
				wEachshell.write(shell_line)
				shell_line = " ".join(["echo Done! >", Eachshell + ".sign\n"])
				wEachshell.write(shell_line)
				wEachshell.close()

				SCRIPT["Phasing"][SampleId] = Eachshell + ":6G"

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

		#run_parallel(PHASEshell, ParalleleNum)
		sys.stderr.write("[ %s ] files generated by 'Reseq Phasing' have been listed in \n\n\t %s \n\n" % (time.asctime(), phaseVcfList))
####################################################### Phasing #############################################################

####################################################### qsub jobs ###########################################################
	split_shell_for_sample = defaultdict(list)
	if len(SCRIPT) > 0:
		basic_qsub = OutputDir + "/qsub.txt"
		wbasic_qsub = open(basic_qsub, 'w')
		if "PRE" in SCRIPT:
			for SampleId in SCRIPT["PRE"].keys():
				sub_scr_list = re.split("\t", SCRIPT["PRE"][SampleId])
				for sub_scr in sub_scr_list:
					wbasic_qsub.write(sub_scr + "\n")

				shell_run =  ";".join(sub_scr_list)
				split_shell_for_sample[SampleId].append(shell_run)
		if "MFQ" in SCRIPT:
			for SampleId in SCRIPT["MFQ"].keys():
				sub_scr_list = re.split("\t", SCRIPT["MFQ"][SampleId])
				if len(sub_scr_list) == 1:
					if ("PRE" in SCRIPT) and (SampleId in SCRIPT["PRE"]):
						pre_scr_list = re.split("\t", SCRIPT["PRE"][SampleId])
						for pre_scr in pre_scr_list:
							wbasic_qsub.write(pre_scr + "\t" + sub_scr_list[0] + "\n")
					else:
						wbasic_qsub.write(sub_scr_list[0] + "\n")
					split_shell_for_sample[SampleId].append(sub_scr_list[0])
				else:
					if ("PRE" in SCRIPT) and (SampleId in SCRIPT["PRE"]):
						pre_scr_list = re.split("\t", SCRIPT["PRE"][SampleId])
						for sub_scr in sub_scr_list:
							for pre_scr in pre_scr_list:
								sub_scr_info = re.split("/", sub_scr)
								pre_scr_info = re.split("/", pre_scr)
								if pre_scr_info[-4] == sub_scr_info[-3]:
									wbasic_qsub.write(pre_scr + "\t" + sub_scr + "\n")
					else:
						for sub_scr in sub_scr_list:
							wbasic_qsub.write(sub_scr + "\n")
					shell_run =  ";".join(sub_scr_list)
					split_shell_for_sample[SampleId].append(shell_run)
		if "ALN" in SCRIPT:
			for SampleId in SCRIPT["ALN"].keys():
				sub_scr_list = re.split("\t", SCRIPT["ALN"][SampleId])
				shell_run =  ";".join(sub_scr_list)
				split_shell_for_sample[SampleId].append(shell_run)
				if ("MFQ" in SCRIPT) and (SampleId in SCRIPT["MFQ"]):
					pre_scr_list = re.split("\t", SCRIPT["MFQ"][SampleId])
					if len(pre_scr_list) == 1:
						for sub_scr in sub_scr_list:
							wbasic_qsub.write(pre_scr_list[0] + "\t" + sub_scr + "\n")
					else:
						for sub_scr in sub_scr_list:
							for pre_scr in pre_scr_list:
								sub_scr_info = re.split("/", sub_scr)
								pre_scr_info = re.split("/", pre_scr)
								if pre_scr_info[-3] == sub_scr_info[-3]:
									wbasic_qsub.write(pre_scr + "\t" + sub_scr + "\n")
				else:
					if ("PRE" in SCRIPT) and (SampleId in SCRIPT["PRE"]):
						pre_scr_list = re.split("\t", SCRIPT["PRE"][SampleId])
						for sub_scr in sub_scr_list:
							for pre_scr in pre_scr_list:
								sub_scr_info = re.split("/", sub_scr)
								pre_scr_info = re.split("/", pre_scr)
								if pre_scr_info[-3] == sub_scr_info[-3]:
									wbasic_qsub.write(pre_scr + "\t" + sub_scr + "\n")
					else:
						for sub_scr in sub_scr_list:
							wbasic_qsub.write(sub_scr + "\n")
		if "MAK" in SCRIPT:
			for SampleId in SCRIPT["MAK"].keys():
				sub_scr_list = re.split("\t", SCRIPT["MAK"][SampleId])
				shell_run =  ";".join(sub_scr_list)
				split_shell_for_sample[SampleId].append(shell_run)
				if ("ALN" in SCRIPT) and (SampleId in SCRIPT["ALN"]):
					pre_scr_list = re.split("\t", SCRIPT["ALN"][SampleId])
					for sub_scr in sub_scr_list:
						for pre_scr in pre_scr_list:
							sub_scr_info = re.split("/", sub_scr)
							pre_scr_info = re.split("/", pre_scr)
							if pre_scr_info[-4] == sub_scr_info[-3]:
								wbasic_qsub.write(pre_scr + "\t" + sub_scr + "\n")
				else:
					for sub_scr in sub_scr_list:
						wbasic_qsub.write(sub_scr + "\n")
		if "BQSR" in SCRIPT:
			for SampleId in SCRIPT["BQSR"].keys():
				sub_scr_list = re.split("\t", SCRIPT["BQSR"][SampleId])
				shell_run =  ";".join(sub_scr_list)
				split_shell_for_sample[SampleId].append(shell_run)
				if ("MAK" in SCRIPT) and (SampleId in SCRIPT["MAK"]):
					pre_scr_list = re.split("\t", SCRIPT["MAK"][SampleId])
					for sub_scr in sub_scr_list:
						for pre_scr in pre_scr_list:
							sub_scr_info = re.split("/", sub_scr)
							pre_scr_info = re.split("/", pre_scr)
							if pre_scr_info[-3] == sub_scr_info[-3]:
								wbasic_qsub.write(pre_scr + "\t" + sub_scr + "\n")
				else:
					for sub_scr in sub_scr_list:
						wbasic_qsub.write(sub_scr + "\n")
		if "STAT" in SCRIPT:
			for SampleId in SCRIPT["STAT"].keys():
				sub_scr_list = re.split("\t", SCRIPT["STAT"][SampleId])
				shell_run =  ";".join(sub_scr_list)
				split_shell_for_sample[SampleId].append(shell_run)
				if ("BQSR" in SCRIPT) and (SampleId in SCRIPT["BQSR"]):
					pre_scr_list = re.split("\t", SCRIPT["BQSR"][SampleId])
					for sub_scr in sub_scr_list:
						for pre_scr in pre_scr_list:
							sub_scr_info = re.split("/", sub_scr)
							pre_scr_info = re.split("/", pre_scr)
							if pre_scr_info[-3] == sub_scr_info[-3]:
								wbasic_qsub.write(pre_scr + "\t" + sub_scr + "\n")
				else:
					for sub_scr in sub_scr_list:
						wbasic_qsub.write(sub_scr + "\n")
		if "MERGE" in SCRIPT:
			for SampleId in SCRIPT["MERGE"].keys():
				sub_scr = SCRIPT["MERGE"][SampleId]
				split_shell_for_sample[SampleId].append(sub_scr)
				if ("STAT" in SCRIPT) and (SampleId in SCRIPT["STAT"]):
					pre_scr_list = re.split("\t", SCRIPT["STAT"][SampleId])
					for pre_scr in pre_scr_list:
						wbasic_qsub.write(pre_scr + "\t" + sub_scr + "\n")
				else:
					wbasic_qsub.write(sub_scr + "\n")
		if "Varcall" in SCRIPT:
			for SampleId in SCRIPT["Varcall"].keys():
				sub_scr_file = SCRIPT["Varcall"][SampleId]
				rsub_scr_file = open(sub_scr_file, 'r')
				shell_run = "A"
				combineshell = "B"
				for sub_scr_file_info in rsub_scr_file:
					(ch, chrshell, combineshell) = re.split("\t", sub_scr_file_info.strip())
					if "MERGE" in SCRIPT and SampleId in SCRIPT["MERGE"]:
						pre_scr = SCRIPT["MERGE"][SampleId]
						wbasic_qsub.write(pre_scr + "\t" + chrshell + "\n")
					if shell_run == "A":
						shell_run = chrshell
					else:
						shell_run = shell_run + ";" + chrshell
					wbasic_qsub.write(chrshell + "\t" + combineshell + "\n")
				split_shell_for_sample[SampleId].append(shell_run)
				split_shell_for_sample[SampleId].append(combineshell)
				rsub_scr_file.close()
		if "FBV" in SCRIPT:
			for SampleId in SCRIPT["FBV"].keys():
				pre_gvcf_shell = dict()
				shell_run = "A"
				combineshell = "B"
				if "Varcall" in SCRIPT and SampleId in SCRIPT["Varcall"]:
					rgvcfshell = open(SCRIPT["Varcall"][SampleId], 'r')
					for gvcfshellinfo in rgvcfshell:
						(ch, chrshell, combineshell) = re.split("\t", gvcfshellinfo.strip())
						pre_gvcf_shell[ch] = chrshell
					rgvcfshell.close()
				sub_scr_file = SCRIPT["FBV"][SampleId]
				rsub_scr_file = open(sub_scr_file, 'r')
				for sub_scr_file_info in rsub_scr_file:
					(ch, chrshell, combineshell) = re.split("\t", sub_scr_file_info.strip())
					if shell_run == "A":
						shell_run = chrshell
					else:
						shell_run = shell_run + ";" + chrshell
					if ch in pre_gvcf_shell:
						wbasic_qsub.write(pre_gvcf_shell[ch] + "\t" + chrshell + "\n")
					wbasic_qsub.write(chrshell + "\t" + combineshell + "\n")
				rsub_scr_file.close()
				split_shell_for_sample[SampleId].append(shell_run)
				split_shell_for_sample[SampleId].append(combineshell)
		if "Phasing" in SCRIPT:
			for SampleId in SCRIPT["Phasing"].keys():
				sub_scr = SCRIPT["Phasing"][SampleId]
				if "Varcall" in SCRIPT and SampleId in SCRIPT["Varcall"]:
					rgvcfshell = open(SCRIPT["Varcall"][SampleId], 'r')
					for gvcfshellinfo in rgvcfshell:
						(ch, chrshell, combineshell) = re.split("\t", gvcfshellinfo.strip())
						wbasic_qsub.write(combineshell + "\t" +  sub_scr + "\n")
					rgvcfshell.close()
				else:
					wbasic_qsub.write(sub_scr + "\n")
				split_shell_for_sample[SampleId].append(sub_scr)
		wbasic_qsub.close()
		
#		qsub_command = " ".join(["perl", os.path.dirname(os.path.abspath(sys.argv[0])) + "/monitor.pl", "-check_in_advance -i", basic_qsub, "-o", OutputDir, "-p", str(ParalleleNum), Qsub_parameter + "\n"])
#		print("\n\n\n %s \n\n\n" % qsub_command)
		print("[ %s ] All shells have been listed in %s \n " % (time.asctime(), basic_qsub))
		print("[ %s ] And shells listed in the second row must run after the shells in the first row!" % time.asctime())

		

####################################################### qsub jobs ###########################################################
