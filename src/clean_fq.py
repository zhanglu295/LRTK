import os, sys, gzip
import getopt
import time
import subprocess
import re
import multiprocessing
import random, string
import pysam
import glob

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
	
	def Barcode(self):
		abs_path = self.get_path('barcode_fa')
		return abs_path

	def Barcode_aln_par(self):
		abs_path = self.get_path('barcode_aln_parameter')
		return abs_path

	def Samtools(self):
		abs_path = self.get_path('samtools')
		return abs_path

	def Fastqc(self):
		abs_path = self.get_path('fastqc')
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

class sequence_quality:
	def countfq(self, quality_info, basic_quality):
		q20 = 0
		q30 = 0
		_quality_list = str(quality_info)
		for qua in _quality_list:
			ordqua = ord(qua) - int(basic_quality)
			if ordqua >= 20:
				q20 += 1
			if ordqua >= 30:
				q30 += 1
		return(q20, q30)

	def check_quality_version(self, quality_info):
		lennum = len(re.compile(r'[a-z]').findall(quality_info))
		qualen = len(quality_info)
		QUALITY = 33
		if lennum / qualen > 0.1:
			QUALITY = 64
		return QUALITY

class extract_barcode:
	def sub_reviese(self, splitfq1, splitFqFile, splitBarcodeFile, splitfq2):
		o_rawfq = open(splitfq1, 'r')
		if splitfq2 != 0:
			o_rawfq2 = open(splitfq2, 'r')
		wsplitFqFile = gzip.open(splitFqFile, 'wb')
		wsplitBarcodeFile = gzip.open(splitBarcodeFile, 'wb')

		while 1:
			rawfq_id1 = o_rawfq.readline()
			rawfq_seq1 = o_rawfq.readline()
			rawfq_str1 = o_rawfq.readline()
			rawfq_qua1 = o_rawfq.readline()

			if splitfq2 == 0:
				rawfq_id2 = o_rawfq.readline()
				rawfq_seq2 = o_rawfq.readline()
				rawfq_str2 = o_rawfq.readline()
				rawfq_qua2 = o_rawfq.readline()
			else:
				rawfq_id2 = o_rawfq2.readline()
				rawfq_seq2 = o_rawfq2.readline()
				rawfq_str2 = o_rawfq2.readline()
				rawfq_qua2 = o_rawfq2.readline()

			if len(rawfq_id1) == 0:
				break
			else:
				barcode_line = rawfq_id1 + rawfq_seq1[0:16] + '\n' + rawfq_str1 + rawfq_qua1[0:16] + "\n"
				fq_line = rawfq_id1 + rawfq_seq1 + rawfq_str1 + rawfq_qua1 + rawfq_id2 + rawfq_seq2 + rawfq_str2 + rawfq_qua2

				wsplitFqFile.write(fq_line.encode())
				wsplitBarcodeFile.write(barcode_line.encode())
		print()
		o_rawfq.close()
		if splitfq2 != 0:
			o_rawfq2.close()
		wsplitFqFile.close()
		wsplitBarcodeFile.close()

	def revise(self, rawfq, fq_prefix_name, outdir, maxsize, _split_done, rawfq2):
		(rawfq_path, rawfq_name) = os.path.split(rawfq)

		if _split_done == 1:
			OutFileList = os.path.join(outdir, "split_fq_barcode.txt")
			return(OutFileList)

		split_shell = outdir + "/split_fq.sh"
		wsplit_shell = open(split_shell, 'w')
		wsplit_shell.write("set -e\n")
		if rawfq.endswith('gz'):
			shell_line = " ".join(["gunzip -c", rawfq, "| split -l", str(maxsize), "-a 4 -", fq_prefix_name + ".1.\n"])
			wsplit_shell.write(shell_line)
			if rawfq2 != 0:
				shell_line = " ".join(["gunzip -c", rawfq2, "| split -l", str(maxsize), "-a 4 -", fq_prefix_name + ".2.\n"])
				wsplit_shell.write(shell_line)
		else:
			shell_line = " ".join(["split -l", str(maxsize), "-a 4", rawfq, fq_prefix_name + ".1.\n"])
			wsplit_shell.write(shell_line)
			if rawfq2 != 0:
				shell_line = " ".join(["split -l", str(maxsize), "-a 4", rawfq2, fq_prefix_name + ".2.\n"])
				wsplit_shell.write(shell_line)
		wsplit_shell.write("echo Done! > " + split_shell + ".sign\n")
		wsplit_shell.close()

		if os.path.isfile(split_shell + ".sign") == False:
			subprocess.call(["sh", split_shell])

		pre = fq_prefix_name + ".1.a*"
		split_fq_file_tmp_list = glob.glob(pre)

		OutFileList = os.path.join(outdir, "split_fq_barcode.txt")
		revise_sign_file = OutFileList + ".sign"
		wOutFileList = open(OutFileList, 'w')
		wrevise_sign_file = open(revise_sign_file, 'w')

		im = 0
		pool = multiprocessing.Pool(processes = process_num)

		for sp in range(0, len(split_fq_file_tmp_list)):
			SplitFqFile = fq_prefix_name + "." + str(sp) + ".fq.gz"
			SplitBarcodeFile = outdir + "/" + b + "." + str(sp) + ".barcode.gz"
			print("[ %s ] split fq and barcode: %s; %s" % (time.asctime(), SplitFqFile, SplitBarcodeFile))
			OutFileListInfo = SplitFqFile + "\t" + SplitBarcodeFile + "\n"
			wOutFileList.write(OutFileListInfo)

			Splitfq1 = split_fq_file_tmp_list[sp]
			Splitfq2 = Splitfq1.replace(".1.a", ".2.a")
			if os.path.isfile(Splitfq2) == False:
				Splitfq2 = 0

			im += 1
			if im <= process_num:
				pool.apply_async(self.sub_reviese, (Splitfq1, SplitFqFile, SplitBarcodeFile, Splitfq2))
			if im == process_num:
				pool.close()
				pool.join()

				im = 0
				pool = multiprocessing.Pool(processes = process_num)
		if im > 0:
			pool.close()
			pool.join()

		wrevise_sign_file.write("done!\n")
		wrevise_sign_file.close()
		wOutFileList.close()
		return(OutFileList)

class modify_barcode:
	def replace_barcode(self, rawfq, barcode, barcode_sam, bwa_path, bwa_parameter, barcode_ref_prefix, new_fq_prefix, samtools_path):
		outdir = os.path.dirname(barcode_sam)
		shelldir = outdir + "/shell"
		check_info(shelldir, "dir")
		tmpprefix = os.path.basename(new_fq_prefix)
		shell = shelldir + "/" + tmpprefix + ".barcode.align.sh"
		logfile = shell + ".log"
#		shell = os.path.join(shelldir, "barcode.align.sh")
		sai = barcode_sam + '.sai'
		while os.path.exists(barcode_sam) == False or os.path.exists(sai) == True:
			oshell = open(shell, 'w')
			sai = barcode_sam + '.sai'
			shell_line = " ".join([bwa_path, 'aln', barcode_ref_prefix, barcode, '-f', sai, bwa_parameter, "2>", logfile, "\n"])
			oshell.write(shell_line)
			shell_line = " ".join([bwa_path, 'samse', barcode_ref_prefix, sai, barcode, '-f', barcode_sam + " 2>>", logfile, "\nrm", sai, "\n"])
			oshell.write(shell_line)
			shell_line = " ".join(["grep -v '^@'", barcode_sam, ">", barcode_sam + ".noheader\n", "mv", barcode_sam + ".noheader", barcode_sam, "\n"])
			oshell.write(shell_line)
			oshell.close()
			print("[ %s ] processing %s\n" % (time.asctime(), shell))
			subprocess.call(["sh", shell])
			print("[ %s ] %s finished!\n" % (time.asctime(), shell))

		if os.path.exists(barcode_sam) and os.path.getsize(barcode_sam):
			pass
		else:
			sys.stderr.write("[ %s ] ERROR: %s failed, error info has been written in %s!\n" % (time.asctime(), shell, logfile))
			sys.exit(-1)

		print("\n[ %s ] modify barcode info of %s\n" % (time.asctime(), rawfq))
		samfile = open(barcode_sam, "r")

		original_barcode_count = 0
		mapped_barcode_count = 0
		corrected_barcode_count = 0

		BXfq = new_fq_prefix + ".BX.fq.gz"
		wBXfq = gzip.open(BXfq, 'wb')
		aBXfq = new_fq_prefix + ".BX_seq.fq.gz"
		waBXfq = gzip.open(aBXfq, 'wb')

		Failedfq = new_fq_prefix + ".failed.fq.gz"
		wFailedfq = gzip.open(Failedfq, 'wb')

		if rawfq.endswith('gz'):
			o_rawfq = gzip.open(rawfq, "rb")
		else:
			o_rawfq = open(rawfq, "r")

		raw_fq1_info_list = list()
		raw_fq2_info_list = list()
		clean_fq1_info_list = list()
		clean_fq2_info_list = list()
		count_info1_list = list()
		count_info2_list = list()

		for n in range(5):
			raw_fq1_info_list.append(0)
			raw_fq2_info_list.append(0)
			clean_fq1_info_list.append(0)
			clean_fq2_info_list.append(0)
			count_info1_list.append(0)
			count_info2_list.append(0)

		QUALITY_version = 33

		while True:
			barcodeinfo = samfile.readline()

			rawfq_id1 = o_rawfq.readline()
			rawfq_seq1 = o_rawfq.readline()
			rawfq_str1 = o_rawfq.readline()
			rawfq_qua1 = o_rawfq.readline()

			rawfq_id2 = o_rawfq.readline()
			rawfq_seq2 = o_rawfq.readline()
			rawfq_str2 = o_rawfq.readline()
			rawfq_qua2 = o_rawfq.readline()

			if len(barcodeinfo) == 0:
				break
			else:
				original_barcode_count += 1
				if rawfq.endswith('gz'):
					rawfq_id1 = rawfq_id1.decode().strip()
					rawfq_seq1 = rawfq_seq1.decode().strip()
					rawfq_str1 = rawfq_str1.decode().strip()
					rawfq_qua1 = rawfq_qua1.decode().strip()

					rawfq_id2 = rawfq_id2.decode().strip()
					rawfq_seq2 = rawfq_seq2.decode().strip()
					rawfq_str2 = rawfq_str2.decode().strip()
					rawfq_qua2 = rawfq_qua2.decode().strip()
				else:
					rawfq_id1 = rawfq_id1.strip()
					rawfq_seq1 = rawfq_seq1.strip()
					rawfq_str1 = rawfq_str1.strip()
					rawfq_qua1 = rawfq_qua1.strip()

					rawfq_id2 = rawfq_id2.strip()
					rawfq_seq2 = rawfq_seq2.strip()
					rawfq_str2 = rawfq_str2.strip()
					rawfq_qua2 = rawfq_qua2.strip()

				barcodeinfo = barcodeinfo.strip()
				barcodeinfo = re.split("\t", barcodeinfo)

				count_info1_list[0] = len(rawfq_seq1) - 23
				count_info2_list[0] = len(rawfq_seq2)
				count_info1_list[1] = 1
				count_info2_list[1] = 1
				count_info1_list[2] = (rawfq_seq1[23:]).count("C") + (rawfq_seq1[23:]).count("G") 
				count_info2_list[2] = rawfq_seq2.count("C") + rawfq_seq2.count("G")

				if count_info1_list[1] == 1:
					QUALITY_version = sequence_quality().check_quality_version(rawfq_qua2)

				(raw1q20, raw1q30) = sequence_quality().countfq(rawfq_qua1[23:], QUALITY_version)
				count_info1_list[3] = raw1q20
				count_info1_list[4] = raw1q30
				(raw2q20, raw2q30) = sequence_quality().countfq(rawfq_qua2, QUALITY_version)
				count_info2_list[3] = raw2q20
				count_info2_list[4] = raw2q30

				for n in range(len(count_info1_list)):
					raw_fq1_info_list[n] += count_info1_list[n]
					raw_fq2_info_list[n] += count_info2_list[n]

				id1 = re.split("\s", rawfq_id1)
				if id1[0] != ('@' + barcodeinfo[0]):
					sys.stderr.write("ERROR - %s and %s do not match!\n" % (id1[0], barcodeinfo[0]))
				else:
					barcode_ori = list(barcodeinfo[9])
					realbarcodequa = str(barcodeinfo[10])
					barcode_err = 0
					if barcodeinfo[2] == "*" or re.search('N', barcodeinfo[9]):
						barcode_err = 1
						if re.search('N', barcodeinfo[9]):
							mapped_barcode_count += 1
						continue
					md = barcodeinfo[17]
					if md.startswith("MD:Z"):
						pass
					else:
						for tmpmd in barcodeinfo[11:]:
							if tmpmd.startswith("MD:Z"):
								md = tmpmd
								
					md = re.split(":", md)
					alter = re.findall(r'\D', md[-1])
					loci = re.findall(r'\d+', md[-1])
					mapped_barcode_count += 1
					if len(alter) == 0:
						corrected_barcode_count += 1
						pass
					elif len(alter) > 1:
						barcode_err = 1
					else:
						corrected_barcode_count += 1
						for i in range(len(alter)):
							ix = int(loci[i]) + i
							barcode_ori[ix] = alter[i]

					if barcode_err == 0:
						for n in range(len(count_info1_list)):
							clean_fq1_info_list[n] += count_info1_list[n]
							clean_fq2_info_list[n] += count_info2_list[n]

						realbarcode = "".join(barcode_ori)

						newfq_id1 = (re.split("\s+", rawfq_id1))[0] + "/1" + "\tBC:Z:" + realbarcode
						newfq_seq1 = rawfq_seq1[23:]
						newfq_qua1 = rawfq_qua1[23:]
						newfq_id2 = (re.split("\s+", rawfq_id2))[0] + "/2" + "\tBC:Z:" + realbarcode
						newfq_seq2 = str(rawfq_seq2)
						newfq_qua2 = str(rawfq_qua2)

						BXout = "\n".join([newfq_id1, newfq_seq1, rawfq_str1, newfq_qua1, newfq_id2, newfq_seq2, rawfq_str2, newfq_qua2, ""])
						wBXfq.write(BXout.encode())

						anewfq_id1 = rawfq_id1
						anewfq_seq1 = realbarcode + rawfq_seq1[16:]
						anewfq_qua1 = str(rawfq_qua1)
						anewfq_id2 = rawfq_id2
						anewfq_seq2 = str(rawfq_seq2)
						anewfq_qua2 = str(rawfq_qua2)
						aBXout = "\n".join([anewfq_id1, anewfq_seq1, rawfq_str1, anewfq_qua1, anewfq_id2, anewfq_seq2, rawfq_str2, anewfq_qua2, ""])
						waBXfq.write(aBXout.encode())

					else:
						Failedout = "\n".join([rawfq_id1, rawfq_seq1, rawfq_str1, rawfq_qua1, rawfq_id2, rawfq_seq2, rawfq_str2, rawfq_qua2, ""])
						wFailedfq.write(Failedout.encode())
		wFailedfq.close()
		wBXfq.close()
		waBXfq.close()

		print("[ %s ] barcode info of %s has been modified\n" % (time.asctime(), rawfq))
		return(raw_fq1_info_list, raw_fq2_info_list, clean_fq1_info_list, clean_fq2_info_list, count_info1_list, count_info2_list, original_barcode_count, mapped_barcode_count, corrected_barcode_count)

class OUTERSOFT:
	def FastQC(self, fq, fastqc_path, statdir):
		tmpshell = os.path.join(os.path.dirname(statdir), "shell/fastqc.sh")
		logfile = tmpshell + ".log"
		wtmpshell = open(tmpshell, 'w')
		shell_line = " ".join([fastqc_path, "-f fastq -o", statdir, fq, "2>", logfile, "&\n"])
		wtmpshell.write(shell_line)
		shell_line = "wait\necho Done\n"
		wtmpshell.write(shell_line)
		wtmpshell.close()

		subprocess.call(["sh", tmpshell])

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
	purification_usage = \
	'''
	modify the barcode sequence, purify and split merged fastq file into two files
	Version: 1.0.0
	Dependents: Python (>=3.0), BWA, SAMtools
	Last Updated Date: 2017-11-14

	Usage: python clean_fq.py <options>

	Basic options:
		-i --input, the input path of compressed or uncompressed fastqs; would be treated as the 1st file of paired fastq files if -2 also in use.
		-2 --input2, the input path of compressed or uncompressed fastqs; would be treated as the 2nd file of paired fastq files if used.
		-o --outputdir, the path of output
		-c --config, the path of configuration file [default: ./config/Basic.config]
		-h --help, print help info
	Advanced options:
		-z --splitsize, the amount of reads of splited fastq, reads num = split_size / 8 [default: 15000000 lines, 1875000 read-pairs, compressed file size: ~300M]
		-n --process_num, the max amout of sub-processing [default: 4]

	'''
	print(purification_usage)

if __name__ == '__main__':
	if len(sys.argv) < 5:
		usage()
		sys.exit(-1)

	inputfq = None
	inputfq2 = 0
	outputdir = None
	ConfigFile = None
	split_size = 15000000
	process_num = 4
	opts, args = getopt.gnu_getopt(sys.argv[1:], 'i:2:o:c:z:n:h', ['input', 'input2', 'outputdir', 'config', 'splitsize', 'process_num', 'help'])
	for o, a in opts:
		if o == '-i' or o == '--input':
			inputfq = a
		if o == '-2' or o == '--input2':
			inputfq2 = a
		if o == '-o' or o == '--outputdir':
			outputdir = a
		if o == '-c' or o == '--config':
			ConfigFile = a
		if o == '-z' or o == '--splitsize':
			split_size = int(a)
		if o == '-n' or o == '--process_num':
			process_num = int(a)
		if o == '-h' or o == '--help':
			usage()
			sys.exit(-1)

	prefix_name = None
	if os.path.isfile(inputfq):
		b = os.path.basename(inputfq)
		if b[(len(b)-9):] == ".fastq.gz":
			b = b[0:(len(b)-9)]
		elif b[(len(b)-6):] == ".fq.gz" or b[(len(b)-6):] == ".fastq":
			b = b[0:(len(b)-6)]
		elif b[(len(b)-3):] == ".fq":
			b = b[0:(len(b)-3)]

		if b.endswith("_1"):
			b = b[0:(len(b)-2)]
		prefix_name = outputdir + "/" + b

		statout = prefix_name + ".fq_statistics.txt"
		newallBXfqfile = prefix_name + ".barcode_id.clean.fq.gz"
		fastqc_file = outputdir + "/fastqc/" + b + ".barcode_id.clean_fastqc.html"
		if os.path.exists(statout) and os.path.getsize(statout) and os.path.exists(newallBXfqfile) and os.path.getsize(newallBXfqfile) and os.path.exists(fastqc_file) and os.path.getsize(fastqc_file):
			sys.stderr.write("[ %s ] fastq pre-processing has been done!\n\n" % time.asctime())
			sys.exit(1)
	else:
		sys.stderr.write("[ %s ] ERROR: %s does not exist!\n" % (time.asctime(), inputfq))
		sys.exit(-1)

	if ConfigFile == None or os.path.exists(ConfigFile) == False:
		sys.stderr.write("configuration file has not been provided or does not exist, Please create it using 'python LRTK-SEQ.py Config'\n")
		sys.exit(-1)

	G = baseinfo()
	G.get_config(ConfigFile)

	R = extract_barcode()
	M = modify_barcode()
	SplitFqBarcodeFileList = outputdir + '/split_fq_barcode.txt'
	revise_sign = SplitFqBarcodeFileList + ".sign"
	split_done = 1
	if os.path.exists(revise_sign) and os.path.getsize(revise_sign) and os.path.exists(SplitFqBarcodeFileList):
		print(SplitFqBarcodeFileList)
		rSplitFqBarcodeFileList = open(SplitFqBarcodeFileList, 'r')
		for splitinfo in rSplitFqBarcodeFileList:
			(SplitFqFile, SplitBarcodeFile) = re.split("\t", splitinfo.strip())
			if os.path.exists(SplitFqFile) == False or os.path.exists(SplitBarcodeFile) == False:
				split_done = 0
		rSplitFqBarcodeFileList.close()
	else:
		split_done = 0
	
	print("[ %s ] extract original barcode info from %s\n" % (time.asctime(), inputfq))
	SplitFqBarcodeFileList = R.revise(inputfq, prefix_name, outputdir, split_size, split_done, inputfq2)
	print("[ %s ] original barcode had been extracted, and raw fq file had been splited, output files has been listed in : %s\n" % (time.asctime(), SplitFqBarcodeFileList))

	SplitCleanFqFileList = outputdir + '/split_clean_fq.txt'
	replace_barcode_sign = SplitCleanFqFileList + ".sign"
	samtools = G.Samtools()

	fq_statistics_result = list()
	BXfqfiletmplist = list()
	BXModifiedfqfiletmplist = list()
	failedfqtmplist = list()
	otherfiletmplist = list()

	if os.path.exists(replace_barcode_sign) and os.path.getsize(replace_barcode_sign) != 0:
		print("[ %s ] barcode info has been modified!\n" % time.asctime())
	else:
		barcoderef = G.Barcode()
		barcoderef_index = barcoderef + ".ann"
		if os.path.isfile(barcoderef_index):
			pass
		else:
			sys.stderr.write("[ %s ] ERROR: index of %s does not exist!\n" % (time.asctime(), barcoderef))
			sys.exit(-1)

		bracode_aln_parameter = G.Barcode_aln_par()
		bwa_path = G.Bwa()
		rSplitFqBarcodeFileList = open(SplitFqBarcodeFileList, 'r')
		im = 0
		pool = multiprocessing.Pool(processes = process_num)

		wSplitCleanFqFileList = open(SplitCleanFqFileList, 'w')
		wreplace_barcode_sign = open(replace_barcode_sign, 'w')
		for splitinfo in rSplitFqBarcodeFileList:
			(split_fq_file, split_barcode_file) = re.split("\t", splitinfo.strip())
			if os.path.exists(split_fq_file) == False or os.path.getsize(split_fq_file) == 0 or os.path.exists(split_barcode_file) == False or os.path.getsize(split_barcode_file) == 0:
				sys.stderr.write("[ %s ] ERROR: %s or %s does not exist or is empty!\n" % (time.asctime(), split_fq_file, split_barcode_file))
				sys.exit(-1)
			barcodesam = (str(split_barcode_file))[0:(len(split_barcode_file)-2)] + "sam"
			otherfiletmplist.append(barcodesam)
			otherfiletmplist.append(split_fq_file)
			otherfiletmplist.append(split_barcode_file)
			newprefix = (str(split_barcode_file))[0:(len(split_barcode_file)-11)]

			im += 1
			if im <= process_num:
				fq_statistics_result.append(pool.apply_async(M.replace_barcode, (split_fq_file, split_barcode_file, barcodesam, bwa_path, bracode_aln_parameter, barcoderef, newprefix, samtools)))
			if im == process_num:
				pool.close()
				pool.join()

				im = 0
				pool = multiprocessing.Pool(processes = process_num)

			BXfqfile = newprefix + ".BX.fq.gz"
			BXfqfiletmplist.append(BXfqfile)
			failedfqfile = newprefix + ".failed.fq.gz"
			failedfqtmplist.append(failedfqfile)

			SplitCleanFqFile = BXfqfile + "\n"
			wSplitCleanFqFileList.write(SplitCleanFqFile)
		wSplitCleanFqFileList.close()

		if im > 0:
			pool.close()
			pool.join()

		statout = prefix_name + ".fq_statistics.txt"
		wstatout = open(statout, 'w')
		raw_fq1_info_list = [0, 0, 0, 0, 0]
		raw_fq2_info_list = [0, 0, 0, 0, 0]
		clean_fq1_info_list = [0, 0, 0, 0, 0]
		clean_fq2_info_list = [0, 0, 0, 0, 0]
		count_info1_list = [0, 0, 0, 0, 0]
		count_info2_list = [0, 0, 0, 0, 0]
		original_barcode_count = 0
		mapped_barcode_count = 0
		corrected_barcode_count = 0
		for res in fq_statistics_result:
			for ai in range(len(raw_fq1_info_list)):
				raw_fq1_info_list[ai] += int(res.get()[0][ai])
				raw_fq2_info_list[ai] += int(res.get()[1][ai])
				clean_fq1_info_list[ai] += int(res.get()[2][ai])
				clean_fq2_info_list[ai] += int(res.get()[3][ai])
				count_info1_list[ai] += int(res.get()[4][ai])
				count_info2_list[ai] += int(res.get()[5][ai])
			original_barcode_count += int(res.get()[6])
			mapped_barcode_count += int(res.get()[7])
			corrected_barcode_count += int(res.get()[8])
		statlog = "#Sequencing quality report\n" + os.path.basename(prefix_name) + "\tYield(bp)\tRead_length(bp)\tGC(%)\tQ20(%)\tQ30(%)\n"
		wstatout.write(statlog)
		raw_fq1_info_list[1] = int(raw_fq1_info_list[0] / raw_fq1_info_list[1])
		raw_fq2_info_list[1] = int(raw_fq2_info_list[0] / raw_fq2_info_list[1])
		clean_fq1_info_list[1] = int(clean_fq1_info_list[0] / clean_fq1_info_list[1])
		clean_fq2_info_list[1] = int(clean_fq2_info_list[0] / clean_fq2_info_list[1])
		for n in range(2, 5):
			raw_fq1_info_list[n] = round(100.00 * raw_fq1_info_list[n] / raw_fq1_info_list[0], 3)
			raw_fq2_info_list[n] = round(100.00 * raw_fq2_info_list[n] / raw_fq2_info_list[0], 3)
			clean_fq1_info_list[n] = round(100.00 * clean_fq1_info_list[n] / clean_fq1_info_list[0], 3)
			clean_fq2_info_list[n] = round(100.00 * clean_fq2_info_list[n] / clean_fq2_info_list[0], 3)
		for n in range(len(raw_fq1_info_list)):
			raw_fq1_info_list[n] = str(raw_fq1_info_list[n])
			raw_fq2_info_list[n] = str(raw_fq2_info_list[n])
			clean_fq1_info_list[n] = str(clean_fq1_info_list[n])
			clean_fq2_info_list[n] = str(clean_fq2_info_list[n])
		statlog = "raw fastq1:\t" + "\t".join(raw_fq1_info_list) + "\n"
		wstatout.write(statlog)
		statlog = "raw fastq2:\t" + "\t".join(raw_fq2_info_list) + "\n"
		wstatout.write(statlog)
		statlog = "clean fastq1:\t" + "\t".join(clean_fq1_info_list) + "\n"
		wstatout.write(statlog)
		statlog = "clean fastq2:\t" + "\t".join(clean_fq2_info_list) + "\n\n"
		wstatout.write(statlog)
		statlog = "#Barcode correction report\namount of barcode in original fastq:\t" + str(original_barcode_count) + "\t100.0%\n"
		wstatout.write(statlog)
		qrate = round(100.0 * mapped_barcode_count / original_barcode_count, 2)
		statlog = "amount of barcode that mapped to the white list:\t" + str(mapped_barcode_count) + "\t" + str(qrate) + "%\n"
		wstatout.write(statlog)
		qrate = round(100.0 * corrected_barcode_count / original_barcode_count, 2)
		statlog = "amount of barcode that perfectly mapped to the white list after corrected:\t" + str(corrected_barcode_count) + "\t" + str(qrate) + "%\n"
		wstatout.write(statlog)
		wstatout.close()

		wreplace_barcode_sign.write("done!\n")
		wreplace_barcode_sign.close()

	tmpdatadir = outputdir + "/tmp"
	check_info(tmpdatadir, "dir")

	newallBXfqfile = prefix_name + ".barcode_id.clean.fq.gz"
	final_fq_list = BXfqfiletmplist
	newallBX_add_fqfile = prefix_name + ".barcode_seq.clean.fq.gz"

	newallfailedfqfile = prefix_name + ".failed.fq.gz"
	tmpshell = newallBXfqfile + ".cat.sh"
	wtmpshell = open(tmpshell, 'w')
	Split_Final_FQfile_List = outputdir + "/Split_final_fq.txt"
	wSplit_Final_FQfile_List = open(Split_Final_FQfile_List, 'w')
	for fi in range(0,len(final_fq_list)):
		wSplit_Final_FQfile_List.write(final_fq_list[fi] + "\n")
		if fi == 0:
			shell_line = " ".join(["cat", final_fq_list[fi], ">", newallBXfqfile + "\n"])
			wtmpshell.write(shell_line)
			seq_fq = (final_fq_list[fi]).replace("BX.fq.gz", "BX_seq.fq.gz")
			shell_line = " ".join(["cat", seq_fq, ">", newallBX_add_fqfile + "\n"])
			wtmpshell.write(shell_line)
			shell_line = " ".join(["cat", failedfqtmplist[fi], ">", newallfailedfqfile + "\n", "mv", failedfqtmplist[fi], tmpdatadir + "\n"])
			wtmpshell.write(shell_line)
		else:
			shell_line = " ".join(["cat", final_fq_list[fi], ">>", newallBXfqfile + "\n"])
			wtmpshell.write(shell_line)
			seq_fq = (final_fq_list[fi]).replace("BX.fq.gz", "BX_seq.fq.gz")
			shell_line = " ".join(["cat", seq_fq, ">>", newallBX_add_fqfile + "\n"])
			wtmpshell.write(shell_line)
			shell_line = " ".join(["cat", failedfqtmplist[fi], ">>", newallfailedfqfile + "\n", "mv", failedfqtmplist[fi], tmpdatadir + "\n"])
			wtmpshell.write(shell_line)
	wSplit_Final_FQfile_List.close()
	othertmpfile = " ".join(otherfiletmplist)
	shell_line = " ".join(["mv", outputdir + "/split*txt*", othertmpfile, tmpdatadir + "\n"])
	wtmpshell.write(shell_line)
	wtmpshell.close()
	subprocess.call(["sh", tmpshell])
#	subprocess.call(["mv", tmpshell, tmpdatadir])

	O = OUTERSOFT()
	print("[ %s ] processing fastq QC ... \n" % time.asctime())
	statdir = outputdir + "/fastqc"
	check_info(statdir, "dir")
	Fastqc_path = G.Fastqc()
	O.FastQC(newallBXfqfile, Fastqc_path, statdir)
	print("[ %s ] fastq QC has finished, results have been writeen to %s \n\n" % (time.asctime(), statdir))

#	subprocess.call(["mv", newallfq1file, newallfq2file, SplitFqBarcodeFileList, tmpdatadir])
#### encode
