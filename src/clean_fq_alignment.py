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

def find_samtools_version(self, samtools_path, outdir = './'):
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
	def revise(self, rawfq, outdir, maxsize, _split_done, rawfq2 = 0):
		(rawfq_path, rawfq_name) = os.path.split(rawfq)
		b = str(rawfq_name)
		if b[(len(b)-9):] == ".fastq.gz":
			b = b[0:(len(b)-9)]
		elif b[(len(b)-6):] == ".fq.gz" or b[(len(b)-6):] == ".fastq":
			b = b[0:(len(b)-6)]
		elif b[(len(b)-3):] == ".fq":
			b = b[0:(len(b)-3)]

		if b.endswith("_1"):
			b = b[0:(len(b)-2)]

		if _split_done == 1:
			fq_prefix_name = outdir + "/" + b
			OutFileList = os.path.join(outdir, "split_fq_barcode.txt")
			return(fq_prefix_name, OutFileList)

		if rawfq.endswith('gz'):
			o_rawfq = gzip.open(rawfq, 'rb')
			if rawfq2 != 0:
				o_rawfq2 = gzip.open(rawfq2, 'rb')
		else:
			o_rawfq = open(rawfq, 'r')
			if rawfq2 != 0:
				o_rawfq2 = open(rawfq2, 'r')

		OutFileList = os.path.join(outdir, "split_fq_barcode.txt")
		revise_sign_file = OutFileList + ".sign"
		wOutFileList = open(OutFileList, 'w')
		wrevise_sign_file = open(revise_sign_file, 'w')
		sp = 0
		SplitFqFile = outdir + "/" + b + "." + str(sp) + ".fq.gz"
		SplitBarcodeFile = outdir + "/" + b + "." + str(sp) + ".barcode.gz"
		print("[ %s ] split fq and barcode: %s; %s" % (time.asctime(), SplitFqFile, SplitBarcodeFile))
		OutFileListInfo = SplitFqFile + "\t" + SplitBarcodeFile + "\n"
		wOutFileList.write(OutFileListInfo)

		wSplitFqFile = gzip.open(SplitFqFile, 'wb')
		wSplitBarcodeFile = gzip.open(SplitBarcodeFile, 'wb')

		SplitSize = 0	
		while True:
			rawfq_id1 = o_rawfq.readline()
			rawfq_seq1 = o_rawfq.readline()
			rawfq_str1 = o_rawfq.readline()
			rawfq_qua1 = o_rawfq.readline()

			if rawfq2 == 0:
				rawfq_id2 = o_rawfq.readline()
				rawfq_seq2 = o_rawfq.readline()
				rawfq_str2 = o_rawfq.readline()
				rawfq_qua2 = o_rawfq.readline()
			else:
				rawfq_id2 = o_rawfq2.readline()
				rawfq_seq2 = o_rawfq2.readline()
				rawfq_str2 = o_rawfq2.readline()
				rawfq_qua2 = o_rawfq2.readline()

			SplitSize += 8

			if len(rawfq_id1) == 0:
				break
			else:
				if rawfq.endswith('gz'):
					rawfq_id1 = rawfq_id1.decode()
					rawfq_seq1 = rawfq_seq1.decode()
					rawfq_str1 = rawfq_str1.decode()
					rawfq_qua1 = rawfq_qua1.decode()
					rawfq_id2 = rawfq_id2.decode()
					rawfq_seq2 = rawfq_seq2.decode()
					rawfq_str2 = rawfq_str2.decode()
					rawfq_qua2 = rawfq_qua2.decode()

				barcode_line = rawfq_id1 + rawfq_seq1[0:16] + '\n' + rawfq_str1 + rawfq_qua1[0:16] + "\n"
				fq_line = rawfq_id1 + rawfq_seq1 + rawfq_str1 + rawfq_qua1 + rawfq_id2 + rawfq_seq2 + rawfq_str2 + rawfq_qua2

				if SplitSize <= maxsize:
					wSplitFqFile.write(fq_line.encode())
					wSplitBarcodeFile.write(barcode_line.encode())
				else:
					wSplitFqFile.close()
					wSplitBarcodeFile.close()
					sp += 1
					SplitSize = 8
					SplitFqFile = outdir + "/" + b + "." + str(sp) + ".fq.gz"
					SplitBarcodeFile = outdir + "/" + b + "." +  str(sp) + ".barcode.gz"
					print("[ %s ] split fq and barcode: %s; %s" % (time.asctime(), SplitFqFile, SplitBarcodeFile))
					OutFileListInfo = SplitFqFile + "\t" + SplitBarcodeFile + "\n"
					wOutFileList.write(OutFileListInfo)

					wSplitFqFile = gzip.open(SplitFqFile, 'wb')
					wSplitBarcodeFile = gzip.open(SplitBarcodeFile, 'wb')

					wSplitFqFile.write(fq_line.encode())
					wSplitBarcodeFile.write(barcode_line.encode())
		print()
		o_rawfq.close()
		wOutFileList.close()
		wSplitFqFile.close()
		wSplitBarcodeFile.close()
		fq_prefix_name = outdir + "/" + b
		wrevise_sign_file.write("done!\n")
		wrevise_sign_file.close()
		return(fq_prefix_name, OutFileList)

class modify_barcode:
	def replace_barcode(self, rawfq, barcode, barcode_sam, bwa_path, bwa_parameter, barcode_ref_prefix, new_fq_prefix, barcode_sequence, samtools_path):
		outdir = os.path.dirname(barcode_sam)
		shelldir = outdir + "/shell"
		check_info(shelldir, "dir")
		tmpprefix = os.path.basename(new_fq_prefix)
		shell = shelldir + "/" + tmpprefix + ".barcode.align.sh"
		logfile = shell + ".log"
#		shell = os.path.join(shelldir, "barcode.align.sh")
		if os.path.exists(barcode_sam) == False or os.path.getsize(barcode_sam) == 0:
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

		if barcode_sequence != None:
			fq2bam = new_fq_prefix + ".fq.sam"
			wfq2bam = open(fq2bam, 'w')
			fq2baminfo = "@HD\tVN:1.3\n@SQ\tSN:chr1\tLN:248956422\n"
			wfq2bam.write(fq2baminfo)

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

						if barcode_sequence != None:
							readid = realbarcode + "_" + ((re.split("\s+", rawfq_id1))[0]).replace('@', '')
							readlength1 = str(len(newfq_seq1)) + "M"
							readlength2 = str(len(newfq_seq2)) + "M"
							fq2baminfo = "\t".join([readid, "*", "*", "*", "*", readlength1, "*", "*", "*", newfq_seq1, newfq_qua1]) + "\n"
							wfq2bam.write(fq2baminfo)
							fq2baminfo = "\t".join([readid, "*", "*", "*", "*", readlength2, "*", "*", "*", newfq_seq2, newfq_qua2]) + "\n"
							wfq2bam.write(fq2baminfo)
					else:
						Failedout = "\n".join([rawfq_id1, rawfq_seq1, rawfq_str1, rawfq_qua1, rawfq_id2, rawfq_seq2, rawfq_str2, rawfq_qua2, ""])
						wFailedfq.write(Failedout.encode())
		wFailedfq.close()
		wBXfq.close()

		if barcode_sequence != None:
			wfq2bam.close()
			shell = shelldir + "/" + tmpprefix + ".sortbam.sh"
			wshell = open(shell, 'w')
			shell_line = " ".join([samtools_path, "view -S -h -b", fq2bam, ">", fq2bam + ".bam"]) + "\n"
			wshell.write(shell_line)
			shell_line = "rm " + fq2bam + "\n"
			wshell.write(shell_line)
			wshell.close()
			subprocess.call(["sh", shell])

		print("[ %s ] barcode info of %s has been modified\n" % (time.asctime(), rawfq))
		return(raw_fq1_info_list, raw_fq2_info_list, clean_fq1_info_list, clean_fq2_info_list, count_info1_list, count_info2_list, original_barcode_count, mapped_barcode_count, corrected_barcode_count)

	def get_new_barcode(self, barcode_num, _start):
		bin_list = (16777216, 4194304, 1048576, 262144, 65536, 16384, 4096, 1024, 256, 64, 16, 4)
		base_dict = {0:"A", 1:"C", 2:"G", 3:"T"}
		n = barcode_num
		plist = list()
		for i in range(0, 12):
			plist.append(0)
		if barcode_num >= bin_list[_start-1]:
			_start = _start - 1
		for b in range(_start, 12):
			com = n // bin_list[b]
			n = n % bin_list[b]
			p = 12 - b
			plist[p] = com
		convert = base_dict[n]
		for r in range(1,12):
			convert = base_dict[plist[r]] + convert
		return(_start, convert)

	def generate_barcode(self, split_fq_sam_file_List, outdir, samtools_path, split_modify_fq_file_List, file_sign, maxsize, barcode_sequence):
		fqsamlist = list()
		rsplit_fq_sam_file_List = open(split_fq_sam_file_List, 'r')
		for split_fq_sam_file_List_info in rsplit_fq_sam_file_List:
			fqsamlist.append((re.split("\t", split_fq_sam_file_List_info.strip()))[1])
			otherfiletmplist.append((re.split("\t", split_fq_sam_file_List_info.strip()))[1])
		mergedFqBam = fqsamlist[0].replace("0.fq.sam.bam", "fq.merged.bam")
		shelldir = outdir + "/shell"
		check_info(shelldir, "dir")
		shell = shelldir + "/merge_fqbam.sh"
		wshell = open(shell, 'w')
		if len(fqsamlist) > 1:
			allfqsam = " ".join(fqsamlist)
			shell_line = " ".join([samtools_path, "merge -f", mergedFqBam, allfqsam]) + "\n"
		else:
			shell_line = " ".join(["ln -sf", fqsamlist[0], mergedFqBam]) + "\n"
		wshell.write(shell_line)
		sv = find_samtools_version(samtools_path, shelldir)
		sortedFqBamprefix = mergedFqBam.replace("fq.merged.bam", "fq.merged.sorted")
		sortedFqBam = sortedFqBamprefix + ".bam"
		if sv == 0:
			shell_line = " ".join([samtools_path, "sort -n -m 1G", mergedFqBam, sortedFqBamprefix]) + "\n"
		else:
			shell_line = " ".join([samtools_path, "sort -n -m 1G -o", sortedFqBamprefix + ".bam", mergedFqBam]) + "\n"
		wshell.write(shell_line)
		shell_line = " ".join([samtools_path, "index", sortedFqBam]) + "\n"
		wshell.write(shell_line)
		wshell.close()
		subprocess.call(["sh", shell])
		new_fq_prefix = sortedFqBam.replace("fq.merged.sorted.bam", "")
		otherfiletmplist.append(mergedFqBam)
		otherfiletmplist.append(sortedFqBam)
		otherfiletmplist.append(sortedFqBam + ".bai")

		wsplit_modify_fq_file_List = open(split_modify_fq_file_List, 'w')
		SplitSize = 0
		rsortedFqBam = pysam.AlignmentFile(sortedFqBam, 'rb')
		readid = "N"
		barcodeid = 0
		barcode_marker = "N"
		new_barcode_sequence = "N"
		start = 11

		s = 0
		split_modify_fq_file = fqsamlist[0].replace("0.fq.sam.bam", "") + str(s) + ".BX.modified.fq.gz"
		wsplit_modify_fq_file_List.write(split_modify_fq_file + "\n")
		wsplit_modify_fq_file = gzip.open(split_modify_fq_file, 'wb')
		for FqBaminfo in rsortedFqBam:
			(real_barcode, real_readid) = re.split("_", FqBaminfo.query_name)
			if FqBaminfo.query_name != readid:
				readid = FqBaminfo.query_name
				if barcode_marker != real_barcode:
					barcode_marker = real_barcode
					(start, new_barcode_suffix) = self.get_new_barcode(barcodeid, start)
					barcodeid += 1
					new_barcode_sequence = barcode_sequence + new_barcode_suffix
				real_readid = '@' + real_readid + "/1\tBC:Z:" + new_barcode_sequence
			else:
				real_readid = '@' + real_readid + "/2\tBC:Z:" + new_barcode_sequence
			SplitSize += 4
			complete_read_info = "\n".join([real_readid, FqBaminfo.query_sequence, "+", pysam.qualities_to_qualitystring(FqBaminfo.query_qualities)]) + "\n"
			if SplitSize > maxsize:
				wsplit_modify_fq_file.close()
				s += 1
				SplitSize = 4
				split_modify_fq_file = fqsamlist[0].replace("0.fq.sam.bam", "") + str(s) + ".BX.modified.fq.gz"
				wsplit_modify_fq_file_List.write(split_modify_fq_file + "\n")
				wsplit_modify_fq_file = gzip.open(split_modify_fq_file, 'wb')
			wsplit_modify_fq_file.write(complete_read_info.encode())
		rsortedFqBam.close()
		wsplit_modify_fq_file_List.close()
		wsplit_modify_fq_file.close()

		wfile_sign = open(file_sign, 'w')
		wfile_sign.write("done!\n")
		wfile_sign.close

		return(split_modify_fq_file_List)

class OUTERSOFT:
	def bwa_fq(self, fq1, outsam, bwa_path, bwa_aln_parameter, bwa_sam_parameter, RGinfo, samtools_path):
		print("[ %s ] fq alignment starts: %s \n" % (time.asctime(), fq1))
		fqprefix = (os.path.basename(fq1)).replace(".gz", "")
		outdir = os.path.dirname(fq1)
		shelldir = os.path.join(outdir, "shell")
		check_info(shelldir, "dir")
		shell = os.path.join(shelldir, fqprefix + ".aln.sh")
		logfile = shell + ".log"
		oshell = open(shell, 'w')
		if RGinfo != 0:
			RGinfo = '"' + RGinfo + '"'
			shell_line = " ".join([bwa_path, 'mem -p -C', '-r', RGinfo, ref_path, fq1, ">", outsam, '2>', logfile + "\n"])
		else:
			shell_line = " ".join([bwa_path, 'mem -p -C', ref_path, fq1, ">", outsam, '2>', logfile + "\n"])
		oshell.write(shell_line)
		outbam = outsam.replace(".sam", ".bam")
		shell_line = " ".join([samtools_path, "view -S -h -F 256 -F 2048 -b", outsam, ">", outbam]) + "\n"
		oshell.write(shell_line)
		oshell.close()

		subprocess.call(["sh", shell])
		print("[ %s ] fq alignment finish: %s \n" % (time.asctime(), outbam))

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

	Usage: python clean_fq_alignment.py <options>

	Basic options:
		-i --input, the input path of compressed or uncompressed fastqs; would be treated as the 1st file of paired fastq files if -2 also in use.
		-2 --input2, the input path of compressed or uncompressed fastqs; would be treated as the 2nd file of paired fastq files if used.
		-o --outputdir, the path of output
		-c --config, the path of configuration file [default: ./config/Basic.config]
		-r --RG, RG info for alignment
		-h --help, print help info
	Advanced options:
		-z --splitsize, the amount of reads of splited fastq, reads num = split_size / 8 [default: 15000000 lines, 1875000 read-pairs, compressed file size: ~300M]
		-b --barcode, the first four bases of barcode sequence [default: None]
		-n --process_num, the max amout of sub-processing [default: 4]

	'''
	print(purification_usage)

if __name__ == '__main__':
	if len(sys.argv) < 5:
		usage()
		sys.exit(-1)

	RGInfo = 0
	inputfq = None
	inputfq2 = 0
	outputdir = None
	ConfigFile = None
	split_size = 15000000
	process_num = 4
	barcode_seq = None
	opts, args = getopt.gnu_getopt(sys.argv[1:], 'i:2:o:c:b:z:n:r:h', ['input', 'input2', 'outputdir', 'config', 'barcode', 'splitsize', 'process_num', 'RG', 'help'])
	for o, a in opts:
		if o == '-i' or o == '--input':
			inputfq = a
		if o == '-2' or o == '--input2':
			inputfq2 = a
		if o == '-o' or o == '--outputdir':
			outputdir = a
		if o == '-c' or o == '--config':
			ConfigFile = a
		if o == '-b' or o == '--barcode':
			barcode_seq = str(a)
		if o == '-z' or o == '--splitsize':
			split_size = int(a)
		if o == '-n' or o == '--process_num':
			process_num = int(a)
		if o == '-r' or o == '--RG':
			RGInfo = a
		if o == '-h' or o == '--help':
			usage()
			sys.exit(-1)

	if os.path.isfile(inputfq):
		pass
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
	(prefix_name, SplitFqBarcodeFileList) = R.revise(inputfq, outputdir, split_size, split_done, inputfq2)
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
				fq_statistics_result.append(pool.apply_async(M.replace_barcode, (split_fq_file, split_barcode_file, barcodesam, bwa_path, bracode_aln_parameter, barcoderef, newprefix, barcode_seq, samtools)))
			if im == process_num:
				pool.close()
				pool.join()

				im = 0
				pool = multiprocessing.Pool(processes = process_num)

			BXfqfile = newprefix + ".BX.fq.gz"
			failedfqfile = newprefix + ".failed.fq.gz"
			failedfqtmplist.append(failedfqfile)
			fq2bamfile = newprefix + ".fq.sam.bam"

			SplitCleanFqFile = BXfqfile + "\n"
			if barcode_seq != None:
				SplitCleanFqFile = "\t".join([BXfqfile, fq2bamfile]) + "\n"
				otherfiletmplist.append(BXfqfile)
			else:
				BXfqfiletmplist.append(BXfqfile)
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

	SplitModifyFqFileList = outputdir + '/split_modify_fq.txt'
	modify_fq_sign = SplitModifyFqFileList + ".sign"
	if barcode_seq != None and (os.path.exists(modify_fq_sign) == False or os.path.getsize(modify_fq_sign) == 0):
		SplitModifyFqFileList = M.generate_barcode(SplitCleanFqFileList, outputdir, samtools, SplitModifyFqFileList, modify_fq_sign, split_size, barcode_seq)
		if os.path.exists(SplitModifyFqFileList) == False or os.path.getsize(SplitModifyFqFileList) == 0:
			print("[ %s ] ERROR: %s does not exist!\n" % (time.asctime(), SplitModifyFqFileList))
			sys.exit(-1)
		rSplitModifyFqFileList = open(SplitModifyFqFileList, 'r')
		for SplitModifyFqFileListinfo in rSplitModifyFqFileList:
			BXModifiedfqfiletmplist.append(SplitModifyFqFileListinfo.strip())
		rSplitModifyFqFileList.close()

	O = OUTERSOFT()
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

	SplitSamFileList = outputdir + '/split_bam.txt'
	aln_sign = SplitSamFileList + ".sign"
	if os.path.exists(aln_sign) == False or os.path.getsize(aln_sign) == 0:
		wSplitSamFileList = open(SplitSamFileList, 'w')
		im = 0
		pool = multiprocessing.Pool(processes = process_num)
		if barcode_seq != None:
			SplitCleanFqFileList = SplitModifyFqFileList
		rSplitCleanFqFileList = open(SplitCleanFqFileList, 'r')
		for SplitCleanFqFileInfo in rSplitCleanFqFileList:
			SplitCleanFqFileInfolist = re.split("\t", SplitCleanFqFileInfo.strip())
			if os.path.exists(SplitCleanFqFileInfolist[0]) == False or os.path.getsize(SplitCleanFqFileInfolist[0]) == 0:
				sys.stderr.write("[ %s ] ERROR: %s dose not exist or is empty!\n" % (time.asctime(), SplitCleanFqFileInfolist[0]))
				sys.exit(-1)
			fqsam = (str(SplitCleanFqFileInfolist[0]))[0:(len(SplitCleanFqFileInfolist[0])-6)] + ".sam"
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
#	subprocess.call(["rm", SplitCleanFqFileList])
	
	SplitBamFilelist = list()
	rSplitSamFileList = open(SplitSamFileList, 'r')
	for SplitSamFileListInfo in rSplitSamFileList:
		SplitBamFilelist.append(SplitSamFileListInfo.strip())
		otherfiletmplist.append(SplitSamFileListInfo.strip())
	allSplitBamFile = " ".join(SplitBamFilelist)
	tmpmergeshell = SplitBamFilelist[0] + ".merge.sh"
	wtmpmergeshell = open(tmpmergeshell, 'w')
	MergeBamFile = prefix_name + ".sorted.bam"
	if os.path.isfile(MergeBamFile):
		subprocess.call(["rm", MergeBamFile])
	if len(SplitBamFilelist) > 1:
		shell_line = " ".join([samtools, "merge -f", MergeBamFile, allSplitBamFile]) + "\n"
	else:
		shell_line = " ".join(["cp", SplitBamFilelist[0], MergeBamFile]) + "\n"
	wtmpmergeshell.write(shell_line)
	wtmpmergeshell.close()
	subprocess.call(["sh", tmpmergeshell])
	otherfiletmplist.append(tmpmergeshell)

	tmpdatadir = outputdir + "/tmp"
	check_info(tmpdatadir, "dir")

	newallBXfqfile = prefix_name + ".BX.modified.fq.gz"
	final_fq_list = BXModifiedfqfiletmplist
	if barcode_seq == None:
		final_fq_list = BXfqfiletmplist
	newallfailedfqfile = prefix_name + ".failed.fq.gz"
	tmpshell = newallBXfqfile + ".cat.sh"
	wtmpshell = open(tmpshell, 'w')
	for fi in range(0,len(final_fq_list)):
		if fi == 0:
			shell_line = " ".join(["cat", final_fq_list[fi], ">", newallBXfqfile + "\n", "mv", final_fq_list[fi], tmpdatadir + "\n"])
			wtmpshell.write(shell_line)
			shell_line = " ".join(["cat", failedfqtmplist[fi], ">", newallfailedfqfile + "\n", "mv", failedfqtmplist[fi], tmpdatadir + "\n"])
			wtmpshell.write(shell_line)
		else:
			shell_line = " ".join(["cat", final_fq_list[fi], ">>", newallBXfqfile + "\n", "mv", final_fq_list[fi], tmpdatadir + "\n"])
			wtmpshell.write(shell_line)
			shell_line = " ".join(["cat", failedfqtmplist[fi], ">>", newallfailedfqfile + "\n", "mv", failedfqtmplist[fi], tmpdatadir + "\n"])
			wtmpshell.write(shell_line)
	othertmpfile = " ".join(otherfiletmplist)
	shell_line = " ".join(["mv", outputdir + "/split*txt*", othertmpfile, tmpdatadir + "\n"])
	wtmpshell.write(shell_line)
	wtmpshell.close()
	subprocess.call(["sh", tmpshell])
	subprocess.call(["mv", tmpshell, tmpdatadir])

	print("[ %s ] processing fastq QC ... \n" % time.asctime())
	statdir = outputdir + "/fastqc"
	check_info(statdir, "dir")
	Fastqc_path = G.Fastqc()
	O.FastQC(newallBXfqfile, Fastqc_path, statdir)
	print("[ %s ] fastq QC has finished, results have been writeen to %s \n\n" % (time.asctime(), statdir))

#	subprocess.call(["mv", newallfq1file, newallfq2file, SplitFqBarcodeFileList, tmpdatadir])
#### encode
