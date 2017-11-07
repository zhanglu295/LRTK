import os, sys, gzip
import getopt
import time
import subprocess
import re
import multiprocessing
import random, string

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
	def revise(self, rawfq, outdir, maxsize):
		if rawfq.endswith('gz'):
			o_rawfq = gzip.open(rawfq, 'rb')
		else:
			o_rawfq = open(rawfq, 'r')

		(rawfq_path, rawfq_name) = os.path.split(rawfq)

		b = str(rawfq_name)
		if b[(len(b)-9):] == ".fastq.gz":
			b = b[0:(len(b)-9)]
		elif b[(len(b)-6):] == ".fq.gz" or b[(len(b)-6):] == ".fastq":
			b = b[0:(len(b)-6)]
		elif b[(len(b)-3):] == ".fq":
			b = b[0:(len(b)-3)]

		OutFileList = os.path.join(outdir, "split_fq_barcode.txt")
		wOutFileList = open(OutFileList, 'w')
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

			rawfq_id2 = o_rawfq.readline()
			rawfq_seq2 = o_rawfq.readline()
			rawfq_str2 = o_rawfq.readline()
			rawfq_qua2 = o_rawfq.readline()

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
		return(fq_prefix_name, OutFileList)

class modify_barcode:
	def replace_barcode(self, rawfq, barcode, barcode_sam, bwa_path, bwa_parameter, barcode_ref_prefix, new_fq_prefix, add):
		outdir = os.path.dirname(barcode_sam)
		shelldir = outdir + "/shell"
		check_info(shelldir)
		tmpprefix = os.path.basename(new_fq_prefix)
		shell = shelldir + "/" + tmpprefix + ".barcode.align.sh"
#		shell = os.path.join(shelldir, "barcode.align.sh")
		oshell = open(shell, 'w')
		sai = barcode_sam + '.sai'
		shell_line = " ".join([bwa_path, 'aln', barcode_ref_prefix, barcode, '-f', sai, bwa_parameter, "\n"])
		oshell.write(shell_line)
		shell_line = " ".join([bwa_path, 'samse', barcode_ref_prefix, sai, barcode, "| grep -v '^@' > " + barcode_sam + "\nrm", sai, "\n"])
		oshell.write(shell_line)
		oshell.close()
		print("[ %s ] processing %s\n" % (time.asctime(), shell))
		subprocess.call(["sh", shell])
		print("[ %s ] %s finished!\n" % (time.asctime(), shell))

		print("\n[ %s ] modify barcode info of %s\n" % (time.asctime(), rawfq))
		samfile = open(barcode_sam, "r")
		newfqfile1 = new_fq_prefix + "_1.fq.gz"
		newfqfile2 = new_fq_prefix + "_2.fq.gz"

		original_barcode_count = 0
		mapped_barcode_count = 0
		corrected_barcode_count = 0

		BXfq = new_fq_prefix + ".BX.fq.gz"
		if add == 1:
			wBXfq = gzip.open(BXfq, 'wb')

		Failedfq = new_fq_prefix + ".failed.fq.gz"
		wFailedfq = gzip.open(Failedfq, 'wb')

		if rawfq.endswith('gz'):
			o_rawfq = gzip.open(rawfq, "rb")
		else:
			o_rawfq = open(rawfq, "r")

		newfq1_out = gzip.open(newfqfile1, 'wb')
		newfq2_out = gzip.open(newfqfile2, 'wb')

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

				id1 = re.split(" ", rawfq_id1)
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
					
						newfq_seq1 = realbarcode + rawfq_seq1[23:]
						newfq_qua1 = realbarcodequa + rawfq_qua1[23:]
						newfq1 = rawfq_id1 + "\n" + newfq_seq1 + "\n" + rawfq_str1 + "\n" + newfq_qua1 + "\n"
						newfq_seq2 = str(realbarcode) + str(rawfq_seq2)
						newfq_qua2 = str(realbarcodequa) + str(rawfq_qua2)
						newfq2 = "\n".join([rawfq_id2, newfq_seq2, rawfq_str2, newfq_qua2, ""])

						newfq1_out.write(newfq1.encode())
						newfq2_out.write(newfq2.encode())

						if add == 1:
							newfq_id1 = rawfq_id1 + " BX:" + realbarcode
							newfq_seq1 = rawfq_seq1[23:]
							newfq_qua1 = rawfq_qua1[23:]
							newfq_id2 = rawfq_id2 + " BX:" + realbarcode
							newfq_seq2 = str(rawfq_seq2)
							newfq_qua2 = str(rawfq_qua2)

							BXout = "\n".join([newfq_id1, newfq_seq1, rawfq_str1, newfq_qua1, newfq_id2, newfq_seq2, rawfq_str2, newfq_qua2, ""])
							wBXfq.write(BXout.encode())
					else:
						Failedout = "\n".join([rawfq_id1, rawfq_seq1, rawfq_str1, rawfq_qua1, rawfq_id2, rawfq_seq2, rawfq_str2, rawfq_qua2, ""])
						wFailedfq.write(Failedout.encode())
		wFailedfq.close()
		newfq1_out.close()
		newfq2_out.close()
		if add == 1:
			wBXfq.close()
		print("[ %s ] barcode info of %s has been modified\n" % (time.asctime(), rawfq))
		return(raw_fq1_info_list, raw_fq2_info_list, clean_fq1_info_list, clean_fq2_info_list, count_info1_list, count_info2_list, original_barcode_count, mapped_barcode_count, corrected_barcode_count)

class OUTERSOFT:
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

	def fill_barcode(self, original_sam, samtools_path):
		print("[ %s ] modify barcode info of %s \n" % (time.asctime(), original_sam))
		original_name = os.path.basename(original_sam)
		new_sam = original_sam.replace('sam', 'new.sam')
		new_bam = new_sam.replace('sam', 'bam')
		sorted_bam = new_sam.replace('new.sam', 'sorted')

		roriginal_sam = open(original_sam, 'r')
		o_new_sam = open(new_sam, 'w')
		convert_strand = {"a":"T", "c":"G", "g":"C", "t":"A"}
		
		while True:
			saminfo = roriginal_sam.readline()
			if len(saminfo) == 0:
				break
			else:
				if saminfo.startswith('@'):
					o_new_sam.write(saminfo)
				else:
					saminfo2 = roriginal_sam.readline()
					saminfo = saminfo.strip()
					saminfo2 = saminfo2.strip()
					saminfolist = re.split('\t', saminfo)
					saminfolist2 = re.split('\t', saminfo2)
					N_info = 0
					for s in range(len(saminfolist)):
						if 'BC:Z' in saminfolist[s]:
							if re.search('N', saminfolist[s], re.IGNORECASE):
								N_info = 1
							if saminfolist[s] != saminfolist2[s]:
								errinfo = saminfolist[0] + "\t" + saminfolist[s] + "\t" + saminfolist2[s] + "\n"
								sys.stderr.write("Warning! The barcode of paired reads are not the same: %s\n" % errinfo)
							BClist = list(saminfolist[s])
							BClist[1] = "X"
							for b in range(5,21):
								if BClist[b].islower():
									BClist[b] = convert_strand[BClist[b]]
							BClist[21] = "-1"
							saminfolist[s] = "".join(BClist[0:22])
							saminfolist2[s] = saminfolist[s]

					if N_info == 0:
						saminfo = "\t".join(saminfolist) + '\n'
						o_new_sam.write(saminfo)
						saminfo2 = "\t".join(saminfolist2) + '\n'
						o_new_sam.write(saminfo2)

		roriginal_sam.close()
		o_new_sam.close()

		shelldir = os.path.dirname(sorted_bam) + "/shell"
		check_info(shelldir)
		shell = shelldir + "/" + os.path.basename(sorted_bam) + ".sh"
		oshell = open(shell, 'w')
		sv = self.find_samtools_version(samtools_path, shelldir)
		if sv == 0:
			shell_line = " ".join([samtools_path, "view -h -S -b", new_sam, ">", new_bam, "\n"])
			oshell.write(shell_line)
			shell_line = " ".join([samtools_path, "sort -m 1G", new_bam, sorted_bam, "\n"])
			oshell.write(shell_line)
		else:
			shell_line = " ".join([samtools_path, "view -h -S -b -o", new_bam, new_sam + "\n"])
			oshell.write(shell_line)
			shell_line = " ".join([samtools_path, "sort -m 2G -o", sorted_bam + ".sorted.bam -T", sorted_bam, new_bam, "\n"])
			oshell.write(shell_line)
			shell_line = " ".join(["mv", sorted_bam + ".sorted.bam", sorted_bam + ".bam\n"])
			oshell.write(shell_line)
		shell_line = " ".join(["rm", original_sam, new_sam, new_bam, "\n"])
		oshell.write(shell_line)
		oshell.close()
		subprocess.call(["sh", shell])
		print("[ %s ] modify barcode info, finished! %s \n" % (time.asctime(), sorted_bam + ".bam"))

	def bwa_fq(self, fq1, fq2, outsam, bwa_path, bwa_aln_parameter, bwa_sam_parameter, RGinfo, samtools_path):
		sai1 = fq1 + '1.sai'
		sai2 = fq2 + '2.sai'

		print("[ %s ] fq alignment starts: %s %s \n" % (time.asctime(), fq1, fq2))
		fqprefix = (os.path.basename(fq1)).replace("_1.fq.gz", "")
		outdir = os.path.dirname(fq1)
		shelldir = os.path.join(outdir, "shell")
		check_info(shelldir)
		shell = os.path.join(shelldir, fqprefix + ".fq.aln.sh")
		oshell = open(shell, 'w')
		shell_line = " ".join([bwa_path, 'aln', ref_path, fq1, bwa_aln_parameter, '-f', sai1, '&\n'])
		oshell.write(shell_line)
		shell_line = " ".join([bwa_path, 'aln', ref_path, fq2, bwa_aln_parameter, '-f', sai2, '&\nwait\n'])
		oshell.write(shell_line)
		if RGinfo != 0:
			RGinfo = '"' + RGinfo + '"'
			shell_line = " ".join([bwa_path, 'sampe', bwa_sam_parameter, '-r', RGinfo, ref_path, sai1, sai2, fq1, fq2, ">", outsam + "\nrm", sai1, sai2, "\n"])
		else:
			shell_line = " ".join([bwa_path, 'sampe', bwa_sam_parameter, ref_path, sai1, sai2, fq1, fq2, ">", outsam + "\nrm", sai1, sai2, "\n"])
		oshell.write(shell_line)
		oshell.close()

		subprocess.call(["sh", shell])
		print("[ %s ] fq alignment finish: %s \n" % (time.asctime(), outsam))

	def FastQC(self, fq_prefix, fastqc_path, statdir):
		fq1 = fq_prefix + "_1.fq.gz"
		fq2 = fq_prefix + "_2.fq.gz"

		tmpshell = os.path.join(statdir, "fastqc.sh")
		wtmpshell = open(tmpshell, 'w')
		shell_line = " ".join([fastqc_path, "-f fastq -o", statdir, fq1, "&\n"])
		wtmpshell.write(shell_line)
		shell_line = " ".join([fastqc_path, "-f fastq -o", statdir, fq2, "&\n"])
		wtmpshell.write(shell_line)
		shell_line = "wait\necho Done\n"
		wtmpshell.write(shell_line)
		wtmpshell.close()

		subprocess.call(["sh", tmpshell])

def check_info(result):
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
	Last Updated Date: 2017-06-01
	Contact: meijp@foxmail.com

	Usage: python fastq_purification.py <options>

	Options:
		-i --input, the input path of compressed or uncompressed fastq files.
		-o --outputdir, the path of output directory
		-c --config, the path of configuration file [default: outdir/config/Basic.config]
		-z --splitsize, the amount of reads of splited fastq, reads num = split_size / 4 [default: 5000000 lines, 1250000 reads]
		-r --RG, RG info for alignment
		-n --process_num, the max amout of sub-processing [default: 4]
		-N --noBX, add BX info to the end of read id or not [default yes]
		-h --help, help info

	'''
	print(purification_usage)

if __name__ == '__main__':
	if len(sys.argv) < 5:
		usage()
		sys.exit(-1)

	RGInfo = 0
	inputfq = None
	outputdir = None
	ConfigFile = None
	addBX = 1
	split_size = 5000000
	process_num = 4
	opts, args = getopt.gnu_getopt(sys.argv[1:], 'i:o:c:z:n:r:hN', ['input', 'outputdir', 'config', 'splitsize', 'process_num', 'RG', 'help', 'noBX'])
	for o, a in opts:
		if o == '-i' or o == '--input':
			inputfq = a
		if o == '-o' or o == '--outputdir':
			outputdir = a
		if o == '-c' or o == '--config':
			ConfigFile = a
		if o == '-z' or o == '--splitsize':
			split_size = int(a)
		if o == '-n' or o == '--process_num':
			process_num = int(a)
		if o == '-r' or o == '--RG':
			RGInfo = a
		if o == '-h' or o == '--help':
			usage()
			sys.exit(-1)
		if o == '-N' or o == '--noBX':
			addBX = 0

	if ConfigFile == None:
		script_abs_path = os.path.abspath(sys.argv[0])
		create_config_py = os.path.join(os.path.dirname(script_abs_path), "create_config.py")
		config_dir = os.path.join(outputdir, "config")
		check_info(config_dir)
		tmpshell = os.path.join(config_dir, "cc.sh")
		wtmpshell = open(tmpshell, 'w')
		shell_line = " ".join(["python", create_config_py, "Basic -o", config_dir, "\n"])
		wtmpshell.write(shell_line)
		wtmpshell.close()
		subprocess.call(["sh", tmpshell])
		subprocess.call(["rm", tmpshell])
		ConfigFile = os.path.join(config_dir, "Basic.config")

	G = baseinfo()
	G.get_config(ConfigFile)

	R = extract_barcode()
	SplitFqBarcodeFileList = outputdir + '/split_fq_barcode.txt'
	print("[ %s ] extract original barcode info from %s\n" % (time.asctime(), inputfq))
	(prefix_name, SplitFqBarcodeFileList) = R.revise(inputfq, outputdir, split_size)
	print("[ %s ] original barcode had been extracted, and raw fq file had been splited, output files has been listed in : %s\n" % (time.asctime(), SplitFqBarcodeFileList))

	barcoderef = G.Barcode()
	bracode_aln_parameter = G.Barcode_aln_par()
	bwa_path = G.Bwa()
	rSplitFqBarcodeFileList = open(SplitFqBarcodeFileList, 'r')
	M = modify_barcode()
	fq_statistics_result = list()
	im = 0
	pool = multiprocessing.Pool(processes = process_num)
	newfq1filetmplist = list()
	newfq2filetmplist = list()
	BXfqfiletmplist = list()
	failedfqtmplist = list()
	otherfiletmplist = list()

	SplitCleanFqFileList = outputdir + '/split_clean_fq.txt'
	wSplitCleanFqFileList = open(SplitCleanFqFileList, 'w')
	for splitinfo in rSplitFqBarcodeFileList:
		(split_fq_file, split_barcode_file) = re.split("\t", splitinfo.strip())
		barcodesam = (str(split_barcode_file))[0:(len(split_barcode_file)-2)] + "sam"
		otherfiletmplist.append(barcodesam)
		otherfiletmplist.append(split_fq_file)
		otherfiletmplist.append(split_barcode_file)
		newprefix = (str(split_barcode_file))[0:(len(split_barcode_file)-11)]

		im += 1
		if im <= process_num:
			fq_statistics_result.append(pool.apply_async(M.replace_barcode, (split_fq_file, split_barcode_file, barcodesam, bwa_path, bracode_aln_parameter, barcoderef, newprefix, addBX,)))
		if im == process_num:
			pool.close()
			pool.join()

			im = 0
			pool = multiprocessing.Pool(processes = process_num)

		newfq1file = newprefix + "_1.fq.gz"
		newfq1filetmplist.append(newfq1file)
		newfq2file = newprefix + "_2.fq.gz"
		newfq2filetmplist.append(newfq2file)
		BXfqfile = newprefix + ".BX.fq.gz"
		BXfqfiletmplist.append(BXfqfile)
		failedfqfile = newprefix + ".failed.fq.gz"
		failedfqtmplist.append(failedfqfile)

		SplitCleanFqFile = "\t".join([newfq1file, newfq2file, BXfqfile]) + "\n"
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

	O = OUTERSOFT()
	ref_path = G.Ref()
	bwa_path = G.Bwa()
	bwa_aln_parameter = G.Fq_aln_par()
	bwa_sam_parameter = G.Fq_sam_par()
	samtools = G.Samtools()

	SplitSamFileList = outputdir + '/split_sam.txt'
	wSplitSamFileList = open(SplitSamFileList, 'w')
	im = 0
	pool = multiprocessing.Pool(processes = process_num)
	rSplitCleanFqFileList = open(SplitCleanFqFileList, 'r')
	for SplitCleanFqFileInfo in rSplitCleanFqFileList:
		SplitCleanFqFileInfolist = re.split("\t", SplitCleanFqFileInfo.strip())
		fqsam = (str(SplitCleanFqFileInfolist[0]))[0:(len(SplitCleanFqFileInfolist[0])-8)] + ".sam"
		SplitSamFileinfo = fqsam + "\n"

		wSplitSamFileList.write(SplitSamFileinfo)

		im += 1
		if im <= process_num:
			pool.apply_async(O.bwa_fq, (SplitCleanFqFileInfolist[0], SplitCleanFqFileInfolist[1], fqsam, bwa_path, bwa_aln_parameter, bwa_sam_parameter, RGInfo, samtools, ))
		if im == process_num:
			pool.close()
			pool.join()

			im = 0
			pool = multiprocessing.Pool(processes = process_num)
	
	if im > 0:
		pool.close()
		pool.join()
	wSplitSamFileList.close()
	subprocess.call(["sh", SplitCleanFqFileList])
	
	rSplitSamFileList = open(SplitSamFileList, 'r')
	im = 0
	pool = multiprocessing.Pool(processes = process_num)
	SplitBamFilelist = list()
	for SplitSamFile in rSplitSamFileList:
		SplitSamFile = SplitSamFile.strip()
		SplitBamFile = SplitSamFile.replace(".sam", ".sorted.bam")
		SplitBamFilelist.append(SplitBamFile)
		im += 1
		if im <= process_num:
			pool.apply_async(O.fill_barcode, (SplitSamFile, samtools, ))
		if im == process_num:
			pool.close()
			pool.join()

			im = 0
			pool = multiprocessing.Pool(processes = process_num)
	if im > 0:
		pool.close()
		pool.join()
#		O.fill_barcode(SplitSamFile, samtools)
	rSplitSamFileList.close()
	subprocess.call(["rm", SplitSamFileList])

	allSplitBamFile = " ".join(SplitBamFilelist)
	tmpmergeshell = SplitBamFilelist[0] + ".merge.sh"
	wtmpmergeshell = open(tmpmergeshell, 'w')
	MergeBamFile = prefix_name + ".sorted.bam"
	if os.path.isfile(MergeBamFile):
		subprocess.call(["rm", MergeBamFile])
	if len(SplitBamFilelist) > 1:
		shell_line = " ".join([samtools, "merge", MergeBamFile, allSplitBamFile]) + "\n"
	else:
		shell_line = " ".join(["cp", SplitBamFilelist[0], MergeBamFile]) + "\n"
	wtmpmergeshell.write(shell_line)
	shell_line = " ".join(["rm", allSplitBamFile]) + "\n"
	wtmpmergeshell.write(shell_line)
	wtmpmergeshell.close()
	subprocess.call(["sh", tmpmergeshell])
	subprocess.call(["rm", tmpmergeshell])

	tmpdatadir = outputdir + "/tmp"
	check_info(tmpdatadir)

	newallfq1file = prefix_name + "_1.fq.gz"
	newallfq2file = prefix_name + "_2.fq.gz"
	newallBXfqfile = prefix_name + ".BX.fq.gz"
	newallfailedfqfile = prefix_name + ".failed.fq.gz"
	tmpshell = newallfq1file + ".cat.sh"
	wtmpshell = open(tmpshell, 'w')
	for fi in range(0,len(newfq1filetmplist)):
		if fi == 0:
			shell_line = " ".join(["cat", newfq1filetmplist[fi], ">", newallfq1file + "\n", "mv", newfq1filetmplist[fi], tmpdatadir + "\n"])
			wtmpshell.write(shell_line)
			shell_line = " ".join(["cat", newfq2filetmplist[fi], ">", newallfq2file + "\n", "mv", newfq2filetmplist[fi], tmpdatadir + "\n"])
			wtmpshell.write(shell_line)
			shell_line = " ".join(["cat", BXfqfiletmplist[fi], ">", newallBXfqfile + "\n", "mv", BXfqfiletmplist[fi], tmpdatadir + "\n"])
			wtmpshell.write(shell_line)
			shell_line = " ".join(["cat", failedfqtmplist[fi], ">", newallfailedfqfile + "\n", "mv", failedfqtmplist[fi], tmpdatadir + "\n"])
			wtmpshell.write(shell_line)
		else:
			shell_line = " ".join(["cat", newfq1filetmplist[fi], ">>", newallfq1file + "\n", "mv", newfq1filetmplist[fi], tmpdatadir + "\n"])
			wtmpshell.write(shell_line)
			shell_line = " ".join(["cat", newfq2filetmplist[fi], ">>", newallfq2file + "\n", "mv", newfq2filetmplist[fi], tmpdatadir + "\n"])
			wtmpshell.write(shell_line)
			shell_line = " ".join(["cat", BXfqfiletmplist[fi], ">>", newallBXfqfile + "\n", "mv", BXfqfiletmplist[fi], tmpdatadir + "\n"])
			wtmpshell.write(shell_line)
			shell_line = " ".join(["cat", failedfqtmplist[fi], ">>", newallfailedfqfile + "\n", "mv", failedfqtmplist[fi], tmpdatadir + "\n"])
			wtmpshell.write(shell_line)
	othertmpfile = " ".join(otherfiletmplist)
	shell_line = " ".join(["mv", othertmpfile, tmpdatadir + "\n"])
	wtmpshell.write(shell_line)
	wtmpshell.close()
	subprocess.call(["sh", tmpshell])
	subprocess.call(["mv", tmpshell, tmpdatadir])

	print("[ %s ] processing fastq QC ... \n" % time.asctime())
	statdir = outputdir + "/fastqc"
	check_info(statdir)
	Fastqc_path = G.Fastqc()
	O.FastQC(prefix_name, Fastqc_path, statdir)
	print("[ %s ] fastq QC has finished, results have been writeen to %s \n\n" % (time.asctime(), statdir))
	subprocess.call(["mv", newallfq1file, newallfq2file, SplitFqBarcodeFileList, tmpdatadir])
#### encode
