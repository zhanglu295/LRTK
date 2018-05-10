import os, sys, re, gzip
import getopt
import time
import subprocess
import multiprocessing
import pysam, random, string

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
	
	def Samtools(self):
		abs_path = self.get_path('samtools')
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

def check_info(result, attribute):
	if attribute == "file":
		if os.path.isfile(result) == False:
			sys.stderr.write("[ %s ] %s does not exist!\n" % (time.asctime(), result))
			sys.exit(-1)
	elif attribute == "dir":
		if os.path.isdir(result) == False:
			os.makedirs(result)

class modify_barcode:
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
	
	def generate_barcode(self, split_fq_sam_file_List, outdir, samtools_path, split_modify_fq_file_List, maxsize, barcode_sequence, sampleid):
		fqsamlist = list()
		rsplit_fq_sam_file_List = open(split_fq_sam_file_List, 'r')
		for split_fq_sam_file_List_info in rsplit_fq_sam_file_List:
			fqsamlist.append(split_fq_sam_file_List_info.strip())
		rsplit_fq_sam_file_List.close()
		mergedFqBam = outdir + "/" + barcode_sequence + ".merged.bam"
		shelldir = outdir + "/shell"
		check_info(shelldir, "dir")
		shell = shelldir + "/" + barcode_sequence + ".merge_fqbam.sh"
		wshell = open(shell, 'w')
		if len(fqsamlist) > 1:
			allfqsam = " ".join(fqsamlist)
			shell_line = " ".join([samtools_path, "merge -f", mergedFqBam, allfqsam]) + "\n"
			wshell.write(shell_line)
			shell_line = "rm " + allfqsam + "\n"
			wshell.write(shell_line)
		else:
			shell_line = " ".join(["cp", fqsamlist[0], mergedFqBam]) + "\n"
			wshell.write(shell_line)
			shell_line = "rm " + fqsamlist[0] + "\n"
			wshell.write(shell_line)
		sv = find_samtools_version(samtools_path, shelldir)
		sortedFqBamprefix = mergedFqBam.replace("merged.bam", "merged.sorted")
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
		subprocess.call(["rm", mergedFqBam])
		new_fq_prefix = sortedFqBam.replace("merged.sorted.bam", "")

		wsplit_modify_fq_file_List = open(split_modify_fq_file_List, 'w')
		SplitSize = 0
		rsortedFqBam = pysam.AlignmentFile(sortedFqBam, 'rb')
		readid = "N"
		barcodeid = 0
		barcode_marker = "N"
		new_barcode_sequence = "N"
		start = 11

		s = 0
		split_modify_fq_file = outdir + "/" + barcode_sequence + "." + str(s) + ".BX.modified.fq.gz"
		split_modify_fq_file_List_info = split_modify_fq_file + "\n"
#		wsplit_modify_fq_file_List.write(sampleid + "\tMerged_libraries\t" + split_modify_fq_file + "\n")
		wsplit_modify_fq_file_List.write(split_modify_fq_file + "\n")
		wsplit_modify_fq_file = gzip.open(split_modify_fq_file, 'wb')
		outFQfile = outdir + "/" + (os.path.basename(merged_fq)).replace("fq.gz", "") + barcode_sequence + ".fq.gz"
		woutFQfile = gzip.open(outFQfile, 'wb')
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
				split_modify_fq_file = outdir + "/" + barcode_sequence + "." + str(s) + ".BX.modified.fq.gz"
				wsplit_modify_fq_file_List.write(split_modify_fq_file + "\n")
				wsplit_modify_fq_file = gzip.open(split_modify_fq_file, 'wb')
			wsplit_modify_fq_file.write(complete_read_info.encode())
			woutFQfile.write(complete_read_info.encode())
		rsortedFqBam.close()
		wsplit_modify_fq_file_List.close()
		wsplit_modify_fq_file.close()
		woutFQfile.close()
		
		subprocess.call(["rm", sortedFqBam, sortedFqBam + ".bai"])

		return(outFQfile)

	def FQ2bam(self, rawfq, outbam, samtools_path):
		print("[ %s ] processing %s\n" % (time.asctime(), rawfq))
		if rawfq.endswith('gz'):
			o_rawfq = gzip.open(rawfq, "rb")
		else:
			o_rawfq = open(rawfq, "r")

		fq2bam = outbam + ".sam"
		wfq2bam = open(fq2bam, 'w')
		fq2baminfo = "@HD\tVN:1.3\n@SQ\tSN:chr1\tLN:248956422\n"
		wfq2bam.write(fq2baminfo)

		while True:
			rawfq_id1 = o_rawfq.readline()
			rawfq_seq1 = o_rawfq.readline()
			rawfq_str1 = o_rawfq.readline()
			rawfq_qua1 = o_rawfq.readline()

			rawfq_id2 = o_rawfq.readline()
			rawfq_seq2 = o_rawfq.readline()
			rawfq_str2 = o_rawfq.readline()
			rawfq_qua2 = o_rawfq.readline()

			if len(rawfq_id1) == 0:
				break
			else:
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

				realbarcode = (re.split(":", (re.split("\s+", rawfq_id1))[-1]))[-1]
				readid = realbarcode + "_" + ((re.split("\s+", rawfq_id1))[0]).replace('@', '')
				readid = (re.split("/", readid))[0]
				readlength1 = str(len(rawfq_seq1)) + "M"
				readlength2 = str(len(rawfq_seq2)) + "M"

				fq2baminfo = "\t".join([readid, "0", "*", "0", "0", readlength1, "*", "0", "0", rawfq_seq1, rawfq_qua1]) + "\n"
				wfq2bam.write(fq2baminfo)
				fq2baminfo = "\t".join([readid, "0", "*", "0", "0", readlength2, "*", "0", "0", rawfq_seq2, rawfq_qua2]) + "\n"
				wfq2bam.write(fq2baminfo)

		wfq2bam.close()
		o_rawfq.close()

		shell_dir = os.path.dirname(outbam) + "/shell"
		check_info(shell_dir, "dir")
		shell = os.path.join(shell_dir, os.path.basename(outbam) + ".sortbam.sh")
		wshell = open(shell, 'w')
		shell_line = " ".join([samtools_path, "view -S -h -b", fq2bam, ">", outbam]) + "\n"
		wshell.write(shell_line)
		shell_line = "rm " + fq2bam + "\n"
		wshell.write(shell_line)
		wshell.close()
		subprocess.call(["sh", shell])

def usage():
	mergefq_usage = \
	'''
	merge fastq files of multiple libraries
	Version: 1.0.0
	Dependents: Python (>=3.0)
	Last Updated Date: 2017-12-26

	Usage: python mergeFQ.py <options>

	Basic options:
		-i --input, list of fastq files of multiple libraries
		-o --output, output file
		-c --config, configuration file
		-n --process_num, the max amout of sub-processing [default: 4]
		-z --splitsize, the amount of reads of splited fastq, reads num = split_size / 8 [default: 15000000 lines, 1875000 read-pairs, compressed file size: ~300M]

	Warnings: only gziped or uncompressed files were supported, and output file would be gziped automatically.

	'''
	print(mergefq_usage)

if __name__ == '__main__':
	if len(sys.argv) < 5:
		usage()
		sys.exit(-1)

	clean_fq_list = None
	merged_fq = None
	process_num = 4
	ConfigFile = None
	splitsize = 15000000

	opts, args = getopt.gnu_getopt(sys.argv[1:], 'i:o:n:c:z:', ['input', 'output', 'process_num', 'config', 'splitsize'])
	for o, a in opts:
		if o == '-i' or o == '--input':
			clean_fq_list = str(a)
		if o == '-o' or o == '--output':
			merged_fq = str(a)
		if o == '-n' or o == '--process_num':
			process_num = int(a)
		if o == '-c' or o == '--config':
			ConfigFile = str(a)
		if o == '-z' or o == '--splitsize':
			splitsize = a
	
	if os.path.isfile(clean_fq_list):
		pass
	else:
		print("[ %s ] ERROR: %s does not exist!\n" % (time.asctime(), clean_fq_list))

	if merged_fq.endswith(".gz") == False:
		merged_fq = merged_fq + ".gz"

	G = baseinfo()
	G.get_config(ConfigFile)
	samtools = G.Samtools()

	M = modify_barcode()

	im = 0
	pool = multiprocessing.Pool(processes = process_num)

	SampleID = None
	i = 0
	rclean_fq_list = open(clean_fq_list, 'r')
	outbamdir = os.path.dirname(merged_fq) + "/Merged_libraries"
	check_info(outbamdir, "dir")

	sign_file = outbamdir + "/shell/cat_FQ.sh.sign"
	if os.path.exist(sign_file) and os.path.getsize(sign_file):
		print("[ %s ] all fastq files have been merged into %s!\n\n" % (time.asctime(), merged_fq))
		sys.exit(1)

	Split_FQ2Bam_list = list()
	barcode_seq_dict = dict()
	check_info(outbamdir, "dir")
	for clean_fq_list_info in rclean_fq_list:
		clean_fq_list_info_List = re.split("\t", clean_fq_list_info.strip())
		SampleID = clean_fq_list_info_List[0]
		LibraryID = clean_fq_list_info_List[1]
		clean_fq_List = clean_fq_list_info_List[2]
		barcode_seq = clean_fq_list_info_List[3]

		if barcode_seq not in barcode_seq_dict:
			if len(barcode_seq_dict) > 0:
				wSplit_FQ2Bam_file.close()
			Split_FQ2Bam_file = outbamdir + "/" + barcode_seq + ".fq2bam.txt"
			barcode_seq_dict[barcode_seq] = Split_FQ2Bam_file
			Split_FQ2Bam_list.append(Split_FQ2Bam_file)
			wSplit_FQ2Bam_file = open(Split_FQ2Bam_file, 'w')

		if os.path.exists(clean_fq_List) == False or os.path.getsize(clean_fq_List) == 0:
			print("[ %s ] ERROR: %s does not exist!\n" % (time.asctime(), clean_fq))
			sys.exit(-1)

		rclean_fq_List = open(clean_fq_List, 'r')
		for clean_fq in rclean_fq_List:
			clean_fq = clean_fq.strip()
			outbam = outbamdir + "/" + os.path.basename(clean_fq)
			if outbam.endswith(".gz"):
				outbam = outbam.replace(".gz", "")
			outbam = outbam + ".fq.bam"
			wSplit_FQ2Bam_file.write(outbam + "\n")

#			M.FQ2bam(clean_fq, outbam, samtools)
			im += 1
			if im <= process_num:
				pool.apply_async(M.FQ2bam, (clean_fq, outbam, samtools))
			if im == process_num:
				pool.close()
				pool.join()

				im = 0
				pool = multiprocessing.Pool(processes = process_num)

	if im > 0:
		pool.close()
		pool.join()
	rclean_fq_list.close()
	wSplit_FQ2Bam_file.close()

	im = 0
	pool = multiprocessing.Pool(processes = process_num)

	outFQdir = outbamdir
	Modified_merged_FQ_Info_list = list()
	Split_Modified_FQ_List = outFQdir + "/Split_Modified_FQ.txt"
	Split_Modified_FQ_List_list = list()
	for barcode_seq in barcode_seq_dict.keys():
		Split_FQ2Bam_file_List = barcode_seq_dict[barcode_seq]

		Split_Modified_FQ_List_each = Split_Modified_FQ_List + "." + barcode_seq
		Split_Modified_FQ_List_list.append(Split_Modified_FQ_List_each)

		Modified_merged_FQ_Info_each = outFQdir + "/" + (os.path.basename(merged_fq)).replace("fq.gz", "") + barcode_seq + ".fq.gz"
		Modified_merged_FQ_Info_list.append(Modified_merged_FQ_Info_each)

#		M.generate_barcode(Split_FQ2Bam_file_List, outFQdir, samtools, Split_Modified_FQ_List_each, int(splitsize), barcode_seq, SampleID)
		im += 1
		if im <= process_num:
			pool.apply_async(M.generate_barcode, (Split_FQ2Bam_file_List, outFQdir, samtools, Split_Modified_FQ_List_each, int(splitsize), barcode_seq, SampleID))
		if im == process_num:
			pool.close()
			pool.join()

			im = 0
			pool = multiprocessing.Pool(processes = process_num)
	if im > 0:
		pool.close()
		pool.join()

	shell_dir = outbamdir + "/shell"
	check_info(shell_dir, "dir")
	shell = shell_dir + "/cat_FQ.sh"
	wshell = open(shell, 'w')
	shell_line = "#!/bin/bash\nset -e\n"
	wshell.write(shell_line)
	for i in range(0,len(Split_Modified_FQ_List_list)):
		if i == 0:
			shell_line = " ".join(["cat", Split_Modified_FQ_List_list[i], ">", Split_Modified_FQ_List + "\nrm", Split_Modified_FQ_List_list[i] + "\n"])
			wshell.write(shell_line)
			shell_line = " ".join(["cat", Modified_merged_FQ_Info_list[i], ">", merged_fq + "\nrm", Modified_merged_FQ_Info_list[i] + "\n"])
			wshell.write(shell_line)
		else:
			shell_line = " ".join(["cat", Split_Modified_FQ_List_list[i], ">>", Split_Modified_FQ_List + "\nrm", Split_Modified_FQ_List_list[i] + "\n"])
			wshell.write(shell_line)
			shell_line = " ".join(["cat", Modified_merged_FQ_Info_list[i], ">>", merged_fq + "\nrm", Modified_merged_FQ_Info_list[i] + "\n"])
			wshell.write(shell_line)
	shell_line = "echo Done > " + shell + ".sign"
	wshell.write(shell_line)
	wshell.close()
	subprocess.call(["sh", shell])

	print("[ %s ] all fastq files have been merged into %s!\n\n" % (time.asctime(), merged_fq))
