import os, sys, gzip
import getopt
import time
import subprocess
import re

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
#				print (name, path)
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

class extract_barcode:
	def revise(self, rawfq, outdir):
		if rawfq.endswith('gz'):
			o_rawfq = gzip.open(rawfq, 'rb')
		else:
			o_rawfq = open(rawfq, 'r')

		(rawfq_path, rawfq_name) = os.path.split(rawfq)

		b = str(rawfq_name)
		if b[(len(b)-8):] == "fastq.gz":
			b = b[0:(len(b)-8)] + "barcode.gz"
		elif b[(len(b)-5):] == "fq.gz":
			b = b[0:(len(b)-5)] + "barcode.gz"
		elif b[(len(b)-5):] == "fastq":
			b = b[0:(len(b)-5)] + "barcode.gz"
		elif b[(len(b)-2):] == "fq":
			b = b[0:(len(b)-2)] + "barcode.gz"

		outfile = os.path.join(outdir, b)
		barcode_out = gzip.open(outfile, 'wb')
	
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
					rawfq_id1 = rawfq_id1.decode()
					rawfq_seq1 = rawfq_seq1.decode()
					rawfq_str1 = rawfq_str1.decode()
					rawfq_qua1 = rawfq_qua1.decode()

				barcode_line = rawfq_id1 + rawfq_seq1[0:16] + '\n' + rawfq_str1 + rawfq_qua1[0:16] + "\n"
				barcode_out.write(barcode_line.encode())

		barcode_out.close()
		return outfile

class modify_barcode:
	def replace_barcode(self, rawfq, barcode_sam, prefix):
		samfile = open(barcode_sam, "r")
		newfqfile1 = prefix + "_1.fq.gz"
		newfqfile2 = prefix + "_2.fq.gz"

		if rawfq.endswith('gz'):
			o_rawfq = gzip.open(rawfq, "rb")
		else:
			o_rawfq = open(rawfq, "r")

		newfq1_out = gzip.open(newfqfile1, 'wb')
		newfq2_out = gzip.open(newfqfile2, 'wb')

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
				if rawfq.endswith('gz'):
					rawfq_id1 = rawfq_id1.decode().strip()
					rawfq_seq1 = rawfq_seq1.decode().strip()
					rawfq_str1 = rawfq_str1.decode().strip()
					rawfq_qua1 = rawfq_qua1.decode().strip()

					rawfq_id2 = rawfq_id2.decode().strip()
					rawfq_seq2 = rawfq_seq2.decode().strip()
					rawfq_str2 = rawfq_str2.decode().strip()
					rawfq_qua2 = rawfq_qua2.decode().strip()
				barcodeinfo = barcodeinfo.strip()
				barcodeinfo = re.split("\t", barcodeinfo)

				id1 = re.split(" ", rawfq_id1)
				if id1[0] != ('@' + barcodeinfo[0]):
					sys.stderr.write("ERROR - %s and %s do not match!\n" % (id1[0], barcodeinfo[0]))
				else:
					barcode_ori = list(barcodeinfo[9])
					realbarcodequa = str(barcodeinfo[10])
					barcode_err = 0
					if barcodeinfo[2] == "*" or re.search('N', barcodeinfo[9]):
						barcode_err = 1
						continue
					for md in barcodeinfo[11:]:
						if md.startswith("MD:Z"):
							md = re.split(":", md)
							alter = re.findall(r'\D', md[-1])
							loci = re.findall(r'\d+', md[-1])
							if len(alter) == 0:
								pass
							elif len(alter) > 1:
								barcode_err = 1
							else:
								for i in range(len(alter)):
									ix = int(loci[i]) + i
									barcode_ori[ix] = alter[i]

					if barcode_err == 0:
						realbarcode = "".join(barcode_ori)
					
						newfq_seq1 = realbarcode + rawfq_seq1[23:]
#					rsq = str(rawfq_seq1[24:])
#					k = "\t".join([str(len(realbarcode)), str(len(rsq)), str(len(newfq_seq1))])
#					print(k)
						newfq_qua1 = realbarcodequa + rawfq_qua1[23:]
						newfq1 = rawfq_id1 + "\n" + newfq_seq1 + "\n" + rawfq_str1 + "\n" + newfq_qua1 + "\n"
						newfq_seq2 = str(realbarcode) + str(rawfq_seq2)
						newfq_qua2 = str(realbarcodequa) + str(rawfq_qua2)
						newfq2 = "\n".join([rawfq_id2, newfq_seq2, rawfq_str2, newfq_qua2, ""])
#					k = "\t".join([str(len(newfq_seq2)), str(len(newfq_qua2))])
#					print(k)

						newfq1_out.write(newfq1.encode())
						newfq2_out.write(newfq2.encode())

		newfq1_out.close()
		newfq2_out.close()

		return(prefix)

	def fill_barcode(self, original_sam, samtools_path, outdir = './'):
		original_name = os.path.basename(original_sam)
		b = str(original_name)

		new_sam = outdir + "/" + b[0:(len(b)-3)] + "new.sam"
		new_bam = outdir + "/" + b[0:(len(b)-3)] + "bam"
		sorted_bam = outdir + "/" + b[0:(len(b)-3)] + "sorted"

		print(" ".join(["modified and sorted bam file:", sorted_bam + '.bam', "\n"]))
		
		o_original_sam = open(original_sam, 'r')
		o_new_sam = open(new_sam, 'w')

		convert_strand = {"a":"T", "c":"G", "g":"C", "t":"A"}

		while True:
			saminfo = o_original_sam.readline()
			if len(saminfo) == 0:
				break
			else:
				if saminfo.startswith('@'):
					o_new_sam.write(saminfo)
				else:
					saminfo2 = o_original_sam.readline()
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

		o_original_sam.close()
		o_new_sam.close()

		shell = sorted_bam + ".sh"
		oshell = open(shell, 'w')
		shell_line = " ".join([samtools_path, "view -h -S -b", new_sam, ">", new_bam, "\n"])
		oshell.write(shell_line)
		shell_line = " ".join([samtools_path, "sort -m 1G", new_bam, sorted_bam, "\n"])
		oshell.write(shell_line)
#		shell_line = " ".join([samtools_path, "view -h", sorted_bam + '.bam > ', sorted_bam + '.sam\n'])
#		oshell.write(shell_line)
		shell_line = " ".join(["#rm", original_sam, new_sam, new_bam, "\n"])
		oshell.write(shell_line)
		oshell.close()

		subprocess.call(["sh", shell])
		return(sorted_bam + '.bam', sorted_bam + '.sam')

class OUTERSOFT:
	def bwa_barcode(self, barcode, barcode_sam, bwa_path, bwa_parameter, ref_path):
		outdir = os.path.dirname(barcode_sam)
		shell = os.path.join(outdir, "barcode.align.sh")
		oshell = open(shell, 'w')
		sai = barcode_sam + '.sai'
		shell_line = " ".join([bwa_path, 'aln', ref_path, barcode, '-f', sai, bwa_parameter, "\n"])
		oshell.write(shell_line)
		shell_line = " ".join([bwa_path, 'samse', ref_path, sai, barcode, "| grep -v '^@' > " + barcode_sam + "\nrm", sai, "\n"])
		oshell.write(shell_line)
		oshell.close()
		print(shell)
		subprocess.call(["sh", shell])
		return(barcode_sam)

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
		-c --config, the path of configuration file [default: outdir/config/QC.config]
		-h --help, help info

	'''
	print(purification_usage)

if __name__ == '__main__':
	if len(sys.argv) < 5:
		usage()
		sys.exit(-1)

	inputfq = None
	outputdir = None
	ConfigFile = None
	opts, args = getopt.gnu_getopt(sys.argv[1:], 'i:o:c:h:', ['input', 'outputdir', 'config', 'help'])
	for o, a in opts:
		if o == '-i' or o == '--input':
			inputfq = a
		if o == '-o' or o == '--outputdir':
			outputdir = a
		if o == '-c' or o == '--config':
			ConfigFile = a
		if o == '-h' or o == '--help':
			usage()
			sys.exit(-1)

	if ConfigFile == None:
		script_abs_path = os.path.abspath(sys.argv[0])
		create_config_py = os.path.join(os.path.dirname(script_abs_path), "create_config.py")
		config_dir = os.path.join(outputdir, "config")
		if os.path.isdir(config_dir):
			pass
		else:
			os.mkdir(config_dir)
		tmpshell = os.path.join(config_dir, "cc.sh")
		wtmpshell = open(tmpshell, 'w')
		shell_line = " ".join(["python", create_config_py, "QC -o", config_dir, "\n"])
		wtmpshell.write(shell_line)
		wtmpshell.close()
		subprocess.call(["sh", tmpshell])
		subprocess.call(["rm", tmpshell])
		ConfigFile = os.path.join(config_dir, "QC.config")

	O = OUTERSOFT()

	G = baseinfo()
	G.get_config(ConfigFile)

	R = extract_barcode()
	sys.stderr.write("[ %s ] extract original barcode info from %s\n" % (time.asctime(), inputfq))
	extract_barcode = R.revise(inputfq, outputdir)
	sys.stderr.write("[ %s ] original barcode had been extracted, output file has been written to: %s\n\n" % (time.asctime(), extract_barcode))

	barcoderef = G.Barcode()
	bracode_aln_parameter = G.Barcode_aln_par()
	bwa_path = G.Bwa()
	barcodesam = (str(extract_barcode))[0:(len(extract_barcode)-2)] + "sam"
	sys.stderr.write("[ %s ] align to the barcode set\n" % time.asctime())
	barcodesam = O.bwa_barcode(extract_barcode, barcodesam, bwa_path, bracode_aln_parameter, barcoderef)
	sys.stderr.write("\n[ %s ] align to the barcode set, output sam file has been written to: %s\n\n" % (time.asctime(), barcodesam))

	M = modify_barcode()
	newprefix = (str(extract_barcode))[0:(len(extract_barcode)-11)]

	sys.stderr.write("[ %s ] replace barcode info\n" % time.asctime())
	newprefix = M.replace_barcode(inputfq, barcodesam, newprefix)
	sys.stderr.write("[ %s ] barcode has been replaced, and new fq files have been written to %s_[1/2].fq.gz\n\n" % (time.asctime(), newprefix))

	sys.stderr.write("[ %s ] processing fastq QC ... \n" % time.asctime())
	statdir = outputdir + "/fastqc"
	if os.path.isdir(statdir):
		pass
	else:
		os.mkdir(statdir)
	Fastqc_path = G.Fastqc()
	O.FastQC(newprefix, Fastqc_path, statdir)
	sys.stderr.write("[ %s ] fastq QC has finished, results have been writeen to %s \n\n" % (time.asctime(), statdir))


#### encode
