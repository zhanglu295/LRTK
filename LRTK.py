import os, re, sys, gzip
import getopt
import time
import subprocess

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
		ABS_PATH = self.configdist[name_variable]
		if ABS_PATH is None:
			print >> sys.stderr, "ERROR - %s - %s is empty!" % (time.asctime(), name_variable)
			sys.exit(-1)
		return ABS_PATH

	def Ref(self):
		abs_path = self.get_path('ref')
		return abs_path
	
	def Bwa(self):
		abs_path = self.get_path('bwa')
		return abs_path
	
	def Barcode(self):
		abs_path = self.get_path('barcode')
		return abs_path

	def Barcode_aln_par(self):
		abs_path = self.get_path('barcode_aln_parameter')
		return abs_path

class extract_barcode:
	def revise(self, rawfq, outdir):
		if rawfq.endswith('gz'):
			o_rawfq = gzip.open(rawfq, 'rb')
		else:
			o_rawfq = open(rawfq, 'r')

		(rawfq_path, rawfq_name) = os.path.split(rawfq)

		b = rawfq_name
		if(rawfq_name.find("fastq")):
	#		rawfq_name.replace("fastq", "barcode", 1)
			strinfo = re.compile('fastq')
			b = strinfo.sub('barcode', rawfq_name)
		else:
	#		rawfq_name.replace("fq", "barcode", 1)
			strinfo = re.compile('fq')
			b = strinfo.sub('barcode', rawfq_name)

		outfile = os.path.join(outdir, b)

		if outfile.endswith('gz'):
			barcode_out = gzip.open(outfile, 'wb')
		else:
			barcode_out = open(outfile, 'w')
	
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
		L1 = LRTK()
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

#					rawfq_id2 = rawfq_id2.decode().strip()
#					rawfq_seq2 = rawfq_seq2.decode().strip()
#					rawfq_str2 = rawfq_str2.decode().strip()
#					rawfq_qua2 = rawfq_qua2.decode().strip()

				barcodeinfo = re.split("\t", barcodeinfo.decode().strip())

				id1 = re.split(" ", rawfq_id1)
				if id1[0] != barcodeinfo[0]:
					print >> sys.stderr, "ERROR - %s and %s do not match!" % (id1[0], barcodeinfo[0])
				else:
					(barcode_ori, realbarcodequa) = (list(barcodeinfo[9]), barcodeinfo[10])
					for md in barcodeinfo[11:]:
						if md.startswith("MD:Z"):
							md = re.split(":", md)
							alter = re.findall(r'\D', md[-1])
							loci = re.findall(r'\d+', md[-1])
							if len(alter) == 0:
								pass
							else:
								for i in range(len(alter)):
									ix = int(loci[i]) + i
									barcode_ori[ix] = alter[i]

					realbarcode = "".join(barcode_ori)
					
					newfq_seq1 = realbarcode + rawfq_seq1[23:]
					newfq_qua1 = realbarcode = rawfq_qua1[23:]
					newfq1 = rawfq_id1 + "\n" + newfq_seq1 + "\n" + rawfq_str1 + "\n" + newfq_qua1 + "\n"
					newfq_seq2 = realbarcode + rawfq_seq2
					newfq_qua2 = realbarcode + rawfq_qua2
					newfq2 = rawfq_id2 + newfq_seq2 + rawfq_str2 + newfq_qua2

					newfq1_out.write(newfq1.encode())
					newfq2_out.write(newfq2.encode())

		newfq1_out.close()
		newfq2_out.close()

class CFCR:
	def calculate(sorted_by_barcode_sam, outfile):
		samfile = open(sorted_by_barcode_bam, 'r')

		allbarcode = dict()
		while True:
			saminfo1 = re.split("\t", samfile.readline())
			saminfo2 = re.split("\t", samfile.readline())

			if len(saminfo1) == 0:
				break
			elif saminfo1[8] == 0:
				pass
			else:
				
				for bc in saminfo1[11:]:
					if bc.startswith("BC"):
						if bc not in allbarcode:
							allbarcode = dict()
							allbarcode[bc] = "\t".join(saminfo1[2,3,7])
						else:
							allbarcode[bc] = allbarcode[bc] + ";" + "\t".join(saminfo1[2,3,7])
		

class OUTERSOFT:
	def bwa_barcode(self, barcode, outfile, bwa_path, bwa_parameter, ref_path):
		outdir = os.path.dirname(outfile)
		shell = os.path.join(outdir, "barcode.align.sh")
		oshell = open(shell, 'w')
		shell_line = " ".join([bwa_path, 'mem', ref_path, barcode, bwa_parameter, "| grep -v '^@SQ'>" + outfile + "\n"])
		oshell.write(shell_line)
		oshell.close()
		print(shell)
		subprocess.call(["sh", shell])
		return(outfile)

	def bwa_fq(fq1, fq2, outprefix, bwa_path, bwa_parameter, samtools_path):
		sai1 = outprefix + '.1.sai'
		sai2 = outprefix + '.2.sai'
		outbam = outprefix + '.bam'
		
		print >> sys.stderr, "fq alignment starts at %s" % (time.asctime())
		subprocess.call([bwa_path, 'aln', ref_path, fq1, bwa_parameter, '-f', sai1])
		subprocess.call([bwa_path, 'aln', ref_path, fq2, bwa_parameter, '-f', sai2])
		subprocess.call([bwa_path, 'sampe', ref_path, sai1, sai2, fq1, fq2, "|", samtools_path + 'view -h -S -b - >' + outbam])
		print >> sys.stderr, "fq alignment ends at %s" % (time.asctime())

	def sort_by_barcode(rawbam, sortedbam, samtools_path, samtools_parameter):
		print >> sys.stderr, "sorting barcode starts at %s" % (time.asctime())
		subprocess.call([samtools_path, 'sort', rawbam + samtools_parameter + '-O bam -o' + sortedbam])
		print >> sys.stderr, "sorting barcode ends at %s" % (time.asctime())

	def sort_by_chr(rawbam, sortedbam, samtools_path, samtools_parameter):
		print >> sys.stderr, "sorting bam starts at %s" % (time.asctime())
		subprocess.call([samtools_path, 'sort', rawbam + samtools_parameter + sortedbam])
		print >> sys.stderr, "sorting bam ends at %s" % (time.asctime())

	def picard_markdup(oribam, markedbam, picard_path, picard_paramter):
		oribam = "I=" + oribam
		markedbam = "O=" + markedbam
		metrics = 'M=' + markedbam
		metrics = metrics.replace('bam', 'metrics')

		print >> sys.stderr, "marking duplication starts at %s" % (time.asctime())
		subprocess.call([picard_path, oribam, markedbam, metrics, picard_paramter])
		print >> sys.stderr, "marking duplication ends at %s" % (time.asctime())

	def gatk_haplotypecaller(bam, outdir, gatk_path, gatk_parameter):

		print >> sys.stderr, "running haplotypecaller starts at %s" % (time.asctime())

		print >> sys.stderr, "running haplotypecaller ends at %s" % (time.asctime())

	def gatk_gvcftovcf(gvcf, vcf, gatk_path, gatk_parameter):

		print >> sys.stderr, "converting vcf starts at %s" % (time.asctime())

		print >> sys.stderr, "converting vcf ends at %s" % (time.asctime())

	def phasing(vcf, phase_path, phase_parameter):

		print >> sys.stderr, "phasing starts at %s" % (time.asctime())

		print >> sys.stderr, "phasing ends at %s" % (time.asctime())

def LRTK_usage():
	usage = \
	'''
	A pipeline to run 10x data
	Version: 0.01
	Date: 2017-06-01
	Author: meijp@foxmail.com
	Usage: python LRTK.py -i config.txt

	Options:
		-c configure file
		-i input fq, compressed and uncompressed fastq files are both availabe
		-o output dir

	'''
	print (usage)
	sys.exit(1)

if __name__ == '__main__':
	if len(sys.argv) <= 1:
		LRTK_usage()
	
	configurationFile = None
	step = '0-'
	qsub_q = None
	qsub_P = None
	qsub = None
	kept = None
	opts, args = getopt.gnu_getopt(sys.argv[1:], 'c:i:o:s:Qkf:')
	for o, a in opts:
		if o == '-c': configurationFile = a
		if o == '-i': inputfq = a
		if o == '-o': outputDir = a
		if o == '-s': step = a
		if o == '-Q': qsub = True
		if o == '-k': kept = True
		if o == '-f': Filter = True

	if configurationFile == None:
		print >> sys.stderr, "config.txt not found!"
		sys.exit(-1)
	elif os.path.isfile(configurationFile) and os.path.exists(configurationFile):
		pass
	else:
		print >> sys.stderr, "config.txt not found!"
		sys.exit(-1)

	if os.path.exists(outputDir):
		pass
	else:
		os.mkdir(outputDir)
	
	if os.path.exists(inputfq):
		pass
	else:
		print >> sys.stderr, "input fq not found!\n"

	_process_ = list()
	if '-' in step:
		sp = re.split('-', step)
		if len(sp) == 1:
			_process_ = range(int(sp[0]), 8)
		elif len(sp) == 2:
			_process_ = range(int(sp[0]), (int(sp[1]) + 1))
		else:
			print >> sys.stderr, "wrong command for -s: only single number(0~8), single number with '-'(0~8-), and two numbers seperated by '-'(0-8) are available"
			sys.exit(-1)
	else:
		if re.match(r'\d+', step):
			_process_.append(int(step))
		else:
			print >> sys.stderr, "wrong command for -s: only single number(0~8), single number with '-'(0~8-), and two numbers seperated by '-'(0-8) are available"
			sys.exit(-1)

	O = OUTERSOFT()
	G = baseinfo()
	G.get_config(configurationFile)

	inputfile = inputfq
	extract_barcode = None
	for p in _process_:
		if p == 0:
			R = extract_barcode()
			sys.stderr.write("\n... ... step 0 ... ... %s \n\n" % time.asctime())
			sys.stderr.write("[ %s ] extract original barcode info from %s\n" % (time.asctime(), inputfq))
			extract_barcode = R.revise(inputfile, outputDir)
			sys.stderr.write("[ %s ] original barcode had been extracted, output file has been written to: %s\n\n" % (time.asctime(), extract_barcode))
		elif p == 1:
			sys.stderr.write("\n... ... step 1 ... ... %s \n\n" % time.asctime())
			barcoderef = G.Barcode()
			bracode_aln_parameter = G.Barcode_aln_par()
			bwa_path = G.Bwa()
			outfile = inputfile.replace("gz", "sam")
			sys.stderr.write("[ %s ] align to the barcode set\n" % time.asctime())
			inputfile = O.bwa_barcode(inputfile, outfile, bwa_path, bracode_aln_parameter, barcoderef)
			sys.stderr.write("\n[ %s ] align to the barcode set, output sam file has been written to: %s\n" % (time.asctime(), inputfile))
		elif p == 2:
			sys.stderr.write("\n... ... step 2 ... ... %s \n\n" % time.asctime())
			M = modify_barcode()
			s
