import os, re, sys, gzip
import pysam

class LRTK:
	def __init__(self):
		self.configdist = {}

	def get_config(self, configFile):
		o_config = open (configFile, 'r')
		while True:
			config_line = o_config.readline().strip()
			if len(config_line) == 0:
				break
			elif config_line.startswith('#'):
				pass
			else:
				(name, path) = re.split("=", config_line)
				name = name.strip()
				path = path.strip()
				path = path.replace("'", "")
				path = path.replace('"', '')

				self.configdist[name] = path
				print (name, path)
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

	def Barcode_aln(self):
		abs_path = self.get_path('barcode_aln')
		return abs_path

class extract_barcode:
	def revise(rawfq, outdir):
		if rawfq.endswith('gz'):
			o_rawfq = gzip.open(rawfq, 'rb')
		else:
			o_rawfq = open(rawfq, 'r')

		(rawfq_path, rawfq_name) = os.path.split(rawfq)
		print(rawfq_name)

		b = rawfq_name
		if(rawfq_name.find("fastq")):
	#		rawfq_name.replace("fastq", "barcode", 1)
			strinfo = re.compile('fastq')
			b = strinfo.sub('barcode', rawfq_name)
		else:
	#		rawfq_name.replace("fq", "barcode", 1)
			strinfo = re.compile('fq')
			b = strinfo.sub('barcode', rawfq_name)

		outfile = os.path.join(rawfq_path, b)
		print(outfile)

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

class modify_barcode:
	def replace_barcode(rawfq, barcode_bam, prefix):
		L1 = LRTK()
		samfile = pysam.AlignmentFile(barcode_bam, "rb")
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
				if id1[0]= != barcodeinfo[0]:
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
								for i in range((len(alter)):
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
	def calculate(sorted_by_barcode_bam, outfile):
		samfile = pysam.AlignmentFile(sorted_by_barcode_bam, 'rb')

		allbarcode = dict()
		while True:
			saminfo1 = re.split("\t", samfile.readline().decode())
			saminfo2 = re.split("\t", samfile.readline().decode())

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
		
		lsls

class OUTERSOFT:
	def bwa_barcode(barcode, outfile, bwa_path, bwa_parameter, ref_path, samtools_path):
		print >> sys.stderr, "barcode alignment starts at %s" % (time.asctime())
		subprocess.call([bwa_path, 'mem', ref_path, barcode, bwa_parameter, "|", samtools_path + 'view -h -S -b - >' + outfile])
		print >> sys.stderr, "barcode alignment ends at %s" % (time.asctime())

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
		-i 

	'''
	print (usage)
	sys.exit(1)

if __name__ == '__main__':
	if len(sys.argv) <= 1:
		LRTK_usage()
	
	configurationFile = None
	qsub_q = None
	qsub_P = None
	qsub = None
	kept = None
	opts, args = getopt.gnu_getopt(sys.argv[1:], 'c:i:o:Qkf:')
	for o, a in opts:
		if o == '-c': configurationFile = a
		if o == '-i': inputDir = a
		if o == '-o': outputDir = a
		if o == '-Q': qsub = True
		if o == '-k': kept = True
		if o == '-f': filter = True

	if os.path.isfile(configurationFile) and os.path.exists(configurationFile):
		pass
	else:
		print >> sys.stderr, "config.txt not found!"
		sys.exit(-1)

	if os.path.exists(outputDir):
		pass
	else:
		os.mkdir(outputDir)
	
	if os.path.exists(inputDir):
		pass
	else:
		print >> sys.stderr, "input dir not found!\n"
	


if len(sys.argv) == 5:
	revise(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
else:
	print("python", sys.argv[0], "10x.fq(fastq)[.gz] outdir bwa_path bwa_parameter")
