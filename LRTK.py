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

	def Fq_aln_par(self):
		abs_path = self.get_path('fq_aln_parameter')
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

	def Samtools(self):
		abs_path = self.get_path('samtools')
		return abs_path

	def Java(self):
		abs_path = self.get_path('java')
		return abs_path

	def Picard(self):
		abs_path = self.get_path('picard')
		return abs_path

	def Picardparameter(self):
		abs_path = self.get_path('picard_parameter')
		return abs_path

	def Dbsnp(self):
		abs_path = self.get_path('dbsnp')
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

				barcode_line = rawfq_id1 + rawfq_seq1[0:15] + '\n' + rawfq_str1 + rawfq_qua1[0:15] + "\n"
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
					newfq_qua1 = realbarcodequa + rawfq_qua1[23:]
					newfq1 = rawfq_id1 + "\n" + newfq_seq1 + "\n" + rawfq_str1 + "\n" + newfq_qua1 + "\n"
					newfq_seq2 = str(realbarcode) + str(rawfq_seq2)
					newfq_qua2 = str(realbarcodequa) + str(rawfq_qua2)
					newfq2 = "\n".join([rawfq_id2, newfq_seq2, rawfq_str2, newfq_qua2, ""])

					newfq1_out.write(newfq1.encode())
					newfq2_out.write(newfq2.encode())

		newfq1_out.close()
		newfq2_out.close()

		return(prefix)

	def fill_barcode(self, original_sam, samtools_path, outdir = './'):
		original_name = os.path.basename(original_sam)

		new_sam = original_name.replace('sam', 'new.sam')
		new_sam = os.path.join(outdir, new_sam)
		new_bam = new_sam.replace('sam', 'bam')
		sorted_bam = new_sam.replace('new.sam', 'sorted')

		print(" ".join(["modified and sorted bam file:", sorted_bam + '.bam', "\n"]))
		
		o_original_sam = open(original_sam, 'r')
		o_new_sam = open(new_sam, 'w')

		while True:
			saminfo = o_original_sam.readline()
			if len(saminfo) == 0:
				break
			else:
				if saminfo.startswith('@'):
					o_new_sam.write(saminfo)
				else:
					saminfo = saminfo.strip()
					saminfolist = re.split('\t', saminfo)
					for s in range(len(saminfolist)):
						if 'BC:Z' in saminfolist[s]:
							saminfolist[s] = saminfolist[s].upper()
					saminfo = "\t".join(saminfolist) + '\n'
					o_new_sam.write(saminfo)

		o_original_sam.close()
		o_new_sam.close()

		shell = sorted_bam + ".sh"
		oshell = open(shell, 'w')
		shell_line = " ".join([samtools_path, "view -h -S -b", new_sam, ">", new_bam, "\n"])
		oshell.write(shell_line)
		shell_line = " ".join([samtools_path, "sort -m 1G", new_bam, sorted_bam, "\n"])
		oshell.write(shell_line)
		shell_line = " ".join([samtools_path, "view -h", sorted_bam + '.bam > ', sorted_bam + '.sam\n'])
		oshell.write(shell_line)
		shell_line = " ".join(["rm", original_sam, new_sam, new_bam, sorted_bam + '.sam', "\n"])
		oshell.write(shell_line)
		oshell.close()

		subprocess.call(["sh", shell])
		return(sorted_bam + '.bam', sorted_bam + '.sam')

class CFCR:
	def split_and_sort_sam(self, sorted_bam, outdir = './'):
		splitdir = os.path.join(outdir, "SamByChr")
		if os.path.exists(splitdir):
			pass
		else:
			os.path.mkdir(splitdir)

		sorted_sam = sorted_bam.replace('bam', 'sam')
		tmpshell = os.path.join(splitdir, "bam2sam.sh")
		wtmpshell = open(tmpshell, 'w')
		shell_line = " ".join([samtools_path, "view -h", sorted_bam, ">", sorted_sam, "\n"])
		wtmpshell.write(shell_line)
		wtmpshell.close()

		chrdict = dict()
		splitsamlist = list()
		headersam = os.path.join(splitdir, 'header.sam')
		outputsam = open(headersam, 'w')
		rsorted_sam = open(sorted_sam, 'r')

		issort = None
		while True:
			saminfo = rsorted_sam.readline()
			if len(saminfo) == 0:
				break
			elif saminfo.startswith('@'):
				if 'coordinate' in saminfo:
					issort = 1
				outputsam.write(saminfo)
			else:
				saminfolist = re.split("\t", saminfo)
				if saminfolist[2] == '*':
					pass
				elif len(chrdict) == 0 or saminfolist[2] not in chrdict:
					if len(chrdict) == 0 and issort == None:
						sys.stderr.write("\n\nWarning: the inputted sam file might have not been sorted %s\n\n" % sorted_sam)
					chrdict[saminfolist[2]] = saminfolist[2]
					outputsam.close()
					chrsam = os.path.join(splitdir, saminfolist[2] + '.sam')
					outputsam = open(chrsam, 'w')
					splitsamlist.append(chrsam)
					sys.stderr.write("[ %s ] split sorted sam file, processing %s: %s ...\n" % (time.asctime(), saminfolist[2], chrsam))
					outputsam.write(saminfo)
				else:
					outputsam.write(saminfo)
		outputsam.close()

		samfilepath = os.path.join(splitdir, "sam.txt")
		wsamfilepath = open(samfilepath, 'w')
		for chrsam in splitsamlist:
			chrsortedsam = chrsam.replace('sam', 'sorted_by_barcode.sam')
			sys.stderr.write("[ %s ] sort sam file, processing %s ...\n" % (time.asctime(), chrsortedsam))
			subprocess.call(["sort -k12", chrsam, ">", chrsortedsam])
			subprocess.call(["rm", chrsam])
			wsamfilepath.write("\n".join([chrsortedsam, ""]))
		wsamfilepath.close()
		return(wsamfilepath)

	def calculate(self, sorted_by_barcode_sam_list, outdir = './'):
		samfile_list = open(sorted_by_barcode_sam, 'r')
		marked_molecule = os.path.join(outdir, "molecule.txt")
		wmarked_molecule = open(marked_molecule, 'w')
		allbarcode = dict()
		mol_id = 1
		for sam in samfile_list:
			sam = sam.strip()
			samfile = open(sam, 'r')
			while True:
				saminfo = samfile.readline().strip()
				saminfolist = re.split("\t", samfile)

				if len(saminfo) == 0:
					break
				elif saminfolist[8] == 0:
					pass

				for bc in saminfolist[11:]:
					if bc.startswith("BC"):
						if bc not in allbarcode:
							if len(allbarcode) == 0:
								pass
							else:
								mol_id = self.FR.(allbarcode[-1], mol_id, wmarked_molecule)
								allbarcode = dict()
							allbarcode[bc] = saminfo
						else:
							allbarcode[bc] = allbarcode[bc] + "\n" + saminfo
			samfile.close()
		mol_id = self.FR.(allbarcode[-1], mol_id, wmarked_molecule)
		wmarked_molecule.close()

	def FR(self, barcodeinfo, molecule_id, writemarker):
		readid = dict()
		onebarcodelist = re.split("\n", barcodeinfo)
		for onebarcode in onebarcodelist:
			eachlist = re.split("\t", onebarcode)
			(eachid, chr_id, start, readlength) = eachlist[0, 2, 3, 5]
			readlength = int(readlength.strip())
			start = int(start)
			end = start + readlength - 1
			reads = onebarcode
			### merge pair-end reads into single fragment
			if eachid not in readid:
				readid[eachid] = [chr_id, start, end, idinfo]
			else:
				if readid[eachid][0] != chr_id:
					sys.stderr.write("%s pair-end reads mapped to different chromosom, would be ignored.\n" % eachid)
					readid.pop(eachid)
				else:
					start = min(readid[eachid][1], start)
					end = max(readid[eachid][2], end)
					dis = end - start + 1
					if dis > 10000:
						sys.stderr.write("%s fragment length > 10kb, would be ignored.\n" % eachid)
						readid.pop(eachid)
					else:
						reads = reads + "\n" + onebarcode
						readid[eachid] = [chr_id, start, end, reads]

		molecule = dict()
		s = int(molecule_id) + 1
		molecule_id = int(molecule_id) + 1
		e = s + 1
		for key, value in readid:
			if molecule_id not in molecule:
				molecule[molecule_id]["chrid"] = value[0]
				molecule[molecule_id]["start"] = value[1]
				molecule[molecule_id]["end"] = value[2]
				molecule[molecule_id]["reads"] = value[3]
			else:
				isnew = 1
				for i in range(s, e):
					if isnew == 0:
						pass
					else:
						if molecule[i]["chrid"] == value[0]:
							dist = molecule[i]["start"] - value[1]
							if dist < 50000:
								isnew = 0
								molecule[i]["start"] = min(molecule[i]["start"], value[1])
								molecule[i]["end"] = max(molecule[i]["end"], value[2])
								molecule[i]["reads"] = molecule[i]["reads"] + "\n" + value[3]
				if isnew == 1:
					molecule_id = molecule_id + 1
					molecule[molecule_id]["chrid"] = value[0]
					molecule[molecule_id]["start"] = value[1]
					molecule[molecule_id]["end"] = value[2]
					molecule[molecule_id]["reads"] = value[3]
					s = s + 1

		for i in range(s, e):
			readlist = re.split("\n", molecule[i]["reads"])
			for rl in readlist:
				outinfo = i + "\t" + rl + "\n"
				writemarker.write(outinfo)

		return(molecule_id)

		

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

	def bwa_fq(self, fq1, fq2, outprefix, bwa_path, bwa_parameter):
		sai1 = outprefix + '.1.sai'
		sai2 = outprefix + '.2.sai'
		outbam = outprefix + '.sam'
		
		sys.stderr.write("[ %s ] fq alignment starts\n" % (time.asctime()))
		outdir = os.path.dirname(fq1)
		shell = os.path.join(outdir, "fq.aln.sh")
		oshell = open(shell, 'w')
		shell_line = " ".join([bwa_path, 'aln', ref_path, fq1, bwa_parameter, '-f', sai1, '&\n'])
		oshell.write(shell_line)
		shell_line = " ".join([bwa_path, 'aln', ref_path, fq2, bwa_parameter, '-f', sai2, '&\nwait\n'])
		oshell.write(shell_line)
		#shell_line = " ".join([bwa_path, 'sampe', ref_path, sai1, sai2, fq1, fq2, "|", samtools_path + ' view -h -S -b - > ' + outbam + "\n"])
		shell_line = " ".join([bwa_path, 'sampe', ref_path, sai1, sai2, fq1, fq2, ">", outbam + "\nrm", sai1, sai2, "\n"])
		oshell.write(shell_line)
		oshell.close()
		print(shell)

		subprocess.call(["sh", shell])
		sys.stderr.write("[ %s ] fq alignment finish\n" % (time.asctime()))
		return(outbam)

	def picard_markdup(self, oribam, markedbam, java_path, picard_path, picard_paramter):
		metrics = markedbam.replace('bam', 'metrics.txt')

		bai = oribam.replace('bam', 'bai')
		isbai = 0
		if os.path.exists(bai):
			isbai = 1
		else:
			bai = oribam + '.bai'
			if os.path.exists(bai):
				isbai = 1
		if isbai:
			pass
		else:
			sys.stderr.write("bai not found, will produce bai file for bam automatically: %s\n" % bai)
			tmpshell = bai + '.sh'
			wtmpshell = open(tmpshell, 'w')
			shell_line = " ".join([java_path, '-jar', picard_path, 'BamIndexStats', "I=" + oribam, "O=" + bai, "\n"])
			wtmpshell.write(shell_line)
			wtmpshell.close()
			subprocess.call(["sh", tmpshell])
			subprocess.call(["rm", tmpshell])

		sys.stderr.write("[ %s ] marking duplication ... \n" % (time.asctime()))
		markedbamdir = os.path.dirname(markedbam)
		tmpshell = os.path.join(markedbamdir, 'mark_duplicate.sh')
		wtmpshell = open(tmpshell, 'w')
		shell_line = " ".join([java_path, '-jar', picard_path, 'MarkDuplicates', "I=" + oribam, "O=" + markedbam, "M=" + metrics, picard_paramter, "\n"])
		wtmpshell.write(shell_line)
		wtmpshell.close()
		subprocess.call(["sh", tmpshell])
		sys.stderr.write("[ %s ] marked bam file has been written to %s\n" % (time.asctime(), markedbam))
		return(markedbam)

	def phasing(self, vcf, phase_path, phase_parameter):

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

	def isset(v):
		try:
			type(eval(v))
		except:
			return 0
		else:
			return 1

	configurationFile = None
	step = '0-'
	inputfq = None
	newprefix = None
	barcodesam = None
	Sortedsam = None
	Sortedbam = None
	Samfilepath = None
	Markedbam = None
	Chrlist = None

	opts, args = getopt.gnu_getopt(sys.argv[1:], 'c:f:o:s:S:e:p:M:P:k:L')
	for o, a in opts:
		if o == '-c': configurationFile = a
		if o == '-f': inputfq = a
		if o == '-o': outputDir = a
		if o == '-s': barcodesam = a
		if o == '-S': step = a
		if o == '-e': extract_barcode = a
		if o == '-p': newprefix = a
		if o == '-m': Sortedsam = a
		if o == '-M': Sortedbam = a
		if o == '-P': Samfilepath = a
		if o == '-k': Markedbam = a
		if o == '-L': Chrlist = a


	if configurationFile == None:
		sys.stderr.write("config.txt not found!\n")
		sys.exit(-1)
	elif os.path.isfile(configurationFile) and os.path.exists(configurationFile):
		pass
	else:
		sys.stderr.write("config.txt not found!\n")
		sys.exit(-1)

	if os.path.exists(outputDir):
		pass
	else:
		os.mkdir(outputDir)

	_process_ = list()
	if '-' in step:
		sp = re.split('-', step)
		if len(sp) == 1:
			_process_ = range(int(sp[0]), 8)
		elif len(sp) == 2:
			_process_ = range(int(sp[0]), (int(sp[1]) + 1))
		else:
			sys.stderr.write("wrong command for -S: only single number(0~8), single number with '-'(0~8-), and two numbers seperated by '-'(0-8) are available\n")
			sys.exit(-1)
	else:
		if re.match(r'\d+', step):
			_process_.append(int(step))
		else:
			sys.stderr.wirte("wrong command for -S: only single number(0~8), single number with '-'(0~8-), and two numbers seperated by '-'(0-8) are available\n")
			sys.exit(-1)

	O = OUTERSOFT()
	G = baseinfo()
	G.get_config(configurationFile)

	for p in _process_:
		if p == 0:
			sys.stderr.write("\n... ... step 0 ... ... %s \n\n" % time.asctime())
			if os.path.exists(inputfq):
				pass
			else:
				sys.stderr.write("input fq not found!\n")
				sys.exit(-1)

			R = extract_barcode()
			sys.stderr.write("[ %s ] extract original barcode info from %s\n" % (time.asctime(), inputfq))
			extract_barcode = R.revise(inputfq, outputDir)
			sys.stderr.write("[ %s ] original barcode had been extracted, output file has been written to: %s\n\n" % (time.asctime(), extract_barcode))
		elif p == 1:
			sys.stderr.write("\n... ... step 1 ... ... %s \n\n" % time.asctime())
			if os.path.exists(extract_barcode):
				pass
			else:
				sys.stderr.write("\nstep 0 is not finished, please run step 0 again or provide extraced bracode info of the original 10x fq using '-e'\n\n")
				sys.exit(-1)

			barcoderef = G.Barcode()
			bracode_aln_parameter = G.Barcode_aln_par()
			bwa_path = G.Bwa()
			barcodesam = extract_barcode.replace("gz", "sam")
			sys.stderr.write("[ %s ] align to the barcode set\n" % time.asctime())
			barcodesam = O.bwa_barcode(extract_barcode, barcodesam, bwa_path, bracode_aln_parameter, barcoderef)
			sys.stderr.write("\n[ %s ] align to the barcode set, output sam file has been written to: %s\n" % (time.asctime(), barcodesam))
		elif p == 2:
			sys.stderr.write("\n... ... step 2 ... ... %s \n\n" % time.asctime())
			if inputfq != None and os.path.exists(inputfq):
				pass
			else:
				sys.stderr.write("\ninput fq not found, please provide the original 10x fq file using '-f'\n\n")
				sys.exit(-1)
			if barcodesam != None and os.path.exists(barcodesam):
				pass
			else:
				sys.stderr.write("\nstep 1 is not finished, please run step 1 again or provide barcode sam file using '-s'\n\n")
				sys.exit(-1)

			M = modify_barcode()
			if isset(newprefix):
				pass
			else:
				newprefix = os.path.basename(inputfq)
				if newprefix.find("fastq"):
					strinfo = re.compile('fastq')
					newprefix = strinfo.sub('new.fq', newprefix)
				else:
					strinfo = re.compile('fq')
					newprefix = strinfo.sub('new.fq', newprefix)
				newprefix = newprefix.replace('.gz', '')
				newprefix = os.path.join(outputDir, newprefix)
			sys.stderr.write("[ %s ] replace barcode info\n" % time.asctime())
			newprefix = M.replace_barcode(inputfq, barcodesam, newprefix)
			sys.stderr.write("[ %s ] barcode has been replaced, and new fq files have been written to %s_[1/2].fq.gz\n" % (time.asctime(), newprefix))
		elif p == 3:
			sys.stderr.write("\n... ... step 3 ... ... %s \n\n" % time.asctime())
			if newprefix == None:
				sys.stderr.write("\nstep 2 is not finished, please run step 2 again or provide the prefix of the paired fq files unsing '-p'\n\n")
				sys.exit(-1)
			else:
				fq1 = newprefix + "_1.fq.gz"
				fq2 = newprefix + "_2.fq.gz"
				if os.path.exists(fq1) and os.path.exists(fq2):
					pass
				else:
					sys.stderr.write("\nstep 2 is not finished, please run step 2 again or provide the prefix of the paired fq files unsing '-p'\n\n")
					sys.exit(-1)

			
			ref_path = G.Ref()
			fq_aln_parameter = G.Fq_aln_par()
			bwa_path = G.Bwa()
			samtools = G.Samtools()

			sys.stderr.write("[ %s ] align to the genome %s\n" % (time.asctime(), ref_path))
			original_sam = O.bwa_fq(fq1, fq2, newprefix, bwa_path, fq_aln_parameter)
			sys.stderr.write("[ %s ] alignment finished. Original sam file has been writeen to %s\n\n" % (time.asctime(), original_sam))

			M = modify_barcode()
			samtools = G.Samtools()
			samdir = os.path.dirname(original_sam)
			(Sortedbam, Sortedsam) = M.fill_barcode(original_sam, samtools, samdir)
		elif p == 4:
			sys.stderr.write("\n... ... step 4 ... ... %s \n\n" % time.asctime())
			if Sortedbam == None:
				sys.stderr.write("\nstep 3 is not finished, please run step 3 again or provide sorted bam files unsing '-M'\n\n")
				sys.exit(-1)
			elif os.path.exists(Sortedbam):
				pass
			else:
				sys.stderr.write("%s is not exists!\n" % (Sortedbam))
				sys.exit(-1)

			javapath = G.Java()
			picardpath = G.Picard()
			picardparameter = G.Picardparameter()
			Markedbam = Sortedbam.replace('bam', 'marked.bam')
			sys.stderr.write("[ %s ] mark duplication ... \n" % time.asctime())
			Markedbam = O.picard_markdup(Sortedbam, Markedbam, javapath, picardpath, picardparameter)
			if os.path.exists(Markedbam):
				sys.stderr.write("[ %s ] duplication reads have been marked, and the output has been written to %s \n" % (time.asctime(), Markedbam))
			else:
				sys.stderr.write("\nERROR: %s has not been created, please check the warning info and try again\n\n" % Markedbam)
				sys.exit(-1)
		elif p == 5:
			sys.stderr.write("\n... ... step 5 ... ... %s \n\n" % time.asctime())
			if Sortedbam == None:
				sys.stderr.write("\nstep 4 is not finished, please run step 4 again or provide marked bam file unsing '-k'\n\n")
				sys.exit(-1)
			elif os.path.exists(Sortedbam):
				pass
			else:
				sys.stderr.write("%s is not exists!\n" % (Sortedbam))
				sys.exit(-1)
			
			C = CFCR()
			samdir = os.path.dirname(Sortedbam)
			sys.stderr.write("[ %s ] split sam file by chr, sort sam file by barcode ... \n" % time.asctime())
			Samfilepath = C.split_and_sort_sam(Sortedbam, samdir)
			sys.stderr.write("[ %s ] splited and sorted sam file list: %s \n\n" % (time.asctime(), Samfilepath))

			sys.stderr.write("[ %s ] calculating CF/CR ... \n" % time.asctime())
			stattxt = C.calculate(Samfilepath, samdir)
			sys.stderr.write("[ %s ] CF/CR has been written to %s \n" % (time.asctime(), stattxt))
		elif p == 6:
			sys.stderr.write("\n... ... stpe 6 ... ... %s \n\n" % time.asctime())
			if Markedbam == None:
				sys.stderr.write("\nstep 5 is not finished, please run step 5 again or provide marked bam files unsing '-k'\n\n")
				sys.exit(-1)
			elif os.path.exists(Markedbam):
				pass
			else:
				sys.stderr.write("%s is not exists!\n" % (Markedbam))
				sys.exit(-1)

			if Chrlist == None:
				sys.stderr.write("\nprovide chr list unsing '-L'\n\n")
				sys.exit(-1)
			elif os.path.exists(Chrlist):
				pass
			else:
				sys.stderr.write("%s is not exists!\n" % Chrlist)
				sys.exists(-1)

			javapath = G.Java()
			gatkpath = G.Gatk()
			ref = G.Ref()
			dbsnp = G.Dbsnp()

			bamdir = os.path.dirname(Markedbam)
			vcfdir = os.path.join(bamdir, "vcf")
			if os.path.exists(vcfdir):
				pass
			else:
				os.mkdir(vcfdir)

			shelldir = op.path.join(vcfdir, "shell")
			if os.path.exists(shelldir):
				pass
			else:
				os.mkdir(vcfdir)

			rChrlist = open(Chrlist, 'r')
			for ch in rChrlist:
				ch = ch.strip()
				chrshell = os.path.join(shelldir, ch + ".vcf.sh")
				chrgvcf = os.path.join(vcfdir, ch + ".gvcf")
				chrvcf = os.path.join(vcfdir, ch + ".vcf")
				wchrshell = open(chrshell, 'w')
				shell_line = " ".join(["set -e\n", javapath, "-Djava.io.tmpdir=" + vcfdir, "-jar", gatkpath, "-T HaplotypeCaller -R", ref, "-I", Markedbam, "-U --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 --dbsnp", dbsnp, "-L", ch, "-o", chrgvcf, "\n"])
				wchrshell.write(shell_line)
				shell_line = " ".join([javapath, "-Xmx2g -jar", gatkpath, "-R", ref, "-T GenotypeGVCFs --variant", chrgvcf, "-o", chrvcf, "--dbsnp", dbsnp, "-allSites -L", ch, "\n"])
				wchrshell.write(shell_line)
				wchrshell.close()

				sys.stderr.write("sh %s\n" % chrshell)

