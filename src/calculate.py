import os, sys, gzip
import getopt
import time
import subprocess
import re
from collections import defaultdict
import random

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

	def Samtools(self):
		abs_path = self.get_path('samtools')
		return abs_path

	def Barcode_index(self):
		abs_path = self.get_path('barcode_index')
		return abs_path

	def Genomesize(self):
		abs_path = self.get_path('genomesize')
		return abs_path

class CFCR:
	def find_samtools_version(self, samtools_path):
		randomstring = ''.join(random.sample(string.ascii_letters + string.digits, 8))
		tmpshell = "~/" + randomstring + ".sh"
		tmplog = "~/" + randomstring + ".log"
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
		subprocess.call(["rm -rf", tmpshell, tmplog])
		print(int(sv))
		return(int(sv))

	def split_and_sort_sam(self, sorted_bam, samtools_path, barcode_index, outdir = './'):
		splitdir = outdir + "/tmp/SamByChr"
		if os.path.exists(splitdir):
			pass
		else:
			os.mkdir(splitdir)

		sorted_sam = sorted_bam
		if sorted_bam.endswith("bam"):
			sorted_sam = sorted_bam.replace('bam', 'sam')
			tmpshell = os.path.join(splitdir, "bam2sam.sh")
			wtmpshell = open(tmpshell, 'w')
			sortshell = tmpshell + ".tmp"
			wsortshell = open(tmpshell + ".tmp", 'w')
			tmpheader = splitdir + "/tmp.header"
			shell_line = " ".join([samtools_path, "view -H", sorted_bam, "| grep coordinate >", tmpheader, "\n"])
			wsortshell.write(shell_line)
			wsortshell.close()
			subprocess.call(["sh", sortshell])
			if os.path.getsize(tmpheader):
				sys.stderr.write("%s has been sorted!\n" % sorted_bam)
			else:
				sys.stderr.write("%s has not been sorted, and would be sorted before calculating CF/CR!\n" % sorted_bam)

				sv = self.find_samtools_version(samtools_path)
				if sv == 0:
					shell_line = " ".join([samtools_path, "sort -m 2G", sorted_bam, sorted_bam + ".sorted\n"])
				else:
					shell_line = " ".join([samtools_path, "sort -m 2G -o", sorted_bam + ".sorted.bam -T", sorted_bam, sorted_bam, "\n"])
				wtmpshell.write(shell_line)
				shell_line = " ".join(["mv", sorted_bam + ".sorted.bam", sorted_bam, "\n"])
				wtmpshell.write(shell_line)
			shell_line = " ".join([samtools_path, "view -h", sorted_bam, ">", sorted_sam, "\n"])
			wtmpshell.write(shell_line)
			wtmpshell.close()
			subprocess.call(["sh", tmpshell])

		chrdict = dict()
		splitsamlist = list()
		headersam = os.path.join(splitdir, 'header.sam')
		outputsam = open(headersam, 'w')
		rsorted_sam = open(sorted_sam, 'r')
		while True:
			saminfo = rsorted_sam.readline()
			if len(saminfo) == 0:
				break
			elif saminfo.startswith('@'):
				outputsam.write(saminfo)
			else:
				saminfo = saminfo.strip()
				saminfolist = re.split("\t", saminfo)
				if saminfolist[2] == '*':
					pass
				elif len(chrdict) == 0 or saminfolist[2] not in chrdict:
					chrdict[saminfolist[2]] = saminfolist[2]
					outputsam.close()
					chrsam = os.path.join(splitdir, saminfolist[2] + '.sam')
					outputsam = open(chrsam, 'w')
					splitsamlist.append(chrsam)
					sys.stderr.write("[ %s ] split sorted sam file, processing %s: %s ...\n" % (time.asctime(), saminfolist[2], chrsam))
					### move the barcode info to the 12th column
					saminfo = self.modify_saminfo(saminfo)
					outputsam.write(saminfo)
				else:
					saminfo = self.modify_saminfo(saminfo)
					outputsam.write(saminfo)
		outputsam.close()
		
		samfilepath = os.path.join(splitdir, "sam.txt")
		wsamfilepath = open(samfilepath, 'w')
		for chrsam in splitsamlist:
			chrsortedsam = chrsam.replace('sam', 'sorted_by_barcode.sam')
			sys.stderr.write("[ %s ] sort sam file, processing %s ...\n" % (time.asctime(), chrsortedsam))
			tmpshell = os.path.join(splitdir, "sort.sh")
			wtmpshell = open(tmpshell, "w")
			shell_line = " ".join(["sort", "-k12,12 -k4n", chrsam, ">", chrsortedsam, "\n"])
			wtmpshell.write(shell_line)
			shell_line = " ".join(["rm", chrsam, "\n"])
			wtmpshell.write(shell_line)
			wtmpshell.close()
			subprocess.call(["sh", tmpshell])
			wsamfilepath.write("\n".join([chrsortedsam, ""]))
		wsamfilepath.close()

		subprocess.call(["rm", sorted_sam])

		return(samfilepath)

	def modify_saminfo(self, orisaminfo):
		saminfolist = re.split("\t", orisaminfo)
		bcindex = 1
		for n in range(len(saminfolist)):
			if saminfolist[n].startswith("BX:Z"):
				bcindex = n
		if bcindex == 11 or bcindex == 1:
			saminfo = "\t".join(saminfolist)
		else:
			newsaminfolist = saminfolist[0:11]
			newsaminfolist.append(saminfolist[bcindex])
			for n in range(11, bcindex):
				newsaminfolist.append(saminfolist[n])
			if bcindex == (len(saminfolist) - 1):
				pass
			else:
				for n in range(bcindex+1, len(saminfolist)):
					newsaminfolist.append(saminfolist[n])
			saminfo = "\t".join(newsaminfolist)
		saminfo = saminfo + "\n"
		return(saminfo)

	def calculate(self, sorted_by_barcode_sam_list, target_size, outdir = './'):
		samfile_list = open(sorted_by_barcode_sam_list, 'r')
		marked_molecule = os.path.join(outdir, "molecule.full.gz")
		wmarked_molecule = gzip.open(marked_molecule, 'wb')
		broken_molecule = os.path.join(outdir, "molecule.broken.txt")
		wbroken_molecule = open(broken_molecule, 'w')

		mol_id = 1
		all_CF = 0
		all_CR = 0
		allbarcode = dict()
		readid_dict = defaultdict(int)
		molecule = defaultdict(list)
		molecule2 = defaultdict(list)
		molecule3 = defaultdict(int)
		molecule4 = defaultdict(list)
		for sam in samfile_list:
			sam = sam.strip()
			sys.stderr.write("[ %s ] processing: %s \n" %(time.asctime(), sam))
			samfile = open(sam, 'r')
			while True:
				saminfo = samfile.readline().strip()
				saminfolist = re.split("\t", saminfo)

				if len(saminfo) == 0:
					break
				elif saminfolist[5] == "*":
					pass
				else:
					### 以barcode 为单位进行处理
					barcode_field = [s for s in saminfolist if "BX:Z:" in s]
					if barcode_field != []:
						bc = barcode_field[0].split(":")[2].split("-")[0]
						start_pos = int(saminfolist[3])
						readlength = len(saminfolist[9])
						end_pos = start_pos + readlength - 1
						chrid = saminfolist[2]
						if bc not in allbarcode:	## 如果barcode不存在，则有两种情况：1. allbarcode是空的   2. allbarcode不是空的，有且只有一个key
							if len(allbarcode) == 0:	## allbarcode是空的，只赋值，不输出
								allbarcode[bc] = bc
								readid_dict[saminfolist[0]] = 1
								molecule[mol_id].append(start_pos)
								molecule2[mol_id].append(end_pos)
								molecule3[mol_id] = readlength
								molecule4[mol_id].append(saminfo)
							else:	## allbarcode不是空的，原来存的信息要输出出来，清空后，重新赋值
								se_num = 0
								pe_num = 0
								for readid_ in readid_dict.keys():
									if readid_dict[readid_] == 2:
										pe_num += 1
									elif readid_dict[readid_] == 1:
										se_num += 1
								if pe_num >= 2 or se_num >= 4 or (pe_num >= 1 and se_num >= 2):
									barcodeinfo = sorted(allbarcode.keys())
									minpos = min(molecule[mol_id])
									maxpos = max(molecule2[mol_id])
									all_CR += molecule3[mol_id]
									small_cf = maxpos - minpos + 1
									all_CF += small_cf
									qn = len(molecule[mol_id])
									for eachinfo in molecule4[mol_id]:
										outmoleculeinfo = "\t".join([chrid, str(minpos), str(maxpos), str(small_cf), str(qn), barcodeinfo[0], str(mol_id), eachinfo]) + "\n"
										wmarked_molecule.write(outmoleculeinfo.encode())
									mol_id += 1
								else:
									for readid_ in readid_dict.keys():
										broken_info = "\t".join([str(readid_), str(readid_dict[readid_])]) + "\n"
										wbroken_molecule.write(broken_info)
								allbarcode = dict()
								readid_dict = defaultdict(int)
								molecule = defaultdict(list)
								molecule2 = defaultdict(list)
								molecule3 = defaultdict(int)
								molecule4 = defaultdict(list)
								allbarcode[bc] = bc
								readid_dict[saminfolist[0]] = 1
								molecule[mol_id].append(start_pos)
								molecule2[mol_id].append(end_pos)
								molecule3[mol_id] = readlength
								molecule4[mol_id].append(saminfo)
						else:
							dist = start_pos - molecule[mol_id][-1]
							if dist < 50000:
								readid_dict[saminfolist[0]] += 1
								molecule[mol_id].append(start_pos)
								molecule2[mol_id].append(end_pos)
								molecule3[mol_id] += readlength
								molecule4[mol_id].append(saminfo)
							else:
								se_num = 0
								pe_num = 0
								for readid_ in readid_dict.keys():
									if readid_dict[readid_] == 2:
										pe_num += 1
									elif readid_dict[readid_] == 1:
										se_num += 1
								if pe_num >= 2 or se_num >= 4 or (pe_num >= 1 and se_num >= 2):
									minpos = min(molecule[mol_id])
									maxpos = max(molecule2[mol_id])
									all_CR += molecule3[mol_id]
									small_cf = maxpos - minpos + 1
									all_CF += small_cf
									qn = len(molecule[mol_id])
									barcodeinfo = sorted(allbarcode.keys())
									for eachinfo in molecule4[mol_id]:
										outmoleculeinfo = "\t".join([chrid, str(minpos), str(maxpos), str(small_cf), str(qn), barcodeinfo[0], str(mol_id), eachinfo]) + "\n"
										wmarked_molecule.write(outmoleculeinfo.encode())
									mol_id += 1
								else:
									for readid_ in readid_dict.keys():
										broken_info = "\t".join([str(readid_), str(readid_dict[readid_])]) + "\n"
										wbroken_molecule.write(broken_info)

								allbarcode = dict()
								readid_dict = defaultdict(int)
								molecule = defaultdict(list)
								molecule2 = defaultdict(list)
								molecule3 = defaultdict(int)
								molecule4 = defaultdict(list)
								allbarcode[bc] = bc
								readid_dict[saminfolist[0]] = 1
								molecule[mol_id].append(start_pos)
								molecule2[mol_id].append(end_pos)
								molecule3[mol_id] = readlength
								molecule4[mol_id].append(saminfo)

			samfile.close()
		wmarked_molecule.close()
		samfile_list.close()
		wbroken_molecule.close()

		cfcr_stat_file = os.path.join(outdir, "CFCR.stat")
		wcfcr_stat_file = open(cfcr_stat_file, 'w')
		stat_line = "Item" + "\t" + "CF/CR length" + "\t" + "genome_size" + "\t" + "CF/CR depth\n"
		wcfcr_stat_file.write(stat_line)
		target_size = int(target_size)
		CF_depth = all_CF / float(target_size)
		CR_depth = all_CR / float(all_CF)

		stat_line = "CF:" + "\t" + str(all_CF) + "\t" + str(target_size) + "\t" + str(CF_depth) + "\n"
		wcfcr_stat_file.write(stat_line)
		stat_line = "CR:" + "\t" + str(all_CR) + "\t" + str(all_CF) + "\t" + str(CR_depth) + "\n"
		wcfcr_stat_file.write(stat_line)
		wcfcr_stat_file.close()
		return(cfcr_stat_file)

def usage():
	calculation_usage = \
	'''
	calculate CF and CR
	Version: 1.0.0
	Dependents: Python (>=3.0), SAMtools
	Last Updated Date: 2017-06-01
	Contact: meijp@foxmail.com

	Usage: python calculate.py <options>

	Options:
		-i --input, path of input bam file
		-o --outputdir, the path of output directory
		-c --config, the path of configuration file [default: outdir/config/QC.config]
		-h --help, help info

	'''
	print(calculation_usage)

if __name__ == '__main__':
	if len(sys.argv) < 5:
		usage()
		sys.exit(-1)

	inputbam = None
	outputdir = None
	ConfigFile = None
	opts, args = getopt.gnu_getopt(sys.argv[1:], 'i:o:c:h:', ['input', 'outputdir', 'config', 'help'])
	for o, a in opts:
		if o == '-i' or o == '--input':
			inputbam = a
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

	G = baseinfo()
	G.get_config(ConfigFile)
	samtools = G.Samtools()
	genomesize = G.Genomesize()
	barcode_index = G.Barcode_index()

	C = CFCR()
	samdir = os.path.dirname(inputbam)
	sys.stderr.write("[ %s ] split sam file by chr, sort sam file by barcode ... \n" % time.asctime())
	Samfilepath = C.split_and_sort_sam(inputbam, samtools, barcode_index, samdir)
	sys.stderr.write("[ %s ] splited and sorted sam file list: %s \n\n" % (time.asctime(), Samfilepath))

	sys.stderr.write("[ %s ] calculating CF/CR ... \n" % time.asctime())
	stattxt = C.calculate(Samfilepath, genomesize, samdir)
	sys.stderr.write("[ %s ] CF/CR has been written to %s \n" % (time.asctime(), stattxt))
