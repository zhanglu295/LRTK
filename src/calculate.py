import os, sys, gzip
import getopt
import time
import subprocess
import re
from collections import defaultdict
import random
import string

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
	def find_samtools_version(self, samtools_path, tmpdir):
		randomstring = ''.join(random.sample(string.ascii_letters + string.digits, 8))
		tmpshell = tmpdir + "/" + randomstring + ".sh"
		tmplog = tmpdir + "/" + randomstring + ".log"
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

				sv = self.find_samtools_version(samtools_path, splitdir)
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

	def calculate(self, sorted_by_barcode_sam_list, target_size, outdir, sorted_bam, samtools_path, molecule_length):
		samfile_list = open(sorted_by_barcode_sam_list, 'r')
		statistics_report_dir = os.path.join(outdir, "statistics")
		if os.path.isdir(statistics_report_dir):
			pass
		else:
			os.makedirs(statistics_report_dir)
		marked_molecule = os.path.join(statistics_report_dir, "molecule.full.gz")
		wmarked_molecule = gzip.open(marked_molecule, 'wb')
		broken_molecule = os.path.join(statistics_report_dir, "molecule.broken.txt")
		wbroken_molecule = open(broken_molecule, 'w')

		mol_id = 1
		all_CF = 0
		all_CR = 0
		allbarcode = dict()
		readid_dict = defaultdict(dict)
		molecule = defaultdict(list)
		molecule2 = defaultdict(list)
		molecule3 = defaultdict(int)
		molecule4 = defaultdict(list)

		molecule_length_distribution = dict()
		barcode_molecule_amount = defaultdict(dict)
		molecule_really_covered = defaultdict(list)
		molecule_coverage_distribution = defaultdict(int)
		molecule_coverage_distribution_sum = 0

		barcode_molecule_amount_file = os.path.join(statistics_report_dir, "barcode_molecule_amount.txt")
		wbarcode_molecule_amount_file = open(barcode_molecule_amount_file, "w")
		wbarcode_molecule_amount_file.write("barcode\tamount_of_molecule\n")

		molecule_really_covered_file = os.path.join(statistics_report_dir, "each_molecule_coverage.txt")
		wmolecule_really_covered_file = open(molecule_really_covered_file, "w")
		wmolecule_really_covered_file.write("molecule_id\tmolecule_length\tcovered\tcoverage\n")

		molecule_coverage_distribution_file = os.path.join(statistics_report_dir, "molecule_coverage_distribution.txt")
		wmolecule_coverage_distribution_file = open(molecule_coverage_distribution_file, "w")
		wmolecule_coverage_distribution_file.write("coverage\tproportion(%)\n")

		filter_molecule_bam_list = list()
		for sam in samfile_list:
			sam = sam.strip()
			sys.stderr.write("[ %s ] processing: %s \n" %(time.asctime(), sam))
			samfile = open(sam, 'r')

			if molecule_length > 0:
				molecule_pass = 0
				filtered_samfile = (sam + ".filtered").replace("sam.filtered", "filtered.sam")
				tmpshell = filtered_samfile + ".sh"
				wtmpshell = open(tmpshell, 'w')
				shell_line = " ".join([samtools_path, "view -H", sorted_bam, "| grep -v coordinate >", filtered_samfile, "\n"])
				wtmpshell.write(shell_line)
				wtmpshell.close()
				subprocess.call(["sh", tmpshell])
				subprocess.call(["rm", tmpshell])

				wfiltered_samfile = open(filtered_samfile, 'a')

			while True:
				saminfo = samfile.readline().strip()
				saminfolist = re.split("\t", saminfo)

				if len(saminfo) == 0:
					break
				elif saminfolist[5] == "*":
					pass
				else:
					barcode_field = [s for s in saminfolist if "BX:Z:" in s]
					if barcode_field != []:
						bc = barcode_field[0].split(":")[2].split("-")[0]
						start_pos = int(saminfolist[3])
						readlength = len(saminfolist[9])
						end_pos = start_pos + readlength - 1
						chrid = saminfolist[2]
						if bc not in allbarcode:	
							if len(allbarcode) == 0:
								allbarcode[bc] = bc
								readid_dict[saminfolist[0]] = 1
								molecule[mol_id].append(start_pos)
								molecule2[mol_id].append(end_pos)
								molecule3[mol_id] = readlength
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
									barcodeinfo = sorted(allbarcode.keys())
									minpos = min(molecule[mol_id])
									maxpos = max(molecule2[mol_id])
									all_CR += molecule3[mol_id]
									small_cf = maxpos - minpos + 1
									all_CF += small_cf
									qn = len(molecule[mol_id])
									for eachinfo in molecule4[mol_id]:
										info_detail_list = re.split("\t", eachinfo)
										if small_cf not in molecule_length_distribution:
											molecule_length_distribution[small_cf] = 1
										else:
											molecule_length_distribution[small_cf] += 1
										outmoleculeinfo = "\t".join([chrid, str(minpos), str(maxpos), str(small_cf), str(qn), barcodeinfo[0], str(mol_id), info_detail_list[0], info_detail_list[1], str(len(info_detail_list[9]))]) + "\n"

										if barcodeinfo[0] not in barcode_molecule_amount:
											barcode_molecule_amount[barcodeinfo[0]][mol_id] = 1
										else:
											if mol_id not in barcode_molecule_amount[barcodeinfo[0]]:
												barcode_molecule_amount[barcodeinfo[0]][mol_id] = 1

										if mol_id not in molecule_really_covered:
											molecule_really_covered[mol_id].append(small_cf)
											molecule_really_covered[mol_id].append(len(info_detail_list[9]))
										else:
											molecule_really_covered[mol_id][1] += len(info_detail_list[9])

										if molecule_length > 0 and molecule_length < small_cf:
											molecule_pass = 1
											pass_filtered_info = eachinfo + "\n"
											wfiltered_samfile.write(pass_filtered_info)
										wmarked_molecule.write(outmoleculeinfo.encode())
									mol_id += 1
								else:
									for readid_ in readid_dict.keys():
										broken_info = "\t".join([str(readid_), str(readid_dict[readid_])]) + "\n"
										wbroken_molecule.write(broken_info)

								for barcode_molecule_barcode in barcode_molecule_amount.keys():
									separate_num = len(barcode_molecule_amount[barcode_molecule_barcode].keys())
									info = barcode_molecule_barcode + "\t" + str(separate_num) + "\n"
									wbarcode_molecule_amount_file.write(info)

								for mol_id_key in molecule_really_covered.keys():
									small_cfcf_rate = 1.0 * molecule_really_covered[mol_id_key][1] / molecule_really_covered[mol_id_key][0]
									info = str(mol_id_key) + "\t" + str(molecule_really_covered[mol_id_key][0]) + "\t" + str(molecule_really_covered[mol_id_key][1]) + "\t" + str(small_cfcf_rate) + "\n"
									wmolecule_really_covered_file.write(info)

									small_cfcf_rate_expand = int(100 * small_cfcf_rate)
									molecule_coverage_distribution[small_cfcf_rate_expand] += 1
									molecule_coverage_distribution_sum += 1

								molecule_really_covered = defaultdict(list)
								barcode_molecule_amount = defaultdict(dict)
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
										info_detail_list = re.split("\t", eachinfo)
										if small_cf not in molecule_length_distribution:
											molecule_length_distribution[small_cf] = 1
										else:
											molecule_length_distribution[small_cf] += 1
										outmoleculeinfo = "\t".join([chrid, str(minpos), str(maxpos), str(small_cf), str(qn), barcodeinfo[0], str(mol_id), info_detail_list[0], info_detail_list[1], str(len(info_detail_list[9]))]) + "\n"

										if barcodeinfo[0] not in barcode_molecule_amount:
											barcode_molecule_amount[barcodeinfo[0]][mol_id] = 1
										else:
											if mol_id not in barcode_molecule_amount[barcodeinfo[0]]:
												barcode_molecule_amount[barcodeinfo[0]][mol_id] = 1

										if mol_id not in molecule_really_covered:
											molecule_really_covered[mol_id].append(small_cf)
											molecule_really_covered[mol_id].append(len(info_detail_list[9]))
										else:
											molecule_really_covered[mol_id][1] += len(info_detail_list[9])

										if molecule_length > 0 and molecule_length < small_cf:
											molecule_pass = 1
											pass_filtered_info = eachinfo + "\n"
											wfiltered_samfile.write(pass_filtered_info)
										wmarked_molecule.write(outmoleculeinfo.encode())
									mol_id += 1
								else:
									for readid_ in readid_dict.keys():
										broken_info = "\t".join([str(readid_), str(readid_dict[readid_])]) + "\n"
										wbroken_molecule.write(broken_info)

								for barcode_molecule_barcode in barcode_molecule_amount.keys():
									separate_num = len(barcode_molecule_amount[barcode_molecule_barcode].keys())
									info = barcode_molecule_barcode + "\t" + str(separate_num) + "\n"
									wbarcode_molecule_amount_file.write(info)
								barcode_molecule_amount = defaultdict(dict)

								for mol_id_key in molecule_really_covered.keys():
									small_cfcf_rate = 1.0 * molecule_really_covered[mol_id_key][1] / molecule_really_covered[mol_id_key][0]
									info = str(mol_id_key) + "\t" + str(molecule_really_covered[mol_id_key][0]) + "\t" + str(molecule_really_covered[mol_id_key][1]) + "\t" + str(small_cfcf_rate) + "\n"
									wmolecule_really_covered_file.write(info)
									small_cfcf_rate_expand = int(100 * small_cfcf_rate)
									molecule_coverage_distribution[small_cfcf_rate_expand] += 1
									molecule_coverage_distribution_sum += 1

								molecule_really_covered = defaultdict(list)
								barcode_molecule_amount = defaultdict(dict)
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

			if molecule_length > 0 and molecule_pass > 0:
				wfiltered_samfile.close()
				tmpshell = filtered_samfile + ".sh"
				tmpdir = os.path.dirname(tmpshell)
				filtered_bamfile = filtered_samfile.replace("filtered.sam", "filtered.bam")
				wtmpshell = open(tmpshell, 'w')
				sv = self.find_samtools_version(samtools_path, tmpdir)
				if sv == 0:
					shell_line = " ".join([samtools_path, "view -h -S", filtered_samfile, "-b |", samtools_path, "sort -m 500M -", filtered_bamfile, "\nmv", filtered_bamfile + ".bam", filtered_bamfile, "\n"])
				else:
					shell_line = " ".join([samtools_path, "view -h -S", filtered_samfile, "-b |", samtools_path, "sort -m 500M -T", filtered_bamfile, "-o", filtered_bamfile,"-\n"])
				wtmpshell.write(shell_line)
				shell_line = " ".join([samtools_path, "index", filtered_bamfile, filtered_bamfile + ".bai\n"])
				wtmpshell.write(shell_line)
				wtmpshell.close()
				subprocess.call(["sh", tmpshell])
				subprocess.call(["rm", tmpshell])

				filter_molecule_bam_list.append(filtered_bamfile)

			samfile.close()
		wmarked_molecule.close()
		samfile_list.close()
		wbroken_molecule.close()

		filter_merge_bam = (sorted_bam + ".filtered").replace("bam.filtered", "filtered.bam")
		if molecule_length > 0 and len(filter_molecule_bam_list) > 0:
			shelldir = os.path.dirname(sorted_bam) + "/shell"
			if os.path.isdir(shelldir):
				pass
			else:
				os.makedirs(shelldir)
			tmpshell = shelldir + "/merge_molecule_filtered_bam.sh"
			wtmpshell = open(tmpshell, 'w')
			sv = self.find_samtools_version(samtools_path, shelldir)
			filtered_bam_name = " ".join(filter_molecule_bam_list)
			if sv == 0:
				if len(filter_molecule_bam_list) > 1:
					shell_line = " ".join([samtools_path, "merge -f", filter_merge_bam, filtered_bam_name, "\n", samtools_path, "index", filter_merge_bam, filter_merge_bam + ".bai\n"])
				else:
					shell_line = " ".join(["ln", "-s", "-f", filter_molecule_bam_list[0], filter_merge_bam + "\n", "ln", "-s", "-f", filter_molecule_bam_list[0] + ".bai", filter_merge_bam + ".bai\n"])
			else:
				shell_line = " ".join([samtools_path, "merge -f", filter_merge_bam, filtered_bam_name])
			wtmpshell.write(shell_line)
			wtmpshell.close()
			subprocess.call(["sh", tmpshell])
		else:
			subprocess.call(["ln", "-s", "-f", sorted_bam, filter_merge_bam])
			old_bai = sorted_bam.replace("marked.bai", "marked.bam.bai")
			new_bai = filter_merge_bam + ".bai"
			subprocess.call(["ln", "-s", "-f", old_bai, new_bai])

		cfcr_stat_file = os.path.join(statistics_report_dir, "CFCR.stat")
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

		molecule_length_distribution_file = os.path.join(statistics_report_dir, "molecule_length_distribution.txt")
		wmolecule_length_distribution_file = open(molecule_length_distribution_file, 'w')
		wmolecule_length_distribution_file.write("#molecule_length\tamount\n")
		for molecule_length_id in sorted(molecule_length_distribution.keys()):
			info = str(molecule_length_id) + "\t" + str(molecule_length_distribution[molecule_length_id]) + "\n"
			wmolecule_length_distribution_file.write(info)
		wmolecule_length_distribution_file.close()
		wbarcode_molecule_amount_file.close()
		wmolecule_really_covered_file.close()

		for n in range(0,max(molecule_coverage_distribution.keys()) + 1):
			rate = round(100.0 * molecule_coverage_distribution[n] / molecule_coverage_distribution_sum, 2)
			info = str(n * 0.01) + "\t" + str(rate) + "\n"
			wmolecule_coverage_distribution_file.write(info)
		wmolecule_coverage_distribution_file.close()

		return(statistics_report_dir)

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
	Molecule_length = 0
	opts, args = getopt.gnu_getopt(sys.argv[1:], 'i:o:m:c:h:', ['input', 'outputdir', 'molecule', 'config', 'help'])
	for o, a in opts:
		if o == '-i' or o == '--input':
			inputbam = a
		if o == '-o' or o == '--outputdir':
			outputdir = a
		if o == '-m' or o == '--molecule':
			Molecule_length = int(a)
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
	report_dir = C.calculate(Samfilepath, genomesize, samdir, inputbam, samtools, Molecule_length)
	sys.stderr.write("[ %s ] files of CF/CR, molecule distribution, molecule coverage and other infos have been placed in directory: %s \n" % (time.asctime(), report_dir))
