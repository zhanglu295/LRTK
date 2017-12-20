import os, sys, gzip
import getopt
import time
import subprocess
import re
from collections import defaultdict
import random
import string
import pysam

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

		tmpshell = os.path.join(splitdir, "check_sort.sh")
		wtmpshell = open(tmpshell, 'w')
		tmpheader = splitdir + "/tmp.header"
		shell_line = " ".join([samtools_path, "view -H", sorted_bam, "| grep coordinate >", tmpheader, "\n"])
		wtmpshell.write(shell_line)
		wtmpshell.close()
		subprocess.call(["sh", tmpshell])
#		subprocess.call([samtools_path, "view -H", sorted_bam, "| grep coordinate >", tmpheader])

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
		subprocess.call(["rm", tmpshell])

		wtmpshell = open(tmpshell, 'w')
		headersam = os.path.join(splitdir, 'header.sam')
		shell_line = " ".join([samtools_path, "view -H", sorted_bam, ">", headersam, "\n"])
		wtmpshell.write(shell_line)
		wtmpshell.close()
		subprocess.call(["sh", tmpshell])

		rsorted_sam = pysam.AlignmentFile(sorted_bam, 'rb')
		duplicated_bam = os.path.join(outdir, "duplicated.bam")
		unmap_bam = os.path.join(outdir, "unmapped.bam")
		additional_bam = os.path.join(outdir, "additional.bam")
		wduplicated_bam = pysam.AlignmentFile(duplicated_bam, 'wb', template = rsorted_sam)
		wunmap_bam = pysam.AlignmentFile(unmap_bam, 'wb', template = rsorted_sam)
		wadditional_bam = pysam.AlignmentFile(additional_bam, 'wb', template = rsorted_sam)

		chrdict = dict()
		splitsamlist = list()
		chrsam = os.path.join(splitdir, "chr1.sam")

		outputsam = pysam.AlignmentFile(chrsam, 'w', template = rsorted_sam)
		for saminfo in rsorted_sam:
			chr_real_name = "chr1"
			if saminfo.reference_id >= 0:
				chr_real_name = rsorted_sam.get_reference_name(saminfo.reference_id)
			if saminfo.flag > 1000:
				wduplicated_bam.write(saminfo)	## duplication reads
			elif chr_real_name == "*":
				wunmap_bam.write(saminfo)	## unmmpaed reads
			elif len(chr_real_name) > 5:
				wadditional_bam.write(saminfo)	## additional chromosomes that not include chr1-22, X, Y, M
			else:
				if chr_real_name not in chrdict:
					chrdict[chr_real_name] = chr_real_name
					outputsam.close()

					newchr = chr_real_name
					chrsam = os.path.join(splitdir, str(newchr) + '.sam')
					outputsam = pysam.AlignmentFile(chrsam, 'w', template = rsorted_sam)
					splitsamlist.append(chrsam)
					sys.stderr.write("[ %s ] split sorted sam file, processing %s: %s ...\n" % (time.asctime(), chr_real_name, chrsam))
					### move the barcode info to the 12th column
				newinfo = pysam.AlignedSegment()
				newinfo = saminfo
				newtags = list()
				BX_tag = ("BX", "1")
				newtags.append(BX_tag)
				for eachtag in newinfo.tags:
					if eachtag[0] == "BX":
						newtags[0] = eachtag
					else:
						newtags.append(eachtag)
				newinfo.tags = newtags
				outputsam.write(newinfo)
		outputsam.close()
		wduplicated_bam.close()
		wunmap_bam.close()
		wadditional_bam.close()
		
		samfilepath = os.path.join(splitdir, "sam.txt")
		wsamfilepath = open(samfilepath, 'w')
		for chrsam in splitsamlist:
			chrsortedsam = chrsam.replace('sam', 'sorted_by_barcode.sam')
			sys.stderr.write("[ %s ] sort sam file, processing %s ...\n" % (time.asctime(), chrsortedsam))
			tmpshell = os.path.join(splitdir, "sort.sh")
			wtmpshell = open(tmpshell, "w")
			shell_line = " ".join(["awk '$1!~/^@/'", chrsam, "|sort", "-k12,12 -k4n", ">", chrsortedsam, "\n"])
			wtmpshell.write(shell_line)
			shell_line = " ".join(["rm", chrsam, "\n"])
			wtmpshell.write(shell_line)
			wtmpshell.close()
			subprocess.call(["sh", tmpshell])
			wsamfilepath.write("\n".join([chrsortedsam, ""]))
		wsamfilepath.close()

#		subprocess.call(["rm", sorted_sam])
		return(samfilepath, headersam)

	def calculate(self, sorted_by_barcode_sam_list, target_size, outdir, sorted_bam, samtools_path, molecule_length, headersamfile):
		samfile_list = open(sorted_by_barcode_sam_list, 'r')
		statistics_report_dir = os.path.join(outdir, "statistics")
		if os.path.isdir(statistics_report_dir):
			pass
		else:
			os.makedirs(statistics_report_dir)
		marked_molecule = os.path.join(statistics_report_dir, "molecule.full.gz")
		wmarked_molecule = gzip.open(marked_molecule, 'wb')
#		broken_molecule = os.path.join(statistics_report_dir, "molecule.broken.txt")
#		wbroken_molecule = open(broken_molecule, 'w')

		mol_id = 1
		all_CF = 0
		all_CR = 0
		allbarcode = dict()
		readid_dict = defaultdict(int)
		molecule = defaultdict(list)
		molecule2 = defaultdict(list)
		molecule3 = defaultdict(int)
		molecule4 = defaultdict(list)

		molecule_length_distribution = dict()
		barcode_molecule_amount = defaultdict(dict)
		molecule_really_covered = defaultdict(list)
		molecule_coverage_distribution = defaultdict(int)
		molecule_coverage_distribution_sum = 0

#		barcode_molecule_amount_file = os.path.join(statistics_report_dir, "barcode_molecule_amount.txt")
#		wbarcode_molecule_amount_file = open(barcode_molecule_amount_file, "w")
#		wbarcode_molecule_amount_file.write("barcode\tamount_of_molecule\n")

		molecule_really_covered_file = os.path.join(statistics_report_dir, "each_molecule_coverage.txt")
		wmolecule_really_covered_file = open(molecule_really_covered_file, "w")
		wmolecule_really_covered_file.write("molecule_id\tmolecule_length\tcovered_bases\tCR\n")

		molecule_coverage_distribution_file = os.path.join(statistics_report_dir, "molecule_coverage_distribution.txt")
		wmolecule_coverage_distribution_file = open(molecule_coverage_distribution_file, "w")
		wmolecule_coverage_distribution_file.write("Coverage\tProportion(%)\n")

		filter_molecule_bam_list = list()
		failed_molecule_bam_list = list()
		for sam in samfile_list:
			sam = sam.strip()
			sys.stderr.write("[ %s ] processing: %s \n" %(time.asctime(), sam))
			samfile = open(sam, 'r')

			tmpsaminfodict = dict()
			if molecule_length > 0:
				molecule_pass = 0
				filtered_samfile = sam.replace("sorted_by_barcode.sam", "sorted_by_barcode.filtered.sam")
				failed_samfile = sam.replace("sorted_by_barcode.sam", "sorted_by_barcode.failed.sam")
				tmpshell = filtered_samfile + ".sh"
				wtmpshell = open(tmpshell, 'w')
				shell_line = " ".join([samtools_path, "view -H", sorted_bam, "| grep -v coordinate >", filtered_samfile, "\n"])
				wtmpshell.write(shell_line)
				shell_line = " ".join([samtools_path, "view -H", sorted_bam, "| grep -v coordinate >", failed_samfile, "\n"])
				wtmpshell.write(shell_line)
				wtmpshell.close()
				subprocess.call(["sh", tmpshell])
				subprocess.call(["rm", tmpshell])

				wfiltered_samfile = open(filtered_samfile, 'a')
				wfailed_samfile = open(failed_samfile, 'a')

			while True:
				saminfo = samfile.readline().strip()
				saminfolist = re.split("\t", saminfo)

				if len(saminfo) == 0:
					break
				elif saminfolist[5] == "*":
					saminfo = saminfo + "\n"
					wfailed_samfile.write(saminfo)
				else:
					if saminfolist[0] in tmpsaminfodict:
						tmpsaminfodict[saminfolist[0]] += "\n" + saminfo
					else:
						tmpsaminfodict[saminfolist[0]] = saminfo
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
											pass_filtered_info_list = info_detail_list[0:11]
											for tag_info in range(11, len(info_detail_list)):
												if info_detail_list[tag_info].startswith("BI:Z:") or info_detail_list[tag_info].startswith("BD:Z:"):
													pass
												else:
													pass_filtered_info_list.append(info_detail_list[tag_info])
											pass_filtered_info = "\t".join(pass_filtered_info_list) + "\n"
											wfiltered_samfile.write(pass_filtered_info)
										wmarked_molecule.write(outmoleculeinfo.encode())
									mol_id += 1
									tmpsaminfodict = dict()
								else:
									for readid_ in tmpsaminfodict.keys():
										broken_info = tmpsaminfodict[readid_] + "\n"
										wfailed_samfile.write(broken_info)
									tmpsaminfodict = dict()

#								for barcode_molecule_barcode in barcode_molecule_amount.keys():
#									separate_num = len(barcode_molecule_amount[barcode_molecule_barcode].keys())
#									info = barcode_molecule_barcode + "\t" + str(separate_num) + "\n"
#									wbarcode_molecule_amount_file.write(info)

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
											pass_filtered_info_list = info_detail_list[0:11]
											for tag_info in range(11, len(info_detail_list)):
												if info_detail_list[tag_info].startswith("BI:Z:") or info_detail_list[tag_info].startswith("BD:Z:"):
													pass
												else:
													pass_filtered_info_list.append(info_detail_list[tag_info])
											pass_filtered_info = "\t".join(pass_filtered_info_list) + "\n"
											wfiltered_samfile.write(pass_filtered_info)
										wmarked_molecule.write(outmoleculeinfo.encode())
									mol_id += 1
									tmpsaminfodict = dict()
								else:
									for readid_ in tmpsaminfodict.keys():
										broken_info = tmpsaminfodict[readid_] + "\n"
										wfailed_samfile.write(broken_info)
									tmpsaminfodict = dict()

#								for barcode_molecule_barcode in barcode_molecule_amount.keys():
#									separate_num = len(barcode_molecule_amount[barcode_molecule_barcode].keys())
#									info = barcode_molecule_barcode + "\t" + str(separate_num) + "\n"
#									wbarcode_molecule_amount_file.write(info)
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

			if molecule_length > 0:
				wfiltered_samfile.close()
				wfailed_samfile.close()
				tmpshell = filtered_samfile + ".sh"
				tmpdir = os.path.dirname(tmpshell)
				filtered_bamfile = filtered_samfile.replace("filtered.sam", "filtered.bam")
				failed_bamfile = failed_samfile.replace("failed.sam", "failed.bam")
				wtmpshell = open(tmpshell, 'w')
				sv = self.find_samtools_version(samtools_path, tmpdir)
				if sv == 0:
					shell_line = " ".join([samtools_path, "view -h -S", filtered_samfile, "-b >", filtered_samfile + ".bam\n", samtools_path, "sort -m 500M", filtered_samfile + ".bam", filtered_bamfile, "\nmv", filtered_bamfile + ".bam", filtered_bamfile, "\n"])
					wtmpshell.write(shell_line)
					shell_line = " ".join([samtools_path, "view -h -S", failed_samfile, "-b >", failed_samfile + ".bam\n",  samtools_path, "sort -m 500M", failed_samfile + ".bam", failed_bamfile, "\nmv", failed_bamfile + ".bam", failed_bamfile, "\n"])
					wtmpshell.write(shell_line)
				else:
					shell_line = " ".join([samtools_path, "view -h -S", filtered_samfile, "-b |", samtools_path, "sort -m 500M -T - -o", filtered_bamfile, "\n"])
					wtmpshell.write(shell_line)
					shell_line = " ".join([samtools_path, "view -h -S", failed_samfile, "-b |", samtools_path, "sort -m 500M -T - -o", failed_bamfile, "\n"])
					wtmpshell.write(shell_line)
				shell_line = " ".join([samtools_path, "index", filtered_bamfile, filtered_bamfile + ".bai\n"])
				wtmpshell.write(shell_line)
				shell_line = " ".join([samtools_path, "index", failed_bamfile, failed_bamfile + ".bai\n"])
				wtmpshell.write(shell_line)
				wtmpshell.close()
				subprocess.call(["sh", tmpshell])
				subprocess.call(["rm", tmpshell])

				filter_molecule_bam_list.append(filtered_bamfile)
				failed_molecule_bam_list.append(failed_bamfile)
			samfile.close()
		wmarked_molecule.close()
		samfile_list.close()
#		wbroken_molecule.close()

		filter_merge_bam = (sorted_bam + ".filtered").replace("bam.filtered", "filtered.bam")
		failed_merge_bam = filter_merge_bam.replace("filtered.bam", "failed.bam")
		if molecule_length > 0 and len(filter_molecule_bam_list) > 0:
			shelldir = os.path.dirname(sorted_bam) + "/shell"
			if os.path.isdir(shelldir):
				pass
			else:
				os.makedirs(shelldir)
			tmpshell = shelldir + "/merge_molecule_filtered_bam.sh"
			wtmpshell = open(tmpshell, 'w')
			if len(filter_molecule_bam_list) > 1:
				filtered_bam_name = " ".join(filter_molecule_bam_list)
				failed_bam_name = " ".join(failed_molecule_bam_list)
				shell_line = " ".join([samtools_path, "merge -f", filter_merge_bam, filtered_bam_name, "\n", samtools_path, "index", filter_merge_bam, filter_merge_bam + ".bai\n"])
				wtmpshell.write(shell_line)
				shell_line = " ".join([samtools_path, "merge -f", failed_merge_bam, failed_bam_name, "\n", samtools_path, "index", failed_merge_bam, failed_merge_bam + ".bai\n"])
				wtmpshell.write(shell_line)
			else:
				shell_line = " ".join(["ln", "-s", "-f", filter_molecule_bam_list[0], filter_merge_bam + "\n", "ln", "-s", "-f", filter_molecule_bam_list[0] + ".bai", filter_merge_bam + ".bai\n"])
				wtmpshell.write(shell_line)
				shell_line = " ".join(["ln", "-s", "-f", failed_molecule_bam_list[0], failed_bam_name + "\n", "ln", "-s", "-f", failed_molecule_bam_list[0] + ".bai", failed_merge_bam + ".bai\n"])
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
		stat_line = "Item" + "\t" + "Total_length_Molecule(CF)/Short_Reads(CR)" + "\t" + "genome_size(CR)/total_length_Molecule(CF)" + "\t" + "CF/CR\n"
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
		wmolecule_length_distribution_file.write("Molecule_length\tNumber\n")
		for molecule_length_id in sorted(molecule_length_distribution.keys()):
			info = str(molecule_length_id) + "\t" + str(molecule_length_distribution[molecule_length_id]) + "\n"
			wmolecule_length_distribution_file.write(info)
		wmolecule_length_distribution_file.close()
#		wbarcode_molecule_amount_file.close()
		wmolecule_really_covered_file.close()

		for n in range(0,max(molecule_coverage_distribution.keys()) + 1):
			rate = round(100.0 * molecule_coverage_distribution[n] / molecule_coverage_distribution_sum, 2)
			info = str(n * 0.01) + "\t" + str(rate) + "\n"
			wmolecule_coverage_distribution_file.write(info)
		wmolecule_coverage_distribution_file.close()

		rmarked_molecule = gzip.open(marked_molecule, 'rb')
		merge_molecule_file = os.path.join(statistics_report_dir, "fragment_info.csv")
		wmerge_molecule_file = open(merge_molecule_file, 'w')
		wmerge_molecule_file.write('Molecule_id'+','+'Baecode_seq'+','+'Chr'+','+'Start_pos'+','+'End_pos'+','+'Molecule_len'+','+'# of reads'+','+'Total_length_reads'+','+'Depth'+'\n')
		molecule_id = 0
		merge_molecule_info_list = list()
		for mm in range(9):
			merge_molecule_info_list.append(0)
		for molecule_read_info in rmarked_molecule:
			molecule_read_info_list = re.split("\t", molecule_read_info.decode().strip())
			if int(molecule_read_info_list[6]) == molecule_id:
				merge_molecule_info_list[7] += int(molecule_read_info_list[9])
			else:
				if molecule_id > 0:
					merge_molecule_info_list[8] = str(round(1.0 * merge_molecule_info_list[7] / int(merge_molecule_info_list[5]), 2))
					merge_molecule_info_list[7] = str(merge_molecule_info_list[7])
					merge_molecule_info = ",".join(merge_molecule_info_list) + "\n"
					wmerge_molecule_file.write(merge_molecule_info)

					merge_molecule_info_list[0] = str(molecule_read_info_list[6])
					merge_molecule_info_list[1] = str(molecule_read_info_list[5])
					merge_molecule_info_list[2] = str(molecule_read_info_list[0])
					merge_molecule_info_list[3] = str(molecule_read_info_list[1])
					merge_molecule_info_list[4] = str(molecule_read_info_list[2])
					merge_molecule_info_list[5] = str(molecule_read_info_list[3])
					merge_molecule_info_list[6] = str(molecule_read_info_list[4])
					merge_molecule_info_list[7] = int(molecule_read_info_list[9])
					molecule_id = int(molecule_read_info_list[6])
				else:
					merge_molecule_info_list[0] = str(molecule_read_info_list[6])
					merge_molecule_info_list[1] = str(molecule_read_info_list[5])
					merge_molecule_info_list[2] = str(molecule_read_info_list[0])
					merge_molecule_info_list[3] = str(molecule_read_info_list[1])
					merge_molecule_info_list[4] = str(molecule_read_info_list[2])
					merge_molecule_info_list[5] = str(molecule_read_info_list[3])
					merge_molecule_info_list[6] = str(molecule_read_info_list[4])
					merge_molecule_info_list[7] = int(molecule_read_info_list[9])
					molecule_id = int(molecule_read_info_list[6])
		merge_molecule_info_list[8] = str(round(1.0 * merge_molecule_info_list[7] / int(merge_molecule_info_list[5]), 2))
		merge_molecule_info_list[7] = str(merge_molecule_info_list[7])
		merge_molecule_info = ",".join(merge_molecule_info_list) + "\n"
		wmerge_molecule_file.write(merge_molecule_info)
		wmerge_molecule_file.close()
		rmarked_molecule.close()

		return(statistics_report_dir)

	def calculate2(self, csvfile, statistics_report_dir):
		CFCR_stat = os.path.join(statistics_report_dir, "CFCR.stat")
		NFP_summary = os.path.join(statistics_report_dir, "NFP_summary.xls")
		MuFL_summary = os.path.join(statistics_report_dir, "Unweight_MuFL.xls")
		rcsvfile = open(csvfile, 'r')
		wCFCR_stat = open(CFCR_stat, 'a')
		wNFP_summary = open(NFP_summary, 'w')
		wMuFL_summary = open(MuFL_summary, 'w')
		wW_MuFL_summary = open(W_MuFL_summary, 'w')
		NFP_dict = dict()
		molecule_len_dict = dict()
		for csvinfo in rcsvfile:
			csvinfolist = re.split(",", csvinfo.strip())
			if csvinfolist[0] != "Molecule_id":
				if csvinfolist[1] not in NFP_dict:
					NFP_dict[csvinfolist[1]] = 1
				else:
					NFP_dict[csvinfolist[1]] += 1

				mst = int(int(csvinfolist[5]) / 1000)
				if mst not in molecule_len_dict:
					molecule_len_dict[mst] = 1
				else:
					molecule_len_dict[mst] +=1
		rcsvfile.close()

		molecule_count = 0
		barcode_count = 0
		sorted_NFP_dict = sorted(NFP_dict.items(), key=lambda NFP_dict:NFP_dict[1], reverse = True)
		for sortedNFP in sorted_NFP_dict:
			info = str(sortedNFP[0]) + "\t" + str(sortedNFP[1]) + "\n"
			wNFP_summary.write(info)
			barcode_count += 1
			molecule_count += int(sortedNFP[1])
		wNFP_summary.close()

		barcode_average_molecule = round(1.0 * molecule_count / barcode_count, 2)
		CFCR_stat_info = "\nNFP:\t" + str(barcode_average_molecule) + "\n"
		wCFCR_stat.write(CFCR_stat_info)

		total_molecule_length = 0
		sorted_molecule_len_dict = sorted(molecule_len_dict.items(), key=lambda molecule_len_dict:molecule_len_dict[0], reverse = False)
		for sortedmoleculelen in sorted_molecule_len_dict:
			info = str(sortedmoleculelen[0]) + "\t" + str(sortedmoleculelen[1]) + "\n"
			total_molecule_length += int(sortedmoleculelen[0]) * int(sortedmoleculelen[1])
			wMuFL_summary.write(info)
		wMuFL_summary.close()
		wW_MuFL_summary.close()

		molecule_length_N50 = 0
		for sortedmoleculelen in sorted_molecule_len_dict:
			if molecule_length_N50 != 1:
				for e in range(0, int(sortedmoleculelen[1])):
					r1 = molecule_length_N50 / total_molecule_length
					r2 = (molecule_length_N50 + int(sortedmoleculelen[0])) / total_molecule_length
					molecule_length_N50 += int(sortedmoleculelen[0])
					if r1 < 0.5 and r2 >= 0.5:
						CFCR_stat_info = "Weighted MuFL:\t" + str(sortedmoleculelen[0]) + "\n"
						wCFCR_stat.write(CFCR_stat_info)
						molecule_length_N50 = 1
		wCFCR_stat.close()
		
		return(statistics_report_dir)

def usage():
	calculation_usage = \
	'''
	calculate linked read statistics
	Version: 1.0.0
	Dependents: Python (>=3.0), SAMtools
	Last Updated Date: 2017-11-15

	Usage: python calculate.py [options]

	Basic Options:
		-i --input, path of input bam
		-o --outputdir, the path of output
		-c --config, the path of configuration file [default: ./config/Basic.config]
		-h --help, print help info
	Advanced Options:
		-m --minlen, minimum molecule length [default: 500bp]

	'''
	print(calculation_usage)

if __name__ == '__main__':
	if len(sys.argv) < 5:
		usage()
		sys.exit(-1)

	inputbam = None
	outputdir = None
	ConfigFile = None
	Molecule_length = 500
	opts, args = getopt.gnu_getopt(sys.argv[1:], 'i:o:m:c:h:', ['input', 'outputdir', 'minlen', 'config', 'help'])
	for o, a in opts:
		if o == '-i' or o == '--input':
			inputbam = a
		if o == '-o' or o == '--outputdir':
			outputdir = a
		if o == '-m' or o == '--minlen':
			Molecule_length = int(a)
		if o == '-c' or o == '--config':
			ConfigFile = a
		if o == '-h' or o == '--help':
			usage()
			sys.exit(-1)

	if ConfigFile == None or os.path.exists(ConfigFile) == False:
		sys.stderr.write("configuration file has not been provided or does not exist, Please create it using 'python LRTK-SEQ.py Config'\n")
		sys.exit(-1)
	
	PLOT_script = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "plot.py")
	if os.path.exists(PLOT_script) == False:
		sys.stderr.write("%s does not exist, the software package might not been downloaded perfectly!" % PLOT_script)
		sys.exit(-1)

	G = baseinfo()
	G.get_config(ConfigFile)
	samtools = G.Samtools()
	genomesize = G.Genomesize()
	barcode_index = G.Barcode_index()

	C = CFCR()
	samdir = os.path.dirname(inputbam)
	sys.stderr.write("[ %s ] split sam file by chromosome, sort sam file by barcode ... \n" % time.asctime())
	(Samfilepath, HeaderSamFile) = C.split_and_sort_sam(inputbam, samtools, barcode_index, samdir)
	sys.stderr.write("[ %s ] splited and sorted sam file list: %s \n\n" % (time.asctime(), Samfilepath))

	sys.stderr.write("[ %s ] calculating CF, CR, NFP and MuFL ... \n" % time.asctime())
	report_dir = samdir + "/statistics"
	report_dir = C.calculate(Samfilepath, genomesize, samdir, inputbam, samtools, Molecule_length, HeaderSamFile)
	sys.stderr.write("[ %s ] statistical results have been generated in directory: %s \n" % (time.asctime(), report_dir))

	stat_csv = report_dir + "/fragment_info.csv"
	report_dir = C.calculate2(stat_csv, report_dir)
	subprocess.call(["python", PLOT_script, stat_csv, report_dir])
	sys.stderr.write("[ %s ] figures have been plotted in %s !\n\n" % (time.asctime(), report_dir))
