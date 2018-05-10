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

	def Sambamba(self):
		abs_path = self.get_path('sambamba')
		return abs_path
	
	def Samtools(self):
		abs_path = self.get_path('samtools')
		return abs_path

	def Barcode_index(self):
		abs_path = self.get_path('barcode_index')
		return abs_path

	def Genomesize(self):
		abs_path = self.get_path('genomesize')
		return abs_path

def check_info(result, attribute):
	if attribute == "file":
		if os.path.isfile(result) == False:
			sys.stderr.write("[ %s ] %s does not exist!\n" % (time.asctime(), result))
			sys.exit(-1)
	elif attribute == "dir":
		if os.path.isdir(result) == False:
			os.makedirs(result)

class CFCR:
	def split_and_sort_sam(self, sorted_bam, sambamba_path, outdir = './'):
		shelldir = os.path.join(outdir, "shell")
		check_info(shelldir, "dir")
		tmpshell = outdir + "/shell/check_sort.sh"
		wtmpshell = open(tmpshell, 'w')
		tmpheader = tmpshell + ".tmp.header"
		shell_line = " ".join([sambamba_path, "view -H", sorted_bam, "| grep queryname >", tmpheader, "\n"])
		wtmpshell.write(shell_line)
		wtmpshell.close()
		subprocess.call(["sh", tmpshell])
		subprocess.call(["rm", tmpshell])

		if os.path.getsize(tmpheader):
			sys.stderr.write("[ %s ] %s has been sorted by reads-name/barcode!\n" % (time.asctime(), sorted_bam))
			if sorted_bam.endswith(".bam"):
				sorted_sam = sorted_bam[0:(len(sorted_bam)-4)] + ".sam"
				bamtosam_shell = shelldir + "/bam2sam.sh"
				wbamtosam_shell = open(bamtosam_shell, "w")
				shell_line = " ".join([sambamba_path, "view -h", sorted_bam, ">", sorted_sam, "\n"])
				wbamtosam_shell.write(shell_line)
				wbamtosam_shell.close()

				subprocess.call(["sh", bamtosam_shell])
			return(sorted_sam)
		else:
			sys.stderr.write("[ %s ] %s has not been sorted by reads-name/barcode, and would be sorted before calculating CF/CR!\n" % (time.asctime(), sorted_bam))

		duplicated_bam = os.path.join(outdir, "duplicated.bam")
		unmap_bam = os.path.join(outdir, "unmapped.bam")
		additional_bam = os.path.join(outdir, "additional.bam")
		barcode_name_bam = os.path.join(outdir, "barcode_name.bam")
		stat_xls = outdir + "/mapping.xls"
		if os.path.exists(stat_xls) and os.path.getsize(stat_xls):
			sys.stderr.write("[ %s ] original bam file has been splited into 4 bam files: \n%s\n%s\n%s\n%s\n\n" % (time.asctime(), duplicated_bam, unmap_bam, additional_bam, barcode_name_bam))
		else:
			wstat_xls = open(stat_xls, 'w')
			total_reads = 0
			unmap_reads = 0
			duplicated_reads = 0
			uniq_map_reads = 0
			rsorted_sam = pysam.AlignmentFile(sorted_bam, 'rb')
			sys.stderr.write("[ %s ] split the original bam file into 4 bam files: \n%s\n%s\n%s\n%s\n\n" % (time.asctime(), duplicated_bam, unmap_bam, additional_bam, barcode_name_bam))
			sys.stderr.write("[ %s ] - Location\t\ttime elapsed(m)\t\tper 1M reads(m)\n" % time.asctime())
			stime = time.time()
			rtime = stime

			wduplicated_bam = pysam.AlignmentFile(duplicated_bam, 'wb', template = rsorted_sam)
			wunmap_bam = pysam.AlignmentFile(unmap_bam, 'wb', template = rsorted_sam)
			wadditional_bam = pysam.AlignmentFile(additional_bam, 'wb', template = rsorted_sam)
			wbarcode_name_bam = pysam.AlignmentFile(barcode_name_bam, 'wb', template = rsorted_sam)
	
			for saminfo in rsorted_sam:
				chr_real_name = "chr1"
				total_reads += 1
				if total_reads > 0 and total_reads % 1000000 == 0:
					etime = time.time()
					elapsed_time = round((etime - stime) / 60, 1)
					per_time = round((etime - rtime) / 60, 1)
					if saminfo.reference_id >= 0:
						epos = saminfo.reference_start
						echr = rsorted_sam.get_reference_name(saminfo.reference_id)
					else:
						epos = "*"
						echr = "*"
					sys.stderr.write("[ %s ] - %s:%s\t\t%s\t\t%s\n" % (time.asctime(), echr, str(epos), str(elapsed_time), str(per_time)))
					rtime = etime
				if saminfo.reference_id >= 0:
					chr_real_name = rsorted_sam.get_reference_name(saminfo.reference_id)
				if saminfo.flag > 1000:
					wduplicated_bam.write(saminfo)	## duplication reads
					duplicated_reads += 1
				elif len(saminfo.cigar) == 0:
					wunmap_bam.write(saminfo)	## unmmpaed reads
					unmap_reads += 1
				elif len(chr_real_name) > 5:
					wadditional_bam.write(saminfo)	## additional chromosomes that not include chr1-22, X, Y, M
				else:
					newinfo = pysam.AlignedSegment()
					newinfo = saminfo
					newtags = list()
					un = 0
					for eachtag in newinfo.tags:
						if eachtag[0] == "BC" or eachtag[0] == "BX":
							#eachtag[1] = eachtag[1].replace("-1", "")
							newinfo.query_name = eachtag[1] + "_" + newinfo.query_name
						if eachtag[0] == "SA":
							un = 1
					if un == 0:
						uniq_map_reads += 1
					wbarcode_name_bam.write(newinfo)
			wduplicated_bam.close()
			wunmap_bam.close()
			wadditional_bam.close()
			wbarcode_name_bam.close()
			sys.stderr.write("[ %s ] Original bam file has been splited into 4 files!\n" % time.asctime())

			mapping_rate = round(100.0 * (total_reads-unmap_reads) /  total_reads, 2)
			duplication_rate = round(100.0 * duplicated_reads / total_reads, 2)
			uniq_map_rate = round(100.0 * uniq_map_reads / total_reads, 2)
			stat_xls_info = "Mapping rate(%):\t" + str(mapping_rate) + "\n"
			wstat_xls.write(stat_xls_info)
			stat_xls_info = "Duplication rate(%):\t" + str(duplication_rate) + "\n"
			wstat_xls.write(stat_xls_info)
			stat_xls_info = "Unique mapping rate(%):\t" + str(uniq_map_rate) + "\n"
			wstat_xls.write(stat_xls_info)
			wstat_xls.close()

		sort_shell = outdir + "/shell/sort_barcode_bam.sh"
		sort_sign = sort_shell + ".sign"
		barcode_name_sorted_bam = os.path.join(outdir, "barcode_name.sorted.bam")
		barcode_name_sorted_sam = os.path.join(outdir, "barcode_name.sorted.sam")
		if os.path.exists(barcode_name_sorted_sam) and os.path.getsize(sort_sign):
			sys.stderr.write("[ %s ] sort bam file by barcode, output: %s\n\n" % (time.asctime(), barcode_name_sorted_bam))
		else:
			sys.stderr.write("[ %s ] sort bam file by barcode, processing ... \n\n\tsh %s \n\n" % (time.asctime(), sort_shell))
			wsort_shell = open(sort_shell, 'w')
			wsort_shell.write("set -e\n")
			shell_line = " ".join([sambamba, "sort -n -t 4 -m 8G --tmpdir=" + outdir, "-o", outdir + "/" + "barcode_name.sorted.bam", barcode_name_bam]) + "\n"
			wsort_shell.write(shell_line)
			shell_line = " ".join([sambamba, "view -h", barcode_name_sorted_bam, ">", barcode_name_sorted_sam, "\n"])
			wsort_shell.write(shell_line)
			shell_line = " ".join(["echo Done! >", sort_sign + "\n"])
			wsort_shell.write(shell_line)
			#shell_line = "rm " + barcode_name_bam + "\n"
			wsort_shell.write(shell_line)
			wsort_shell.close()
			subprocess.call(["sh", sort_shell])
			sys.stderr.write("[ %s ] sort bam file by barcode, done: %s\n\n" % (time.asctime(), barcode_name_sorted_sam))

		return(barcode_name_sorted_sam)

#	def barcode_calculate(self, molecule_id, rbam, wbam, wcsv, mol_length, baminfo_list):
	def barcode_calculate(self, molecule_id, wcsv, mol_length, baminfo_list):
		readid_dict = defaultdict(list)
		each_barcode_sequence = "A"
		for baminfo in baminfo_list:
			baminfolist = re.split("\t", baminfo.strip())
			each_barcode_sequence = (re.split("_", baminfolist[0]))[0]
			readid = (re.split("_", baminfolist[0]))[1]
#			print("original\t" + baminfolist[0])
			if readid not in readid_dict:
				readid_dict[readid].append(baminfolist[2])
				readid_dict[readid].append(baminfo)
			else:
				if baminfolist[2] != readid_dict[readid][0]:
					del readid_dict[readid]
				else:
					readid_dict[readid].append(baminfo)

		chrinfo_dict = defaultdict(dict)
		for readid in readid_dict.keys():
#			print("readid\t" + readid_dict[readid][0])
			chrid = readid_dict[readid][0]
			readlist1 = re.split("\t", readid_dict[readid][1])
			pos = int(readlist1[3])
			if pos not in chrinfo_dict[chrid]:
				chrinfo_dict[chrid][pos] = readid_dict[readid][1:]
			else:
				pz = pos - 200
				while(pz < pos+200):
					if pz not in chrinfo_dict[chrid]:
						chrinfo_dict[chrid][pz] = readid_dict[readid][1:]
						pz = pos + 500
					else:
						pz += 1

		CR_dict = dict()
		molecule_dict = defaultdict(dict)
		mi = 0
		for chrinfo in chrinfo_dict.keys():
			each_chr_dict = chrinfo_dict[chrinfo]
			sorted_dict = sorted(each_chr_dict.items(), key=lambda each_chr_dict:each_chr_dict[0], reverse = False)
			for pos in sorted_dict:
#				print(str(pos[0]) + "\t" + str(pos[1]))
#				sys.exit(-1)
				baminfolist = re.split("\t", pos[1][0])
				baminfo2list = list()
				readid = baminfolist[0]
				small_cr = len(baminfolist[9])
				start = int(baminfolist[3])
				end = start + small_cr
				if len(pos[1]) == 2:
					baminfo2list = re.split("\t", pos[1][1])
					small_cr += len(baminfo2list[9])
					start = min(start, int(baminfo2list[3]))
					end = max(end, int(baminfo2list[3]) + len(baminfo2list[9]))
				CR_dict[readid] = small_cr
				small_cf = end - start + 1
				if small_cf < 50000:
					if mi not in molecule_dict:
						molecule_dict[mi]["start"] = start
						molecule_dict[mi]["end"] = end
						molecule_dict[mi]["info"] = list()
						molecule_dict[mi]["info"].append(pos[1][0])
						if len(baminfo2list) != 0:
							molecule_dict[mi]["info"].append(pos[1][1])
					else:
						if (molecule_dict[mi]["end"] < start and (start - molecule_dict[mi]["end"] > 50000)) or (molecule_dict[mi]["start"] > end and (molecule_dict[mi]["start"] - end > 50000)):
							mi += 1
							molecule_dict[mi]["start"] = start
							molecule_dict[mi]["end"] = end
							molecule_dict[mi]["info"] = list()
							molecule_dict[mi]["info"].append(pos[1][0])
							if len(baminfo2list) != 0:
								molecule_dict[mi]["info"].append(pos[1][1])
						else:
							start = min(start, molecule_dict[mi]["start"])
							end = max(end, molecule_dict[mi]["end"])
						molecule_dict[mi]["start"] = start
						molecule_dict[mi]["end"] = end
#						molecule_dict[mi]["info"] = list()
						molecule_dict[mi]["info"].append(pos[1][0])
						if len(baminfo2list) != 0:
							molecule_dict[mi]["info"].append(pos[1][1])

		#each_barcode_sequence
		molecule_num_in_each_barcode = len(molecule_dict)
		'''
		molecule_num_in_each_barcode = 0
		for realmi in molecule_dict.keys():
			Molecule_len = molecule_dict[realmi]["end"] - molecule_dict[realmi]["start"] + 1
			if Molecule_len >= mol_length:
				molecule_num_in_each_barcode += 1
		'''
		clean_baminfo_list = list()

		if molecule_num_in_each_barcode >= 1:
			for realmi in molecule_dict.keys():
#				print(str(realmi) + "\t" + str(molecule_dict[realmi]["info"]))
				Molecule_len = molecule_dict[realmi]["end"] - molecule_dict[realmi]["start"] + 1
				if Molecule_len >= mol_length and len(molecule_dict[realmi]["info"]) > 1:
					molecule_id += 1
#					Chrid = rbam.get_reference_name(molecule_dict[realmi]["info"][0].reference_id)
					Chrid = (re.split("\t", molecule_dict[realmi]["info"][0]))[2]
					Start_pos = str(molecule_dict[realmi]["start"])
					End_pos = str(molecule_dict[realmi]["end"])
					Molecule_len = molecule_dict[realmi]["end"] - molecule_dict[realmi]["start"] + 1
					reads_num = len(molecule_dict[realmi]["info"])
					Total_length_reads = 0
					for i in range(0,len(molecule_dict[realmi]["info"])):
						baminfo = molecule_dict[realmi]["info"][i].strip()
						tmplist = re.split("_", baminfo)
						for i in range(1, len(tmplist)):
							if i == 1:
								baminfo = tmplist[i]
							else:
								baminfo = baminfo + "_" + tmplist[i]
						MI_tag = "MI:i:" + str(molecule_id)
						baminfo = baminfo + "\t" + MI_tag + "\n"
						clean_baminfo_list.append(baminfo)

						baminfo_list = re.split("\t", baminfo)
						if len(baminfo_list) < 9:
							print(baminfo)
						Total_length_reads += len(baminfo_list[9])
					Depth = round(1.0 * Total_length_reads / Molecule_len, 2)
					csvinfo = ",".join([str(molecule_id), each_barcode_sequence, Chrid, Start_pos, End_pos, str(Molecule_len), str(reads_num), str(Total_length_reads), str(Depth)]) + "\n"
					wcsv.write(csvinfo)

		return(molecule_id, clean_baminfo_list)

	def calculate(self, sorted_by_barcode_sam, outdir, sambamba_path, molecule_length, samtools_path):
		statistics_report_dir = outdir
		check_info(statistics_report_dir, "dir")

		merge_molecule_file = os.path.join(statistics_report_dir, "fragment_info.csv")
		wmerge_molecule_file = open(merge_molecule_file, 'w')
		wmerge_molecule_file.write('Molecule_id'+','+'Barcode_seq'+','+'Chr'+','+'Start_pos'+','+'End_pos'+','+'Molecule_len'+','+'# of reads'+','+'Total_length_reads'+','+'Depth'+'\n')

		rsorted_by_barcode_sam = open(sorted_by_barcode_sam, 'r')
		outsam = outdir + "/molecule_filtered.sam"
		woutsam = open(outsam, 'w')

		barcode_sequence = "A"
		mol_id = 0
		barcode_baminfo_list = list()
		p_total_reads = 0
		sys.stderr.write("[ %s ] - %s processing ... \n" % (time.asctime(), sorted_by_barcode_sam))
		sys.stderr.write("[ %s ] - Location\t\ttime elapsed(m)\t\tper 1M reads(m)\n" % time.asctime())
		stime = time.time()
		rtime = stime

		for saminfo in rsorted_by_barcode_sam:
			if saminfo.startswith("@"):
				woutsam.write(saminfo)
			else:
				p_total_reads += 1
				saminfolist = re.split("\t", saminfo.strip())
				new_barcode = (re.split("_", saminfolist[0]))[0]
				if p_total_reads > 0 and p_total_reads % 1000000 == 0:
					etime = time.time()
					elapsed_time = round((etime - stime) / 60, 1)
					per_time = round((etime - rtime) / 60, 1)
					sys.stderr.write("[ %s ] - %s\t\t%s\t\t%s\n" % (time.asctime(), new_barcode, str(elapsed_time), str(per_time)))
		#		if new_barcode == "AAACAAGTGCAGCACC":
		#			sys.exit(-1)
				if new_barcode != barcode_sequence:
					if barcode_sequence != "A":
						#mol_id = self.barcode_calculate(mol_id, rsorted_by_barcode_bam, woutbam, wmerge_molecule_file, molecule_length, barcode_baminfo_list)
						(mol_id, Clean_BamInfo_list) = self.barcode_calculate(mol_id, wmerge_molecule_file, molecule_length, barcode_baminfo_list)
						for Clean_BamInfo in Clean_BamInfo_list:
							woutsam.write(Clean_BamInfo)
					barcode_baminfo_list = list()
					barcode_sequence = new_barcode
					barcode_baminfo_list.append(saminfo)
				else:
					barcode_baminfo_list.append(saminfo)
			#mol_id = self.barcode_calculate(mol_id, rsorted_by_barcode_bam, woutbam, wmerge_molecule_file, molecule_length, barcode_baminfo_list)
		(mol_id, Clean_BamInfo_list) = self.barcode_calculate(mol_id, wmerge_molecule_file, molecule_length, barcode_baminfo_list)
		for Clean_BamInfo in Clean_BamInfo_list:
			woutsam.write(Clean_BamInfo)

		wmerge_molecule_file.close()
		rsorted_by_barcode_sam.close()
		woutsam.close()

		outsortbam = outsam.replace("molecule_filtered.sam", "molecule_filtered.sorted.bam")
		tmpshell = outdir + "/sort_final_bam.sh"
		wtmpshell = open(tmpshell, 'w')
		wtmpshell.write("set -e\n")
		shell_line = " ".join([samtools_path, "sort -m 5G --threads 4", "-o", outsortbam, outsam + "\n"])
		wtmpshell.write(shell_line)
		shell_line = " ".join([samtools_path, "index", outsortbam]) + "\n"
		wtmpshell.write(shell_line)
		wtmpshell.close()
		subprocess.call(["sh", tmpshell])
#		subprocess.call(["rm", tmpshell])
#		subprocess.call(["rm", outsam])

		return(outsortbam, merge_molecule_file)

	def calculate2(self, Filter_Sorted_Bam, csvfile, statistics_report_dir, all_genome_size, d_bed, sambamba_path, samtools_path):
		depth_shell = statistics_report_dir + "/depth.sh"
		depth_file = statistics_report_dir + "/all_depth.txt"
		wdepth_shell = open(depth_shell, 'w')
#		shell_line = " ".join([sambamba_path, "depth base -L", d_bed, Filter_Sorted_Bam, ">", depth_file]) + "\n"
		shell_line = " ".join([samtools_path, "depth -b", d_bed, Filter_Sorted_Bam, ">", depth_file]) + "\n"
		wdepth_shell.write(shell_line)
		wdepth_shell.close()
		subprocess.call(["sh", depth_shell])

		chr_size_dict = dict()
		chr_size_dict["All"] = 0
		rd_bed = open(d_bed, 'r')
		for bed in rd_bed:
#			bed_list = re.split("[:-]", bed.strip())
			bed_list = re.split("\t", bed.strip())
			chr_size_dict["All"] += int(bed_list[2]) - int(bed_list[1])
			if bed_list[0] not in chr_size_dict:
				chr_size_dict[bed_list[0]] = int(bed_list[2]) - int(bed_list[1])
			else:
				chr_size_dict[bed_list[0]] += int(bed_list[2]) - int(bed_list[1])
		rd_bed.close()

		depth_dict = defaultdict(dict)
		cov = (0, 1, 2, 4, 10, 20, 30, 50)
		for _chr_ in chr_size_dict.keys():
			for eachcov in cov:
				depth_dict[_chr_][eachcov] = 0

		rdepth_file = open(depth_file, 'r')
		for depthinfo in rdepth_file:
			depth_list = re.split("\t", depthinfo.strip())
			depth_list[2] = int(depth_list[2])
			for eachcov in cov:
				if depth_list[2] >= eachcov:
					if eachcov == 0:
						depth_dict["All"][eachcov] += depth_list[2]
						depth_dict[depth_list[0]][eachcov] += depth_list[2]
					else:
						depth_dict["All"][eachcov] += 1
						depth_dict[depth_list[0]][eachcov] += 1
		rdepth_file.close()
		subprocess.call(["rm", depth_file])

		summary_xls = statistics_report_dir + "/summary.xls"
		wsummary_xls = open(summary_xls, 'w')
		info = "\t".join(["Chr", "size(bp)", "average_depth", "coverage>=1X(%)", "coverage>=2X(%)", "coverage>=4X(%)", "coverage>=10X(%)", "coverage>=20X(%)", "coverage>=30X(%)", "coverage>=50X(%)",])
		wsummary_xls.write(info + "\n")
		for c in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, "X", "Y", "M", "MT", "All"):
			c = str(c)
			if c not in chr_size_dict:
				c = "chr" + c
			if c in chr_size_dict:
				info = c + "\t" + str(chr_size_dict[c])
				average_depth = round(1.0 * depth_dict[c][0] / chr_size_dict[c], 2)
				info = info + "\t" + str(average_depth)
				for i in range(1, len(cov)):
					coverage = round(100.0 * depth_dict[c][cov[i]] / chr_size_dict[c], 2)
					info = info + "\t" + str(coverage)
				wsummary_xls.write(info + "\n")
		wsummary_xls.write("\n\n\n######################################################################################\n")
#		wsummary_xls.close()

		NFP_summary = os.path.join(statistics_report_dir, "NFP_summary.xls")
#		W_MuFL_summary = os.path.join(statistics_report_dir, "weight_MuFL.xls")
		MuFL_summary = os.path.join(statistics_report_dir, "Unweight_MuFL.xls")
		rcsvfile = open(csvfile, 'r')
		wNFP_summary = open(NFP_summary, 'w')
		wMuFL_summary = open(MuFL_summary, 'w')
#		wW_MuFL_summary = open(W_MuFL_summary, 'w')
		all_CF = 0
		all_CR = 0
		#all_reads_len = 0
		NFP_dict = dict()
		molecule_len_dict = dict()
		for csvinfo in rcsvfile:
			csvinfolist = re.split(",", csvinfo.strip())
			if csvinfolist[0] != "Molecule_id":
				all_CF += int(csvinfolist[5])
				all_CR += int(csvinfolist[7])
				#all_reads_len += int()
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
		
		CF_depth = round(1.0 * all_CF / int(all_genome_size), 2)
		#CFCR_stat_info = "CF:\t" + str(all_CF) + "\t" + str(all_genome_size) + "\t" + str(CF_depth) + "\n"
		CFCR_stat_info = "CF:\t" + str(CF_depth) + "\n"
		wsummary_xls.write(CFCR_stat_info)
		CR_depth = round(1.0 * all_CR / all_CF)
		#CFCR_stat_info = "CR:\t" + str(all_CR) + "\t" + str(all_CF) + "\t" + str(CR_depth) + "\n"
		CFCR_stat_info = "CR:\t" + str(CR_depth) + "\n"
		wsummary_xls.write(CFCR_stat_info)
		C = round(1.0 * CF_depth * CR_depth, 2)
		CFCR_stat_info = "C:\t" + str(C) + "\n\n"
		wsummary_xls.write(CFCR_stat_info)

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
		wsummary_xls.write(CFCR_stat_info)

		total_molecule_length = 0
		total_molecule_num = 0
		wMuFL_summary.write("Molecule_len(kb)\tNums\n")
		sorted_molecule_len_dict = sorted(molecule_len_dict.items(), key=lambda molecule_len_dict:molecule_len_dict[0], reverse = False)
		for sortedmoleculelen in sorted_molecule_len_dict:
			info = str(sortedmoleculelen[0]) + "\t" + str(sortedmoleculelen[1]) + "\n"
			total_molecule_length += (int(sortedmoleculelen[0]) * int(sortedmoleculelen[1]))
			total_molecule_num += int(sortedmoleculelen[1])
			wMuFL_summary.write(info)
		wMuFL_summary.close()
#		wW_MuFL_summary.close()
		Un_MuFL = round((total_molecule_length / total_molecule_num), 2)
		CFCR_stat_info = "Unweighted MuFL(Kb):\t" + str(Un_MuFL) + "\n"
		wsummary_xls.write(CFCR_stat_info)

		molecule_length_N50 = 0
		for sortedmoleculelen in sorted_molecule_len_dict:
			ml = int(sortedmoleculelen[0]) * int(sortedmoleculelen[1])
			r1 = 1.0 * molecule_length_N50 / total_molecule_length
			r2 = 1.0 * (molecule_length_N50 + ml) / total_molecule_length
			if r1 < 0.5 and r2 >= 0.5:
				CFCR_stat_info = "Weighted MuFL(Kb):\t" + str(sortedmoleculelen[0]) + "\n"
				wsummary_xls.write(CFCR_stat_info)
			molecule_length_N50 += ml
			'''
			if molecule_length_N50 != 10000:
				for e in range(0, int(sortedmoleculelen[1])):
					r1 = molecule_length_N50 / total_molecule_length
					r2 = (molecule_length_N50 + int(sortedmoleculelen[0])) / total_molecule_length
					molecule_length_N50 += int(sortedmoleculelen[0])
					if r1 < 0.5 and r2 >= 0.5:
						CFCR_stat_info = "Weighted MuFL(Kb):\t" + str(sortedmoleculelen[0]) + "\n"
						wsummary_xls.write(CFCR_stat_info)
						molecule_length_N50 = 10000
			'''
#		wsummary_xls.close()

		mapping_xls = statistics_report_dir + "/mapping.xls"  # mapping.xls
		rmapping_xls = open(mapping_xls, 'r')
		for xlsinfo in rmapping_xls:
			wsummary_xls.write(xlsinfo)
		rmapping_xls.close()
		wsummary_xls.close()
#		subprocess.call(["cat", mapping_xls, ">>", summary_xls])
#		subprocess.call(["rm", mapping_xls])
		
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
	sambamba = G.Sambamba()
	genomesize = G.Genomesize()
	barcode_index = G.Barcode_index()
	samtools = G.Samtools()

	intervals_0 = os.path.dirname(ConfigFile) + "/original.intervals"

	bed_0 = outputdir + "/0based.bed"
	rintervals_0 = open(intervals_0, 'r')
	wbed_0 = open(bed_0, 'w')
	for inte in rintervals_0:
		inlist = re.split("[:-]", inte.strip())
		outin = "\t".join(inlist) + "\n"
		wbed_0.write(outin)
	rintervals_0.close()
	wbed_0.close()

	C = CFCR()
	samdir = outputdir
	#sys.stderr.write("[ %s ] sort bam file by barcode ... \n" % time.asctime())
	modified_bam = C.split_and_sort_sam(inputbam, sambamba, samdir)
	sys.stderr.write("[ %s ] sorted bam file: %s \n\n" % (time.asctime(), modified_bam))

	sys.stderr.write("[ %s ] calculating CF, CR, NFP and MuFL ... \n" % time.asctime())
	report_dir = outputdir
	stat_csv = report_dir + "/fragment_info.csv"
	filtered_bam = report_dir + "/molecule_filtered.sorted.bam"
	(filtered_bam, stat_csv) = C.calculate(modified_bam, samdir, sambamba, Molecule_length, samtools)    #stat_csv = report_dir + "/fragment_info.csv"
	sys.stderr.write("[ %s ] filtered bam file: %s \n" % (time.asctime(), filtered_bam))

##	sys.stderr.write("[ %s ] calculating depth ... \n" % (time.asctime()))
##	report_dir = C.calculate2(filtered_bam, stat_csv, report_dir, genomesize, bed_0, sambamba, samtools)
##	sys.stderr.write("[ %s ] statistical results have been generated in directory: %s \n" % (time.asctime(), report_dir))
#	subprocess.call(["python", PLOT_script, stat_csv, report_dir])
##	sys.stderr.write("[ %s ] figures have been plotted in %s !\n\n" % (time.asctime(), report_dir))
