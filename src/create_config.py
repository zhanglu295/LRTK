import os, sys
import getopt
import time
import re
import glob

class create_config:
	def make_intervals(self, original_region, size, outputdir, dataset_dir):
		dataset_dir = str(dataset_dir)
		fai = dataset_dir + '/GATK_bundle/Homo_sapiens_assembly38.fasta.fai'
		maxsize = dict()
		rfai = open(fai, 'r')
		for fai_info in rfai:
			fai_info_list = re.split("\t", fai_info.strip())
			maxsize[fai_info_list[0]] = int(fai_info_list[1])
		rfai.close()

		roriginal_region = open(original_region, 'r')
		based_bed = outputdir + "/original.intervals"
		wbased_bed = open(based_bed, 'w')
		out_intervals = outputdir + "/extend.intervals"
		wout_intervals = open(out_intervals, 'w')
		chr_id = None
		ori_min = None
		ori_max = 0
		ext_min = None
		ext_max = 0
		for bedinfo in roriginal_region:
			(chrid, start, end) = re.split("\t", bedinfo.strip())
			start = int(start)
			end = int(end)
			if len(chrid) < 6:
				bed_0 = chrid + ":" + str(start) + "-" + str(end) + "\n"
				wbased_bed.write(bed_0)
			if chr_id == None:
				chr_id = chrid
				ori_min = start
				ori_max = end
				if start - size < 0:
					ext_min = 1
				else:
					ext_min = start - size
				ext_max = end + size
			else:
				if chr_id != chrid:
					if maxsize[chr_id] <= ext_max:
						ext_max = maxsize[chr_id]
					intervals = str(chr_id) + ":" + str(ext_min) + "-" + str(ext_max) + "\n"
					wout_intervals.write(intervals)
					chr_id = chrid
					ori_min = start
					ori_max = end
					if start - size < 0:
						ext_min = 1
					else:
						ext_min = start - size
					ext_max = end + size
				else:
					if ori_min > start:
						sys.stderr.write("ERROR: the input bed file was not sorted!\n")
						sys.exit(-1)
					else:
						if start - size > ext_max:
							intervals = str(chr_id) + ":" + str(ext_min) + "-" + str(ext_max) + "\n"
							wout_intervals.write(intervals)
							ori_min = start
							ori_max = end
							ext_min = start - size
							ext_max = end + size
						else:
							ori_min = start
							ori_max = end
							ext_max = end + size
		if maxsize[chr_id] <= ext_max:
			ext_max = maxsize[chr_id]
		intervals = str(chr_id) + ":" + str(ext_min) + "-" + str(ext_max-size) + "\n"
		wout_intervals.write(intervals)
		wout_intervals.close()
		roriginal_region.close()
		wbased_bed.close()

		return(out_intervals)

	def find_dbsnp(self, dataset_dir):
		vcflist = glob.glob(dataset_dir + "/GATK_bundle/dbsnp_*.hg38.vcf*")
		vcf_version = 0
		vcf = None
		for vcffile in vcflist:
			if vcffile.endswith("idx") or vcffile.endswith("tbi"):
				pass
			else:
				vr = (re.split("\.", ((re.split("/",vcffile))[-1])))[0]
				vr = int(vr.replace("dbsnp_", ""))
				if vr > vcf_version:
					vcf_version = vr
					vcf = vcffile
		return(vcf)

	def Basic_config(self, Basic_config_file, script_dir, dataset_dir):
		wBasic_config_file = open(Basic_config_file, 'w')
		script_abs_path = os.path.abspath(sys.argv[0])
		bwa = script_dir + '/bwa'
		samtools = script_dir + '/samtools'
		sambamba = script_dir + '/sambamba'
		picard = script_dir + '/picard.jar'
		fastqc = script_dir + '/FastQC/fastqc'
		gatk = script_dir + '/GenomeAnalysisTK.jar'
		java = script_dir + '/java'
		config_line = "\n".join(["###basic software", "bwa = " + bwa, "sambamba = " + sambamba, "samtools = " + samtools, "picard = " + picard, "java = " + java, "", ""])
		wBasic_config_file.write(config_line)
		barcode_fa = dataset_dir + '/barcode_whitelist.fasta'
		config_line = "\n".join(["###barcode modification", "barcode_fa = " + barcode_fa, "barcode_aln_parameter = -l 5", "fastqc = " + fastqc, "", ""])
		wBasic_config_file.write(config_line)
		genome_fa = dataset_dir + '/GATK_bundle/Homo_sapiens_assembly38.fasta'
		config_line = "\n".join(["###alignment", "ref = " + genome_fa, "fq_aln_parameter = -B 16", "fq_sampe_parameter = -a 1000", "", ""])
		wBasic_config_file.write(config_line)
		config_line = "\n".join(["###merge multiple bam files belong to the same files", "picard_merge_parameter = ASSUME_SORTED=true USE_THREADING=true VALIDATION_STRINGENCY=LENIENT", "", ""])
		wBasic_config_file.write(config_line)
		config_line = "\n".join(["###mark duplication", "picard_mark_parameter = BARCODE_TAG=BX REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true", "", ""])
		wBasic_config_file.write(config_line)
		config_line = "\n".join(["###statistics", "barcode_index = -k12", "genomesize = 2923732648", "", ""])
		wBasic_config_file.write(config_line)
#		dbsnp_vcf = dataset_dir + '/GATK_bundle/dbsnp_146.hg38.vcf'
		dbsnp_vcf = self.find_dbsnp(dataset_dir)
		gold_indel_vcf = dataset_dir + '/GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
		hapmap = dataset_dir + '/GATK_bundle/hapmap_3.3.hg38.vcf.gz'
		omni = dataset_dir + '/GATK_bundle/1000G_omni2.5.hg38.vcf.gz'
		mills = dataset_dir + '/GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
		high_confidence = dataset_dir + '/GATK_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz'
		config_line = "\n".join(["###BQSR", "gatk = " + gatk, "dbsnp = " + dbsnp_vcf, "gold_indel = " + gold_indel_vcf, "hapmap = " + hapmap, "high_confidence = " + high_confidence, "omni = " + omni, "mills = " + mills, "intervals = " + extended_intervals, "", ""])
		wBasic_config_file.write(config_line)
		wBasic_config_file.close()

	def Reseq_config(self, Reseq_config_file, script_dir, dataset_dir):
		wReseq_config_file = open(Reseq_config_file, 'w')
		script_dir = str(script_dir)
		dataset_dir = str(dataset_dir)
		genome_fa = dataset_dir + '/GATK_bundle/Homo_sapiens_assembly38.fasta'
		gatk = script_dir + '/GenomeAnalysisTK.jar'
		sambamba = script_dir + '/sambamba'
#		dbsnp = dataset_dir + '/GATK_bundle/dbsnp_146.hg38.vcf'
		dbsnp_vcf = self.find_dbsnp(dataset_dir)
		hapmap = dataset_dir + '/GATK_bundle/hapmap_3.3.hg38.vcf.gz'
		omni = dataset_dir + '/GATK_bundle/1000G_omni2.5.hg38.vcf.gz'
		mills = dataset_dir + '/GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
		high_confidence = dataset_dir + '/GATK_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz'
		java = script_dir + '/java'
		config_line = "\n".join(["###basic software", "java = " + java, "", ""])
		wReseq_config_file.write(config_line)
		config_line = "\n".join(["###variation call", "gatk = " + gatk, "sambamba = " + sambamba, "ref = " + genome_fa, "dbsnp = " + dbsnp_vcf, "hapmap = " + hapmap, "high_confidence = " + high_confidence, "omni = " + omni, "mills = " + mills, "HaplotypeCaller_par = --variant_index_type LINEAR --variant_index_parameter 128000", "GenotypeGVCFs_par = None", "intervals = " + extended_intervals, "", ""])
		wReseq_config_file.write(config_line)
		config_line = "\n".join(["###SV call", "", ""])
		wReseq_config_file.write(config_line)
		Hapcut2 = script_dir + '/HAPCUT2'
		ExtractHAIRS = script_dir + '/extractHAIRS'
		fgbio = script_dir + '/fgbio-0.4.0-SNAPSHOT.jar'
		config_line = "\n".join(["###phasing", "extractHAIRS = " + ExtractHAIRS, "HAPCUT2 = " + Hapcut2, "fgbio = " + fgbio, "", ""])
		wReseq_config_file.write(config_line)
		wReseq_config_file.close()

		config_dir = os.path.dirname(Reseq_config_file)
		chrlistfile = os.path.join(config_dir, "chrlist.txt")
		wchrlistfile = open(chrlistfile, 'w')
		genomedictfile = dataset_dir + "/GATK_bundle/Homo_sapiens_assembly38.dict"
		rgenomedictfile = open(genomedictfile, 'r')
		for chrinfo in rgenomedictfile:
			chrinfo = chrinfo.strip()
			if chrinfo.startswith('@SQ'):
				chrinfolist = re.split("\t", chrinfo)
				chrid = chrinfolist[1].replace("SN:", "") + "\n"
				wchrlistfile.write(chrid)
		wchrlistfile.close()
		rgenomedictfile.close()
def usage():
	create_usage = \
	'''
	generate configuration file for LRTK-SEQ
	Version: 1.0.0
	Dependents: Python (>=3.0)
	Last Updated Date: 2017-11-14

	Usage: python create_config.py <command> [options]

	Command:
		Basic       generate configuration file for pre-processing
		Reseq       generate configuration file for resequencing
		all         generate  pre-processing and resequencing configuration files
	Basic options:
		-o --outputdir, output directory.
	Advanced options:
		-s --softwarepath, software directory [default: ./bin]
		-d --datasetpath, dataset directory [default: ./dataset]
		-h --help, print help info

	'''
	print(create_usage)

if __name__ == '__main__':
	if len(sys.argv) < 4:
		usage()
		sys.exit(-1)
	else:
		command = sys.argv[1]
		outputdir = None
		SoftwarePathDir = os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + "/bin"
		DatasetPathDir = os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + "/dataset"
		target_bed = DatasetPathDir + "/GATK_bundle/nonN.bed"
		extend_size = 500

		commandlist = ["Basic", "Reseq", "Denovo", "all"]
		if command not in commandlist:
			sys.stderr.write(" %s is not supported in command!\n\n" % command)
			usage()
			sys.exit(-1)

		if os.path.exists(target_bed) == False or os.path.isfile(target_bed) == False:
			reffa = DatasetPathDir + "/GATK_bundle/Homo_sapiens_assembly38.fasta"
			find_nonN_src = os.path.dirname(os.path.abspath(sys.argv[0])) + "/find_nonN_region.py"
			print("[ %s ] ERROR: %s does not exist, you can create it by using the command:\n" % (time.asctime(), target_bed))
			print("python %s %s %s\n" % (find_nonN_src, reffa, target_bed))
			sys.exit()

		opts, args = getopt.gnu_getopt(sys.argv[2:], 'o:s:d:b:e:h:', ['outputdir', 'softwarepath', 'datasetpath', 'bed', 'extend', 'help'])
		for o, a in opts:
			if o == '-o' or o == '--outputdir':
				outputdir = str(a)
			if o == '-s' or o == '--softwarepath':
				SoftwarePathDir = str(a)
			if o == '-d' or o == '--datasetpath':
				DatasetPathDir = str(a)
			if o == '-b' or o == '--bed':
				target_bed = str(a)
			if o == '-e' or o == '--extend':
				extend_size = int(a)
			if o == '-h' or o == '--help':
				usage()
				sys.exit(-1)

		C = create_config()
		extended_intervals = C.make_intervals(target_bed, extend_size, outputdir, DatasetPathDir)
		if command == "Basic" or command == "all":
			config_file = os.path.join(outputdir, "Basic.config")
			C.Basic_config(config_file, SoftwarePathDir, DatasetPathDir)
			sys.stderr.write("[ %s ] configuration file for `Basic` has been generated: %s \n" % (time.asctime(), config_file))

		if command == "Reseq" or command == "all":
			config_file = os.path.join(outputdir, "Reseq.config")
			C.Reseq_config(config_file, SoftwarePathDir, DatasetPathDir)
			sys.stderr.write("[ %s ] configuration file for `Reseq` has been generated: %s \n" % (time.asctime(), config_file))

#		if command == "Denovo" or command == "all":
#			config_file = os.path.join(outputdir, "Denovo.config")
#			C.Denovo_config(config_file, SoftwarePathDir, DatasetPathDir)
#			sys.stderr.write("[ %s ] configuration file for de novo assembly has been created: %s \n" % (time.asctime(), config_file))
