import os, sys
import getopt
import time
import re

class create_config:
	def Basic_config(self, Basic_config_file, script_dir, dataset_dir):
		wBasic_config_file = open(Basic_config_file, 'w')
		script_abs_path = os.path.abspath(sys.argv[0])
		bwa = script_dir + '/bwa'
		samtools = script_dir + '/samtools'
		picard = script_dir + '/picard.jar'
		fastqc = script_dir + '/FastQC/fastqc'
		java = 'java'
		config_line = "\n".join(["###basic software", "bwa = " + bwa, "samtools = " + samtools, "picard = " + picard, "java = " + java, "", ""])
		wBasic_config_file.write(config_line)
		barcode_fa = dataset_dir + '/barcode/barcode.fa'
		config_line = "\n".join(["###barcode modification", "barcode_fa = " + barcode_fa, "barcode_aln_parameter = -l 5", "fastqc = " + fastqc, "", ""])
		wBasic_config_file.write(config_line)
		genome_fa = dataset_dir + '/genome/genome.fa'
		config_line = "\n".join(["###alignment", "ref = " + genome_fa, "fq_aln_parameter = -B 16", "fq_sampe_parameter = -a 1000", "", ""])
		wBasic_config_file.write(config_line)
		config_line = "\n".join(["###merge multiple bam files belong to the same files", "picard_merge_parameter = ASSUME_SORTED=true USE_THREADING=true VALIDATION_STRINGENCY=LENIENT", "", ""])
		wBasic_config_file.write(config_line)
		config_line = "\n".join(["###mark duplication", "picard_mark_parameter = BARCODE_TAG=BX REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true", "", ""])
		wBasic_config_file.write(config_line)
		config_line = "\n".join(["###statistics", "barcode_index = -k12", "genomesize = 2861343787", "", ""])
		wBasic_config_file.write(config_line)
		wBasic_config_file.close()

	def Reseq_config(self, Reseq_config_file, script_dir, dataset_dir):
		wReseq_config_file = open(Reseq_config_file, 'w')
		script_dir = str(script_dir)
		dataset_dir = str(dataset_dir)
		genome_fa = dataset_dir + '/genome/genome.fa'
		gatk = script_dir + '/gatk-package-4.beta.1-local.jar'
		dbsnp = dataset_dir + '/dbsnp137.vcf'
		java = 'java'
		config_line = "\n".join(["###basic software", "java = " + java, "", ""])
		wReseq_config_file.write(config_line)
		config_line = "\n".join(["###variation call", "gatk = " + gatk, "ref = " + genome_fa, "dbsnp = " + dbsnp, "HaplotypeCaller_par = --variant_index_type LINEAR --variant_index_parameter 128000", "GenotypeGVCFs_par = None", "", ""])
		wReseq_config_file.write(config_line)
		config_line = "\n".join(["###SV call", "", ""])
		wReseq_config_file.write(config_line)
		Hapcut2 = script_dir + '/HAPCUT2'
		ExtractHAIRS = script_dir + '/extractHAIRS'
		fgbio = script_dir + '/fabio-0.2.0-SNAPSHOT.jar'
		config_line = "\n".join(["###phasing", "extractHAIRS = " + ExtractHAIRS, "HAPCUT2 = " + Hapcut2, "fgbio = " + fgbio, "", ""])
		wReseq_config_file.write(config_line)
		wReseq_config_file.close()

		config_dir = os.path.dirname(Reseq_config_file)
		chrlistfile = os.path.join(config_dir, "chrlist.txt")
		wchrlistfile = open(chrlistfile, 'w')
		genomedictfile = dataset_dir + "/genome/genome.dict"
		rgenomedictfile = open(genomedictfile, 'r')
		for chrinfo in rgenomedictfile:
			chrinfo = chrinfo.strip()
			if chrinfo.startswith('@SQ'):
				chrinfolist = re.split("\t", chrinfo)
				chrid = chrinfolist[1].replace("SN:", "") + "\n"
				wchrlistfile.write(chrid)
		wchrlistfile.close()
		rgenomedictfile.close()

	def Denovo_config(self, Denovo_config_file, script_dir, dataset_dir):
		wDenovo_config_file = open(Denovo_config_file, 'w')
		script_abs_path = os.path.abspath(sys.argv[0])

		config_line = "\n".join(["###denovo assemble", "", ""])
		wDenovo_config_file.write(config_line)

		wDenovo_config_file.close()

def usage():
	create_usage = \
	'''
	create configuration file for LRTK
	Version: 1.0.0
	Dependents: Python (>=3.0)
	Last Updated Date: 2017-06-01
	Contact: meijp@foxmail.com

	Usage: python create_config.py <command> [options]

	Command:
		Basic       create configuration file for basic processing
		Reseq       create configuration file for resequencing
		Denovo      create configuration file for de novo assembly
		all         create all configuration files

	Options:
		-o --outputdir, the path of output directory
		-s --softwarepath, the path of directory where all software were installed
		-d --datasetpath, dataset directory
		-h --help, help info

	'''
	print(create_usage)

if __name__ == '__main__':
	if len(sys.argv) < 4:
		usage()
		sys.exit(-1)
	else:
		command = sys.argv[1]
		outputdir = None
		SoftwarePathDir = None
		DatasetPathDir = None

		commandlist = ["Basic", "Reseq", "Denovo", "all"]
		if command not in commandlist:
			sys.stderr.write(" %s is not supported in command!\n\n" % command)
			usage()
			sys.exit(-1)

		opts, args = getopt.gnu_getopt(sys.argv[2:], 'o:s:d:h:', ['outputdir', 'softwarepath', 'datasetpath', 'help'])
		for o, a in opts:
			if o == '-o' or o == '--outputdir':
				outputdir = str(a)
			if o == '-s' or o == '--softwarepath':
				SoftwarePathDir = str(a)
			if o == '-d' or o == '--datasetpath':
				DatasetPathDir = str(a)
			if o == '-h' or o == '--help':
				usage()
				sys.exit(-1)

		C = create_config()
		if command == "Basic" or command == "all":
			config_file = os.path.join(outputdir, "Basic.config")
			C.Basic_config(config_file, SoftwarePathDir, DatasetPathDir)
			sys.stderr.write("[ %s ] configuration file for quality control has been created: %s \n" % (time.asctime(), config_file))
		if command == "Reseq" or command == "all":
			config_file = os.path.join(outputdir, "Reseq.config")
			C.Reseq_config(config_file, SoftwarePathDir, DatasetPathDir)
			sys.stderr.write("[ %s ] configuration file for resequencing has been created: %s \n" % (time.asctime(), config_file))
		if command == "Denovo" or command == "all":
			config_file = os.path.join(outputdir, "Denovo.config")
			C.Denovo_config(config_file, SoftwarePathDir, DatasetPathDir)
			sys.stderr.write("[ %s ] configuration file for de novo assembly has been created: %s \n" % (time.asctime(), config_file))
