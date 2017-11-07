import os, sys, gzip
import getopt
import time
import subprocess
import re
import random, string

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

	def Ref(self):
		abs_path = self.get_path('ref')
		return abs_path

	def Fq_aln_par(self):
		abs_path = self.get_path('fq_aln_parameter')
		return abs_path

	def Fq_sam_par(self):
		abs_path = self.get_path('fq_sampe_parameter')
		return abs_path

	def Bwa(self):
		abs_path = self.get_path('bwa')
		return abs_path

	def Samtools(self):
		abs_path = self.get_path('samtools')
		return abs_path

class modify_barcode:
	def find_samtools_version(self, samtools_path, outdir = './'):
		randomstring = ''.join(random.sample(string.ascii_letters + string.digits, 8))
		print(randomstring)
		tmpshell = outdir + "/" + randomstring + ".sh"
		tmplog = outdir + "/" + randomstring + ".log"
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

	def fill_barcode(self, original_sam, samtools_path, outdir = './'):
		original_name = os.path.basename(original_sam)

		new_sam = original_name.replace('sam', 'new.sam')
		new_sam = os.path.join(outdir, new_sam)
		new_bam = new_sam.replace('sam', 'bam')
		sorted_bam = new_sam.replace('new.sam', 'sorted')

		print(" ".join(["modified and sorted bam file:", sorted_bam, ".bam\n"]))
		
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

		shelldir = os.path.dirname(sorted_bam) + "/shell"
		check_info(shelldir)
		shell = shelldir + "/" + os.path.basename(sorted_bam) + ".sh"
		oshell = open(shell, 'w')
		sv = self.find_samtools_version(samtools_path, shelldir)
		if sv == 0:
			shell_line = " ".join([samtools_path, "view -h -S -b", new_sam, ">", new_bam, "\n"])
			oshell.write(shell_line)
			shell_line = " ".join([samtools_path, "sort -m 1G", new_bam, sorted_bam, "\n"])
			oshell.write(shell_line)
		else:
			shell_line = " ".join([samtools_path, "view -h -S -b -o", new_bam, new_sam + "\n"])
			oshell.write(shell_line)
			shell_line = " ".join([samtools_path, "sort -m 2G -o", sorted_bam + ".sorted.bam -T", sorted_bam, new_bam, "\n"])
			oshell.write(shell_line)
			shell_line = " ".join(["mv", sorted_bam + ".sorted.bam", sorted_bam + ".bam\n"])
			oshell.write(shell_line)
		shell_line = " ".join(["#rm", original_sam, new_sam, new_bam, "\n"])
		oshell.write(shell_line)
		oshell.close()

		subprocess.call(["sh", shell])
		return(sorted_bam + '.bam', sorted_bam + '.sam')

def check_info(result):
	if os.path.isdir(result):
		pass
	else:
		os.makedirs(result)

class OUTERSOFT:
	def bwa_fq(self, fq1, fq2, outprefix, bwa_path, bwa_aln_parameter, bwa_sam_parameter, RGinfo):
		sai1 = outprefix + '.1.sai'
		sai2 = outprefix + '.2.sai'
		outsam = outprefix + '.sam'
		
		sys.stderr.write("[ %s ] fq alignment starts\n" % (time.asctime()))
		outdir = os.path.dirname(fq1)
		shelldir = os.path.join(outdir, "shell")
		check_info(shelldir)
		shell = os.path.join(shelldir, "fq.aln.sh")
		oshell = open(shell, 'w')
		shell_line = " ".join([bwa_path, 'aln', ref_path, fq1, bwa_aln_parameter, '-f', sai1, '&\n'])
		oshell.write(shell_line)
		shell_line = " ".join([bwa_path, 'aln', ref_path, fq2, bwa_aln_parameter, '-f', sai2, '&\nwait\n'])
		oshell.write(shell_line)
		#shell_line = " ".join([bwa_path, 'sampe', ref_path, sai1, sai2, fq1, fq2, "|", samtools_path + ' view -h -S -b - > ' + outbam + "\n"])
		if RGinfo != 0:
			RGinfo = '"' + RGinfo + '"'
			shell_line = " ".join([bwa_path, 'sampe', bwa_sam_parameter, '-r', RGinfo, ref_path, sai1, sai2, fq1, fq2, ">", outsam + "\nrm", sai1, sai2, "\n"])
		else:
			shell_line = " ".join([bwa_path, 'sampe', bwa_sam_parameter, ref_path, sai1, sai2, fq1, fq2, ">", outsam + "\nrm", sai1, sai2, "\n"])
		oshell.write(shell_line)
		oshell.close()
		print(shell)

		subprocess.call(["sh", shell])
		sys.stderr.write("[ %s ] fq alignment finish\n" % (time.asctime()))
		return(outsam)

def usage():
	alignment_usage = \
	'''
	alignment and calculate
	Version: 1.0.0
	Dependents: Python (>=3.0), BWA, SAMtools, Picard (>=2.0)
	Last Updated Date: 2017-06-01
	Contact: meijp@foxmail.com

	Usage: python alignment.py <options>

	Options:
		-i --input, prefix of pair-end fq file, input_1.fq.gz and input_2.fq.gz must be valid
		-o --outputdir, the path of output directory
		-r --RG, RG info for alignment
		-c --config, the path of configuration file [default: outdir/config/Basic.config]
		-h --help, help info

	'''
	print(alignment_usage)

if __name__ == '__main__':
	if len(sys.argv) < 5:
		usage()
		sys.exit(-1)

	inputfq_prefix = None
	outputdir = None
	RGInfo = 0
	ConfigFile = None
	opts, args = getopt.gnu_getopt(sys.argv[1:], 'i:o:c:r:h:', ['input', 'outputdir', 'config', 'RG', 'help'])
	for o, a in opts:
		if o == '-i' or o == '--input':
			inputfq_prefix = a
		if o == '-o' or o == '--outputdir':
			outputdir = a
		if o == '-c' or o == '--config':
			ConfigFile = a
		if o == '-r' or o == '--RG':
			RGInfo = a
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
		shell_line = " ".join(["python", create_config_py, "Basic -o", config_dir, "\n"])
		wtmpshell.write(shell_line)
		wtmpshell.close()
		subprocess.call(["sh", tmpshell])
		subprocess.call(["rm", tmpshell])
		ConfigFile = os.path.join(config_dir, "Basic.config")

	if RGInfo == 0:
		sys.stderr.write("\n\nWarnings: RG info was not defined, this might cause unpredictable error in subsequent blocks!\n\n")

	O = OUTERSOFT()

	G = baseinfo()
	G.get_config(ConfigFile)

	ref_path = G.Ref()
	bwa_path = G.Bwa()
	bwa_aln_parameter = G.Fq_aln_par()
	bwa_sam_parameter = G.Fq_sam_par()
	samtools = G.Samtools()

	fq1 = inputfq_prefix + "_1.fq.gz"
	fq2 = inputfq_prefix + "_2.fq.gz"
	newprefix = inputfq_prefix
	sys.stderr.write("[ %s ] align to the genome %s\n" % (time.asctime(), ref_path))
	original_sam = O.bwa_fq(fq1, fq2, newprefix, bwa_path, bwa_aln_parameter, bwa_sam_parameter, RGInfo)
	sys.stderr.write("[ %s ] alignment finished. Original sam file has been writeen to %s\n\n" % (time.asctime(), original_sam))

	M = modify_barcode()
	samdir = os.path.dirname(original_sam)
	sys.stderr.write("[ %s ] modify barcode info of original sam file %s\n" % (time.asctime(), original_sam))
	(Sortedbam, Sortedsam) = M.fill_barcode(original_sam, samtools, samdir)
	sys.stderr.write("[ %s ] barcode info has been modified, new bam file has been written to %s\n\n" % (time.asctime(), Sortedbam))

	tmpdir = samdir + "/tmp"
	if os.path.isdir(tmpdir):
		pass
	else:
		os.mkdir(tmpdir)
	mvshell = tmpdir + "/tmp.mv.2.sh"
	wmvshell = open(mvshell, 'w')
	shell_line = " ".join(["mv", fq1, fq2, tmpdir, "\n"])
	wmvshell.write(shell_line)
	shell_line = " ".join(["mv", original_sam, tmpdir, "\n"])
	wmvshell.write(shell_line)
	shell_line = " ".join(["mv", inputfq_prefix + ".new.sam", inputfq_prefix + ".new.bam", tmpdir, "\n"])
	wmvshell.write(shell_line)
	wmvshell.close()
	subprocess.call(["sh", mvshell])
