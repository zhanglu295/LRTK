import os, sys, gzip
import getopt
import time
import subprocess
import re

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
	def split_and_sort_sam(self, sorted_bam, samtools_path, barcode_index, outdir = './'):
		splitdir = os.path.join(outdir, "SamByChr")
		if os.path.exists(splitdir):
			pass
		else:
			os.mkdir(splitdir)

		sorted_sam = sorted_bam.replace('bam', 'sam')
		tmpshell = os.path.join(splitdir, "bam2sam.sh")
		wtmpshell = open(tmpshell, 'w')
		shell_line = " ".join([samtools_path, "view -h", sorted_bam, ">", sorted_sam, "\n"])
		wtmpshell.write(shell_line)
		wtmpshell.close()
		subprocess.call(["sh", tmpshell])

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
				saminfo = saminfo.strip()
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
			shell_line = " ".join(["sort", "-k12", chrsam, ">", chrsortedsam, "\n"])
			wtmpshell.write(shell_line)
			shell_line = " ".join(["rm", chrsam, "\n"])
			wtmpshell.write(shell_line)
			wtmpshell.close()
			subprocess.call(["sh", tmpshell])
			wsamfilepath.write("\n".join([chrsortedsam, ""]))
		wsamfilepath.close()
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
		marked_molecule = os.path.join(outdir, "molecule.txt.gz")
		wmarked_molecule = gzip.open(marked_molecule, 'wb')
		allbarcode = dict()
		mol_id = 0
		barid = None
		all_CF = 0
		all_CR = 0
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
#					sys.stderr.write("%s : unmapped reads would be ignored!\n" % saminfolist[0])
					pass
				else:
					for bc in saminfolist[11:]:
						if bc.startswith("BX"):
							barid = bc
							if bc not in allbarcode:
								if len(allbarcode) == 0:
									pass
								else:
									onekey = sorted(allbarcode)
#								print(allbarcode[onekey[0]])
									(mol_id, all_CF, all_CR) = self.FR(allbarcode[onekey[0]], mol_id, wmarked_molecule, barid, all_CF, all_CR)
									allbarcode = dict()
								allbarcode[bc] = saminfo
							else:
								allbarcode[bc] = allbarcode[bc] + "\n" + saminfo
			samfile.close()
			onekey = sorted(allbarcode)
		if len(onekey) > 0:
			(mol_id, all_CF, all_CR) = self.FR(allbarcode[onekey[0]], mol_id, wmarked_molecule, barid, all_CF, all_CR)
		wmarked_molecule.close()

		cfcr_stat_file = os.path.join(outdir, "CFCR.stat")
		wcfcr_stat_file = open(cfcr_stat_file, 'w')
		stat_line = "Item" + "\t" + "CF/CR length" + "\t" + "genome_size" + "\t" + "CF/CR depth\n"
		wcfcr_stat_file.write(stat_line)
		target_size = int(target_size)
		CF_depth = all_CF / float(target_size)
		CR_depth = all_CR / float(target_size)
#		k1 = " ".join([str(all_CF), str(target_size), str(CF_depth)])
#		print(k1)
#		k2 = " ".join([str(all_CR), str(target_size), str(CR_depth)])
#		print(k2)
		stat_line = "CF:" + "\t" + str(all_CF) + "\t" + str(target_size) + "\t" + str(CF_depth) + "\n"
		wcfcr_stat_file.write(stat_line)
		stat_line = "CR:" + "\t" + str(all_CR) + "\t" + str(target_size) + "\t" + str(CR_depth) + "\n"
		wcfcr_stat_file.write(stat_line)
		wcfcr_stat_file.close()
		return(cfcr_stat_file)

	def FR(self, barcodeinfo, molecule_id, writemarker, barcodeid, cf, cr):
		readid = dict()
		onebarcodelist = re.split("\n", barcodeinfo)
		for onebarcode in onebarcodelist:
			eachlist = re.split("\t", onebarcode)
			eachid = eachlist[0]
			chr_id = eachlist[2]
			start = eachlist[3]
			readlength = len(eachlist[9])
			start = int(start)
			end = start + readlength - 1
			reads = onebarcode
			### merge pair-end reads into single fragment
			if eachid not in readid:
				readid[eachid] = [chr_id, start, end, onebarcode]
				cr = cr + readlength
			else:
				if readid[eachid][0] != chr_id:
#					sys.stderr.write("%s pair-end reads mapped to different chromosom, would be ignored.\n" % eachid)
					readid.pop(eachid)
					cr = cr - readlength
				else:
					start = min(readid[eachid][1], start)
					end = max(readid[eachid][2], end)
					dis = end - start + 1
					cr = cr + readlength
					if dis > 10000:
#						sys.stderr.write("%s fragment length > 10kb, would be ignored.\n" % eachid)
						readid.pop(eachid)
					else:
						reads = readid[eachid][3] + "\n" + reads
						readid[eachid] = [chr_id, start, end, reads]

		molecule = dict()
		s = int(molecule_id) + 1
		molecule_id = int(molecule_id) + 1
		e = s + 1
		for key in readid:
			value = readid[key]
			if molecule_id not in molecule:
				molecule[molecule_id] = list()
				molecule[molecule_id].append(value[0])	## chr id
				molecule[molecule_id].append(value[1])  ## start
				molecule[molecule_id].append(value[2])  ## end
				molecule[molecule_id].append(value[3])  ## reads info
			else:
				isnew = 1
				for i in range(s, e):
					if isnew == 0:
						pass
					else:
						if molecule[i][0] == value[0]:
							dist = molecule[i][1] - value[1]
							if dist < 50000:
								isnew = 0
								molecule[i][1] = min(molecule[i][1], value[1])
								molecule[i][2] = max(molecule[i][2], value[2])
								molecule[i][3] = molecule[i][3] + "\n" + value[3]
				if isnew == 1:
					molecule_id = molecule_id + 1
					molecule[molecule_id] = list()
					molecule[molecule_id].append(value[0])
					molecule[molecule_id].append(value[1])
					molecule[molecule_id].append(value[2])
					molecule[molecule_id].append(value[3])
#					molecule[molecule_id]["chrid"] = value[0]
#					molecule[molecule_id]["start"] = value[1]
#					molecule[molecule_id]["end"] = value[2]
#					molecule[molecule_id]["reads"] = value[3]
					s = s + 1

		for i in range(s, e):
#			readlist = re.split("\n", molecule[i]["reads"])
			if i in molecule:
#				print(molecule[i][3])
				readlist = re.split("\n", molecule[i][3])
				cf = cf + molecule[i][2] - molecule[i][1] + 1
				for rl in readlist:
					outinfo = str(i) + "\t" + barcodeid + "\t" + rl + "\n"
					writemarker.write(outinfo.encode())

		return(molecule_id, cf, cr)

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
