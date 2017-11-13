import os, sys
import re

if len(sys.argv) != 3:
	print("python %s fa nonN.bed" % sys.argv[0])
	sys.exit(-1)

fa = sys.argv[1]
out = sys.argv[2]

rfa = open(fa, 'r')
wout = open(out, 'w')

chrid = None
start = 0
end = 0
last_str = None
all_length = 0

for fainfo in rfa:
	fainfo = fainfo.strip()
	if re.search(">", fainfo):
		if start > 0:
			outinfo = str(chrid) + "\t" + str(start) + "\t" + str(end) + "\n"
			wout.write(outinfo)
		fainfolist = re.split("\s", fainfo)
		chrid = fainfolist[0].replace(">", "")
		start = 0
		all_length = 0
	else:
		fainfo = str(fainfo)
		N = fainfo.count("N")
		falength = len(fainfo)
		if start == 0:
			if N > 0:
				if N < falength:
					start = all_length + N + 1
				else:
					start = 0
			else:
				start = all_length + 1
		else:
			if N > 0:
				end = all_length + falength - N
				outinfo = str(chrid) + "\t" + str(start) + "\t" + str(end) + "\n"
				wout.write(outinfo)
				start = 0
			else:
				end = all_length + falength
		all_length += falength
if start > 0:
	outinfo = str(chrid) + "\t" + str(start) + "\t" + str(end) + "\n"
	wout.write(outinfo)
wout.close()
rfa.close()
