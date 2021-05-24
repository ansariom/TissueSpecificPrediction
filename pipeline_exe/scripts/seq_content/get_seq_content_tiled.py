import numpy as np
import sys

args = sys.argv[1:] 
#fname = "../ibdc_roe-only/peaks_3000_region_oneLine.fa"
fname = args[4]

half = 3000
tile = int(args[0])
nucs_up = int(args[1])
nucs_down = int(args[2])
nwins = (nucs_down + nucs_up) / tile
print("num wins = " + str(nwins))
at_content = True
gc_content = False
outfile = args[3]

def getSeqContent(line, chars):
    sum = 0.0
    for ch in chars:
	sum = sum + line.count(ch)
    return sum/tile

print "start loading file .."
with open(fname) as f:
    content = f.readlines()
print "file loaded!"

content = [x.strip() for x in content]
print len(content)

i = 0
all_arr = []

header = ","

for w in range(0, nwins):
    if at_content:
	header = header + ",TAContent" + str(w+1)
    if gc_content:
	header = header + ",GCContent" + str(w+1)
print "header created!"

while i < len(content):
    l = content[i]
    arr = []
    if l.startswith(">"):
	seqName = l[1:] + "_0"
	print seqName
	i = i + 1
	l = content[i]
	i = i + 1
	start = half - nucs_up
    	to = start + tile + 1
    	arr.append(seqName)
	end = half + nucs_down + 2
    	while to < end:
	    if at_content:
	    	feature = getSeqContent(l[start:to], 'AT')
	    	arr.append(feature)
	    if gc_content:
		feature = getSeqContent(l[start:to], 'GC')
		arr.append(feature)
	    start = to
            to = to + tile
	all_arr.append(np.asanyarray(arr))

all_arr = np.asanyarray(all_arr)
np.savetxt(outfile, all_arr, delimiter = ",", fmt = "%s", header = header)


