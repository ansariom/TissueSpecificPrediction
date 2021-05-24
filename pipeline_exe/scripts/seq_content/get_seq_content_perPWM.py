import numpy as np
import sys

fname = "../ibdc_roe-only/peaks_3000_region_oneLine.fa"
indir = "seq_content_tests/"
suffix = "pwms_set1"
pwm_infile = indir + suffix + ".txt"

half = 3000
#tile = 20
nucs_up = 200
nucs_down = 50
at_content = True
gc_content = False
outfile = sys.argv[1]

def load_matrices(infile):
    mtx_dict = {}
    mtx_file = open(infile)
    for line in iter(mtx_file):
	line = line.rstrip()
        if (line.startswith(">")):
            pwm_name = line[2:]
            mtx_arr = []
            continue
        if(line.startswith("=")):
            line_parts = line.split("\t")
            mtx_arr.append(line_parts[1:])
        else:
            arr = np.array(mtx_arr, dtype=float)
            mtx_dict[pwm_name] = arr
    return mtx_dict

def getSeqContent(line, chars, w):
    sum = 0.0
    for ch in chars:
	sum = sum + line.count(ch)
    return sum/w

with open(fname) as f:
    content = f.readlines()

content = [x.strip() for x in content]
i = 0
all_arr = []

header = ","

pwms = load_matrices(pwm_infile)
w = -1 * nucs_up
e = nucs_down + 1
while w < e:
    for pwm in iter(pwms.keys()):
	if at_content:
	    header = header + "," + pwm +  "_TAContent_" + str(w)
	if gc_content:
	    header = header + "," + pwm + "_GCContent_" + str(w)
    w = w + 1

while i < len(content):
    l = content[i]
    arr = []
    if l.startswith(">"):
	seqName = l[1:]
	i = i + 1
	l = content[i]
	i = i + 1
	start = half - nucs_up
    	to = start
    	arr.append(seqName)
	end = half + nucs_down + 1
    	while start < end:
	    for pwm in iter(pwms.keys()):
		mtx = pwms[pwm]
		to = to + mtx.shape[0] + 1
		if at_content:
	    	    feature = getSeqContent(l[start:to], 'AT', mtx.shape[0])
	    	    arr.append(feature)
	        if gc_content:
		    feature = getSeqContent(l[start:to], 'GC', mtx.shape[0])
		    arr.append(feature)
	    start = start + 1
	    to = start
	all_arr.append(np.asanyarray(arr))

all_arr = np.asanyarray(all_arr)
np.savetxt(outfile, all_arr, delimiter = ",", fmt = "%s", header = header)


