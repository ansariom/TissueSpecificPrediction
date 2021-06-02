
import sys
infasta = sys.argv[1]
outfile = sys.argv[2]

d = 80

with open(infasta) as f:
	seqs = f.readlines()

fw = open(outfile, 'w')

l = []
for i in range(1,len(seqs),2):
	l.append(seqs[i-1])
	s = seqs[i].upper()
	idx = 0
	while(idx < len(s) - 1):
		l.append(s[idx:idx+d])
		if (idx < (len(s) - 1)):
			l.append('\n')
		idx = idx + d

fw.writelines(l)
fw.close()




