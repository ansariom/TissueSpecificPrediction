'''
Created on Dec 9, 2017

@author: mitra
'''

import numpy as np
import math
import sys
from numpy import dtype
bg_M0 = [0.3486, 0.1599,    0.1498,    0.3417]
bg_M1 = [0.4020, 0.1501, 0.1419, 0.3060,
0.3708, 0.1786, 0.1252, 0.3254,
0.3595, 0.1434, 0.1693, 0.3279,
0.2788, 0.1685, 0.1607, 0.3920]

bg_M1 = np.asarray(bg_M1).reshape(4,4)
print(bg_M1)

di_nt= "ACGT"

#pwm_file = "/Users/mitra/Downloads/pwms_peat_core_agris_huges_dec2017.mat"
pwm_file = sys.argv[1]
outfile = sys.argv[2]

line = ""
pwm_name = ""
s = []
freq = []

#out = open("../max_loglik_scores_for_pwms.txt", "w")
out = open(outfile, "w")

f = open(pwm_file, "r")
for line in iter(f):
    line = line.rstrip()
    if line.startswith(">"):
        
        num_log = 0
        dnum_log = 0
        if len(s) > 0:
            for i in xrange(len(s)):
                if (i == 0):
                    bg_freq = bg_M0[s[i]]
                else:
                    bg_freq = bg_M1[s[i-1],s[i]]
                num_log = num_log + math.log(freq[i])
                dnum_log = dnum_log + math.log(bg_freq)
        loglik = num_log - dnum_log
        # print(pwm_name + "\t" + str(loglik))
        out.write(pwm_name + "\t" + str(loglik) + "\n")
        pwm_name = line[2:]
        s = []
        freq = []
    elif line.startswith("="):
        current_pos = np.asarray(line[2:].split('\t'), dtype = float)
        freq.append(np.amax(current_pos))
        s.append(np.argmax(current_pos))
        
# the last one
if len(s) > 0:
    for i in xrange(len(s)):
        if (i == 0):
            bg_freq = bg_M0[s[i]]
        else:
            bg_freq = bg_M1[s[i-1],s[i]]
        num_log = num_log + math.log(freq[i])
        dnum_log = dnum_log + math.log(bg_freq)
    loglik = num_log - dnum_log
    # print(pwm_name + "\t" + str(loglik))
    out.write(pwm_name + "\t" + str(loglik) + "\n")
out.close()
f.close()
