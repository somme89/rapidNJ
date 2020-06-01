import sys
import os
import string
import random
import time

if len(sys.argv) > 3 or len(sys.argv) < 2:
    print "USAGE: out_dir [-g]"
    sys.exit()
out_dir = sys.argv[1]
useGaps = 0
if len(sys.argv) == 3:
    if sys.argv[2] == "-g":
        useGaps = 1
    else:
        print "USAGE: out_dir [-g]"
        sys.exit()

seq_lengths = range(100000,3100000,100000)
seq_lengths.append(10000)
seq_lengths.append(50000)
seq_count = 100

for length in seq_lengths:
    print length
    if(useGaps):
        os.system("bin/sim_seq " + str(seq_count) + " " + str(length) + " d s " + " -g > " + out_dir + "/seq_" + str(seq_count) + "_" + str(length) + ".sth")
    else:
        os.system("bin/sim_seq " + str(seq_count) + " " + str(length) + " d s " + " > " + out_dir + "/seq_" + str(seq_count) + "_" + str(length) + ".sth")


seq_lengths = range(10000,200000,10000)
seq_lengths.append(100)
seq_lengths.append(1000)
seq_lengths.append(5000)
seq_count = 1000

for length in seq_lengths:
    print length
    if(useGaps):
        os.system("bin/sim_seq " + str(seq_count) + " " + str(length) + " d s " + " -g > " + out_dir + "/seq_" + str(seq_count) + "_" + str(length) + ".sth")
    else:
        os.system("bin/sim_seq " + str(seq_count) + " " + str(length) + " d s " + " > " + out_dir + "/seq_" + str(seq_count) + "_" + str(length) + ".sth")

seq_lengths = range(5000,100000,5000)
seq_lengths.append(50)
seq_lengths.append(500)
seq_lengths.append(2500)
seq_count = 2000

for length in seq_lengths:
    print length
    if(useGaps):
        os.system("bin/sim_seq " + str(seq_count) + " " + str(length) + " d s " + " -g > " + out_dir + "/seq_" + str(seq_count) + "_" + str(length) + ".sth")
    else:
        os.system("bin/sim_seq " + str(seq_count) + " " + str(length) + " d s " + " > " + out_dir + "/seq_" + str(seq_count) + "_" + str(length) + ".sth")

seq_lengths = range(10000,200000,10000)
seq_lengths.append(100)
seq_lengths.append(1000)
seq_lengths.append(5000)
seq_lengths = [x/3 for x in seq_lengths]
seq_count = 3000

for length in seq_lengths:
    print length
    if(useGaps):
        os.system("bin/sim_seq " + str(seq_count) + " " + str(length) + " d s " + " -g > " + out_dir + "/seq_" + str(seq_count) + "_" + str(length) + ".sth")
    else:
        os.system("bin/sim_seq " + str(seq_count) + " " + str(length) + " d s " + " > " + out_dir + "/seq_" + str(seq_count) + "_" + str(length) + ".sth")
