import sys
import os

counter = 0
for line in open(sys.argv[1]):
	if line.startswith(">"):
		print ">seq_" + str(counter)
		counter += 1
	else:
		print line,