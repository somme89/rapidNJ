import os
import sys

for line in open(sys.argv[1]):
	line = line.replace("'","")
	print line,