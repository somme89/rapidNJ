import os
import sys
from subprocess import PIPE, Popen, STDOUT

cmd = './make_phylip_trees.sh'
p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
output = p.stdout.read()

cmd = '../x64/Release/bin/rapidnj_x64.exe LGT_70_short_names.fa -a jc -x LGT_70_short_names_rnj_tree.nwk'
p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
output = p.stdout.read()

cmd = 'python strip_quotes.py LGT_70_short_names_rnj_tree.nwk > LGT_70_short_names_rnj_tree_no_quotes.nwk '
p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
output = p.stdout.read()

cmd = '../tqdist/quartet_dist.exe LGT_70_short_names_rnj_tree_no_quotes.nwk LGT_70_short_names_phylip_tree.nwk'
p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
output = p.stdout.read()
if int(output) != 0:
	print "ERROR: LG_70 tree distance was not 0 but", output
	sys.exit(1)
else:
	print "LG_70 SUCCESS!"

os.remove("LGT_70_short_names_rnj_tree.nwk")
os.remove("LGT_70_short_names_rnj_tree_no_quotes.nwk")
os.remove("LGT_70_short_names_phylip_tree.nwk")
os.remove("LGT_70_short_names.fa")
os.remove("LGT_70_short_names.phylip")

cmd = '../x64/Release/bin/rapidnj_x64.exe PF15224_full_short_names.fa -x PF15224_full_short_names_rnj_tree.nwk'
p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
output = p.stdout.read()

cmd = 'python strip_quotes.py PF15224_full_short_names_rnj_tree.nwk > PF15224_full_short_names_rnj_tree_no_quotes.nwk '
p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
output = p.stdout.read()

cmd = '../tqdist/quartet_dist.exe PF15224_full_short_names_rnj_tree_no_quotes.nwk PF15224_full_short_names_phylip_tree.nwk'
p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
output = p.stdout.read()
if int(output) != 112:
	print "ERROR: PF15224 tree distance was not 112 but", output
	sys.exit(1)
else:
	print "PF15224 SUCCESS!"

os.remove("PF15224_full_short_names_rnj_tree.nwk")
os.remove("PF15224_full_short_names_rnj_tree_no_quotes.nwk")
os.remove("PF15224_full_short_names_phylip_tree.nwk")
os.remove("PF15224_full_short_names.fa")
os.remove("PF15224_full_short_names.phylip")

os.remove("outfile")
os.remove("infile")
os.remove("input")
