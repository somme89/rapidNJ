#!/bin/bash
python shorten_tree_names.py LGT_70.fa > LGT_70_short_names.fa
./fasta2phylip.pl LGT_70_short_names.fa > LGT_70_short_names.phylip
cp LGT_70_short_names.phylip infile
rm outfile
echo -e "D" > input
echo -e "D" >> input
echo -e "Y" >> input
../phylip_package/dnadist.exe < input
mv outfile infile
rm outtree
echo -e "Y" > input
../phylip_package/neighbor.exe < input
mv outtree LGT_70_short_names_phylip_tree.nwk

python shorten_tree_names.py PF15224_full.fa > PF15224_full_short_names.fa
./fasta2phylip.pl PF15224_full_short_names.fa > PF15224_full_short_names.phylip
cp PF15224_full_short_names.phylip infile
rm outfile
echo -e "P" > input
echo -e "P" >> input
echo -e "P" >> input
echo -e "Y" >> input
../phylip_package/protdist.exe < input
mv outfile infile
rm outtree
echo -e "Y" > input
../phylip_package/neighbor.exe < input
mv outtree PF15224_full_short_names_phylip_tree.nwk