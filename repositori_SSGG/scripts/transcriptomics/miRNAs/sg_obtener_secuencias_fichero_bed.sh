#!/bin/bash

i=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M)

for a in ${i[@]}
do
	awk '{if ($1=="chr'$a'") print $4"\t"$2"\t"$3"\t"$6}' $1 > $a\_$1_tmp
	sg_extract_seq.pl $a\_$1_tmp /data/results/Solid0065/referencia_genoma_humano/GRch37.58/chr$a > $a\_$1\_seq.fasta_tmp
done
cat *_seq.fasta_tmp > $1\_final.fasta
rm *_tmp
