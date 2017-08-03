#/bin/bash

i=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M)

for j in ${i[@]}
do
    sed 's/chr//' $1 | awk '{if ($1=="'$j'") print "chr"$0}' > $4/_tmp_target_$j
	sed 's/chr//' $2 | awk '{if ($1=="'$j'") print "chr"$0}' > $4/_tmp_indels_$j
	sg_extraer_indels_en_rango.pl $4/_tmp_target_$j $4/_tmp_indels_$j > $4/_tmp_indels-on-target_$j
done
head -n 1 $2 > $4\/$3
cat $4/_tmp_indels-on-target_* >> $4\/$3
rm $4/_tmp_target* $4/_tmp_indels*
