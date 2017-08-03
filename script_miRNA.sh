tmap mapall -i fastq -f /home/jordi/Feina/crotalus/Neonato/assembly/blastresult/blast/phosphosdiest.class.txt.fa -r ../../over_neonato.fastq -a 3 stage1 map1 stage2 map2 map3 mapvsw > resultats_miRNA.sam
samtools view -bS resultats_miRNA.sam > resultats_miRNA.bam
samtools sort resultats_miRNA.bam resultats_miRNA.sorted
samtools index resultats_miRNA.sorted.bam
