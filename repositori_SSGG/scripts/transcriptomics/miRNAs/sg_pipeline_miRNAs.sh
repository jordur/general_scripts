#!/bin/bash
#Input: 1) $1=nombre muestra 2)#$2=fichero fastq (ruta completa)
#1) Eliminar los nucleótidos del adaptador en el extremo 3'.
#/share/apps/Python-2.6/python /share/apps/cutadapt-0.9.3/cutadapt -c -e 0.12 -a 330201030313112312 --bwa -x TRIM: $1_F3.csfasta $1_F3_QV.qual > $1_sinAdaptador.fastq
#2) Mapear las lecturas con BWA.
# -l = longitud del seed (son miRNAs, el seed tiene que ser pequenyo).
# -k = numero de mismatches en el seed (pongo 2 porque trabajamos con color space y para que existe 1 polimorfismo hemos de ver 2 mismatches)
#/share/apps/bwa-0.5.9/bwa aln -c -l 16 -k 2 -t 8 -n 4 /data/results/Solid0065/referencia_genoma_humano/GRch37.58/hg19.fa $2 > $1.sai
#/share/apps/bwa-0.5.9/bwa samse -f $1.sam /data/results/Solid0065/referencia_genoma_humano/GRch37.58/hg19.fa $1.sai $2
#3) Eliminar las secuencias que mapean en multiples sitios en funcion de los valores de mapeo (PhredQV=20) y del flag XA (alineamientos secundarios de las lecturas).
#samtools view -S -q 20 $1.sam > $1_Q20.sam
#grep -v "XA:" $1_Q20.sam > $1_Q20_sin_XA.sam
awk '{if ((/XM:i:0/) || (/XM:i:1/) || (/XM:i:2/) || (/XM:i:3/))  print $0}' ../$1_Q20_sin_XA.sam > $1_Q20_sin_XA_3mm.sam_tmp
cat /data/results/Solid0065/BF41_microRNAs_FM/secondary/15/OK/cab $1_Q20_sin_XA_3mm.sam_tmp > $1_Q20_sin_XA_3mm.sam 
samtools view -S -b $1_Q20_sin_XA_3mm.sam > $1_Q20_sin_XA_3mm.bam
coverageBed -abam $1_Q20_sin_XA_3mm.bam -b /data/results/Solid0065/BF41_microRNAs_FM/datos_partida/repeats_UCSC.bed > $1_repeats
coverageBed -abam $1_Q20_sin_XA_3mm.bam -b /data/results/Solid0065/BF41_microRNAs_FM/datos_partida/Homo_sapiens.bed -s > $1_genes
coverageBed -abam $1_Q20_sin_XA_3mm.bam -b /data/results/Solid0065/BF41_microRNAs_FM/datos_partida/Pfam.bed -s > $1_Pfam
coverageBed -abam $1_Q20_sin_XA_3mm.bam -b /data/results/Solid0065/BF41_microRNAs_FM/datos_partida/hsa_mod.gff -s > $1_miRNAs_conocidos
awk '{if ($11>20) print $10}' $1_miRNAs_conocidos | sed 's/ID="//' | sed 's/";//' > $1_nombres_miRNAs
#/data/results/Solid0065/BF41_microRNAs_FM/scripts/sg_extraer_miRNAs_en_muestra.pl $1_nombres_miRNAs /data/results/Solid0065/BF41_microRNAs_FM/datos_partida/DNA_humano_miRNAs_maduros.fasta > maduros
#/data/results/Solid0065/BF41_microRNAs_FM/scripts/sg_extraer_miRNAs_en_muestra.pl $1_nombres_miRNAs /data/results/Solid0065/BF41_microRNAs_FM/datos_partida/DNA_humano_miRNA_star.fasta > star
a
#4) Convertir el fichero sam en arf
bwa_sam_converter.pl $1_Q20_sin_XA_3mm.sam
#5) Contabilizar el número de veces que aparece cada miRNA.
collapse_reads.pl reads.fa hsa > reads_collapsed.fa
#6) Eliminar secuencias menores a 17nt
/share/apps/mirdeep2/fastaparse.pl reads_collapsed.fa -a 17 > reads_collapsed_17nt.fa
#7) Ejecutar wrapper de miRDeep
miRDeep2.pl reads_collapsed_17nt.fa /data/results/Solid0065/referencia_genoma_humano/GRch37.58/hg19.fa reads_vs_genome.arf DNA_humano_miRNA_maduros.fasta /data/results/Solid0065/BF41_microRNAs_FM/datos_partida/DNA_miRNAs_cercanos_humano_maduro.fasta /data/results/Solid0065/BF41_microRNAs_FM/datos_partida/DNA_miRNAs_horquillas_humano.fasta -t Human -s /data/results/Solid0065/BF41_microRNAs_FM/datos_partida/DNA_humano_miRNA_star.fasta 2> report.log
