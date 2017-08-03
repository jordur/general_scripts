SAM and BAM filtering one-liners  

@author: David Fredman, david.fredmanAAAAAA@gmail.com (sans poly-A tail)  
@dependencies: http://sourceforge.net/projects/bamtools/ and http://samtools.sourceforge.net/

Please extend with additional/faster/better solutions via a pull request!


BWA mapping (using piping for minimal disk I/O)

	bwa aln -t 8 targetGenome.fa reads.fastq | bwa samse targetGenome.fa - reads.fastq\
	| samtools view -bt targetGenome.fa - | samtools sort - reads.bwa.targetGenome
	
	samtools index reads.bwa.targetGenome.bam

Count number of records (unmapped reads + each aligned location per mapped read) in a bam file:

	samtools view -c filename.bam

Count with flagstat for additional information:

	samtools flagstat filename.bam

Count the number of alignments (reads mapping to multiple locations counted multiple times)

	samtools view -F 0x04 -c filename.bam

Count number of mapped reads (not mapped locations) for left and right mate in read pairs

	samtools view -F 0x40 filename.bam | cut -f1 | sort | uniq | wc -l
	samtools view -f 0x40 -F 0x4 filename.bam | cut -f1 | sort | uniq | wc -l #left mate
	samtools view -f 0x80 -F 0x4 filename.bam | cut -f1 | sort | uniq  | wc -l #right mate

Remove unmapped reads, keep the mapped reads:

	samtools view -F 0x04 -b in.bam > out.aligned.bam

Count UNmapped reads:

	samtools view -f4 -c in.bam

Require minimum mapping quality (to retain reliably mapped reads):

	samtools view -q 30 -b in.bam > aligned_reads.q30.bam
	samtools view -q 30 -c in.bam #to count alignments with score >30

Require match to be on the sense strand of the reference (samtools flag)

	samtools view -F 16

Require match to be on antisense strand (samtools flag)

	samtools view -f 16

Require at least N matches at the start of the read:

	$N=6
	samtools view in.bam \
	| perl -lane 'next unless $F[5] =~ /^(\d+)M/;print if $1 >= $N;'

Filter by number of mismatches in BWA generated output, use BWA-specific flag:

Tag	Meaning  
NM     Edit distance  
MD     Mismatching positions/bases  
AS     Alignment score  
BC     Barcode sequence  
X0     Number of best hits  
X1     Number of suboptimal hits found by BWA  
XN     Number of ambiguous bases in the reference  
XM     Number of mismatches in the alignment  
XO     Number of gap opens  
XG     Number of gap extentions  
XT     Type: Unique/Repeat/N/Mate-sw  
XA     Alternative hits; format: (chr,pos,CIGAR,NM;)*  
XS     Suboptimal alignment score  
XF     Support from forward/reverse alignment  
XE     Number of supporting seeds  

To keep only reads that map without any mismatches:

	bamtools filter -tag XM:0 -in reads.bam -out reads.noMismatch.bam

Retain only uniquely mapping reads (reads with a single unambigous mapping location):

If BWA was used it is possible to use the BWA XT flag value U for unique (analogously, R is for repeat). I did not find a simple way to do this with samtools or bamtools, so grep to the rescue:

	samtools view reads.bam | grep 'XT:A:U' | samtools view -bS -T referenceSequence.fa - > reads.uniqueMap.bam

However, the concept of "uniquely mapping" is not the cleanest idea - in most scenarios any given read could be placed elsewhere although it may be a lower scoring alignment. Thus, you could instead filter based on mapping quality, to retain the "reliably mapped" reads. Different mappers have different scoring models. As a rule of thumb, min values of 5 or 10 will work well. If you used bowtie/bowtie2, try:

	samtools view -b -q 10 foo.bam > foo.filtered.bam
	