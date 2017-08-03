#!/usr/bin/perl -w
# 2011/09/08 - árbol sw.
use File::ReadBackwards;
use File::Basename;

sub usage {
	print "\nHow to use: sg_Qvalues_SOLiD.pl <csfasta> <qual>\n\n";
	print "Description: this perl function calculates the quality values of a .csfasta file and its associated .qual file. For computing reasons only a maximum of 10M sequences should be taken from the original .csfasta and .qual files.\n\n";
	print "INPUT: <csfasta> and <qual> are the .csfasta and .qual files from a given strand.\n\n";
	print "OUTPUT: <base_file_name_Qstats.tsv> is a two tab-separated file containing the stats of each read position and the stats of each sequence\n";
	print "Please notice that the ouput file will be placed in the same directory than the input files, and with the same file extension\n";
	exit(1);
}

if(scalar(@ARGV) != 2){
    usage();
}

my $max_number_seqs=10000000; #Maximum number of sequences to be taken from the original .csfasta and .qual files
my $max_QV=50; #Maximum quality value of the .qual file

my $line_file; #auxiliar variable for printing lines to file substituting dots with commas
my $csfasta_name=$ARGV[0];
my $qual_name=$ARGV[1];

# The number of reads in the csfasta and qual files will be obtained and compared, so as to determine if both numbers of reads match together:
my $csfasta_seqs=`grep -c ">" $csfasta_name 2>&1`;
$csfasta_seqs=~ s/\r|\n//g; # Number of sequences in the csfasta file
my $qual_seqs=`grep -c ">" $qual_name 2>&1`;
$qual_seqs=~ s/\r|\n//g; # Number of sequences in the qual file
my $csfasta_reads=`grep -v "#" $csfasta_name | head -2 | tail -1 | wc | awk '{print \$3-2}' 2>&1`;
$csfasta_reads=~ s/\r|\n//g; # Read length in the csfasta file
print "Number of sequences in csfasta file: $csfasta_seqs, in qual file: $qual_seqs \nNumber of reads per sequence: $csfasta_reads \n";

if ( $csfasta_seqs ne $qual_seqs ){
	print "The .csfasta and .qual number of sequences read are different!!! The programm will exit!!!";
	exit;
}

my $reduce_seqs=0;

# A maximum of 10M sequences will be considered for the calculation of statistics:
($basename_csfasta, $dir1, $ext1) = fileparse($ARGV[0], qr/\.[^.]*/);
($basename_qual, $dir2, $ext2) = fileparse($ARGV[1], qr/\.[^.]*/);
if ( $csfasta_seqs > $max_number_seqs ) {
	my $cmd=`sg_randomLineFileAdvanced.pl $csfasta_name $qual_name $max_number_seqs 10Mseqs_`;
	($basename_csfasta, $dir1, $ext1) = fileparse($ARGV[0], qr/\.[^.]*/);
	($basename_qual, $dir2, $ext2) = fileparse($ARGV[1], qr/\.[^.]*/);
	$csfasta_name=$dir1 . "10Mseqs_" . $basename_csfasta . $ext1;
	$qual_name=$dir2 . "10Mseqs_" . $basename_qual . $ext2;
	$csfasta_seqs=$max_number_seqs;
	$qual_seqs=$max_number_seqs;
	$reduce_seqs=1;
}

# input and output files get opened
open(csfasta_file, "<", $csfasta_name) || die("Could not open file!"); #$#ARGV]);
open(qual_file, "<", $qual_name) || die("Could not open file!");

($base, $dir) = fileparse($ARGV[0], qr/\.[^.]*/);
open(qstats_file, ">", $dir . $base . "_Qstats.tsv" ) || die("Could not open file!");
open(log_file, ">", $dir . $base . "_log.txt") || die("Could not open file!"); #Creation of output log file

# Creation of temporary directory:
my $cmd1=`mkdir $dir/temp_dir_$base`;

# Initialization of statistical variables
for ($i=0;$i<$csfasta_reads;$i++) {
	$position_quals_mean[$i]=0; # Mean values of qualities of the correctly read positions of the sequences
	$position_quals_max[$i]=0; # Maximum values of qualities of the correctly read positions of the sequences
	$position_quals_min[$i]=0; # Minimum values of qualities of the correctly read positions of the sequences
	$position_quals_median[$i]=0; # Median values of qualities of the correctly read positions of the sequences
	$position_quals_percent20[$i]=0; # 20-percentil of qualities of the correctly read positions of the sequences
	$position_quals_percent80[$i]=0; # 80-percentil of qualities of the correctly read positions of the sequences
	$position_quals_0_mean[$i]=0; # Mean values of qualities of the correctly read 0-base-pairs of the sequences
	$position_quals_1_mean[$i]=0; # Mean values of qualities of the correctly read 1-base-pairs of the sequences
	$position_quals_2_mean[$i]=0; # Mean values of qualities of the correctly read 2-base-pairs of the sequences
	$position_quals_3_mean[$i]=0; # Mean values of qualities of the correctly read 3-base-pairs of the sequences
	$position_quals_0_count[$i]=0; # Count of 0-base-pairs per position
	$position_quals_1_count[$i]=0; # Count of 1-base-pairs per position
	$position_quals_2_count[$i]=0; # Count of 2-base-pairs per position
	$position_quals_3_count[$i]=0; # Count of 3-base-pairs per position
	$position_bases_read[$i]=0; # Number of correctly read base-pairs per position
	$aux=$i+1;
	# Temporary files for the quality values of each position are created:
	my $cmd=`grep -v ">" $qual_name | awk '{print \$$aux}' > $dir/temp_dir_$base/temp_$i`;
}

for ($i=0;$i<$max_QV;$i++) {
	$per_seq_QVs[$i]=0; # array for storing the quality scores per sequence
	$per_pos_QVs[$i]=0; # array for storing the quality scores per position
}

my $seq_number=0; # Counter for the sequence number

# Each line of the qual file will be read and analyzed:
while ( $line_qual=<qual_file> ) {
	# The corresponding line in the csfasta file is read:
	$line_csfasta=<csfasta_file>;
	# It will be detected if the csfasta file is shorter than the qual file:
	if ( $line_csfasta eq "" ){
		print "The csfasta file $basename_csfasta has a shorter number of sequences than the qual file $basename_qual!!!!\n";
		exit;
	}
	# Creation of the vectors with qual values and base values:
	my @qual_values = split(' ',$line_qual);
	my @base_values = split(//,$line_csfasta);
	
	# Only the lines of the qual and csfasta files containing the bases and qual values will be processed:
	if ( substr($qual_values[0],0,1) ne ">" ) {
		
		$seq_number++;
		
		# Seqs with different number of base-pairs won't be considered:
		if ( (($#qual_values+1) != $csfasta_reads) || (($#base_values-1) != $csfasta_reads ) ) {
			print "File $basename_csfasta: the sequence ",$seq_number," has a different number of base-pairs and will not be considered !!!!\n";
		}
		else {
			# Initialization of sequence dependent variables:
			my $position=0;
			my $seqs_quals_mean=0;
			my $seqs_quals_mean_all=0;
			my $seq_bases_read=0;
	
			foreach my $val (@qual_values) {
				if ($seq_number == 1) {
					$position_quals_min[$position]=$val;
				}
				
				# The previous version of the program used the commented lines, so that no temporary disk storage was needed:
				#$seqs_quals_mean[$seq_number-1]=$seqs_quals_mean[$seq_number-1]+$val;
				#$seqs_quals_values[$position][$seq_number-1]=$val;
				$seqs_quals_mean_all=$seqs_quals_mean_all+$val;
	
				#print qstats_file "seq: $seq_number \t position: $position val: $val\t mean: $position_quals_mean[$position] \t max: $position_quals_max[$position] \t min: $position_quals_min[$position] \t mean0: $position_quals_0_mean[$position] \t base_value: $base_values[$position+1]\n";
				if ($csfasta_reads > $position) {
					
					# When the position/base-pair is correctly read, then the different scores are calculated:
					if ( $base_values[$position+1] ne "." ) {
						
						# Maximum and minimum values get stored
						if( $val > $position_quals_max[$position] ) {
							$position_quals_max[$position]=$val;
						}
						if ( $val < $position_quals_min[$position] ) {
							$position_quals_min[$position]=$val;
						}

						$seq_bases_read++;
						$position_bases_read[$position]++;
						$position_quals_mean[$position]=$position_quals_mean[$position]+$val;
						$seqs_quals_mean=$seqs_quals_mean+$val;
						
						# Histogram of QVs per position:
						$per_pos_QVs[int($val)]++;

						# if the base-pair corresponds to 0:
						if ( $base_values[$position+1] eq "0" ) {
							$position_quals_0_mean[$position]=$position_quals_0_mean[$position] + $val;
							$position_quals_0_count[$position]++;
						}
						# if the base-pair corresponds to 1:
						if ( $base_values[$position+1] eq "1" ) {
							$position_quals_1_mean[$position]=$position_quals_1_mean[$position]+$val;
							$position_quals_1_count[$position]++;
						}
						# if the base-pair corresponds to 2:
						if ( $base_values[$position+1] eq "2" ) {
							$position_quals_2_mean[$position]=$position_quals_2_mean[$position]+$val;
							$position_quals_2_count[$position]++;
						}
						# if the base-pair corresponds to 3:
						if ( $base_values[$position+1] eq "3" ) {
							$position_quals_3_mean[$position]=$position_quals_3_mean[$position]+$val;
							$position_quals_3_count[$position]++;
						}
					}
				}
				else {
					print "Inconsistent quality file (some sequences are longer than ",$csfasta_reads," bases!!!!\n";
					exit;
				}
				$position++;
			}
			# Histogram of QVs per sequence:
			my $seqs_mean;
			if ( $seq_bases_read == 0 ) {
				# In case that no base-pair could be identified as a color, the mean quality of the sequence will be -1:
				$seqs_mean=-1;
			} else {
				$seqs_mean=$seqs_quals_mean / $seq_bases_read;
				if ( $seqs_mean > 33 ) {
					print log_file "Seq Number $seq_number\n";
					print log_file $line_qual;
				}
			}
			if ( int($seqs_mean) > 33 ) {
				print log_file "Seq Number $seq_number\n";
				print log_file $line_qual;
			}
			if ( $seqs_mean ne -1) {
				$per_seq_QVs[int($seqs_mean)]++;
			} else {
				print log_file "Seq Number $seq_number\n";
				print log_file $line_qual;
			}
		}
	}
}

my $total_0_bases=0;
my $total_1_bases=0;
my $total_2_bases=0;
my $total_3_bases=0;

# Process of position temp files to obtain median values:
for ($i=0;$i<$csfasta_reads;$i++) {
	my $cmd=`sort -n $dir/temp_dir_$base/temp_$i > $dir/temp_dir_$base/temp_sorted_$i`;
	$median_line=int($qual_seqs/2);
	$percent20_line=int($qual_seqs*0.2);
	$percent80_line=int($qual_seqs*0.8);
	$position_quals_median[$i]=`head -$median_line $dir/temp_dir_$base/temp_sorted_$i | tail -1 2>&1`;
	$position_quals_median[$i]=~ s/\r|\n//g;
	$position_quals_percent20[$i]=`head -$percent20_line $dir/temp_dir_$base/temp_sorted_$i | tail -1 2>&1`;
	$position_quals_percent20[$i]=~ s/\r|\n//g;
	$position_quals_percent80[$i]=`head -$percent80_line $dir/temp_dir_$base/temp_sorted_$i | tail -1 2>&1`;
	$position_quals_percent80[$i]=~ s/\r|\n//g;
	$position_quals_mean_all[$i]=`awk '{sum=sum+\$1} END {print sum}' $dir/temp_dir_$base/temp_sorted_$i`;
	$position_quals_mean_all[$i]=~ s/\r|\n//g;
	$position_quals_mean_all[$i]=$position_quals_mean_all[$i]/$qual_seqs;
	$total_0_bases=$total_0_bases + $position_quals_0_count[$i];
	$total_1_bases=$total_1_bases + $position_quals_1_count[$i];
	$total_2_bases=$total_2_bases + $position_quals_2_count[$i];
	$total_3_bases=$total_3_bases + $position_quals_3_count[$i];
}

my $cmd=`rm -rf $dir/temp_dir_$base`;

if ( $reduce_seqs == 1) {
	$cmd=`if [ -e $dir1/10Mseqs_$basename_csfasta$ext1 ]; then rm $dir1/10Mseqs_$basename_csfasta$ext1; fi`;
	$cmd=`if [ -e $dir2/10Mseqs_$basename_qual$ext2 ]; then rm $dir2/10Mseqs_$basename_qual$ext2; fi`;
}

# Results get printed to output_file:
print qstats_file "base position\tmean of all quality values\tbases correctly read\tmean of read quality values\tmaximum quality value\tminimum quality value\tmedian\t0-color calls mean QV\t1-color calls mean QV\t2-color calls mean QV\t3-color calls mean QV\t0-color calls (%)\t1-color calls (%)\t2-color calls (%)\t3-color calls (%)\t20-percentil\t80-percentil\n";
for ($i=0;$i<$csfasta_reads;$i++) {	
	if ($position_bases_read[$i] == 0) {
		$position_quals_mean[$i]=0;
		$position_quals_0_mean[$i]=0;
		$position_quals_1_mean[$i]=0;
		$position_quals_2_mean[$i]=0;
		$position_quals_3_mean[$i]=0;
		$aux_pos_quals_0_count=0;
		$aux_pos_quals_1_count=0;
		$aux_pos_quals_2_count=0;
		$aux_pos_quals_3_count=0;
	}
	else {
		$position_quals_mean[$i]=$position_quals_mean[$i] / $position_bases_read[$i];
		$position_quals_0_mean[$i]=$position_quals_0_mean[$i] / $position_bases_read[$i];
		$position_quals_1_mean[$i]=$position_quals_1_mean[$i] / $position_bases_read[$i];
		$position_quals_2_mean[$i]=$position_quals_2_mean[$i] / $position_bases_read[$i];
		$position_quals_3_mean[$i]=$position_quals_3_mean[$i] / $position_bases_read[$i];
		$aux_pos_quals_0_count=$position_quals_0_count[$i]/$position_bases_read[$i]*100;
		$aux_pos_quals_1_count=$position_quals_1_count[$i]/$position_bases_read[$i]*100;
		$aux_pos_quals_2_count=$position_quals_2_count[$i]/$position_bases_read[$i]*100;
		$aux_pos_quals_3_count=$position_quals_3_count[$i]/$position_bases_read[$i]*100;
	}

	$line_file=$i."\t".$position_quals_mean_all[$i]."\t".$position_bases_read[$i]."\t".$position_quals_mean[$i]."\t".$position_quals_max[$i]."\t".$position_quals_min[$i]."\t".$position_quals_median[$i]."\t".$position_quals_0_mean[$i]."\t".$position_quals_1_mean[$i]."\t".$position_quals_2_mean[$i]."\t".$position_quals_3_mean[$i]."\t".$aux_pos_quals_0_count . "\t".$aux_pos_quals_1_count . "\t".$aux_pos_quals_2_count . "\t".$aux_pos_quals_3_count . "\t".$position_quals_percent20[$i]."\t".$position_quals_percent80[$i]."\n";

	# Dots get changed into commas, for compatibility with spanish-excel:
	$line_file=~ s/\./,/g;
	print qstats_file $line_file;
}

# Header of the sequence-related values will be printed:
print qstats_file "\nQV\tper position QV scores\tsequences count per mean sequence quality\n";
for ($i=0;$i<$max_QV;$i++){
	print qstats_file $i,"\t",$per_pos_QVs[$i],"\t",$per_seq_QVs[$i],"\n";
}
#print seqstats_file "sequence\tmean of quality values\tbases correctly read\tmean of quality values of bases correctly read\n";
#print seqstats_file $seq_number,"\t",$seqs_quals_mean_all / $csfasta_reads,"\t",$seq_bases_read,"\t",$seqs_quals_mean / $seq_bases_read,"\n";

print "Number of sequences in csfasta file $basename_csfasta: \t$csfasta_seqs, in qual file $basename_qual: \t$qual_seqs \tNumber of reads per sequence: \t$csfasta_reads \n";

# Opened files get closed
close(csfasta_file);
close(qual_file);
close(qstats_file);

