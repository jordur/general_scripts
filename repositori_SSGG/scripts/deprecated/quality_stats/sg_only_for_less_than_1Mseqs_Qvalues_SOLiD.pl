#!/usr/bin/perl -w
# 2011/07/13 - árbol sw.
use File::ReadBackwards;
use File::Basename;

sub usage {
	print "\nHow to use: sg_Qvalues_SOLiD.pl <csfasta> <qual>\n\n";
	print "Description: this perl function calculates the quality values of a .csfasta file and its associated .qual file. For computing reasons only a maximum of 10M sequences should be taken from the original .csfasta and .qual files.\n\n";
	print "INPUT: <csfasta> and <qual> are the .csfasta and .qual files from a given strand.\n\n";
	print "OUTPUT: <base_file_name_Qstats.tsv> and <base_file_name_Seqstats.tsv> are two tab-separated files containing the stats of each read position (_Qstats.tsv) and the stats of each sequence (_Seqstats.tsv)\n";
	print "Please notice that the ouput files will be placed in the same directory than the input files, and with the same file extension\n";
	exit(1);
}

sub system_bash {
	my @args = ( "bash", "-c", shift );
	system(@args);
}

if(scalar(@ARGV) != 2){
    usage();
}

# input and output files get opened
open(csfasta_file, "<", $ARGV[0]) || die("Could not open file!"); #$#ARGV]);
open(qual_file, "<", $ARGV[1]) || die("Could not open file!");
my $csfasta_name=$ARGV[0];
my $qual_name=$ARGV[1];

($base, $dir, $ext) = fileparse($ARGV[0], qr/\.[^.]*/);
open(qstats_file, ">", $dir . $base . "_Qstats.tsv" ) || die("Could not open file!");
open(seqstats_file, ">", $dir . $base . "_Seqstats.tsv" ) || die("Could not open file!");

# The number of reads in the csfasta and qual files will be obtained and compared, so as to determine if both numbers of reads match together:

my $csfasta_seqs=`grep -c ">" $csfasta_name 2>&1`;
$csfasta_seqs=~ s/\r|\n//g;
my $qual_seqs=`grep -c ">" $qual_name 2>&1`;
$qual_seqs=~ s/\r|\n//g;
my $csfasta_reads=`head -2 $csfasta_name | tail -1 | wc | awk '{print \$3-2}' 2>&1`;
$csfasta_reads=~ s/\r|\n//g;
print "Number of sequences in csfasta file: $csfasta_seqs, in qual file: $qual_seqs \nNumber of reads per sequence: $csfasta_reads \n";

if ( $csfasta_seqs ne $qual_seqs ){
	print "The .csfasta and .qual number of sequences read are different!!! The programm will exit!!!";
	exit;
}

# Variables definition and initialization:
for ($i=0;$i<$qual_seqs;$i++) {
	$seqs_quals_mean[$i]=0;
}

for ($i=0;$i<$csfasta_reads;$i++) {
	$position_quals_mean[$i]=0;
	$position_quals_max[$i]=0;
	$position_quals_min[$i]=0;
	$position_quals_0_mean[$i]=0;
	$position_quals_1_mean[$i]=0;
	$position_quals_2_mean[$i]=0;
	$position_quals_3_mean[$i]=0;
	for ($j=0;$j<$qual_seqs;$j++) {
		$seqs_quals_values[$i][$j]=0;
	}
}

my $seq_number=0;

while ( $line_qual=<qual_file> ) {
	$line_csfasta=<csfasta_file>;
	if ( $line_csfasta eq "" ){
		print "The csfasta file has a shorter number of sequences than the qual file!!!!";
		exit;
	}
	my @qual_values = split(' ',$line_qual);
	my @base_values = split(//,$line_csfasta);
	if ( substr($qual_values[0],0,1) ne ">" ) {
		$seq_number++;
		my $position=0;

		foreach my $val (@qual_values) {
			if ($seq_number == 1) {
				$position_quals_min[$position]=$val;
			}
			
			$seqs_quals_mean[$seq_number-1]=$seqs_quals_mean[$seq_number-1]+$val;
			$seqs_quals_values[$position][$seq_number-1]=$val;

			if ($csfasta_reads > $position) {
				$position_quals_mean[$position]=$position_quals_mean[$position]+$val;
				if ( $val > $position_quals_max[$position] ) {
					$position_quals_max[$position]=$val;
				}
				if ( $val < $position_quals_min[$position] ) {
					 $position_quals_min[$position]=$val;
				}

				# if the base-pair corresponds to 0:
				if ( $base_values[$position+1] eq "0" ) {
					$position_quals_0_mean[$position]=$position_quals_0_mean[$position] + $val;
				}
				# if the base-pair corresponds to 1:
				if ( $base_values[$position+1] eq "1" ) {
					$position_quals_1_mean[$position]=$position_quals_1_mean[$position]+$val;
				}
				# if the base-pair corresponds to 2:
				if ( $base_values[$position+1] eq "2" ) {
					$position_quals_2_mean[$position]=$position_quals_2_mean[$position]+$val;
				}
				# if the base-pair corresponds to 3:
				if ( $base_values[$position+1] eq "3" ) {
					$position_quals_3_mean[$position]=$position_quals_3_mean[$position]+$val;
				}				
			}
			else {
				print "Inconsistent quality file (some sequences are longer than ",$csfasta_reads," bases!!!!";
				exit;
			}
			$position++;
		}
	}
}

# Results get printed to output_file:
print qstats_file "base position\tmean of quality values\tmaximum quality value\tminimum quality value\tmedian\t0-basepairs mean QV\t1-basepairs mean QV\t2-basepairs mean QV\t3-basepairs mean QV\n";
for ($i=0;$i<$csfasta_reads;$i++) {	
	# The array of quality values per position and sequence has to be sorted so that the "median" value can be obtained:
	$qual_values_unsorted=\@{$seqs_quals_values[$i]};
	@qual_values_sorted=sort { $a <=> $b } @{$qual_values_unsorted};
	$position_quals_mean[$i]=$position_quals_mean[$i] / $qual_seqs;
	$position_quals_0_mean[$i]=$position_quals_0_mean[$i] / $qual_seqs;
	$position_quals_1_mean[$i]=$position_quals_1_mean[$i] / $qual_seqs;
	$position_quals_2_mean[$i]=$position_quals_2_mean[$i] / $qual_seqs;
	$position_quals_3_mean[$i]=$position_quals_3_mean[$i] / $qual_seqs;
	print qstats_file $i,"\t",$position_quals_mean[$i],"\t",$position_quals_max[$i],"\t",$position_quals_min[$i],"\t",$qual_values_sorted[$qual_seqs/2],"\t",$position_quals_0_mean[$i],"\t",$position_quals_1_mean[$i],"\t",$position_quals_2_mean[$i],"\t",$position_quals_3_mean[$i],"\n";
}

print seqstats_file "sequence\tmean of quality values\n";
$i=0;
foreach my $val (@seqs_quals_mean) {
	print seqstats_file $i,"\t",$val/$csfasta_reads,"\n";
	$i++;
}

# Opened files get closed
close(csfasta_file);
close(qual_file);
close(qstats_file);
close(seqstats_file);

