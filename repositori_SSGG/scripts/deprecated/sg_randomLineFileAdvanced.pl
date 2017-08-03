#!/usr/bin/perl -w
# 2011/09/16 - árbol sw.
#use strict;
use File::ReadBackwards;
use File::Basename;

sub usage {
	print"\n\nINPUTS:\n\t\tI:File in .csfasta format\n\t\tII:File in .qual format\n\t\tIII:Number of random lines\n\t\tIV:output file prefix\n\t\tV (optional):File 2 in .csfasta format\n\t\tVI (optional):File 2 in .qual format\n\nNote 1: the resultant files will be output to the directories where the input files are in, adding the given prefix at the beginning ot their names.\n\nNote 2: parameters V and VI need to be added in case of mate-pair reads or pair-end reads.\n\n";
    exit(1);
}
if( (scalar(@ARGV) != 6) && (scalar(@ARGV) != 4) ){
    usage();
}

# variables definitions
my @numberLine=();

my $rand=0;
my $entry=0;
my $count=0;
my $countLine=0;
my $exit=0;

my $line="";
my $line2="";
my $line3="";
my $line4="";

my $CS_OUT;
my $QV_OUT;
my $CS_OUT2;
my $QV_OUT2;
my $FILE;
my $FILE2;
my $FILE3;
my $FILE4;

my $input=$ARGV[0];
my $inputQval=$ARGV[1];
my $numberRand=$ARGV[2];

# Names, paths and extensions of the firsts csfasta and qual files are obtained:
($fasta1_base, $fasta1_dir, $fasta1_ext) = fileparse($ARGV[0], qr/\.[^.]*/);
($qual1_base, $qual1_dir, $qual1_ext) = fileparse($ARGV[1], qr/\.[^.]*/);

# Output files names for firsts csfasta and qual files are obtained:
my $outputCsfasta=$fasta1_dir . $ARGV[3] . $fasta1_base . $fasta1_ext;
my $outputQval=$qual1_dir . $ARGV[3] . $qual1_base . $qual1_ext;

# The output files get opened:
open($CS_OUT,'>',$outputCsfasta) || die("Could not open file $outputCsfasta!");
open($QV_OUT,'>',$outputQval) || die("Could not open file $outputQval!");

# The number of reads in the csfasta and qual files will be obtained and compared, so as to determine if both numbers of reads match together:
#my $csfasta_seqs=`grep -c ">" $input 2>&1`;
#$csfasta_seqs=~ s/\r|\n//g; # Number of sequences in the csfasta file
#my $qual_seqs=`grep -c ">" $inputQval 2>&1`;
#$qual_seqs=~ s/\r|\n//g; # Number of sequences in the qual file
#my $csfasta2_seqs=`grep -c ">" $input2 2>&1`;
#$csfasta2_seqs=~ s/\r|\n//g; # Number of sequences in the csfasta file
#my $qual2_seqs=`grep -c ">" $inputQval2 2>&1`;
#$qual2_seqs=~ s/\r|\n//g; # Number of sequences in the qual file

#print "Number of sequences in csfasta file 1: $csfasta_seqs, in qual file 1: $qual_seqs \n";
#print "Number of sequences in csfasta file 2: $csfasta2_seqs, in qual file 2: $qual2_seqs \n";

#if ( $csfasta_seqs <> $csfasta2_seqs || $csfasta_seqs <> $qual_seqs || $csfasta_seqs <> $qual2_seqs ) {
#	print "The numbers of sequences in csfasta & qual files are different!!!"
#	exit (1);
#}

# generate the vector of random lines to get from files (0 for lines which will be left out, 1 for lines which will be included in the new files)

open(FILE,'<',$input);

while($line=<FILE>) {
	chomp($line);
	
	if($count<$numberRand) {
	
		if($line=~/>/) {
			$rand=int(rand(2));

			if($rand==1) {
				print $CS_OUT "$line\n";
				$numberLine[$countLine]=1;
				$entry=1;
			}
			else {
				$numberLine[$countLine]=0;
				$entry=0;
			}
		}
		else {
			if($entry==1) {
				print $CS_OUT "$line\n";
				$entry=0;
				$count=$count+1;
			}
		}
	}
	else {
		last;
	}
	$countLine=$countLine+1;
}

close(FILE);

$countLine=0;
$entry=0;
$exit=0;

open(FILE2,'<',$inputQval) || die;
while($line2=<FILE2>)
{
	chomp($line2);
	if($exit<$numberRand) {
		if($line2=~/>/) {
			if($numberLine[$countLine]==1) {
				print $QV_OUT "$line2\n";
				$entry=1;
			}
			else {
				$entry=0;
			}
		}		
		else {
			if($entry==1) {
				print $QV_OUT "$line2\n";
				$exit=$exit+1;
				$entry=0;
			}
		}
	}
	else
	{
		last;
	}
	$countLine=$countLine+1;
}
close(FILE2);	


# if the seconds csfasta and qual files are also entered as parameters, then they get processed:
if ( scalar(@ARGV) == 6 ){
	
	my $input2=$ARGV[4];
	my $inputQval2=$ARGV[5];

	# Names, paths and extensions of the seconds csfasta and qual files are obtained:
	($fasta2_base, $fasta2_dir, $fasta2_ext) = fileparse($ARGV[4], qr/\.[^.]*/);
	($qual2_base, $qual2_dir, $qual2_ext) = fileparse($ARGV[5], qr/\.[^.]*/);
	my $outputCsfasta2=$fasta2_dir . $ARGV[3] . $fasta2_base . $fasta2_ext;
	my $outputQval2=$qual2_dir . $ARGV[3] . $qual2_base . $qual2_ext;
	open($CS_OUT2,'>',$outputCsfasta2) || die("Could not open file $outputCsfasta2!");
	open($QV_OUT2,'>',$outputQval2) || die("Could not open file $outputQval2!");

	$countLine=0;
	$entry=0;
	$exit=0;

	open(FILE3,'<',$input2) || die;
	while($line3=<FILE3>) {
		chomp($line3);

		if($exit<$numberRand) {
			if($line3=~/>/) {
				if($numberLine[$countLine]==1) {
					print $CS_OUT2 "$line3\n";
					$entry=1;
				}
				else {
					$entry=0;
				}
			}
			else {
				if($entry==1) {
					print $CS_OUT2 "$line3\n";
					$exit=$exit+1;
					$entry=0;
				}
			}
		}
		else {
			last;
		}
		$countLine=$countLine+1;
	}

	close(FILE3);

	$countLine=0;
	$entry=0;
	$exit=0;

	open(FILE4,'<',$inputQval2) || die;
	while($line4=<FILE4>) {
		chomp($line4);
		if($exit<$numberRand) {
			if($line4=~/>/) {
				if($numberLine[$countLine] == 1) {
					print $QV_OUT2 "$line4\n";
                	$entry=1;
				}
				else {
					$entry=0;
				}
			}
			else {
				if($entry==1) {
					print $QV_OUT2 "$line4\n";
					$exit=$exit+1;
					$entry=0;
				}
			}
		}
		else
		{
			last;
		}
                $countLine=$countLine+1;
	}

	close(FILE4);
}
