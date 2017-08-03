#!/usr/bin/perl -w
# 2011/07/13 - árbol sw.
use File::ReadBackwards;
use File::Basename;

sub usage {
	print "\nHow to use: sg_extrat_qual_seqs.pl <fasta_seqs_file> <qual_file>\n\n";
	print "Description: this perl function extract the qualities corresponding to the sequences included in the input fasta file from the input qual file. Please notice that both files must be sorted by read tags.\n\n";
	print "INPUT: <fasta_seqs_file> is the fasta (or csfasta) file containing the desired sequences whose qualities have to be extracted. <qual_file> is the file which contains all the qualitites relative to the fasta file.\n\n";
	print "OUTPUT: <qual_seqs_file> is the file containing only the qualities corresponding to the sequences in the fasta file.\n";
	print "Please notice that the ouput files will be placed in the same directory than the input files, and with the same file extension\n";
	exit(1);
}

if(scalar(@ARGV) != 2){
    usage();
}

# input and output files get opened
open(fasta_seqs_file, "<", $ARGV[0]) || die("Could not open file!\n"); #$#ARGV]);
open(qual_file, "<", $ARGV[1]) || die("Could not open file!\n");

($base, $dir) = fileparse($ARGV[0], qr/\.[^.]*/); #($base, $dir, $ext) = fileparse($ARGV[0], qr/\.[^.]*/);
my $qual_seqs_name=$dir . $base . "_QV.qual";
if ( -e $qual_seqs_name ) {
	print "File ",$qual_seqs_name," already exists! Do you want to overwrite it (y/n)?";
	$answer=<STDIN>;
	chomp ($answer);
	if ($answer eq "y"){
		open(qual_seqs_file, ">", $dir . $base . "_QV.qual") || die("Could not open file!\n");
	}
	else {
		print "Execution aborted!!";
		exit 1;
	}
}
else {
	open(qual_seqs_file, ">", $dir . $base . "_QV.qual") || die("Could not open file!\n");
}

# Odd lines contain the read tags, while even lines contain the sequences (colorspace-dinucleotides, nucleotides or quality values)
# Read tags have the following structure: >34_385_6571_F3
# The last 3 characters need to be removed for the subsequent string comparison, when the comparison is performed between mate-pair or pair-end reads.

# A read tag from file_1 will be read, and string-compared with the current read tag from file_2. If the read tag from file_1 is "smaller", then it will be removed (together with its sequence), and the next read tag from file_1 will be retrieved. If it's "greater", the read tag from file_2 (together with its sequence) will be removed, and the next read tag from file_2 will be retrieved. If both are equal, then the read tags from both files (together with its sequences) get stored in their respective files.

$lines_processed=0;
$flag_file_to_read=3; # value 3 for reading both files, 2 for reading file_2, 1 for reading file_1, 0 for reading no file, -1 for finishing (once the end of one of the files has been reached)

while ( $flag_file_to_read > 0 ) {
	
	# The next tag and sequence from fasta file will be read
	if ($flag_file_to_read==1 || $flag_file_to_read==3) {
		$fasta_line_tag=<fasta_seqs_file>;

		# Once a new line of the file has been read, it will be splitted in fields separated by '_'. It will be detected if the line is a "read tag" (it begins with ">").
		if ( $fasta_line_tag ) {
			@fasta_tag_aux=(split /,/, $fasta_line_tag);
			@fasta_tag=(split /_/, $fasta_tag_aux[0]);
			if (substr( $fasta_line_tag, 0 ,1) eq ">") { #$fasta_tag[0], 0, 1) eq ">") {
				$fasta_tag[0]=~ s/>//g;
				# Sequence associated to the tag will be retrieved:
				$fasta_line_seq=<fasta_seqs_file>;
				if ( $fasta_line_seq ) {
					if (substr($fasta_line_seq, 0 ,1) eq ">") { #$fasta_seq[0], 0, 1) eq ">") {
						print "ERROR: Problem in fasta file $ARGV[0]!! No sequence associated to the tag was found!!\n";
						exit 1;
					}
					else {
						if ( $flag_file_to_read==3 ) {
							$flag_file_to_read=2;
						}
						else{
							$flag_file_to_read=0;
						}
					}
				}
			}
			else {
				print "WARNING: line $fasta_line_tag from $ARGV[0] won't be considered!!! \n";
			}
		}
		else {
			$flag_file_to_read=-1;
		}
	}
	
	# The next tag and sequence from qual file will be retrieved
	if ($flag_file_to_read==2 || $flag_file_to_read==3) {
		$qual_line_tag=<qual_file>;

		# It will be displayed the millions of lines processed
		$lines_processed=$lines_processed+1;
		if ( ($lines_processed % 10000000) == 0 ) {
			print "Info: Lines processed from Qual file $ARGV[1]: ", $lines_processed, "\n";
		}

		# Once a new line of the file has been read, it will be splitted in fields separated by '_'. It will be detected if the line is a "read tag" (it begins with ">").
		if ( $qual_line_tag ) {
			@qual_tag=(split /_/, $qual_line_tag);
			if (substr($qual_line_tag, 0, 1) eq ">") { #$qual_tag[0], 0, 1) eq ">") {
				$qual_tag[0]=~ s/>//g;
				# Sequence associated to the tag will be retrieved:
				$qual_line_seq=<qual_file>;
				if ( $qual_line_seq ) {
					if (substr($qual_line_seq, 0 ,1) eq ">") { #$fasta_seq[0], 0, 1) eq ">") {
						print "ERROR: Problem in qual file $ARGV[1]!! No sequence associated to the tag was found!!\n";
						exit 1;
					}
					else {
						if ( $flag_file_to_read==3 ) {
							$flag_file_to_read=1;
						}
						else {
							$flag_file_to_read=0;
						}
					}
				}
			}
			else {
				print "WARNING: line $qual_line_tag won't be considered from qual file $ARGV[1]!! \n";
			}
		}
		else {
			print "WARNING: The qual file $ARGV[1] has reached its end!!\n";
			$flag_file_to_read=-1;
		}
	}
	
	# If no file has reached its end, then the read tags will be compaired:
	if ( $flag_file_to_read == 0 ) {
		$i=0;
		# The different numbers from the read tags get compaired one by one:
		while (( $fasta_tag[$i] eq $qual_tag[$i] ) && ( $i < (scalar @fasta_tag-1) )) {
			$i++;
		}
		if ( $i < (scalar @fasta_tag-1) ) {
			if ( $fasta_tag[$i] < $qual_tag[$i] ) {
				print "WARNING: The tag ID $fasta_tag_aux[0] wasn't found in the qual file $ARGV[1]!!!\n";
				$flag_file_to_read=1;
			}
			else {
				$flag_file_to_read=2;
			}
		}
		else {
			# Read tag match! Read tag and sequence will be written to the qual_seqs_file
			print qual_seqs_file $fasta_line_tag;
			print qual_seqs_file $qual_line_seq;
			$flag_file_to_read=3;
		}
	}
}
		
# Opened files get closed
close(fasta_seqs_file);
close(qual_seqs_file);
close(qual_file);
