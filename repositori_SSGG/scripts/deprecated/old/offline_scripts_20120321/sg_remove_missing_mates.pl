#!/usr/bin/perl -w
# 2011/07/13 - árbol sw.
use File::ReadBackwards;
use File::Basename;

sub usage {
        print "\nHow to use: sg_remove_missing_mates.pl <file_1> <file_2>\n\n";
	print "Description: this perl function removes the reads from files 1 and 2, whose mate-pair reads are missing in the other file. Please notice that both files must be sorted by read tags. It is also important to carry out this script by hand on all the sample files, that is, on the .fasta/.csfasta and .qual associated files.\n\n";
        print "INPUT: <file_1> and <file_2> are both .fasta, .csfasta or .qual, with the usual file structure >read_tag followed by sequence or quality values.\n\n";
        print "OUTPUT: <file_1_without_mates> and <file_2_without_mates>\n";
	print "Please notice that the ouput files will be placed in the same directory than the input files, and with the same file extension\n";
        exit(1);
}

if(scalar(@ARGV) != 2){
    usage();
}

# input and output files get opened
open(file_1, "<", $ARGV[0]) || die("Could not open file!"); #$#ARGV]);
open(file_2, "<", $ARGV[1]) || die("Could not open file!");

($base, $dir, $ext) = fileparse($ARGV[0], qr/\.[^.]*/);
open(output_1, ">", $dir . $base . "_only_mates" . $ext) || die("Could not open file!");
($base, $dir, $ext) = fileparse($ARGV[1], qr/\.[^.]*/);
open(output_2, ">", $dir . $base . "_only_mates" . $ext) || die("Could not open file!");

# Odd lines contain the read tags, while even lines contain the sequences (colorspace-dinucleotides, nucleotides or quality values)
# Read tags have the following structure: >34_385_6571_F3
# The last 3 characters need to be removed for the subsequent string comparison.

# A read tag from file_1 will be read, and string-compared with the current read tag from file_2. If the read tag from file_1 is "smaller", then it will be removed (together with its sequence), and the next read tag from file_1 will be retrieved. If it's "greater", the read tag from file_2 (together with its sequence) will be removed, and the next read tag from file_2 will be retrieved. If both are equal, then the read tags from both files (together with its sequences) get stored in their respective files.

$lines_processed=0;
$flag_file_to_read=3; # value 3 for reading both files, 2 for reading file_2, 1 for reading file_1 and 0 for reading no file

while ( $flag_file_to_read ) {

	# The next read with read tag from file_1 will be searched
	while ( $flag_file_to_read==1 || $flag_file_to_read==3) {
		$line_file_1=<file_1>;

		# It will be displayed the millions of lines processed
		$lines_processed=$lines_processed+1;
	        if ( ($lines_processed % 1000000) == 0 ) {
               		print "Lines processed: ", $lines_processed, "\n";
	        }

		# Once a new line of the file has been read, it will be splitted in fields separated by '_'. It will be detected if the line is a "read tag" (it begins with ">").
		if ( $line_file_1 ) {
       	        	@read_tag_1=(split /_/, $line_file_1);
                        if (substr($read_tag_1[0], 0, 1) eq ">") {
				$read_tag_1[0]=~ s/>//g;
                                if ( $flag_file_to_read==1 ) {
                                        $flag_file_to_read=4;
                                }
                                else {
                                        $flag_file_to_read=2;
                                }
                        }
        	}
        	else {
       	        	$flag_file_to_read=0;
        	}
	}
	
	# The next read with read tag from file_2 will be searched
	while ( $flag_file_to_read==2 || $flag_file_to_read==3) {
		$line_file_2=<file_2>;

		# Once a new line of the file has been read, it will be splitted in fields separated by '_'. It will be detected if the line is a "read tag" (it begins with ">").
		if ( $line_file_2 ) {
			@read_tag_2=(split /_/, $line_file_2);
			#$read_tag_2=substr($line_file_2, 0, - 4);
			if (substr($read_tag_2[0], 0, 1) eq ">") {
				$read_tag_2[0]=~ s/>//g;
				if ( $flag_file_to_read==2 ) {
					$flag_file_to_read=4;
				}
				else {
					$flag_file_to_read=1;
				}
			}
		}
		else {
       	                $flag_file_to_read=0;
               	}
	}
	
	# If no file has reached its end, then the read tags will be compaired:
	if ( $flag_file_to_read ) {
		$i=0;
		# The different numbers from the read tags get compaired one by one:
		while (( $read_tag_1[$i] eq $read_tag_2[$i] ) && ( $i < (scalar @read_tag_1-1) )) {
			$i++;
		}
		if ( $i < (scalar @read_tag_1-1) ) {
			if ( $read_tag_1[$i] < $read_tag_2[$i] ) {
				$flag_file_to_read=1;
			}
			else {
				$flag_file_to_read=2;
			}
		}
		else {
			# Read tag match! Read tags and sequences will be stored in both files
			
			$aux_1=<file_1>;
			$aux_2=<file_2>;
			if ( $aux_1 && $aux_2 ){
				if ( (substr($aux_1, 0, 1) ne ">") && (substr($aux_2, 0, 1) ne ">") ) {
					print output_1 $line_file_1;
					print output_1 $aux_1;
					print output_2 $line_file_2;
					print output_2 $aux_2;
					$flag_file_to_read=3;
				}
				elsif ( substr($aux_1, 0, 1) eq ">" ) {
					$flag_file_to_read=1;
				}
				elsif ( substr($aux_2, 0, 1) ne ">" ){
					$flag_file_to_read=2;
				}
				else {
					$flag_file_to_read=3;
				}
			}
			else {
				$flag_file_to_read=0;
			}
		}
	}
}

# Opened files get closed
close(file_1);
close(file_2);
close(output_1);
close(output_2);

