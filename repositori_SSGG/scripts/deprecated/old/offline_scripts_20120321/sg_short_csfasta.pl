#!/usr/bin/perl -w
# 2011/10/14 - árbol sw.

use File::ReadBackwards;
use File::Basename;
use lib '/share/apps/bioperl';
use SOLiD_Reads;

sub usage {
	print "\nHow to use: sg_short_SOLiD_fasta.pl <fasta_seqs_file>\n\n";
	print "Description: this perl function shorts regarding the TAG_IDs a csfasta/qual file.\n\n";
	print "INPUT:\n";
	print "<fasta/qual file> is the csfasta/qual file containing the desired sequences/qualities to get shorted after their TAG_IDs.\n";
	print "<parallel_threads> is the number of threads that will run in parallel\n\n";
	print "OUTPUT: <fasta/qual shorted_file> file containing the shorted TAG_IDs and sequences/qualities corresponding to the original fasta/qual file.\n";
	print "Please notice that the ouput files will be placed in the same directory than the input files, and with the same file extension\n";
	exit(1);
}


sub next_read {
# retrieves the next read (tag_id and sequence/qual values) from the qual/fasta file
	my ( $file, $filehandle ) = @_;

	my $tag=<$filehandle>;
	# Once a new line of the file has been read, it will be splitted in fields separated by '_'. It will be detected if the line is a "read tag" (it begins with ">").
	if ( $tag ) {
		if (substr( $tag, 0 ,1) eq ">") {
			# Sequence associated to the tag will be retrieved:
			my $seq=<$filehandle>;
			if ( $seq ) {
				if (substr($seq, 0 ,1) eq ">") {
					print "ERROR: Problem in fasta file $file!! No sequence associated to the tag was found!!\n";
					exit 2;
				} 
				else {
					# If all was ok, then the Reads_SOLiD object will be returned:
					my $read=SOLiD_Reads::create_from_tag_seq($tag, $seq);
					return (2,$read);
				}
			} 
			else {
				print "WARNING: Last read doesn't contain a sequence associated to the tag in fasta file $file!!\n";
				return (1,undef);
			}
		} 
		else {
			print "WARNING: Incoherent line $tag in file $file!! Should be an ID tag!!\n";
			return (1,undef);
		}
	} 
	else {
		print "Info: Last line of file $file was reached.\n";
		return (0,undef);
	}
}


sub short_reads {
# shorts "read" into the "reads" array
# grabs N reads from file and gets them shorted into an array of SOLiD_Reads objetcs
	# Definitions
	my ( $read, @reads )= @_;

	if ( $read->reads_compare( $reads[0] ) == -1 ) {
		unshift @reads, $read;
	}
	elsif ( $read->reads_compare ( $reads[$#reads] ) == 1 ) {
		push @reads, $read;
	}
	else {
		push @reads, $read;
		my $flag_minor=0;
		$j=$#reads;
		# The new element gets placed at the end of the array, and will be swapped with the previous elements until the next previous element is determined to be smaller
		while ( $flag_minor==0 && $j > 0 ) {
			if ( $reads[$j]->reads_compare ( $reads[$j-1] ) == -1 ) {
				($reads[$j], $reads[$j-1])=($reads[$j-1], $reads[$j]);
				$j--;
			}
			else {
				$flag_minor=1;
			}
		}
	}

	return @reads;
}


# ---------------------------------------
# --------------- main ------------------
# ---------------------------------------

if(scalar(@ARGV) != 1){
    usage();
}

# Definitions
$fasta_file_name=$ARGV[0];
#$parallel_threads=$ARGV[1];

my $max_seqs=50000; # Should be something like 200000~50000000 (40MB~10GB of RAM)!! # Maximum number of sequences to be read and stored in RAM. Each sequence will contain the bead TAG_ID (panel_xpixel_ypixel_tagtype) and the sequence/quality values
my $max_files=100; # Should be something like 100~1000!! # Maximum number of files to get opened simultaneously. It has to be considered that a read (TAG_ID and sequence) from each file will be kept in RAM storage

my $temp_file; # name of the temporary subfile to be created

my @filereads; # array for storing the reads from files to be shorted
my @filehandles; # array for storing the handles to the temporary files

# input and output files get opened
open(fasta_file, "<", $fasta_file_name) || die("Could not open file!\n"); #$#ARGV]);

($base, $dir, $ext) = fileparse($ARGV[0], qr/\.[^.]*/); #($base, $dir, $ext) = fileparse($ARGV[0], qr/\.[^.]*/);
my $shorted_file_name=$dir . $base . "_shorted" . $ext;

# Odd lines contain the read tags (in fasta/csfasta/qual files), while even lines contain the sequences/qualities (colorspace-dinucleotides, nucleotides or quality values)
# Read tags have the following structure: >34_385_6571_F3 (bead TAG_ID: panel_xpixel_ypixel_tagtype)
# The last 3 characters need to be removed for the subsequent string comparison, when the comparison is performed between mate-pair or pair-end reads.

my $j=0; # auxiliar variable for counting the number of auxiliar files created
my $lines=0; # variable containing the number of lines processed from input file
do {
	do {
		my $i=0; # auxiliar variable for counting the number of reads added to a file
		my @reads; # array for storing the reads that get shorted
		do {
			do {
				($status_general,$read) = next_read($fasta_file_name, *fasta_file);
			} while ( $status_general == 1 );
		
			if ( $status_general != 0 ) {
				$i++;
				$lines++;
				if ( $i == 1) {
					push @reads, $read;
				}
				else {
					@reads=short_reads($read,@reads);
				}
				# Info report: each 10k lines a log will be output
				if ( ($lines % 10000) == 0 ) {
					print "Info: Lines processed from file $fasta_file_name: ", $lines, "\n";
				}
			}
		} while ( $i<$max_seqs && $status_general != 0 );
		if ( $i > 0 ){
			$j++;
			$temp_file=$dir . $base . "_temp" . $j;
	
			local *FILE;
			open(FILE,">",$temp_file) || die;
			#push the typeglobe to the end of the filehandles array
			push @filehandles, *FILE;
			push @filenames, $temp_file;
			
			foreach (@reads) {
				print {$filehandles[$#filehandles]} $_->return_file_lines;
			}

			#push @filereads, $reads[0];
			close $filehandles[$#filehandles];
		}
	} while ( $j<$max_files && $status_general != 0 );
	
	# The output files get opened for reading and the firsts reads stored in the $filreads array
	for (my $k=0; $k<$#filehandles+1; $k++) {
		open($filehandles[$k], "<", $filenames[$k]);
		do {
			($status, $filereads[$k])=next_read($filenames[$k],$filehandles[$k]);
		} while ($status==1);
		if ($status==0) {
			close $filehandles[$k];
            # the file has to be deleted
            unlink $filenames[$k];
            splice @filehandles, $k, 1;
            splice @filenames, $k, 1;
            splice @filereads, $k, 1;
		}
	}
	
	# The auxiliar temp file get opened
	local *FILE;
	$aux_file=*FILE;
	$temp_file=$dir . $base . "_temp";
	open (FILE, ">", $temp_file);
	
	# The pre-shorted files will get shorted into the results shorted file:
	my $min_read;
	my $min;
	while ($#filehandles > 0) {
		$min_read=$filereads[0];
		$min=0;
		for ($k=1; $k<$#filehandles+1; $k++) {
			if ($min_read->reads_compare($filereads[$k]) == 1) {
				$min_read=$filereads[$k];
				$min=$k;
			}
		}
		print $aux_file $filereads[$min]->return_file_lines;
		do {
			($status,$filereads[$min])=next_read($filenames[$min],$filehandles[$min]);
		} while ($status == 1);
		if ( $status == 0 ) {
			close $filehandles[$min];
			# the file has to be deleted
			unlink $filenames[$min];
			splice @filehandles, $min, 1;
			splice @filenames, $min, 1;
			splice @filereads, $min, 1;
		}
	}

	# Append the remaining file contents to the auxiliar shorted file
	do {
		print $aux_file $filereads[0]->return_file_lines;
		do {
			($status, $filereads[0])=next_read($filenames[0],$filehandles[0]);
		} while ($status == 1);
	} while ($status != 0);
	close $filehandles[0];

	# the file has to be deleted
	unlink $filenames[0];
	close $aux_file;

	# Rename the auxiliar file to the 1st temp file
	rename($temp_file, $dir . $base . "_temp1");
	local *FILE;
	$filehandles[0]=*FILE;
	$filenames[0]=$dir . $base . "_temp1";
	$j=1;
} while ( $status_general != 0 );

# close and rename resulting file
rename($dir . $base . "_temp1", $shorted_file_name);

# Opened files get closed
close(fasta_file);
