#!/usr/bin/perl
# This program is based on one written by
# Kay Hofmann, PhD                    
# Bioinformatics Group  
# MEMOREC Stoffel GmbH
# Stoeckheimer Weg 1
# D50829 Koeln/Germany
# which used command line input of parameters and, by default, output the
# results to the screen (but this could be redirected to a single output
# file using the unix > command.
# It was modified by DAJ to:
# (1) prompt user for the splitting parameters
# (2) output to a user-specified directory
# (3) output as a seperate, FASTA format file for each subsequence
# introductory blurb, first clearing a bit of screen
print "\n";
print "\n";
print "\n";
print "\n";
print "\n";
print "\n";
print "********************************************************************************\n";
print "*****************program to split a sequence into subsequences******************\n";
print "********************************************************************************\n";
print "\n";
print "This program is designed to take a user-specified sequence and convert it into\n";
print "subsequences, of user-specified size and overlap.\n";
print "\n";
print "\n";
print "Acceptable input formats are: FASTA, EMBL, GCG.\n";
print "The output format is fasta.\n";
print "\n";
print "\n";
###################################################################
# determine users requirements
###################################################################
# asking for the input file
print "For this program to run, the sequence file to be processed must be in the\n";
print "current working directory. Ie. in:\n";
system pwd;
print "\n";
print "Please enter the name of the sequence file to be processed:\n";
$infile = <STDIN>;
chomp $infile;
# checking to see if the file exists
until (-f $infile)  { # the -f bit means check if there is a file called ...
	print "Sorry, that file does not appear to exist, please try again.\n";
	$infile = <STDIN>;
	chomp $infile;
}
# creating a directory for the output files
print "\n";
print "What would you like to call the directory that the subsequence output files\n";
print "will be stored in? It will be created as a NEW subdirectory of the current\n";
print "working directory.\n";
$output_directory = <STDIN>;
chomp $output_directory;
# checking if this suggested directory already exists, rejecting until an unused name is supplied 
while (-d $output_directory) { # the -d bit means check if there is a directory called ...
	print "Sorry, that directory already exists, please enter another name.\n";
	$output_directory = <STDIN>;
	chomp $output_directory;
}
# once an acceptable output directory name is entered,, create a directory called it
mkdir ($output_directory, 0777) || die "sorry system cannot create output directory $output_directory";
opendir (OUTPUTDIR, $output_directory) || die "sorry, system cannot open output directory $output_directory";
# now to request the sequence splitting parameters
print "\n";
print "What size do you want the subsequences to be?\n";
$split = <STDIN>;
chomp $split;
until ($split =~ /\d/) { # checking that it is a valid answer
	print "sorry that is not an acceptable answer, please try again\n";
	$split = <STDIN>;
	chomp $split;
} 
print "\n";
print "What size do you want the overlap beteween subsequences to be?\n";
$overlap = <STDIN>;
chomp $overlap;
until ($overlap =~ /\d/) { # checking that it is a valid answer
	print "sorry that is not an acceptable answer, please try again\n";
	$overlap = <STDIN>;
	chomp $overlap;
} 
# printing confirmation of files and parameters
print "\n";
print "Confirmation - the program is working with file: $infile\n";
print "and will output to directory: $output_directory\n";
print "\n";
print "subsequences will be $split long\n";
print "and overlap by $overlap.\n";
print "\n";
# now need to open file $infile for reading from, and determine the sequence file format
open (INFILE, $infile) || die "sorry, system can't open $infile for reading from";
$line=<INFILE>; # grab the first line
if ($line =~ /^\>/) {
	$format="FASTA";
	print "input sequence format is fasta\n";
}
elsif ($line =~ /^ID\s/){
	$format="EMBL";
	print "input sequence format is EMBL\n";}
else {
	while (!($line =~ /\.\./) && ($line=<INFILE>)) {
	}
	if ($line=~/\.\./) {
		$format="GCG";
		print "input sequence format is GCG\n";
	}
	else {
		die "Sorry, sequence format is not recognized\n";
	}
}
close (INFILE);
###################################################################
# the main program
###################################################################
&getseq; # go to subroutine getseq
$seqlen = length($Seq); # determine length of sequence
print "sequence length = $seqlen\n";
print "\n";
$filecount = 1; # set marker for number of files created
$pos1 = 1; # set position of start of subsequence
$seqleft = $seqlen; # set counter for sequence remaining to be processed
# when dividing the sequence into subsequences, the header of the
#subsequence files will contain a counter for the part of the sequence
#contained in the file (from $filecount), plus the start and end points
#of the subsequence
#The second line will contain the size of the subsequence. Because the
#size of the sequences will rarely be an exact multiple of the
#subsequence size, it is necessary to divide the parsing of the sequence
#into 2 bits:
#firstly process all of the subsequences that can be made of size
# = $split
#then  process the final bit left over (< $split in size).
#Similarly, within each subsequence, when it is parsed, it is
#reformatted into 60 base long lines. Because the length of the
#subsequence will rarely be an exact multiple of the line length, it is
#necessary to divide the parsing of the subsequence into 2 bits:
#firstly process all of the output lines that can be made of length 60
#then  process the final bit left over (<60 in size).
#processing all of the subsequences that can be made of size = $split
while ($seqleft > $split) {
	$pos2 = $pos1+$split-1;
	$outfile = $output_directory."/".$infile."_".$pos1."to".$pos2;
	open (OUTFILE, ">>$outfile") || die "sorry, system can't open outfile $outfile for writing to";
	print "outfile is called $outfile\n";
	print OUTFILE ">".$SeqID."_".$filecount." ".$pos1." to ".$pos2." ". $seqDE."\n";
	print OUTFILE $split."\n";
	$part = substr($Seq,($pos1 - 1),$split);
	$partlength = $split;
	$partpos = 0;
	# outputting all lines of length 60
	while ($partlength >60) { # the default line length
		$subpart = substr($part, $partpos, 60);
		print OUTFILE $subpart."\n";
		$partpos = $partpos + 60;
		$partlength = $partlength - 60;
	}
	#outputting the final line of length <60
	$subpart = substr($part, $partpos);
	print OUTFILE $subpart."\n";
	$filecount = $filecount + 1;
	$pos1 = $pos2 - $overlap + 1;
	$seqleft = $seqleft - $split + $overlap;
	print "sequence remaining = $seqleft\n";
	print "\n";
	close OUTFILE;
}

#processing the final subsequence of size < $split
$outfile = $output_directory."/".$infile."_".$pos1."to".$seqlen;
open (OUTFILE, ">>$outfile") || die "sorry, system can't open outfile for writing to";
print OUTFILE ">".$SeqID."_".$filecount." ".$pos1." to ".$seqlen." ". $seqDE."\n";
$part = substr($Seq, ($pos1 -1));
$finalpart = length($part);
print OUTFILE $finalpart."\n";
$partlength = $finalpart;
$partpos = 0;
# outputting all lines of length 60
while ($partlength >60) { # the default line length
	$subpart = substr($part, $partpos, 60);
	print OUTFILE $subpart."\n";
	$partpos = $partpos + 60;
	$partlength = $partlength - 60;
}
#outputting the final line of length <60
$subpart = substr($part, $partpos);
print OUTFILE $subpart."\n";
print "outfile is called $outfile\n";
print "sequence remaining = 0\n";
close (OUTFILE);
close (INFILE);
closedir (OUTPUTDIR);
print "\n";
print "\n";
print "*****FINISHED PROCESSING*****\n";
print "$filecount files created\n";
# end of program
##################################################################
# subroutine to get sequence in any format and put it into Seq, SeqID and SeqDE
##################################################################
sub getseq {
	open (INFILE, $infile) || die "sorry, system can't open $infile for reading from";
	if ($format eq "FASTA") {
		&getfasta;
		$Seq = $fastabuffer;
		$SeqID = $FastaID;
		$SeqDE = $FastaDE;
	}
	elsif ($format eq "EMBL") {
		&getembl;
		$Seq = $emblbuffer;
		$SeqID = $EmblID;
		$SeqDE = $EmblDE;
	}
	else {
		&getgcg;
		$Seq = $gcgbuffer;
		$SeqID = $gcgid;
		$SeqDE = $gcgde;
	}
}
###################################################################
# subroutine to read Pearson/FASTA format sequence (not to be called externally) 
###################################################################
sub getfasta {
	$fastabuffer="";
	$FastaID="";
	$FastaDE="";
	$line="";
	until (($fastaline =~ /^\>/) || eof(INFILE)) {
		$fastaline=<INFILE>;
	}
	if ($fastaline =~/^\>(\S+)\s(.*)$/) {
		$FastaID = $1;
		$FastaDE = $2;
	}
	while ($FastaID =~/^(.*\|)(.+)$/) {
		$FastaID = $2;
		}
	until (($line =~ /^\>/) || eof(INFILE)) {
		$line = <INFILE>;
		if (!($line =~ /^\>/)) {
			$fastabuffer .= $line
		}
	}
	if ($line =~ /^\>/) {
		$fastaline = $line;
	}
	else {
		$fastaline = "";
	}
	$fastabuffer =~ tr/a-z/A-Z/;
	$fastabuffer =~ s/[^A-Z]//g;
}
###################################################################
# subroutine to read EMBL/Swissprot format sequence (not to be called externally) 
###################################################################
sub getembl {
	$emblbuffer = "";
	$EmblID = "";
	$EmblDE = "";
	$line = "";
	until (($line =~ /^ID\s/) || eof(INFILE)) {
		$line=<INFILE>;
	}
	if ($line =~/^ID\s+(\w+).*$/) {
		$EmblID=$1;
	}
	until (($line =~ /^SQ\s/) || eof(INFILE)) {
		$line=<INFILE>;
		if ($line =~ /^DE\s+(.*)/) {
			if($EmblDE) {
				$EmblDE .=" ";
			}
			$EmblDE .= $1;
		}
	}
	if ($line =~ /^SQ\s/) {
		until (($line =~ /^\/\//) || eof(INFILE)) {
			$line = <INFILE>;
			if (!($line =~ /^\/\//)) {
				$emblbuffer .= $line;
			}
		}
	}
	$emblbuffer =~ tr/a-z/A-Z/;
	$emblbuffer =~ s/[^A-Z]//g;
}
###################################################################
# subroutine to read GCG format sequence (not to be called externally) 
###################################################################
sub getgcg {
	$gcgbuffer = "";
	$gcgid = "";
	$line = "";
	until (($line =~ /\.\./) || eof(INFILE)) {
		$line = <INFILE>;
	}
	if ($line =~/^(\w*).*\.\./) {
		$gcgid=$1;
	}
	until (eof(INFILE)) {
		$line = <INFILE>;
		$gcgbuffer .= $line;
	}
	$gcgbuffer =~ tr/a-z/A-Z/;
	$gcgbuffer =~ s/[^A-Z]//g;
}

