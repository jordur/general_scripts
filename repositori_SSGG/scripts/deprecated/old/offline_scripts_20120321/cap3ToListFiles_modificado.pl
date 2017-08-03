#!/usr/bin/perl

use strict;
use warnings;

#parse a cap3 output file and create two files with info about which sequences are in 
#which contig
#first file will have format of one line per sequence:  contigName  seqName
#second file will have format of one line per contig:   contigName  seqName1,seqName2,seqName3
#the latter file is primarily used as input for a script that creates genelists for GeneSpring

#here we want to process all the .ace files in a directory.
#I then  need to concatenate all these to get my final lists.
#This produces ludicrous numbers of output files because I'm sending the job off to a cluster and can't 
#have all the processes trying to access a single file...there must be a smarter way to do this. 


my $outputDir = ".";
my $multiIDLine = 0; #use this if you are parsing out an id from a longer line.
#Check out the line below that deals with which elements to keep. Line starting if(multiIDLine)
#######################################

if ($#ARGV != 0) { die "usage: < cap3ToListFiles.pl inputCap3File >\n\n"; }

my $file = $ARGV[0];

my %contHash;
my $contName;

#split filename on . and take first part only
$file =~ s|\.\/||;
my @nameBits = split(/\./, $file);

#my $outfile1 = $outputDir . "/" . $nameBits[0] . "_outfile1.txt";
my $outfile2 = $outputDir . "/" . $nameBits[0] . "_outfile2.xls";
#print "outfile1 is $outfile1 and outfile2 is $outfile2\n";
open (CAP3FILE, $file) or die "I could not open $file: $! \n\n";


while ( <CAP3FILE> ) {
    chomp($_);
    s/\x0D//;  #get rid of ^M if its there.
#search for next occurrence of "CO Contig"
	if ( /CO (Contig\d*) / )  
	{ 
		$contName = $1;
		$contHash{$contName} = [];  
	}
	#so contName should change everytime we have a new contig to deal with
	#because its a global variable, it should remain the hashname until it changes again when
	#a new CO line is found	
	#if ( /^AF( .+) \w \-?\d+/ ) 
	if ( /^AF( .+) \w \-?\d+$/ )
	{
		my $names = $1;
		# #remove the .ab1 name from the filenames if it appears
		$names =~ s/.ab1//g;
		
		if ($multiIDLine) 
		{
			my @seqIds = split(/\|/, $names);
			push (@{$contHash{$contName}},$seqIds[4]);
		}

		else { push(@{$contHash{$contName}},$names) } ;
		
	}
}

close CAP3FILE;
#`print to file
#open (OUTFILE1, ">$outfile1") or warn "I could not open outfile1.\n";
open (OUTFILE2, ">$outfile2") or warn "I could not open outfile2.\n";

my $firstTouch = 1;
my $contador=0;
print OUTFILE2 "contig","\t","num_secuencias","\t","seqs_in_contig\n";

for my $contig ( sort (keys %contHash) ) 
{
	unless ($firstTouch) {  print OUTFILE2 "\n"; }
	$firstTouch = 0;
        my @seqNames = @{$contHash{$contig}};
       
	print OUTFILE2 "$contig","\t",$#seqNames+1,"\t";

	my $seqNum =  0;	
	foreach my $sname (@seqNames) 
	{	
		  #$contador ++;$sname+1\n";
       # print OUTFILE1 "$sname\t$nameBits[0]\t$contig\n";
	print OUTFILE2 "$sname";
	if ($seqNum < ($#seqNames)) { print OUTFILE2 ","; }
	$seqNum++;
	}
}
	
#close OUTFILE1;
close OUTFILE2;

print "Done $file\n";

