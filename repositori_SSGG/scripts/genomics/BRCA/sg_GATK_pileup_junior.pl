#!/usr/bin/perl -w

use warnings;
use strict;
use Getopt::Long;
use Cwd;
use File::Copy;
use Statistics::Descriptive qw(:all);
use Data::Dumper;
use Tie::File;

sub getstat ($$$);

my $input_id;
my $output_id;

GetOptions (	
				"i=s"	=> \$input_id,	# input list file
				"o=s"	=> \$output_id, # output name
				
);

print "\n\n=================================================================
Obtaining Per base Pileup v 1.0\n";

if (!$input_id) {
	die "
		Options:
	-i GATK pileup file name (vcf from --output_mode EMIT_ALL_SITES)
	-o Output name path (optional)
\n".localtime()."\n=================================================================\n\n";
}

print "\nStarting process... ".localtime()."\n\n";

if (not defined $output_id) {
	my @name = split (/\./, $input_id);
	$output_id = $name[0];
}

open (VCF, "< $input_id") or die "Can´t opne $input_id\n";
open (PIL, "> $output_id.pileup") or die "Can´t create $output_id.pileup\n";
open (IGV, "> $output_id.igv") or die "Can´t create $output_id.igv\n";
print IGV "chr\tstart\tend\tprobe\t$output_id\n";

while (my $line = <VCF>) {
	chomp $line;
	if ($line =~ m/^#/) {
		next;
	}
	
	my (@data) = split(/\s+/, $line);
	my ($chr, $position, $ref, $depth);
	$chr = $data[0];
	$position = $data[1];
	$ref = $data[3];
	
	if ($data[9] eq "./.") {
		$depth = 0;
	}
	else {
		my (@fields) = split(/:/, $data[9]);
		my (@depths) = split(/,/, $fields[1]);
		if ($depths[1]) {
			$depth = $depths[0] + $depths[1];
		}		
		else {
			$depth = $depths[0];
		}
	}
	
	print PIL "$chr\t$position\t$ref\t$depth\n";
	
	my ($start, $end);
	
	$start = $position;
	$end = $position;
	
	print IGV "$chr\t$start\t$end\t",$chr."_".$position,"\t$depth\n";
}

close IGV;
close PIL;
close VCF;

print "Process finished... ".localtime()."\n";

exit;



