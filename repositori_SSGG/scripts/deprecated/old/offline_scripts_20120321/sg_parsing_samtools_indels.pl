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

my $pileup_id;
my $vcf_id;
my $output_id;

GetOptions (	
				"p=s"	=> \$pileup_id,	# input list file
				"v=s"	=> \$vcf_id,  # probes file		
);


if (!$pileup_id || !$vcf_id) {
	die "
\n\n=================================================================
	Normalisation step 1 script v 1.0
			Options:
		-p samtools indels pileup file
		-v VCF file
n".localtime()."\n=================================================================\n\n";
}

open (VCF, "< $vcf_id") or die "Can´t open $vcf_id\n";

while (my $line = <VCF>) {
	chomp $line;
	if ($line =~ m/^#/) {
		print "$line\n";
		next;
	}	
	else {
		my @data = split (/\s+/, $line);
		if ($data[7] =~ m/indel/i) {
				my $command = '$1 == "'.$data[0].'" && $2 == '.$data[1];
				my @info =  split (/\s+/,`awk \'$command\' $pileup_id `);
				my $freq;
				if ($info[8] eq '*') {
					my $rate = $info[11] / $info[7];
					$freq = sprintf("%.2f",$rate);
				}
				elsif ($info[9] eq '*') {
					my $rate = $info[10] / $info[7];
					$freq = sprintf("%.2f",$rate);
				}
				else {
					$freq = 0;
				}
				print "$data[0]\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$data[5]\t$data[6]\t";
				my @col8 = split(/;/, $data[7]);
				if ($freq > 0.85) {
					print "AC=2;AF=1.00;AN=1;";
				}	
				else {
					print "AC=1;AF=0.50;AN=2;";
				}	
				my $sb = "NA";
				if ($col8[7]){
					my @SB = split(/,/, $col8[7]);
					$SB[0] =~ s/PV4=//g;
					$sb = "pval($SB[0])";
				}
				print "$col8[1];Dels=$freq;HRun=NA;HaplotypeScore=NA;$col8[5];MQ0=0;QD=NA;SB=$sb;sumGLbyD=NA;SF=NA\tGT:DP:GQ:PL\t";
				my @flag = split (/:/, $data[9]);
				print "$flag[0]:$info[7]:$flag[2]:$flag[1]\n";
		}
		else {
			print "$line\n";
		}
	}
}
close VCF;

exit;


