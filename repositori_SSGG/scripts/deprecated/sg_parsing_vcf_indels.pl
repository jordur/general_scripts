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

my $vcf_id;
my $output_id;

GetOptions (	
				"v=s"	=> \$vcf_id,  # probes file		
);


if (!$vcf_id) {
	die "
\n\n=================================================================
	Parsing script v 1.0
			Options:
		-v VCF file
\n".localtime()."\n=================================================================\n\n";
}

open (VCF, "< $vcf_id") or die "Can´t open $vcf_id\n";

while (my $line = <VCF>) {
	chomp $line;
	if ($line =~ m/^#/) {
		print "$line\n";
		next;
	}	
	else {
		my ($ref, $var, $pos);
		my @data = split (/\s+/, $line);

		if ($data[4] =~ m/\,/) {
			my @alleles = split (/\,/, $data[4]);
			my $allele;
			foreach my $al (@alleles) {
				if (defined $allele) {
					if (length($al) < length($allele)){
						$allele = $al;
					}
				}
				else {
					$allele = $al;
				}
			}
			$data[4] = $allele;
		}		

		if (length($data[3]) > length($data[4])) { #deletion
			$var = '-';
			$ref = $data[3];
			$ref =~ s/$data[4]//;
			$data[3] = $ref;
			$data[4] = $var;
			$pos = $data[1] + 1;
			$data[1] = $pos;
			for (my $i = 0; $i < scalar(@data); $i++) {
				print "$data[$i]\t";		
			}
			print "\n";
		}
		elsif (length($data[3]) < length($data[4])) { #insertion
			$ref = '-';
			$var = $data[4];
			$var =~ s/$data[3]//;
			$data[3] = $ref;
			$data[4] = $var;
			$pos = $data[1] + 1;
			$data[1] = $pos;
			for (my $i = 0; $i < scalar(@data); $i++) {
				print "$data[$i]\t";		
			}
			print "\n";
		}
		else {
			print "$line\n";
		}
	}
}
close VCF;

exit;


