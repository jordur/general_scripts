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
my $junior_flag;

GetOptions (	
				"p=s"	=> \$pileup_id,	# input list file
				"v=s"	=> \$vcf_id,  # probes file		
				"j"		=> \$junior_flag, #flag for junior analysis
);


if (!$pileup_id || !$vcf_id) {
	die "
\n\n=================================================================
	Normalisation step 1 script v 1.0
			Options:
		-p GATK indels pileup file name
		-v VCF file name
		-j for Junior analysis
n".localtime()."\n=================================================================\n\n";
}

my %Depth_changes;

open (VCF, "< $vcf_id") or die "Can´t open $vcf_id\n";

my ($print_line, $last_chr, $last_pos, $last_type, $last_length);

while (my $line = <VCF>) {
	chomp $line;
	my @data = split (/\s+/, $line);
	if ($line =~ m/^#/ || length($data[3]) == length($data[4])) {
		print "$line\n";
		next;
	}	
	else {
		my $type;
		my $chr = $data[0];
		my $pos = $data[1];
		if (length($data[3]) > length($data[4])) {
			$type = 'del';
		}
		elsif (length($data[3]) < length($data[4])) {
			$type = 'ins';
		}
		
		my $length = abs(length($data[3]) - length($data[4]));
		
		if ($print_line) {
			if ($chr eq $last_chr && ($pos - $last_pos) <= abs(length($data[3]) - length($data[4])) && $type eq $last_type){
				my @last_data = split(/\s+/, $print_line);
				my @data = split(/\s+/, $line);
				my ($depth, $var, $freq, $last_depth, $last_var, $last_freq);
				my @col10 = split(/:/, $data[9]);
				my @col8 = split (/;/, $data[7]);
				$depth = $col10[1];
				if ($depth == "0") { $depth = 1;}
				if ($col8[4] =~ m/ds/i) {
					$col8[5] =~ s/Dels=//;
					$freq = $col8[5];
				}
				else {
					$col8[4] =~ s/Dels=//;
					$freq = $col8[4];
				}
				if ($freq eq '') { $freq = 0;}
				$var = sprintf("%.0f",($freq * $depth)/ (1 - $freq));
				@col10 = split(/:/, $last_data[9]);
				@col8 = split (/;/, $last_data[7]);
				$last_depth = $col10[1];
				if ($col8[4] =~ m/ds/i) {
					$col8[5] =~ s/Dels=//;
					$last_freq = $col8[5];
				}
				else {
					$col8[4] =~ s/Dels=//;
					$last_freq = $col8[4];
				}
				if ($last_freq eq '') { $last_freq = 0;}
				$last_var = sprintf("%.0f",($last_freq * $last_depth)/ (1 - $last_freq));
				$last_depth -= $var;
				$depth += $var;
				my $new_freq;
				if(($last_depth + $last_var + $var) == 0) {
					$new_freq = 0;
				}
				else {
					$new_freq = sprintf("%.2f",(($var+$last_var)/($last_depth + $last_var + $var)));
				}
				
				for (my $i = $pos + 1; $i< ($pos + $length) + 1; $i++) {
					$Depth_changes{$chr}{$i}{depth} = $depth;
				}
				for (my $i = $last_pos + 1; $i< ($last_pos + $last_length) + 1; $i++) {
					$Depth_changes{$chr}{$i}{depth} = $last_depth;
				}

				if ($last_freq > 0.90) {
					$col8[0] = "AC=2";
					$col8[1] = "AF=1.00";
					$col8[2] = "AN=1";
				}
				else {
					$col8[0] = "AC=1";
					$col8[1] = "AF=0.50";
					$col8[2] = "AN=2";
				}
				
				if ($col8[4] =~ m/ds/i) {
					$col8[5] = "Dels=$new_freq";;
				}
				else {
					$col8[4] = "Dels=$new_freq";
				}
				$col10[1] = $last_depth;
				
				$last_data[9] = join(":", @col10);
				$last_data[7] = join(";", @col8);
				$print_line = join("\t", @last_data);
				
			}
			else {
				print "$print_line\n";
				$print_line = $line;
				$last_chr = $chr;
				$last_pos = $pos;
				$last_type = $type;
				$last_length = $length;
			}
		}
		else {
			$print_line = $line;
			$last_chr = $chr;
			$last_pos = $pos;
			$last_type = $type;
			$last_length = $length;
		}
	}
}
close VCF;

open (PIL, "< $pileup_id") or die "Can´t open $pileup_id\n";
my @name = split (/\./, ($pileup_id));
pop (@name);
my $outname = join (".", @name);

open (OUT, "> $outname.vcf") or die "Can´t create $outname.vcf\n";

while (my $line = <PIL>) {
	chomp $line;
	if ($line =~ m/^#/) {
		print OUT "$line\n";
		next;
	}	
	else {
		my @data = split (/\s+/, $line);
		if ($junior_flag) {
			if ($data[0] =~ m/BRCA2/) {
				$data[0] = "chr13";
				$data[1] += 32889616;
			}
			elsif ($data[0] =~ m/BRCA1/) {
				$data[0] = "chr17";
#				my $tmp_pos = 41277501 - $data[1];
				$data[1] += 41196311;
			}
		}

		my $chr = $data[0];
		my $pos = $data[1];
		
		my $print_line;
		
		if ($Depth_changes{$chr}{$pos}) {
			my @col10 = split (/:/, $data[9]);
			$col10[1] = $Depth_changes{$chr}{$pos}{depth};
			$data[9] = join (":", @col10);
		}

		$print_line = join ("\t", @data);
		print OUT "$print_line\n";
	}
}
close OUT;
close PIL;

exit;


