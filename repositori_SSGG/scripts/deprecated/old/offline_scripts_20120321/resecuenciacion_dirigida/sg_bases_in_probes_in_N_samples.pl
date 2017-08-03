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

my $samples_id;
my $probe_id;
my $output_id;

GetOptions (	
				"i=s"	=> \$samples_id,	# input list file
				"b=s"	=> \$probe_id,  # probes file
				"o=s"	=> \$output_id, # output name
				
);

print "\n\n=================================================================
Normalisation step 1 script v 1.0\n";

if (!$samples_id || !$probe_id) {
	die "
			Options:
		-i Input list file name (Chr Pos Sample1 Sample2 ...)
		-b Probes file (Chr Start End Gene Name)
		-o Output name path (optional)
\n".localtime()."\n=================================================================\n\n";
}

print "\n Starting process... ".localtime()."\n\n";

if (not defined $output_id) {
	my @name = split (/\./, $samples_id);
	$output_id = $name[0];
}

my %Global_depth;
my $header = `grep -i -m 1 ^#chr $samples_id`;
chomp $header;
my @header = split (/\s+/, $header);
my %Samples;
my $samples_number = 0;
my @samples;

for (my $i = 2; $i < (scalar(@header)); $i++) {
	$Samples{$samples_number}{name} = $header[$i];
	$samples_number++;
	push (@samples, $header[$i]);
}
print "Number of samples = $samples_number \n Samples = @samples\n";

undef @header;

open (BED, "< $probe_id") or die "Can´t open $probe_id\n";

my $samples;
my $genes;
my $runs;
my %runs;

open (OUT, "> $output_id.byprobes") or die "Can´t create $output_id.byprobes\n";
print OUT "Gene\tProbe";
foreach my $sample (@samples) {
	print OUT "\t$sample";
}
print OUT "\n";

while (my $probe = <BED>) {
	chomp $probe;
	my ($chr, $start, $end, $gene, $name) = split (/\s+/, $probe);
	if ($start =~ m/start/i) {
		next;
	}
	
	my %samples_depth;
	
	print OUT "$gene\t$name";
	my $probe_depth;
    my $command = '($1 == "'.$chr.'" || $1 == "chr'.$chr.'") && $2 >= '.$start.' && $2 <= '.$end;
	my $length = $end - $start;
	
	my @data = `awk \'$command\' $samples_id`; 
	chomp @data;
	
	foreach my $line (@data) {
		chomp $line;
		my @fields = split (/\s+/, $line);
		
		my ($chrom) = shift (@fields);
		my ($position) = shift (@fields);
		
		for (my $i = 0; $i < $samples_number; $i++) {
			if (defined $samples_depth{$i}{probe_depth}) {
				$samples_depth{$i}{probe_depth} .= " $fields[$i]";
			}
			else {
				$samples_depth{$i}{probe_depth} = $fields[$i];
			}
		}
	}
	
	for (my $i = 0; $i < $samples_number; $i++) {
		my $median;
		my $mean;
		getstat($samples_depth{$i}{probe_depth}, \$mean, \$median);
		print OUT "\t$median";
	}
	print OUT "\n";
}
close OUT;
close BED;

print "Process finished... ".localtime()."\n";
exit;

sub getstat ($$$) {
	my ($data, $mean, $median) = @_;
	my $stat = Statistics::Descriptive::Full->new();
	my @Depth = split (/\s+/, $data);
	$stat->add_data(@Depth);
	$$median = $stat->median();
	$$mean = $stat->mean();
	#print "$$mean\t$$median\n";
}
