#!/usr/bin/perl

#=================================================#
# Getting coverage from global file to probe file #
#=================================================#


use warnings;
use strict;
use Getopt::Long;
use Cwd;
use File::Copy;
#use Parallel::ForkManager;
use Statistics::Descriptive qw(:all);
use Data::Dumper;
use Tie::File;

my $output_id;
my $coverage_id;
my $bed_id;

GetOptions (	"o=s"	=> \$output_id,	# output file root
				"i=s" 	=> \$coverage_id, # coverage file root
				"b=s"   => \$bed_id, #bed file
);

print "\n\n=================================================================
Getting global coverage stat script v 1.0\n";

if (not defined $coverage_id or not defined $bed_id) {
	die "
			Options:
		-i Input Pileup file name (chr position base depth)
		-b BED file (chr start end gene_name exon_number) 
		-o Output file name root (optional)
\n".localtime()."\n=================================================================\n\n";
}

print "\n Starting process... ".localtime()."\n\n";

if (not defined $output_id) {
	my @name = split (/\./, $coverage_id);
	$output_id = $name[0]."_median";
}

open (BED, "< $bed_id") or die "Can´t open $bed_id\n";
open (PROBE, "> $output_id.probes") or die "Can't open $output_id.probes\n";
open (PIL, "> $output_id.allbases_pileup") or die "Can't open $output_id.allbases_pileup\n";
my %Global_depth;

print "Analising positions...".localtime()."\n";

while (my $line = <BED>) {
	
	my ($chr, $start, $end, $gene, $exon) = split (/\s+/, $line);
	if ($chr =~ m/Chr/) {
		next;
	}
	$chr =~ s/chr//g;
	
	my @probe_depth;
	#print "Analising $line...".localtime()."\n";
	my $length = $end - $start;
	
	for (my $i = 0; $i <= $length; $i++) {
		my $depth;
		my $pos = $start + $i;
		my $col = '$4';
		#print "grep -m 1 chr$chr\[\[:space:\]\]$pos\[\[:space:\]\] $coverage_id | awk \'\{print $col\}\'\n";
		
		my $data = `grep -m 1 chr$chr\[\[:space:\]\]$pos\[\[:space:\]\] $coverage_id | awk \'\{print $col\}\'`;
		chomp $data;
		
		$depth = ($data || '0');
		
		push (@probe_depth, $depth);
		
		if (defined $Global_depth{$gene}{depth}) {
			$Global_depth{$gene}{depth} .= "$depth ";	
		}
		else {
			$Global_depth{$gene}{depth} = "$depth ";
		}
		
		if (defined $Global_depth{depth}){
			$Global_depth{depth}.= "$depth ";
		}	
		else {
			$Global_depth{depth} = "$depth ";
		}
		print PIL "$chr\t$pos\t$depth\n";
	}
	
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@probe_depth);
	my $median = $stat->median();
	my $mean = $stat->mean();
	
	print PROBE "$chr\t$start\t$end\t$gene\t$exon\t$median\t$mean\n";

}

close BED;
close PROBE;

open (GENE, "> $output_id.genes") or die "Can't open $output_id.genes\n";
open (BED, "< $bed_id") or die "Can´t open $bed_id\n";

my $last_gene = 'NA';

print "Printing values...".localtime()."\n";

while (my $line = <BED>) {
	my ($chr, $start, $end, $gene, $exon) = split (/\s+/, $line);
	if ($chr =~ m/Chr/) {
		next;
	}
	if ($gene ne $last_gene) {
		my $stat = Statistics::Descriptive::Full->new();
		my @Depth = split (/\s+/, $Global_depth{$gene}{depth});
		$stat->add_data(@Depth);
		my $median = $stat->median();
		my $mean = $stat->mean();
		print GENE "$gene\t$median\t$mean\n";
		$last_gene = $gene;	
	}
}

my $stat = Statistics::Descriptive::Full->new();
my @Depth = split (/\s+/, $Global_depth{depth});
$stat->add_data(@Depth);
my $median = $stat->median();
my $mean = $stat->mean();
print GENE "Total_sample\t$median\t$mean\n";

print "Process Finished...".localtime()."\n";

close GENE;
close BED;

exit;