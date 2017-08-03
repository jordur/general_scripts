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

my $list_id;
my $probe_id;
my $output_id;

GetOptions (	
				"i=s"	=> \$list_id,	# input list file
				"b=s"	=> \$probe_id,  # probes file
				"o=s"	=> \$output_id, # output name
				
);

print "\n\n=================================================================
Normalisation step 1 script v 1.0\n";

if (!$list_id || !$probe_id) {
	die "
			Options:
		-i Input list file name (Sample_file_path Sample_name Run_id)
		-b Probes file (Chr Start End Name)
		-o Output name path (optional)
\n".localtime()."\n=================================================================\n\n";
}

print "\n Starting process... ".localtime()."\n\n";

if (not defined $output_id) {
	my @name = split (/\./, $list_id);
	$output_id = $name[0];
}

my %Global_depth;

open (BED, "< $probe_id") or die "Can´t open $probe_id\n";
open (LIST, "< $list_id") or die "Can´t open $list_id\n";
my @Data = <LIST>;
close LIST;

my $samples;
my $genes;
my $runs;
my %runs;

while (my $probe = <BED>) {
	chomp $probe;
	my ($chr, $start, $end, $gene, $name) = split (/\s+/, $probe);
	if ($chr =~ m/chr/i) {
		next;
	}
	if (!$genes) {
		$genes = $gene;
	}
	elsif ($genes && $genes !~ m/$gene/i) {
		$genes .= " $gene";
	}
	
	my $probe_depth;
	my $col = '$3';
    my $command = '($1 == '.$chr.' || $1 == "chr'.$chr.'") && $2 >= '.$start.'&& $2 <= '.$end;
	my $length = $end - $start;
	
	foreach my $sample_line (@Data) {
		chomp $sample_line;
		my ($file, $sample, $run) = split (/\s+/, $sample_line);
		
		if ($runs{$run}{samples} && $runs{$run}{samples} !~ m/$sample/i ) {
			$runs{$run}{samples} .= " $sample";
		}
		elsif (!$runs{$run}{samples}) {
			$runs{$run}{samples} = $sample;
		}
		#print "Analising $sample...".localtime()."\n";
		if ($samples && $samples !~ m/$sample/i) {
			$samples .= " $sample";
		}
		elsif (!$samples) {
			$samples = "$sample";
		}
		
		if ($runs && $runs !~ m/$run/i) {
			$runs .= " $run";
		}
		elsif (!$runs) {
			$runs = "$run";
		}

		my $probe_depth;
		my @data = `awk \'$command {print $col}\' $file`;
		chomp @data;
		for (my $i = 0; $i <= $length; $i++) {
			my $depth = 0;
			
			if (defined $data[$i]) {
				$depth = $data[$i];
			}
			
			my $pos = $start + $i;
			#print "$chr\t$pos\t$depth\n";
			if ($probe_depth) {
				$probe_depth .= " $depth";
			}		
			else {
				$probe_depth = "$depth";
			}
			
			# Global data
			if (defined $Global_depth{global}{depth}) {
				$Global_depth{global}{depth} .= " $depth";	
			}
			else {
				$Global_depth{global}{depth} = "$depth";
			}
			if (defined $Global_depth{global}{$gene}{depth}) {
				$Global_depth{global}{$gene}{depth} .= " $depth";	
			}
			else {
				$Global_depth{global}{$gene}{depth} = "$depth";
			}
			if (defined $Global_depth{run}{$run}{depth}) {
				$Global_depth{run}{$run}{depth} .= " $depth";	
			}
			else {
				$Global_depth{run}{$run}{depth} = "$depth";
			}
			if (defined $Global_depth{run}{$run}{$gene}{depth}) {
				$Global_depth{run}{$run}{$gene}{depth} .= " $depth";	
			}
			else {
				$Global_depth{run}{$run}{$gene}{depth} = "$depth";
			}
			
			if (defined $Global_depth{sample}{$sample}{depth}) {
				$Global_depth{sample}{$sample}{depth} .= " $depth";	
			}
			else {
				$Global_depth{sample}{$sample}{depth} = "$depth";
			}
			if (defined $Global_depth{sample}{$sample}{$gene}{depth}) {
				$Global_depth{sample}{$sample}{$gene}{depth} .= " $depth";	
			}
			else {
				$Global_depth{sample}{$sample}{$gene}{depth} = "$depth";
			}

			#Probe data
			if (defined $Global_depth{probe}{global}{$chr}{$start}{depth}) {
				$Global_depth{probe}{global}{$chr}{$start}{depth} .= " $depth";	
			}
			else {
				$Global_depth{probe}{global}{$chr}{$start}{depth} = "$depth";
			}
			if (defined $Global_depth{probe}{$run}{$chr}{$start}{depth}) {
				$Global_depth{probe}{$run}{$chr}{$start}{depth} .= " $depth";	
			}
			else {
				$Global_depth{probe}{$run}{$chr}{$start}{depth} = "$depth";
			}
			
		}
		my $median;
		my $mean;
		getstat($probe_depth, \$mean, \$median);
		$Global_depth{probe}{$sample}{$chr}{$start}{median} = $median;	
		$Global_depth{probe}{$sample}{$chr}{$start}{mean} = $mean;
	}
}
close BED;

open (GRUN, "> $output_id\_allruns.genes") or die "Can´t create $output_id\__allruns.genes\n";
print GRUN "Run\tGene\tMedian\tMean\n";

foreach my $run ((split(/\s+/, $runs))) {
	foreach my $gene ((split(/\s+/, $genes))) {
		my $median;
		my $mean;
		getstat($Global_depth{run}{$run}{$gene}{depth}, \$mean, \$median);
		print GRUN "$run\t$gene\t$median\t$mean\n";
		$Global_depth{run}{$run}{$gene}{median} = $median;
		$Global_depth{run}{$run}{$gene}{mean} = $mean;
		undef $Global_depth{run}{$run}{$gene}{depth};
	}
	my $median;
	my $mean;
	getstat($Global_depth{run}{$run}{depth}, \$mean, \$median);
	print GRUN "$run\tBoth\t$median\t$mean\n";
	$Global_depth{run}{$run}{median} = $median;
	$Global_depth{run}{$run}{mean} = $mean;
	#undef $Global_depth{run}{$run}{depth};
}

foreach my $gene ((split(/\s+/, $genes))) {
	my $median;
	my $mean;
	getstat($Global_depth{global}{$gene}{depth}, \$mean, \$median);
	print GRUN "All\t$gene\t$median\t$mean\n";
	$Global_depth{global}{$gene}{median} = $median;
	$Global_depth{global}{$gene}{mean} = $mean;
	undef $Global_depth{global}{$gene}{depth};
}

my $G_median;
my $G_mean;
getstat($Global_depth{global}{depth}, \$G_mean, \$G_median);
print GRUN "All\tBoth\t$G_median\t$G_mean\n";
$Global_depth{global}{mean} = $G_mean;
$Global_depth{global}{median} = $G_median;
#undef $Global_depth{global}{depth};
close GRUN;

foreach my $run ((split(/\s+/, $runs))) {
	open (RUN, "> $run\_allsamples.genes") or die "Can´t create $run\_allsamples.genes\n";
	print RUN "Sample\tGene\tMedian\tMean\n";
	foreach my $sample ((split(/\s+/, $runs{$run}{samples}))) {
		foreach my $gene ((split(/\s+/, $genes))) {
			my $median;
			my $mean;
			getstat($Global_depth{sample}{$sample}{$gene}{depth}, \$mean, \$median);
			print RUN "$sample\t$gene\t$median\t$mean\n";
			undef $Global_depth{sample}{$sample}{$gene}{depth};		
		}
		my $median;
		my $mean;
		getstat($Global_depth{sample}{$sample}{depth}, \$mean, \$median);
		print RUN "$sample\tBoth\t$median\t$mean\n";
		$Global_depth{sample}{$sample}{median} = $median;
		$Global_depth{sample}{$sample}{mean} = $mean;
	#	undef $Global_depth{sample}{$sample}{depth};		
	}
	close RUN;
}

open (BED, "< $probe_id") or die "Can´t open $probe_id\n";
open (GNORM, "> $output_id\_global.normalised") or die "Can´t create $output_id\_global.normalised\n";

print GNORM "probe";
foreach my $sample ((split(/\s+/, $samples))) {
	print GNORM "\t$sample";
}
print GNORM "\n";

while (my $probe = <BED>) {
	chomp $probe;
	my ($chr, $start, $end, $gene, $name) = split (/\s+/, $probe);
	if ($chr =~ m/chr/i) {
		next;
	}
	print GNORM "$name";
	my $glob_mean;
	my $glob_median;
	getstat($Global_depth{probe}{global}{$chr}{$start}{depth}, \$glob_mean, \$glob_median);
	my $glob_rate =  $glob_median / $Global_depth{global}{median};
	foreach my $sample ((split(/\s+/, $samples))) {
		if ($Global_depth{sample}{$sample}{median} != 0) {
			my $sample_rate = $Global_depth{probe}{$sample}{$chr}{$start}{median} / $Global_depth{sample}{$sample}{median};
			my $norm_rate = $sample_rate / $glob_rate;
			print GNORM "\t$norm_rate";
		}
		else {
			print GNORM "\tN/A";
		}
	}
	print GNORM "\n";
}
close GNORM;
close BED;

foreach my $run ((split(/\s+/, $runs))) {
	open (RNORM, "> $run\_run.normalised") or die "Can´t create $run\_run.normalised\n";
	print RNORM "probe";
	foreach my $sample ((split(/\s+/, $runs{$run}{samples}))) {
		print RNORM "\t$sample";
	}
	print RNORM "\n";
	
	open (BED, "< $probe_id") or die "Can´t open $probe_id\n";
	while (my $probe = <BED>) {
		chomp $probe;
		my ($chr, $start, $end, $gene, $name) = split (/\s+/, $probe);
		if ($chr =~ m/chr/i) {
			next;
		}
		print RNORM "$name";
		my $run_mean;
		my $run_median;
		getstat($Global_depth{probe}{$run}{$chr}{$start}{depth}, \$run_mean, \$run_median);
		my $run_rate =  $run_median / $Global_depth{run}{$run}{median};
		foreach my $sample ((split(/\s+/, $runs{$run}{samples}))) {
			if ($Global_depth{sample}{$sample}{median} != 0 && $run_rate != 0) {
				my $sample_rate = $Global_depth{probe}{$sample}{$chr}{$start}{median} / $Global_depth{sample}{$sample}{median};
				my $norm_rate = $sample_rate / $run_rate;
				print RNORM "\t$norm_rate";
			}
			else {
				print RNORM "\tN/A";
			}
		}
		print RNORM "\n";
	}
	close RNORM;	
	close BED;
}

foreach my $sample ((split(/\s+/, $samples))) {
	open (PROBE, "> $sample\_stat.probes") or die "Can't open $sample\_stat.probes\n";
	open (BED, "< $probe_id") or die "Can´t open $probe_id\n";
	while (my $probe = <BED>) {
		chomp $probe;
		my ($chr, $start, $end, $gene, $name) = split (/\s+/, $probe);
		if ($chr =~ m/chr/i) {
			next;
		}
		print PROBE 	"$chr\t$start\t$end\t$gene\t$name\t$Global_depth{probe}{$sample}{$chr}{$start}{median}\t$Global_depth{probe}{$sample}{$chr}{$start}{mean}\n";	
	}
	
	close PROBE;
	close BED;
}

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
