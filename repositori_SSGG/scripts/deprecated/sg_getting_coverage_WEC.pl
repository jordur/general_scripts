#!/usr/bin/perl

#=================================================#
# Getting coverage from global file to probe file #
#=================================================#


use warnings;
use strict;
use Getopt::Long;
use Parallel::ForkManager;
use Cwd;
use File::Copy;
use Statistics::Descriptive qw(:all);
use Data::Dumper;
use Tie::File;

sub process ($$$);
sub getcov ($$$$);
sub getprobes ($$$$);

my $output_id;
my $probe_id;
my $coverage_id;
my $threads = 2;

GetOptions (	"o=s"	=> \$output_id,	# output file root
				"p=s"	=> \$probe_id,	# probes file
				"c=s" 	=> \$coverage_id, # coverage file root
				"t=i"   => \$threads, #number of processors
);

print "\n\n=================================================================
Getting coverage from global file to probe file script v 1.0\n";

if (not defined $probe_id or not defined $coverage_id) {
	die "
			Options:
		-p Probes file name
		-c Coverage file name 
		-o Output file name root (optional)
		-t Number of processors (default = 2)
\n".localtime()."\n=================================================================\n\n";
}

print "\n Starting process... ".localtime()."\n\n";

if (not defined  $output_id ) {
	$output_id = $coverage_id;
}

$threads = 8 if ($threads > 8);

 my $pm = new Parallel::ForkManager($threads); # Run in parallel with 8 threads max
 $pm->run_on_finish( sub { my ($pid, $exit_code, $ident) = @_; } );
 $pm->run_on_start( sub { my ($pid,$ident)=@_; } );
 $pm->run_on_wait( sub {}, 0.5 );

my @global_depth;

foreach my $chr (1..24) {
	my $pid = $pm->start($chr) and next;
	if ($chr == 23) { $chr = 'X';}
	elsif ($chr == 24) { $chr = 'Y';}
	my %Regions;
	my @Probe_Starts;
	my @CNV_Starts;
	getprobes($chr, \%Regions, \@Probe_Starts, \@CNV_Starts);
	getcov($chr, \%Regions, \@Probe_Starts, \@CNV_Starts);
	process ($chr, \%Regions, \@Probe_Starts);
	$pm->finish;
}
$pm->wait_all_children;

sub process ($$$) {
  	my ($chr_id, $regions, $starts) = @_;
  	my $chr = $chr_id;
  	if ($chr_id eq 'X') { $chr = 23;}
  	elsif ($chr_id eq 'Y') { $chr = 24;}
	open (STAT, "> $output_id.stat.$chr_id") or die "Can't create $output_id.stat.$chr_id\n";
	
	print "Writing output file for chr$chr_id... ".localtime()."\n";
	
	PROBE:foreach my $start (@{$starts}) {
#		print "$start\n";
		my @Depth;
		foreach my $depth (keys %{$$regions{$chr}{probes}{$start}{Depth}}) {
			for (my $i = 0; $i < $$regions{$chr}{probes}{$start}{Depth}{$depth}{value}; $i++) {
				push (@Depth, $depth);
			}
		}
		my $length = $$regions{$chr}{probes}{$start}{end} - $$regions{$chr}{probes}{$start}{start};
		my $diff = $length - scalar (@Depth);
		my $coverage = scalar (@Depth) / $length;

		if ($diff < -1) {
			print "Something weird happens in chr$chr_id:$start-$$regions{$chr}{probes}{$start}{end}:\n$diff -> $length - ".scalar (@Depth)."\n";
		}
		elsif ($diff > 0) {
			for (my $j = 0; $j < $diff; $j++) {
				push (@Depth, 0);
			}
		}
		my $stat = Statistics::Descriptive::Full->new();
		$stat->add_data(@Depth);
		my $median = $stat->median();
		my @mean = split (/\./, $stat->mean());
		my $max = $stat->max();
		my $min = $stat->min();
		if (not defined $max and not defined $min) {
			$max = $min = 0;
		}
		print STAT "$chr_id\t$$regions{$chr}{probes}{$start}{start}\t$$regions{$chr}{probes}{$start}{end}\t$median\t$mean[0]\t$min\t$max\t$coverage\n";
	}
	close STAT;
}		

sub getcov ($$$$) {
  	my ($chr_id, $regions, $starts) = @_;
  	my $chr = $chr_id;
  	if ($chr_id eq 'X') { $chr = 23;}
  	elsif ($chr_id eq 'Y') { $chr = 24;}
  	print "Getting coverage from $coverage_id.$chr_id... ".localtime()."\n";
	open (OUT, "> $output_id.cov.$chr_id") or die "Can't create $output_id.$chr_id\n";
	open (COV, "< $coverage_id.$chr_id") or die "Can't open $coverage_id.$chr_id \n";
	my @data = @{$starts};
LINE:while (my $line = <COV>) {
		chomp $line;
		my ($chrom, $pos, $base, $depth) = split (/\s+/, $line);
		PROBE:foreach my $start (@data) {
			my $end = $$regions{$chr}{probes}{$start}{end};
			if ($pos < $start) {
				next LINE;
			}
			elsif ($pos > $end) {
				shift @data;
				next PROBE;
			}
			elsif ($pos >= $start && $pos <= $end) {
				print OUT "$line\n";
				if (not defined $$regions{$chr}{probes}{$start}{Depth}{$depth}{value}) {
					$$regions{$chr}{probes}{$start}{Depth}{$depth}{value} = 1;
					last PROBE;
					
				}
				else {
					$$regions{$chr}{probes}{$start}{Depth}{$depth}{value}++;
					last PROBE;
				}	
			}
		}
	}
	close COV;
	close OUT;
}

sub getprobes ($$$$) {
  	my ($chr_id, $regions, $probe_starts, $cnv_starts) = @_;
  	my $chr = $chr_id;
  	if ($chr_id eq 'X') { $chr = 23;}
  	elsif ($chr_id eq 'Y') { $chr = 24;}
  	print "Getting probes for chr$chr_id... ".localtime()."\n";
	my @regs = `grep ^chr$chr_id\[[:space:]] $probe_id`;
#	print "Total number of probes: ".scalar(@regs)."\n";
	foreach my $probe (@regs) {
		chomp $probe;
		my ($chr, $start, $end, $name) = split (/\s+/, $probe);
		if (not defined $name) { $name = "chr$chr:$start-$end"};
		my @name = split (/\|/, $name);
		my $probe_name = pop @name;
#		print "$name\n";
		$chr =~ s/chr//;
		if ($chr eq 'X') { $chr = 23;}
		elsif ($chr eq 'Y') { $chr = 24;}
		if (defined $$regions{$chr}{probes}{$start}) {
			next;
		}
		else {			
	#		print "$probe_name\n";
			$$regions{$chr}{probes}{$start}{start} = $start;
			$$regions{$chr}{probes}{$start}{end} = $end;
			$$regions{$chr}{probes}{$start}{name} = $probe_name;
			push (@{$probe_starts}, $start);
		}
	}
#	print "Total number of probes after analysis: ".scalar(@{$starts})."\n";	
}

exit;

