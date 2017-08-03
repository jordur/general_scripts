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

sub process ($$);

my $input_id;
my $probe_id;
my $threads = 2;
my $output_id;

GetOptions (	"o=s"	=> \$output_id,	# output file root
				"p=s"	=> \$probe_id,	# probes file
				"t=i"   => \$threads, #number of processors
				"i=s"	=> \$input_id, # rate thresholds
);

print "\n\n=================================================================
Getting reads on target using bamtools script v 1.0\n";

if (not defined $probe_id or not defined $input_id) {
	die "
			Options:
		-p Probes file name
		-i Input BAM file name 
		-o Output file name root (optional)
		-t Number of processors (default = 2)
\n".localtime()."\n=================================================================\n\n";
}

print "\n Starting process... ".localtime()."\n\n";

if (not defined  $output_id ) {
	my @name = split (/\./,$input_id); 
	$output_id = $name[0];
}

$threads = 8 if ($threads > 8);

 my $pm = new Parallel::ForkManager($threads); # Run in parallel with 8 threads max
 $pm->run_on_finish( sub { my ($pid, $exit_code, $ident) = @_; } );
 $pm->run_on_start( sub { my ($pid,$ident)=@_; } );
 $pm->run_on_wait( sub {}, 0.5 );

open (OUT, "> $output_id.stat.all") or die "Can't create $output_id.stat\n";

foreach my $chr_id (1..24) {
	my $pid = $pm->start($chr_id) and next;
	my $chr = $chr_id;
	if ($chr_id == 23) { $chr = 'X';}
	elsif ($chr_id == 24) { $chr = 'Y';}
	my @chr_reads;

	# chr_reads must be initialized!!
	for (my $i=1;$i<25;$i++){
		$chr_reads[$i]=0;
	}
	process($chr, \@chr_reads);
	print OUT "$chr\t$chr_reads[$chr_id]\n";
	$pm->finish;
}
$pm->wait_all_children;
close OUT;

sub process ($$) {
  	my ($chr_id, $chr_reads) = @_;
  	my $chr = $chr_id;
  	
  	print "Getting probes for chr$chr_id... ".localtime()."\n";
	my @regs = `grep ^chr$chr_id\[[:space:]] $probe_id`;
	open (PROB, "> $output_id.probe_stat.$chr_id") or die "ERROR: Can't create $output_id.probe_stat.$chr_id\n";
	foreach my $probe (@regs) {
		chomp $probe;
		my ($chr, $start, $end, $name) = split (/\s+/, $probe);
		$chr =~ s/chr//;
		my $region = "chr".$chr.":".$start."..".$end;
		my $reads = `bamtools count -in $input_id -region $region`;
		chomp $reads;
		if ($chr eq 'X') { $chr = 23;}
		elsif ($chr eq 'Y') { $chr = 24;}
		$$chr_reads[$chr] += $reads;
		print PROB "$chr\t$start\t$end\t$reads\n";
	}
	print "Analysis finished for chr$chr_id\n";
	close PROB;
}

exit;

