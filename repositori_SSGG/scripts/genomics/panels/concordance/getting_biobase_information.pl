#!/usr/bin/perl -w

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

sub process ($$$$);
sub getsnps ($$$$);
sub getprobes ($$$);
sub getting_cov ($$);

my $output_id;
my $biobase_id;
my $probe_id;
my $threads = 2;
my $coverage_id;

GetOptions (	"b=s"	=> \$biobase_id, # input file
				"o=s"	=> \$output_id,	# output file root
				"p=s"	=> \$probe_id,	# probes file
				"t=i"   => \$threads, #number of processors		
				"c=s"	=> \$coverage_id #HapMap samples ID
);

print "\n\n=================================================================
Getting variant from global file to probe file script v 1.0\n";

if (not defined $probe_id or not defined $biobase_id or not defined $coverage_id) {
	die "
			Options:
		-b Biobase file root name (for example: /path_to/bioBaseSortUniques)
		-p Probes file name
		-c Coverage file root name
		-o Output file name root (optional)
		-t Number of processors (default = 2)
\n".localtime()."\n=================================================================\n\n";
}

print "\n Starting process... ".localtime()."\n\n";

if (not defined $output_id) {
	my @name = split (/\./, $coverage_id);
	$output_id = $name[0];
}

my $pm = new Parallel::ForkManager($threads); # Run in parallel 

$pm->run_on_finish( sub { my ($pid, $exit_code, $ident) = @_; } );
$pm->run_on_start( sub { my ($pid,$ident)=@_; } );
$pm->run_on_wait( sub {}, 0.5 );

open (SUM, "> $output_id\_biobase.sum") or die "Can't create $output_id\_biobase.sum\n";
print SUM "#Chr\tcovered\tnon_covered\n";
close SUM;

foreach my $chr (1..24) {
	my $pid = $pm->start($chr) and next;
	if ($chr == 23) { $chr = 'X';}
	elsif ($chr == 24) { $chr = 'Y';}
	my %Regions;
	my @Probe_Starts;
	my %Variants;
	getprobes($chr, \%Regions, \@Probe_Starts);
	getsnps($chr, \%Regions, \@Probe_Starts, \%Variants);
	my ($covered, $non_covered) = (0,0);
	process ($chr, \%Variants, \$covered, \$non_covered);
	open (SUM, ">> $output_id\_biobase.sum") or die "Can't create $output_id\_biobase.sum\n";
	print SUM $chr,"\t",
			  $covered ,"\t",
			  $non_covered ,"\n";
	close SUM;
	$pm->finish;
}
$pm->wait_all_children;

sub process ($$$$) {
  	my ($chr_id, $variants, $covered, $non_covered) = @_;
  	my $chr = $chr_id;
  	if ($chr_id eq 'X') { $chr = 23;}
  	elsif ($chr_id eq 'Y') { $chr = 24;}
	open (COV, "> $output_id\_biobase_covered.chr$chr_id") or die "Can't create $output_id\_covered.chr$chr_id\n";
	open (NON, "> $output_id\_biobase_non_covered.chr$chr_id") or die "Can't create $output_id\non_covered.chr$chr_id\n";
	print COV "chr\tposition\tgene\tvar_id\tdepth\n";
	print NON "chr\tposition\tgene\tvar_id\tdepth\n";
	print "Writing output files for chr$chr_id... ".localtime()."\n";
	
	foreach my $pos (sort {$$variants{$chr}{$a}<=>$$variants{$chr}{$b}} keys %{$$variants{$chr}}) {
		foreach my $name (sort {$$variants{$chr}{$pos}{$a} cmp $$variants{$chr}{$pos}{$b}} keys %{$$variants{$chr}{$pos}}) {
			my $cov = (&getting_cov ($chr_id, $pos) || '0');
			if ($cov < 10) {
				print NON "$chr_id\t$pos\t$$variants{$chr}{$pos}{$name}{gene}\t$name\t$cov\n";
				$$non_covered++;
			}
			else {
				print COV "$chr_id\t$pos\t$$variants{$chr}{$pos}{$name}{gene}\t$name\t$cov\n";
				$$covered++;
			}
		}
	}
	close COV;
	close NON;	
	close SUM;
}		

sub getting_cov ($$) {
	my ($chr, $pos) = @_;
	my $col = '$4';
	my $cov = `grep -m 1 "chr$chr\[\[:space:\]\]$pos\[\[:space:\]\]" $coverage_id.chr$chr.pileup | awk \'{print $col}\'`;
	chomp $cov;
	#print "$$chr\t$pos\t$cov\n";
	return $cov;
}

sub getsnps ($$$$) {
  	my ($chr_id, $regions, $starts, $variants) = @_;
  	my $chr = $chr_id;
  	if ($chr_id eq 'X') { $chr = 23;}
  	elsif ($chr_id eq 'Y') { $chr = 24;}
  	print "Getting variants for chr$chr_id... ".localtime()."\n";
	open (BIO, "< $biobase_id\_chr$chr_id.txt") or die "Can't open $biobase_id\_chr$chr_id.txt\n";
	my @data = @{$$starts[$chr]};
INDEL:while (my $variant = <BIO>) {
		if ($variant !~ m/^\d+/) { next;}
		chomp $variant;
		my ($chrom, $name, $gene, $hg18, $hg19) = split (/\t/, $variant);
		$hg19 =~ s/:/ /g;
		$hg19 =~ s/-/ /g;
		my $pos = (split (/\s+/, $hg19))[1];
		if ($chr eq 'X') { $chr = 23;}
		elsif ($chr eq 'Y') { $chr = 24;}		
	PROBE:foreach my $start (@data) {
			if ($pos >= $$regions{$chr}{probes}{$start}{start} and $pos <= $$regions{$chr}{probes}{$start}{end}){
				$$variants{$chr}{$pos}{$name}{gene} = $gene;
				last PROBE;
			}
 		}
	}
	close BIO;
}

sub getprobes ($$$) {
  	my ($chr_id, $regions, $probe_starts) = @_;
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
			push (@{$$probe_starts[$chr]}, $start);
		}
	}
#	print "Total number of probes after analysis: ".scalar(@{$starts})."\n";	
}

print "Process finished...".localtime()."\n";

exit;