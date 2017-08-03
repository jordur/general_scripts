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

sub process ($);
sub getsnps ($$$$);
sub getprobes ($$$);
sub getpos ($$);
sub getfreq ($);
sub getgenotype ($$$);

my $output_id;
my $input_id;
my $probe_id;
my $threads = 2;
my $samples_id;

GetOptions (	"i=s"	=> \$input_id, # input file
				"o=s"	=> \$output_id,	# output file root
				"p=s"	=> \$probe_id,	# probes file
				"t=i"   => \$threads, #number of processors		
				"n=s"	=> \$samples_id #HapMap samples ID
);

print "\n\n=================================================================
Getting variant from global file to probe file script v 1.0\n";

if (not defined $probe_id or not defined $input_id or not defined $samples_id) {
	die "
			Options:
		-i Variant input file root (example ALL_phase1_release_v3.20101123.snps_indels_svs.genotypes)
		-p Probes file name
		-n HapMap IDs (comma separated)
		-o Output file name root (optional)
		-t Number of processors (default = 2)
\n".localtime()."\n=================================================================\n\n";
}

print "\n Starting process... ".localtime()."\n\n";

if (not defined $output_id) {
	my @name = split (/\./, $input_id);
	$output_id = $name[0]."_final";
}

$threads = 8 if ($threads > 8);

 my $pm = new Parallel::ForkManager($threads); # Run in parallel with 8 threads max
 $pm->run_on_finish( sub { my ($pid, $exit_code, $ident) = @_; } );
 $pm->run_on_start( sub { my ($pid,$ident)=@_; } );
 $pm->run_on_wait( sub {}, 0.5 );

foreach my $chr (1..23) {
	my $pid = $pm->start($chr) and next;
	if ($chr == 23) { $chr = 'X';}
	elsif ($chr == 24) { $chr = 'Y';}
	my %Regions;
	my @Probe_Starts;
	my %Indels;
	getprobes($chr, \%Regions, \@Probe_Starts);
	#getsnps($chr, \%Regions, \@Probe_Starts, \%Indels);
	process ($chr);
	$pm->finish;
}
$pm->wait_all_children;

sub process ($) {
  	my ($chr_id) = @_;
	my @samples = split (/,/, $samples_id);
	foreach my $sample (@samples) {
		open (OUT, "> ${sample}_${output_id}.chr$chr_id.snv") or die "Can't create ${sample}_${output_id}.chr$chr_id.snv\n";
		print "Writing output file for chr$chr_id... ".localtime()."\n";
		print OUT "#CHR\tPOSITION\tREF_AL\tVAL_AL\tVAL_AL_FREQ\tSNP_ID\t";
		my $first_line = `gunzip -c $input_id.chr$chr_id.vcf.gz | grep -m 1 \"^#CHROM\"`;
		chomp $first_line;
		print OUT "$sample\n";
		my $pos = &getpos ($sample, $first_line);
		open (TMP, "< $output_id.$chr_id.tmp") or die "Can't open $output_id.$chr_id.tmp\n";
		while (my $line = <TMP>) {
			my @data = split(/\s+/, $line);
			my $allele_freq = &getfreq($data[7]);
			my $genotype = &getgenotype($data[3], $data[4], $data[$pos]);
			print OUT "$data[0]\t$data[1]\t$data[3]\t$data[4]\t$allele_freq\t$data[2]\t$genotype\n";
		}
		close OUT;
	}
}		

sub getgenotype ($$$) {
	my ($ref, $var, $info) = @_;
	my $geno = (split (/:/, $info))[0];
	$geno =~ s/\|/\//g;
	$geno =~ s/0/$ref/g;
	$geno =~ s/1/$var/g;
	#print "$info\t$geno\n";
	return $geno;
}

sub getfreq ($) {
	my $info = shift;
	my @data = split (/;/, $info);
	foreach my $field (@data) {
		if ($field =~ m/AF=/ && $field !~ m/LD/) {
			$field =~ s/AF=//g;
			#print "$field - $info\n";
			return $field;
		}
	}
}

sub getpos ($$) {
	my ($id, $line) = @_;
	my @data = split(/\s+/,$line);
	
	for (my $i = 0; $i <= scalar(@data); $i++){
		if ($data[$i] eq $id) {
			return $i;
			last;
		}
	}
}

sub getsnps ($$$$) {
  	my ($chr_id, $regions, $starts, $indels) = @_;
  	my $chr = $chr_id;
  	if ($chr_id eq 'X') { $chr = 23;}
  	elsif ($chr_id eq 'Y') { $chr = 24;}
  	print "Getting variants for chr$chr_id... ".localtime()."\n";
	open (VCF, "gunzip -c $input_id.chr$chr_id.vcf.gz|") or die "Can't open $input_id.chr$chr_id.vcf.gz\n";
	open (OUT, "> $output_id.$chr_id.tmp") or die "Can't create $output_id.$chr_id.tmp\n";
	my @data = @{$$starts[$chr]};
INDEL:while (my $variant = <VCF>) {
		if ($variant =~ m/^#/) { next;}
		chomp $variant;
		my ($chr, $pos) = split (/\s+/, $variant);
		if ($chr eq 'X') { $chr = 23;}
		elsif ($chr eq 'Y') { $chr = 24;}		
	PROBE:foreach my $start (@data) {
			if ($pos > $$regions{$chr}{probes}{$start}{start} and $pos < $$regions{$chr}{probes}{$start}{end}){
				print OUT "$variant\n";
				last PROBE;
			}
 		}
	}
	close VCF;
	close OUT;
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

exit;

