#!/usr/bin/env perl
########################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Template script for developing Perl scripts
# @Author: Arbol
# @Contributors: 
########################################################

# -----------------------
# -- Include libraries --
# -----------------------

use warnings;
use strict;
use File::Basename;
use Getopt::Long;
use Cwd;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use file_handle;

# ---------------------------
# -- functions definitions --
# ---------------------------


# ---------------------------------
# -- global & default parameters --
# ---------------------------------
my ($f3fasta,$f3qual,$f5fasta,$f5qual,$r1fastq,$r2fastq);
my $output = cwd() . "/output_file.tsv";

# --------------------------------
# -- Input parameters & options --
# --------------------------------
(my $basename, my $dir, my $ext) = fileparse($0, qr/\.[^.]*/);

GetOptions (
				"f3f=s"	=> \$f3fasta,
				"f3q=s"	=> \$f3qual,
				"f5f=s"	=> \$f5fasta,
				"f5q=s"	=> \$f5qual,
				"r1=s"	=> \$r1fastq,
				"r2=s"	=> \$r2fastq,
				"o=s"	=> \$output,
);

print "\n=================================================================
$basename: Script for obtaining sequencing lane metrics v0.1\n";

if ((not defined($f3fasta) or not defined($f3qual) or not defined($f5fasta) or not defined($f5qual)) and (not defined($r1fastq) or not defined($r2fastq))) {
	die "
			Parameters:
	For fasta+qual files:
		-f3f  F3 fasta or csfasta file
		-f3q  F3 qual file
		-f5f  F5 fasta or csfasta file
		-f5q  F5 qual file
	For fastq files:
		-r1   Read group 1 fastq.gz file (gzipped fastq file)
		-r2   Read group 2 fastq.gz file (gzipped fastq file)
	General parameters:
		-o    Output name (optional)
\n".localtime()."\n=================================================================\n\n";
}

&main ();
print  "\nINFO: Process successfully finished... ".localtime()."\n";

# ---------------
# -- functions --
# ---------------
sub main (){

	my $path = `pwd`;
	chomp $path;
	my (%f3fmetrics,%f5fmetrics,%f3qmetrics,%f5qmetrics,%r1metrics,%r2metrics);
	
	# Print output header
	my $out = &get_out_file_handle($output);
	print $out "File_name\tRead_length\tNumber_of_reads\n";
	
	# Open files and read their metrics
	if (defined($f3fasta)){
		%f3fmetrics = GetMetrics($f3fasta,"fasta");
		print $out $f3fasta."\t".$f3fmetrics{read_length}."\t".$f3fmetrics{reads}."\n";
	}
	if (defined($f3qual)){
		%f3qmetrics = GetMetrics($f3qual,"qual");
		print $out $f3qual."\t".$f3qmetrics{read_length}."\t".$f3qmetrics{reads}."\n";
	}
	if (defined($f5fasta)){
		%f5fmetrics = GetMetrics($f5fasta,"fasta");
		print $out $f5fasta."\t".$f5fmetrics{read_length}."\t".$f5fmetrics{reads}."\n";
	}
	if (defined($f5qual)){
		%f5qmetrics = GetMetrics($f5qual,"qual");
		print $out $f5qual."\t".$f5qmetrics{read_length}."\t".$f5qmetrics{reads}."\n";
	}
	if (defined($r1fastq)){
		%r1metrics = GetMetrics($r1fastq,"fastq");
		print $out $r1fastq.": sequences"."\t".$r1metrics{seqs}{read_length}."\t".$r1metrics{seqs}{reads}."\n";
		print $out $r1fastq.": qualities"."\t".$r1metrics{quals}{read_length}."\t".$r1metrics{quals}{reads}."\n";
	}
	if (defined($r2fastq)){
		%r2metrics = GetMetrics($r2fastq,"fastq");
		print $out $r2fastq.": sequences"."\t".$r2metrics{seqs}{read_length}."\t".$r2metrics{seqs}{reads}."\n";
		print $out $r2fastq.": qualities"."\t".$r2metrics{quals}{read_length}."\t".$r2metrics{quals}{reads}."\n";
	}
	
	if (defined($f3fasta) and defined($f3qual)){
		if ($f3fmetrics{reads} != $f3qmetrics{reads}){
			die "ERROR: Not concordant number of reads in F3 fasta/csfasta file (".$f3fmetrics{reads}." reads) and F3 qual file (".$f3qmetrics{reads}." reads)\n";
		}
	}
	if (defined($f5fasta) and defined($f5qual)){
		if ($f5fmetrics{reads} != $f5qmetrics{reads}){
			die "ERROR: Not concordant number of reads in F5 fasta/csfasta file (".$f5fmetrics{reads}." reads) and F5 qual file (".$f5qmetrics{reads}." reads)\n";
		}
	}
	if (defined($r1fastq) and defined($r2fastq)){
		if ($r1metrics{seqs}{reads} != $r1metrics{quals}{reads}){
			die "ERROR: Not concordant number of reads in R1 fastq file (".$r1metrics{seqs}{reads}." sequences and ".$r1metrics{quals}{reads}." quality values)\n";
		}
		if ($r2metrics{seqs}{reads} != $r2metrics{quals}{reads}){
			die "ERROR: Not concordant number of reads in R1 fastq file (".$r2metrics{seqs}{reads}." sequences and ".$r2metrics{quals}{reads}." quality values)\n";
		}
	}
	
	# Close output file handle
	$out->close();
}

sub GetMetrics(){
	my ($file,$type) = @_;
	my %metrics;
	my $fh;
	
	# Open file for reading
	# In case of fastq, it's assumed that is a gzipped fastq file (.fast.gz)
	if ($type eq "fastq"){
		$fh = IO::Uncompress::Gunzip->new( $file )
    		or die "ERROR: IO::Uncompress::Gunzip failed: $GunzipError\n";
	} else {
		$fh = &get_in_file_handle($file);
	}
	
	# Check each line
	my $reads = 0; # Number of reads in file
	my $fastq_info = "";
	my $fastq_qualreads =0; # Number of qual reads in fastq files only
	my $nts = 0; # Number of nucleotides per read
	my $fastq_qualnts = 0; # Number of qual values per read in case of fastq files
	my $current_read;
	my $all_nts; # total number of nucleotides in the file
	my $all_qual_vals; # total number of quality values in the file
	my $unconstant = 0; # for checking if all reads have same length
	while (<$fh>){
		my $line = $_;
		chomp ($line);
		if ($type eq "fasta" or $type eq "qual"){
			if (substr($line,0,1) eq ">"){
				$current_read = $line;
				$reads++;
			} else {
				my $current_nts;
				if ($type eq "fasta"){
					$current_nts = length($line)-1;
					$all_nts += $current_nts;
				} elsif ($type eq "qual"){
					my @quals = split(" ",$line);
					$current_nts = $#quals+1;
					$all_qual_vals += $current_nts;
				}
				if ($nts == 0){
					$nts = $current_nts;
				} elsif ($current_nts != $nts){
					$unconstant = 1;
				}
			}
		} elsif ($type eq "fastq"){
			$current_read = $line;
			$reads++;
			my $line = <$fh>;
			chomp ($line);
			my $current_nts;
			$current_nts = length($line)-1;
			$all_nts += $current_nts;
			if ($nts == 0){
				$nts = $current_nts;
			} elsif ($current_nts != $nts){
				$unconstant = 1;
			}
			$line = <$fh>;
			chomp ($line);
			$fastq_qualreads++;
			$line = <$fh>;
			chomp ($line);
			$current_nts = length($line)-1;
			$all_qual_vals += $current_nts;
			if ($fastq_qualnts == 0){
				$fastq_qualnts = $current_nts;
			} elsif ($current_nts != $fastq_qualnts){
				$unconstant = 1;
			}
		} else {
			die "ERROR: Invalid file type: $type!!\n";
		}
	}
	$fh->close();
	if ($unconstant){
		print "WARNING: Not constant number of nucleotides per read in file $file\n";
		$nts = $all_nts / $reads;
		if ($fastq_qualreads != 0){
			$fastq_qualnts = $all_qual_vals / $fastq_qualreads;
		}
	}
	
	if ($type eq "fastq"){
		$metrics{seqs}{read_length} = $nts;
		$metrics{quals}{read_length} = $fastq_qualnts;
		$metrics{seqs}{reads} = $reads;
		$metrics{quals}{reads} = $fastq_qualreads;
	} else {
		$metrics{read_length} = $nts;
		$metrics{reads} = $reads;
	}
	return %metrics;
}

exit;