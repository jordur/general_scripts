#!/usr/bin/env perl
########################################################
# Copyright 2012 - Sistemas Genomicos S.L.
# @Desc: Script to Parsing Somatic GATK in target regions
# @Author: Arbol & JM Rosa
# @Contributors: 
########################################################

use warnings;
use strict;
use Getopt::Long;
use Parallel::ForkManager;
use Cwd;
use File::Copy;
use File::Path;
use Statistics::Descriptive qw(:all);
use Data::Dumper;
use Tie::File;
use genomics::genotype;
use Scalar::Util qw(looks_like_number);

# Global parameters
my $F3_var_id = "";
my $F3_gatk_id = "";
my $pair_var_id;
my $pair_gatk_id;
my $F3_gff_id = "";
my $pair_gff_id = "";
my $output_id;

# Script parameters
GetOptions (	
				"vf=s"	=> \$F3_var_id,	
				"vp=s"	=> \$pair_var_id,	
				"gf=s"	=> \$F3_gatk_id,	
				"gp=s"	=> \$pair_gatk_id,
				"bf=s"	=> \$F3_gff_id,
				"bp=s"	=> \$pair_gff_id,					
				"o=s"	=> \$output_id,	
);

print "\n\n=================================================================
Collecting Somatic Vars from VarScan, GATK and BioScope Gff in target regions :: v 1.1\n";

if (not defined $F3_var_id or not defined $F3_gatk_id or not defined $pair_var_id or not defined $pair_gatk_id or not defined $F3_gff_id or not defined $pair_gff_id) {
	die "
			Options:
		-vf VarScan F3 file (optional)
		-vp VarScan pair file (mandatory)
		-gf GATK F3 file (optional)
		-gp GATK pair file (mandatory)
		-bf Bioscope GFF F3 file (optional)
		-bp Bioscope GFF pair file (optional)
		-o Output name (optional)
\n".localtime()."\n=================================================================\n\n";
}

&main();

sub main(){
	my $path = `pwd`;
	my %Variants;      # Variable for storing variants from all callers
	my %config;        # Variable for storing samples order and columns
	my $output_config; # Variable for storing output format
	chomp $path;
	 
	if (not defined $output_id) {
		$output_id = "varscan_gatk_gff_variants";
	}
	 
	# Open log, output and error files
	open (LOG, ">$output_id.log") or die "Can't create log file\n";
	print LOG "sg_collecting_varscan_gatk_somatic.pl -vf $F3_var_id -vp $pair_var_id -gf $F3_gatk_id -gp $pair_gatk_id -bf $F3_gff_id -bp $pair_gff_id -o $output_id in $path\n";
	print LOG "\nStarting process... ".localtime()."\n";
	open (OUT, ">$output_id.collected.vcf") or die "Can't create $output_id.scollected.vcf\n";
	open (ERR, ">$output_id.error") or die "Can't create $output_id.error\n";
	
	# Print header
	print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDATA\n";
	
	# Samples order in GATK Pair-Ends will be taken for output
	open (my $PGATK, "< $pair_gatk_id") or die "Can't open $pair_gatk_id\n";
	my $line = <$PGATK>;
	close $PGATK;
	$output_config = GetConfig($line);
	print OUT "#Samples: ";
	my @sample_ids = &get_ids($$output_config{samples});
	for (my $i=0;$i<=$#sample_ids;$i++){
		foreach(@sample_ids){
			if ($$output_config{samples}{$_}{order} == $i){
				print OUT "S" . ($i+1) . "=" . $_; 
			}
		}
		if ($i != $#sample_ids){
			print OUT ";";
		} else {
			print OUT "\n";
		}
	}
		
	# Check the different variant calling files. The earlier one variant calling file appears in the checking,
	# the higher its priority is
	open ($PGATK, "< $pair_gatk_id") or die "Can't open $pair_gatk_id\n";
	&RetrieveVarsInformation($PGATK, \%config, \%Variants, "gatk", "pair", $output_config);
	close $PGATK;
	
	open (my $VARSC, "< $pair_var_id") or die "Can't create $pair_var_id\n";
	&RetrieveVarsInformation($VARSC, \%config, \%Variants, "varscan", "pair", $output_config);
	close $VARSC;
	
	if ($F3_gatk_id ne ""){
		open (my $F3GATK, "< $F3_gatk_id") or die "Can't create $F3_gatk_id\n";
		&RetrieveVarsInformation ($F3GATK, \%config, \%Variants,"gatk","F3", $output_config);
		close $F3GATK;
	}
	
	if ($F3_var_id ne ""){
		open (my $F3VARSC, "< $F3_var_id") or die "Can't open $F3_var_id\n";
		&RetrieveVarsInformation ($F3VARSC, \%config, \%Variants,"varscan","F3", $output_config);
		close $F3VARSC;
	}
	
	if ($pair_gff_id ne ""){
		open (my $PAIRGFF, "< $pair_gff_id") or die "Can't open $pair_gff_id\n";
		&RetrieveVarsInformation ($PAIRGFF, \%config, \%Variants,"gff","pair", $output_config);
		close $PAIRGFF;
	}
	
	if ($F3_gff_id ne ""){
		open (my $F3GFF, "< $F3_gff_id") or die "Can't open $F3_gff_id\n";
		&RetrieveVarsInformation ($F3GFF, \%config, \%Variants,"gff","F3", $output_config);
		close $F3GFF;
	}
	
	# Sort and output variants
	foreach my $chr (sort {$a <=> $b} keys %Variants) {
		my $chrom;
		if ($chr eq "23"){
			$chrom = "chrX";
		} elsif ($chr eq "24"){
			$chrom = "chrY";
		} elsif ($chr eq "25"){
			$chrom = "chrM";
		} else {
			$chrom = "chr" . $chr;
		}

		# La linea comentada estaba anteriormente y no consigue ordenar efectivamente el hash
		#foreach my $position (sort {$Variants{$chr}{$a} <=>$Variants{$chr}{$b}} keys %{$Variants{$chr}}) {
		foreach my $position (sort {$a <=> $b} keys %{$Variants{$chr}}) {
			foreach my $ref (keys %{$Variants{$chr}{$position}}) {
				foreach my $var (keys %{$Variants{$chr}{$position}{$ref}}) {
					#my @data = split (/\s+/, $Variants{$chr}{$position}{$ref}{$var}{string});
					print OUT "$chrom\t$position\t.\t$ref\t$var\t1000\t.\t";
					my $sb = $Variants{$chr}{$position}{$ref}{$var}{SB};
					my $info;
					for(my $i=0;$i<=$#sample_ids;$i++){
						my $sample;
						foreach (@sample_ids){
							if ($$output_config{samples}{$_}{order} == $i){
								$sample = $_;
							}
						}
						$info .= "S".($i+1)."D=".$Variants{$chr}{$position}{$ref}{$var}{samples}{$sample}{depth}.";";
						$info .= "S".($i+1)."F=".$Variants{$chr}{$position}{$ref}{$var}{samples}{$sample}{freq}.";";
						$info .= "S".($i+1)."G=".$Variants{$chr}{$position}{$ref}{$var}{samples}{$sample}{geno}.";";
					}
					print OUT $info."SB=$sb\tGT:GK:VS:GF:PA:F3\t0/1:$Variants{$chr}{$position}{$ref}{$var}{gatk}:$Variants{$chr}{$position}{$ref}{$var}{varscan}:$Variants{$chr}{$position}{$ref}{$var}{gff}:$Variants{$chr}{$position}{$ref}{$var}{pair}:$Variants{$chr}{$position}{$ref}{$var}{F3}\n";
				}
			}
		}
	}
	
	close OUT;
	close ERR;
	
	print LOG "\nProcess finished... ".localtime()."\n";
	close LOG;
	
	print  "\nProcess finished... ".localtime()."\n";
}

sub CallerPriority(){
# Function that returns the priority of the used variant caller
	my ($variant_caller,$type) = @_;
	
	if ($variant_caller eq "gatk"){
		if ($type eq "pair"){
			return 6;
		} else {
			return 4;
		}
	} elsif ($variant_caller eq "varscan"){
		if ($type eq "pair"){
			return 5;
		} else {
			return 3;
		}
	} elsif ($type eq "pair"){
		return 2;
	} else {
		return 1;
	}
}

sub GetVariantCallPriority(){
# Function that returns the highest priority of the variant callers of a chosen variant
	my ($variant) = @_;
	my $priority = 0;
	
	if ($$variant{gatk} != 0){
		if ($$variant{pair} != 0){
			$priority = 6;
		} else {
			if ($priority < 4){
				$priority = 4;
			}
		}
	} elsif ($$variant{varscan} != 0){
		if ($$variant{pair} != 0){
			$priority = 5;
		} else {
			$priority = 3;
		}
	} elsif ($$variant{pair} != 0){
		$priority = 2;
	} else {
		$priority = 1;
	}
	
	return $priority;
}

sub ResetPreeminentCaller(){
# Function that resets variant-calls for a given position (variant)
	my ($variant) = @_;
	
	if ($$variant{gatk} != 0){
		$$variant{gatk} = 1;
	}
	if ($$variant{varscan} != 0){
		$$variant{varscan} = 1;
	}
	if ($$variant{gff} != 0){
		$$variant{gff} = 1;
	}
	if ($$variant{pair} != 0){
		$$variant{pair} = 1;
	}
	if ($$variant{F3} != 0){
		$$variant{F3} = 1;
	}
}

sub StoreVariant(){
# Function that inits a variants hash
	my ($Variants,$current_variant,$variant_caller,$type) = @_;
	
	my $chr = $$current_variant{chr};
	my $position = $$current_variant{position};
	my $ref = $$current_variant{ref};
	my $var = $$current_variant{var};
	
	if (not defined($$Variants{$chr}{$position}{$ref}{$var}{gatk})){
		$$Variants{$chr}{$position}{$ref}{$var}{gatk} = 0;
	}
	if (not defined($$Variants{$chr}{$position}{$ref}{$var}{varscan})){
		$$Variants{$chr}{$position}{$ref}{$var}{varscan} = 0;
	}
	if (not defined($$Variants{$chr}{$position}{$ref}{$var}{gff})){
		$$Variants{$chr}{$position}{$ref}{$var}{gff} = 0;
	}
	if (not defined($$Variants{$chr}{$position}{$ref}{$var}{pair})){
		$$Variants{$chr}{$position}{$ref}{$var}{pair} = 0;
	}
	if (not defined($$Variants{$chr}{$position}{$ref}{$var}{F3})){
		$$Variants{$chr}{$position}{$ref}{$var}{F3} = 0;
	}
	
	if (defined ($$Variants{$chr}{$position}{$ref}{$var}{SB})){
		# Strand bias value of the most valuable caller will remain
		if (&CallerPriority($variant_caller,$type) > &GetVariantCallPriority($$Variants{$chr}{$position}{$ref}{$var})){
			$$Variants{$chr}{$position}{$ref}{$var}{SB} = $$current_variant{SB};
		}
		# There is a change of preeminent caller...
		&ResetPreeminentCaller($$Variants{$chr}{$position}{$ref}{$var});		
	} else {
		$$Variants{$chr}{$position}{$ref}{$var}{SB} = $$current_variant{SB};
	}
	$$Variants{$chr}{$position}{$ref}{$var}{$variant_caller} = 2;
	$$Variants{$chr}{$position}{$ref}{$var}{$type} = 2;
	$$Variants{$chr}{$position}{$ref}{$var}{samples} = $$current_variant{samples};
}

sub ObtainFreqPos(){
# Function that obtains the relative position of a value in an array considering its value
	my ($array,$pos) = @_;
	
	my $position = 0;
	for (my $i=0;$i<=$#{$array};$i++){
		if($pos != $i and $$array[$pos] > $$array[$i]){
			$position++;
		}
	}
	return $position;
}

sub CheckVariantEquivalence(){
# Function that decides wether the information of a variant caller for a variant (variant) must be substituted with the information of another caller (current variant)
# Used criteria are: relative values of frequencies among samples and depth must be over half the depth of the preeminent caller
	my ($variant,$current_var) = @_;
	
	my @samples = get_ids($$current_var{samples});
	my @freqs_existent;
	my @freqs_current;
	my @pos_existent;
	my @pos_current;
	
	# Create vector containing the allelic frequencies for the variants
	foreach my $sample (@samples){
		push(@freqs_existent,$$variant{samples}{$sample}{freq});
		push(@freqs_current,$$current_var{samples}{$sample}{freq});
		
		# Check that variant depth is over half the depth of the preeminent caller for every sample
		my $var_depth_normalized;
		if (Scalar::Util::looks_like_number($$variant{'samples'}{$sample}{'depth'})){
			$var_depth_normalized = $$variant{'samples'}{$sample}{'depth'};
		} else {
			$var_depth_normalized = 0;
		}
		my $current_var_depth_normalized;
		if (Scalar::Util::looks_like_number($$current_var{'samples'}{$sample}{'depth'})){
			$current_var_depth_normalized = $$current_var{'samples'}{$sample}{'depth'};
		} else {
			$current_var_depth_normalized = 0;
		}
				
		if ($current_var_depth_normalized < $var_depth_normalized/2){
			return 0;
		}
	}
	
	# Create vector containing relative position of the allelic frequencies among samples
	for (my $i=0;$i<=$#samples;$i++){
		push(@pos_existent,&ObtainFreqPos(\@freqs_existent,$i));
		push(@pos_current,&ObtainFreqPos(\@freqs_current,$i));	
	}
		
	for (my $i=0;$i<=$#pos_existent;$i++){
		if ($pos_existent[$i] != $pos_current[$i]){
			return 0;
		}
	}

	# If all criteria are satisfied, return OK=1
	return 1;
}

sub GetSamplesInformation(){
# Function for comparing current variant information with variant information from current chosen variant caller
# Finally, only information from chosen variant caller will remain in the output file
	my ($config, $variants, $line, $variant_caller, $type) = @_;
	
	my %current_var = &GetFields($config,$line);
	my $chr = $current_var{chr};
	my $position = $current_var{position};
	my $ref = $current_var{ref};
	my $var = $current_var{var};
	
	# If a variant is already existent (already called with a preeminent caller) and its predicted genotype for a sample is Homo_ref or P_Homo_ref (freq < 0.12), 
	# then it will be checked if current genotype for the sample in other variant calling strategies has a greater freq. If so, information from this non-preeminent
	# caller is included in the collected vcf. If not, information from preeminent caller remains for the variant
	if (defined($$variants{$chr}{$position}{$ref}{$var})){
		# Variant has been detected by caller, although it's still unknown its preevalency
		if ($$variants{$chr}{$position}{$ref}{$var}{$variant_caller} == 0){
			$$variants{$chr}{$position}{$ref}{$var}{$variant_caller} = 1;
		}
		if ($$variants{$chr}{$position}{$ref}{$var}{$type} == 0){
			$$variants{$chr}{$position}{$ref}{$var}{$type} = 1;
		}
		
		# Check preevalency of the variant (by checking all samples)
		my @samples = &get_ids($$variants{$chr}{$position}{$ref}{$var}{samples});
		foreach (@samples){
			my $sample = $_;
			if ($$variants{$chr}{$position}{$ref}{$var}{samples}{$sample}{freq} < $current_var{samples}{$sample}{freq} and $$variants{$chr}{$position}{$ref}{$var}{samples}{$sample}{freq} < 0.12){
				if (&CheckVariantEquivalence($$variants{$chr}{$position}{$ref}{$var},\%current_var)){
					&StoreVariant($variants,\%current_var,$variant_caller,$type);
				}
			}
		}
	} else {
		&StoreVariant($variants,\%current_var,$variant_caller,$type);
	}
	#if $$variants{}
}

sub GetConfig(){
# Function that returns the configuration of an snvs file
	my $header = shift;
	
	my %config;
	my @data = split(/\s+/, $header);
	my $nsamp = 0;

	# Static values
	$config{chr} = 0;
	$config{position} = 1;
	$config{ref} = 2;
	$config{var} = 3;
	$config{type} = 4;
	$config{SB} = $#data;
	
	# Sample values
	for(my $i=0;$i<=$#data;$i++){
		my $tag;
		if (length($data[$i]) > 5){
			$tag = substr($data[$i],length($data[$i])-6,6)	
		} else {
			$tag = "";
		}
		if($tag eq "_depth"){
			my $sample_name = substr($data[$i],0,length($data[$i])-6);
			$config{samples}{$sample_name}{depth} = $i;
			$config{samples}{$sample_name}{geno} = $i + 1;
			$config{samples}{$sample_name}{freq} = $i + 2;
			$config{samples}{$sample_name}{order} = $nsamp;
			$nsamp ++;
			$i = $i + 2;
		}
	}
	return \%config;
}

sub GetFields(){
# Function for getting fields from a snvs file line, once its configuration has been checked
	my ($config, $line) = @_;
	
	my @data = split (/\s+/, $line);
	
	my %variants;
	my @variant_ids = get_ids($config);
	foreach (@variant_ids){
		my $id = $_;
		if ($id eq "chr"){
			my $chr = $data[0];
			$chr =~ s/chr//;
			$chr =~ s/X/23/;
			$chr =~ s/Y/24/;
			$chr =~ s/M/25/;
			if ($chr =~ m/\\0\\0\\/) {
				print "$line\n";
			}
			$variants{$id} = $chr;
		} elsif ($id ne "samples"){
			$variants{$id} = $data[$$config{$id}];
		} else {
			my @sample_ids = get_ids($$config{samples});
			foreach (@sample_ids){
				my $sample_id = $_;
				$variants{samples}{$sample_id}{depth} = $data[$$config{samples}{$sample_id}{depth}];
				$variants{samples}{$sample_id}{geno} = $data[$$config{samples}{$sample_id}{geno}];
				my $freq;
				if ($data[$$config{samples}{$sample_id}{freq}] eq "-"){
					$freq = 0;
				} else {
					$freq = $data[$$config{samples}{$sample_id}{freq}];
				}
				$variants{samples}{$sample_id}{freq} = $freq;
				# Genotype is once again predicted here. In future versions, a more sophisticated prediction must be carried out
				my @tmp = split(",", genotype->get_genotype ($data[$$config{samples}{$sample_id}{depth}],$data[$$config{samples}{$sample_id}{freq}],$data[$$config{ref}],$data[$$config{var}]));
				$variants{samples}{$sample_id}{geno} = $tmp[0];
			}
		}
	}
	
	return %variants;
}

sub get_ids ()
{
	my ($self) = @_;
	return keys (%{$self});
}

sub CompareConfigs()
# Function for comparing 2 files configurations, checking sample names consistency
{
	my ($conf_input, $conf_output) = @_;
	my @input_samples = &get_ids($$conf_input{samples});
	foreach (@input_samples){
		my $input_sample = $_;
		my $found = 0;
		my @output_samples = &get_ids($$conf_output{samples});
		foreach (@output_samples){
			if ($_ eq $input_sample){
				return 0;
			}
		}
	}
	return 1;
}

sub RetrieveVarsInformation(){
# Function that returns a hash with the information of the variants of an snvs file
	my ($file_handle, $config, $Variants, $variant_caller, $type, $output_config) = @_;
	while (my $line = <$file_handle>){
		chomp $line;
		if ($line =~ m/^#Chr/) {
			# Get configuration of current variants file
			$config = &GetConfig($line);
			
			# If current variants file configuration is different from output configuration, exit the programm
			if (&CompareConfigs($config,$output_config)){
				print Dumper($config);
				print Dumper($output_config);
				die ("ERROR: Inconsistent information in $variant_caller - $type!!!\n"); 
			}
		}
		else {
			&GetSamplesInformation($config,$Variants,$line,$variant_caller,$type);
		}
	} 
}

exit;