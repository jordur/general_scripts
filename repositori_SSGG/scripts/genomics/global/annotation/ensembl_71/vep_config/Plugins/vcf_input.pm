=head1 LICENSE
                                                                                                                     
 Copyright (c) 2013 Sistemas Genómicos S.L - All rights reserved.                                                                      
                                                                                                        
 This software is distributed under a modified Apache license.                                                                                                                      
                                                                                                                     
=head1 CONTACT                                                                                                       

 Guillermo Marco Puche <guillermo.marco@sistemasgenomicos.com>
    
=cut

=head1 NAME

 vcf_input

=head1 SYNOPSIS

 VCF input file must be compressed and idexed with tabix before using the plugin:
 bgzip -c example.vcf > example.vcf.gz
 tabix -p vcf test.vcf.gz

 mv vcf_input.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin vcf_input,/path/to/file.vcf.gz

=head1 DESCRIPTION
 
 A VEP plugin that extracts the 8th column form VCF input file. (INFO)
 It extracts the corresponding 8th column corresponding to our own weird VCF
 format at writes to VEP output as "SAMPLES" column.

=cut

package vcf_input;

use strict;
use warnings;
use Data::Dumper;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub new {
    my $class = shift;
    
    my $self = $class->SUPER::new(@_);

	my $file = $self->params->[0];
	$self->{file} = $file;
	
	
	# test tabix
  	die "ERROR: tabix does not seem to be in your path\n" unless `which tabix 2>&1` =~ /tabix$/;
  	
  	# check files exist
	die "ERROR: Input VCF compressed & indexed file: $file not found\n" unless -e $file;
    die "ERROR: Tabix index file $file\.tbi not found - perhaps you need to create it first?\n" unless -e $file.'.tbi';
	

    return $self;
}

sub version {
    return '71';
}

sub feature_types {
    return ['Feature', 'Intergenic'];
}

sub get_header_info {
    return {
    	Chr => "Chromosome",
    	Pos => "Position",
    	Ref_Allele => "Reference Allele",
    	Var_Allele => "Variation Allele",
        Samples => "Columns containing Samples info: Genotype,Depth,Var/Depth, BIAS of VCF weird format.",
        Strand_Bias => "Strand Bias for each variant."
    };
}

sub run {
	
	my $self = shift;
    my $vf = shift;
    my $line_hash = shift;
    
    ## Def vars
	my ($samples,@samples,$sample_names,@sample_names,@samples_values,%res);
	my (@sample_space, @sample_comma, @sample_equal);
	my ($info, $chr, $pos, $ref_allele, $var_allele);
	my (@names, $strand_bias);
	my $pos_string;
    
    
    ## Def Handlers
    my $config = $self->{config};
    my $ind_cols = $config->{ind_cols};
    my $line = $vf->{base_variation_feature_overlap}->{base_variation_feature}->{_line};
    
    my @split_line = split ("\t", $line);
	#print "\nCHR: $split_line[0]\tPOS: $split_line[1]\tREF_ALLELE: $split_line[3]\tVAR:ALLELE: $split_line[4]\n";
		
	#We store the samples values information on $samples
	$samples = $split_line[7];
	$res{"Chr"} = $split_line[0];
	$res{"Pos"} = $split_line[1];
	$res{"Ref_Allele"} = $split_line[3];
	$res{"Var_Allele"} = $split_line[4];
	
	## Variable & handlers declaration
	#my $tva = @_;
	#my $vf = $tva->variation_feature;
	
	## Call get_samples_name() subroutine and obtain the samples name !
	## Remove *.gz extension from input file.
	my $vcf = substr($self->{file}, 0, -3);
	@names = &get_samples_name($vcf);

	#We split S1D=12;S1G=...
	@samples = split(";", $samples);
	
	#We obtain Strand Bias
	$strand_bias = (split("=",$samples[-1]))[1];
	
	#Reformat Strand Bias..
	#if ($strand_bias eq "---") {$strand_bias = "-";}
	 
	#We add Strand Bias to the hash %res
	$res{"Strand_Bias"} = $strand_bias; 
	#We remove Strand Bias from array
	splice (@samples, -1);

	#We store on @samples_values array the numeric value only.
	foreach my $sample_value (@samples)
	{
		push (@samples_values, (split("=", $sample_value))[1]);
	}
	
	#We obtain the number of samples.(will be used when Strand Bias will be fixed)
	my $num_muestras = scalar(@names);
	
	#For each sample we asign the correspoding value to the %res hash
	#Keys are sample column names ie: $col_name."_Depth
	#Values are obtained from @samples_values array.
	foreach my $sample_name (@names){
		
		#We asign the values to the hash and then we delete
		#the values from the @samples_values array.
		
		foreach my $i (@samples_values) {
			$res{$sample_name."_Depth"} = $samples_values[0];
			$res{$sample_name."_Var/Depth"} = $samples_values[1];
			$res{$sample_name."_Genotype"} = $samples_values[2];
			last;
			}
			
		#We remove the already asigned values.
		shift @samples_values for 1..3;
	}
	
    return { %res };
}

sub get_samples_name {
	my $vcf_input = $_[0];
	
	my (@names,@sample_space, @sample_comma, @sample_equal);
		
	#Open VCF file
	open VCF_INPUT, $vcf_input or die $!;
	
	#Extract the first non header line and get the samples
	#name and the INFO column to count the samples.
	while (<VCF_INPUT>){
		if ($_ =~ /^#Samples/) {
			@sample_space = (split(" ",$_))[1];
		}
	}
	
	close VCF_INPUT;
	
	foreach my $i (@sample_space){
		@sample_comma = split(";",$i);
	}
	foreach my $j (@sample_comma){
		push (@sample_equal, (split("=", $j))[1]);
	}

	@names = @sample_equal;
	return @names;
}
1;