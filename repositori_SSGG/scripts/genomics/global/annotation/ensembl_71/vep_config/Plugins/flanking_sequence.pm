=head1 LICENSE
                                                                                                                                                                                       
                                                                                                                     
 This software is distributed under a modified Apache license.                                                                                                                      
                                                                                                                     
=head1 CONTACT                                                                                                       

 Guillermo Marco Puche <guillermo.marco@sistemasgenomicos.com>
 
=cut

=head1 NAME

 flanking_sequence

=head1 SYNOPSIS

perl variant_effect_predictor.pl -i variations.vcf --plugin flanking_sequence

=head1 DESCRIPTION
 
 Obtains the flanking sequence on F5-F3 strands around the variation location
 (50 bp each side)

=cut

package flanking_sequence;

use strict;
use warnings;

use Data::Dumper;
use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub new {
    my $class = shift;
    
    my $self = $class->SUPER::new(@_);
	

    return $self;
}

sub version {
    return '71';
}

sub feature_types {
    return ['Feature'];
    #return ['Feature', 'Intergenic'];
}


sub get_header_info {
    return {
        Flanking_sequence => "Column containing the flanking sequence. 5' to 3' 50 bases on each side.",
    };
}



sub run {
	
	## Objetos que se heredan
	my $self = shift;
	#my $tva = shift;
    my $vf = shift;
	
	## Def Handlers
    my $config = $self->{config};
    #my $ind_cols = $config->{ind_cols};
    my $line = $vf->{base_variation_feature_overlap}->{base_variation_feature}->{_line};
	
	## Variables
	my ($five_prime_seq, $three_prime_seq, $ref_allele, $var_allele) = "";
	
	## Split line
    my @split_line = split ("\t", $line);

	## Store the allele values information	
	$ref_allele = $split_line[3];
	$var_allele = $split_line[4];

	my $vf2 = $vf->variation_feature;

	## Will McLaren fix
	if(!defined($vf2->{slice})) {
		my $sa = $config->{sa} || Bio::EnsEMBL::Registry->get_adaptor($config->{species}, 'core', 'slice');
		$vf2->{slice} = $sa->fetch_by_region('chromosome', $vf2->{chr});
	}
	
	if (defined $vf2->feature_Slice->expand(50, -1)->seq){
		$five_prime_seq = $vf2->feature_Slice->expand(50, -1)->seq;	
	}
	
	if (defined $vf2->feature_Slice->expand(-1, 50)->seq){
		$three_prime_seq = $vf2->feature_Slice->expand(-1, 50)->seq;
	}

	my $flanking_sequence = $five_prime_seq . "\[$ref_allele\/$var_allele\]" .$three_prime_seq;
	
	    return { "Flanking_sequence" => $flanking_sequence, };
}

1;