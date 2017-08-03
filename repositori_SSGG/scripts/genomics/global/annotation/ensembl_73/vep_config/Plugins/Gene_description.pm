=head1 LICENSE
                                                                                                                     
 Copyright (c) 2013 Sistemas Genómicos S.L - All rights reserved.                                                                      
                                                                                                        
 This software is distributed under a modified Apache license.                                                                                                                      
                                                                                                                     
=head1 CONTACT                                                                                                       

 Guillermo Marco Puche <guillermo.marco@sistemasgenomicos.com>
    
=cut

=head1 NAME

 vcf_input

=head1 SYNOPSIS



=head1 DESCRIPTION
 
 If gene exists plugin adds Gene_description column to output.

=cut

package Gene_description;

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
	
    return $self;
}

sub version {
    return '73';
}

sub feature_types {
    return ['Transcript'];
}

sub get_header_info {
    return {
    	Gene_description => "Gene description.",
    };
}

sub run {
	
	my $self = shift;
	my $tva = shift;
	
	## Def vars
	my %res;
	my $gene;
	my $gene_desc = "-";
	
	## Def Handlers
    my $config = $self->{config};
    
    my $tv = $tva->transcript_variation;
    my $tr = $tv->transcript;
    my $transcript_id = $tr->stable_id();
    $gene = $config->{ga}->fetch_by_transcript_stable_id($transcript_id);
    $gene_desc = $gene->description;
    
    if (defined $gene_desc){
#    	my $gene_id = $gene->external_name;
    	$gene_desc = $gene->description;
    	
    	$gene_desc =~ s/; /\&/g;
		$gene_desc =~ s/;/\&/g;
    	
    	$gene_desc =~ s/\s+/_/g;
		$gene_desc =~ s/\[.*\]//g;
		$gene_desc =~ s/_$//;
    }
    
    $res{"Gene_description"} = $gene_desc;
    
    return { %res };
}
1;