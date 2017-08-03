=head1 LICENSE                                                                                                                                                                           
                                                                                                                
 This software is distributed under a modified Apache license.                                                                                                                      
                                                                                                                     
=head1 CONTACT                                                                                                       

 Guillermo Marco Puche <guillermo.marco@sistemasgenomicos.com>
 
=cut

=head1 NAME

 Intepro

=head1 SYNOPSIS

perl variant_effect_predictor.pl --plugin Interpro2

=head1 DESCRIPTION
 
 Obtains the Interpro AC & Description. 

=cut

package Interpro2;

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
	return ['Transcript']
}

sub get_header_info {
    return {
        InterPro_IDs => "Interpro AC.",
        InterPro_descriptions => "Interpro Description.",
    };
}

sub format_coords {
	my ($start, $end) = @_;
	
	if(!defined($start)) {
		return '-';
	}
	elsif(!defined($end)) {
		return $start;
	}
	elsif($start == $end) {
		return $start;
	}
	elsif($start > $end) {
		return $end.'-'.$start;
	}
	else {
		return $start.'-'.$end;
	}
}

sub run {
		my %interpro;
		my ($self, $tva) = @_;
		my ($interpro_ac, $idesc);
		my @interpro_data;
		$interpro_data[0]="";
		$interpro_data[1]="";
	
		my $tr = $tva->transcript;
		my $translation = $tr->translation;
		
		if ($translation){
			my $pfeatures = $translation->get_all_ProteinFeatures();
			while ( my $pf = shift @{$pfeatures} ) {
				
				if ($pf->interpro_ac() && $pf->idesc()) {	
					$interpro_ac = $pf->interpro_ac;
					$idesc = $pf->idesc;
					$idesc =~s/[\$#@~!\+\-&*()\[\]]+//g;
							
					if (!$interpro_data[0] && !$interpro_data[1]) {
							$interpro_data[0] = $interpro_ac;
							$interpro_data[1] = $idesc;
					}
					
					elsif ($interpro_data[0] =~ /$interpro_ac/ || $interpro_data[1] =~ /$idesc/) {
						next;
					}
					
					else {
						$interpro_data[0] .= ','.$interpro_ac;
						$interpro_data[1] .= ','.$idesc;
					}
							
					$interpro{"InterPro_IDs"} = $interpro_data[0];
					$interpro{"InterPro_descriptions"} = $interpro_data[1];
				}
			}
		}
	return { %interpro, };
}

1;
