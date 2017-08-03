=head1 LICENSE                                                                                                                                                                           
                                                                                                                
 This software is distributed under a modified Apache license.                                                                                                                      
                                                                                                                     
=head1 CONTACT                                                                                                       

 Guillermo Marco Puche <guillermo.marco@sistemasgenomicos.com>
 
=cut

=head1 NAME

 Intepro

=head1 SYNOPSIS

perl variant_effect_predictor.pl --plugin interpro

=head1 DESCRIPTION
 
 Obtains the Interpro AC & Description. 

=cut

package Interpro;

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

#sub feature_types {
#    return ['Transcript'];
#}

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
				
			my $aa_coordinate=&format_coords($translation->start, $translation->end);		
			
			if ($aa_coordinate ne '-') {
				my $aa_start = $aa_coordinate;
				my $aa_end = $aa_coordinate;
				
					if ($aa_coordinate =~ /(\d+)-(\d+)/) {
						my @pos = split (/-/, $aa_coordinate);
						$aa_start = $pos[0];
						$aa_end = $pos[1];	
						}
				foreach my $pf(@{$translation->get_all_ProteinFeatures}) {
				#while ( my $pf = shift @{$translation->get_all_ProteinFeatures} ) {
					
					if(($aa_start >=$pf->start && $aa_start<=$pf->end) || ($aa_end >=$pf->start && $aa_end<=$pf->end)){
						if ($pf->interpro_ac() && $pf->idesc()) {	
							$interpro_ac = $pf->interpro_ac;
							$idesc = $pf->idesc;
							$idesc =~s/[\$#@~!\+\-&*()\[\]]+//g;

							#Added chomp fix
							#Quantifier follows nothing in regex; marked by <-- HERE in m/(+ <-- HERE )RNA_virus_helicase_core_dom/
							#\Q$replace\E
							#chomp($interpro_ac);
							#chomp($idesc);
							
							if (!$interpro_data[0] && !$interpro_data[1]) {
							#if ($interpro_data[0] eq "" && $interpro_data[1] eq "") {
										$interpro_data[0] = $interpro_ac;
										$interpro_data[1] = $idesc;
							}
							#elsif ($interpro_data[0] eq '$interpro_ac' || $interpro_data[1] eq '$idesc') {
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
		}
	}	
	return { %interpro, };
}

1;
