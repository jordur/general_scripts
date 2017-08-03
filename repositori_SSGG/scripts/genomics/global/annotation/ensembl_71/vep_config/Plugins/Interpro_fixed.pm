=head1 LICENSE                                                                                                                                                                                                                                                                                          
 This software is distributed under a modified Apache license.                                                                                                                      
                                                                                                                     
=head1 CONTACT                                                                                                       

 Guillermo Marco Puche <guillermo.marco@sistemasgenomicos.com>
 
=cut

=head1 NAME

 Intepro

=head1 SYNOPSIS

perl variant_effect_predictor.pl --plugin Interpro_fixed

=head1 DESCRIPTION
 
 Obtains the Interpro AC & Description. Fixed 28/06/2013. Now reports are correct ! 

=cut

package Interpro_fixed;

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
    return ['Transcript'];
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
	my $self = shift;
	my $tva = shift;

	my %interpro;
	my $Transcript="";
   	my $var_feature_vep = $tva->variation_feature;
   	my $tr = $tva->transcript;
   	
   	if ($tr){
   		$Transcript = $tr->stable_id;	
   	}
    	
   	my $config = $self->{config};
    	
	my ($interpro_ac, $idesc);
	my @interpro_data;
		
	$interpro_data[0]="";
	$interpro_data[1]="";
		
	my $new_vf = $tva->variation_feature;

	foreach my $tv(@{$new_vf->get_all_TranscriptVariations}) {
		my $transcript_id=$tv->transcript->stable_id;
		#print "\n$Transcript \t $transcript_id\n";
		#print "$transcript_id - $Transcript\n";
		if ($transcript_id && $transcript_id eq $Transcript) {
			my $slice_adaptor = $config->{reg}->get_adaptor(($config->{species}), 'core', 'Transcript');
			my $tr = $slice_adaptor->fetch_by_stable_id($Transcript);
			my $translation = $tr->translation;

			if($translation) {
				my $pfeatures = $translation->get_all_ProteinFeatures();
				my $aa_coordinate=&format_coords($tv->translation_start, $tv->translation_end);
				
				#print "AA coordinate: $aa_coordinate\n";
				
				if ($aa_coordinate ne '-') {
					my $aa_start = $aa_coordinate;
					my $aa_end = $aa_coordinate;
					
					if ($aa_coordinate =~ /(\d+)-(\d+)/) {
						my @pos = split (/-/, $aa_coordinate);
						$aa_start = $pos[0];
						$aa_end = $pos[1];	
					}
					while ( my $pfeature = shift @{$pfeatures} ) {
						my $logic_name = $pfeature->analysis()->logic_name();
						#print "if(($aa_start >=$pfeature->start() && $aa_start<=$pfeature->end()) || ($aa_end >=$pfeature->start() && $aa_end<=$pfeature->end()))\n";
						if(($aa_start >=$pfeature->start() && $aa_start<=$pfeature->end()) || ($aa_end >=$pfeature->start() && $aa_end<=$pfeature->end())){					
							if ($pfeature->interpro_ac() && $pfeature->idesc()) {
								my $ipro_id = $pfeature->interpro_ac();
								my $ipro_desc = $pfeature->idesc();
								
								#if ($ipro_desc){
									#$ipro_desc =~ s/[\$#@~!\+\-&*()\[\]]+//g;
								#}
								
								if (!$interpro_data[0] && !$interpro_data[1]) {
									$interpro_data[0] = $ipro_id;
									$interpro_data[1] = $ipro_desc;
									
								}
								elsif ($interpro_data[0] =~ /$ipro_id/ || $interpro_data[1] =~ /\Q$ipro_desc\E/) {
									next;
								}
								else {
									$interpro_data[0] .= ','.$ipro_id;
									$interpro_data[1] .= ','.$ipro_desc;
								}
								
								$interpro{"InterPro_IDs"} = $interpro_data[0];
								$interpro{"InterPro_descriptions"} = $interpro_data[1];
							}
						}
					}
				}
			}
		}
	}
		
	
	return { %interpro, };
}

1;
