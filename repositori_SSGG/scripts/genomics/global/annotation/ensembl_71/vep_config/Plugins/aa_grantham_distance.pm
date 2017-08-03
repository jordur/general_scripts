=head1 LICENSE                                                                                                                                                                          
                                                                                                                     
 This software is distributed under a modified Apache license.                                                                                                                      
                                                                                                                     
=head1 CONTACT                                                                                                       

 Guillermo Marco Puche <guillermo.marco@sistemasgenomicos.com>
 JM.Rosa aka Coyote por sus dos funciones get_aa_distance() y get_distance()
 
=cut

=head1 NAME

 aa_grantham_distance

=head1 SYNOPSIS

perl variant_effect_predictor.pl -i variations.vcf --plugin aa_grantham_distance

=head1 DESCRIPTION
 
 Obtains the AA Grantham Distance for each HGVSp

=cut

package aa_grantham_distance;

use strict;
use warnings;

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
    return ['Feature','Intergenic'];
}

sub variant_feature_types {
    return ['VariationFeature'];
}

sub get_header_info {
    return {
        AA_GRANTHAM_DISTANCE => "Column containing the AA GRANTHAM DISTANCE obtained from HGVSc.",
    };
}

sub run {
	my ($self, $tva) = @_;
	my $aa_grantham_distance = '-';
	
	my $vf = $tva->variation_feature;

	#Intentamos obtener el valor del HGSVp
	my $hgvsp = $tva->hgvs_protein;
	
	#Si existe HGVSp parseamos el formato para obtener la distancia.
	if(defined $hgvsp and $hgvsp !~ '(p.=)'){

		$aa_grantham_distance = &get_aa_distance ($hgvsp);
	}
    
	    return { "AA_GRANTHAM_DISTANCE" => $aa_grantham_distance, };
}

sub get_aa_distance () {
	my ($hgvsp) = @_;	
	my @aminoacids = split (/:p\./, $hgvsp);
	my @aa = split (/\d+/, $aminoacids[1]);
	my $amino_cond="Gly|Ala|Val|Leu|Ile|Met|Phe|Trp|Pro|Ser|Thr|Cys|Tyr|Asn|Gln|Asp|Glu|Lys|Arg|His";
	my $dist = '-';
	
	if(($aa[0] =~ /^($amino_cond)$/i) && ($aa[1] =~ m/^($amino_cond)$/i)){
		$dist = &get_distance ($aa[0], $aa[1]);
	}
	return $dist;
}

sub get_distance () {
	my ($aa1, $aa2, $hgvsp) = @_;
	my %Distances;
	$Distances{Arg}{Ser}=110;
	$Distances{Leu}{Ser}=145;
	$Distances{Leu}{Arg}=102;
	$Distances{Pro}{Ser}=74;
	$Distances{Pro}{Arg}=103;
	$Distances{Pro}{Leu}=98;
	$Distances{Thr}{Ser}=58;
	$Distances{Thr}{Arg}=71;
	$Distances{Thr}{Leu}=92;
	$Distances{Thr}{Pro}=38;
	$Distances{Ala}{Ser}=99;
	$Distances{Ala}{Arg}=112;
	$Distances{Ala}{Leu}=96;
	$Distances{Ala}{Pro}=27;
	$Distances{Ala}{Thr}=58;
	$Distances{Val}{Ser}=124;
	$Distances{Val}{Arg}=96;
	$Distances{Val}{Leu}=32;
	$Distances{Val}{Pro}=68;
	$Distances{Val}{Thr}=69;
	$Distances{Val}{Ala}=64;
	$Distances{Gly}{Ser}=56;
	$Distances{Gly}{Arg}=125;
	$Distances{Gly}{Leu}=138;
	$Distances{Gly}{Pro}=42;
	$Distances{Gly}{Thr}=59;
	$Distances{Gly}{Ala}=60;
	$Distances{Gly}{Val}=109;
	$Distances{Ile}{Ser}=142;
	$Distances{Ile}{Arg}=97;
	$Distances{Ile}{Leu}=5;
	$Distances{Ile}{Pro}=95;
	$Distances{Ile}{Thr}=89;
	$Distances{Ile}{Ala}=94;
	$Distances{Ile}{Val}=29;
	$Distances{Ile}{Gly}=135;
	$Distances{Phe}{Ser}=155;
	$Distances{Phe}{Arg}=97;
	$Distances{Phe}{Leu}=22;
	$Distances{Phe}{Pro}=114;
	$Distances{Phe}{Thr}=103;
	$Distances{Phe}{Ala}=113;
	$Distances{Phe}{Val}=50;
	$Distances{Phe}{Gly}=153;
	$Distances{Phe}{Ile}=21;
	$Distances{Tyr}{Ser}=144;
	$Distances{Tyr}{Arg}=77;
	$Distances{Tyr}{Leu}=36;
	$Distances{Tyr}{Pro}=110;
	$Distances{Tyr}{Thr}=92;
	$Distances{Tyr}{Ala}=112;
	$Distances{Tyr}{Val}=55;
	$Distances{Tyr}{Gly}=147;
	$Distances{Tyr}{Ile}=33;
	$Distances{Tyr}{Phe}=22;
	$Distances{Cys}{Ser}=112;
	$Distances{Cys}{Arg}=180;
	$Distances{Cys}{Leu}=198;
	$Distances{Cys}{Pro}=169;
	$Distances{Cys}{Thr}=149;
	$Distances{Cys}{Ala}=195;
	$Distances{Cys}{Val}=192;
	$Distances{Cys}{Gly}=159;
	$Distances{Cys}{Ile}=198;
	$Distances{Cys}{Phe}=205;
	$Distances{Cys}{Tyr}=194;
	$Distances{His}{Ser}=89;
	$Distances{His}{Arg}=29;
	$Distances{His}{Leu}=99;
	$Distances{His}{Pro}=77;
	$Distances{His}{Thr}=47;
	$Distances{His}{Ala}=86;
	$Distances{His}{Val}=84;
	$Distances{His}{Gly}=98;
	$Distances{His}{Ile}=94;
	$Distances{His}{Phe}=100;
	$Distances{His}{Tyr}=83;
	$Distances{His}{Cys}=174;
	$Distances{Gln}{Ser}=68;
	$Distances{Gln}{Arg}=43;
	$Distances{Gln}{Leu}=113;
	$Distances{Gln}{Pro}=76;
	$Distances{Gln}{Thr}=42;
	$Distances{Gln}{Ala}=91;
	$Distances{Gln}{Val}=96;
	$Distances{Gln}{Gly}=87;
	$Distances{Gln}{Ile}=109;
	$Distances{Gln}{Phe}=116;
	$Distances{Gln}{Tyr}=99;
	$Distances{Gln}{Cys}=154;
	$Distances{Gln}{His}=24;
	$Distances{Asn}{Ser}=46;
	$Distances{Asn}{Arg}=86;
	$Distances{Asn}{Leu}=153;
	$Distances{Asn}{Pro}=91;
	$Distances{Asn}{Thr}=65;
	$Distances{Asn}{Ala}=111;
	$Distances{Asn}{Val}=133;
	$Distances{Asn}{Gly}=80;
	$Distances{Asn}{Ile}=149;
	$Distances{Asn}{Phe}=158;
	$Distances{Asn}{Tyr}=143;
	$Distances{Asn}{Cys}=139;
	$Distances{Asn}{His}=68;
	$Distances{Asn}{Gln}=46;
	$Distances{Lys}{Ser}=121;
	$Distances{Lys}{Arg}=26;
	$Distances{Lys}{Leu}=107;
	$Distances{Lys}{Pro}=103;
	$Distances{Lys}{Thr}=78;
	$Distances{Lys}{Ala}=106;
	$Distances{Lys}{Val}=97;
	$Distances{Lys}{Gly}=127;
	$Distances{Lys}{Ile}=102;
	$Distances{Lys}{Phe}=102;
	$Distances{Lys}{Tyr}=85;
	$Distances{Lys}{Cys}=202;
	$Distances{Lys}{His}=32;
	$Distances{Lys}{Gln}=53;
	$Distances{Lys}{Asn}=94;
	$Distances{Asp}{Ser}=65;
	$Distances{Asp}{Arg}=96;
	$Distances{Asp}{Leu}=172;
	$Distances{Asp}{Pro}=108;
	$Distances{Asp}{Thr}=85;
	$Distances{Asp}{Ala}=126;
	$Distances{Asp}{Val}=152;
	$Distances{Asp}{Gly}=94;
	$Distances{Asp}{Ile}=168;
	$Distances{Asp}{Phe}=177;
	$Distances{Asp}{Tyr}=160;
	$Distances{Asp}{Cys}=154;
	$Distances{Asp}{His}=81;
	$Distances{Asp}{Gln}=61;
	$Distances{Asp}{Asn}=23;
	$Distances{Asp}{Lys}=101;
	$Distances{Glu}{Ser}=80;
	$Distances{Glu}{Arg}=54;
	$Distances{Glu}{Leu}=138;
	$Distances{Glu}{Pro}=93;
	$Distances{Glu}{Thr}=65;
	$Distances{Glu}{Ala}=107;
	$Distances{Glu}{Val}=121;
	$Distances{Glu}{Gly}=98;
	$Distances{Glu}{Ile}=134;
	$Distances{Glu}{Phe}=140;
	$Distances{Glu}{Tyr}=122;
	$Distances{Glu}{Cys}=170;
	$Distances{Glu}{His}=40;
	$Distances{Glu}{Gln}=29;
	$Distances{Glu}{Asn}=42;
	$Distances{Glu}{Lys}=56;
	$Distances{Glu}{Asp}=45;
	$Distances{Met}{Ser}=135;
	$Distances{Met}{Arg}=91;
	$Distances{Met}{Leu}=15;
	$Distances{Met}{Pro}=87;
	$Distances{Met}{Thr}=81;
	$Distances{Met}{Ala}=84;
	$Distances{Met}{Val}=21;
	$Distances{Met}{Gly}=127;
	$Distances{Met}{Ile}=10;
	$Distances{Met}{Phe}=28;
	$Distances{Met}{Tyr}=36;
	$Distances{Met}{Cys}=196;
	$Distances{Met}{His}=87;
	$Distances{Met}{Gln}=101;
	$Distances{Met}{Asn}=142;
	$Distances{Met}{Lys}=95;
	$Distances{Met}{Asp}=160;
	$Distances{Met}{Glu}=126;
	$Distances{Trp}{Ser}=177;
	$Distances{Trp}{Arg}=101;
	$Distances{Trp}{Leu}=61;
	$Distances{Trp}{Pro}=147;
	$Distances{Trp}{Thr}=128;
	$Distances{Trp}{Ala}=148;
	$Distances{Trp}{Val}=88;
	$Distances{Trp}{Gly}=184;
	$Distances{Trp}{Ile}=61;
	$Distances{Trp}{Phe}=40;
	$Distances{Trp}{Tyr}=37;
	$Distances{Trp}{Cys}=215;
	$Distances{Trp}{His}=115;
	$Distances{Trp}{Gln}=130;
	$Distances{Trp}{Asn}=174;
	$Distances{Trp}{Lys}=110;
	$Distances{Trp}{Asp}=181;
	$Distances{Trp}{Glu}=152;
	$Distances{Trp}{Met}=67;
	$Distances{Ser}{Arg} = 110;
	$Distances{Ser}{Leu} = 145;
	$Distances{Arg}{Leu} = 102;
	$Distances{Ser}{Pro} = 74;
	$Distances{Arg}{Pro} = 103;
	$Distances{Leu}{Pro} = 98;
	$Distances{Ser}{Thr} = 58;
	$Distances{Arg}{Thr} = 71;
	$Distances{Leu}{Thr} = 92;
	$Distances{Pro}{Thr} = 38;
	$Distances{Ser}{Ala} = 99;
	$Distances{Arg}{Ala} = 112;
	$Distances{Leu}{Ala} = 96;
	$Distances{Pro}{Ala} = 27;
	$Distances{Thr}{Ala} = 58;
	$Distances{Ser}{Val} = 124;
	$Distances{Arg}{Val} = 96;
	$Distances{Leu}{Val} = 32;
	$Distances{Pro}{Val} = 68;
	$Distances{Thr}{Val} = 69;
	$Distances{Ala}{Val} = 64;
	$Distances{Ser}{Gly} = 56;
	$Distances{Arg}{Gly} = 125;
	$Distances{Leu}{Gly} = 138;
	$Distances{Pro}{Gly} = 42;
	$Distances{Thr}{Gly} = 59;
	$Distances{Ala}{Gly} = 60;
	$Distances{Val}{Gly} = 109;
	$Distances{Ser}{Ile} = 142;
	$Distances{Arg}{Ile} = 97;
	$Distances{Leu}{Ile} = 5;
	$Distances{Pro}{Ile} = 95;
	$Distances{Thr}{Ile} = 89;
	$Distances{Ala}{Ile} = 94;
	$Distances{Val}{Ile} = 29;
	$Distances{Gly}{Ile} = 135;
	$Distances{Ser}{Phe} = 155;
	$Distances{Arg}{Phe} = 97;
	$Distances{Leu}{Phe} = 22;
	$Distances{Pro}{Phe} = 114;
	$Distances{Thr}{Phe} = 103;
	$Distances{Ala}{Phe} = 113;
	$Distances{Val}{Phe} = 50;
	$Distances{Gly}{Phe} = 153;
	$Distances{Ile}{Phe} = 21;
	$Distances{Ser}{Tyr} = 144;
	$Distances{Arg}{Tyr} = 77;
	$Distances{Leu}{Tyr} = 36;
	$Distances{Pro}{Tyr} = 110;
	$Distances{Thr}{Tyr} = 92;
	$Distances{Ala}{Tyr} = 112;
	$Distances{Val}{Tyr} = 55;
	$Distances{Gly}{Tyr} = 147;
	$Distances{Ile}{Tyr} = 33;
	$Distances{Phe}{Tyr} = 22;
	$Distances{Ser}{Cys} = 112;
	$Distances{Arg}{Cys} = 180;
	$Distances{Leu}{Cys} = 198;
	$Distances{Pro}{Cys} = 169;
	$Distances{Thr}{Cys} = 149;
	$Distances{Ala}{Cys} = 195;
	$Distances{Val}{Cys} = 192;
	$Distances{Gly}{Cys} = 159;
	$Distances{Ile}{Cys} = 198;
	$Distances{Phe}{Cys} = 205;
	$Distances{Tyr}{Cys} = 194;
	$Distances{Ser}{His} = 89;
	$Distances{Arg}{His} = 29;
	$Distances{Leu}{His} = 99;
	$Distances{Pro}{His} = 77;
	$Distances{Thr}{His} = 47;
	$Distances{Ala}{His} = 86;
	$Distances{Val}{His} = 84;
	$Distances{Gly}{His} = 98;
	$Distances{Ile}{His} = 94;
	$Distances{Phe}{His} = 100;
	$Distances{Tyr}{His} = 83;
	$Distances{Cys}{His} = 174;
	$Distances{Ser}{Gln} = 68;
	$Distances{Arg}{Gln} = 43;
	$Distances{Leu}{Gln} = 113;
	$Distances{Pro}{Gln} = 76;
	$Distances{Thr}{Gln} = 42;
	$Distances{Ala}{Gln} = 91;
	$Distances{Val}{Gln} = 96;
	$Distances{Gly}{Gln} = 87;
	$Distances{Ile}{Gln} = 109;
	$Distances{Phe}{Gln} = 116;
	$Distances{Tyr}{Gln} = 99;
	$Distances{Cys}{Gln} = 154;
	$Distances{His}{Gln} = 24;
	$Distances{Ser}{Asn} = 46;
	$Distances{Arg}{Asn} = 86;
	$Distances{Leu}{Asn} = 153;
	$Distances{Pro}{Asn} = 91;
	$Distances{Thr}{Asn} = 65;
	$Distances{Ala}{Asn} = 111;
	$Distances{Val}{Asn} = 133;
	$Distances{Gly}{Asn} = 80;
	$Distances{Ile}{Asn} = 149;
	$Distances{Phe}{Asn} = 158;
	$Distances{Tyr}{Asn} = 143;
	$Distances{Cys}{Asn} = 139;
	$Distances{His}{Asn} = 68;
	$Distances{Gln}{Asn} = 46;
	$Distances{Ser}{Lys} = 121;
	$Distances{Arg}{Lys} = 26;
	$Distances{Leu}{Lys} = 107;
	$Distances{Pro}{Lys} = 103;
	$Distances{Thr}{Lys} = 78;
	$Distances{Ala}{Lys} = 106;
	$Distances{Val}{Lys} = 97;
	$Distances{Gly}{Lys} = 127;
	$Distances{Ile}{Lys} = 102;
	$Distances{Phe}{Lys} = 102;
	$Distances{Tyr}{Lys} = 85;
	$Distances{Cys}{Lys} = 202;
	$Distances{His}{Lys} = 32;
	$Distances{Gln}{Lys} = 53;
	$Distances{Asn}{Lys} = 94;
	$Distances{Ser}{Asp} = 65;
	$Distances{Arg}{Asp} = 96;
	$Distances{Leu}{Asp} = 172;
	$Distances{Pro}{Asp} = 108;
	$Distances{Thr}{Asp} = 85;
	$Distances{Ala}{Asp} = 126;
	$Distances{Val}{Asp} = 152;
	$Distances{Gly}{Asp} = 94;
	$Distances{Ile}{Asp} = 168;
	$Distances{Phe}{Asp} = 177;
	$Distances{Tyr}{Asp} = 160;
	$Distances{Cys}{Asp} = 154;
	$Distances{His}{Asp} = 81;
	$Distances{Gln}{Asp} = 61;
	$Distances{Asn}{Asp} = 23;
	$Distances{Lys}{Asp} = 101;
	$Distances{Ser}{Glu} = 80;
	$Distances{Arg}{Glu} = 54;
	$Distances{Leu}{Glu} = 138;
	$Distances{Pro}{Glu} = 93;
	$Distances{Thr}{Glu} = 65;
	$Distances{Ala}{Glu} = 107;
	$Distances{Val}{Glu} = 121;
	$Distances{Gly}{Glu} = 98;
	$Distances{Ile}{Glu} = 134;
	$Distances{Phe}{Glu} = 140;
	$Distances{Tyr}{Glu} = 122;
	$Distances{Cys}{Glu} = 170;
	$Distances{His}{Glu} = 40;
	$Distances{Gln}{Glu} = 29;
	$Distances{Asn}{Glu} = 42;
	$Distances{Lys}{Glu} = 56;
	$Distances{Asp}{Glu} = 45;
	$Distances{Ser}{Met} = 135;
	$Distances{Arg}{Met} = 91;
	$Distances{Leu}{Met} = 15;
	$Distances{Pro}{Met} = 87;
	$Distances{Thr}{Met} = 81;
	$Distances{Ala}{Met} = 84;
	$Distances{Val}{Met} = 21;
	$Distances{Gly}{Met} = 127;
	$Distances{Ile}{Met} = 10;
	$Distances{Phe}{Met} = 28;
	$Distances{Tyr}{Met} = 36;
	$Distances{Cys}{Met} = 196;
	$Distances{His}{Met} = 87;
	$Distances{Gln}{Met} = 101;
	$Distances{Asn}{Met} = 142;
	$Distances{Lys}{Met} = 95;
	$Distances{Asp}{Met} = 160;
	$Distances{Glu}{Met} = 126;
	$Distances{Ser}{Trp} = 177;
	$Distances{Arg}{Trp} = 101;
	$Distances{Leu}{Trp} = 61;
	$Distances{Pro}{Trp} = 147;
	$Distances{Thr}{Trp} = 128;
	$Distances{Ala}{Trp} = 148;
	$Distances{Val}{Trp} = 88;
	$Distances{Gly}{Trp} = 184;
	$Distances{Ile}{Trp} = 61;
	$Distances{Phe}{Trp} = 40;
	$Distances{Tyr}{Trp} = 37;
	$Distances{Cys}{Trp} = 215;
	$Distances{His}{Trp} = 115;
	$Distances{Gln}{Trp} = 130;
	$Distances{Asn}{Trp} = 174;
	$Distances{Lys}{Trp} = 110;
	$Distances{Asp}{Trp} = 181;
	$Distances{Glu}{Trp} = 152;
	$Distances{Met}{Trp} = 67;
	
	my $dist;
	my $distance;
	
	$dist = $Distances{$aa1}{$aa2};	
	if ($dist < 70) {
		$distance = 'small_('.$dist.')';
	}
	elsif ($dist > 140) {
		$distance = 'large_('.$dist.')';
	}
	else {
		$distance = 'moderate_('.$dist.')';
	}
	
	return $distance;
}

1;
