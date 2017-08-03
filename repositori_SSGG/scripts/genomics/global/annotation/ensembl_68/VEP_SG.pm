=head1 LICENSE

 Copyright (c) 1999-2012 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

 Please email comments or questions to jm.rosa@sistemasgenomicos.com

=cut

# EnsEMBL module for Bio::EnsEMBL::Variation::Utils::Sequence
#
#

=head1 NAME

VEP_SG - Methods used by the Variant Effect Predictor

=head1 SYNOPSIS

  VEP_SG qw(configure);

  my $config = configure();

=head1 METHODS

=cut


use strict;
use warnings;

package VEP_SG;

# module list
use Getopt::Long;
use FileHandle;
use File::Path qw(make_path);
use Storable qw(nstore_fd fd_retrieve freeze thaw);
use Scalar::Util qw(weaken);
use Digest::MD5 qw(md5_hex);

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::VariationFeature;
use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor;
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(MAX_DISTANCE_FROM_TRANSCRIPT overlap);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(unambiguity_code);
use Bio::EnsEMBL::Variation::Utils::EnsEMBL2GFF3;
use Bio::EnsEMBL::Variation::StructuralVariationFeature;
use Bio::EnsEMBL::Variation::DBSQL::StructuralVariationFeatureAdaptor;
use Bio::EnsEMBL::Variation::TranscriptStructuralVariation;

# we need to manually include all these modules for caching to work
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::ProteinFeature;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DBSQL::GeneAdaptor;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::DBSQL::TranslationAdaptor;
use Bio::EnsEMBL::DBSQL::TranscriptAdaptor;
use Bio::EnsEMBL::DBSQL::MetaContainer;
use Bio::EnsEMBL::DBSQL::CoordSystemAdaptor;

use Exporter;
use vars qw(@ISA @EXPORT_OK);
@ISA = qw(Exporter);

# open socket pairs for cross-process comms
# use Socket;
# socketpair(CHILD, PARENT, AF_UNIX, SOCK_STREAM, PF_UNSPEC) or die "ERROR: Failed to open socketpair: $!";
# CHILD->autoflush(1);
# PARENT->autoflush(1);

@EXPORT_OK = qw(
    &parse_line
    &vf_to_consequences
    &validate_vf
    &read_cache_info
    &dump_adaptor_cache
    &load_dumped_adaptor_cache
    &load_dumped_variation_cache
    &get_all_consequences
    &get_slice
    &build_slice_cache
    &build_full_cache
    &regions_from_hash
    &get_time
    &debug
    &convert_to_vcf
    &progress
    &end_progress
    @REG_FEAT_TYPES
    @OUTPUT_COLS
    @VEP_WEB_CONFIG
    %FILTER_SHORTCUTS
);

our @OUTPUT_COLS = qw(
    Uploaded_variation
    Location
    Allele
    Gene
    Feature
    Feature_type
    Consequence
    cDNA_position
    CDS_position
    Protein_position
    Amino_acids
    Codons
    Existing_variation
    Extra
);

our @REG_FEAT_TYPES = qw(
    RegulatoryFeature
    MotifFeature
);

our @VEP_WEB_CONFIG = qw(
    format
    check_existing
    coding_only
    core_type
    hgnc
    protein
    hgvs
    terms
    check_frequency
    freq_filter
    freq_gt_lt
    freq_freq
    freq_pop
    sift
    polyphen
    regulatory
);

our %FILTER_SHORTCUTS = (
    upstream => {
        '5KB_upstream_variant' => 1,
        '2KB_upstream_variant' => 1,
    },
    downstream => {
        '5KB_downstream_variant'  => 1,
        '2KB_downstream_variant'  => 1,
        '500B_downstream_variant' => 1,
    },
    utr => {
        '5_prime_UTR_variant' => 1,
        '3_prime_UTR_variant' => 1,
    },
    splice => {
        splice_donor_variant    => 1,
        splice_acceptor_variant => 1,
        splice_region_variant   => 1,
    },
    coding_change => {
        stop_lost            => 1,
        stop_gained          => 1,
        missense_variant     => 1,
        frameshift_variant   => 1,
        inframe_insertion    => 1,
        inframe_deletion     => 1,
    },
    regulatory => {
        regulatory_region_variant => 1,
        TF_binding_site_variant   => 1,
    },
);

# parses a line of input, returns VF object(s)
sub parse_line {
    my $config = shift;
    my $line   = shift;
    
	my @data = split (/\s+/, $line);
	
	my ($chr,$start,$end, $ref, $var);
	$chr = $data[0];
	$chr =~ s/chr//;
	$start =  $end = $data[1];
	$ref = $data[3];
	$var = $data[4];
	
	#is deletion
	if (length ($ref) > length ($var)) {
		$start++;
		$end += (length ($ref) - 1);
		$ref =~ s/^$var//;
		$var = '-';
	}
	
	#is insertion
	if (length ($ref) < length ($var)) {
		$start++;
		$var =~ s/^$ref//;
		$ref = '-';
	}
	
	my $vf = Bio::EnsEMBL::Variation::VariationFeature->new_fast({
		start          => $start,
		end            => $end,
		allele_string  => $ref.'/'.$var,
		strand         => 1,
		map_weight     => 1,
		adaptor        => $config->{vfa},
		chr            => $chr,
	});
    
    return $vf;
}

# takes a variation feature and returns ready to print consequence information
sub vf_to_consequences {
    my $config = shift;
    my $vf = shift;
    
	my @return;
	push (@return, \&check_existing ($config, $vf)) if $config->{check_existing};
	
	return @return;

}


sub check_existing {
	my ($config, $vf) = @_;
	my $vfs = $config->{vfa}->fetch_all_by_Slice($vf->{slice});
	my @Vfs;
	my $n = 0;
	foreach my $var (@{$vfs}) {
		if (&checking_variant ($vf, $var) == 1 && !$vf->{existing}) {
			#print "$var->{variation_name}\n";
			$vf->{existing} = $var->{variation_name};
			$vf->{gmaf} = 1;
		}
		elsif (&checking_variant ($vf, $var) == 1 && $vf->{existing}) {
			$vf->{existing} .= ",$var->{variation_name}";
			$vf->{gmaf} .= ",1";
		}
		else {
			next;
		}
	}
	return $vf;
}

#Checking position and alleles  	print "Variation: ", $var->variation_name, " with alleles ", $var->allele_string, 
	       # " in chromosome ", $vf->{slice}->seq_region_name, " and position ", $var->start,"-",$var->end,"\n";	

sub checking_variant {
	my ($vf, $var) = @_;
	my $i = $vf->{start} - ($vf->{slice}->{start} + $var->{start} - 1);
	if ($i == 0) {
		return 1 if &check_alleles($vf, $var) == 1;
		return 0 if &check_alleles($vf, $var) != 1;
	}
	else {
		return 0;
	}
}

sub check_alleles {
	my ($vf, $var) = @_;
	print "$vf->{allele_string} -> $var->{allele_string}\n";
	return 1 if $vf->{allele_string} eq $var->{allele_string};
	return 0 if $vf->{allele_string} ne $var->{allele_string};
}

# gets a slice from the slice adaptor
sub get_slice {
    my $config = shift;
    my $chr = shift;
    my $start = shift;
    my $end = shift;
   
    my $slice;
    
    # first try to get a chromosome
    eval { $slice = $config->{sa}->fetch_by_region('chromosome',$chr,$start,$end,1); };
    
    # if failed, try to get any seq region
    if(!defined($slice)) {
        $slice = $config->{sa}->fetch_by_region(undef, $chr,$start,$end,1);
    }
    
    return $slice;
}

# DEBUG AND STATUS METHODS
##########################

# gets time
sub get_time() {
    my @time = localtime(time());

    # increment the month (Jan = 0)
    $time[4]++;

    # add leading zeroes as required
    for my $i(0..4) {
        $time[$i] = "0".$time[$i] if $time[$i] < 10;
    }

    # put the components together in a string
    my $time =
        ($time[5] + 1900)."-".
        $time[4]."-".
        $time[3]." ".
        $time[2].":".
        $time[1].":".
        $time[0];

    return $time;
}

# prints debug output with time
sub debug {
    my $text = (@_ ? (join "", @_) : "No message");
    my $time = get_time;
    
    print $time." - ".$text.($text =~ /\n$/ ? "" : "\n");
}

1;
