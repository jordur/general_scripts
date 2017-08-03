#!/usr/bin/perl -w

=head1 LICENSE

  This script was developed by Sistemas Genomicos SL

=head1 CONTACT

  Please email comments or questions to:
  jm.rosa@sistemasgenomicos.com

=cut

=head1 NAME

Customized Variant Effect Predictor - a script to predict the consequences of genomic variants

Version 1.0

by Juan Manuel Rosa (jm.rosa@sistemasgenomicos.com)
=cut

use strict;
use Getopt::Long;
use FileHandle;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor;
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(MAX_DISTANCE_FROM_TRANSCRIPT);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(unambiguity_code);
use Bio::EnsEMBL::Utils::Exception qw(throw);	
use Bio::EnsEMBL::Transcript();
use Bio::EnsEMBL::ProteinFeature;
use Bio::EnsEMBL::OntologyTerm;

# configure from command line opts
my $config = &configure(scalar @ARGV);

# run the main sub routine
&main($config);

	# this is the main sub-routine - it needs the configured $config hash
sub main {
	my $config = shift;
	
	debug("Starting...") if defined $config->{verbose};

	my %transcripts;
	
	if (defined $config->{bed_file}) {
		get_transcript_translation($config, \%transcripts);
	}
	
	my $species = $config->{species};
	my $input = $config->{input_file};
	my $output = $config->{tmp_file};
	my $host = $config->{host};
	my $port = $config->{port};
	my $user = $config->{user};
	my $password = $config->{password};
	my $pass = '';
	if ($password) {
		$pass = "--password $password";
	}
	
	
	if ($species eq 'human' || $species eq 'homo_sapiens') {
		my $command = "variant_effect_predictor_64.pl -i $input --format vcf -o $output --sift b --polyphen b --condel b --regulatory --hgnc  --hgvs --no_intergenic --force_overwrite --species $species --host $host --port $port --user $user $pass --check_existing";
		
		print "Launching command:\n$command\n".localtime()."\n";
	
		#`variant_effect_predictor_64.pl -i $input --format vcf -o $output --sift b --polyphen b --condel b --regulatory --hgnc  --hgvs 	--no_intergenic --force_overwrite --species $species --host $host --port $port --user $user $pass --check_existing > ensembl64_tmp.log`;
	}
	else {
		my $command = "variant_effect_predictor_64.pl -i $input --format vcf -o $output --regulatory --hgnc  --hgvs --no_intergenic --force_overwrite --species $species --host $host --port $port --user $user $pass --check_existing";
		
		print "Launching command:\n$command\n".localtime()."\n";
	
		#`variant_effect_predictor_64.pl -i $input --format vcf -o $output --regulatory --hgnc  --hgvs 	--no_intergenic --force_overwrite --species $species --host $host --port $port --user $user $pass --check_existing  > ensembl64_tmp.log`;
	}
	
	print "Completing information... ".localtime()."\n";
	
	my $line_number = 0;
	my $vf_count;
	my $tmp_file_handle = $config->{tmp_file_handle}; 
	my $last_number = 0;
	my %populations;

	my $out_file_handle = $config->{out_file_handle};
	
# read the file
	while (my $line = <$tmp_file_handle>) {
		chomp $line;
		$line_number++;

		# header line?
		next if ($line =~/^\#/);

		if (defined $config->{bed_file}) {
			my $index;
			checking_transcript(\%transcripts, $line, \$index);
			print "$index\n";
			if ($index == 0) {
				next;
			}
			else {
				print_output($config, $line, \%populations, \$last_number, $out_file_handle, \%transcripts);
			}
		}
		else {
			print_output($config, $line, \%populations, \$last_number, $out_file_handle, \%transcripts);
		}
		
		$vf_count++;
		debug("Processed $vf_count variants") if $vf_count =~ /0$/ && defined($config->{verbose});
	}
	print $out_file_handle "\n";
	
	foreach my $pop (keys %populations) {
		print $out_file_handle "(",$populations{$pop}{number},")\t$pop\n",
	}
	
	debug("Finished!") if defined $config->{verbose};
}

sub checking_transcript () {
	my ($transcripts, $line, $value) = @_;
	my $trans_id = (split (/\s+/, $line))[4];
	if (defined $$transcripts{$trans_id}{name}) {
		$$value = 1;
	}
	else {
		$$value = 0;
	}
	print "$$value\t";
}

sub print_output {
	my ($config, $line, $populations, $last_number, $out_file_handle, $transcripts) = @_;
	
	my ($Uploaded_variation, $Location, $Allele, $Gene, $Transcript, $Feature_type, $Consequence, $cDNA_position, $CDS_position,$Protein_position, $Amino_acids, $Codons, $Existing_variation, $Extra) = split (/\s+/, $line);
	
	my $trans_id;
	
	if (defined $config->{bed_file}) {
		if (defined $$transcripts{$Transcript}) {
			$trans_id = $$transcripts{$Transcript}{name};
		}
	}
		
	my @rs = split (/,/, $Existing_variation);
	my ($genotype, @maf, %HGNC, $condel, $polyphen, $sift, $hgvsc, $conservation, @coord, $end, @interpro, @sequence, $sequence);
	
	foreach my $snp (@rs) {
		$genotype = (split/_/, $Uploaded_variation)[2];
		$maf[3] = 1;
		if ($snp =~ m/rs/i) {
			get_rs_data($config, $snp, $genotype, \@maf, $populations, $last_number);
		}
		if ($maf[3] == 0) {
			next;
		}
		else {
			$Location =~ s/\d+://;
			my @extra = split (/;/, $Extra);
			
			foreach my $field (@extra) {
				if ($field =~ m/HGVSc/) {
					$hgvsc = $field;
					$hgvsc =~ s/HGVSc=//;
				}
				elsif ($field =~ m/PolyPhen/) {
					$polyphen = $field;
					$polyphen =~ s/PolyPhen=//;
				}
				elsif ($field =~ m/Condel/) {
					$condel = $field;
					$condel =~ s/Condel=//;
				} 
				elsif ($field =~ m/SIFT/) {
					$sift = $field;
					$sift =~ s/SIFT=//;
				} 
			}
			
			undef @extra;
			
			if ($Transcript =~ m/ENST/ && $Extra =~ m/HGNC/) {
				get_hgnc($config, $Transcript, \%HGNC);
			}
			@coord = split (/\_/, $Uploaded_variation);
			$end = $coord[1] + (length ($coord[2]) - 3);
			get_conservation_score ($coord[0], $coord[1], $end, $conservation);
			get_interprot ($config, $Transcript, \@interpro, $coord[0], $coord[1], $end, $genotype);
			get_sequence ($config, $coord[0], $coord[1], $genotype, \@sequence);
			$sequence  = $sequence[0]."[".$genotype."]".$sequence[1];
						
			print $out_file_handle $Uploaded_variation,"\t",
		#	print $Uploaded_variation,"\t",
				$Location,"\t",
				"Depth\t",
				"Homozygosity\t",
				($conservation || "-"),"\t",
				($snp || '-'),"\t",
				($maf[0] || '-'),"\t",
				($HGNC{hgnc_id} || $Gene),"\t",
				($trans_id || $Transcript),"\t",			
				"(hgvsc_biobase || \$hgvsc)","\t",
				($Consequence || '-'),"\t",
				"(hgvsc_protein || '-')","\t",
				"(aa_conservation || '-')",
				"(aa_distance || '-')",
				($interpro[0] || '-'),"\t",
				($condel || '-'),"\t",
				($sift || '-'),"\t",
				($polyphen || '-'),"\t",
				"(biobase_desc} || '-')","\t",
				($interpro[1] || '-'),"\t",
				($maf[1] || '-'),"\t",
				$sequence,"\t",
			"\n";
		}
	}
}


sub get_sequence () {
	my ($config, $chr, $position, $genotype, $sequence) = @_;
	my ($start, $end, @inter);
	$start = $position - 50;
	$inter[0] = $position - 1;
	$inter[1] = $position + 1;
	$end = $position + 50;
	
	my @alleles = split (/\//, $genotype);
	if (length($alleles[0]) > length($alleles[1])) {
		$inter[1] += (length($alleles[0]) - length($alleles[1]));
		$end = $inter[1] + 49;
#		print "$start\t$inter[0]\t$inter[1]\t$end\n";
	}
	my $slice =  $config->{sa}->fetch_by_region('chromosome', $chr, $start, $inter[0]);
	$$sequence[0] = $slice->seq();
	$slice =  $config->{sa}->fetch_by_region('chromosome', $chr, $inter[1], $end);
	$$sequence[1] = $slice->seq();
	return $sequence;
}

sub get_interprot () {
	my ($config, $Transcript, $array, $chr, $start, $end, $allele_string) = @_;
	my $slice = $config->{sa}->fetch_by_region('chromosome', $chr);
	my $vfa = $config->{reg}->get_adaptor($config->{species}, 'variation', 'variationfeature');
	my $new_vf = Bio::EnsEMBL::Variation::VariationFeature->new(
	  -start => $start,
	  -end => $end,
	  -strand => 1,
	  -slice => $slice,           # the variation must be attached to a slice
	  -adaptor => $vfa,
	  -allele_string => $allele_string,
	);
	foreach my $tv(@{$new_vf->get_all_TranscriptVariations}) {
		my $transcript_id=$tv->transcript->stable_id;
		if ($transcript_id && $transcript_id eq $Transcript) {
			my $slice_adaptor = $config->{reg}->get_adaptor(($config->{species}), 'core', 'Transcript');
			my $tr = $slice_adaptor->fetch_by_stable_id($Transcript);
			my $translation = $tr->translation;

			if($translation) {
				my $pfeatures = $translation->get_all_ProteinFeatures();
				my $aa_coordinate=&format_coords($tv->translation_start, $tv->translation_end);
				
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
						if(($aa_start >=$pfeature->start() && $aa_start<=$pfeature->end()) || ($aa_end >=$pfeature->start() && $aa_end<=$pfeature->end())){					
							if ($pfeature->interpro_ac() && $pfeature->idesc()) {
								my $ipro_id = $pfeature->interpro_ac();
								my $ipro_desc = $pfeature->idesc();
								if (!$$array[0] && !$$array[1]) {
									$$array[0] = $ipro_id;
									$$array[1] = $ipro_desc;
								}
								elsif ($$array[0] =~ /$ipro_id/ || $$array[1] =~ /$ipro_desc/) {
									next;
								}
								else {
									$$array[0] .= ','.$ipro_id;
									$$array[1] .= ','.$ipro_desc;
								}
							}
						}
					}
				}
			}
		}
	}
	return $array;			
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

sub get_rs_data() {

	my ($config, $snp, $genotype, $array, $populations, $last_number) = @_;
	my @snps = split(/,/, $snp);
	my $var = $config->{va}->fetch_by_name($snps[0]);
	my $alleles = $var->get_all_Alleles();
	my $ancestral = ($var->ancestral_allele() || '-');
	my @tmp_alleles = split (/\//, $genotype);
		
	my @samp_alleles;
	if (length ($tmp_alleles[0]) != length ($tmp_alleles[1])){
		if (length($tmp_alleles[0]) > length($tmp_alleles[1])) {
			my $length = length($tmp_alleles[1]);
			$tmp_alleles[1] = '-';
			my @bases = split ("", $tmp_alleles[0]);
			for (my $i = 0; $i < $length; $i++) {
				shift(@bases);	
			}
			$tmp_alleles[0] = join("", @bases);
		}
		elsif (length($tmp_alleles[1]) > length($tmp_alleles[0])) {
			my $length = length($tmp_alleles[0]);
			$tmp_alleles[0] = '-';
			my @bases = split ("", $tmp_alleles[1]);
			for (my $i = 0; $i < $length; $i++) {
				shift(@bases);	
			}
			$tmp_alleles[1] = join("", @bases);
		}
	}
	
	my @final_alleles;
	my $value = 0;
	my $comp = 0;
	my $times = 0;
	foreach my $allele (@{$alleles}) {
		my $allele_string = $allele->allele;
		if ($times == 0) {
			if ($allele_string eq $tmp_alleles[0] || $allele_string eq $tmp_alleles[1]) {
				push (@final_alleles, $allele);
				$value++;
			}
			elsif ($allele_string eq &get_complementary($tmp_alleles[0])|| $allele_string eq &get_complementary($tmp_alleles[1])) {
				push (@final_alleles, $allele);
				$value++;
				$comp++;
			}
			$times++;
		}
		else {
			if ($comp == 0 && $value > 0) {
				if ($allele_string eq $tmp_alleles[0] || $allele_string eq $tmp_alleles[1]) {
					push (@final_alleles, $allele);
					$value++;
				}
			}
			elsif ($comp > 0 && $value > 0) {
				if ($allele_string eq &get_complementary($tmp_alleles[0])|| $allele_string eq &get_complementary($tmp_alleles[1])) {
					push (@final_alleles, $allele);
					$value++;
				}
			}
			else {
				if ($allele_string eq $tmp_alleles[0] || $allele_string eq $tmp_alleles[1]|| $allele_string eq &get_complementary($tmp_alleles[0])|| $allele_string eq &get_complementary($tmp_alleles[1]))  {
					push (@final_alleles, $allele);
					$value++;
					$comp++
				}	
			}
		}
	}
	
	if ($value > 1) {
		my %alleles;

		foreach my $allele (@final_alleles) {
			get_maf_pop($allele, \%alleles, $last_number, $populations)
		}
		$$array[0] = '';
		$$array[1] = '';
		foreach my $print_allele (keys %alleles) {
			$$array[0] .= "$print_allele/";
			$$array[1] .= "$alleles{$print_allele}{pop}/";
		}
		$$array[0] =~ s/\/$//;
		$$array[1] =~ s/\/$//;
		$$array[3] = 1;
		return $array;
	}
	else {
		$$array[3] = 0;
		return $array;
	}
}

sub get_maf_pop() {
	my ($allele, $alleles, $last_number, $populations) = @_;
	my $allele_string = $allele->allele;
	my $string;
	
	if (!$allele->frequency || $allele->frequency < 0.05) {
		my $freq = $allele->frequency;
		$string = $allele->allele."_".($freq || "NAF");
		my $popu = $allele->population;
		my $name = $popu->{name};
		if (defined $name && !$$populations{$name}) {
			$$last_number++;
			$$populations{$name}{number} = $$last_number;
		}
		
		foreach my $known_pop (keys %{$populations}) {
			my $pop_num = $$populations{$known_pop}{number};
			if (defined $name && $known_pop eq $name && $$alleles{$string} && $$alleles{$string}{pop} !~ m/\($pop_num\)/) {
				$$alleles{$string}{pop} .= ",($pop_num)";
			}
			
			elsif (defined $name && $known_pop eq $name && !$$alleles{$string}{pop}) {
					$$alleles{$string}{pop} = "($pop_num)";
			}
		}
		return $alleles;
	}
}

sub get_complementary () {
	my $bases = shift;
	my $complementary = $bases;
	$complementary =~ s/actg/tgac/ig;
	return $complementary;
}

sub get_hgnc () {
	my $config = shift;
	my $transcript = shift;
	my $hgnc = shift;

	my $gene = $config->{ga}->fetch_by_transcript_stable_id($transcript);

	$$hgnc{hgnc_id} = $gene->external_name;

	my $HGNC_Desc = $gene->description;
	
	$HGNC_Desc =~ s/\s+/_/g;
	$HGNC_Desc =~ s/\[.*\]//g;
	$HGNC_Desc =~ s/_$//;
	$$hgnc{hgnc_desc} = $HGNC_Desc;
	
	return $hgnc;
	
}
 
sub get_conservation_score {
	my $contador=0;
	my ($var1,$var2,$var3, $conservation)=@_;
	$conservation = '';

	my $nslice = $config->{sa}->fetch_by_region('toplevel', $_[0], $_[1], $_[2]+1);
	my $display_size = ($nslice->end) - $nslice->start + 1;
	my $scores = $config->{cs_adaptor}->fetch_all_by_MethodLinkSpeciesSet_Slice($config ->{mlss}, $nslice, $display_size);

	foreach my $score (@$scores) {
		if (defined $score->diff_score) {
				$conservation .= sprintf ("%.4f", $score->diff_score).",";
		}	
	}
	$conservation =~ s/,$//;
	return $conservation;
}

sub get_transcript_translation () {
	my ($config, $transcripts) = @_;
 	my $bed_file_handle = $config->{bed_file_handle}; 
	while (my $line = <$bed_file_handle>) {
		chomp $line;
		my ($gene_id, $transcript_id, $NM_id) = split (/\s+/, $line);
		if (!defined $$transcripts{$transcript_id}{name}) {
			$$transcripts{$transcript_id}{name} = $NM_id;
			$$transcripts{$transcript_id}{gene} = $gene_id;
		}
	}
	return $transcripts;
}

sub configure {
	my $args = shift;
	
	my $config = {};
	
	GetOptions(
		$config,
		'help',
		
		# input options,
		'config=s',
		'input_file=s',
		'format=s',
		'bed_file=s',
		
		# DB options
		'species=s',
		'registry=s',
		'host=s',
		'user=s',
		'port=s',
		'password=s',
		'db_version=i',
		'genomes',
		
		# runtime options
		'most_severe',
		'buffer_size=i',
		'chunk_size=s',
		'check_ref',
		'check_existing=i',
		'failed=i',
		'whole_genome',
		'tmp_dir=s',
		'gp',
		
		# output options
		'output_file=s',
		'terms=s',
		'verbose',
		'quiet',
		'coding_only',
		'protein',
		'hgnc',
		'hgvs',
		'sift=s',
		'polyphen=s',
		'condel=s',
	);
	
	# print usage message if requested or no args supplied
	if(defined($config->{help}) || !$args) {
		&usage;
		exit(0);
	}
	
	# summarise options if verbose
	if(defined $config->{verbose}) {
		my $header =<<INTRO;
#---------------------------------------------#
# CUSTOMIZED ENSEMBL VARIANT EFFECT PREDICTOR #
#---------------------------------------------#

version 1.0

By sistemas Genomicos

Configuration options:

INTRO
		print $header;
		
		my $max_length = (sort {$a <=> $b} map {length($_)} keys %$config)[-1];
		
		foreach my $key(sort keys %$config) {
			print $key.(' ' x (($max_length - length($key)) + 4)).$config->{$key}."\n";
		}
		
		print "\n".("-" x 20)."\n\n";
	}
	
	print "Creating tmp file ensembl64.tmp... ".localtime()."\n";
	`touch ensembl64.tmp`;

	# set defaults
	
	$config->{species} ||= "homo_sapiens";
	$config->{host}    ||= 'ensembldb.ensembl.org';
	$config->{port}    ||= 5306;
	$config->{user}         ||= 'anonymous';
	$config->{output_file}  ||= "variant_effect_output.txt";
	$config->{tmp_file}      = "ensembl64.tmp";
	
	# connect to databases
	$config->{reg} = &connect_to_dbs($config);
	
	# get input file handles
	$config->{in_file_handle} = &get_in_file_handle($config);
	
	# get bed file handles
	$config->{bed_file_handle} = &get_bed_file_handle($config);	
	
	# get tmp file handles
	$config->{tmp_file_handle} = &get_tmp_file_handle($config);	
	
	# configure output file
	$config->{out_file_handle} = &get_out_file_handle($config);
	
	#getting adaptors
	
	$config->{va} = $config->{reg}->get_adaptor($config->{species}, 'variation', 'variation');
	$config->{sa} = $config->{reg}->get_adaptor($config->{species}, 'core', 'slice');
	$config->{mlss_adaptor} = $config->{reg}->get_adaptor("Multi", "compara", "MethodLinkSpeciesSet");
	$config ->{mlss} = $config->{mlss_adaptor}->fetch_by_method_link_type_species_set_name("GERP_CONSERVATION_SCORE", "mammals");
	$config->{cs_adaptor} = $config->{reg}->get_adaptor("Multi", 'compara', 'ConservationScore');
	$config->{ga} = $config->{reg}->get_adaptor($config->{species}, 'core', 'gene');
	$config->{tva} = $config->{reg}->get_adaptor($config->{species}, 'variation', 'transcriptvariation');
	
	return $config;
}

sub connect_to_dbs {
    my $config = shift;
    
    # get registry
    my $reg = 'Bio::EnsEMBL::Registry';
	$reg->load_registry_from_db(
		-host       => $config->{host},
		-user       => $config->{user},
		-pass       => $config->{password},
		-port       => $config->{port},
		-db_version => $config->{db_version},
		-species    => $config->{species} =~ /^[a-z]+\_[a-z]+/i ? $config->{species} : undef,
		-verbose    => $config->{verbose},
		-no_cache   => $config->{no_slice_cache},
	);
         
	eval { $reg->set_reconnect_when_lost() };
    
    return $reg;
}

sub usage {
	my $usage =<<END;
#---------------------------------------------#
# CUSTOMIZED ENSEMBL VARIANT EFFECT PREDICTOR #
#---------------------------------------------#

version 1.0

By Sistemas Genomicos

Usage:
perl variant_effect_predictor_Ensembl64.pl [arguments]

Options
=======

--help                 Display this message and quit

-i | --input_file      Input file - vcf format.

-o | --output_file     Output file. Write to STDOUT by specifying -o STDOUT - this
                       will force --quiet [default: "variant_effect_output.txt"]
                       
--species [species]    Species to use [default: "human"]
                                       
--host                 Manually define database host [default: "ensembldb.ensembl.org"]
--user 		           Database username [default: "anonymous"]
--port                 Database port [default: 5306]
--password             Database password [default: no password]

--bed_file			   BED file with transcripts ID

END

	print $usage;
}
# gets file handle for input
sub get_in_file_handle {
    my $config = shift;

    # define the filehandle to read input from
    my $in_file_handle = new FileHandle;
    
    if(defined($config->{input_file})) {
        
        # check defined input file exists
        die("ERROR: Could not find input file ", $config->{input_file}, "\n") unless -e $config->{input_file};
        
        if($config->{input_file} =~ /\.gz$/){
            $in_file_handle->open($config->{compress}." ". $config->{input_file} . " | " ) or die("ERROR: Could not read from input file ", $config->{input_file}, "\n");
        }
        else {
            $in_file_handle->open( $config->{input_file} ) or die("ERROR: Could not read from input file ", $config->{input_file}, "\n");
        }
    }
    
    # no file specified - try to read data off command line
    else {
        $in_file_handle = 'STDIN';
        debug("Reading input from STDIN (or maybe you forgot to specify an input file?)...") unless defined $config->{quiet};
    }
    
    return $in_file_handle;
}

# gets temp file handle for input
sub get_tmp_file_handle {
    my $config = shift;

    # define the filehandle to read input from
    my $tmp_file_handle = new FileHandle;
    
    if(defined($config->{tmp_file})) {
        
        # check defined input file exists
        die("ERROR: Could not find input file ", $config->{tmp_file}, "\n") unless -e $config->{tmp_file};
        
        if($config->{tmp_file} =~ /\.gz$/){
            $tmp_file_handle->open($config->{compress}." ". $config->{tmp_file} . " | " ) or die("ERROR: Could not read from input file ", $config->{tmp_file}, "\n");
        }
        else {
            $tmp_file_handle->open( $config->{tmp_file} ) or die("ERROR: Could not read from input file ", $config->{tmp_file}, "\n");
        }
    }
    
    # no file specified - try to read data off command line
    else {
        $tmp_file_handle = 'STDIN';
        debug("Reading input from STDIN (or maybe you forgot to specify an input file?)...") unless defined $config->{quiet};
    }
    
    return $tmp_file_handle;
}

sub get_bed_file_handle {
    my $config = shift;

    # define the filehandle to read input from
    my $bed_file_handle = new FileHandle;
    
    if(defined($config->{bed_file})) {
        
        # check defined input file exists
        die("ERROR: Could not find input file ", $config->{bed_file}, "\n") unless -e $config->{bed_file};
        
        if($config->{bed_file} =~ /\.gz$/){
            $bed_file_handle->open($config->{compress}." ". $config->{bed_file} . " | " ) or die("ERROR: Could not read from input file ", $config->{bed_file}, "\n");
        }
        else {
            $bed_file_handle->open( $config->{bed_file} ) or die("ERROR: Could not read from input file ", $config->{bed_file}, "\n");
        }
    }
    
    # no file specified - try to read data off command line
    else {
        $bed_file_handle = 'STDIN';
        debug("Reading input from STDIN (or maybe you forgot to specify an input file?)...") unless defined $config->{quiet};
    }
    
    return $bed_file_handle;
}

sub get_out_file_handle {
	my $config = shift;
	
	# define filehandle to write to
	my $out_file_handle = new FileHandle;
	$out_file_handle->open(">".$config->{output_file}) or die("ERROR: Could not write to output file ", $config->{output_file}, "\n");
	
	# make header
	my $time = &get_time;
	my $core_mca = $config->{reg}->get_adaptor($config->{species}, 'core', 'metacontainer');
	my $db_string = $core_mca->dbc->dbname." on ".$core_mca->dbc->host if defined $core_mca;
	my $version_string =
		"Using API version ".$config->{reg}->software_version.
		", DB version ".(defined $core_mca && $core_mca->get_schema_version ? $core_mca->get_schema_version : '?');
	
	my $header =<<HEAD;
## CUSTOMIZED ENSEMBL VARIANT EFFECT PREDICTOR v1.0
## Output produced at $time
## Connected to $db_string
## $version_string
## Extra column keys:
## HGNC     : HGNC gene identifier
## ENSP     : Ensembl protein identifer
## HGVSc    : HGVS coding sequence name
## HGVSp    : HGVS protein sequence name
## SIFT     : SIFT prediction
## PolyPhen : PolyPhen prediction
## Condel   : Condel SIFT/PolyPhen consensus prediction
## InterPro_ID :Protein Domain ID from InterPro database
## InterPro_Desc :Protein Domain ID description InterPro database
## MAF	:Minor Allele Frequency < 0.05 
## Pop :Population in which MAF < 0.05

HEAD
	
	# add headers
	print $out_file_handle $header;
	
	# add column headers
	print $out_file_handle join "\t", qw(
		#Uploaded_variation
		Location
		Ensembl_ID
		Transcript
		HGNC_ID
		HGNC_Description
		Consequence
		cDNA_position
		CDS_position
		Protein_position
		Amino_acids
		Codons
		Existing_variation
		HGVSc
		Conservation_value
		Condel_value
		SIFT_value
		Polyphen_value
		InterPro_ID
		InterPro_Desc
		MAF
		Pop
	);

	print $out_file_handle "\n";
	
	return $out_file_handle;
}

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
