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
use Bio::EnsEMBL::DBSQL::OntologyDBAdaptor;

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
		my $command = "variant_effect_predictor_64.pl -i $input --format pileup -o $output --sift b --polyphen b --condel b --regulatory --hgnc  --hgvs --no_intergenic --force_overwrite --species $species --host $host --port $port --user $user $pass --check_existing";
		
		print "Launching command:\n$command\n".localtime()."\n";
	
		`variant_effect_predictor_64.pl -i $input --format guess -o $output --sift b --polyphen b --condel b --regulatory --hgnc  --hgvs 	--no_intergenic --force_overwrite --species $species --host $host --port $port --user $user $pass --check_existing > ensembl64_tmp.log`;
	}
	else {
		my $command = "variant_effect_predictor_64.pl -i $input --format pileup -o $output --regulatory --hgnc  --hgvs --no_intergenic --force_overwrite --species $species --host $host --port $port --user $user $pass --check_existing";
		
		print "Launching command:\n$command\n".localtime()."\n";
	
		`variant_effect_predictor_64.pl -i $input --format guess -o $output --regulatory --hgnc  --hgvs 	--no_intergenic --force_overwrite --species $species --host $host --port $port --user $user $pass --check_existing  > ensembl64_tmp.log`;
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

		my $index;
		if (defined $config->{bed_file}) {
			checking_transcript(\%transcripts, $line, \$index);
			#print "$index\n";
			if ($index == 0) {
				next;
			}
			else {
				print_output($config, $line, \%populations, \$last_number, $out_file_handle, \%transcripts, $index);
			}
		}
		else {
			$index = 0;
			print_output($config, $line, \%populations, \$last_number, $out_file_handle, \%transcripts, $index);
		}
		
		$vf_count++;
		debug("Processed $vf_count variants") if $vf_count =~ /0$/ && defined($config->{verbose});
	}
	print $out_file_handle "\n\n#Code\tPop_Name\n";
	
	foreach my $number (sort {$a <=> $b} keys %{$populations{numbers}}) {
		print $out_file_handle "$number\t$populations{numbers}{$number}{name}\n",
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
	#print "$$value\t";
}

sub print_output {
	my ($config, $line, $populations, $last_number, $out_file_handle, $transcripts, $index) = @_;
	
	my ($Uploaded_variation, $Location, $Allele, $Gene, $Transcript, $Feature_type, $Consequence, $cDNA_position, $CDS_position,$Protein_position, $Amino_acids, $Codons, $Existing_variation, $Extra) = split (/\s+/, $line);
	
	my $trans_id;
	
	if (defined $config->{bed_file}) {
		if (defined $$transcripts{$Transcript}) {
			$trans_id = $$transcripts{$Transcript}{name};
		}
	}
		
	my @rs = split (/,/, $Existing_variation);
	my ($genotype, @maf, %HGNC, $condel, $polyphen, $sift, $hgvsc, $hgvsp, @coord, $end, @interpro, @sequence, $sequence);
	
	foreach my $snp (@rs) {
		$genotype = (split/_/, $Uploaded_variation)[2];
		$maf[3] = 1;
		if ($snp =~ m/rs/i) {
			get_rs_data($config, $snp, $genotype, \@maf, $populations, $last_number);
		}
		if ($maf[3] == 0) {
			$snp = '-';
			$maf[0] = '-';
			$maf[1] = '-';
		}
		$Location =~ s/\d+://;
		my @extra = split (/;/, $Extra);
		my $NM_id;
		my $prot_id;
		foreach my $field (@extra) {
			if ($field =~ m/HGVSc/) {
				if ($index == 1) {
					$field =~ s/HGVSc=//;
					$hgvsc = $field;
					my @data = split (/.\d+:/, $field);
					$NM_id = $$transcripts{$data[0]}{name};
					$prot_id = $$transcripts{$data[0]}{prot};
					$hgvsc = $NM_id.":".$data[1];
				}
				else {
					$hgvsc = $field;
					$hgvsc =~ s/HGVSc=//;
				}
			}
			elsif ($field =~ m/HGVSp/) {
				if ($index == 1) {
					if ($field !~ m/ENST/) {
						$field =~ s/HGVSp=//;
						$hgvsp = $field;
						my @data = split (/.\d+:/, $field);
						$hgvsp = $prot_id.":".$data[1];
					}
				}
				else {
					$hgvsp = $field;
					$hgvsp =~ s/HGVSp=//;
				}
			}
			elsif ($field =~ m/PolyPhen/) {
				$polyphen = $field;
				$polyphen =~ s/PolyPhen=//;
				$polyphen =~ s/\(/_\(/;
			}
			elsif ($field =~ m/Condel/) {
				$condel = $field;
				$condel =~ s/Condel=//;
				$condel =~ s/\(/_\(/;
			} 
			elsif ($field =~ m/SIFT/) {
				$sift = $field;
				$sift =~ s/SIFT=//;
				$sift =~ s/\(/_\(/;
			} 
		}
		
		undef @extra;
		
		if ($Transcript =~ m/ENST/ && $Extra =~ m/HGNC/) {
			get_hgnc($config, $Transcript, \%HGNC);
		}
		
		my $aa_distance;
		if (defined $hgvsp) {
			$aa_distance = &get_aa_distance ($hgvsp);
		}
		
		@coord = split (/\_/, $Uploaded_variation);
		$end = $coord[1] + (length ($coord[2]) - 3);
		my $conservation = &get_conservation_score ($coord[0], $coord[1], $end);
		get_interprot ($config, $Transcript, \@interpro, $coord[0], $coord[1], $end, $genotype);
		get_sequence ($config, $coord[0], $coord[1], $genotype, \@sequence);
		$sequence  = $sequence[0]."[".$genotype."]".$sequence[1];
		my ($depth, $rate, @Info);
		get_depth($config, $Uploaded_variation, \$depth, \$rate, \@coord, \@Info);
		my $SB = &get_SB($config, $Uploaded_variation);
		my $aa_conservation= &get_aa_conservation($config);
		my $biobase_desc = &get_biobase_info($coord[0], $coord[1]);
		
		
		
		print $out_file_handle ($HGNC{hgnc_id} || $Gene),"\t",
			$Info[0],"\t",
			$Info[1],"\t",
			$Info[2],"\t",
			$Info[3],"\t",
			$Info[4],"\t",
			($depth || '-'),"\t",
			($rate || '-'),"\t",
			($SB || '-'),"\t",
			($conservation || "-"),"\t",
			($snp || '-'),"\t",
			($maf[0] || '-'),"\t",
			($trans_id || $Transcript),"\t",			
			($hgvsc || '-'),"\t",				
			($Consequence || '-'),"\t",
			($hgvsp || '-'),"\t",
			($aa_conservation || '-'),"\t",
			($aa_distance || '-'),"\t",
			($interpro[0] || '-'),"\t",
			($condel || '-'),"\t",
			($sift || '-'),"\t",
			($polyphen || '-'),"\t",
			($biobase_desc || "-\t-"),"\t",
			($interpro[1] || '-'),"\t",
			($maf[1] || '-'),"\t",
			$sequence,"\t",
		"\n";
	}
}

#FUNCION A ACTUALIZAR
sub get_biobase_info () {
	my ($chr, $position) = @_;
	my $bio_file_handle = &get_bio_file_handle($chr);
	my $desc = "-\t-";
	while (my $line = <$bio_file_handle>) {
		chomp $line;
		my @data = split (/\t/, $line);
		my ($chr, $start, $end) = split (/[:-]/, $data[4]);
		if (!$end or $end eq '' or $end eq '+') { $end = $start + 1;}
		if ($position >= $start && $position <= $end) {
			
			$desc = $data[1]."/".$data[9]."/".$data[7];
			$desc =~ s/\s+/_/g;

			if($data[11]!~/A-Z]/)
			{
				$desc .= "\thttp://www.ncbi.nlm.nih.gov/pubmed/".$data[11];
			}
			else
			{
				$desc .= "\t".$data[11];
			}
			
			last;
		}
	}
	return $desc;
}
#######################
	
sub get_aa_conservation () {
	my ($config) = @_;
	my $aa_conservation = 'NA';
	return $aa_conservation;
}

sub get_aa_distance () {
	my ($hgvsp) = @_;	
	my @aminoacids = split (/:p\./, $hgvsp);
	my @aa = split (/\d+/, $aminoacids[1]);
	$aa[1] =~ s/fs.*//g;
	my $dist = '-';
	if ($aa[1] !~ m/X/) {
		$dist = &get_distance ($aa[0], $aa[1]);
	}
	return $dist;
}

sub get_SB () {
	my ($config, $variation) = @_;	
	my $match = 'SB=';
	my $value = &analysing_field_8 ($config, $variation, $match);
	return $value;
}

sub get_depth() {
	my ($config, $variation, $depth, $rate, $coord, $info) = @_;
	my $column = '$10';
	my $line = &using_awk_vcf ($config, $variation, $column);
	my @data = split (/:/, $line);
	my @Depth = split (/,/, $data[1]);
	my $ratio;
	if (not defined $Depth[1]) {
		my $match = 'Dels=';
		my $value = &analysing_field_8 ($config, $variation, $match);
		$ratio = $value;
		$match = 'DP=';
		$value = &analysing_field_8 ($config, $variation, $match);
		$$depth = $value;
	}
	else {
		$$depth = $Depth[0] + $Depth[1]; 
		if ($$depth > 0 ) {
			$ratio = sprintf("%.2f", ($Depth[1] / $$depth));
		}
		else {
			$ratio = 0;
		}
	}
	
	$$info[0] = $$coord[0];
	$$info[1] = $$coord[1];
	my $ref = (split (/\//, $$coord[2]))[0];
	my $var = (split (/\//, $$coord[2]))[1];
	$$info[2] = $ref;
	$$info[3] = $var;
	
	if ($ratio == 0) {
		$$rate = "NE_(".$ratio.")";
		$$info[4] = $ref."/".$var;
	}
	elsif ($ratio > 0 && $ratio <= 0.12) {
		$$rate = "P_Homo_ref_(".$ratio.")";
		$$info[4] = $ref."/".$ref;
	}
	elsif ($ratio > 0.12 && $ratio < 0.35) {
		$$rate = "UNC_FP_(".$ratio.")";
		$$info[4] = $ref."/".$var;
	}
	elsif ($ratio >= 0.35 && $ratio < 0.65) {
		$$rate = "P_Hetero_(".$ratio.")";
		$$info[4] = $ref."/".$var;
	}
	elsif ($ratio > 0.65 && $ratio < 0.85) {
		$$rate = "UNC_Het_(".$ratio.")";
		$$info[4] = $ref."/".$var;
	}
	else {
		$$rate = "P_Homo_var_(".$ratio.")";
		$$info[4] = $var."/".$var;
	}
}

sub analysing_field_8 (){
	my ($config, $variation, $match) = @_;
	my $column = '$8';
	my $value;
	my $line = &using_awk_vcf ($config, $variation,$column);
	my @fields =  split (/;/, $line);
	foreach my $field (@fields) {
		if ($field =~ m/$match/) {
			$field =~ s/$match//;
			$value = $field;
		}
	}
	return $value;
}

sub using_awk_vcf () {
	my ($config, $variation, $column) = @_;
	my @loc = split (/_/, $variation);
	my @alleles = split (/\//, $loc[2]);
	my $command = 'awk \'(('.'$1 == '.$loc[0].') || ('.'$1 == "chr'.$loc[0].'")) && ('.'$2 == '.$loc[1].') && ('.'$4 == "'.$alleles[0].'") && ('.'$5 == "'.$alleles[1].'") {print '.$column.'}\' '.$config->{parsed_file};
	my $line = `$command`;
	return $line;
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
			get_maf_pop($allele, \%alleles, $last_number, $populations);
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
		my $freq;
		if ($allele->frequency) {
			$freq = sprintf ("%.4f", $allele->frequency);
		}
		$string = $allele->allele."_".($freq || "NAF");
		my $popu = $allele->population;
		my $name = $popu->{name};
		if (defined $name && !$$populations{pops}{$name}) {
			$$last_number++;
			$$populations{numbers}{$$last_number}{name} = $name;
			$$populations{pops}{$name}{number} = $$last_number;
		}
		
		foreach my $known_pop (keys %{$$populations{pops}}) {
			my $pop_num = $$populations{pops}{$known_pop}{number};
			if (defined $name && $known_pop eq $name && $$alleles{$string} && $$alleles{$string}{pop} !~ m/$pop_num/) {
				$$alleles{$string}{pop} .= ",$pop_num";
			}
			
			elsif (defined $name && $known_pop eq $name && !$$alleles{$string}{pop}) {
					$$alleles{$string}{pop} = "$pop_num";
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
}
 
sub get_conservation_score {
	my $contador=0;
	my ($var1,$var2,$var3)=@_;
	my $conservation = '';

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
		my ($gene_id, $transcript_id, $NM_id, $prot_id) = split (/\s+/, $line);
		if (!defined $$transcripts{$transcript_id}{name}) {
			$$transcripts{$transcript_id}{name} = $NM_id;
			$$transcripts{$transcript_id}{gene} = $gene_id;
			$$transcripts{$transcript_id}{prot} = $prot_id;
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
		'parsed_file=s',
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
--parsed_file			Parsed-indel input file - vcf format.

-o | --output_file     Output file. Write to STDOUT by specifying -o STDOUT - this
                       will force --quiet [default: "variant_effect_output.txt"]
                       
--species [species]    Species to use [default: "human"]
                                       
--host                 Manually define database host [default: "ensembldb.ensembl.org"]
--user 		          Database username [default: "anonymous"]
--port                 Database port [default: 5306]
--password             Database password [default: no password]

--bed_file		   BED file with transcripts ID

END

	print $usage;
}

# gets file handle for input
sub get_bio_file_handle {
	my $chr = shift;
	my $bio_file_handle = new FileHandle;
	my $file = "/share/references/biobase/bioBaseSortUniques_chr$chr.txt";
	$bio_file_handle->open($file) or die("ERROR: Could not read from input file: $file\n");
	return $bio_file_handle;
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
	
	# add column headers
	print $out_file_handle join "\t", qw(
		#HGNC_ID
		Chromosome
		Position
		Ref_Allele
		Var_Allele
		Genotype
		Depth
		Ratio_Var/Depth
		Strand_Bias
		Base_Conservation_Score
		Existing_Variation
		MAF
		Transcript
		HGVSc
		Consequence
		HGVSp
		AA_Conservation_Score
		AA_Distance_Grantham
		InterPro_ID
		Condel_Value
		SIFT_Value
		Polyphen_Value
		Biobase_Description
		Biobase_Link
		InterPro_Desc
		Populations_MAF
		Flanking_Sequence
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

sub get_distance () {
	my ($aa1, $aa2) = @_;
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

	my $dist = $Distances{$aa1}{$aa2};
	my $distance;
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
