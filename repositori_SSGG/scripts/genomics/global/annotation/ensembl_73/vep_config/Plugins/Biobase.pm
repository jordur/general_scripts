=head1 LICENSE
                                                                                                                     
 Copyright (c) 2013 Sistemas Genómicos S.L - All rights reserved.                                                                      
                                                                                                                     
 This software is distributed under a modified Apache license.                                                                                                                      
                                                                                                                     
=head1 CONTACT                                                                                                       

 Guillermo Marco Puche <guillermo.marco@sistemasgenomicos.com>
    
=cut

=head1 NAME

 biobase

=head1 SYNOPSIS

 Gets variant information from HGMD.

=head1 DESCRIPTION
 
 Returns two fields. HGMD_info and Related_Publication.

=cut

package Biobase;

use strict;
#use warnings;

use FileHandle;
use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Data::Dumper;

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
    return '73';
}

sub feature_types {
    return ['Feature','Intergenic'];
}

sub get_header_info {
    return { 
    	HGMD_info => "HGMD database description.",
    	Related_publication => "HGMD database link."
    };
}

sub run {
	
	### Def Handler's vars
	my $self = shift;
    my $vf = shift;
    my $line_hash = shift;
    my $gene;
	my $gene_desc = "-";
    
	### Def Handlers
    my $config = $self->{config};
    my $ind_cols = $config->{ind_cols};
    my $line = $vf->{base_variation_feature_overlap}->{base_variation_feature}->{_line};
    
    
    
    ### Def vars
	my ($chr, $pos, $ref_allele, $var_allele, $variation);
	
	## Split VCF input line and get information
	my @split_line = split ("\t", $line);
	$chr = $split_line[0];		
	$pos = $split_line[1];
	$ref_allele = $split_line[3];
	# $var_allele must be obtained from currently parse line hash, to keep functionality in case of multiple alternative alleles 
	my $var_allele = $line_hash->{Allele};
	
	## Get Biobase file
	my $bio_file_handle = &get_bio_file_handle($chr);
	my $desc = "-\t-";
	
	# Load Biobase information for chromosome into memory
	my @f=<$bio_file_handle>;

	# Variables for working with changes and their reverse-complementaries
	#my $ref_allele="";
    my $ref_allele2="";
    #my $var_allele="";
    my $var_allele2="";

	my $once = 0;
	
	## Transcription and adaptation of code of sg_extraerBiobaseInformation.pl... Uuuufffff... 
	## I must keep in mind, to never again undergo such a task...
	$ref_allele2 = &complementario($ref_allele);
	$var_allele2 = &complementario($var_allele);
	
	my ($position, $positionA, $positionB);
	my $count = 0;

	my $SNVs = "NO";
	my $entry = 0;
	#@annotationExplode=split("\t",$line);
	
	# Check if the variant is a SNV or indel
	if(length($ref_allele) == length ($var_allele)) {
		# SNV
		$position=$pos;
		$SNVs="YES";
	} else {
		# indel
		$positionA = $pos; 
		# Previous annotated files contained indel length
		$positionB = $pos + abs(length($ref_allele) - length ($var_allele));
		$SNVs="NO";
	}
	#$line=~s/SNV\t-/SNV\t0/;
	for(my $k=$count;$k<scalar(@f);$k++) {
		chomp($f[$k]);

		my $biostart=0;
		my $bioend=0;
		my @explode1=split("\t",$f[$k]);
		my @explode2=split(":",$explode1[4]);
		if($explode2[1] !~/-/) {
			$biostart=$explode2[1];
			$bioend=$explode2[1];
		} else {
			my @explode3=split("-",$explode2[1]);
			$biostart=$explode3[0];
			$bioend=$explode3[1];
		}
		my $strand=$explode2[2];
		if($SNVs eq "YES") {
			if($position >= $biostart && $position <= $bioend) {
				if($explode1[8] eq "-" || $explode1[8] eq "null") { ##ESA COLUMNA ESTA VACIA
					$entry=1;
					$count=$k;
					
					if($once==0){
						$once=1;
					}
					if($explode1[11]!~/[A-Z]/) {##PUBMED IDENTIFICATION
						$desc = $explode1[1]."/".$explode1[9]."/".$explode1[17]."/".$explode1[7];
						$desc =~ s/\s+/_/g;
						$desc .= "\thttp://www.ncbi.nlm.nih.gov/pubmed/".$explode1[11];
					} else {
						$desc = $explode1[1]."/".$explode1[9]."/".$explode1[17]."/".$explode1[7];
						$desc =~ s/\s+/_/g;
						$desc .= "\t".$explode1[11];
					}
				} else {
					$explode1[8]=~s/\[/\t/;
					$explode1[8]=~s/\]/\t/;
					$explode1[8]=~s/\//-/;
					my @SNVspecific = split("\t",$explode1[8]);
					if($SNVspecific[1]=~/-/) {
						my @SNVspecific2=split("-",$SNVspecific[1]);
						if($strand eq "+") {
							if($ref_allele eq $SNVspecific2[0] && $var_allele eq $SNVspecific2[1]) {
								$entry=1;
								$count=$k;

								if($once==0) {
									$once=1;
								}
								if($explode1[11]!~/A-Z]/){ ## PUBMED IDENTIFICATION
									$desc = $explode1[1]."/".$explode1[9]."/".$explode1[17]."/".$explode1[7];
									$desc =~ s/\s+/_/g;
									$desc .= "\thttp://www.ncbi.nlm.nih.gov/pubmed/".$explode1[11];
								} else {
									$desc = $explode1[1]."/".$explode1[9]."/".$explode1[17]."/".$explode1[7];
									$desc =~ s/\s+/_/g;
									$desc .= "\t".$explode1[11];
								}
							} else {
								###BECAUSE BIOBASE DATABASE AREN'T STANDARD HGVS, IF SOME BASE IS EQUAL TO REFERENCE IS OK 
								if($ref_allele  eq $SNVspecific2[1] && $var_allele eq $SNVspecific2[0]) {
									$entry=1;
									$count=$k;
									if($once==0) {
										$once=1;
									} 
									if($explode1[11]!~/A-Z]/) {
										## PUBMED IDENTIFICATION
										$desc = $explode1[1]."/".$explode1[9]."/".$explode1[17]."/".$explode1[7];
										$desc =~ s/\s+/_/g;
										$desc .= "\thttp://www.ncbi.nlm.nih.gov/pubmed/".$explode1[11];
									} else {
										$desc = $explode1[1]."/".$explode1[9]."/".$explode1[17]."/".$explode1[7];
										$desc =~ s/\s+/_/g;
										$desc .= "\t".$explode1[11];
									}
								}
							}
						} else {
							if($ref_allele2 eq $SNVspecific2[0] && $var_allele2 eq $SNVspecific2[1]) {
								#COMPLEMENTARY STRAND
								$entry=1;
								$count=$k;
								if($once==0) {
									$once=1;
								}
								### OJOOOO!!! la linea siguiente era erronea en el codigo fuente de sg_extraerBiobaseInformation
								## se hacia referencia a la variable "explode", que no aparece definida en ningun caso. He creido 
								## conveniente sustituirla por "explode1", aunque el codigo entero queda puesto en duda hasta su 
								## revision exhaustiva. 
								## Hacer notar que este tipo de errores se evitan simplemente incluyendo la libreria strict
								if($explode1[11]!~/A-Z]/) {
									$desc = $explode1[1]."/".$explode1[9]."/".$explode1[17]."/".$explode1[7];
									$desc =~ s/\s+/_/g;
									$desc .= "\thttp://www.ncbi.nlm.nih.gov/pubmed/".$explode1[11];
								} else {
									$desc = $explode1[1]."/".$explode1[9]."/".$explode1[17]."/".$explode1[7];
									$desc =~ s/\s+/_/g;
									$desc .= "\t".$explode1[11];
								}
							} else {
								###BECAUSE BIOBASE DATABASE AREN'T STANDARD HGVS, IF SOME BASE IS EQUAL TO REFERENCE IS OK
								if($ref_allele2 eq $SNVspecific2[1] && $var_allele2 eq $SNVspecific2[0]) {
									$entry=1;
									$count=$k;
									
									if($once==0) {
										$once=1;
									}
									# Aqui tambien se hacia uso de $explode, sin que se hubiera definido previamente
									if($explode1[11]!~/A-Z]/) {
										$desc = $explode1[1]."/".$explode1[9]."/".$explode1[17]."/".$explode1[7];
										$desc =~ s/\s+/_/g;
										$desc .= "\thttp://www.ncbi.nlm.nih.gov/pubmed/".$explode1[11];
									} else {
										$desc = $explode1[1]."/".$explode1[9]."/".$explode1[17]."/".$explode1[7];
										$desc =~ s/\s+/_/g;
										$desc .= "\t".$explode1[11];
									}
								}
							}
						}
					}
				}
			}
		} else {
			my $entryIndels=0;
			if($explode1[scalar(@explode1)-1]=~/Indels/) {
				my $start_indel=0;
				my $end_indel=0;
				if($explode2[1]=~/-/) {
					my @explode3=split("-",$explode2[1]);
					$start_indel=$explode3[0];
					$end_indel=$explode3[1];
				} else {
					$start_indel=$explode2[1]; ## OFTEN DELECTIONS
				}
				
				if($end_indel!=0) {
					if($positionA ==  $start_indel || $positionA == $start_indel -1 ) {
						$entryIndels=1;
					}
				} else {
					## SOMETIMES THE BIOBASE DELECTION WITH ONLY START COORDINATE HAVE REFERENCE TO NEXT NUCLEOTIDE
					if($start_indel == $positionA || $start_indel -1 == $positionA) {
						$entryIndels=1;
					}
				}
			}
			if($explode1[18] eq "Indels" && $entryIndels==1) {
				if($explode1[8] eq "-" || $explode1[8] eq "null") { ##COLUMNA VACIA
					$entry=1;
					$count=$k;
					if($once==0) {
						$once=1;
					}
					# Aqui tambien se hacia uso de $explode, sin que se hubiera definido previamente
					if($explode1[11]!~/A-Z]/) {
						$desc = $explode1[1]."/".$explode1[9]."/".$explode1[17]."/".$explode1[7];
						$desc =~ s/\s+/_/g;
						$desc .= "\thttp://www.ncbi.nlm.nih.gov/pubmed/".$explode1[11];
					} else {
						$desc = $explode1[1]."/".$explode1[9]."/".$explode1[17]."/".$explode1[7];
						$desc =~ s/\s+/_/g;
						$desc .= "\t".$explode1[11];
					}
				} else {
					$explode1[8]=~s/\[/\t/;
					$explode1[8]=~s/\]/\t/;
					my @indelSpecific_tmp=split("\t",$explode1[8]);
					my @indelSpecific=split("/",$indelSpecific_tmp[1]);
					my $len_0;
					my $len_1;
					if($indelSpecific[0] eq "-") {
						$len_0=0;
					} else {
						$len_0=length($indelSpecific[0]);
					}

					if($indelSpecific[1] eq "-") {
						$len_1=0;
					} else {
						$len_1=length($indelSpecific[1]);
					}
					$indelSpecific[0]=~s/-//;
					$indelSpecific[1]=~s/-//;
					if(($len_0-$len_1)!=0) {
						my $indel_type2=$indelSpecific[0];
						$indel_type2=~s/$indelSpecific[1]//;
						if(($len_0-$len_1)>0) {
							$variation=$ref_allele;
							$variation=~s/$var_allele//;
							if($strand eq "-") {
								$variation=&complementario($variation);
							}
							if($indel_type2 eq $variation) {
								$entry=1;
								$count=$k;
								
								if($once==0) {
									$once=1;
								}
								# Aqui tambien se hacia uso de $explode, sin que se hubiera definido previamente
								if($explode1[11]!~/A-Z]/) {
									$desc = $explode1[1]."/".$explode1[9]."/".$explode1[17]."/".$explode1[7];
									$desc =~ s/\s+/_/g;
									$desc .= "\thttp://www.ncbi.nlm.nih.gov/pubmed/".$explode1[11];
								} else {
									$desc = $explode1[1]."/".$explode1[9]."/".$explode1[17]."/".$explode1[7];
									$desc =~ s/\s+/_/g;
									$desc .= "\t".$explode1[11];
								}
							}
						}
						if($len_0-$len_1 < 0) {
							$indel_type2==$indelSpecific[1];
							$indel_type2=~s/$indelSpecific[0]//;
							$variation=$var_allele;
							$variation=~s/$ref_allele//;
							if($strand eq "-") {
								$variation= &complementario($variation);
							}
							if($indel_type2 eq $variation) {
								$entry=1;
								$count=$k;
								if($once==0) {
									$once=1;
								}
								# Aqui tambien se hacia uso de $explode, sin que se hubiera definido previamente
								if($explode1[11]!~/A-Z]/) {
									$desc = $explode1[1]."/".$explode1[9]."/".$explode1[17]."/".$explode1[7];
									$desc =~ s/\s+/_/g;
									$desc .= "\thttp://www.ncbi.nlm.nih.gov/pubmed/".$explode1[11];
								}
								else {
									$desc = $explode1[1]."/".$explode1[9]."/".$explode1[17]."/".$explode1[7];
									$desc =~ s/\s+/_/g;
									$desc .= "\t".$explode1[11];
								}
							}
						}
					}
					else {
						$entry=1;
						$count=$k;
						if($once==0) {
							$once=1;
						}
						# Aqui tambien se hacia uso de $explode, sin que se hubiera definido previamente
						if($explode1[11]!~/A-Z]/) {
							$desc = $explode1[1]."/".$explode1[9]."/".$explode1[17]."/".$explode1[7];
							$desc =~ s/\s+/_/g;
							$desc .= "\thttp://www.ncbi.nlm.nih.gov/pubmed/".$explode1[11];
						} else {
							$desc = $explode1[1]."/".$explode1[9]."/".$explode1[17]."/".$explode1[7];
							$desc =~ s/\s+/_/g;
							$desc .= "\t".$explode1[11];
						}
					}
				}
			}	
		}
	}
	#We split output to return the different colums
	#I don't like chorizo
	my @desc = split(/\s/, $desc);
	my $hgmd_desc .= $desc[0];
	my $hgmd_link .= $desc[1];
	return { 
		"HGMD_info" => $hgmd_desc,
		"Related_publication" => $hgmd_link
	};
}

sub get_bio_file_handle {
	my $chr = shift;
	my $bio_file_handle = new FileHandle;
	my $file = $ENV{'BIOBASE_PATH'};
	my $biobase = $file."/bioBaseSortUniques_".$chr.".txt";
	
	$bio_file_handle->open($biobase) or die("ERROR: Could not read from input file: $file\n");
	return $bio_file_handle;
}

sub complementario
# Function that returns the reverse-complementary of a given position nucleotide
{
	my $x=shift;
	my $y="";
	
	my @split=split("",$x);

	for(my $ii=0;$ii<scalar(@split);$ii++)
	{
		if($split[$ii] eq "A")
		{
			$y="T".$y;
		}	
		if($split[$ii] eq "T")
		{
			$y="A".$y;
		}
		if($split[$ii] eq "C")
		{
			$y="G".$y;
		}
		if($split[$ii] eq "G")
		{
			$y="C".$y;
		}
	}
	return($y);
}

1;