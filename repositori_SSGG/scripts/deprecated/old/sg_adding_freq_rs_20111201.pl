#!/usr/bin/perl

use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Data::Dumper;

my ($input, $output) = @ARGV;
if (not defined $input or not defined $output) {
	die "sg_adding_freq_rs.pl <input> <output>\n";
}

print "Process started...".localtime()."\n";

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

my %changes;
my %Complementary;

# Describe complementary bases
$Complementary{'A'} = 'T';
$Complementary{'T'} = 'A';
$Complementary{'C'} = 'G';
$Complementary{'G'} = 'C';
$Complementary{'-'} = '-';

open (FH, "< $input") or die "Can´t open $input\n";
open (FREQ, "> $output\_with_freq.xls") or die "Can´t create $output\_with_freq.xls\n";

while (my $line = <FH>) {
	chomp $line;
	my @data = split (/\s+/, $line);
	my $snp = $data[8];

	if ($data[1] =~ m/chr/i) {
		$line .= "MAF_dbSNP\tPop";
		print FREQ $line,"\n";
		next;
	}
		
	if ($snp eq '-') {
		$line .= "-\t-";
		print FREQ $line,"\n";
		next;
	}
	
	my $var_adaptor = $registry->get_adaptor('human', 'variation', 'variation');
	my $var = $var_adaptor->fetch_by_name($snp);
	
	#print $snp,"\t",$var->ancestral_allele(),"\n";
	my $ancestral = ($var->ancestral_allele() || '-');
	
	my $alleles = $var->get_all_Alleles();
	
	my $freq;
	my $pop;
	
	if ($data[5] eq 'indel') {
		if (length($data[3]) > length($data[4])) {
			$data[4] = '-';
			my @bases = split ("", $data[3]);
			shift(@bases);
			$data[3] = join("", @bases);
		}
		elsif (length($data[4]) > length($data[3])) {
			$data[3] = '-';
			my @bases = split ("", $data[4]);
			shift(@bases);
			$data[4] = join("", @bases);
		}
		foreach my $allele (@{$alleles}) {
			my $allele_string = $allele->allele;
			my $string;
			if (!$allele->frequency) {
				my $popu = $allele->population;
				$string = $allele->allele."_NAF";
				if ($freq) {
					if ($freq !~ m/$string/i) {
						$freq .= "/".$string;
					}
				}
				else {
					$freq = $string;
				}
				if ($pop) {
					if ($popu && $pop !~ m/$popu->{name}/i) {
						$pop .= "/".($popu->{name} || '-');
					}					
				}
				else {
					$pop = ($popu->{name} || '-');
				}
			}
			elsif ($allele->frequency && $allele->frequency < 0.05) {
				my $popu = $allele->population;
				$string = $allele->allele."_".$allele->frequency;
				if ($freq) {
					if ($freq !~ m/$string/i) {
						$freq .= "/".$string;
					}
				}
				else {
					$freq = $string;
				}
				if ($pop) {
					if ($popu && $pop !~ m/$popu->{name}/i) {
						$pop .= "/".($popu->{name} || '-');
					}					
				}
				else {
					$pop = ($popu->{name} || '-');
				}
			}
		}
	}
	
	else {
		foreach my $allele (@{$alleles}) {
			my $allele_string = $allele->allele;
			my $string;
			
			#print Dumper $allele,"\n";
			if ($ancestral ne '-') {
				#allele = variant and ancestral allele = reference
				
				if (($allele_string eq $data[4] && $ancestral eq $data[3]) || ($allele_string eq $Complementary{$data[4]} && $ancestral eq $Complementary{$data[3]})) {

					if (!$allele->frequency) {
						my $popu = $allele->population;
						$string = $allele->allele."_NAF";
						if ($freq) {
							if ($freq !~ m/$string/i) {
								$freq .= "/".$string;
							}
						}
						else {
							$freq = $string;
						}
						if ($pop) {
							if ($popu && $pop !~ m/$popu->{name}/i) {
								$pop .= "/".($popu->{name} || '-');
							}					
						}
						else {
							$pop = ($popu->{name} || '-');
						}
					}
					elsif ($allele->frequency && $allele->frequency < 0.05) {
						my $popu = $allele->population;
						$string = $allele->allele."_".$allele->frequency;
						if ($freq) {
							if ($freq !~ m/$string/i) {
								$freq .= "/".$string;
							}
						}
						else {
							$freq = $string;
						}
						if ($pop) {
							if ($popu && $pop !~ m/$popu->{name}/i) {
								$pop .= "/".($popu->{name} || '-');
							}					
						}
						else {
							$pop = ($popu->{name} || '-');
						}
					}
				}

				#allele = reference and ancestral allele = variant
				
				elsif (($allele_string eq $data[3] && $ancestral eq $data[4]) || ($allele_string eq $Complementary{$data[3]} && $ancestral eq $Complementary{$data[4]})) {			
					
					if (!$allele->frequency) {
						my $popu = $allele->population;
						$string = $allele->allele."_NAF";
						if ($freq) {
							if ($freq !~ m/$string/i) {
								$freq .= "/".$string;
							}
						}
						else {
							$freq = $string;
						}
						if ($pop) {
							if ($popu && $pop !~ m/$popu->{name}/i) {
								$pop .= "/".($popu->{name} || '-');
							}					
						}
						else {
							$pop = ($popu->{name} || '-');
						}
					}
					elsif ($allele->frequency && $allele->frequency < 0.05) {
						my $popu = $allele->population;
						$string = $allele->allele."_".$allele->frequency;
						if ($freq) {
							if ($freq !~ m/$string/i) {
								$freq .= "/".$string;
							}
						}
						else {
							$freq = $string;
						}
						if ($pop) {
							if ($popu && $pop !~ m/$popu->{name}/i) {
								$pop .= "/".($popu->{name} || '-');
							}					
						}
						else {
							$pop = ($popu->{name} || '-');
						}
					}
				}
				elsif ($allele_string eq $ancestral && $data[7] !~ m/hetero/i) {
					if (!$allele->frequency) {
						my $popu = $allele->population;
						$string = $allele->allele."_NAF";
						if ($freq) {
							if ($freq !~ m/$string/i) {
								$freq .= "/".$string;
							}
						}
						else {
							$freq = $string;
						}
						if ($pop) {
							if ($popu && $pop !~ m/$popu->{name}/i) {
								$pop .= "/".($popu->{name} || '-');
							}					
						}
						else {
							$pop = ($popu->{name} || '-');
						}
					}
					elsif ($allele->frequency && $allele->frequency < 0.05) {
						my $popu = $allele->population;
						$string = $allele->allele."_".$allele->frequency;
						if ($freq) {
							if ($freq !~ m/$string/i) {
								$freq .= "/".$string;
							}
						}
						else {
							$freq = $string;
						}
						if ($pop) {
							if ($popu && $pop !~ m/$popu->{name}/i) {
								$pop .= "/".($popu->{name} || '-');
							}					
						}
						else {
							$pop = ($popu->{name} || '-');
						}
					}
				}
			}
			else {
				if (!$allele->frequency) {
					my $popu = $allele->population;
					$string = $allele->allele."_NAF";
					if ($freq) {
						if ($freq !~ m/$string/i) {
							$freq .= "/".$string;
						}
					}
					else {
						$freq = $string;
					}
					if ($pop) {
						if ($popu && $pop !~ m/$popu->{name}/i) {
							$pop .= "/".($popu->{name} || '-');
						}					
					}
					else {
						$pop = ($popu->{name} || '-');
					}
				}
				elsif ($allele->frequency && $allele->frequency < 0.05) {
					my $popu = $allele->population;
					$string = $allele->allele."_".$allele->frequency;
					if ($freq) {
						if ($freq !~ m/$string/i) {
							$freq .= "/".$string;
						}
					}
					else {
						$freq = $string;
					}
					if ($pop) {
						if ($popu && $pop !~ m/$popu->{name}/i) {
							$pop .= "/".($popu->{name} || '-');
						}					
					}
					else {
						$pop = ($popu->{name} || '-');
					}
				}
			}
		}
	}
	
	$line .= "\t".($freq || '-')."\t".($pop || '-');
	
	print FREQ $line,"\n";

}
print "Process finished... ".localtime()."\n";

close FH;
close FREQ;
exit;
