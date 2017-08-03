#!usr/bin/perl -w
use strict;
use Bio::Seq;
use Bio::SeqIO;

open OUT , ">outfile.fasta";

sub codon2aa{
my($codon)=@_;
$codon= uc $codon;
my(%genetic_code) = (
 	'TCA'=>'S', #SERINE
	'TCC'=>'S', #SERINE
	'TCG'=>'S',  #SERINE
	'TCT'=>'S', #SERINE
	'TTC'=>'F', #PHENYLALANINE
	'TTT'=>'F', #PHENYLALANINE
	'TTA'=>'L', #LEUCINE
	'TTG'=>'L', #LEUCINE
	'TAC'=>'Y', #TYROSINE
	'TAT'=>'Y', #TYROSINE
	'TAA'=>'*', #STOP
	'TAG'=>'*', #STOP
	'TGC'=>'C', #CYSTEINE
	'TGT'=>'C', #CYSTEINE
	'TGA'=>'*', #STOP
	'TGG'=>'W', #TRYPTOPHAN
	'CTA'=>'L', #LEUCINE
	'CTC'=>'L', #LEUCINE
	'CTG'=>'L', #LEUCINE
	'CTT'=>'L', #LEUCINE
	'CCA'=>'P', #PROLINE
	'CAT'=>'H', #HISTIDINE
	'CAA'=>'Q', #GLUTAMINE
	'CAG'=>'Q', #GLUTAMINE
	'CGA'=>'R', #ARGININE
	'CGC'=>'R', #ARGININE
	'CGG'=>'R', #ARGININE
	'CGT'=>'R', #ARGININE
	'ATA'=>'I', #ISOLEUCINE
	'ATC'=>'I', #ISOLEUCINE
	'ATT'=>'I', #ISOLEUCINE
	'ATG'=>'M', #METHIONINE
	'ACA'=>'T', #THREONINE
	'ACC'=>'T', #THREONINE
	'ACG'=>'T', #THREONINE
	'ACT'=>'T', #THREONINE
	'AAC'=>'N', #ASPARAGINE
	'AAT'=>'N', #ASPARAGINE
	'AAA'=>'K', #LYSINE
	'AAG'=>'K', #LYSINE
	'AGC'=>'S', #SERINE
	'AGT'=>'S', #SERINE
	'AGA'=>'R', #ARGININE
	'AGG'=>'R', #ARGININE
	'CCC'=>'P', #PROLINE
	'CCG'=>'P', #PROLINE
	'CCT'=>'P', #PROLINE
	'CAC'=>'H', #HISTIDINE
	'GTA'=>'V', #VALINE
	'GTC'=>'V', #VALINE
	'GTG'=>'V', #VALINE
	'GTT'=>'V', #VALINE
	'GCA'=>'A', #ALANINE
	'GCC'=>'A', #ALANINE
	'GCG'=>'A', #ALANINE
	'GCT'=>'A', #ALANINE
	'GAC'=>'D', #ASPARTIC ACID
	'GAT'=>'D', #ASPARTIC ACID
	'GAA'=>'E', #GLUTAMIC ACID
	'GAG'=>'E', #GLUTAMIC ACID
	'GGA'=>'G', #GLYCINE
	'GGC'=>'G', #GLYCINE
	'GGG'=>'G', #GLYCINE
	'GGT'=>'G', #GLYCINE
);


if(exists $genetic_code{$codon}){
    return $genetic_code{$codon};
    }else{
    print "";
     }
}

print "Enter your fasta file name:\n";
my $file = <STDIN>;

my $in  = Bio::SeqIO->new(-file => $file , '-format' => 'Fasta')  or die "Failed to open input file: $!";

while ( my $seq = $in->next_seq() ) {
# 	print OUT ">",$seq->id(),"\n";
    	my $dna = $seq->seq();
# print "Enter a DNA sequence:\n";
# chomp(my $dna=<STDIN>);

	my $protein1;
	my $protein2;
	my $protein3;
	my $codon;
	
		for(my $i=0; $i<(length($dna)-2); $i+=3){
			$codon = substr($dna,$i,3);
			$protein1 .= codon2aa($codon);
		}
		      if ($protein1 =~ s/[0-9]//g){
# 		      print $protein1,"\n";

		      print OUT ">",$seq->id(),".1\n",$protein1,"\n";
# 		      print OUT "\n";
		      }else{
		      print OUT ">",$seq->id(),".1\n",$protein1,"\n";
# 		      print OUT "\n";
		      }
	substr($dna,0,1)="";
# 
# 
		for(my $i=0; $i<(length($dna)-2); $i+=3){
			$codon = substr($dna,$i,3);
			$protein2 .= codon2aa($codon);
		}
		      if ($protein2 =~ s/[0-9]//g){

		      print OUT ">",$seq->id(),".2\n",$protein2,"\n";
# 		      print OUT "\n";
		      }else{
		      print OUT ">",$seq->id(),".2\n",$protein2,"\n";
# 		      print OUT "\n";
		      }	

	substr($dna,0,1)="";

		for(my $i=0; $i<(length($dna)-2); $i+=3){
			$codon = substr($dna,$i,3);
			$protein3 .= codon2aa($codon);
		}
		      if ($protein3 =~ s/[0-9]//g){

		      print OUT ">",$seq->id(),".3\n",$protein3,"\n";
# 		      print OUT "\n";
		      }else{
		      print OUT ">",$seq->id(),".3\n",$protein3,"\n";
# 		      print OUT "\n";
		      }

  }
close OUT;