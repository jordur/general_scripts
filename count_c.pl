#! usr/bin/perl -w
use strict;
use Bio::SeqIO;

# open FILE,"<translation.txt";
open OUT, ">gt_0.1_translation.txt";
my $in = Bio::SeqIO->new(-file => "<translation.txt"  , '-format' => 'Fasta') or die "Failed to open input file: $!";

while ( my $seq = $in->next_seq() ) {	
# 	print ">",$seq->id(),"\n";
	my $prot = $seq->seq();
	my@aa = split ('',$prot);
# 	print $aa[0],"\n";
# 	print scalar@aa,"\n";
	my $countC=0;
	foreach (@aa){
		if (/C/){
			++$countC;
# 			print $countC,"\t",scalar@aa,"\n";
			my $num = $countC/scalar@aa;
		if ($num gt 0.1){
			print OUT ">",$seq->id(),"\n",$prot,"\n";
			}
		}
	}



}
close OUT;


