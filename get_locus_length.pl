#!/usr/bin/perl -w

use strict;
use Bio::DB::GenBank;
use Bio::SeqIO;


# my $out = Bio::SeqIO->new ( "-file"   => ">prot.fa",
# 			    "-format" => "Fasta" );

open OUT, ">lengh.txt";
my $gb = Bio::DB::GenBank->new();
my @gi = <>;
foreach my $a (@gi) {
		
# 	my $gi = <STDIN>;
	chomp($a);
	
	print "Fetching gi '", $a, "' ... ";
	my $seq = $gb->get_Seq_by_gi($a); # GI Number
	print "done\n. ";
# 	print $seq,"\n";
	my @features = $seq->get_SeqFeatures(); # just top level
  foreach my $feat ( @features ) {
		if( $feat->primary_tag eq 'source' ) {
# 			print $feat,"\n";
			 my $len = $feat->length;
			print OUT $len,"\n";
# 			my ($prot_id) = $feat->get_tag_values('protein_id');
# 			
# 			print "Fetching prot id ", $prot_id, " ... ";
# 			my $prot_obj = $gb->get_Seq_by_acc( $prot_id );
# 			print "done. ";
# 			
# 			print "Writting fasta ... ";
# 			$out->write_seq( $prot_obj );
# 			
# 			print "done\n";
		}
  }
}
close OUT;
__END__

