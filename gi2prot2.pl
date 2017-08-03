#!/usr/bin/perl -w

use strict;
use Bio::DB::GenBank;
use Bio::SeqIO;


my $out = Bio::SeqIO->new ( "-file"   => ">ref.fasta",
			    "-format" => "Fasta" );
my $a = "Waiting for reference sequence...\n";
	print $a;
my $gb = Bio::DB::GenBank->new();

	my $gi = <STDIN>;
	chomp($gi);
	
	print "Fetching gi ", $gi, " ... ";
	my $seq = $gb->get_Seq_by_gi($gi); # GI Number
	print "done. ";
	
	my @features = $seq->get_SeqFeatures(); # just top level
  foreach my $feat ( @features ) {
		if( $feat->primary_tag eq 'CDS' ) {
			my ($prot_id) = $feat->get_tag_values('protein_id');
			
			print "Fetching prot id ", $prot_id, " ... ";
			my $prot_obj = $gb->get_Seq_by_acc( $prot_id );
			print "done. ";
			
			print "Writting fasta ... ";
			$out->write_seq( $prot_obj );
			
			print "done\n";
		}
  }


__END__

=head1 NAME - myGiTest

=head1 AUTHOR

Javier Santoyo B<email> jsantoyo@cipf.es

=head1 DATE

July  9, 2009

