#! /usr/bin/perl -w
use strict;
use Bio::DB::GenBank;
use Bio::AnnotatableI;

open INPUT, "< uniq_gi";
# my $line = <STDIN>;
my $gb = Bio::DB::GenBank->new();
#         print "Waiting for reference id...\n";
#         print "Enter file name:\n";
	my @file=<INPUT>;
foreach my $gi(@file){

#       print "Fetching gi ", $gi, " ... ";
        my $seq = $gb->get_Seq_by_gi($gi); # GI Number
#       print "done. ";

        my @features = $seq->get_SeqFeatures(); # just top level
if(my $feat->has_tag('protein_id') { 
	foreach $feat ( @features ) {
                if( $feat->primary_tag eq 'CDS' ) {
                        my ($prot_id) = $feat->get_tag_values('protein_id');

#                       print "Fetching prot id ", $prot_id, " ... ";
                        my $prot_obj = $gb->get_Seq_by_acc( $prot_id );
#                       print "done. ";

#                       print "Writting fasta ... ";

#                         $out->write_seq( $prot_obj );

#                       print "done\n";

                        my $secuencia=$prot_obj->seq();
                        open REF,">>ref.fasta";
#                         chomp(my $ref=<REF>);
                        print  REF ">",$prot_id,"\n";
                        print REF "$secuencia\n";


                }
	}
   }
}
close REF;
close INPUT;
