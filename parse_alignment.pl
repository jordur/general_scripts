#!usr/bin/perl -w
use strict;
use Bio::AlignIO;
use Bio::Align::AlignI;
use Bio::SimpleAlign ;


my $inputfilename = "l-amino_mafft";
    my $in  = Bio::AlignIO->new(-file => $inputfilename , '-format' => 'clustalw');
#     my $out = Bio::AlignIO->new(-file => ">out.l-amino" , '-format' => 'fasta');


#     while ( my $aln = $in->next_aln() ) {
#         $out->write_aln($aln);
#     }

if( my $aln = $in->next_aln ){
# 	$aln->uniq_seq();
# 	print $aln,"\n";
# 	
	
#  for my $s ($aln->each_seq ) {
#    if( $s->id eq $ref_seq ) {
# 	my $threshold=$aln->no_sequences;
# 	my $seq= $aln->consensus_string(200/$threshold)."\n";
# 

# 	foreach  my $seq ( $aln->each_seq() ) {
# 		print $seq,"\n";
		
		
		
			
# 	}
	


# 	my $pos_start=1;
# 	my $pos_end=516;
#      my $col_start = $aln->column_from_residue_number($s->id, $pos_start);
#      my $col_end = $aln->column_from_residue_number($s->id, $pos_end);
#      print "grabbing columns $col_start .. $col_end\n";
#      my $piece = $aln->slice($col_start, $col_end);
#      $out->write_aln($piece);
#      last; # all done
#    		}
#  	}
# }

 foreach my $seq ($aln->each_seq) {
      my $res = $seq->subseq($pos, $pos);
      $count{$res}++;
  }
  foreach $res (keys %count) {
      printf "Res: %s  Count: %2d\n", $res, $count{$res};
	  }
}


