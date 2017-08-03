#!/usr/bin/perl -w
use Bio::SearchIO;

my $in = new Bio::SearchIO(-format => 'blast',-file   => @ARGV) ;
print "Query_name\tQuery_length\tHit_description\tHit_length\t#hsps\tHit_frame\tHit_start\tQuery_start\tHit_end\tQuery_end\tIdentical_residues\thsp_length\tBit_score\tE-value\t","%","_conserved\tAccession_ID\n";
while( my $result = $in->next_result ) 
{
	while( my $hit = $result->next_hit )
	{
    			while( my $hsp = $hit->next_hsp ) 
			{
          			my $blastframe = ($hsp->query->frame + 1) * $hsp->query->strand;
				print $result->query_name,"\t",
				$result->query_length,"\t",
				$hit->description,"\t",
				$hit->length,"\t",                        
				$hit->num_hsps,"\t",
				$blastframe,"\t",
				$hsp->start('hit'),"\t",
				$hsp->start('query'),"\t",
                        	$hsp->end('hit'),"\t",
				$hsp->end('query'),"\t",
				$hsp->num_conserved,"\t",
				$hsp->length,"\t",
				$hit->bits,"\t",
				$hsp->evalue,"\t",
				($hsp->frac_conserved)*100,"\t",  
                        	$hit->accession,"\n",
      			}
    	}
}


