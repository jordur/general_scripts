#!/usr/bin/perl -w
use Bio::SearchIO;

my $in = new Bio::SearchIO(-format => 'blast',-file   => @ARGV) ;
print "query_name\tquery_length\thit description\thit length\t#hsps\tstrand hit\tstart hit\tstart query\tend hit\tend query\tnumero nt iguales\thsp length\tbit score\tE-valor\tFracc_conserved\taccession\n";
while( my $result = $in->next_result ) 
{
	while( my $hit = $result->next_hit )
	{
    			while( my $hsp = $hit->next_hsp ) 
			{
          			print $result->query_name,"\t",
                       		$result->query_length,"\t",
				$hit->description,"\t",
				$hit->length,"\t",                        
				$hit->num_hsps,"\t", 
				$hsp->strand('hit'),"\t",
				#$hsp->strand('query'),"\t",                  
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


