#!usr/bin/perl -w
use Bio::SearchIO::blasttable;
use Bio::Search::Result::GenericResult;
use Bio::Search::Hit::GenericHit;

open FILE, "<est_results.txt";
open OUT, ">est_best_hits.txt";
my @est=<FILE>;
my $file="est_results.txt";
#~ print @est,"\n";



my $parser = Bio::SearchIO->new(-file   => $file,
                                 -format => 'blasttable');
				 #~ print $parser;

foreach my $blast(@est){
	while ($blast =~ /# BLASTN/){
	my $result = $parser->next_result;
#~ print $result;
	my $id = $result->query_name();
	my @hits = $result->hits;
	#~ print $hits[0],"\n";
	foreach my $hit ($hits[0]){
		my $hit_names=$hit->name();
		print OUT $id,"\t",$hit_names,"\n";
	}
			
			}
		
	}	
close FILE;
close OUT;