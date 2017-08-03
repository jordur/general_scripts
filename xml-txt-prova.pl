#!/usr/bin/perl -w
use strict;
use Bio::Perl;
use Bio::SearchIO;
use Bio::SeqFeature::SimilarityPair;

my $blast_report = "l-amino";

# Get the report
my $searchio = new Bio::SearchIO (-format => 'blastxml',
				  -file   => $blast_report);



while( my $result = $searchio->next_result ) {
    
    while(my $hit = $result->next_hit){
	my  $hit_name = $hit->name ;
      	my $strand = $hit->strand('hit');

	$strand =~ s/1//;

	my %mis_hsps;
        my %mis_from;
        my %mis_to;
	
	my $count = 0;

	#if ($hit_name eq "uaccno=FGSMDPN09FM4Z2") {
	    while(my $hsp = $hit->next_hsp){
		$count++;
		my $key = $count ."|". $strand ."". ($hsp->frame + 1);
                my $from = $count;
		my $feat  = new Bio::SeqFeature::FeaturePair(-feature1 =>  $hsp->query,
							     -feature2 => $hsp->hit);

                my $start = $hsp->start('hit');
                my $end = $hsp->end('hit');
		print $start . " ".  $end . "\n";
                $mis_from{$start} = $end;
                
                #$mis_from{$count} = $start;
                #$mis_to{$count} = $end;
		#my $evalue = $hsp->evalue;
		#$mis_hsps{$key} = $hsp->evalue;
                #print $hit_name ." ". $key . " ". $evalue . "\n";
		}
	    
            #my $best_from = pop @mis_from;
            my @sorted_start = sort{$mis_from{$b} <=> $mis_from{$a} }keys%mis_from;
            my @sorted_end = sort{$mis_from{$b} <=> $mis_from{$a} }values%mis_from;
	    #my @sorted_hsps = sort{$mis_hsps{$b} <=> $mis_hsps{$a} }keys %mis_hsps;
	    #my @sorted_from = sort{ values %mis_from;
            #my @sorted_to = sort{$mis_to{$b} <=> $mis_to{$a} }values %mis_to;
            #print @sorted_from;
	    my $first_start = pop @sorted_start;
            my $first_end = shift @sorted_end;
	    #my $best_from = pop @sorted_from;
            #my $best_to = pop @sorted_to;
            #print $best_from;
	    #$best_frame =~ s/[\d+]\|//;
	    #print $hit_name ." ". $first_start ." ". $first_end ."\n";
	}
    }
#}



#my $feat  = new Bio::SeqFeature::FeaturePair(-feature1 =>  $hsp->query, -feature2 => $hsp->hit,);
