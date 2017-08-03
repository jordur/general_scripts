#!usr/bin/perl 
use strict;
use Bio::AlignIO;
use Bio::Align::AlignI;
use Bio::SimpleAlign ;
use Bio::SeqIO;

#  open COUNT, ">count_nucleotides.txt";
print "Input file name:\n";
my $inputfilename=<STDIN>;
	my $in  = Bio::SeqIO->new(-file => $inputfilename , '-format' => 'fasta');
		while ( my $seq = $in->next_seq() ) {
			my $len = $seq->length(),"\n";
			my $res=$seq->subseq(1,$len),"\n";
			my @res=split '',$res;
#			if (scalar@res < 60){
#				print ">",$res,"\n";
#				}
			print scalar@res,"\n";
#        print "Sequence ",$seq->id," first 10 bases ",$seq->subseq(1,10),"\n";
    }

#~ close COUNT;	


