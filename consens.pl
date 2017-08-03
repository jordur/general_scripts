#!usr/bin/perl -w
use strict;
use Bio::Align::AlignI;
use Bio::AlignIO;

my @files= qw/reprolysin_cobalt/;
foreach (@files){
###########Parseo MSA i agafa consens d'una posició representada com a mínim n cops ##########
#~ open OUTFILE, ">>/home/jordi/Feina/prova/rm/class_script/final/toxines_senseNA/consensus/$_.cons.fasta";
# print "Enter protein name:\n";
my $prot=$_;
# print $_;
# print "Enter snake name\n";
# my $snake=<STDIN>;
my $str = Bio::AlignIO->new('-file' => $prot);
my $aln = $str->next_aln();
my $depth = $aln->no_sequences;
#~ print $depth;
my $result= $aln->consensus_string(400/$depth),"\n";
$result =~ s/\?/X/g;
print ">\n",$result,"\n";
#~ close OUTFILE;

}


