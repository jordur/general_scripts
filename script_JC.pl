#! /usr/bin/perl -w
use strict;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::Tools::Run::Alignment::Clustalw;

print "Enter sequence name:\n";
my $input = <STDIN>;
chomp($input);
my $in  = Bio::AlignIO->new(-file => $input,
                         -format => 'clustalw');
while ( my $aln = $in->next_aln ) { 
#   my  $out->write_aln($aln); 
#   print $out,"\n";
   my $percent_ident = $aln->percentage_identity;
   print $percent_ident,"\n";

}

