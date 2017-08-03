#! /usr/bin/perl -w
use strict;
use Bio::TreeIO;
use strict;
my ($filein,$fileout) = @ARGV;
my ($format,$oformat) = qw(newick nexus);
my $in = Bio::TreeIO->new(-file => $filein, -format => $format);
my $out= Bio::TreeIO->new(-format => $oformat, -file => ">$fileout");

while( my $t = $in->next_tree ) {
  $out->write_tree($t);
}
