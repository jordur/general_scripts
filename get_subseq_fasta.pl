#!/usr/bin/perl -w
use strict;

use Bio::SeqIO; 
use Bio::Seq;
use Bio::DB::GenBank;
use Bio::PrimarySeqI;
use File::Basename;
use Getopt::Long;
use Cwd;

# ---------------------------------
# -- global & default parameters --
# ---------------------------------
my $input;
my $start;
my $end;
my $output;

# --------------------------------
# -- Input parameters & options --
# --------------------------------
(my $basename, my $dir, my $ext) = fileparse($0, qr/\.[^.]*/);

GetOptions (
                                "i=s"   => \$input,
                                "start=s"   => \$start,
                                "end=s"	=> \$end,
                                "output=s" => \$output,
);

print "\n=================================================================
$basename: Obtaining a subset of a fasta sequence\n";

if (not defined $input) {
        die "
                Options:
                -i input fasta file
                -start Start position
                -end End position
                -ouput Output file 
\n".localtime()."\n=================================================================\n\n";
}

print  "\nProcess started at... ".localtime()."\n";
&main ();
print  "\nProcess finished... ".localtime()."\n";


# ---------------
# -- functions --
# ---------------
sub main (){
	open OUT, ">".$output;
	my $in = Bio::SeqIO->new (-file =>  $input, -format => "Fasta");
	my $id;
	my $desc;
	my $subseq;
	while (my $seq = $in->next_seq) { 
		$id = $seq->primary_id;
		$desc = $seq -> description;
		#$subseq = $in->subseq(10,40);
		#chomp($id);
		#chomp($desc);
		print OUT ">",$id,"\t",$desc,"\n",$seq->subseq($start,$end),"\n";
		print ">",$id,"\t",$desc,"\n",$seq->subseq($start,$end),"\n";
	}
}

