#! /usr/bin/perl -w

# -----------------------
# -- Include libraries --
# -----------------------

use strict;
use Bio::SeqIO;
use Bio::SeqIO::qual;
use Bio::Seq::PrimaryQual;
use Cwd;
use Getopt::Long;

# --------------------------------- 
# -- global & default parameters --
# ---------------------------------

my $lineA = $ARGV[0];
my $lineB = $ARGV[1];
my $sample = $ARGV[2];
my $name = $ARGV[3];
my $seqname = $ARGV[4]; 
my $lineC = $ARGV[5];
my $path = getcwd();
my $batch = $ARGV[6];

open OUT1 , ">".$path."/".$batch."/".$sample."/results/assembly/".$sample."-rename-contigs.fa";
open OUT2 , ">".$path."/".$batch."/".$sample."/results/assembly/".$sample."-rename-singlets.fa";
open OUT3 , ">".$path."/".$batch."/".$sample."/results/assembly/".$sample."-rename-contigs.qual";

# --------------------------------
# -- Input parameters & options --
# --------------------------------

if (@ARGV != 7){
    die "
        Check your options!
\n".localtime()."\n=================================================================\n\n";

}

# ---------------
# -- functions --
# --------------- 


#### Parse contigs CAP3 file ###
my $input  = Bio::SeqIO->new(
				 -format => 'fasta',
				 -file   => $lineA,
				 );
    
my $count=0;
while ( my $seq = $input->next_seq() ) {
    $count++;
    my $id = $seq->id;
    if ($id =~ /(\w)+/){
	my $sequence = $seq->seq;
    print OUT1 ">ci|0000000".$count."|vnm|".$sample."|".$seqname." ".$name."\n".$sequence."\n";
    }
    
}
    

##### Parse CAP3 singlets file ###
my $input2  = Bio::SeqIO->new(
			      -format => 'fasta',
			      -file   => $lineB,
			      );
$count=0;
while ( my $seq = $input2->next_seq() ) {
   my $id = $seq->id;
   if ($id =~ /(\w)+/){
	my $sequence = $seq->seq;
	$count++;
	print OUT2  ">si|0000000".$count."|vnm|".$sample."|".$seqname." ".$name."\n".$sequence."\n";
   }
   
}

#### Parse contigs qual CAP3 file ###   
my $input3  = Bio::SeqIO->new(
			      -format => 'qual',
			      -file   => $lineC,
			      );
$count=0;
while ( my $seq = $input3->next_seq()) {
   $count++;
   my $id = $seq->id;
   my @quals = @{$seq->qual()};
   if ($id =~ /(\w)+/){
	print OUT3 ">ci|0000000".$count."|vnm|".$sample."|".$seqname." ".$name."\n"."@quals\n";
   }
   
}
