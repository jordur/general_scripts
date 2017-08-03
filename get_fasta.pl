#!usr/bin/env perl -w
use strict;
use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Align::AlignI;
use Bio::SimpleAlign ;

#### Parsing file #####
#open OUT, ">repro_UTR_target_scan_format"; 
#open IN, $ARGV[0];
#print "Enter input fasta file name:\n";
my $in1 = $ARGV[0];
my $in2;
if ($in1 =~/.\/(.*).collapsed.orf/){
    $in2 = $1;
}

my $id2;
#my $in1 = <IN>;
#print "Quina serp és?:";
#print "Quina toxina és?:";
#my $serp =<STDIN>;
#chomp($serp);
#open OUT, ">".$serp."_UTR_target_scan_format"; 
#print "Adult(1) o Neonato(2)?:";
#my $dev = <STDIN>;
#chomp($dev);

#print "@identificadors\n";
#    chomp($identificador);
#    print $identificador."\n";
my $in = Bio::SeqIO->new( -file => $in1, '-format' => 'Fasta')  or die "Failed to open input file: $!";
my $count=1;
while ( my $seq = $in->next_seq() ) {
#    foreach my $identificador (@identificadors){
#	chomp(@identificadors);
#        print $identificador."\n";
    my $dna = $seq->seq();
    my $id = $seq ->id();
#    print $id,"\n";
    if ($id =~ /lcl\|(.+)/){
	$id2 = $1;
	#print $id2,"\n";
	if($id2 =~ /(Locus.+)\/\d+/){
	    #print $id2,"\n";
	    $id2 = $1;
	    #print $id2,"\n";
	}
    }
#	print $in2,"\n";
 	open OUT,">".$in2."_".$id2."_".$count.".fasta";
#	if ( $id !~ $identificador){
#    	if ($id ~ /@identificador/){    
	print OUT ">".$id2."\n".$dna."\n";
#	    print $identificador." found!\n";
#	}else{
#	    print OUT ">".${serp}."_".$id."\n".$dna."\n";
#	    print OUT $id."\t".$dna."\t".$dev."\n";
#	    print OUT $id."\t".$dev."\t".$dna."\n";
	#print OUT ">".$id."\n".$dna."\n";
#	}
	#    }
	$count+=1;
	close OUT;
}
#close IN;
