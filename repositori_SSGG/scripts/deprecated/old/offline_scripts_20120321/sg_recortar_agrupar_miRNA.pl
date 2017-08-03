#!/usr/bin/perl

sub usage {
    print "\nCOMO SE USA:\n";
    print "EL ORDEN ES MUY IMPORTANTE \n\t ";
    print "PRIMERO FICHERO .ma  \n\t ";
    print "SEGUNDO FICHERO .gff con las secuencias en base-space> \n";
    exit(1);
}
if(scalar(@ARGV) == 0){
    usage();
}

open(MA,"<",@ARGV[0]);
open(GFF,"<",@ARGV[1]);
open(FASTAMIRNA,">", 'miRNA.fasta') or die "miRNA.fasta.\n";
my @ma;
my @gff;
my @ma_linea;
my @gff_linea;
while  (my $linea=<MA>)
{
  	chomp($linea);
  	push (@ma,$linea);
}

for (my $a=0; $a<=$#ma; $a=$a+2)
{
	@longitud = split (/\./,$ma[$a]);
	push (@ma_linea,$longitud[2]);
	$trans = $longitud[0];
	$trans =~s />//;
	$trans =~s /,.*$//;
	push (@nombre,$trans);
}

close (MA);


while ( my $linea2=<GFF>)
{
  chomp($linea2);
  push (@gff,$linea2);
}

close (GFF);


for (my $i=0; $i<=$#gff;$i++)
{

	@gff_linea = split (/\t/,$gff[$i]);
	for (my $u=0;$u<=$#nombre;$u++)
	{
		if ($gff_linea[0] eq  $nombre[$u])
		{ 
			$gff_linea[8]=~s/b=//;
        		$gff_linea[8]=~s/;.*$//;
			$miRNA= substr ($gff_linea[8],0,$ma_linea[$u]);
			print FASTAMIRNA $nombre[$u],"\t",$miRNA,"\n";
			 	
		}
	}

}

close (FASTAMIRNA);
