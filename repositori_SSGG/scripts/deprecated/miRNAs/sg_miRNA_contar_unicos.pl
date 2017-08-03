#!/usr/bin/perl
sub usage {
    print "\nCOMO SE USA:\n";
    print "EL ORDEN ES MUY IMPORTANTE \n\t ";
    print "Primero, fichero con 2 columnas, 1 con nombre secuencia, 2 secuencia  \n\t ";
    print "Segundo, fichero con 1 columna con las secuencias unicas en base-space \n";
    print "\n";
    exit(1);
}
if(scalar(@ARGV) == 0){
    usage();
}

open(FASTA,"<",@ARGV[0]);
open(SECUENCIA,"<",@ARGV[1]);
#open(FASTAMIRNA,">", 'miRNA.fasta') or die "miRNA.fasta.\n";
my @fasta;
my @secuencia;
#my @fasta_linea;
#my @sequencia_linea;
while  (my $fasta_linea=<FASTA>)
{
  	chomp($fasta_linea);
	#$fasta_linea =~ tr/agtcn/AGTCN/; 
  	push (@fasta,$fasta_linea);
}

close (FASTA);


while ( my $secuencia_linea=<SECUENCIA>)
{
  chomp($secuencia_linea);
  #$secuencia_linea=~ tr/agtcn/AGTCN/;
  push (@secuencia,$secuencia_linea);
}

close (SECUENCIA);

$contador=0;
for (my $i=0; $i<=$#secuencia;$i++)
{
	for (my $a=0; $a<=$#fasta;$a++)
	{	
		@valor = split (/\t/,$fasta[$a]);	
		
		if ($valor[1] eq $secuencia[$i])
		{
			$contador=$contador+1;
			#splice (@valor,$a,1);
			#print @valor,"\n";
		}
		
	}
	print  ">miRNA",$i,"_x",$contador,"\n";
	print $secuencia[$i],"\n";
	$contador=0;
	#print @valor,"\n";
}

#close (FASTAMIRNA);
