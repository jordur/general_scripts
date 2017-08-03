#!/usr/bin/perl

sub usage
{
    print "\nCOMO SE USA:\n";
    print " Input: 1) fichero de blastclust  2) fichero con secuencias en fasta en base-space\n";
    print " Output: fichero con el conteo de los miRNAs. Una linea por cada miRNA con coverage > 10x  \n";
    print "\n";
    exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}

open(BLASTCLUST,"<",@ARGV[0]);
open(FASTA,"<",@ARGV[1]);
open (CONTEO, ">", 'conteo.tab') or die "No puedo abrir conteo.tab.\n";

my @blastclust;
my @split_blastclust;
my $lineabc;
my $lineafasta;
my @split_por_espacio;
$contado=0;

while  (my $lineafasta=<FASTA>)
{
        chomp($lineafasta);
        push (@fasta,$lineafasta);
}

while  (my $lineabc=<BLASTCLUST>)
{
        chomp($lineabc);
	#@split_por_espacio = split (/\s+/, $lineabc);
	#print $split_por_espacio[0],"\t";
	for (my $i=0; $i<length ($lineabc); $i++)
	{
		$caracter = substr($lineabc,$i,2);
		if ($caracter eq "F3")
		{
			$count++;
		}
	}
	if ($count >1)
	{
		@split_por_espacio = split (/\s+/, $lineabc);
		#print $split_por_espacio[0],"\t";
		#print $count,"\t";
		for (my $a=0; $a<=$#fasta;$a=$a+2)
        	{
                	$cabecera = $fasta[$a];
			$cabecera =~ s/>//;
			$cabecera1 = $cabecera;
                	$cabecera1 =~ s/,.*$//;
			if ($cabecera1 eq $split_por_espacio[0])
			{
				#print $split_por_espacio[0],"\t";
				print $cabecera,"\t";
				print $count,"\t";
				print $fasta[$a+1],"\n";
			}
			#else
			#{
				#splice (@fasta,$a);
				#splice (@fasta,$a+1);
			#}
        	}
		#print $split_por_espacio[0],"\t";
		#print $count,"\t";
	}
	$count = 0;
}


close (BLASTCLUST);
close (FASTA);

#for ($i =0;$i<=$#doscolumnas ; $i++)
#{
#	@split_dos = split (/\t/, $doscolumnas[$i]);
#	print FASTA ">".$split_dos[0],"\n";
#	print FASTA $split_dos[1],"\n";
#}

close (CONTEO);
