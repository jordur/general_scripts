#!/usr/bin/perl

sub usage
{
    print "\nPARA QUE SIRVE: compara fichero del conteo de miRBase con el fichero de miRNAs maduros\n";
    print " Input: 1) fichero conteo de miRBase 2) fichero con miRNAs maduros ej: mir2<TAB>chr<TAB>inicio<TAB>final<TAB>orientacion\n";
    print " Output: fichero tabulado para utilizar con cal_miRNA_extract_seq.pl mir2<TAB>inicio<TAB>final<TAB>orientacion  \n\n ";
    exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}

open(CONTEO,"<",@ARGV[0]);
open(MIRBASE,"<",@ARGV[1]);

my @conteo;
my @mirbase2;

while  (my $linea2=<CONTEO>)
{
        chomp($linea2);
        push (@conteo,$linea2);
}

close (CONTEO);


while  (my $linea=<MIRBASE>)
{
        chomp($linea);
        push (@mirbase2, $linea);
}

close (MIRBASE);

my $contador = 0;
for (my $a=0; $a<=$#mirbase2 ;$a++)
{
	@mirbase= split ("\t",$mirbase2[$a]);	
	print $mirbase[0],"\t",$mirbase[2],"\t",$mirbase[3],"\t",$mirbase[4],"\t";
	for (my $i=0;$i<=$#conteo;$i++)
	{
		@conteo_split= split ("\t",$conteo[$i]);
		@conteo_split2 = split (/\"/,$conteo[$i]);
		if (($mirbase[0] eq $conteo_split2[3]) && ((abs($conteo_split[3]-$mirbase[2])) <= 10))
		{
			$contador = $contador+$conteo_split[7];
		}
	}
	print $contador,"\n";
	$contador = 0;	
}
