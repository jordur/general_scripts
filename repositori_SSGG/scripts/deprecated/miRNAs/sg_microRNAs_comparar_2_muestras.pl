#!/usr/bin/perl
sub usage {
    print "\nPARA QUE SIRVE ESTO: Compara dos ficheros tabulados con el conteo de los miRNAs de dos muestras diferentes \n siempre que el coverage entre las dos sea de 10x o más. Ej: OrR1 Vrs OrR2\n";
    print "Input: 1)y 2) ficheros .ma.tab de las dos muestras que se quieren comparar 3) nombre del cromosoma\n";
    print "Output: 1)Fichero con miRNAs de la primera muestra 2)Fichero con miRNAs de la segunda muestra 3)Fichero con miRNAs comunes a las dos  \n";
    print "\n";
    exit(1);
}
if(scalar(@ARGV) == 0){
    usage();
}

open(TAB,"<",@ARGV[0]);
open(TAB2,"<",@ARGV[1]);
my $nombre_chr=@ARGV[2];
open(ENUNO,">","$nombre_chr.UNO.ma.tab") or die "No puedo abrir $nombre_chr.UNO.ma.tab.\n";
open(ENDOS,">","$nombre_chr.DOS.ma.tab") or die "No puedo abrir $nombre_chr.DOS.ma.tab.\n";
open(ENAMBAS,">","$nombre_chr.COMUNES.ma.tab") or die "No puedo abrir $nombre_chr.COMUNES.ma.tab.\n";


my @tab;
my $contador=0;
while  (my $tab_linea=<TAB>)
{
  	chomp($tab_linea);
	push (@tab, $tab_linea);
}

while (my $tab_linea2=<TAB2>)
{
	chomp($tab_linea2);
	push (@tab2,$tab_linea2);
}

close (TAB);
close (TAB2);

my $contador=1;
my $numero_miRNA=0;

for (my $a=0; $a<=$#tab2;)
{
	@split_tab_linea2=split(/\t/,$tab2[$a]);
	@cogerx=split("_x",$split_tab_linea2[0]);
	#$contador=$contador+1;
	#print join ("\n",@cogerx),"\n";
	#print "es2  ",$a,"\t",$split_tab_linea2[1],"\t",$split_tab_linea2[2],"\t",$split_tab_linea2[3],"\n";
	for (my $o=0; $o<=$#tab;$o++)
	{
        	@split_tab = split (/\t/,$tab[$o]);
		@cogerx2=split("x",$split_tab[0]);
		 #print "es1  ",$o,"\t",$split_tab[1],"\t",$split_tab[2],"\t",$split_tab[3],"\n";
		if(($split_tab_linea2[1] eq $split_tab[1]) && ($split_tab_linea2[2] eq $split_tab[2]) && ($split_tab_linea2[3] eq $split_tab[3]))
		{
			$contador=$contador+1;
			$numero_miRNA=$numero_miRNA+1;
			
			if (($cogerx[1]+$cogerx2[1])>=10)
			{
				print ENAMBAS "dme-mir-SG",$numero_miRNA,"_",$split_tab[3],$split_tab[1],"_",$split_tab[2],"_x",$cogerx[1]+$cogerx2[1],"\t",$split_tab[1],"\t",$split_tab[2],"\t",$split_tab[3],"\n";
			}
		}
	}
	
	if (($contador == 1) && ($cogerx[1] >=10) )
	{
		print ENDOS $tab2[$a],"\n";
	}
	$contador=1;
	$a++;
}

$contar=1;
for (my $e=0; $e<=$#tab;$e++)
{
	@split_viceversa = split (/\t/,$tab[$e]);
	@cogerx3=split("x",$split_viceversa[0]);
	for (my $u=0; $u<=$#tab2;$u++)
	{
		@split_viceversa2 = split (/\t/,$tab2[$u]);
		if(($split_viceversa[1] eq $split_viceversa2[1]) && ($split_viceversa[1] eq $split_viceversa2[1]) && ($split_viceversa[2] eq $split_viceversa2[2]) && ($split_viceversa[3] eq $split_viceversa2[3]))
                {
			$contar=$contar+1;
				
      

                }
	}
		
	if (($contar==1) && ($cogerx3[1] >=10))
	{
		print ENUNO $tab[$e],"\n";	
	}
	$contar=1;

}

close (ENUNO);
close (ENDOS);
close (ENAMBAS);
