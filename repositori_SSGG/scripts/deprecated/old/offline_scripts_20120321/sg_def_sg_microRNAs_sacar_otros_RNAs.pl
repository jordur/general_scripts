#!/usr/bin/perl
sub usage
{
    print "\nPARA QUE SIRVE: Saca los grupos de lecturas que coinciden con ncRNAs o cualquier otro fichero con formato nombre<TAB>chr<TAB>inicio<TAB>final<TAB>orientacion\n";
    print " Input: 1) fichero tabulado del conteo de miRNAs 2) base de datos tabulada ej: mir2<TAB>chr<TAB>inicio<TAB>final<TAB>orientacion  3) nombre del cromosoma y tipo de molecula que buscajos Ej: CHR1.ncRNAs\n";
    exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}

open(CONTEO,"<",@ARGV[0]);
open(MIRBASE,"<",@ARGV[1]);
my $nombre_chr=@ARGV[2];
open(SENSE,">","$nombre_chr.conocidos.sense.ma.tab") or die "No puedo abrir $nombre_chr.conocidos.sense.ma.tab.\n";
open(ANTISENSE,">","$nombre_chr.conocidos.antisense.ma.tab") or die "No puedo abrir $nombre_chr.conocidos.antisense.ma.tab.\n";
open(NUEVOS,">","$nombre_chr.nuevos.ma.tab") or die "No puedo abrir $nombre_chr.nuevos.ma.tab.\n";


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
my $contador_antisense=0;
for (my $a=0; $a<=$#mirbase2 ;$a++)
{
	@mirbase= split ("\t",$mirbase2[$a]);	
	for (my $i=0;$i<=$#conteo;$i++)
	{
		@conteo_split= split ("\t",$conteo[$i]);
		@conteo_split2 = split ("_x",$conteo_split[0]);
		#print $mirbase[4],"\n";#" igual ",$mirbase[4],"\n";
		if ((($conteo_split[1]>=$mirbase[2]) && ($conteo_split[2]<=$mirbase[3]) && ($conteo_split[3] eq $mirbase[4])))
		{
			print SENSE $conteo_split[0],"\t",$conteo_split[1],"\t",$conteo_split[2],"\t",$conteo_split2[1],"\n"; #mirbase[0]."_".$mirbase[1]."_".$mirbase[4]."_".$mirbase[2]."_".$mirbase[3],"\t",$conteo_split[1],"\t",$conteo_split[2],"\t",$conteo_split2[1],"\n";		
			$contador=$contador+$conteo_split2[1];
			#print $conteo[$i],"\n";
			push (@comparar,$conteo[$i]);
			
		}
		elsif ((($conteo_split[1]>=$mirbase[2]) && ($conteo_split[2]<=$mirbase[3]) && ($conteo_split[3] ne $mirbase[4])))
                {
                        print ANTISENSE $conteo_split[0],"\t",$conteo_split[1],"\t",$conteo_split[2],"\t",$conteo_split2[1],"\n";
                        $contador_antisense=$contador_antisense+$conteo_split2[1];
			push (@comparar,$conteo[$i]);
                }

	}
	if ($contador_antisense>0)
	{
		print ANTISENSE $mirbase[0]."_".$mirbase[1]."_".$mirbase[4]."_".$mirbase[2]."_".$mirbase[3],"=",$contador_antisense,"\n";
		if ($contador >0)
		{
			print SENSE $mirbase[0]."_".$mirbase[1]."_".$mirbase[4]."_".$mirbase[2]."_".$mirbase[3],"=",$contador,"\n";
		}
	}
	elsif ($contador >0)
	{
		print SENSE $mirbase[0]."_".$mirbase[1]."_".$mirbase[4]."_".$mirbase[2]."_".$mirbase[3],"=",$contador,"\n";
	}
	$contador =0;
	$contador_antisense=0;

}

#print @comparar,"\n";
my $siEsta=0;
for (my $u;$u<=$#conteo;$u++)
{
	#print $conteo[$u],"\n";
	
	for (my $mo;$mo<=$#comparar;$mo++)
	{
		if ($conteo[$u] eq $comparar[$mo])
		{
			$siEsta=$siEsta+1;
			#print "\n";	
		}
	}
	#print $siEsta,"\n";
	if ($siEsta == 0)
	{
		print NUEVOS $conteo[$u],"\n";
	}
#	print $si_esta,"\n";
	$siEsta=0;

}

close (SENSE);
close (ANTISENSE);
close (NUEVOS);
