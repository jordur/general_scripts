#!/usr/bin/perl -w
use File::ReadBackwards;

sub usage {
        print "\nCOMO SE USA: \n\n";
	print "OUTPUT: matriz_snps.txt\n\n";
        exit(1);
}

if(scalar(@ARGV) == 0){
    usage();
}

open(POSICIONES,$ARGV[$#ARGV]);
	
	while (my $posiciones=<POSICIONES>)
	{
		chomp($posiciones);
		push (@posiciones_referencia,$posiciones);
	}

close(POSICIONES);
#print  "------",join ("\t",@posiciones_referencia),"\n";


$str = "comp";
for($j = 0; $j <= $#ARGV-1; $j++)
{
	
	open (SALIDA,">",$ARGV[$j].$str);
	#print SALIDA join ("\n",@posiciones_referencia),"\n";
	open (FICHERO,"<",$ARGV[$j]) or die ("ack - $!");
#	open (SALIDA,">",$ARGV[$j].$str);
	# print "----ABRO FICHERO----\n";
	my $ultima_linea_fichero  = File::ReadBackwards->new($ARGV[$j])->readline;
	if (!$ultima_linea_fichero)
	{
	#	print "FICHERO NO EXISTE\n";
		for ($p=0; $p <= $#posiciones_referencia ; $p++)
		{
			print SALIDA "0\n";#$posiciones_referencia[$p],"--0\n";
		}	
	}
	else
	{
	chomp ($ultima_linea_fichero);
	my $val=0;
	while ($posiciones_muestra = <FICHERO>)
	{
		chomp ($posiciones_muestra);
		for ($mo=$val; $mo <= $#posiciones_referencia ; $mo++)
		{
			chomp ($posiciones_referencia[$mo]);
	#		print "LO QUE ENTRA-REF-MUESTRA--> ",$posiciones_referencia[$mo],"\t",$posiciones_muestra,"\n";
			if($posiciones_referencia[$mo] eq $posiciones_muestra)
			{
				if (($mo < $#posiciones_referencia) && ($posiciones_referencia[$mo] eq $ultima_linea_fichero))
				{
					$resto=$#posiciones_referencia-$mo;
					for ($quedan =0; $quedan <=$resto; $quedan++)
					{
	#					print "ULTIMAS-->",$posiciones_referencia[$mo+$quedan],"--0\n";
						print SALIDA "0\n";
					}
					last;
					
				}
				else
				{
	#				print "ES IGUAL-->",$posiciones_referencia[$mo],"--1\n";
					print SALIDA "1\n";
					$val=$mo+1;
					last;
				}
			}
			elsif ($posiciones_muestra > $posiciones_referencia[$mo])
			{
				$num=$posiciones_muestra-$posiciones_referencia[$mo]-1;
				for ($o=0; $o <= $num; $o++)
				{
	#				print "MUESTRA MAYOR----> ",$posiciones_referencia[$mo+$o],"--0\n";		
					print SALIDA "0\n";
				}
	#			print "MUESTRA MAYOR----> ",$posiciones_referencia[$mo+$num+1],"--1\n";
				print SALIDA "1\n";
				$val=$mo+$num+2;
				last;
			}
		}
	}
	}
	close(FICHERO);
	close(SALIDA);
		
}










