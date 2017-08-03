#!/usr/bin/perl
sub usage
{
    print "\nPARA QUE VALE:\n";
    print "Input: Fichero con los datos de genotipado de Anita B \n ";
    print "Output: Fichero con los SNPs separados \n";
    exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}


open(GRUPOS,"<",@ARGV[0]);
while (my $linea = <GRUPOS>)
{
	chomp ($linea);
	push (@array, $linea);
}
close (GRUPOS);

my $valor = `awk '{print \$1}' @ARGV[0] | sort -u  2>&1`;
my $count=0;
my @split_fichero= split ("\n",$valor);
my @split_split= split ("\t",$split_fichero[0]),"\n";



sub imprimo ()
{
	for ($u=0; $u<=$#split_fichero; $u++)
	{
		
		@split_split= split ("\t",$split_fichero[$u]),"\n";
		open (OUT, ">",'ini'.$split_split[0].'.tmp')or die "No puedo abrir $split_split[0].tmp \n";
		my $fichero_tmp='ini'.$split_split[0].'.tmp';
		for ($a=0; $a<= $#array ; $a++)
		{
			@split_linea_anterior= split ("\t", $array[$a]);
			if ($split_linea_anterior[0] eq $split_split[0])
			{
				print OUT $array[$a],"\n";
			}
					
		}
		close (OUT);
	}


}

$countG=0;
$countS=0;
$countN=0;
$numG=0;
$numS=0;
$numN=0;

sub individuals ()
{
        $name=`ls -l | grep ini | sed 's/.*ini/ini/' 2>&1`;
	chomp($name);
	@split_name= split ("\n",$name);
	$genotypes= `awk '{print \$3}' @ARGV[0] | sort -u  2>&1`;
	chomp($genotypes);
	@split_genotypes =  split ("\n", $genotypes);
	for ($e=0; $e<=$#split_name; $e++)
	{
		print "\#\#\#\#\#\#\#\#\n$split_name[$e]\t",$e+1,"\t\t\n";
		my @tic;
	        open(TMP,"<",$split_name[$e]);
        	while (my $toe=<TMP>)
        	{
        		chomp($toe);
	        	push (@tic,$toe);
        	}
        	close(TMP);
		for ($v=0; $v<=$#tic; $v++)
		{
			@tic_split=split ("\t", $tic[$v]);
			$tic_split[1]=~s/\x0D//;
			$nombre=substr ($tic_split[1],0,1);
			if ($nombre eq "S")
			{
				$numS=$numS+1;
			}
			elsif ($nombre eq "N")
			{
				$numN=$numN+1;
			}
			elsif ($nombre eq "G")
			{
				$numG=$numG+1;
			}
		}

		for (my $i=0; $i<=$#split_genotypes ; $i++)
		{
			$split_genotypes[$i]=~s/\x0D//;
			for (my $uno=0; $uno<=$#tic ; $uno++)
 	             	{
               		        chomp ($tic[$uno]);
				@split_tic = split ("\t", $tic[$uno]);
				$split_tic[1]=~s/\x0D//;
	                        $tipo_individuo=substr ($split_tic[1],0,1);
				chomp($split_tic[2]);
				$split_tic[2]=~s/\x0D//;
				
				if (($split_genotypes[$i] eq $split_tic[2]))
				{
					
				        if ($tipo_individuo eq "S")
        				{
				                $countS=$countS+1;
        				}
				        elsif ($tipo_individuo eq "N")
        				{
				                $countN=$countN+1;
        				}
				        elsif ($tipo_individuo eq "G")
        				{
				                $countG=$countG+1;
        				}
				}
			}
					print "genotipo_buscado_$split_genotypes[$i]\tindividuos_S_totales_en_ese_SNP\t$numS\t individuos_con_ese_genotipo_$split_genotypes[$i]\t$countS\n";
					print "genotipo_buscado_$split_genotypes[$i]\tindividuos_N_totales_en_ese_SNP\t$numN\t individuos_con_ese_genotipo_$split_genotypes[$i]\t$countN\n";
					print "genotipo_buscado_$split_genotypes[$i]\tindividuos_G_totales_en_ese_SNP\t$numG\t individuos_con_ese_genotipo_$split_genotypes[$i]\t$countG\n";
			print "------","\n";
			$countS=0;
			$countN=0;
			$countG=0;
		}
		$numS=0;
		$numN=0;
		$numG=0;
	}
}


imprimo();
individuals();

