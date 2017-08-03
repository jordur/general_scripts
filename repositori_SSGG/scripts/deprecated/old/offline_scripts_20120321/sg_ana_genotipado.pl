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
$countC=0;
$countH=0;
$countN=0;

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
		print "\#\#\#\#\#\#\#\#\nFICHERO$e\n";
		my @tic;
	        open(TMP,"<",$split_name[$e]);
        	while (my $toe=<TMP>)
        	{
        		chomp($toe);
	        	push (@tic,$toe);
        	}
        	close(TMP);
		for (my $i=0; $i<=$#split_genotypes ; $i++)
		{
			print "genotipo_buscado ", $split_genotypes[$i],"\n";
			$split_genotypes[$i]=~s/\x0D//;
			for (my $uno=0; $uno<=$#tic ; $uno++)
 	             	{
               		        chomp ($tic[$uno]);
				@split_tic = split ("\t", $tic[$uno]);
	                        $tipo_individuo=substr ($split_tic[1],0,1);
				#print "genotipo_buscar $split_genotypes[$i] genotipo_individuo $split_tic[2]\n";
#				print "genotipo_buscarrrrrrrrrrrrrrrrrrrrrrrrrrrrr  $split_genotypes[$i] genotipo_individuo $split_tic[2]\n";
#$split_genotypes[$genos] ,"\tgenotipo_individuo", $split_tic[2],"\n";
				chomp($split_tic[2]);
				$split_tic[2]=~s/\x0D//;
				#print $split_tic[2],"\n";
				# print "genotipo_buscarrrrrrrrr $split_genotypes[$i] genotipo_individuo $split_tic[2]\n";	
				if (($tipo_individuo eq "C") && ($split_genotypes[$i] eq $split_tic[2]))
				{
					$countC=$countC+1;
					#push (@C, $split_genotypes[$genos]);
					
					#print "individuosC\t genotipos_iguales $split_tic[2] = $countC \n";
						
				}
			
				elsif (($tipo_individuo eq "H") && ($split_genotypes[$i] eq $split_tic[2]))
				{
					$countH=$countH+1;
					#print "individuosH\t genotipos iguales $split_tic[2] = $countH\n";
				}
				elsif  (($tipo_individuo eq "N") && ($split_genotypes[$i] eq $split_tic[2]))
				{
					$countN=$countN+1;
					#print "individuosN\t genotipos iguales $split_tic[2]= $countN\n";
				}
				elsif  (($tipo_individuo eq "G") && ($split_genotypes[$i] eq $split_tic[2]))
				{
					$countG=$countG+1;
					#print "individuosG\t genotipos iguales $split_tic[2]= $countG\n";
				}
			}
			if ($countC!=0)
			{
				print "individuosC\t genotipos_iguales $split_genotypes[$i] = $countC \n";
			}
			elsif ($countH!=0)
			{
				print "individuosH\t genotipos_iguales $split_genotypes[$i] = $countH\n";
			}
			elsif ($countN!=0)
			{
				print "individuosN\t genotipos_iguales $split_genotypes[$i] = $countN\n";
			}
			elsif ($countG!=0)
			{
				print "individuosG\t genotipos_iguales $split_genotypes[$i] = $countG\n";
			}
			print "------","\n";
			$countC=0;
			$countH=0;
			$countN=0;
			$countG=0;
		}
	}
}



imprimo();
individuals();

