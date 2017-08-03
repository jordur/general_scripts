#!/usr/bin/perl
sub usage
{	
	print "\nPARA QUE SIRVE: Parsea la salida de miRDeep\n";
    	print "Input: 1) Fichero *.precursores que se obtiene con miRDeepl.pl 2) nombre del cromosoma\n";
	print "Output: 1) Fichero tabulado con la información de los miRNAs predichos 2) Fichero con los miRNAs incluidos en cada miRNA detectado para hacer el recuento real de los miRNAs maduros nuevos\n\n\n";
    exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}

open(MIRDEEP,"<",@ARGV[0]);
my $nombre_chr=@ARGV[1];
open(PARSEADO,">","$nombre_chr.predicciones.parseado") or die "No puedo abrir $nombre_chr.predicciones.parseado\n";
open(GRUPOS,">","$nombre_chr.predicciones.grupos") or die "No puedo abrir $nombre_chr.predicciones.grupos\n";
#open(MIMADUROS,">","$nombre_chr.predicciones.miRNAs_maduros.tab") or die "No puedo abrir $nombre_chr.predicciones.miRNAs_maduros.tab\n";
my $inicio;
my @array;
while  (my $linea=<MIRDEEP>)
{
        chomp($linea);
 	push (@array,$linea);       
}

close (MIRDEEP);

$coverage=0;
$contador =1;
$uno=0;
print PARSEADO "#score\t","loop_beg\t","loop_end\t","loop_seq\t","mature_beg\t","mature_end\t","mature_query\t","mature_chr_beg\t","mature_chr_end\t","mature_strand\t","mature_seq\t","pre_seq\t","chr\t","pri_seq\t","star_beg\t","star_end\t","star_seq\t","groups_in\t","coverage\n";
for (my $i=0; $i<=$#array ; $i++)
{
	@split_linea= split ("  \t",$array[$i]);
	@split_linea2 = split ("\t",$array[$i]);
	if ($split_linea2[0] eq "score")
	{
		print PARSEADO $split_linea2[1],"\t";
	}
	elsif ($split_linea[0] eq "loop_beg")
	{
		print PARSEADO $split_linea[1],"\t";		
	}
	elsif ($split_linea[0] eq "loop_end")
	{
		print PARSEADO $split_linea[1],"\t";		
	}
	elsif ($split_linea[0] eq "loop_seq")
	{
		print PARSEADO $split_linea[1],"\t";		
	}
	elsif ($split_linea[0] eq "mature_beg")
	{
		print PARSEADO $split_linea[1],"\t";		
	}
	elsif ($split_linea[0] eq "mature_end")
	{
		print PARSEADO $split_linea[1],"\t";		
	}
	elsif ($split_linea[0] eq "mature_query")
	{
		@split2=split ("_",$split_linea[1]);
        	$signo_mature= substr($split2[4],0,1);
        	if ($signo_mature eq "-")
		{
			$valor1= substr($split2[4],1,length ($split2[4]));
			print PARSEADO "dme-mir-grupo_",$contador,"_",$split2[4],"_",$split2[5],"\t",$valor1,"\t",$split2[5],"\t",$signo_mature,"\t";
		}
		else
		{
			print PARSEADO "dme-mir-grupo_",$contador,"_",$split2[4],"_",$split2[5],"\t",$split2[4],"\t",$split2[5],"\t","+","\t";
		}
	}
	elsif ($split_linea[0] eq "mature_seq")
	{
		print PARSEADO $split_linea[1],"\t";		
	}
	elsif ($split_linea[0] eq "pre_seq")
	{
		print PARSEADO $split_linea[1],"\t";		
	}
	elsif ($split_linea[0] eq "pri_id")
        {
                @split_pri= split ("_",$split_linea[1]);
		print PARSEADO $split_pri[0],"\t";
        }
	elsif ($split_linea[0] eq "pri_seq")
	{
		print PARSEADO $split_linea[1],"\t";		
	}
	elsif ($split_linea[0] eq "star_beg")
	{
		print PARSEADO $split_linea[1],"\t";
	}
	elsif ($split_linea[0] eq "star_end")
	{
		print PARSEADO $split_linea[1],"\t";		
	}
	elsif ($split_linea[0] eq "star_seq")
	{
		print PARSEADO $split_linea[1],"\t";		
	}
	elsif ($array[$i] =~ /Plus/)
	{
		@split3=split ("_",$split_linea2[0]);
		#print "---------->", $split3[4],"\n";
		$signo= substr($split3[4],0,1);
	        $valor_inicio= substr($split3[4],1,length ($split3[4]));
		if ($signo_mature eq "-")
                {
			#$valor_inicio= substr($split3[4],1,length ($split3[4]));
                	print PARSEADO $split_linea2[0],",";
               		print GRUPOS $split_linea2[0],"\t",$valor_inicio,"\t",$split3[4],"\t",$signo,"\n";
               		@cov = split ("_x",$split_linea[0]);
               		$coverage = $coverage+$cov[1];
		}
		else
		{
			
			print PARSEADO $split_linea2[0],",";
			print GRUPOS $split_linea2[0],"\t",$valor_inicio,"\t",$split3[4],"\t","+","\n";
			@cov2 = split ("_x",$split_linea[0]);
			$coverage = $coverage+$cov2[1];
		}
	}
	elsif ($array[$i] =~ /^$/)
	{
		$uno=$uno+1; #Si hay una línea en blanco
		if ($uno == 2 || $i == $#array)#Si hay una línea en blanco o hemos llegado al final del fichero
		{
			print PARSEADO "\t",$coverage,"\n";
			$coverage=0;
			$contador++;
			$uno=0;
		}

	}
}	
close (PARSEADO);
close (GRUPOS);
#close (MIMADUROS);
