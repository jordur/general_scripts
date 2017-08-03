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
        	$signo= substr($split2[1],0,1);
        	$valor1= substr($split2[1],1,length ($split2[1]));
		print PARSEADO "dme-mir-grupo_",$contador,"_",$split2[1],"_",$split2[2],"\t",$valor1,"\t",$split2[2],"\t",$signo,"\t";
		#print MIMADUROS "dme-mir-grupo_",$contador,"_",$split2[1],"_",$split2[2],"\t",$valor1,"\t",$split2[2],"\t",$signo,"\n";		
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
	elsif ($inicio = substr ($split_linea2[0],0,7) eq "dme-mir")
	{
		@split3=split ("_",$split_linea2[0]);
                $signo2= substr($split3[1],0,1);
                $valor2= substr($split3[1],1,length ($split3[1]));
		print PARSEADO $split_linea2[0],",";
		print GRUPOS $split_linea2[0],"\t",$valor2,"\t",$split3[2],"\t",$signo2,"\n";
		@cov = split ("_x",$split_linea[0]);
		$coverage = $coverage+$cov[1];
	}
	elsif ($array[$i] =~ /^$/)
	{
		$uno=$uno+1;
		if ($uno == 2 || $i == $#array)
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
